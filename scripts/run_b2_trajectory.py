#!/usr/bin/env python3
"""Small H24 B2 time loop with restartable local-memory controller state.

This entrypoint exposes the released kernels and checkpoint state. It does not
reproduce the paper's full observer suite or the N20 moving-window continuation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src/python"))
from b2_api import step as b2_step


def sha_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def array_sha(a: np.ndarray) -> str:
    return sha_bytes(np.ascontiguousarray(a).tobytes())


def load_field(path: Path) -> np.ndarray:
    if path.suffix == ".npy":
        field = np.load(path, allow_pickle=False)
    else:
        with np.load(path, allow_pickle=False) as data:
            field = data["phi"]
    field = np.ascontiguousarray(field, dtype=np.float64)
    if field.ndim != 3 or field.shape[0] != 3:
        raise ValueError("input field must have shape (3, ny, nx)")
    return field


def stats(phi: np.ndarray, recipe: dict, step_no: int, factor: float) -> dict:
    total = np.sum(phi, axis=0)
    return {
        "step": int(step_no),
        "physical_time_W2": float(step_no * recipe["dt_W2"] / factor),
        "raw_time_W2": float(step_no * recipe["dt_W2"]),
        "field_raw_sha256": array_sha(phi),
        "phi_min": float(np.min(phi)),
        "phi_max": float(np.max(phi)),
        "partition_linf": float(np.max(np.abs(total - 1.0))),
        "phase_integral_W2": (np.sum(phi, axis=(1, 2)) * recipe["dx_W"] ** 2).tolist(),
    }


def atomic_checkpoint(path: Path, phi: np.ndarray, suppression: np.ndarray,
                      metadata: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    with tmp.open("wb") as stream:
        np.savez_compressed(
            stream,
            phi=np.ascontiguousarray(phi, dtype=np.float64),
            suppression=np.ascontiguousarray(suppression, dtype=np.uint8),
            metadata=np.asarray(json.dumps(metadata, sort_keys=True)),
        )
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(tmp, path)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    source = ap.add_mutually_exclusive_group(required=True)
    source.add_argument("--input", type=Path, help="H24 .npz (phi) or .npy initial field")
    source.add_argument("--resume", type=Path, help="checkpoint emitted by this script")
    ap.add_argument("--recipe", type=Path, default=ROOT / "configs/b2_h24.json")
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--stencil", choices=("d5", "d9"), default="d9")
    ap.add_argument("--route", choices=("native", "halo", "local"), default="halo")
    ap.add_argument("--steps", type=int, required=True, help="additional Euler steps to run")
    ap.add_argument("--observe-every", type=int, default=100)
    ap.add_argument("--checkpoint-every", type=int, default=1000)
    ap.add_argument("--library", type=Path, help="Override the built shared library")
    args = ap.parse_args()

    if args.steps < 1 or args.observe_every < 1 or args.checkpoint_every < 1:
        ap.error("steps and output intervals must be positive")
    recipe_path = args.recipe.resolve()
    recipe_bytes = recipe_path.read_bytes()
    recipe = json.loads(recipe_bytes)
    recipe_sha = sha_bytes(recipe_bytes)
    if recipe.get("schema") != "b2-h24-single-step-recipe-v1":
        ap.error("the trajectory entrypoint is limited to the fixed H24 recipe")
    if args.route == "local" and args.stencil != "d9":
        ap.error("the retained local-memory trajectory uses D9 only")

    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    observation_path = output_dir / "observations.jsonl"
    suffix = ".dylib" if sys.platform == "darwin" else ".so"
    kernel_name = "local_memory" if args.route == "local" else args.route
    library = args.library or ROOT / "build/kernels" / f"b2_{args.stencil}_{kernel_name}{suffix}"
    if not library.is_file():
        raise FileNotFoundError(f"kernel library not found: {library}; run scripts/build_kernels.py first")

    factor = float(recipe.get("planar_factor_a", 1.0))
    cadence = float(recipe.get("local_controller_cadence_W2", 0.0002))
    cursor = 1
    step_no = 0
    previous_suppression = np.zeros((2400, 100), dtype=np.uint8)
    controller_log = []
    if args.input is not None:
        if args.route == "local":
            from local_memory_logic import next_suppression
        phi = load_field(args.input)
        if phi.shape != (3, 2400, 100):
            ap.error("only the fixed H24 domain (3, 2400, 100) is supported here")
        initial_hash = array_sha(phi)
        expected_hash = recipe.get("seed_sha256_raw")
        if expected_hash and initial_hash != expected_hash:
            ap.error("initial field raw hash differs from the selected frozen recipe")
        if args.route == "local":
            previous_suppression, meta = next_suppression(
                phi, None, int(recipe["radius_cells"]), evaluate_overlap=True)
            previous_suppression = np.ascontiguousarray(previous_suppression, dtype=np.uint8)
            controller_log.append({"step": 0, "event": "initial_evaluation", **{
                key: value for key, value in meta.items() if key != "overlap_centers"}})
        if any(output_dir.iterdir()):
            raise FileExistsError(f"fresh run requires an empty output directory: {output_dir}")
        run_meta = {
            "schema": "b2-h24-trajectory-checkpoint-v1",
            "recipe_sha256": recipe_sha,
            "route": args.route,
            "stencil": args.stencil,
            "initial_field_raw_sha256": initial_hash,
            "global_step": 0,
            "next_controller_index": cursor,
            "controller_events": controller_log,
        }
        with observation_path.open("w", encoding="utf-8") as stream:
            row = stats(phi, recipe, 0, factor)
            row.update({"route": args.route, "stencil": args.stencil,
                        "suppressed_support_nodes": int(previous_suppression.sum()) if args.route == "local" else None})
            stream.write(json.dumps(row, sort_keys=True) + "\n")
    else:
        with np.load(args.resume, allow_pickle=False) as data:
            phi = np.ascontiguousarray(data["phi"], dtype=np.float64)
            previous_suppression = np.ascontiguousarray(data["suppression"], dtype=np.uint8)
            run_meta = json.loads(str(data["metadata"].item()))
        if run_meta.get("schema") != "b2-h24-trajectory-checkpoint-v1":
            ap.error("checkpoint schema is not supported")
        if run_meta.get("recipe_sha256") != recipe_sha or run_meta.get("route") != args.route or run_meta.get("stencil") != args.stencil:
            ap.error("checkpoint recipe/route/stencil does not match the requested run")
        if phi.shape != (3, 2400, 100) or previous_suppression.shape != (2400, 100):
            ap.error("checkpoint dimensions do not match the H24 grid")
        step_no = int(run_meta["global_step"])
        cursor = int(run_meta["next_controller_index"])
        if not observation_path.exists():
            raise FileNotFoundError(f"resume requires the original observations file in {output_dir}")

    for _ in range(args.steps):
        suppression = previous_suppression if args.route == "local" else None
        phi, kernel_summary = b2_step(
            library, phi, recipe, stencil=args.stencil, route=args.route,
            suppression=suppression,
        )
        step_no += 1
        event_rows = []
        if args.route == "local":
            from local_memory_logic import next_suppression
            # At the first integration state on/after each physical target,
            # evaluate new overlaps before the following Euler update.
            target_step = math.ceil(cursor * cadence * factor / recipe["dt_W2"] - 1e-12)
            evaluate = step_no >= target_step
            previous_suppression, meta = next_suppression(
                phi, previous_suppression, int(recipe["radius_cells"]),
                evaluate_overlap=evaluate,
            )
            previous_suppression = np.ascontiguousarray(previous_suppression, dtype=np.uint8)
            if evaluate:
                event = {"step": step_no, "event": "controller_evaluation", "target_index": cursor,
                         "target_time_W2": cursor * cadence, "actual_time_W2": step_no * recipe["dt_W2"] / factor}
                event.update({key: value for key, value in meta.items() if key != "overlap_centers"})
                run_meta.setdefault("controller_events", []).append(event)
                event_rows.append(event)
                while step_no >= target_step:
                    cursor += 1
                    target_step = math.ceil(cursor * cadence * factor / recipe["dt_W2"] - 1e-12)
            run_meta["next_controller_index"] = cursor
        run_meta["global_step"] = step_no

        # The final state is always observed; intermediate rows follow the
        # requested cadence. A resumed run appends to the same JSONL file.
        if step_no % args.observe_every == 0 or _ == args.steps - 1:
            row = stats(phi, recipe, step_no, factor)
            row.update({
                "route": args.route,
                "stencil": args.stencil,
                "suppressed_support_nodes": int(previous_suppression.sum()) if args.route == "local" else None,
                "controller_evaluations": event_rows,
                "kernel_summary": kernel_summary,
            })
            with observation_path.open("a", encoding="utf-8") as stream:
                stream.write(json.dumps(row, sort_keys=True) + "\n")

        if step_no % args.checkpoint_every == 0:
            atomic_checkpoint(output_dir / f"checkpoint_step{step_no:09d}.npz", phi, previous_suppression, run_meta)

    final_checkpoint = output_dir / "checkpoint_latest.npz"
    atomic_checkpoint(final_checkpoint, phi, previous_suppression, run_meta)
    print(json.dumps({
        "status": "completed",
        "route": args.route,
        "stencil": args.stencil,
        "global_step": step_no,
        "next_controller_index": cursor,
        "field_raw_sha256": array_sha(phi),
        "observations": str(observation_path),
        "checkpoint": str(final_checkpoint),
    }, indent=2))


if __name__ == "__main__":
    main()
