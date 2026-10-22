#!/usr/bin/env python3
"""Bounded Common-launched qualification for the released numerical sources."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def sha(a: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()


def run(command: list[str]) -> None:
    subprocess.run(command, check=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--openphase-root", type=Path, required=True)
    ap.add_argument("--fftw-prefix", type=Path)
    ap.add_argument("--output-dir", type=Path, default=ROOT / "build/qualification/release_smoke")
    args = ap.parse_args()
    outdir = args.output_dir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    if any(outdir.iterdir()):
        raise FileExistsError(f"qualification output directory must be empty: {outdir}")

    run([sys.executable, str(ROOT / "scripts/build_kernels.py")])
    kernel_report = outdir / "kernel_qualification.json"
    run([sys.executable, str(ROOT / "tests/qualify_kernels.py"), "--output", str(kernel_report)])

    cli_dir = outdir / "b2_step_cli"
    run([sys.executable, str(ROOT / "tests/qualify_b2_step_cli.py"),
         "--output-dir", str(cli_dir)])
    cli_report = cli_dir / "qualification.json"

    openphase_report = outdir / "openphase"
    openphase_command = [sys.executable, str(ROOT / "tests/run_openphase_smoke.py"),
                         "--openphase-root", str(args.openphase_root.resolve()),
                         "--output-dir", str(openphase_report)]
    if args.fftw_prefix:
        openphase_command += ["--fftw-prefix", str(args.fftw_prefix.resolve())]
    run(openphase_command)

    runner = str(ROOT / "scripts/run_b2_trajectory.py")
    seed = str(ROOT / "inputs/h24_cold_initial_field.npz")
    recipe = str(ROOT / "configs/b2_h24.json")
    split_dir = outdir / "b2_split"
    common = [sys.executable, runner, "--recipe", recipe, "--stencil", "d9",
              "--route", "local", "--output-dir", str(split_dir),
              "--observe-every", "255", "--checkpoint-every", "510"]
    run(common + ["--input", seed, "--steps", "510"])
    run(common + ["--resume", str(split_dir / "checkpoint_latest.npz"), "--steps", "5"])

    continuous_dir = outdir / "b2_continuous"
    continuous = [sys.executable, runner, "--recipe", recipe, "--stencil", "d9",
                  "--route", "local", "--output-dir", str(continuous_dir),
                  "--observe-every", "255", "--checkpoint-every", "515",
                  "--input", seed, "--steps", "515"]
    run(continuous)

    with np.load(split_dir / "checkpoint_latest.npz", allow_pickle=False) as split:
        split_phi = np.asarray(split["phi"])
        split_mask = np.asarray(split["suppression"])
        split_meta = json.loads(str(split["metadata"].item()))
    with np.load(continuous_dir / "checkpoint_latest.npz", allow_pickle=False) as full:
        full_phi = np.asarray(full["phi"])
        full_mask = np.asarray(full["suppression"])
        full_meta = json.loads(str(full["metadata"].item()))
    if split_meta["global_step"] != 515 or split_meta["next_controller_index"] != 2:
        raise AssertionError("split checkpoint did not preserve the first controller cursor")
    if full_meta["global_step"] != 515 or full_meta["next_controller_index"] != 2:
        raise AssertionError("continuous checkpoint did not advance the first controller cursor")
    if split_meta["controller_events"][-1].get("step") != 515:
        raise AssertionError("controller evaluation is absent at the first scheduled state")
    if not np.array_equal(split_phi, full_phi) or not np.array_equal(split_mask, full_mask):
        raise AssertionError("split/resumed local-memory path differs from continuous path")

    result = {
        "schema": "bounded-portable-release-qualification-v1",
        "status": "PASS",
        "scope": "Seven source kernel smoke; four two-step OpenPhase caller checks; exact 515-step local-memory restart check. No full-horizon production run.",
        "kernel_qualification": {"path": "kernel_qualification.json",
                                 "sha256": hashlib.sha256(kernel_report.read_bytes()).hexdigest()},
        "b2_step_cli": {
            "path": "b2_step_cli/qualification.json",
            "sha256": hashlib.sha256(cli_report.read_bytes()).hexdigest(),
        },
        "openphase_qualification": {"path": "openphase/qualification.json",
                                    "sha256": hashlib.sha256((openphase_report / "qualification.json").read_bytes()).hexdigest()},
        "b2_local_memory_restart": {
            "steps": 515,
            "first_controller_step": 515,
            "next_controller_index": 2,
            "split_field_raw_sha256": sha(split_phi),
            "continuous_field_raw_sha256": sha(full_phi),
            "suppression_masks_equal": True,
            "field_states_equal": True,
        },
    }
    result_path = outdir / "qualification.json"
    result_path.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": "PASS", "qualification": "qualification.json",
                      "b2_field_raw_sha256": sha(split_phi)}, indent=2))


if __name__ == "__main__":
    main()
