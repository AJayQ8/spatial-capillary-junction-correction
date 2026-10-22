#!/usr/bin/env python3
"""Qualify the frozen H24/N20 one-step CLI recipes and input matching."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def raw_sha(array: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output-dir", type=Path,
                    default=ROOT / "build/qualification/b2_step_cli")
    args = ap.parse_args()
    outdir = args.output_dir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    if any(outdir.iterdir()):
        raise FileExistsError(f"CLI qualification directory must be empty: {outdir}")

    subprocess.run([sys.executable, str(ROOT / "scripts/build_kernels.py")], check=True)
    n20_recipe = json.loads((ROOT / "configs/b2_n20_warm.json").read_text())
    if (n20_recipe["dx_W"], n20_recipe["eta_W"], n20_recipe["dt_W2"]) != (0.005, 0.1, 1.0e-7):
        raise AssertionError("N20 frozen recipe no longer matches its retained grid, eta and timestep")

    cases = [
        ("h24_d9_halo", "configs/b2_h24.json", "inputs/h24_cold_initial_field.npz",
         "d9", "halo", (3, 2400, 100)),
        ("n20_d5_halo", "configs/b2_n20_warm.json", "inputs/n20_warm_time1_phi.npy",
         "d5", "halo", (3, 800, 200)),
    ]
    records = []
    for label, recipe_rel, input_rel, stencil, route, shape in cases:
        output = outdir / f"{label}.npz"
        completed = subprocess.run([
            sys.executable, str(ROOT / "scripts/run_b2_step.py"),
            "--recipe", str(ROOT / recipe_rel), "--input", str(ROOT / input_rel),
            "--stencil", stencil, "--route", route, "--output", str(output),
        ], check=True, capture_output=True, text=True)
        report = json.loads(completed.stdout)
        with np.load(output, allow_pickle=False) as archive:
            result = np.ascontiguousarray(archive["phi"], dtype=np.float64)
        if tuple(result.shape) != shape or raw_sha(result) != report["output_raw_sha256"]:
            raise AssertionError(f"one-step CLI result mismatch for {label}")
        records.append({"case": label, "recipe": recipe_rel, "input": input_rel,
                        "input_raw_sha256": report["input_raw_sha256"],
                        "output_raw_sha256": report["output_raw_sha256"],
                        "shape": list(result.shape)})

    bad_output = outdir / "wrong_recipe.npz"
    rejected = subprocess.run([
        sys.executable, str(ROOT / "scripts/run_b2_step.py"),
        "--recipe", str(ROOT / "configs/b2_h24.json"),
        "--input", str(ROOT / "inputs/n20_warm_time1_phi.npy"),
        "--stencil", "d5", "--route", "halo", "--output", str(bad_output),
    ], check=False, capture_output=True, text=True)
    rejection = rejected.stderr.lower()
    clear_mismatch = ("requires an input field with shape" in rejection
                      or "raw hash" in rejection)
    if rejected.returncode == 0 or bad_output.exists() or not clear_mismatch:
        raise AssertionError("H24 recipe was not clearly rejected for the N20 warm field")

    receipt = {
        "schema": "b2-one-step-cli-qualification-v1",
        "status": "PASS",
        "cases": records,
        "wrong_h24_recipe_rejected_for_n20": True,
        "scope": "Two one-step CLI checks and one expected recipe/field mismatch rejection; no trajectory run.",
    }
    (outdir / "qualification.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"status": "PASS", "cases": len(records),
                      "wrong_recipe_rejected": True, "output": "qualification.json"}, indent=2))


if __name__ == "__main__":
    main()
