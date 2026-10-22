#!/usr/bin/env python3
"""Run one retained H24 OpenPhase route with its pinned caller and recipe."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import subprocess

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ROUTES = {
    "native-full": ("openphase_h24_native_full_dt", "configs/openphase_h24_native_full_dt.opi",
                    "configs/openphase_h24_full_dt.json", "inputs/schedule_full_dt.tsv", "src/openphase_h24_native_v004.cpp", False),
    "halo-full": ("openphase_h24_halo", "configs/openphase_h24_halo_full_dt.opi",
                  "configs/openphase_h24_full_dt.json", "inputs/schedule_full_dt.tsv", "src/openphase_h24_halo_v001.cpp", True),
    "native-half": ("openphase_h24_native_halfdt", "configs/openphase_h24_native_halfdt.opi",
                    "configs/openphase_h24_halfdt.json", "inputs/schedule_half_dt.tsv", "src/openphase_h24_native_halfdt_v001.cpp", False),
    "halo-half": ("openphase_h24_halo", "configs/openphase_h24_halo_halfdt.opi",
                  "configs/openphase_h24_halfdt.json", "inputs/schedule_half_dt.tsv", "src/openphase_h24_halo_v001.cpp", True),
}


def sha_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha_file(path: Path) -> str:
    return sha_bytes(path.read_bytes())


def replace_steps(opi: str, steps: int) -> str:
    if steps < 1:
        raise ValueError("steps must be positive")
    keys = {"$nSteps", "$FTime", "$STime", "$tRstrt"}
    output = []
    for line in opi.splitlines():
        if any(line.startswith(key) for key in keys):
            head, sep, tail = line.partition(":")
            if not sep:
                raise ValueError(f"malformed OPI field: {line}")
            value, *rest = tail.strip().split(maxsplit=1)
            suffix = (" " + rest[0]) if rest else ""
            line = head + ": " + str(steps) + suffix
        output.append(line)
    return "\n".join(output) + "\n"


def short_schedule(steps: int) -> str:
    return f"kind\tstep\ttag\nS\t0\tscalar\nS\t{steps}\tscalar\n"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--route", choices=tuple(ROUTES), required=True)
    ap.add_argument("--steps", default="full",
                    help="full uses the frozen schedule; a positive integer creates a short run")
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--binary-dir", type=Path, default=ROOT / "build/openphase")
    ap.add_argument("--max-seconds", type=float, default=86400.0)
    args = ap.parse_args()
    if args.max_seconds <= 0:
        ap.error("--max-seconds must be positive")
    try:
        steps = None if args.steps == "full" else int(args.steps)
    except ValueError:
        ap.error("--steps must be 'full' or a positive integer")
    if steps is not None and steps < 1:
        ap.error("--steps must be positive")

    program, opi_rel, recipe_rel, schedule_rel, source_rel, is_halo = ROUTES[args.route]
    binary = args.binary_dir.resolve() / program
    if not binary.is_file():
        raise FileNotFoundError(f"built caller missing: {binary}; use scripts/build_openphase_callers.py")
    opi_src, schedule_src = ROOT / opi_rel, ROOT / schedule_rel
    recipe = json.loads((ROOT / recipe_rel).read_text())
    seed_path = ROOT / "inputs/h24_cold_initial_field.npz"
    with np.load(seed_path, allow_pickle=False) as archive:
        phi = np.ascontiguousarray(archive["phi"], dtype=np.float64)
    input_hash = sha_bytes(phi.tobytes())
    if input_hash != recipe["seed_sha256_raw"]:
        raise ValueError("H24 initial field raw hash differs from frozen recipe")

    outdir = args.output_dir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    if any(outdir.iterdir()):
        raise FileExistsError(f"output directory must be empty: {outdir}")
    raw_path = outdir / "initial_field.raw"
    raw_path.write_bytes(phi.tobytes())
    opi_text = opi_src.read_text()
    if steps is not None:
        opi_text = replace_steps(opi_text, steps)
        schedule_text = short_schedule(steps)
        schedule_hash = sha_bytes(schedule_text.encode())
    else:
        schedule_text = schedule_src.read_text()
        schedule_hash = sha_file(schedule_src)
    opi_path = outdir / "ProjectInput.opi"
    schedule_path = outdir / "schedule.tsv"
    opi_path.write_text(opi_text)
    schedule_path.write_text(schedule_text)

    summary = outdir / "summary.tsv"
    fields = outdir / "fields.oph2df1"
    tags = outdir / "field_tags.tsv"
    halo = outdir / "halo_metrics.tsv"
    argv = [str(binary), "native", str(opi_path), str(raw_path), str(summary),
            str(fields), str(tags), str(schedule_path), str(recipe["factor_a"]),
            str(recipe["area_t99_short_W2"]), str(recipe["area_t99_long_W2"]),
            str(args.max_seconds), "0"]
    if is_halo:
        argv += [str(recipe["chi"]), str(recipe["radius_cells"]), "0", str(halo)]
    completed = subprocess.run(argv, check=False, text=True, capture_output=True)
    (outdir / "stdout.log").write_text(completed.stdout)
    (outdir / "stderr.log").write_text(completed.stderr)
    if completed.returncode != 0:
        raise subprocess.CalledProcessError(completed.returncode, argv,
                                            output=completed.stdout, stderr=completed.stderr)

    scalar_rows = []
    with summary.open(newline="") as stream:
        scalar_rows = list(csv.DictReader(stream, delimiter="\t"))
    final_step = int(float(scalar_rows[-1]["step"])) if scalar_rows else -1
    expected_step = steps if steps is not None else int(max(
        int(row.split("\t")[1]) for row in schedule_text.splitlines()[1:] if row.startswith("S\t")))
    partial = (summary.with_suffix(summary.suffix + ".status").exists() or final_step < expected_step)
    record = {
        "schema": "openphase-route-run-v1",
        "status": "PARTIAL" if partial else "COMPLETED",
        "route": args.route,
        "caller_source": source_rel,
        "caller_source_sha256": sha_file(ROOT / source_rel),
        "compiled_binary_sha256": sha_file(binary),
        "upstream_revision": recipe["upstream_revision"],
        "recipe": recipe_rel,
        "recipe_sha256": sha_file(ROOT / recipe_rel),
        "opi_sha256": sha_bytes(opi_text.encode()),
        "schedule_sha256": schedule_hash,
        "input_raw_sha256": input_hash,
        "raw_timestep_W2": recipe["raw_dt_W2"],
        "expected_terminal_step": expected_step,
        "observed_terminal_step": final_step,
        "scalar_rows": len(scalar_rows),
        "outputs": {p.name: sha_file(p) for p in (summary, fields, tags) if p.exists()},
    }
    if is_halo and halo.exists():
        record["outputs"][halo.name] = sha_file(halo)
    (outdir / "run_manifest.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps({"status": record["status"], "route": args.route,
                      "terminal_step": final_step, "manifest": "run_manifest.json"}, indent=2))
    if partial:
        raise SystemExit(3)


if __name__ == "__main__":
    main()
