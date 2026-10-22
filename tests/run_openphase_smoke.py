#!/usr/bin/env python3
"""Build and advance each released H24 OpenPhase caller by two steps."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--openphase-root", type=Path, required=True)
    ap.add_argument("--fftw-prefix", type=Path)
    ap.add_argument("--output-dir", type=Path, default=ROOT / "build/qualification/openphase_smoke")
    args = ap.parse_args()
    outdir = args.output_dir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    if any(outdir.iterdir()):
        raise FileExistsError(f"smoke output directory must be empty: {outdir}")
    binary_dir = outdir / "callers"
    build = [sys.executable, str(ROOT / "scripts/build_openphase_callers.py"),
             "--openphase-root", str(args.openphase_root.resolve()),
             "--output-dir", str(binary_dir)]
    if args.fftw_prefix:
        build += ["--fftw-prefix", str(args.fftw_prefix.resolve())]
    subprocess.run(build, check=True)
    records = []
    for route in ("native-full", "halo-full", "native-half", "halo-half"):
        route_dir = outdir / route
        command = [sys.executable, str(ROOT / "scripts/run_openphase_route.py"),
                   "--route", route, "--steps", "2", "--binary-dir", str(binary_dir),
                   "--output-dir", str(route_dir)]
        subprocess.run(command, check=True)
        manifest = json.loads((route_dir / "run_manifest.json").read_text())
        if manifest["status"] != "COMPLETED" or manifest["observed_terminal_step"] != 2:
            raise AssertionError(f"OpenPhase smoke route incomplete: {route}")
        records.append({"route": route, "manifest_sha256": sha(route_dir / "run_manifest.json"),
                        "input_raw_sha256": manifest["input_raw_sha256"],
                        "final_step": manifest["observed_terminal_step"],
                        "summary_sha256": manifest["outputs"]["summary.tsv"]})
    result = {"schema": "openphase-two-step-smoke-v1", "status": "PASS",
              "scope": "Two Euler updates per caller; not a trajectory reproduction.",
              "upstream_revision": "f201751168f55eb5bd9539aa925fb8d61e5c3649",
              "routes": records}
    (outdir / "qualification.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": "PASS", "route_count": len(records),
                      "qualification": str(outdir / "qualification.json")}, indent=2))


if __name__ == "__main__":
    main()
