#!/usr/bin/env python3
"""Build the pinned OpenPhase caller programs against an existing upstream build."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
UPSTREAM_REVISION = "f201751168f55eb5bd9539aa925fb8d61e5c3649"
CALLERS = {
    "openphase_h24_native_full_dt": "src/openphase_h24_native_v004.cpp",
    "openphase_h24_halo": "src/openphase_h24_halo_v001.cpp",
    "openphase_h24_native_halfdt": "src/openphase_h24_native_halfdt_v001.cpp",
}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--openphase-root", type=Path, required=True,
                    help="Existing ICAMS/OpenPhase checkout built at the pinned revision")
    ap.add_argument("--fftw-prefix", type=Path, default=os.environ.get("FFTW_PREFIX"),
                    help="Optional FFTW prefix containing lib/libfftw3")
    ap.add_argument("--output-dir", type=Path, default=ROOT / "build/openphase")
    args = ap.parse_args()
    upstream = args.openphase_root.resolve()
    revision = subprocess.check_output(
        ["git", "-C", str(upstream), "rev-parse", "HEAD"], text=True).strip()
    if revision != UPSTREAM_REVISION:
        raise SystemExit(f"OpenPhase revision mismatch: expected {UPSTREAM_REVISION}, got {revision}")
    required = [upstream / "include", upstream / "external", upstream / "lib/libOpenPhase.a",
                upstream / "external/miniz.o", upstream / "external/WinBase64/libbase64.o"]
    missing = [str(p) for p in required if not p.exists()]
    if missing:
        raise SystemExit("OpenPhase build is incomplete; missing: " + ", ".join(missing))
    cxx = os.environ.get("CXX", "clang++").split()
    if not cxx:
        raise SystemExit("CXX is empty")
    outdir = args.output_dir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    records = []
    for name, source_rel in CALLERS.items():
        source = ROOT / source_rel
        output = outdir / name
        command = [*cxx, "-std=c++17", "-O3", "-DNDEBUG", "-DACADEMIC",
                   "-Wno-unknown-pragmas", "-I", str(upstream / "include"),
                   "-I", str(upstream / "external"), str(source), "-o", str(output),
                   "-L", str(upstream / "lib"), "-lOpenPhase",
                   str(upstream / "external/miniz.o"),
                   str(upstream / "external/WinBase64/libbase64.o")]
        if args.fftw_prefix:
            fftw = args.fftw_prefix.resolve()
            command += ["-L", str(fftw / "lib"), "-lfftw3"]
        else:
            command += ["-lfftw3"]
        subprocess.run(command, check=True)
        records.append({"program": name, "source": source_rel,
                        "source_sha256": sha(source), "binary_sha256": sha(output)})
    manifest = {"schema": "openphase-caller-build-v1",
                "upstream": "ICAMS/OpenPhase", "upstream_revision": revision,
                "callers": records}
    (outdir / "build_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"status": "built", "count": len(records),
                      "output_dir": str(outdir), "upstream_revision": revision}, indent=2))


if __name__ == "__main__":
    main()
