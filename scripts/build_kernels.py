#!/usr/bin/env python3
"""Compile the exact released B2 and Grose C kernels into build/kernels."""
from __future__ import annotations
import argparse, hashlib, json, os, pathlib, platform, shlex, subprocess

ROOT = pathlib.Path(__file__).resolve().parents[1]
SOURCES = sorted((ROOT / "src/kernels").glob("*.c"))

def sha(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()

def compile_one(src: pathlib.Path, outdir: pathlib.Path, cc: str) -> dict:
    suffix = ".dylib" if platform.system() == "Darwin" else ".so"
    out = outdir / (src.stem + suffix)
    cmd = shlex.split(cc) + ["-O3", "-std=gnu11", "-fPIC", "-dynamiclib" if suffix == ".dylib" else "-shared", str(src), "-o", str(out), "-lm"]
    subprocess.run(cmd, check=True)
    return {"source": src.relative_to(ROOT).as_posix(), "source_sha256": sha(src),
            "library": out.relative_to(ROOT).as_posix(), "library_sha256": sha(out),
            "command": cmd}

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cc", default=os.environ.get("CC", "cc"))
    ap.add_argument("--out", default=str(ROOT / "build/kernels"))
    args = ap.parse_args()
    outdir = pathlib.Path(args.out).resolve(); outdir.mkdir(parents=True, exist_ok=True)
    entries = [compile_one(src, outdir, args.cc) for src in SOURCES]
    manifest = {"schema": "portable-kernel-build-v1", "compiler": args.cc, "kernels": entries}
    (outdir / "build_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"built": len(entries), "manifest": str(outdir / "build_manifest.json")}, indent=2))

if __name__ == "__main__":
    main()
