#!/usr/bin/env python3
"""Verify release file bytes using the distributed SHA-256 manifest."""
from __future__ import annotations
import hashlib
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]

def main() -> int:
    manifest = ROOT / "provenance" / "SHA256SUMS"
    if not manifest.is_file():
        print("Missing provenance/SHA256SUMS", file=sys.stderr)
        return 2
    failed = []
    checked = 0
    for line in manifest.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        digest, name = line.split("  ", 1)
        path = ROOT / name
        if not path.resolve().is_relative_to(ROOT.resolve()):
            raise ValueError(f"Path outside repository: {name}")
        if not path.is_file():
            failed.append(f"missing: {name}")
            continue
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        checked += 1
        if actual != digest:
            failed.append(f"hash mismatch: {name}")
    if failed:
        print("\n".join(failed), file=sys.stderr)
        return 1
    print(f"PASS: {checked} release files match SHA-256 manifest")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
