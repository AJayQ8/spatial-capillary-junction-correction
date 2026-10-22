#!/usr/bin/env python3
"""Public one-step kernel smoke; optionally compare to original private sources.

Without --reference-root this is self-contained and checks the distributed
sources on the distributed initial fields. With that option it additionally
compiles the retained authoring-tree sources and requires exact source and
output identity. Neither mode runs a production trajectory.
"""
from __future__ import annotations

import argparse
import ctypes
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src/python"))
from b2_api import step as b2_step

DP = ctypes.POINTER(ctypes.c_double)


def sha_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def raw_sha(a: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()


def load_field(path: Path) -> np.ndarray:
    if path.suffix == ".npy":
        field = np.load(path, allow_pickle=False)
    else:
        with np.load(path, allow_pickle=False) as archive:
            field = archive["phi"]
    field = np.ascontiguousarray(field, dtype=np.float64)
    if field.ndim != 3 or field.shape[0] != 3:
        raise ValueError(f"expected (3,ny,nx) input field: {path}")
    return field


def compile_source(source: Path, output: Path) -> None:
    suffix = "-dynamiclib" if sys.platform == "darwin" else "-shared"
    cmd = [*os.environ.get("CC", "cc").split(), "-O3", "-std=gnu11", "-fPIC",
           suffix, str(source), "-o", str(output), "-lm"]
    subprocess.run(cmd, check=True, capture_output=True, text=True)


def local_mask(phi: np.ndarray, radius: int) -> tuple[np.ndarray, dict]:
    from local_memory_logic import next_suppression
    mask, summary = next_suppression(phi, None, radius, evaluate_overlap=True)
    return np.ascontiguousarray(mask, dtype=np.uint8), summary


def invariant_summary(phi: np.ndarray) -> dict:
    if not np.isfinite(phi).all():
        raise AssertionError("step output contains non-finite phase values")
    partition = float(np.max(np.abs(phi.sum(axis=0) - 1.0)))
    if partition > 3e-12:
        raise AssertionError(f"projected phase sum defect {partition:g}")
    return {"partition_linf": partition, "phi_min": float(phi.min()), "phi_max": float(phi.max())}


def grose_step(library: Path, phi: np.ndarray, recipe: dict) -> np.ndarray:
    ny, nx = phi.shape[1:]
    sigma = np.ascontiguousarray(recipe["sigma"], dtype=np.float64).ravel()
    mobility = np.ascontiguousarray(recipe["mobility"], dtype=np.float64).ravel()
    mu0 = np.ascontiguousarray(recipe["mu0"], dtype=np.float64).ravel()
    lib = ctypes.CDLL(str(library))
    fn = lib.grose_step_v001
    fn.argtypes = [DP, DP, ctypes.c_int, ctypes.c_int, ctypes.c_double,
                   ctypes.c_double, ctypes.c_double, ctypes.c_double,
                   DP, DP, DP, DP]
    fn.restype = ctypes.c_int
    contiguous = lambda a: np.ascontiguousarray(a).ctypes.data_as(DP)
    out = np.empty_like(phi)
    stats = np.zeros(8, dtype=np.float64)
    status = fn(contiguous(phi), contiguous(out), ny, nx,
                recipe["h_W"], recipe["eta_W"], recipe["dt_W2"],
                recipe["amplitude"], contiguous(sigma), contiguous(mobility),
                contiguous(mu0), contiguous(stats))
    if status:
        raise RuntimeError(f"Grose kernel returned status {status}")
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reference-root", type=Path,
                    help="Optional original AJ Physics source root; not needed by public smoke mode")
    ap.add_argument("--output", type=Path,
                    default=ROOT / "build/qualification/kernel_qualification.json")
    args = ap.parse_args()
    reference_root = args.reference_root.resolve() if args.reference_root else None

    h24_recipe = json.loads((ROOT / "configs/b2_h24.json").read_text())
    n20_recipe = json.loads((ROOT / "configs/b2_n20_warm.json").read_text())
    grose_recipe = json.loads((ROOT / "configs/grose_h24.json").read_text())
    h24 = load_field(ROOT / "inputs/h24_cold_initial_field.npz")
    n20 = load_field(ROOT / "inputs/n20_warm_time1_phi.npy")
    if raw_sha(h24) != h24_recipe["seed_sha256_raw"]:
        raise AssertionError("H24 frozen field hash does not match its recipe")
    if raw_sha(n20) != n20_recipe["initial_field"]["raw_sha256"]:
        raise AssertionError("N20 frozen field hash does not match its recipe")

    source_pairs = {
        "b2_d5_native": (ROOT / "src/kernels/b2_d5_native.c", "bench/b2_fast.c"),
        "b2_d5_halo": (ROOT / "src/kernels/b2_d5_halo.c", "bench/junction_capillary_halo_probe_v001/capillary_halo_fast.c"),
        "b2_d5_local_memory": (ROOT / "src/kernels/b2_d5_local_memory.c", "bench/grain_overlap_local_memory_v001/capillary_halo_local.c"),
        "b2_d9_native": (ROOT / "src/kernels/b2_d9_native.c", "bench/calibration_control_results/isotropic_stencil_diagnosis_v001/b2_nine_point_full_v001.c"),
        "b2_d9_halo": (ROOT / "src/kernels/b2_d9_halo.c", "bench/grain_applicability_diagnosis_v001/nine_point_eta008_v001/capillary_halo_nine_point.c"),
        "b2_d9_local_memory": (ROOT / "src/kernels/b2_d9_local_memory.c", "bench/grain_applicability_diagnosis_v001/nine_point_eta008_v001/cold_fixed_recipe_v001/capillary_halo_local_nine_point.c"),
        "grose_full_kernel": (ROOT / "src/kernels/grose_full_kernel.c", "bench/grose_complete_comparison_v001/grose_full_kernel_v002.c"),
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    checks = []
    suffix = ".dylib" if sys.platform == "darwin" else ".so"
    with tempfile.TemporaryDirectory(prefix="scjc-qualify-") as temp_dir:
        tmp = Path(temp_dir)
        public_libs: dict[str, Path] = {}
        reference_libs: dict[str, Path] = {}
        source_records = []
        for name, (public_source, relative_reference) in source_pairs.items():
            public_source_hash = sha_file(public_source)
            public_lib = tmp / f"{name}_public{suffix}"
            compile_source(public_source, public_lib)
            public_libs[name] = public_lib
            record = {"kernel": name, "release_source_sha256": public_source_hash,
                      "compiled_library_sha256": sha_file(public_lib)}
            if reference_root is not None:
                original = reference_root / relative_reference
                if not original.is_file():
                    raise FileNotFoundError(f"reference source missing: {relative_reference}")
                original_hash = sha_file(original)
                if public_source_hash != original_hash:
                    raise AssertionError(f"{name} release source differs from retained source")
                ref_lib = tmp / f"{name}_reference{suffix}"
                compile_source(original, ref_lib)
                reference_libs[name] = ref_lib
                record.update({"retained_source_sha256": original_hash,
                               "source_byte_match": True,
                               "retained_compiled_library_sha256": sha_file(ref_lib)})
            source_records.append(record)

        # Native and whole-halo D5/D9 kernels use the same H24 field and recipe.
        for stencil in ("d5", "d9"):
            outcomes = {}
            for route in ("native", "halo"):
                name = f"b2_{stencil}_{route}"
                output, witness = b2_step(public_libs[name], h24, h24_recipe,
                                          stencil=stencil, route=route)
                invariant = invariant_summary(output)
                row = {"check": f"{name}_public_smoke", "input_raw_sha256": raw_sha(h24),
                       "output_raw_sha256": raw_sha(output), **invariant, "kernel": witness}
                if reference_root is not None:
                    ref_output, _ = b2_step(reference_libs[name], h24, h24_recipe,
                                            stencil=stencil, route=route)
                    defect = float(np.max(np.abs(output - ref_output)))
                    if not np.array_equal(output, ref_output):
                        raise AssertionError(f"{name} retained-source output differs by {defect:g}")
                    row.update({"retained_output_raw_sha256": raw_sha(ref_output),
                                "max_abs_difference": defect, "exact_output_match": True})
                checks.append(row)
                outcomes[route] = output

            zero = dict(h24_recipe)
            zero["lambda_halo"] = 0.0
            native_zero, _ = b2_step(public_libs[f"b2_{stencil}_native"], h24, zero,
                                     stencil=stencil, route="native")
            halo_zero, _ = b2_step(public_libs[f"b2_{stencil}_halo"], h24, zero,
                                   stencil=stencil, route="halo")
            defect = float(np.max(np.abs(native_zero - halo_zero)))
            if defect > 2e-13:
                raise AssertionError(f"{stencil} chi=0 halo/native identity defect {defect:g}")
            checks.append({"check": f"{stencil}_chi0_halo_equals_native",
                           "max_abs_difference": defect})

        # Local-memory state update uses the exact N20 warm recipe, not H24 settings.
        mask, mask_meta = local_mask(n20, int(n20_recipe["radius_cells"]))
        for stencil in ("d5", "d9"):
            name = f"b2_{stencil}_local_memory"
            output, witness = b2_step(public_libs[name], n20, n20_recipe,
                                      stencil=stencil, route="local", suppression=mask)
            invariant = invariant_summary(output)
            row = {"check": f"{name}_n20_local_state_smoke",
                   "input_raw_sha256": raw_sha(n20), "output_raw_sha256": raw_sha(output),
                   "suppressed_nodes": int(mask.sum()), "suppression": {
                       key: value for key, value in mask_meta.items() if key != "overlap_centers"},
                   **invariant, "kernel": witness}
            if reference_root is not None:
                ref_output, _ = b2_step(reference_libs[name], n20, n20_recipe,
                                        stencil=stencil, route="local", suppression=mask)
                defect = float(np.max(np.abs(output - ref_output)))
                if not np.array_equal(output, ref_output):
                    raise AssertionError(f"{name} retained-source output differs by {defect:g}")
                row.update({"retained_output_raw_sha256": raw_sha(ref_output),
                            "max_abs_difference": defect, "exact_output_match": True})
            checks.append(row)

        grose_name = "grose_full_kernel"
        grose_out = grose_step(public_libs[grose_name], h24, grose_recipe)
        grose_inv = invariant_summary(grose_out)
        grose_row = {"check": "grose_nonzero_step_smoke", "input_raw_sha256": raw_sha(h24),
                     "output_raw_sha256": raw_sha(grose_out), **grose_inv}
        if reference_root is not None:
            reference_out = grose_step(reference_libs[grose_name], h24, grose_recipe)
            defect = float(np.max(np.abs(grose_out - reference_out)))
            if not np.array_equal(grose_out, reference_out):
                raise AssertionError(f"Grose retained-source output differs by {defect:g}")
            grose_row.update({"retained_output_raw_sha256": raw_sha(reference_out),
                              "max_abs_difference": defect, "exact_output_match": True})
        checks.append(grose_row)

    result = {
        "schema": "portable-kernel-qualification-v2",
        "status": "PASS",
        "mode": "public plus retained-source comparison" if reference_root is not None else "public self-contained smoke",
        "note": "One-step operator checks only; no full production trajectory was rerun.",
        "inputs": {
            "h24": {"shape": list(h24.shape), "raw_sha256": raw_sha(h24)},
            "n20_warm": {"shape": list(n20.shape), "raw_sha256": raw_sha(n20),
                         "source_time_W2": n20_recipe["initial_field"]["source_time_W2"],
                         "lab_y_origin_W": n20_recipe["initial_field"]["lab_y_origin_W"]},
        },
        "sources": source_records,
        "checks": checks,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": "PASS", "mode": result["mode"],
                      "check_count": len(checks), "output": str(args.output)}, indent=2))


if __name__ == "__main__":
    main()
