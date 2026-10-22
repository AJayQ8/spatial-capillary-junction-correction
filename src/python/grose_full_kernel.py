"""ctypes wrapper for the new full-field Grose implementation.

The wrapper does not import the scalar oracle. It owns reusable NumPy work
arrays, keeps the force update and profile correction as separate calls, and
exposes selected correction witnesses for independent equation review.
"""
from __future__ import annotations

import ctypes
from pathlib import Path

import numpy as np


_DP = ctypes.POINTER(ctypes.c_double)
_BP = ctypes.POINTER(ctypes.c_uint8)
_IP = ctypes.POINTER(ctypes.c_int)


def digest_array(a: np.ndarray) -> str:
    import hashlib
    return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()


class GroseFullKernel:
    """Full-field kernels with method-local reusable output buffers.

    Returned arrays are borrowed buffers and are overwritten by the next call
    to the same method. A method rejects an input that overlaps one of its
    output buffers, so it cannot accidentally turn a read into an in-place
    update. The intended force/profile sequence alternates ``out_step`` and
    ``out_profile``; those buffers do not overlap.
    """

    def __init__(self, library: str | Path, shape: tuple[int, int, int]):
        self.library_path = Path(library).resolve()
        self.lib = ctypes.CDLL(str(self.library_path))
        self.ny, self.nx = map(int, shape[1:])
        if shape[0] != 3 or min(self.ny, self.nx) < 2:
            raise ValueError("shape must be (3, ny>=2, nx>=2)")
        n = self.ny * self.nx
        self.out_step = np.empty(shape, dtype=np.float64)
        self.rhs = np.empty(shape, dtype=np.float64)
        self.lap = np.empty(shape, dtype=np.float64)
        self.eligible = np.empty(shape, dtype=np.uint8)
        self.out_profile = np.empty(shape, dtype=np.float64)
        self.psi = np.empty((2, self.ny, self.nx), dtype=np.float64)
        self.icf = np.empty((2, self.ny, self.nx), dtype=np.uint8)
        self.normal_x = np.empty((2, self.ny, self.nx), dtype=np.float64)
        self.normal_y = np.empty_like(self.normal_x)
        self.ext_x = np.empty_like(self.normal_x)
        self.ext_y = np.empty_like(self.normal_x)
        self.ext_weight = np.empty_like(self.normal_x)
        self.psi_corrected = np.empty_like(self.normal_x)
        self.anchor_x = np.empty_like(self.normal_x)
        self.anchor_y = np.empty_like(self.normal_x)
        self.fraction_x = np.empty_like(self.normal_x)
        self.fraction_y = np.empty_like(self.normal_x)
        self.anchor_valid = np.empty((2, self.ny, self.nx), dtype=np.uint8)
        self._n = n

        self._rhs_fn = self.lib.grose_rhs_field_v001
        self._rhs_fn.argtypes = [_DP, _DP, _DP, _BP, ctypes.c_int, ctypes.c_int,
            ctypes.c_double, ctypes.c_double, ctypes.c_double, _DP, _DP, _DP]
        self._rhs_fn.restype = ctypes.c_int

        self._step_fn = self.lib.grose_step_v001
        self._step_fn.argtypes = [_DP, _DP, ctypes.c_int, ctypes.c_int,
            ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
            _DP, _DP, _DP, _DP]
        self._step_fn.restype = ctypes.c_int

        self._profile_fn = self.lib.grose_profile_correct_v001
        self._profile_fn.argtypes = [_DP, _DP, ctypes.c_int, ctypes.c_int,
            ctypes.c_double, ctypes.c_double, _DP, _BP, _DP, _DP, _DP, _DP,
            _DP, _DP, _DP, _DP, _DP, _DP, _DP, _BP, _IP, ctypes.c_int, _DP, _DP]
        self._profile_fn.restype = ctypes.c_int

    @staticmethod
    def _ptr(a: np.ndarray):
        return np.ascontiguousarray(a, dtype=np.float64).ctypes.data_as(_DP)

    def _check_phi(self, phi: np.ndarray) -> np.ndarray:
        p = np.ascontiguousarray(phi, dtype=np.float64)
        if p.shape != (3, self.ny, self.nx) or not np.isfinite(p).all():
            raise ValueError("phi must be finite and match the kernel's (3,ny,nx) grid")
        return p

    @staticmethod
    def _reject_output_alias(p: np.ndarray, outputs: tuple[tuple[str, np.ndarray], ...],
                             method: str) -> None:
        for name, output in outputs:
            if np.shares_memory(p, output):
                raise ValueError(f"{method} input aliases its reusable {name} output buffer")

    def evaluate(self, phi: np.ndarray, *, h: float, eta: float,
                 amplitude: float, sigma: np.ndarray, mobility: np.ndarray,
                 mu0: np.ndarray) -> dict:
        p = self._check_phi(phi)
        self._reject_output_alias(p, (("rhs", self.rhs), ("laplacian", self.lap),
                                      ("eligible", self.eligible)), "evaluate")
        s, m, mu = (np.ascontiguousarray(v, dtype=np.float64) for v in (sigma, mobility, mu0))
        if any(a.shape != (3, 3) for a in (s, m, mu)):
            raise ValueError("pair matrices must have shape (3,3)")
        code = self._rhs_fn(self._ptr(p), self._ptr(self.rhs), self._ptr(self.lap),
            self.eligible.ctypes.data_as(_BP), self.ny, self.nx, float(h),
            float(eta), float(amplitude), self._ptr(s), self._ptr(m), self._ptr(mu))
        if code:
            raise RuntimeError(f"grose_rhs_field_v001 returned {code}")
        return {"rhs": self.rhs, "laplacian": self.lap, "eligible": self.eligible}

    def advance(self, phi: np.ndarray, *, h: float, eta: float, dt: float,
                amplitude: float, sigma: np.ndarray, mobility: np.ndarray,
                mu0: np.ndarray) -> tuple[np.ndarray, dict]:
        p = self._check_phi(phi)
        self._reject_output_alias(p, (("step", self.out_step),), "advance")
        s, m, mu = (np.ascontiguousarray(v, dtype=np.float64) for v in (sigma, mobility, mu0))
        stats = np.zeros(8, dtype=np.float64)
        code = self._step_fn(self._ptr(p), self._ptr(self.out_step), self.ny, self.nx,
            float(h), float(eta), float(dt), float(amplitude), self._ptr(s), self._ptr(m),
            self._ptr(mu), self._ptr(stats))
        if code:
            raise RuntimeError(f"grose_step_v001 returned {code}")
        summary = {"eliminated_phases": float(stats[0]), "crossing_events": float(stats[1]),
            "max_local_substeps": int(stats[2]), "max_roundoff_negative": float(stats[3]),
            "max_partition_cleanup": float(stats[4]), "partition_linf": float(stats[5])}
        return self.out_step, summary

    def correct_profiles(self, phi: np.ndarray, *, h: float, eta: float,
                         watch_yx: tuple[tuple[int, int], ...] = ()) -> tuple[np.ndarray, dict, list]:
        p = self._check_phi(phi)
        self._reject_output_alias(p, (("profile", self.out_profile),), "correct_profiles")
        coords = np.asarray(watch_yx, dtype=np.int32)
        if coords.size == 0:
            coords = np.empty((0, 2), dtype=np.int32)
        if coords.ndim != 2 or coords.shape[1] != 2:
            raise ValueError("watch_yx must be a sequence of (y,x) pairs")
        debug = np.full((len(coords), 2, 20), np.nan, dtype=np.float64)
        stats = np.zeros(16, dtype=np.float64)
        code = self._profile_fn(self._ptr(p), self._ptr(self.out_profile), self.ny, self.nx,
            float(h), float(eta), self._ptr(self.psi),
            self.icf.ctypes.data_as(_BP), self._ptr(self.normal_x), self._ptr(self.normal_y),
            self._ptr(self.ext_x), self._ptr(self.ext_y), self._ptr(self.ext_weight),
            self._ptr(self.psi_corrected), self._ptr(self.anchor_x), self._ptr(self.anchor_y),
            self._ptr(self.fraction_x), self._ptr(self.fraction_y),
            self.anchor_valid.ctypes.data_as(_BP),
            coords.ctypes.data_as(_IP), len(coords), self._ptr(debug), self._ptr(stats))
        if code:
            raise RuntimeError(f"grose_profile_correct_v001 returned {code}")
        rows = []
        keys = ("group", "icf", "extension_valid", "anchor_found", "target_psi",
            "anchor_x_grid", "anchor_y_grid", "anchor_psi", "signed_distance_cells",
            "theta_anchor", "theta_raw", "theta_clipped", "candidate_psi",
            "delta_psi", "beta1_fraction_at_anchor", "beta2_fraction_at_anchor",
            "corrected_psi", "anchor_interpolated", "fraction_sample_x_grid",
            "fraction_sample_y_grid")
        for idx, (y, x) in enumerate(coords.tolist()):
            for group in range(2):
                row = {"y": int(y), "x": int(x), "group": group}
                values = debug[idx, group]
                for k, name in enumerate(keys):
                    v = float(values[k])
                    row[name] = v if np.isfinite(v) else None
                rows.append(row)
        summary = {
            "edge_candidates_alpha": int(stats[0]), "edge_candidates_beta": int(stats[1]),
            "interior_alpha": int(stats[2]), "interior_beta": int(stats[3]),
            "corrected_alpha": int(stats[4]), "corrected_beta": int(stats[5]),
            "no_anchor_alpha": int(stats[6]), "no_anchor_beta": int(stats[7]),
            "max_abs_profile_delta_alpha": float(stats[8]),
            "max_abs_profile_delta_beta": float(stats[9]),
            "normalization_fallbacks": int(stats[10]), "new_beta_without_anchor": int(stats[11]),
            "nearest_bulk_interior_sources_alpha": int(stats[12]),
            "nearest_bulk_interior_sources_beta": int(stats[13]),
        }
        return self.out_profile, summary, rows
