"""Small ctypes interface for the released, pinned B2 C kernels."""
from __future__ import annotations

import ctypes
from pathlib import Path
import numpy as np

DP = ctypes.POINTER(ctypes.c_double)
BP = ctypes.POINTER(ctypes.c_uint8)
IP = ctypes.POINTER(ctypes.c_int)


def _p(a, typ=DP):
    return np.ascontiguousarray(a).ctypes.data_as(typ)


def step(library: str | Path, phi: np.ndarray, recipe: dict, *, stencil: str,
         route: str, suppression: np.ndarray | None = None) -> tuple[np.ndarray, dict]:
    """Run one C Euler step. Returns a copied field and compact kernel witnesses."""
    p = np.ascontiguousarray(phi, dtype=np.float64)
    if p.ndim != 3 or p.shape[0] != 3 or min(p.shape[1:]) < 2:
        raise ValueError("phi must have shape (3, ny>=2, nx>=2)")
    ny, nx = map(int, p.shape[1:])
    n = ny * nx
    out = np.empty_like(p); rhs = np.empty_like(p); lap = np.empty_like(p)
    active = np.empty(p.shape, dtype=np.uint8); flux = np.zeros(3, dtype=np.float64)
    lib = ctypes.CDLL(str(Path(library).resolve()))
    gamma = np.ascontiguousarray(recipe["gamma"], dtype=np.float64).ravel()
    delta = np.ascontiguousarray(recipe["delta_f"], dtype=np.float64).ravel()
    mobility = np.ascontiguousarray(recipe["mobility"], dtype=np.float64).ravel()
    scalar = (ctypes.c_int, ctypes.c_int, ctypes.c_double, ctypes.c_double,
              ctypes.c_double, ctypes.c_double)
    if route == "native":
        fn = lib.b2_step
        fn.argtypes = [DP, DP, DP, DP, BP, DP, DP, DP, DP, *scalar]
        fn.restype = ctypes.c_int
        code = fn(_p(p), _p(out), _p(rhs), _p(lap), _p(active, BP), _p(gamma),
                  _p(delta), _p(mobility), _p(flux), ny, nx,
                  float(recipe["dx_W"]), float(recipe["eta_W"]),
                  float(recipe["dt_W2"]), float(recipe["epsilon_factor"]))
        summary = {"flux": flux.tolist()}
    elif route in ("halo", "local"):
        summary_arr = np.zeros(3, dtype=np.float64)
        w0 = np.empty((ny, nx), dtype=np.float64)
        horizontal = np.empty_like(w0); halo = np.empty_like(w0)
        radius = int(recipe["radius_cells"])
        deque = np.empty(max(ny, nx) + 2 * radius + 8, dtype=np.int32)
        fn = lib.b2_capillary_halo_step
        if route == "local":
            if suppression is None:
                raise ValueError("local route requires a 2D uint8 suppression mask")
            mask = np.ascontiguousarray(suppression, dtype=np.uint8)
            if mask.shape != (ny, nx):
                raise ValueError("suppression mask shape mismatch")
            fn.argtypes = [DP, BP, DP, DP, DP, BP, DP, DP, DP, DP, DP, DP,
                           DP, DP, IP, *scalar, ctypes.c_double, ctypes.c_int]
            args = [_p(p), _p(mask, BP)]
        else:
            fn.argtypes = [DP, DP, DP, DP, BP, DP, DP, DP, DP, DP, DP,
                           DP, DP, IP, *scalar, ctypes.c_double, ctypes.c_int]
            args = [_p(p)]
        fn.restype = ctypes.c_int
        args += [_p(out), _p(rhs), _p(lap), _p(active, BP), _p(gamma),
                 _p(delta), _p(mobility), _p(flux), _p(summary_arr), _p(w0),
                 _p(horizontal), _p(halo), _p(deque, IP), ny, nx,
                 float(recipe["dx_W"]), float(recipe["eta_W"]),
                 float(recipe["dt_W2"]), float(recipe["epsilon_factor"]),
                 float(recipe["lambda_halo"]), radius]
        code = fn(*args)
        summary = {"flux": flux.tolist(), "halo_max": float(summary_arr[0]),
                   "halo_positive_nodes": int(summary_arr[1]),
                   "halo_integral_W2": float(summary_arr[2])}
    else:
        raise ValueError("route must be native, halo or local")
    if code:
        raise RuntimeError(f"{stencil}/{route} kernel returned status {code}")
    summary.update({"stencil": stencil, "route": route, "grid": [ny, nx]})
    return out.copy(), summary
