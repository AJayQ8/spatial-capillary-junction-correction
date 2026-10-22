#!/usr/bin/env python3
"""Apply one Grose-Asimow kernel update to a supplied field."""
from __future__ import annotations
import argparse, ctypes, json, pathlib, sys, hashlib
import numpy as np
ROOT=pathlib.Path(__file__).resolve().parents[1]
DP=ctypes.POINTER(ctypes.c_double)
def ptr(a): return np.ascontiguousarray(a,dtype=np.float64).ctypes.data_as(DP)
def sha(a): return hashlib.sha256(np.ascontiguousarray(a,dtype=np.float64).tobytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--recipe',default=str(ROOT/'configs/grose_h24.json'))
    ap.add_argument('--input',required=True); ap.add_argument('--output',required=True)
    ap.add_argument('--library',default=None); args=ap.parse_args()
    r=json.loads(pathlib.Path(args.recipe).read_text()); pth=pathlib.Path(args.input)
    if pth.suffix=='.npy': phi=np.load(pth,allow_pickle=False)
    else:
        with np.load(pth,allow_pickle=False) as z: phi=z['phi']
    phi=np.ascontiguousarray(phi,dtype=np.float64)
    if phi.ndim!=3 or phi.shape[0]!=3: raise ValueError('input must have shape (3,ny,nx)')
    suffix='.dylib' if sys.platform=='darwin' else '.so'
    libpath=pathlib.Path(args.library) if args.library else ROOT/'build/kernels'/f'grose_full_kernel{suffix}'
    lib=ctypes.CDLL(str(libpath.resolve())); fn=lib.grose_step_v001
    fn.argtypes=[DP,DP,ctypes.c_int,ctypes.c_int,ctypes.c_double,ctypes.c_double,
        ctypes.c_double,ctypes.c_double,DP,DP,DP,DP]; fn.restype=ctypes.c_int
    out=np.empty_like(phi); stats=np.zeros(8,dtype=np.float64)
    sigma=np.ascontiguousarray(r['sigma'],dtype=np.float64).ravel()
    mobility=np.ascontiguousarray(r['mobility'],dtype=np.float64).ravel()
    mu0=np.ascontiguousarray(r['mu0'],dtype=np.float64).ravel()
    code=fn(ptr(phi),ptr(out),phi.shape[1],phi.shape[2],r['h_W'],r['eta_W'],
        r['dt_W2'],r['amplitude'],ptr(sigma),ptr(mobility),ptr(mu0),ptr(stats))
    if code: raise RuntimeError(f'grose_step_v001 returned {code}')
    dest=pathlib.Path(args.output); dest.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(dest,phi=out)
    report={'schema':'single-kernel-step-result-v1','method':'Grose-Asimow source-based C kernel',
        'input_raw_sha256':sha(phi),'output_raw_sha256':sha(out),'shape':list(phi.shape),
        'dt_W2':r['dt_W2'],'stats':{'eliminated_phases':int(stats[0]),'crossing_events':int(stats[1]),
        'max_local_substeps':int(stats[2]),'max_roundoff_negative':float(stats[3]),
        'max_partition_cleanup':float(stats[4]),'partition_linf':float(stats[5])}}
    dest.with_suffix(dest.suffix+'.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
if __name__=='__main__': main()
