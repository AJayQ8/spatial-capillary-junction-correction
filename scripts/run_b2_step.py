#!/usr/bin/env python3
"""Apply one released B2 Euler update; this is not a trajectory runner."""
from __future__ import annotations
import argparse, hashlib, json, pathlib, sys
import numpy as np
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src/python"))
from b2_api import step

def raw_hash(a): return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--recipe", default=str(ROOT/"configs/b2_h24.json"))
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--stencil", choices=("d5","d9"), default="d9")
    ap.add_argument("--route", choices=("native","halo","local"), default="halo")
    ap.add_argument("--library", help="Override build/kernels library path")
    args=ap.parse_args()
    recipe=json.loads(pathlib.Path(args.recipe).read_text())
    inp=pathlib.Path(args.input)
    if inp.suffix==".npy": phi=np.load(inp,allow_pickle=False)
    else:
        with np.load(inp,allow_pickle=False) as z: phi=z["phi"]
    phi=np.ascontiguousarray(phi,dtype=np.float64)
    schema=recipe.get("schema")
    if schema=="b2-h24-single-step-recipe-v1":
        expected_shape=(3,2400,100)
        expected_hash=recipe.get("seed_sha256_raw")
    elif schema=="b2-n20-warm-continuation-recipe-v1":
        expected_shape=(3,800,200)
        expected_hash=recipe.get("initial_field",{}).get("raw_sha256")
    else:
        expected_shape=None
        expected_hash=None
    if expected_shape is not None and phi.shape!=expected_shape:
        ap.error(f"recipe {schema} requires an input field with shape {expected_shape}")
    input_hash=raw_hash(phi)
    if expected_hash and input_hash!=expected_hash:
        ap.error("input field raw hash does not match the selected frozen recipe; choose the matching recipe")
    if args.route=="local":
        from local_memory_logic import next_suppression
        mask,_=next_suppression(phi,None,int(recipe["radius_cells"]),evaluate_overlap=True)
        mask=np.ascontiguousarray(mask,dtype=np.uint8)
    else: mask=None
    suffix=".dylib" if sys.platform=="darwin" else ".so"
    lib=pathlib.Path(args.library) if args.library else ROOT/"build/kernels"/f"b2_{args.stencil}_{args.route if args.route != 'local' else 'local_memory'}{suffix}"
    # Native and halo source names use b2_d9_native, b2_d9_halo, etc.
    if not lib.exists():
        expected=ROOT/"build/kernels"/f"b2_{args.stencil}_{'local_memory' if args.route=='local' else args.route}{suffix}"
        lib=expected
    out,summary=step(lib,phi,recipe,stencil=args.stencil,route=args.route,suppression=mask)
    output=pathlib.Path(args.output); output.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(output,phi=out)
    summary.update({"input_raw_sha256":input_hash,"output_raw_sha256":raw_hash(out),
                    "output":output.name,"input_shape":list(phi.shape),
                    "suppressed_nodes":int(mask.sum()) if mask is not None else None})
    output.with_suffix(output.suffix+".json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2))
if __name__=="__main__": main()
