# Numerical sources and simulation commands

The release contains the exact seven B2/Grose C kernels used by the retained
comparisons, compact H24 and N20 initial fields, four frozen OpenPhase input
files, two OpenPhase observation schedules, and the three authored OpenPhase
callers. The paper's numerical results and their exported data remain the
authoritative record. The commands below expose selected kernels and allow new
runs; the short release qualification is not a fresh full-horizon reproduction.

## B2 and Grose source kernels

Install Python 3.10+, NumPy and SciPy from `requirements.txt`, plus a C compiler.
From the repository root, build the C libraries:

```sh
python3 scripts/build_kernels.py
```

The distributed kernels cover D5/D9 native B2, whole-halo B2, local-memory B2,
and the project-authored Grose–Asimow reconstruction. The self-contained public
smoke test compiles those sources in a temporary directory and checks one step
on the packaged H24/N20 fields, including the zero-strength halo/native
identity:

```sh
python3 tests/qualify_kernels.py
```

This check needs no authoring worktree. It tests the distributed code and input
files; it does not compare against a separately fetched source archive. The
original N20 field is the time-one continuation state, not a cold start. Its
recipe records the source time, raw field hash and laboratory-coordinate origin.
The one-step CLI can compare native and halo updates on that exact input:

```sh
python3 scripts/run_b2_step.py \
  --recipe configs/b2_n20_warm.json \
  --input inputs/n20_warm_time1_phi.npy \
  --stencil d5 --route halo \
  --output build/n20_first_step.npz
```

The CLI checks each frozen recipe's grid and raw initial-field hash, so the H24
parameters cannot be applied silently to the N20 field (or vice versa).

The released trajectory loop is deliberately limited to the fixed H24 grid and
the B2 kernels. It writes observations and atomic checkpoints. For the local
memory route, each checkpoint retains the suppression mask and next controller
index as well as the phase field:

```sh
python3 scripts/run_b2_trajectory.py \
  --input inputs/h24_cold_initial_field.npz \
  --route local --stencil d9 --steps 1000 \
  --output-dir build/h24_local_example
python3 scripts/run_b2_trajectory.py \
  --resume build/h24_local_example/checkpoint_latest.npz \
  --route local --stencil d9 --steps 1000 \
  --output-dir build/h24_local_example
```

The H24 loop records simple field/invariant summaries; it does not reproduce the
paper's complete geometric profile, event, shape or signed-speed observer
pipeline. The N20 continuation additionally recenters its moving window by
integer rows, so the general H24 loop must not be applied to N20. The release
provides the N20 initial field and a one-step operator path, not a moving-window
N20 production driver. The Grose entrypoint is likewise a one-step kernel
adapter; it is a reconstruction from the stated convention, not the authors'
executable or a full production comparator.

## OpenPhase transfer callers

The OpenPhase source itself is not vendored. Obtain the public ICAMS/OpenPhase
repository at the pinned revision `f201751168f55eb5bd9539aa925fb8d61e5c3649`,
install its FFTW/compiler dependencies, and build the library following its
upstream instructions. The retained build used the serial static academic
configuration:

```sh
git clone https://github.com/ICAMS/OpenPhase.git OpenPhase
git -C OpenPhase checkout f201751168f55eb5bd9539aa925fb8d61e5c3649
cd OpenPhase
make -j1 CXX=clang++ SETTINGS="no-openmp static academic"
cd ..
```

On macOS, the pinned upstream metadata source uses Linux's `/proc/self/exe`.
The release includes a build-only compatibility preinclude for that call. Use
it when building the upstream library on macOS, leaving upstream source files
unchanged:

```sh
make -j1 \
  CXX="clang++ -include /absolute/path/to/release/src/openphase_macos_metadata_compat_v001.h" \
  SETTINGS="no-openmp static academic"
```

Then build the retained callers against that checkout (set `FFTW_PREFIX` if
FFTW is outside standard compiler paths):

```sh
python3 scripts/build_openphase_callers.py \
  --openphase-root ./OpenPhase --fftw-prefix /path/to/fftw
```

Caller identity matters. `openphase_h24_native_v004.cpp` is the native baseline;
it contains an optional regularized branch but does not implement the spatial
correction. `openphase_h24_halo_v001.cpp` applies the whole-halo correction and
is used at both raw timestep levels. The `native_halfdt` caller is a separate
half-step native source with a positive-finite timestep admission guard. The
frozen full- and half-step OPI files, schedules, H24 input and recipes select
the relevant combinations. The corrected transfer must use the halo caller.

Run a full saved schedule by selecting one route. This can require substantial
local compute; the release qualification uses only two steps per caller:

```sh
python3 scripts/run_openphase_route.py \
  --route halo-half --steps full \
  --output-dir build/openphase_halo_half
```

Use `--steps 2` for a two-update check. The half-step experiment is a fixed-grid
two-level timestep sensitivity check; it does not estimate an asymptotic
convergence order. No full-horizon OpenPhase run is repeated by the short
qualification command below.

## Bounded release qualification

The maintainer qualification builds the C kernels, exercises public one-step
updates, checks restart equivalence through the first H24 local-memory
controller evaluation, and compiles/runs each OpenPhase caller for two Euler
steps. The OpenPhase test needs an existing checkout at the exact upstream
revision and the FFTW prefix used by that build:

```sh
python3 tests/run_release_smoke.py \
  --openphase-root ./OpenPhase --fftw-prefix /path/to/fftw \
  --output-dir build/qualification/release_smoke
```

This qualification is short and bounded; it does not execute the paper's
full-horizon cases, establish a cost frontier, or validate numerical behavior
beyond the listed steps and stored observers. `provenance/source_map.json`
records release-relative source and input hashes. OpenPhase remains under its
upstream GPL-3.0 license; preserve its notices when redistributing a combined
build.
