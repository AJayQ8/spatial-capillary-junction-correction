# Reproducibility levels

The release separates the following operations.

1. **Inspect retained results.** Open `paper/manuscript.pdf` and
   `paper/supplement.pdf`; inspect the ten CSV tables in `data/figure_data/`.
2. **Rebuild figures.** `scripts/build_figures.py` consumes those frozen CSVs;
   it does not solve a phase-field equation or refit measurements. The five
   result/graphical-abstract PNGs reproduced byte-for-byte in the retained
   Matplotlib 3.9.4 environment. `provenance/figure_rebuild_check.json` records
   this check. PDF serialization or fonts may vary across environments.
   `scripts/build_method_schematic.py` separately builds the illustrative
   schematic; it is not simulation data.
3. **Build and exercise numerical kernels.** The short qualification described
   in `docs/SIMULATIONS.md` builds the portable sources and checks updates
   against retained source implementations. These checks validate packaging and
   short-step equivalence; they do not rerun or independently validate a full
   production trajectory.
4. **Run the benchmark recipes.** The supplied source, inputs and settings can
   be used for new runs. Production horizons require substantial local compute.
   Follow `docs/SIMULATIONS.md` for commands, upstream dependencies and output
   conventions. Do not equate a truncated demonstration with a full paper run.

`provenance/source_map.json` maps exported numerical files to retained source
identifiers and hashes. `provenance/SHA256SUMS` authenticates release file bytes;
run `python scripts/verify_release.py` to check it. Hash agreement establishes
file identity, not scientific correctness.

## What the release does not claim

The comparisons cover the smooth driven-junction and two-pair consumption cases
in the paper. They do not establish a universal correction, a causal angle law,
a full continuum-limit study, or an optimized accuracy-cost frontier. The late
OpenPhase long-pair shape miss remains in the data. Empty-boundary observations
are recorded as such. The Grose-Asimow route is a declared reconstruction, not
an author-supplied implementation.

The full internal execution history is preserved with the research project.
This public package contains the compact inputs and reproducibility material,
not every checkpoint, machine-specific operational receipt or interrupted run.
The initial GitHub release has no invented DOI; a later Zenodo record can cite
this exact tag.

## Rebuild the manuscript

The retained sources use Elsevier's CAS single-column class. With Tectonic
0.17.0 (or a suitable TeX distribution) available:

```sh
cd paper
tectonic --keep-logs main.tex
tectonic --keep-logs supplement.tex
```

This creates `main.pdf` and `supplement.pdf`; the retained main article is named
`manuscript.pdf` in the repository. Copy the source directory to a build folder
first if you want to preserve the checked-in supplement PDF. A clean anonymous
GitHub clone rebuilt both documents with matching extracted page text and page
counts (11 and 6); `provenance/paper_rebuild_check.json` records that check.
The retained class emits an overfull-box warning at `maketitle`; no change to
manuscript layout was made for this software release.
