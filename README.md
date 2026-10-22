# Spatial capillary correction for driven junctions

Code, benchmark inputs and figure data for **A spatial capillary correction for driven junction motion in double-obstacle multiphase-field simulations**, by Ayas Alzanki.

[Manuscript](paper/manuscript.pdf) · [Supplement](paper/supplement.pdf) · [Simulation guide](docs/SIMULATIONS.md)

The method extends a ternary-phase signal into a neighborhood of each junction and amplifies the capillary contribution without changing the bulk-driving term. The examples cover driven fronts, interacting junction pairs, a reconstruction of the Grose–Asimow method, and an implementation in OpenPhase using the same correction strength and relative radius.

![Spatial correction schematic](figures/method_schematic.png)

## Rebuild the figures

Python 3.10 or newer is required. From the repository root:

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python scripts/verify_release.py
python scripts/build_figures.py
python scripts/build_method_schematic.py --output-dir build/figures
```

The scripts read the included CSV files and write to `build/figures/`. The method schematic is an illustration. It does not represent simulation data.

## Run the calculations

Start with the [simulation guide](docs/SIMULATIONS.md). To compile and check the B2 and reconstructed Grose kernels:

```sh
python scripts/build_kernels.py
python tests/qualify_kernels.py
```

The B2 trajectory driver supports the fixed-domain H24 case. The N20 moving-window driver and complete Grose measurement pipeline are not included; their kernels, initial fields and reported measurements are supplied. OpenPhase run commands, inputs and upstream dependencies are described in the guide.

The release tests check kernel updates and selected restart behavior. They do not repeat the full benchmark calculations. See [Reproducibility](docs/REPRODUCIBILITY.md) for details.

## Contents

| Directory | Contents |
| --- | --- |
| `src/` | B2 and Grose kernels, Python helpers and OpenPhase callers. |
| `scripts/` | Build and run commands, figure scripts and checksum verification. |
| `configs/` | Benchmark parameters and observation settings. |
| `inputs/` | Initial fields and OpenPhase observation schedules. |
| `data/figure_data/` | Ten CSV datasets used in the figures and tables. |
| `figures/` | Figures and graphical abstract. |
| `paper/` | Manuscript, supplement and editable LaTeX sources. |
| `tests/` | Software checks. |
| `provenance/` | File hashes, source mapping and test records. |

## Scope

The paper demonstrates accuracy improvements in the tested driven-junction settings. It does not provide a universal parameter rule or an optimized accuracy–cost comparison. The data include a late OpenPhase shape error outside the selected tolerance and observations with no remaining boundary. The Grose–Asimow comparison uses our reconstruction of the published equations. The local-memory rule is an optional, history-dependent extension.

## Citation and license

Use [CITATION.cff](CITATION.cff) to cite version `1.0.1`, and cite the accompanying manuscript when using its results. A Zenodo DOI will be added when the archive is published.

Original material is licensed under GPL-3.0-only. The Elsevier LaTeX files retain their LPPL licenses; other dependencies retain their own licenses. See [LICENSE](LICENSE) and [THIRD_PARTY.md](THIRD_PARTY.md).
