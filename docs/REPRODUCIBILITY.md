# Reproducibility

## Data and figures

The ten CSV files in `data/figure_data/` contain the measurements plotted in the paper. `scripts/build_figures.py` rebuilds the result figures; `scripts/build_method_schematic.py` creates the explanatory schematic. The figure scripts do not run simulations or refit measurements.

The five result and graphical-abstract PNGs reproduced byte-for-byte with Matplotlib 3.9.4. PDF serialization and fonts can differ between environments.

## Simulations

The [simulation guide](SIMULATIONS.md) lists the supported run commands, supplied inputs and external dependencies. The public B2 trajectory driver covers fixed-domain H24. N20 moving-window and complete Grose measurement drivers are not included. The OpenPhase callers can run the supplied full observation schedules.

The release checks cover kernel updates, input/configuration checks, short OpenPhase runs and one B2 restart comparison. They do not validate every controller transition or repeat a complete production trajectory. The test records in `provenance/` describe the coverage. Large historical checkpoint archives are not included.

## File checks

Run `python scripts/verify_release.py` to compare files with `provenance/SHA256SUMS`. `provenance/source_map.json` records source and input identities. Checksums establish file identity, not the correctness of the numerical method.

## Build the documents

The sources use Elsevier's CAS single-column class. Copy `paper/` to a build directory and run:

```sh
tectonic --keep-logs main.tex
tectonic --keep-logs supplement.tex
```

These commands create `main.pdf` and `supplement.pdf`. The checked-in main article is named `manuscript.pdf`. Tectonic 0.17.0 was used for the release builds.
