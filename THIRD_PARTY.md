# License and attribution notes

Original research code and accompanying original materials in this repository
are distributed under the GNU General Public License, version 3 only, in
[LICENSE](LICENSE). Copyright 2026 Ayas Alzanki. Third-party files retain their
own notices and licenses as described below.

## OpenPhase

The OpenPhase transfer uses [ICAMS/OpenPhase](https://github.com/ICAMS/OpenPhase),
pinned to commit `f201751168f55eb5bd9539aa925fb8d61e5c3649`.
OpenPhase is distributed under GPL-3.0; its source and libraries are not vendored
here. The repository provides the research caller and build instructions. Users
must obtain and build that upstream dependency separately. Preserve upstream
copyright and license notices in any redistributed combined work.

## Elsevier manuscript styles

The files `paper/cas-sc.cls`, `paper/cas-common.sty`, `paper/README`,
`paper/manifest.txt`, and the CAS thumbnail images originate in Elsevier's CAS
bundle (class version 2.4, 2024-05-04). They retain their upstream LaTeX Project
Public License 1.3c-or-later notices. The two class/style files were checked
byte-for-byte against the official CTAN bundle when preparing this release;
they are unmodified. This directory contains only the subset needed by the
manuscript, not a complete CAS distribution. The complete unmodified bundle is
available from [CTAN](https://ctan.org/pkg/els-cas-templates). Formatting changes
for the article are made in the manuscript source.

`paper/elsarticle-num.bst` is Elsevier's numeric bibliography style, version 2.1;
its header specifies LPPL 1.3-or-later and retains Elsevier's copyright notice.
These third-party style files are excluded from the repository's default GPL
license. The relevant license is available from
[the LaTeX Project](https://www.latex-project.org/lppl/).

## Scientific precedents

The benchmark and B2 force decomposition follow Hoffrogge et al.; see the full
bibliography in `paper/references.bib`. The Grose-Asimow comparator is a
project-authored reconstruction from the published equations, not code supplied
or endorsed by those authors. See `docs/SIMULATIONS.md` and Supplementary Section
S7 for its declared conventions and limitations.

Python dependencies (NumPy, SciPy, Matplotlib, and development tools) are installed
separately and retain their respective licenses. No credentials, prebuilt
OpenPhase libraries, or private laboratory runtime are distributed.
