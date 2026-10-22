# Zenodo archive handoff

The GitHub repository and versioned release are the source snapshot for the
archive. Zenodo publication is pending; this file does not assert that a DOI
has been reserved or a record published.

For the initial record, use the complete source archive attached to release
`v1.0.0`, including code, configurations, compact inputs, data and manuscript
sources. Preserve its SHA-256 checksum. The GitHub automatically generated
source ZIP is also a complete snapshot of the tagged repository.

Suggested metadata, to be checked in the Zenodo form:

- Title: **Spatial capillary correction for driven multiphase-field junctions:
  code, benchmarks and figure data**.
- Creator: **Ayas Alzanki**, ORCID **0009-0004-3121-3778**.
- Affiliation: **Energy and Building Research Center, Kuwait Institute for
  Scientific Research (KISR)**.
- Version: **1.0.0**; resource type: **Software**.
- Description: use the scope summary in `README.md`, including the declared
  Grose–Asimow reconstruction and separately obtained OpenPhase dependency.
- License: original material **GPL-3.0-only**, with the retained third-party
  LPPL notices documented in `THIRD_PARTY.md`.
- Related source: <https://github.com/AJayQ8/spatial-capillary-junction-correction/releases/tag/v1.0.0>.

Once the record is available, insert its actual DOI in `CITATION.cff`, this
repository's README and the manuscript data-availability statement. Preserve
the original release rather than silently replacing its files. A DOI for the
archive must not be described as a DOI for a peer-reviewed journal article.
