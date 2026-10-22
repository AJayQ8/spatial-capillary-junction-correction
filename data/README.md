# Exported measurements

The ten CSV files in `figure_data/` are the frozen measurements used by the
accompanying manuscript. Lengths are normalized by strip width W and physical
time by W²/(M₀γ₀); `raw_dt` is the raw solver timestep. See the manuscript and
supplement for the force and clock conventions.

| CSV | Content |
| --- | --- |
| `n20_terminal_fixed92_profiles.csv` | Aligned terminal profiles at the same 92 physical transverse positions. |
| `n20_shape_history.csv` | Developed-front aligned RMS histories. |
| `n20_route_metrics.csv` | Developed-front speed and terminal shape summaries. |
| `h24_shape_samples.csv` | Short/long consumption-pair shape observations for B2 and reconstructed Grose routes. |
| `h24_speed_t99_metrics.csv` | Consumption speeds and disappearance-time upper errors. |
| `halo_parameter_robustness.csv` | Four finite strength/radius perturbations and the baseline. |
| `openphase_halo_shape_samples.csv` | OpenPhase native/corrected shape observations. |
| `openphase_route_metrics.csv` | OpenPhase speed, event, and shape summaries. |
| `openphase_half_timestep_summary.csv` | Original/half-step OpenPhase summaries with clock and retained archive identifiers. |
| `outer_profile_geometry.csv` | Saved-field outer-profile angle descriptors and signed speed errors. |

The shape lower/upper columns bracket segment-sampling readout, not continuum
model error. An empty boundary is a separate outcome, not zero shape error.
Most speed-error columns are nonnegative magnitudes; columns explicitly named
`signed_speed_error_percent` use (v − v*)/|v*|. Angles describe outer-profile
fits, not measured core/Herring angles or a causal selection law. A blank CSV
cell means unavailable/not applicable, not zero.

`receipt_id` and `archive_sha256` are retained provenance identifiers; they are
not URLs or requirements for running the public scripts. The complete internal
checkpoint archives are not included. Compact initial fields are supplied
separately in `inputs/`.
