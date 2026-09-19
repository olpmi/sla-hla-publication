# Source data notes

What each archived input under `data/source_data/` is, and which published output
it feeds. Historical source filenames are kept as they are, so a name may not
match the figure number it now supplies; the mappings below are authoritative.

## Classifier

Canonical classifier predictions and validation tables are the archived
five-seed results. The prediction tables themselves are in
`data/archived_predictions/`.

## Eplet positions and agreement

`figureS8_S13_UPDATED_*` supplies the locus-resolved eplet profiles shown as
**Figures S8 and S9**. `figureS_identity_eplet_agreement_corrected_p_plotted_values.csv`
supplies **Figure S10**.

`eplet_statistics_v6.csv` is the authoritative current statistical table (and
publication Table S3); its permutation p-values are computed from the archived
integer counts. `revised_eplet_statistics_corrected_p.csv` is the preserved
earlier rounded evidence — the v6 file avoids double rounding for displayed
p-values. `revised_eplet_statistics.csv` is a current alias for references in the
retained Table S6 annotations.

`Figure4_plotted_values.csv` records the mismatch-rate conversion and the
bootstrap interval sign inversion.

Statistical definitions for these tables are in
`EPLET_STATISTIC_DEFINITIONS.md`.

## Structural panels and annotations

`structural_panel_manifest_v6.csv` describes the current **Figures 5 and 6**,
with revised label counts; `structural_panel_manifest_original_measurements.csv`
preserves the original measurement manifest.

`structural_labels_v6.csv` and `structural_label_render_audit.csv` give the exact
displayed letters, coordinates and image anchors.
`structural_candidate_annotations_v6.csv` holds the full supplied annotation
candidate pool, with canonical pair-map agreement and label selections.
**Table S5 contains displayed positions only.**

CE RMSDs and label placement are reverified in
`structural_render_report_v6.json`.

Model paths recorded inside the original manifest identify source provenance and
are not paths in this repository. The packaged model files are in
`data/structures/models/` (24 models) and `data/structures/reference/` (3
experimental reference structures); match them by basename. The rendered PyMOL
panels are in `data/structural_panels/`.
