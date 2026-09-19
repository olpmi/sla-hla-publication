# v6 source data

Canonical classifier predictions and validation tables are preserved from the v5 evidence set. The residue-replacement data are held separately under repository_handoff/ and are not publication supplemental data.

Historical source filenames retain their provenance. figureS8_S13_UPDATED_* supplies the locus-resolved eplet profiles now shown as Figures S8 and S9 (v5 S15/S16). figureS_identity_eplet_agreement_corrected_p_plotted_values.csv supplies Figure S10 (v5 S17). Figure_number_crosswalk.csv gives the publication mapping.

The authoritative current statistical table is eplet_statistics_v6.csv (and publication Table S3), whose permutation p-values are computed from the archived integer counts. revised_eplet_statistics_corrected_p.csv is the preserved earlier rounded evidence; the v6 file avoids double rounding for displayed p-values. revised_eplet_statistics.csv is a current alias for references in retained Table S6 annotations. Figure4_plotted_values.csv records the mismatch-rate conversion and bootstrap interval sign inversion.

structural_panel_manifest_v6.csv describes current Figures 5/6, with revised label counts. The original measurement manifest remains preserved in the source files. structural_labels_v6.csv and structural_label_render_audit.csv provide exact displayed letters/coordinates and image anchors. structural_candidate_annotations_v6.csv contains the full supplied annotation candidate pool, with canonical pair-map agreement and label selections. Table S5 contains displayed positions only. Raw models and the original annotation table are in structural_inputs/.

Historical model paths in the original manifest identify source provenance; use the basename under structural_inputs/models/ for packaged model files. CE RMSDs and placement are reverified in structural_render_report_v6.json. No inference, classifier retraining, or antibody-binding analysis was run for v6.
