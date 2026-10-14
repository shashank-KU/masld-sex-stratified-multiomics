# common/

`setup.R` is sourced by the first chunk of every `<Folder>/<Folder>.Rmd`: packages, paths (`data/`, `results/`), FDR thresholds, helpers, Source Data writer and the step-to-step workspace cache.

Run order (each script loads the previous step, or runs it first if needed):
Table1_cohort -> Fig1_hepatic_transcriptome -> Hepatic_lipidome -> Hepatic_polar_metabolome -> Fig2_serum_metabolome -> Fig3_external_validation -> Fig4_paired_liver_serum -> Fig5_pPRS -> Fig6_pathway_map -> Supp_EPoS_serum_validation -> Assemble_supplementary_tables

Settings (optional environment variables): MASLD_DATA_DIR, MASLD_RESULTS_DIR, MASLD_SUPP_XLSX.
