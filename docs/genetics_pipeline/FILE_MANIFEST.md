# File Manifest

Provenance for the genetics and MR code in this repository. Every file was
taken from the `regenie_pipeline` working directory on the HMS O2 cluster and
then curated for release: the published set covers the **full cohort** only,
for **BMI**, **HbA1c**, and **TRIG_HDL_RATIO**. Scripts were edited to drop
other traits and the glycemic strata, so some differ from their working copies.

Paths below are relative to `/n/groups/patel/sivateja/regenie_pipeline/`.

## Split-sample REGENIE GWAS

- `gwas_regenie/scripts/combine_regenie_results.sh` ← `scripts/combine_regenie_results.sh`
- `gwas_regenie/scripts/convert_to_fuma_format.py` ← `scripts/convert_to_fuma_format.py`
- `gwas_regenie/scripts/filter_pheno_to_genotype.sh` ← `scripts/filter_pheno_to_genotype.sh`
- `gwas_regenie/scripts/fix_duplicates.sh` ← `scripts/fix_duplicates.sh`
- `gwas_regenie/scripts/fix_eids_rds.R` ← `scripts/fix_eids_rds.R`
- `gwas_regenie/scripts/generate_fuma_p0.1_files.sh` ← `scripts/generate_fuma_p0.1_files.sh`
- `gwas_regenie/scripts/generate_inputs.R` ← `scripts/generate_inputs.R`
- `gwas_regenie/scripts/generate_inputs_memory_efficient.R` ← `scripts/generate_inputs_memory_efficient.R`
- `gwas_regenie/scripts/generate_regenie_inputs_prespec7_protein_bmi.R` ← `scripts/generate_regenie_inputs_prespec7_protein_bmi.R`
- `gwas_regenie/scripts/pheno_list_full_per_chr.txt` ← `scripts/pheno_list_full_per_chr.txt`
- `gwas_regenie/scripts/pheno_list_hallmarks.txt` ← `scripts/pheno_list_hallmarks.txt`
- `gwas_regenie/scripts/pheno_list_prespec7_protein_bmi_full.txt` ← `scripts/pheno_list_prespec7_protein_bmi_full.txt`
- `gwas_regenie/scripts/pheno_list_proteins.txt` ← `scripts/pheno_list_proteins.txt`
- `gwas_regenie/scripts/prefilter_snps.sh` ← `scripts/prefilter_snps.sh`
- `gwas_regenie/scripts/preflight_checks.R` ← `scripts/preflight_checks.R`
- `gwas_regenie/scripts/process_hallmarks_fuma.sh` ← `scripts/process_hallmarks_fuma.sh`
- `gwas_regenie/scripts/read_regenie_summary_stats.R` ← `scripts/read_regenie_summary_stats.R`
- `gwas_regenie/scripts/read_regenie_summary_stats.py` ← `scripts/read_regenie_summary_stats.py`
- `gwas_regenie/scripts/regenie_to_ldsc_munge_input.py` ← `scripts/regenie_to_ldsc_munge_input.py`
- `gwas_regenie/scripts/resume_regenie.sh` ← `scripts/resume_regenie.sh`
- `gwas_regenie/scripts/run_generate_inputs.sh` ← `scripts/run_generate_inputs.sh`
- `gwas_regenie/scripts/run_regenie_array.sh` ← `scripts/run_regenie_array.sh`
- `gwas_regenie/scripts/run_regenie_array_per_chr.sh` ← `scripts/run_regenie_array_per_chr.sh`
- `gwas_regenie/scripts/run_regenie_test.sh` ← `scripts/run_regenie_test.sh`
- `gwas_regenie/scripts/run_test_regenie.sh` ← `scripts/run_test_regenie.sh`
- `gwas_regenie/scripts/step1_concatenate_hallmarks.sh` ← `scripts/step1_concatenate_hallmarks.sh`
- `gwas_regenie/scripts/step2_convert_to_fuma.sh` ← `scripts/step2_convert_to_fuma.sh`
- `gwas_regenie/scripts/step3_filter_p0.1.sh` ← `scripts/step3_filter_p0.1.sh`
- `gwas_regenie/scripts/submit_full_per_chr.sh` ← `scripts/submit_full_per_chr.sh`
- `gwas_regenie/scripts/submit_hallmarks_full.sh` ← `scripts/submit_hallmarks_all_strata.sh` (renamed)
- `gwas_regenie/scripts/submit_proteins_full.sh` ← `scripts/submit_proteins_all_strata.sh` (renamed)
- `gwas_regenie/scripts/submit_regenie.sh` ← `scripts/submit_regenie.sh`
- `gwas_regenie/scripts/submit_regenie_prespec7_proteins_full.sh` ← `scripts/submit_regenie_prespec7_proteins_full.sh`
- `gwas_regenie/scripts/submit_test_bmi.sh` ← `scripts/submit_test_bmi.sh`
- `gwas_regenie/scripts/test_bmi_minimal.sh` ← `scripts/test_bmi_minimal.sh`
- `gwas_regenie/scripts/test_bmi_quick.sh` ← `scripts/test_bmi_quick.sh`
- `gwas_regenie/scripts/test_bmi_quick_subset.sh` ← `scripts/test_bmi_quick_subset.sh`
- `gwas_regenie/scripts/test_bmi_tiny.sh` ← `scripts/test_bmi_tiny.sh`
- `gwas_regenie/scripts/test_bmi_with_filter.sh` ← `scripts/test_bmi_with_filter.sh`
- `gwas_regenie/scripts/test_full_stratum.sh` ← `scripts/test_full_stratum.sh`
- `gwas_regenie/scripts/test_regenie.sh` ← `scripts/test_regenie.sh`
- `gwas_regenie/scripts/test_single_protein.sh` ← `scripts/test_single_protein.sh`

## GWAS pipeline documentation

- `gwas_regenie/docs/README.md` ← `README.md`
- `gwas_regenie/docs/RESUME_GUIDE.md` ← `RESUME_GUIDE.md`
- `gwas_regenie/docs/RUN_PLAN.md` ← `RUN_PLAN.md`
- `gwas_regenie/docs/SUMMARY_STATS_PATHS.md` ← `SUMMARY_STATS_PATHS.md`
- `gwas_regenie/docs/TEST_RUN_GUIDE.md` ← `TEST_RUN_GUIDE.md`

## Two-sample and bidirectional MR

- `mr_twosample/merge_cis_pqtl_coloc.R` ← `merge_cis_pqtl_coloc.R`

- `mr_twosample/scripts/bidirectional_mr_bmi_to_proteins.R` ← `scripts/bidirectional_mr_bmi_to_proteins.R`
- `mr_twosample/scripts/bidirectional_mr_decode_trig_hdl.R` ← `scripts/bidirectional_mr_decode_trig_hdl.R`
- `mr_twosample/scripts/bidirectional_mr_lep_bmi.R` ← `scripts/bidirectional_mr_lep_bmi.R`
- `mr_twosample/scripts/combine_cis_trans_results.R` ← `scripts/combine_cis_trans_results.R`
- `mr_twosample/scripts/combine_reverse_mr_r2clump.R` ← `scripts/combine_reverse_mr_r2clump.R`
- `mr_twosample/scripts/generate_bidirectional_mr_tables.R` ← `scripts/generate_bidirectional_mr_tables.R`
- `mr_twosample/scripts/lepr_bmi_adjusted_for_lep.R` ← `scripts/lepr_bmi_adjusted_for_lep.R`
- `mr_twosample/scripts/merge_trig_hdl_bidirectional.R` ← `scripts/merge_trig_hdl_bidirectional.R`
- `mr_twosample/scripts/mr_decode_lep_bmi.R` ← `scripts/mr_decode_lep_bmi.R`
- `mr_twosample/scripts/mr_decode_lepr_bmi.R` ← `scripts/mr_decode_lepr_bmi.R`
- `mr_twosample/scripts/mr_summary_table.R` ← `scripts/mr_summary_table.R`
- `mr_twosample/scripts/reverse_mr_bmi_to_all_proteins.R` ← `scripts/reverse_mr_bmi_to_all_proteins.R`
- `mr_twosample/scripts/reverse_mr_phenotype_to_all_proteins.R` ← `scripts/reverse_mr_phenotype_to_all_proteins.R`
- `mr_twosample/scripts/reverse_mr_r2clump_sensitivity.R` ← `scripts/reverse_mr_r2clump_sensitivity.R`
- `mr_twosample/scripts/reverse_observational_pewas.R` ← `scripts/reverse_observational_pewas.R`
- `mr_twosample/scripts/submit_all_proteins_cis_trans.sh` ← `scripts/submit_all_proteins_cis_trans.sh`
- `mr_twosample/scripts/submit_cis_trans_chunked.sh` ← `scripts/submit_cis_trans_chunked.sh`
- `mr_twosample/scripts/submit_decode_lep_bmi.sh` ← `scripts/submit_decode_lep_bmi.sh`
- `mr_twosample/scripts/submit_leptin_cis_trans.sh` ← `scripts/submit_leptin_cis_trans.sh`
- `mr_twosample/scripts/submit_reverse_mr_r2clump.sh` ← `scripts/submit_reverse_mr_r2clump.sh`
- `mr_twosample/scripts/submit_reverse_pewas.sh` ← `scripts/submit_reverse_pewas.sh`
- `mr_twosample/scripts/twosampleMR_all_proteins_cis_trans.R` ← `scripts/twosampleMR_all_proteins_cis_trans.R`
- `mr_twosample/scripts/twosampleMR_cis_trans_chunked.R` ← `scripts/twosampleMR_cis_trans_chunked.R`
- `mr_twosample/scripts/twosampleMR_leptin_cis_trans.R` ← `scripts/twosampleMR_leptin_cis_trans.R`
- `mr_twosample/scripts/twosampleMR_proteins_bmi.R` ← `scripts/twosampleMR_proteins_bmi.R`
- `mr_twosample/scripts/twosampleMR_proteins_hba1c.R` ← `scripts/twosampleMR_proteins_hba1c.R`
- `mr_twosample/scripts/twosampleMR_proteins_trig_hdl_ratio_full.R` ← `scripts/twosampleMR_proteins_trig_hdl_ratio_full.R`

## Figures and summary tables

- `figures_tables/run_variance_decomposition.sh` ← `run_variance_decomposition.sh`
- `figures_tables/validation_bmi_cis_mr_concordant_plot.R` ← `validation_bmi_cis_mr_concordant_plot.R`
- `figures_tables/validation_bmi_cis_mr_plot.R` ← `validation_bmi_cis_mr_plot.R`
- `figures_tables/validation_hba1c_cis_mr_plot.R` ← `validation_hba1c_cis_mr_plot.R`
- `figures_tables/variance_decomposition_bmi.R` ← `variance_decomposition_bmi.R`

- `figures_tables/scripts/create_master_integrated_table.R` ← `scripts/create_master_integrated_table.R`
- `figures_tables/scripts/fig2_enr_atlas.R` ← `scripts/fig2_enr_atlas.R`
- `figures_tables/scripts/fig_bidirectional_mr_summary.R` ← `scripts/fig_bidirectional_mr_summary.R`
- `figures_tables/scripts/fig_concordance_hba1c_forward.R` ← `scripts/fig_concordance_hba1c_forward.R`
- `figures_tables/scripts/fig_forward_obs_vs_step.R` ← `scripts/fig_forward_obs_vs_step.R`
- `figures_tables/scripts/fig_fwd_vs_rev_mr_scatterplots.R` ← `scripts/fig_fwd_vs_rev_mr_scatterplots.R`
- `figures_tables/scripts/fig_reverse_pewas_summary.R` ← `scripts/fig_reverse_pewas_summary.R`
- `figures_tables/scripts/funnel_table_obs_step_mr.R` ← `scripts/funnel_table_obs_step_mr.R`
- `figures_tables/scripts/funnel_table_obs_step_mr_global_bonferroni.R` ← `scripts/funnel_table_obs_step_mr_global_bonferroni.R`
- `figures_tables/scripts/reproduce_triangulation_numbers.R` ← `scripts/reproduce_triangulation_numbers.R`
- `figures_tables/scripts/reproduce_triangulation_v2.R` ← `scripts/reproduce_triangulation_v2.R`
- `figures_tables/scripts/supplemental_tables_manuscript_and_hba1c_forward.R` ← `scripts/supplemental_tables_manuscript_and_hba1c_forward.R`
- `figures_tables/scripts/validation_concordance_all_phenotypes.R` ← `scripts/validation_concordance_all_phenotypes.R`

## Genetics pipeline documentation

- `docs/genetics_pipeline/FILE_MANIFEST.md` ← `scripts/FILE_MANIFEST.md`
- `docs/genetics_pipeline/WALKTHROUGH_SPLIT_SAMPLE_MR.md` ← `WALKTHROUGH_SPLIT_SAMPLE_MR.md`
