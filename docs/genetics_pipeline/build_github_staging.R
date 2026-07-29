#!/usr/bin/env Rscript
# Build MAP_D_github_staging_OPTION_A/ for merging into stejat98/MAP-D (Option A monorepo)

root <- "/n/groups/patel/sivateja/regenie_pipeline"
out <- file.path(root, "MAP_D_github_staging_OPTION_A")

# Repository scope: BMI, HbA1c, and TRIG_HDL_RATIO only.
out_of_scope_patterns <- c(
  "alst", "_cad", "cad_", "ckd", "nafld", "_hf_", "scz",
  "motrpac", "drugbank", "phewas_", "safety", "plm", "lasso_quadrant",
  "prediabetes_vs_t2d", "benchmark_timing"
)

in_scope <- function(name) {
  lname <- tolower(name)
  !any(vapply(out_of_scope_patterns, function(p) grepl(p, lname, fixed = TRUE), NA))
}

classify <- function(name) {
  if (name %in% c(
    "submit_leptin_cis_trans.sh", "submit_all_proteins_cis_trans.sh",
    "submit_cis_trans_chunked.sh", "submit_decode_lep_bmi.sh",
    "submit_twosampleMR_strata.sh", "submit_reverse_pewas.sh",
    "submit_reverse_mr_r2clump.sh"
  )) return("mr_twosample")

  mr_prefixes <- c(
    "bidirectional_mr_", "twosampleMR_", "reverse_mr_",
    "reverse_observational_pewas", "mr_decode_", "mr_lepr_",
    "mr_summary_table", "generate_MR_", "generate_bidirectional_mr_tables",
    "merge_trig_hdl_bidirectional", "combine_cis_trans_results",
    "combine_reverse_mr_",
    "lepr_bmi_adjusted", "compare_lepr"
  )
  if (any(vapply(mr_prefixes, function(p) startsWith(name, p), NA))) return("mr_twosample")

  fig_prefixes <- c(
    "fig", "funnel_table", "supplemental_tables", "validation_",
    "reproduce_triangulation", "create_master_integrated"
  )
  if (any(vapply(fig_prefixes, function(p) startsWith(name, p), NA))) return("figures_tables")

  "gwas_regenie"
}

safe_copy <- function(from, to) {
  dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
  invisible(file.copy(from, to, overwrite = TRUE))
}

if (dir.exists(out)) unlink(out, recursive = TRUE)
dir.create(out, recursive = TRUE)

scripts_dir <- file.path(root, "scripts")
files <- list.files(scripts_dir, full.names = FALSE)
counts <- c(gwas_regenie = 0L, mr_twosample = 0L, figures_tables = 0L)
manifest <- list()

for (fn in files) {
  fp <- file.path(scripts_dir, fn)
  if (!file.info(fp)$isdir && in_scope(fn)) {
    bucket <- classify(fn)
    rel <- if (bucket == "gwas_regenie") {
      file.path("gwas_regenie", "scripts", fn)
    } else {
      file.path(bucket, "scripts", fn)
    }
    safe_copy(fp, file.path(out, rel))
    counts[bucket] <- counts[bucket] + 1L
    manifest[[length(manifest) + 1L]] <- c(rel, fp)
  }
}

root_extras <- list(
  c("merge_cis_pqtl_coloc.R", "mr_twosample"),
  c("variance_decomposition_bmi.R", "figures_tables"),
  c("run_variance_decomposition.sh", "figures_tables"),
  c("validation_bmi_cis_mr_concordant_plot.R", "figures_tables"),
  c("validation_bmi_cis_mr_plot.R", "figures_tables"),
  c("validation_hba1c_cis_mr_plot.R", "figures_tables")
)
for (row in root_extras) {
  fn <- row[1]
  bucket <- row[2]
  p <- file.path(root, fn)
  if (file.exists(p)) {
    dest <- file.path(out, bucket, fn)
    safe_copy(p, dest)
    manifest[[length(manifest) + 1L]] <- c(file.path(bucket, fn), p)
  }
}

doc_names <- c(
  "README.md", "RUN_PLAN.md", "TEST_RUN_GUIDE.md", "RESUME_GUIDE.md",
  "MONITOR_JOB.md", "STATUS_UPDATE.md", "SUMMARY_STATS_PATHS.md",
  "RESUME_HALLMARKS_FUMA.md"
)
for (name in doc_names) {
  p <- file.path(root, name)
  if (file.exists(p)) {
    rel <- file.path("gwas_regenie", "docs", name)
    safe_copy(p, file.path(out, rel))
    manifest[[length(manifest) + 1L]] <- c(rel, p)
  }
}

gitignore_txt <- c(
  "# --- Suggested additions when merging MAP-D monorepo genetics/MR folders ---",
  "# Large / sensitive outputs (keep local only)",
  "gwas_regenie/results/",
  "gwas_regenie/inputs/",
  "mr_twosample/results/",
  "figures_tables/results/",
  "**/slurm/*.out",
  "**/slurm/*.err",
  "*.regenie.gz",
  "*.bed", "*.bim", "*.fam",
  "# UK Biobank individual-level exports (never commit)",
  "*.RDS", "*.rds", "*.fst",
  "olink_eids*.RDS",
  "# Local secrets (OpenGWAS JWT, etc.)",
  ".secrets/", ""
)
writeLines(gitignore_txt, file.path(out, "SUGGESTED_ROOT_GITIGNORE_APPEND.txt"))

walkthrough <- file.path(root, "WALKTHROUGH_SPLIT_SAMPLE_MR.md")
if (file.exists(walkthrough)) {
  rel <- file.path("docs", "genetics_pipeline", basename(walkthrough))
  safe_copy(walkthrough, file.path(out, rel))
  manifest[[length(manifest) + 1L]] <- c(rel, walkthrough)
}

man_lines <- c(
  "# Staged files for MAP-D Option A monorepo", "",
  sprintf("- **gwas_regenie/scripts/**: %d files", counts["gwas_regenie"]),
  sprintf("- **mr_twosample/scripts/**: %d files", counts["mr_twosample"]),
  sprintf("- **figures_tables/scripts/**: %d files", counts["figures_tables"]),
  "", "## All files", ""
)
mf <- do.call(rbind, manifest)
o <- order(mf[, 1])
for (i in o) {
  man_lines <- c(man_lines, sprintf("- `%s` \u2190 `%s`", mf[i, 1], mf[i, 2]))
}
writeLines(man_lines, file.path(out, "FILE_MANIFEST.md"))

merge_md <- c(
  "# Merging this staging folder into [MAP-D](https://github.com/stejat98/MAP-D)",
  "",
  "This directory was generated from `regenie_pipeline` as **Option A**: add three top-level folders to the existing MAP-D repo.",
  "",
  "**Regenerate this staging folder:** from `regenie_pipeline`, run `module load gcc/14.2.0 R/4.4.2 && Rscript build_github_staging.R`.",
  "",
  "## Layout after merge",
  "",
  "Place contents at the **root** of `MAP-D` (alongside existing `src/`, `main_analysis.R`, etc.):",
  "",
  "| Staged folder | Target on GitHub |",
  "|---------------|------------------|",
  "| `gwas_regenie/` | `gwas_regenie/` |",
  "| `mr_twosample/` | `mr_twosample/` |",
  "| `figures_tables/` | `figures_tables/` |",
  "",
  "## Steps",
  "",
  "1. Clone MAP-D: `git clone https://github.com/stejat98/MAP-D.git && cd MAP-D`",
  "2. Copy (or merge) each staged folder into the repo root.",
  "3. Append `SUGGESTED_ROOT_GITIGNORE_APPEND.txt` to `.gitignore` (edit paths if needed).",
  "4. Refactor hardcoded cluster paths in R/Python scripts to `config.R` or environment variables before pushing.",
  "5. Add `renv` / `sessionInfo()` and document REGENIE + TwoSampleMR versions in README.",
  "",
  "## Not included (by design)",
  "",
  "- `results/`, `inputs/`, `test_results/`, `filtered_snps/`, `py_pkgs/`",
  "- UK Biobank individual-level data",
  "- SLURM `.out` / `.err` logs",
  "- Analyses outside the BMI / HbA1c / TRIG_HDL_RATIO scope",
  "- LDSC tooling: use the official [LDSC](https://github.com/bulik/ldsc) repository with `regenie_to_ldsc_munge_input.py`",
  ""
)
writeLines(merge_md, file.path(out, "MAP_D_MONOREPO_MERGE.md"))

message("Created: ", out)
message("Counts: ", paste(names(counts), counts, sep = "=", collapse = ", "))
