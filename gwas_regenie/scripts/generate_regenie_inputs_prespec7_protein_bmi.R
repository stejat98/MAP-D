#!/usr/bin/env Rscript
# REGENIE inputs for pre-specified Olink proteins in *proteomics-only* EIDs, **BMI in covariates**
# (apples-to-apples with BMI-adjusted ALST hallmarks GWAS for two-sample MR).
#
# Sample split (unchanged from main pipeline):
#   — protein phenotypes: eid ∈ Olink    (this script)
#   — ALST hallmark:      eid ∉ Olink    (independent individuals for MR)
#
# Olink panel note: **IGF1** and **CRP** have no matching `GENE;` row in coding143 — not in this panel;
#   use external pQTL (e.g. deCODE) for MR for those analytes.
#
# Genes with UKB Olink columns (protein_x.<coding>) here:
#   ACAN, FGF21, GDF15, GHR, IGFBP2, IGFBP3, LEP
#
# Resume: by default, skips a phenotype if pheno.txt + covar.txt already exist
#   (avoids re-reading the large main_df and overwriting good outputs). Override with:
#   REGENIE_PRESPEC_FORCE=1  # rewrite all

suppressPackageStartupMessages({
  library(data.table)
})

PIPE <- Sys.getenv("REGENIE_PIPELINE", "/n/groups/patel/sivateja/regenie_pipeline")
MAIN_RDS <- Sys.getenv(
  "MAPD_DATA_RDS",
  "/n/groups/patel/sivateja/UKB/PEWAS_results/data_plus_GLP_complications_glycemic_status_HbA1c_adjusted.RDS"
)
PROTEOMICS_EIDS_RDS <- "/n/groups/patel/sivateja/olink_eids_for_proteins_gwas.RDS"
ADJ_RDATA <- "/home/st320/UKB_PEWAS/adjustments_survival_analysis.Rdata"
OUTPUT_BASE <- file.path(PIPE, "inputs")
STRATUM_NAME <- "full"

# coding column index in main_df = protein_x.<coding> (must exist)
PRESPEC_EXPOSURE_COLS <- c(
  "protein_x.12",   # ACAN
  "protein_x.1039", # FGF21
  "protein_x.1137", # GDF15
  "protein_x.1156", # GHR
  "protein_x.1344", # IGFBP2
  "protein_x.1345", # IGFBP3
  "protein_x.1572"  # LEP
)

`%notin%` <- Negate(`%in%`)

regenie_prespec_force <- function() {
  v <- tolower(Sys.getenv("REGENIE_PRESPEC_FORCE", "0"))
  v %in% c("1", "true", "yes")
}

out_root <- file.path(OUTPUT_BASE, STRATUM_NAME, "proteins_only")
list_file <- file.path(PIPE, "scripts", "pheno_list_prespec7_protein_bmi_full.txt")

dir_nonempty_files <- function(d1, d2) {
  s1 <- tryCatch(file.info(d1)$size, error = function(e) NA_real_)
  s2 <- tryCatch(file.info(d2)$size, error = function(e) NA_real_)
  isTRUE(s1 > 0) && isTRUE(s2 > 0)
}

if (!regenie_prespec_force()) {
  have <- vapply(PRESPEC_EXPOSURE_COLS, function(phen) {
    od <- file.path(out_root, phen)
    pe <- file.path(od, "pheno.txt")
    ce <- file.path(od, "covar.txt")
    file.exists(pe) && file.exists(ce) && dir_nonempty_files(pe, ce)
  }, logical(1L))
  if (all(have)) {
    message(
      "All ", length(PRESPEC_EXPOSURE_COLS), " prespec7 directories present; ",
      "skipping main_df load (set REGENIE_PRESPEC_FORCE=1 to rebuild from RDS)."
    )
    lines <- vapply(PRESPEC_EXPOSURE_COLS, function(phen) {
      if (file.exists(file.path(out_root, phen, "pheno.txt"))) {
        return(paste("full", "proteins_only", phen, sep = "|"))
      }
      NA_character_
    }, character(1L))
    lines <- lines[!is.na(lines)]
    writeLines(lines, list_file)
    message("Pheno list (n=", length(lines), "): ", list_file)
    q("no", status = 0, runLast = FALSE)
  }
}

e0 <- new.env()
load(ADJ_RDATA, envir = e0)
if (exists("adjustments", e0, inherits = FALSE)) {
  ad <- e0$adjustments
} else {
  h <- names(e0)[grepl("adjust", names(e0), ignore.case = TRUE)][1L]
  ad <- e0[[h]]
}
ad <- unique(c(as.character(ad), "x.738"))
ad <- as.character(ad)

cat("Loading main_df...\n")
main_df <- readRDS(MAIN_RDS)
if (!is.data.table(main_df)) main_df <- as.data.table(main_df)

e1 <- new.env()
load(PROTEOMICS_EIDS_RDS, envir = e1)
if (exists("olink_eids_for_proteins_gwas", e1, inherits = FALSE)) {
  olink_eids <- e1$olink_eids_for_proteins_gwas
} else {
  n1 <- names(e1)[grepl("eid|olink", names(e1), ignore.case = TRUE)][1L]
  olink_eids <- e1[[n1]]
}

covariate_cols <- intersect(ad, names(main_df))
if (!"BMI" %in% names(main_df)) stop("BMI must be in main_df for this analysis.")
covariate_cols <- unique(c(covariate_cols, "BMI"))
covariate_cols <- covariate_cols[covariate_cols %in% names(main_df)]
cat("Covariates (incl. BMI): n =", length(covariate_cols), "\n")

df <- main_df
cat("Stratum full N =", nrow(df), "\n")
df_p <- df[eid %in% olink_eids]
cat("proteins_only N =", nrow(df_p), "\n")

create_regenie_inputs_protein_bmi <- function(df_subset, phenotype_name, output_dir) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  if (!phenotype_name %in% names(df_subset)) {
    message("  SKIP: ", phenotype_name, " not in data")
    return(FALSE)
  }
  cat("  ", basename(output_dir), ": BMI in covar (same set as other PEWAS + ALST order)\n", sep = "")
  covar_cols_this <- covariate_cols
  pd <- df_subset[, .(eid, phenotype = get(phenotype_name))]
  pd <- pd[!is.na(phenotype)]
  if (nrow(pd) < 100L) return(FALSE)
  ve <- pd$eid
  covd <- df_subset[eid %in% ve, c("eid", covar_cols_this), with = FALSE]
  vcc <- c()
  for (col in covar_cols_this) {
    if (col %in% names(covd)) {
      cd <- covd[[col]]
      if (!all(is.na(cd)) && length(unique(cd[!is.na(cd)])) > 1) vcc <- c(vcc, col)
    }
  }
  if (length(vcc) > 0) covd <- covd[, c("eid", vcc), with = FALSE] else covd <- covd[, .(eid)]
  md <- merge(pd, covd, by = "eid", all.x = TRUE)
  if (any(duplicated(md$eid))) md <- md[!duplicated(eid)]
  pf <- file.path(output_dir, "pheno.txt")
  po <- md[, .(FID = eid, IID = eid, p = phenotype)]
  setnames(po, "p", phenotype_name)
  fwrite(po, file = pf, sep = " ", quote = FALSE, na = "NA")
  cc <- setdiff(colnames(md), c("eid", "phenotype"))
  if (length(cc) > 0) {
    co <- md[, c("eid", cc), with = FALSE]
    co[, FID := eid]
    co[, IID := eid]
    setcolorder(co, c("FID", "IID", cc))
    co[, eid := NULL]
  } else {
    co <- md[, .(FID = eid, IID = eid)]
  }
  cf <- file.path(output_dir, "covar.txt")
  fwrite(co, file = cf, sep = " ", quote = FALSE, na = "NA")
  message("    Wrote ", nrow(co), " rows, ", ncol(co) - 2L, " covariates + BMI")
  TRUE
}

ok <- 0L
for (phen in PRESPEC_EXPOSURE_COLS) {
  od <- file.path(out_root, phen)
  pf_ex <- file.path(od, "pheno.txt")
  cf_ex <- file.path(od, "covar.txt")
  if (!regenie_prespec_force() && file.exists(pf_ex) && file.exists(cf_ex) && dir_nonempty_files(pf_ex, cf_ex)) {
    message("  RESUME: using existing ", phen, " (set REGENIE_PRESPEC_FORCE=1 to rewrite)")
    ok <- ok + 1L
    next
  }
  if (create_regenie_inputs_protein_bmi(df_p, phen, od)) ok <- ok + 1L
}
cat("\nWrote", ok, "of", length(PRESPEC_EXPOSURE_COLS), "phenotype directories under\n  ", out_root, "\n", sep = "")

lines <- vapply(PRESPEC_EXPOSURE_COLS, function(phen) {
  if (file.exists(file.path(out_root, phen, "pheno.txt"))) {
    return(paste("full", "proteins_only", phen, sep = "|"))
  }
  NA_character_
}, character(1L))
lines <- lines[!is.na(lines)]
writeLines(lines, list_file)
message("Pheno list (n=", length(lines), "): ", list_file)
