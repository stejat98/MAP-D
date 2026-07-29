#!/usr/bin/env Rscript
# FIG2-style atlas: UK Biobank Olink only (no STEP trial data).
# Panel A: Spearman correlation of UKB forward beta vectors (one full proteome per
#           trait×stratum file). Tables are full_joined by Exposure; cor() uses
#           pairwise.complete.obs so each rho uses the overlap for that pair only
#           (not intersection across all nine columns).
# Panel B: LDSC genetic correlations (full glycemic stratum, held-out GWAS sample).
# Values are read from Supplementary_Table_LDSC_genetic_correlations_pairwise.csv (same as S23).
.libPaths(c("/n/groups/patel/sivateja/R_libs", .libPaths()))
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

pewas_dir <- "/n/groups/patel/sivateja/UKB/PEWAS_results"
out_pdf <- "/n/groups/patel/sivateja/regenie_pipeline/FIG2_ENR_ATLAS.pdf"
out_png <- "/n/groups/patel/sivateja/regenie_pipeline/FIG2_ENR_ATLAS.png"
out_b_pdf <- "/n/groups/patel/sivateja/regenie_pipeline/FIG2_LDSC_genetic_correlations_standalone.pdf"
out_b_png <- "/n/groups/patel/sivateja/regenie_pipeline/FIG2_LDSC_genetic_correlations_standalone.png"
ldsc_pairwise <- "/n/groups/patel/sivateja/regenie_pipeline/results/supplemental_tables/Supplementary_Table_LDSC_genetic_correlations_pairwise.csv"

strata_rds <- c("non_T2D", "prediabetes", "T2D")
strata_lab <- c("Normoglycemic", "Prediabetes", "T2D")
suffix <- "adj_fasting_time.RDS"

traits <- c("BMI", "HbA1c", "TRIG_HDL_RATIO")
trait_short <- c("BMI", "HbA1c", "TG/HDL")

pull_forward_betas <- function(phenotype, stratum) {
  fn <- file.path(
    pewas_dir,
    paste0(phenotype, "_Linear_regression_proteomic_lm_results_", stratum, "_", suffix)
  )
  if (!file.exists(fn)) stop("Missing file: ", fn)
  d <- readRDS(fn)
  d <- d[grepl("protein_x", d$term, fixed = TRUE), , drop = FALSE]
  if (nrow(d) == 0) stop("No protein rows in ", fn)
  d %>%
    dplyr::select(Exposure, estimate) %>%
    distinct(Exposure, .keep_all = TRUE)
}

cat("UKB forward PEWAS → proteomic correlation matrix (Panel A)\n")

col_dfs <- list()
for (ti in seq_along(traits)) {
  pheno <- traits[ti]
  tlab <- trait_short[ti]
  for (si in seq_along(strata_rds)) {
    st <- strata_rds[si]
    slab <- strata_lab[si]
    colnm <- paste0(tlab, " | ", slab)
    df <- pull_forward_betas(pheno, st)
    names(df)[2] <- colnm
    col_dfs[[colnm]] <- df
  }
}

wide <- col_dfs[[1]]
for (k in seq_len(length(col_dfs))[-1]) {
  nm <- names(col_dfs)[k]
  wide <- full_join(wide, col_dfs[[nm]], by = "Exposure")
}

X <- as.matrix(wide[, -1, drop = FALSE])
M_prot <- cor(X, method = "spearman", use = "pairwise.complete.obs")
diag(M_prot) <- 1

n_per_col <- colSums(!is.na(X))
Ma <- !is.na(X)
overlap_counts <- crossprod(Ma) # [i,j] = proteins with beta in both column i and j
min_pair_n <- min(overlap_counts)
cat(sprintf("  Proteins per column (non-NA beta): %s\n", paste(range(n_per_col), collapse = "-")))
cat(sprintf("  Smallest pairwise overlap for Spearman rho: %d proteins\n", min_pair_n))

# Panel A: within each glycemic stratum, 3×3 Spearman rho among BMI, HbA1c, TG/HDL (faceted)
prot_long <- lapply(strata_lab, function(slab) {
  idx <- paste0(trait_short, " | ", slab)
  subm <- M_prot[idx, idx, drop = FALSE]
  dimnames(subm) <- list(trait_short, trait_short)
  subm %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Var1") %>%
    pivot_longer(-Var1, names_to = "Var2", values_to = "value") %>%
    mutate(subgroup = slab)
}) %>%
  bind_rows() %>%
  mutate(
    subgroup = factor(subgroup, levels = strata_lab),
    Var1 = factor(Var1, levels = trait_short),
    Var2 = factor(Var2, levels = trait_short)
  )

trait_short_ldsc <- function(x) ifelse(x == "TG_HDL_ratio", "TG/HDL", as.character(x))

rg_pairs <- read.csv(ldsc_pairwise, check.names = FALSE) %>%
  mutate(
    pair_label = paste0(trait_short_ldsc(trait_A_short), "-", trait_short_ldsc(trait_B_short)),
    pair_label = factor(
      pair_label,
      levels = c("BMI-HbA1c", "BMI-TG/HDL", "HbA1c-TG/HDL")
    ),
    label_x = pmax(rg_95CI_upper, genetic_correlation_rg, na.rm = TRUE) + 0.012
  )

heat_theme <- function() {
  theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 8),
      axis.text.y = element_text(color = "black", size = 8),
      legend.position = "right"
    )
}

pA <- ggplot(prot_long, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "grey92", linewidth = 0.35) +
  geom_text(aes(label = sprintf("%.2f", value)), size = 3, color = "grey15") +
  scale_fill_gradient(
    low = "#5E35B1",
    high = "#F9A825",
    limits = c(0, 1),
    oob = scales::squish,
    name = expression("Spearman" ~ rho)
  ) +
  coord_fixed() +
  facet_wrap(~subgroup, nrow = 1, scales = "fixed") +
  labs(x = NULL, y = NULL, tag = "A") +
  heat_theme() +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "grey95", color = NA),
    legend.position = "right",
    plot.tag = element_text(face = "bold")
  )

bar_theme <- function() {
  theme_minimal(base_size = 11) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(color = "black", size = 9),
      axis.text.y = element_text(color = "black", size = 9)
    )
}

xmax_b <- max(rg_pairs$label_x, na.rm = TRUE) * 1.06

pB <- ggplot(rg_pairs, aes(x = genetic_correlation_rg, y = pair_label)) +
  geom_col(fill = "black", orientation = "y", width = 0.62) +
  geom_errorbar(
    aes(xmin = rg_95CI_lower, xmax = rg_95CI_upper),
    orientation = "y",
    width = 0.22,
    linewidth = 0.35,
    color = "grey35"
  ) +
  geom_text(
    aes(x = label_x, label = sprintf("%.4f", genetic_correlation_rg)),
    hjust = 0,
    size = 3,
    color = "grey15"
  ) +
  scale_x_continuous(
    limits = c(0, xmax_b),
    expand = expansion(mult = c(0.01, 0.02)),
    name = expression("LDSC" ~ r[g])
  ) +
  scale_y_discrete(limits = rev(levels(rg_pairs$pair_label))) +
  labs(x = NULL, y = NULL, tag = "B") +
  bar_theme() +
  theme(plot.tag = element_text(face = "bold"))

combined <- pA / pB +
  plot_layout(heights = c(1.05, 0.38))

ggsave(out_pdf, combined, width = 12, height = 11.5, dpi = 300)
ggsave(out_png, combined, width = 12, height = 11.5, dpi = 300)

# Standalone Panel B (LDSC r_g) for Illustrator — same geometry as composite panel B, no panel letter
pB_standalone <- pB + labs(tag = NULL) + theme(plot.tag = element_blank())
ggsave(out_b_pdf, pB_standalone, width = 6, height = 3.25, dpi = 300)
ggsave(out_b_png, pB_standalone, width = 6, height = 3.25, dpi = 300)

n_snps_txt <- paste(sort(unique(rg_pairs$n_SNPs_LDSC)), collapse = ", ")
caption_a <- sprintf(
  paste0(
    "Panel A. Within each glycemic stratum (normoglycemia, prediabetes, type 2 diabetes), ",
    "Spearman correlation (rho) between UK Biobank Olink forward association beta vectors for ",
    "body mass index (BMI), glycated hemoglobin (HbA1c), and triglyceride-to-HDL ratio (TG/HDL); ",
    "each vector is the set of protein-level estimates from phenotype ~ scaled protein models ",
    "(fasting time adjusted). Tables were full-joined on protein identifier; correlations use ",
    "pairwise complete observations across proteins. Smallest pairwise overlap N = %d; ",
    "non-missing betas per column N = %d-%d."
  ),
  min_pair_n,
  min(n_per_col),
  max(n_per_col)
)
caption_b <- paste0(
  "Panel B. Linkage disequilibrium score regression (LDSC) genetic correlations (r_g) between ",
  "GWAS for BMI, HbA1c, and TG/HDL in the UK Biobank held-out (non-Olink) sample, full glycemic ",
  "stratum (REGENIE). Point estimates and 95% confidence intervals are taken from ",
  "Supplementary_Table_LDSC_genetic_correlations_pairwise.csv (HapMap3 intersection); ",
  "N SNPs per LDSC run: ", n_snps_txt, ". Bars are annotated at the right with r_g."
)

cat(
  "\n",
  "================================================================================\n",
  "Figure caption (for manuscript; not on figure)\n",
  "================================================================================\n\n",
  caption_a, "\n\n",
  caption_b, "\n\n",
  "================================================================================\n",
  sep = ""
)

cat(
  "\nWrote:\n  ", out_pdf, "\n  ", out_png, "\n",
  "  ", out_b_pdf, "\n  ", out_b_png, "\n",
  sep = ""
)
