# =============================================================================
# plots.R  ---  MAP-D figure builders (pure functions of data + params)
#
# Every function takes the mapd data list (or a component) and returns a plotly
# or ggplot object. Kept separate from app.R so figures can be rendered / QA'd
# headlessly without the Shiny runtime.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(plotly); library(stringr)
})

# ---- palette (Okabe-Ito, colourblind-safe) ---------------------------------
HALLMARK_COL <- c("BMI" = "#E69F00", "HbA1c" = "#0072B2", "TG/HDL" = "#009E73")
STAGE_COL    <- c("Normoglycemic" = "#56B4E9", "Prediabetes" = "#E69F00", "T2D" = "#D55E00")
QUAD_COL <- c(
  "Reversion (UKB+ / GLP1RA−)"  = "#009E73",  # green  : up in UKB, down on drug
  "Reversion (UKB− / GLP1RA+)"  = "#D55E00",  # red    : down in UKB, up on drug
  "Persistent (UKB+ / GLP1RA+)"      = "#E69F00",  # orange : up in both
  "Persistent (UKB− / GLP1RA−)" = "#0072B2" # blue : down in both
)
CAD_COL <- c("Higher risk (OR>1)" = "#D55E00", "Lower risk (OR<1)" = "#0072B2")

BONF_P <- 0.05 / 26307   # manuscript Bonferroni threshold (1.9e-6)

# ---- shared theme ----------------------------------------------------------
theme_mapd <- function(base = 13) {
  theme_minimal(base_size = base, base_family = "Helvetica") +
    theme(
      plot.title    = element_text(face = "bold", size = base + 2),
      plot.subtitle = element_text(colour = "grey35"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey92"),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      axis.title = element_text(face = "bold")
    )
}

# floor p-values so -log10 is finite / plottable
safe_neglog10 <- function(p) {
  p <- as.numeric(p)
  p[is.na(p)] <- 1
  p <- pmax(p, 1e-300)
  -log10(p)
}

empty_plot <- function(msg = "No data for this selection") {
  plot_ly() |>
    add_annotations(text = msg, x = 0.5, y = 0.5, showarrow = FALSE,
                    font = list(size = 16, color = "grey40")) |>
    layout(xaxis = list(visible = FALSE), yaxis = list(visible = FALSE))
}

as_plotly <- function(p, height = NULL, tooltip = "text") {
  ggplotly(p, tooltip = tooltip, height = height) |>
    layout(legend = list(orientation = "h", y = -0.2)) |>
    config(displaylogo = FALSE,
           modeBarButtonsToRemove = c("lasso2d", "select2d", "autoScale2d"))
}

# =============================================================================
# ATLAS  (Fig 2) -- proteome-wide volcano per hallmark
#   direction: "reverse" (Hallmark -> Protein) or "forward" (Protein -> Hallmark)
# =============================================================================
plot_atlas_volcano <- function(master, hallmark, direction = "reverse",
                               label_top = 12) {
  d <- master |> filter(Hallmark == hallmark)
  if (direction == "reverse") {
    d$beta <- d$Rev_Obs_Beta; d$p <- d$Rev_Obs_P
    xlab <- "Effect size  (β per SD of hallmark)"
    sub  <- "Hallmark → Protein  (reverse observational)"
  } else {
    d$beta <- d$Fwd_Obs_Beta; d$p <- d$Fwd_Obs_P
    xlab <- "Effect size  (β per SD of protein)"
    sub  <- "Protein → Hallmark  (forward observational)"
  }
  d <- d |> filter(!is.na(beta), !is.na(p))
  if (!nrow(d)) return(empty_plot())
  d$nlp <- safe_neglog10(d$p)
  d$Significant <- d$p < BONF_P
  d$dir <- ifelse(d$beta > 0, "Positive", "Negative")
  d$grp <- ifelse(!d$Significant, "Not significant",
                  ifelse(d$beta > 0, "Positive (Bonferroni)", "Negative (Bonferroni)"))
  cols <- c("Not significant" = "grey78",
            "Positive (Bonferroni)" = unname(HALLMARK_COL[hallmark]),
            "Negative (Bonferroni)" = "#555555")
  d$tt <- paste0("<b>", d$Gene, "</b><br>β = ", signif(d$beta, 3),
                 "<br>P = ", signif(d$p, 3))
  top <- d |> filter(Significant) |> slice_max(abs(beta), n = label_top)
  p <- ggplot(d, aes(beta, nlp, colour = grp, text = tt)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_hline(yintercept = -log10(BONF_P), linetype = "dashed", colour = "grey40") +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey60") +
    scale_colour_manual(values = cols, name = NULL) +
    labs(title = paste0(hallmark, " atlas"), subtitle = sub,
         x = xlab, y = expression(-log[10](P))) +
    theme_mapd()
  if (nrow(top)) p <- p + ggrepel_or_text(top)
  as_plotly(p)
}

# text labels without a hard ggrepel dependency
ggrepel_or_text <- function(df) {
  geom_text(data = df, aes(label = Gene), size = 3, colour = "grey15",
            vjust = -0.6, show.legend = FALSE)
}

# ---- Top proteins ranking (lollipop) ---------------------------------------
plot_top_proteins <- function(master, hallmark, direction = "reverse", n = 15) {
  d <- master |> filter(Hallmark == hallmark)
  if (direction == "reverse") { d$beta <- d$Rev_Obs_Beta; d$p <- d$Rev_Obs_P }
  else                        { d$beta <- d$Fwd_Obs_Beta; d$p <- d$Fwd_Obs_P }
  d <- d |> filter(!is.na(beta), p < BONF_P) |> slice_max(abs(beta), n = n)
  if (!nrow(d)) return(empty_plot("No Bonferroni-significant proteins"))
  d$Gene <- factor(d$Gene, levels = d$Gene[order(d$beta)])
  d$tt <- paste0("<b>", d$Gene, "</b><br>β = ", signif(d$beta, 3),
                 "<br>P = ", signif(d$p, 3))
  p <- ggplot(d, aes(beta, Gene, text = tt)) +
    geom_segment(aes(x = 0, xend = beta, y = Gene, yend = Gene), colour = "grey75") +
    geom_point(aes(colour = beta > 0), size = 3) +
    scale_colour_manual(values = c("TRUE" = unname(HALLMARK_COL[hallmark]),
                                   "FALSE" = "#555555"), guide = "none") +
    geom_vline(xintercept = 0, colour = "grey60") +
    labs(title = paste0("Top ", n, " ", hallmark, " associations"),
         x = "Effect size (β)", y = NULL) +
    theme_mapd()
  as_plotly(p, height = max(360, n * 26))
}

# =============================================================================
# ARCHITECTURE (Fig 3A) -- cross-hallmark correlation heatmap
# =============================================================================
plot_corr_heatmap <- function(corr_df, title = "Cross-hallmark correlation") {
  d <- corr_df
  d$lab <- sprintf("%.2f", d$r)
  p <- ggplot(d, aes(Hallmark_A, Hallmark_B, fill = r, text = paste0(
        Hallmark_A, " vs ", Hallmark_B, "<br>r = ", round(r, 3)))) +
    geom_tile(colour = "white", linewidth = 1) +
    geom_text(aes(label = lab), size = 4.2, colour = "grey10") +
    scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00",
                         midpoint = 0, limits = c(-1, 1), name = "Pearson r") +
    labs(title = title, x = NULL, y = NULL) +
    coord_equal() + theme_mapd() +
    theme(panel.grid = element_blank())
  as_plotly(p)
}

# =============================================================================
# TRIANGULATION (Fig 4) -- UKB reverse-obs beta vs STEP semaglutide effect
#   MR-significant proteins drawn as squares. Coloured by reversion/persistence.
# =============================================================================
plot_triangulation <- function(master, hallmark, trial = NULL, only_trial_sig = TRUE) {
  if (is.null(trial)) trial <- if (hallmark == "HbA1c") "STEP2" else "STEP1"
  d <- master |> filter(Hallmark == hallmark)
  d$step  <- if (trial == "STEP2") d$Estimate_STEP2 else d$Estimate_STEP1
  d$stepp <- if (trial == "STEP2") d$Pval_STEP2 else d$Pval_STEP1
  d$ukb   <- d$Rev_Obs_Beta
  d <- d |> filter(!is.na(step), !is.na(ukb))
  if (only_trial_sig) d <- d |> filter(stepp < 0.05)
  if (!nrow(d)) return(empty_plot())
  quad <- function(u, s) dplyr::case_when(
    u > 0 & s < 0 ~ "Reversion (UKB+ / GLP1RA−)",
    u < 0 & s > 0 ~ "Reversion (UKB− / GLP1RA+)",
    u > 0 & s > 0 ~ "Persistent (UKB+ / GLP1RA+)",
    TRUE          ~ "Persistent (UKB− / GLP1RA−)")
  d$Quadrant <- quad(d$ukb, d$step)
  d$MR <- ifelse(d$Rev_MR_Sig %in% TRUE | d$Fwd_MR_Sig %in% TRUE,
                 "MR-supported", "Not MR-significant")
  d$tt <- paste0("<b>", d$Gene, "</b><br>UKB β = ", signif(d$ukb, 3),
                 "<br>", trial, " effect = ", signif(d$step, 3),
                 "<br>", d$Quadrant,
                 ifelse(d$MR == "MR-supported", "<br>■ MR-supported", ""))
  p <- ggplot(d, aes(ukb, step, colour = Quadrant, shape = MR, text = tt)) +
    geom_hline(yintercept = 0, colour = "grey55") +
    geom_vline(xintercept = 0, colour = "grey55") +
    geom_point(alpha = 0.8, size = 2) +
    scale_colour_manual(values = QUAD_COL, name = NULL) +
    scale_shape_manual(values = c("MR-supported" = 15, "Not MR-significant" = 16), name = NULL) +
    labs(title = paste0(hallmark, "  ×  ", trial, " semaglutide trial"),
         subtitle = "Squares = genetic (MR) support · opposite signs = reversion, same = persistent",
         x = "UK Biobank effect  (hallmark → protein β)",
         y = paste0(trial, " semaglutide effect (log2 FC)")) +
    theme_mapd()
  as_plotly(p)
}

# quadrant counts for the triangulation view
triangulation_counts <- function(master, hallmark, trial = NULL, only_trial_sig = TRUE) {
  if (is.null(trial)) trial <- if (hallmark == "HbA1c") "STEP2" else "STEP1"
  d <- master |> filter(Hallmark == hallmark)
  d$step <- if (trial == "STEP2") d$Estimate_STEP2 else d$Estimate_STEP1
  d$stepp <- if (trial == "STEP2") d$Pval_STEP2 else d$Pval_STEP1
  d$ukb  <- d$Rev_Obs_Beta
  d <- d |> filter(!is.na(step), !is.na(ukb))
  if (only_trial_sig) d <- d |> filter(stepp < 0.05)
  if (!nrow(d)) return(data.frame())
  d$class <- ifelse(sign(d$ukb) != sign(d$step), "Reversion", "Persistent")
  d |> count(class, name = "n")
}

# =============================================================================
# BIDIRECTIONAL MR (genetics) -- forest of significant instruments per hallmark
# =============================================================================
plot_mr_forest <- function(master, hallmark, max_n = 30) {
  d <- master |> filter(Hallmark == hallmark)
  fwd <- d |> filter(Fwd_MR_Sig %in% TRUE) |>
    transmute(Gene, beta = Fwd_MR_Beta, se = Fwd_MR_SE, p = Fwd_MR_P,
              Direction = "Protein-first (protein → hallmark)")
  rev <- d |> filter(Rev_MR_Sig %in% TRUE) |>
    transmute(Gene, beta = Rev_MR_Beta, se = Rev_MR_SE, p = Rev_MR_P,
              Direction = "Hallmark-first (hallmark → protein)")
  dd <- bind_rows(fwd, rev) |> filter(!is.na(beta))
  if (!nrow(dd)) return(empty_plot("No genome-wide MR-significant proteins"))
  dd <- dd |> group_by(Direction) |> slice_max(abs(beta / pmax(se, 1e-9)), n = max_n) |> ungroup()
  dd$lo <- dd$beta - 1.96 * dd$se; dd$hi <- dd$beta + 1.96 * dd$se
  dd$lab <- factor(dd$Gene, levels = unique(dd$Gene[order(dd$Direction, dd$beta)]))
  dd$tt <- paste0("<b>", dd$Gene, "</b><br>", dd$Direction,
                  "<br>β = ", signif(dd$beta, 3), "<br>P = ", signif(dd$p, 3))
  p <- ggplot(dd, aes(beta, lab, colour = Direction, text = tt)) +
    geom_vline(xintercept = 0, colour = "grey60") +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, alpha = 0.6) +
    geom_point(size = 2.4) +
    scale_colour_manual(values = c("Protein-first (protein → hallmark)" = "#CC79A7",
                                   "Hallmark-first (hallmark → protein)" = "#0072B2"),
                        name = NULL) +
    labs(title = paste0(hallmark, ": genetically-supported (MR) proteins"),
         x = "MR causal estimate (β)", y = NULL) +
    theme_mapd() + theme(legend.position = "top")
  as_plotly(p, height = max(380, nrow(dd) * 22))
}

# causal-architecture summary bars across hallmarks
plot_causal_architecture <- function(mr_summary) {
  d <- mr_summary |>
    select(Hallmark, `Protein-first (protein → hallmark)` = Forward_MR_sig,
           `Hallmark-first (hallmark → protein)` = Reverse_MR_sig) |>
    pivot_longer(-Hallmark, names_to = "Direction", values_to = "n")
  d$Hallmark <- factor(d$Hallmark, levels = c("BMI", "HbA1c", "TG/HDL"))
  d$tt <- paste0("<b>", d$Hallmark, "</b><br>", d$Direction, ": ", d$n, " proteins")
  p <- ggplot(d, aes(Hallmark, n, fill = Direction, text = tt)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.65) +
    geom_text(aes(label = n), position = position_dodge(width = 0.7),
              vjust = -0.4, size = 3.6, show.legend = FALSE) +
    scale_fill_manual(values = c("Protein-first (protein → hallmark)" = "#CC79A7",
                                 "Hallmark-first (hallmark → protein)" = "#0072B2"), name = NULL) +
    labs(title = "Causal architecture differs by hallmark",
         subtitle = "BMI is hallmark-first dominant · HbA1c is protein-first dominant · TG/HDL is bidirectional",
         x = NULL, y = "MR-significant proteins") +
    theme_mapd()
  as_plotly(p)
}

# MR volcano for a single direction
plot_mr_volcano <- function(master, hallmark, direction = "reverse") {
  d <- master |> filter(Hallmark == hallmark)
  if (direction == "reverse") { d$beta <- d$Rev_MR_Beta; d$p <- d$Rev_MR_P; d$sig <- d$Rev_MR_Sig
    sub <- "Hallmark-first MR (hallmark → protein)" }
  else { d$beta <- d$Fwd_MR_Beta; d$p <- d$Fwd_MR_P; d$sig <- d$Fwd_MR_Sig
    sub <- "Protein-first MR (protein → hallmark)" }
  d <- d |> filter(!is.na(beta), !is.na(p))
  if (!nrow(d)) return(empty_plot("No MR estimates in this direction"))
  d$nlp <- safe_neglog10(d$p)
  d$grp <- ifelse(d$sig %in% TRUE, "MR-significant", "Not significant")
  d$tt <- paste0("<b>", d$Gene, "</b><br>β = ", signif(d$beta, 3),
                 "<br>P = ", signif(d$p, 3))
  top <- d |> filter(grp == "MR-significant") |> slice_max(nlp, n = 12)
  p <- ggplot(d, aes(beta, nlp, colour = grp, text = tt)) +
    geom_point(alpha = 0.65, size = 1.6) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey45") +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey60") +
    scale_colour_manual(values = c("MR-significant" = "#D55E00",
                                   "Not significant" = "grey78"), name = NULL) +
    { if (nrow(top)) geom_text(data = top, aes(label = Gene), size = 3,
                               vjust = -0.6, colour = "grey15", show.legend = FALSE) } +
    labs(title = paste0(hallmark, " MR volcano"), subtitle = sub,
         x = "MR causal estimate (β)", y = expression(-log[10](P))) +
    theme_mapd()
  as_plotly(p)
}

# =============================================================================
# EVIDENCE FUNNEL -- winnowing of proteins across evidence layers
# =============================================================================
plot_funnel <- function(funnel, hallmark) {
  d <- funnel |> filter(Hallmark == hallmark, stage_num <= 4, !is.na(n))
  if (!nrow(d)) return(empty_plot())
  d$stage_short <- factor(d$stage_short, levels = unique(d$stage_short[order(d$stage_num)]))
  d$tt <- paste0("<b>", d$stage_short, "</b><br>", d$Subgroup, ": ", d$n, " proteins")
  p <- ggplot(d, aes(stage_short, n, fill = Subgroup, text = tt)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = n), position = position_dodge(width = 0.8),
              vjust = -0.35, size = 3, show.legend = FALSE) +
    scale_fill_manual(values = STAGE_COL, name = "Glycemic stage") +
    labs(title = paste0(hallmark, ": evidence funnel"),
         subtitle = "Observational → trial → genetics → concordant",
         x = NULL, y = "Proteins retained") +
    theme_mapd() +
    theme(axis.text.x = element_text(angle = 18, hjust = 1))
  as_plotly(p)
}

# =============================================================================
# INCIDENT CAD (Fig 5)
# =============================================================================
plot_cad_volcano <- function(cad, label_top = 12) {
  d <- cad |> filter(!is.na(estimate), !is.na(FDR))
  d$grp <- ifelse(!d$Sig, "Not significant (FDR≥0.05)", d$Direction)
  cols <- c("Not significant (FDR≥0.05)" = "grey80", CAD_COL)
  d$tt <- paste0("<b>", d$Gene, "</b><br>OR = ", round(d$OR, 3),
                 "<br>FDR = ", signif(d$FDR, 3))
  top <- d |> filter(Sig) |> slice_max(abs(estimate), n = label_top)
  p <- ggplot(d, aes(estimate, neglog10FDR, colour = grp, size = neglog10FDR, text = tt)) +
    geom_point(alpha = 0.7) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey45") +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey60") +
    scale_colour_manual(values = cols, name = NULL) +
    scale_size_continuous(range = c(1, 5), guide = "none") +
    { if (nrow(top)) geom_text(data = top, aes(label = Gene), size = 3,
                               vjust = -0.7, colour = "grey15", show.legend = FALSE) } +
    labs(title = "Incident coronary artery disease",
         subtitle = "Association of trial-triangulated proteins with incident CAD (UK Biobank)",
         x = "log odds ratio per SD protein", y = expression(-log[10](FDR))) +
    theme_mapd()
  as_plotly(p)
}

# =============================================================================
# GTEx tissue enrichment
# =============================================================================
plot_gtex <- function(gtex, hallmark, stage, top = 15) {
  d <- gtex |> filter(Hallmark == hallmark, Subgroup == stage)
  if (!nrow(d)) return(empty_plot())
  d <- d |> count(Term, name = "n") |> slice_max(n, n = top)
  d$Term <- factor(d$Term, levels = d$Term[order(d$n)])
  d$tt <- paste0("<b>", d$Term, "</b><br>", d$n, " proteins")
  p <- ggplot(d, aes(n, Term, text = tt)) +
    geom_col(fill = unname(HALLMARK_COL[hallmark]), width = 0.7) +
    labs(title = paste0(hallmark, " · ", stage),
         subtitle = "Tissue-enriched sources (GTEx)", x = "Proteins", y = NULL) +
    theme_mapd()
  as_plotly(p, height = max(340, nrow(d) * 26))
}

# =============================================================================
# DrugBank directional bar
# =============================================================================
plot_drugbank <- function(drugbank, top = 20) {
  d <- drugbank |> distinct(gene_symbol, OR, Direction, action_needed) |>
    slice_max(abs(log(OR)), n = top)
  d$gene_symbol <- factor(d$gene_symbol, levels = d$gene_symbol[order(d$OR)])
  d$tt <- paste0("<b>", d$gene_symbol, "</b><br>OR = ", round(d$OR, 3),
                 "<br>action: ", d$action_needed)
  p <- ggplot(d, aes(OR, gene_symbol, fill = Direction, text = tt)) +
    geom_col(width = 0.7) +
    geom_vline(xintercept = 1, colour = "grey55") +
    scale_fill_manual(values = CAD_COL, name = NULL) +
    labs(title = "Druggable CAD-associated proteins",
         subtitle = "Directionally matched approved drugs (DrugBank)",
         x = "Odds ratio (incident CAD)", y = NULL) +
    theme_mapd()
  as_plotly(p, height = max(360, nrow(d) * 26))
}
