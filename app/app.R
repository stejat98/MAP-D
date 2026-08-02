# =============================================================================
# MAP-D : Metabolic Atlas of the Proteome in Diabetes
# Interactive companion to Tangirala, Isaac, ... Tierney & Patel (Nature Metabolism)
#
# Guides the reader through the manuscript's triangulation narrative (Fig 1 -> 5):
#   observational atlas  ->  semaglutide trials  ->  bidirectional MR genetics
#   ->  evidence funnel  ->  incident CAD  ->  drug targeting  ->  protein explorer
#
# Data: data/mapd.rds  (built from SUPP_TABLE_1 by build_data.R -- single source
# of truth). Rebuild with:  Rscript build_data.R
# =============================================================================

suppressPackageStartupMessages({
  library(shiny); library(bslib); library(dplyr); library(tidyr)
  library(plotly); library(DT); library(stringr)
})

source("plots.R")
m <- readRDS("data/mapd.rds")
HM <- m$meta$hallmarks          # BMI, HbA1c, TG/HDL
SG <- m$meta$stages             # Normoglycemic, Prediabetes, T2D

# ---- small UI helpers ------------------------------------------------------
guide_box <- function(...) {
  div(class = "guide-box",
      div(class = "guide-icon", "\U0001F9ED"),   # compass
      div(...))
}
hallmark_pick <- function(id, sel = "BMI")
  radioButtons(id, "Hallmark", choices = HM, selected = sel, inline = TRUE)
dir_pick <- function(id, sel = "reverse")
  radioButtons(id, "Direction",
               choices = c("Hallmark → Protein (reverse)" = "reverse",
                           "Protein → Hallmark (forward)" = "forward"),
               selected = sel)

# =============================================================================
# THEME
# =============================================================================
theme <- bs_theme(
  version = 5,
  bg = "#ffffff", fg = "#1c2733",
  primary = "#0072B2", secondary = "#009E73",
  base_font  = font_google("Inter"),
  heading_font = font_google("Inter"),
  "navbar-bg" = "#12263a"
) |>
  bs_add_rules("
    .guide-box{display:flex;gap:14px;background:#eef4fb;border-left:5px solid #0072B2;
      border-radius:10px;padding:14px 18px;margin:6px 0 20px 0;font-size:0.95rem;color:#2a3b4d;}
    .guide-icon{font-size:1.5rem;line-height:1.2;}
    .hero{background:linear-gradient(120deg,#12263a 0%,#0072B2 60%,#009E73 120%);
      color:#fff;padding:34px 40px;border-radius:16px;margin-bottom:22px;}
    .hero h1{font-weight:800;margin:0;font-size:2.3rem;letter-spacing:-.5px;}
    .hero p{opacity:.92;margin:.5rem 0 0 0;max-width:960px;}
    .card{border:none;box-shadow:0 2px 14px rgba(20,40,70,.08);border-radius:14px;}
    .arch-card{border-radius:12px;padding:16px 18px;height:100%;}
    .section-note{color:#5a6b7b;font-size:.9rem;}
  ")

# =============================================================================
# UI
# =============================================================================
ui <- page_navbar(
  title = span(strong("MAP-D"), span(" · Metabolic Atlas of the Proteome in Diabetes",
                                     style = "font-weight:400;opacity:.85;")),
  theme = theme, id = "nav", fillable = FALSE,

  # ---------------------------------------------------------------- OVERVIEW
  nav_panel(
    "Overview", icon = icon("compass"),
    div(class = "hero",
        h1("Mapping the proteome from normoglycemia to type 2 diabetes"),
        p("MAP-D integrates cross-sectional proteomics, semaglutide trials, and genetics ",
          "(Mendelian randomization) to triangulate the circulating proteins that track ",
          "three cardiometabolic hallmarks — adiposity (BMI), glycemia (HbA1c), and ",
          "insulin resistance (TG/HDL) — across glycemic stages, and asks which are ",
          "reversed by GLP-1 receptor agonists versus persistently dysregulated.")),
    layout_column_wrap(
      width = 1/4, class = "mb-2",
      value_box("Proteins", m$meta$n_proteins, showcase = icon("dna"),
                theme = "primary", p("Olink Explore 3072")),
      value_box("Hallmarks × stages", "3 × 3", showcase = icon("layer-group"),
                theme = "secondary", p("BMI · HbA1c · TG/HDL")),
      value_box("CAD-tested proteins", m$meta$n_cad_tested, showcase = icon("heart-pulse"),
                theme = "primary", p(paste0(m$meta$n_cad_sig, " FDR-significant"))),
      value_box("Druggable hits", m$meta$n_drug_pairs, showcase = icon("pills"),
                theme = "secondary", p(paste0(m$meta$n_drugs, " approved drugs")))
    ),
    layout_column_wrap(
      width = 1/2,
      card(card_header("How to read this atlas"),
        card_body(
          p(strong("The central idea is triangulation."),
            " A protein is most credible when three independent lines of evidence agree:"),
          tags$ol(
            tags$li(strong("Observational"), " — the protein correlates with a hallmark in UK Biobank (the ", strong("Atlas"), " tab)."),
            tags$li(strong("Interventional"), " — semaglutide changes it in the STEP trials (the ", strong("Triangulation & Trials"), " tab)."),
            tags$li(strong("Genetic"), " — Mendelian randomization supports a causal direction (the ", strong("Bidirectional MR"), " tab).")),
          p("The ", strong("Evidence Funnel"), " shows how proteins narrow as evidence is stacked; ",
            strong("Persistent → CAD"), " and ", strong("Drug Targeting"),
            " follow the therapeutically persistent proteins to disease risk and to approved drugs. ",
            "Use ", strong("Protein Explorer"), " to pull every layer of evidence for any single protein."))),
      card(card_header("Causal architecture at a glance"),
        card_body(
          plotlyOutput("ov_causal", height = "260px"),
          p(class = "section-note",
            "Adiposity predominantly reshapes the proteome (reverse-dominant); glycemia is ",
            "protein-driven (forward-dominant); insulin resistance shows bidirectional feedback.")))
    ),
    card(card_header("Study design — analysis layers (Supplementary Table 1)"),
         card_body(DTOutput("ov_specs")))
  ),

  # ---------------------------------------------------------------- ATLAS
  nav_panel(
    "Atlas", icon = icon("map"),
    guide_box(strong("The observational atlas (Fig 2–3)."),
      " Each point is one of 2,923 proteins. The volcano shows effect size vs significance for a hallmark; ",
      "toggle the direction of the model. The architecture panel shows how the three hallmarks' ",
      "proteomic signatures relate — observationally and genetically."),
    layout_sidebar(
      sidebar = sidebar(width = 300,
        hallmark_pick("atlas_hm"), dir_pick("atlas_dir"),
        sliderInput("atlas_topn", "Top proteins to rank", 5, 40, 15, 5),
        hr(), p(class = "section-note",
          "Dashed line = Bonferroni threshold (p < 1.9×10⁻⁶). Hover any point for the gene, β, and P.")),
      layout_column_wrap(width = 1/2,
        card(full_screen = TRUE, card_header("Volcano"), plotlyOutput("atlas_volcano", height = "460px")),
        card(full_screen = TRUE, card_header("Top-ranked proteins"), plotlyOutput("atlas_top", height = "460px")))
    ),
    card(card_header("Cross-hallmark correlation (proteome-wide)"),
         plotlyOutput("atlas_corr", height = "420px"),
         card_footer(class = "section-note",
           "Pearson correlation of hallmark→protein effect sizes across all 2,923 proteins."))
  ),

  # ---------------------------------------------------------- TRIANGULATION
  nav_panel(
    "Triangulation & Trials", icon = icon("arrows-to-circle"),
    guide_box(strong("Do observation and therapy agree? (Fig 4)."),
      " X = the protein’s hallmark association in UK Biobank; Y = its change under semaglutide in the STEP trial. ",
      strong("Opposite signs = “reversion”"), " (therapy normalizes it); ", strong("same signs = “persistent”"),
      " (therapy does not). ", strong("Squares"), " mark proteins with genetic (MR) support."),
    layout_sidebar(
      sidebar = sidebar(width = 300,
        hallmark_pick("tri_hm"),
        radioButtons("tri_trial", "Semaglutide trial",
                     c("STEP 1 (obesity, non-T2D)" = "STEP1", "STEP 2 (T2D)" = "STEP2"),
                     selected = "STEP1"),
        checkboxInput("tri_sig", "Only proteins modulated by semaglutide (trial P < 0.05)", TRUE),
        hr(), htmlOutput("tri_counts")),
      card(full_screen = TRUE, card_header("UK Biobank vs semaglutide"),
           plotlyOutput("tri_plot", height = "560px"))
    )
  ),

  # ------------------------------------------------------------ MR GENETICS
  nav_panel(
    "Bidirectional MR", icon = icon("dna"),
    guide_box(strong("Genetics resolves cause from consequence."),
      " Split-sample Mendelian randomization tests two directions: ",
      strong("protein-first"), " (genetically-predicted protein → hallmark) and ",
      strong("hallmark-first"), " (genetically-predicted hallmark → protein). ",
      "The balance of the two reveals each hallmark’s causal architecture."),
    layout_column_wrap(width = 1/3,
      div(class = "arch-card", style = "background:#fdf1e3;",
          h5("BMI — adiposity"), p(strong("Hallmark-first dominant.")),
          p(class="section-note","106 hallmark-first proteins (reshaped by adiposity) vs 10 protein-first — the body’s fat mass broadly remodels the proteome.")),
      div(class = "arch-card", style = "background:#e7f0fb;",
          h5("HbA1c — glycemia"), p(strong("Protein-first dominant.")),
          p(class="section-note","19 protein-first proteins causally influence glycemia and none are altered by it — upstream proteins drive glucose control.")),
      div(class = "arch-card", style = "background:#e6f5ef;",
          h5("TG/HDL — insulin resistance"), p(strong("Bidirectional.")),
          p(class="section-note","17 protein-first and 34 hallmark-first signals, with proteins (e.g. SHBG, PLTP) showing feedback in both directions."))),
    layout_sidebar(
      sidebar = sidebar(width = 300,
        hallmark_pick("mr_hm"),
        radioButtons("mr_dir", "Volcano direction",
                     c("Hallmark-first (hallmark → protein)" = "reverse",
                       "Protein-first (protein → hallmark)" = "forward"), selected = "reverse"),
        hr(), p(class = "section-note",
          "Forest shows genome-wide MR-significant proteins with 95% CI. Volcano shows all tested proteins in the chosen direction.")),
      layout_column_wrap(width = 1/2,
        card(full_screen = TRUE, card_header("MR-supported proteins (forest)"),
             plotlyOutput("mr_forest", height = "520px")),
        card(full_screen = TRUE, card_header("MR volcano"),
             plotlyOutput("mr_volcano", height = "520px")))
    ),
    card(card_header("Causal architecture across hallmarks"),
         plotlyOutput("mr_causal", height = "300px"))
  ),

  # --------------------------------------------------------- EVIDENCE FUNNEL
  nav_panel(
    "Evidence Funnel", icon = icon("filter"),
    guide_box(strong("Stacking evidence winnows candidates."),
      " From all Bonferroni-significant observational associations, each additional requirement — ",
      "trial modulation, then genetic (MR) support, then directional concordance — retains fewer, ",
      "higher-confidence proteins. Bars are split by glycemic stage."),
    layout_sidebar(
      sidebar = sidebar(width = 280, hallmark_pick("fn_hm", "BMI"),
        p(class = "section-note", "Only BMI and HbA1c carry the full funnel in the workbook.")),
      card(full_screen = TRUE, card_header("Evidence funnel"),
           plotlyOutput("fn_plot", height = "480px"))),
    card(card_header("Funnel counts"), card_body(DTOutput("fn_table")))
  ),

  # -------------------------------------------------------------- CAD
  nav_panel(
    "Persistent → CAD", icon = icon("heart-pulse"),
    guide_box(strong("Do persistent proteins carry residual risk? (Fig 5)."),
      " Each trial-triangulated protein is tested for association with incident coronary artery ",
      "disease in UK Biobank. Points right of centre (OR > 1) mark higher risk. ",
      "The strongest signals (e.g. SCARB2, TAFA5, PLTP, LRRN1) index inflammatory and stress pathways."),
    layout_sidebar(
      sidebar = sidebar(width = 280,
        sliderInput("cad_label", "Label top proteins", 0, 25, 12, 1),
        p(class = "section-note", paste0(m$meta$n_cad_sig, " of ",
          m$meta$n_cad_tested, " proteins are FDR-significant."))),
      card(full_screen = TRUE, card_header("Incident CAD volcano"),
           plotlyOutput("cad_plot", height = "520px"))),
    card(card_header("CAD associations"), card_body(DTOutput("cad_table")))
  ),

  # -------------------------------------------------------- DRUG TARGETING
  nav_panel(
    "Drug Targeting", icon = icon("pills"),
    guide_box(strong("From risk protein to approved drug."),
      " CAD-associated proteins are matched to approved drugs with a directionally appropriate ",
      "action — ", strong("inhibitors/antagonists"), " for higher-risk (OR>1) proteins, ",
      strong("agonists/activators"), " for protective (OR<1) proteins — nominating repurposing candidates."),
    layout_column_wrap(width = 1/1,
      card(full_screen = TRUE, card_header("Druggable CAD proteins"),
           plotlyOutput("drug_plot", height = "480px"))),
    card(card_header("Drug–protein matches (DrugBank)"), card_body(DTOutput("drug_table")))
  ),

  # ------------------------------------------------------- PROTEIN EXPLORER
  nav_panel(
    "Protein Explorer", icon = icon("magnifying-glass-chart"),
    guide_box(strong("Every line of evidence for one protein."),
      " Search a gene symbol to assemble its full evidence card: observational associations across ",
      "hallmarks, MR causal estimates, semaglutide trial effects, incident-CAD risk, drug matches, and tissue source."),
    layout_sidebar(
      sidebar = sidebar(width = 320,
        selectizeInput("pe_gene", "Protein / gene", choices = NULL,
                       options = list(placeholder = "Type a gene, e.g. LEP, GHR, SHBG, EGFR")),
        uiOutput("pe_headline")),
      card(card_header("Observational & genetic evidence, by hallmark"),
           card_body(DTOutput("pe_master"))),
      layout_column_wrap(width = 1/2,
        card(card_header("Incident CAD"), card_body(uiOutput("pe_cad"))),
        card(card_header("Drug matches"), card_body(DTOutput("pe_drug")))),
      card(card_header("Tissue sources (GTEx)"), card_body(uiOutput("pe_gtex")))
    )
  ),

  # ---------------------------------------------------------- DATA & METHODS
  nav_panel(
    "Data & Methods", icon = icon("book"),
    layout_column_wrap(width = 1/2,
      card(card_header("Download the atlas"),
        card_body(
          p("All figures are computed from the manuscript supplementary workbook."),
          downloadButton("dl_master", "Integrated evidence (obs + MR + trials)", class = "btn-primary mb-2"), br(),
          downloadButton("dl_cad", "Incident-CAD associations", class = "btn-primary mb-2"), br(),
          downloadButton("dl_drug", "DrugBank matches", class = "btn-primary mb-2"), br(),
          downloadButton("dl_gtex", "GTEx tissue enrichment", class = "btn-primary"))),
      card(card_header("Resources & citation"),
        card_body(
          p(strong("Live atlas: "), a("chiragjp-mapd.share.connect.posit.cloud",
              href = "https://chiragjp-mapd.share.connect.posit.cloud/", target = "_blank")),
          p(strong("Summary statistics: "), a("figshare (MAP-D atlas)",
              href = "https://figshare.com/articles/dataset/MAP-D_protein-wide_association_summary_statistics/30007306?file=67080170", target = "_blank")),
          p(strong("GWAS summary statistics: "), a("figshare (BMI/HbA1c/TG-HDL)",
              href = "https://doi.org/10.6084/m9.figshare.32048550", target = "_blank")),
          p(strong("Code: "), a("github.com/stejat98/MAP-D", href = "https://github.com/stejat98/MAP-D", target = "_blank")),
          hr(),
          p(class = "section-note",
            "Tangirala, Isaac, Gehad, ... Tierney & Patel. Assessing Type 2 Diabetes and GLP-1 ",
            "agonist response trajectories with a proteomic atlas of disease progression.")))),
    card(card_header("Glossary"),
      card_body(
        tags$dl(
          tags$dt("Hallmark → Protein (reverse)"), tags$dd("Protein modeled as outcome of the hallmark; β = protein change per SD of hallmark. Comparable to trial effects."),
          tags$dt("Protein → Hallmark (forward)"), tags$dd("Hallmark modeled as outcome; identifies proteins whose levels may influence the hallmark."),
          tags$dt("Reversion vs Persistent"), tags$dd("Whether semaglutide reverses (opposite sign) or maintains (same sign) the observational direction."),
          tags$dt("MR (Mendelian randomization)"), tags$dd("Uses genetic instruments to infer causal direction, independent of drug exposure."))))
  )
)

# =============================================================================
# SERVER
# =============================================================================
server <- function(input, output, session) {

  updateSelectizeInput(session, "pe_gene", choices = m$search_index,
                       selected = "LEP", server = TRUE)

  dt <- function(df, ...) datatable(df, rownames = FALSE, filter = "top",
      options = list(pageLength = 10, scrollX = TRUE, dom = "tip"), ...)

  # clean scientific formatting for p-values / betas (keeps NA and 0-underflow tidy)
  fmt_p <- function(x) {
    x <- as.numeric(x)
    ifelse(is.na(x), "—",
      ifelse(x == 0, "<1e-300", formatC(x, format = "g", digits = 3)))
  }
  fmt_n <- function(x, d = 3) {
    x <- as.numeric(x)
    ifelse(is.na(x), "—", formatC(signif(x, d), format = "g", digits = d))
  }

  # ---- Overview ----
  output$ov_causal <- renderPlotly(plot_causal_architecture(m$mr_summary))
  output$ov_specs  <- renderDT(datatable(m$specs, rownames = FALSE,
      options = list(pageLength = 5, scrollX = TRUE, dom = "t")))

  # ---- Atlas ----
  output$atlas_volcano <- renderPlotly(plot_atlas_volcano(m$master, input$atlas_hm, input$atlas_dir))
  output$atlas_top     <- renderPlotly(plot_top_proteins(m$master, input$atlas_hm, input$atlas_dir, input$atlas_topn))
  output$atlas_corr    <- renderPlotly(plot_corr_heatmap(m$hallmark_corr_all))

  # ---- Triangulation ----
  output$tri_plot <- renderPlotly(
    plot_triangulation(m$master, input$tri_hm, input$tri_trial, input$tri_sig))
  output$tri_counts <- renderUI({
    cts <- triangulation_counts(m$master, input$tri_hm, input$tri_trial, input$tri_sig)
    if (!nrow(cts)) return(p("No proteins for this selection."))
    rev_n <- cts$n[cts$class == "Reversion"]; per_n <- cts$n[cts$class == "Persistent"]
    tagList(h6("Protein counts"),
      p(span(style = "color:#009E73;font-weight:700;", "Reversion: "),
        if (length(rev_n)) rev_n else 0),
      p(span(style = "color:#E69F00;font-weight:700;", "Persistent: "),
        if (length(per_n)) per_n else 0))
  })

  # ---- MR ----
  output$mr_forest  <- renderPlotly(plot_mr_forest(m$master, input$mr_hm))
  output$mr_volcano <- renderPlotly(plot_mr_volcano(m$master, input$mr_hm, input$mr_dir))
  output$mr_causal  <- renderPlotly(plot_causal_architecture(m$mr_summary))

  # ---- Funnel ----
  output$fn_plot  <- renderPlotly(plot_funnel(m$funnel, input$fn_hm))
  output$fn_table <- renderDT({
    d <- m$funnel |> filter(Hallmark == input$fn_hm) |>
      select(Stage = stage_short, `Glycemic stage` = Subgroup, Proteins = value) |>
      arrange(`Glycemic stage`)
    dt(d)
  })

  # ---- CAD ----
  output$cad_plot  <- renderPlotly(plot_cad_volcano(m$cad, input$cad_label))
  output$cad_table <- renderDT({
    d <- m$cad |> transmute(Gene, Protein = Desc, OR = round(OR, 3),
        `OR 95% CI` = paste0(round(OR_lo, 2), "–", round(OR_hi, 2)),
        P = fmt_p(p.value), FDR = fmt_p(FDR),
        Direction, `FDR<0.05` = Sig, N = SampleSize)
    dt(d)
  })

  # ---- Drugs ----
  output$drug_plot  <- renderPlotly(plot_drugbank(m$drugbank))
  output$drug_table <- renderDT({
    d <- m$drugbank |> transmute(Gene = gene_symbol, Protein = protein_name,
        OR = round(OR, 3), Direction, `Drug` = drug_name,
        `Drug action` = actions, `Action needed` = action_needed,
        DrugBank = drugbank_id, FDR = fmt_p(FDR))
    dt(d)
  })

  # ---- Protein explorer ----
  pe <- reactive({ req(input$pe_gene); input$pe_gene })

  output$pe_headline <- renderUI({
    g <- pe(); mm <- m$master |> filter(Gene == g)
    if (!nrow(mm)) return(p("No record."))
    cc <- unique(na.omit(mm$Causal_Class)); cc <- cc[cc != "Neither"]
    cadrow <- m$cad |> filter(Gene == g)
    tagList(
      h4(g), p(class = "section-note", mm$Gene[1]),
      if (length(cc)) p(strong("Causal class: "), paste(cc, collapse = ", ")) else NULL,
      if (nrow(cadrow)) p(strong("Incident CAD OR: "),
        sprintf("%.2f (FDR %.2g)", cadrow$OR[1], cadrow$FDR[1])) else NULL
    )
  })

  output$pe_master <- renderDT({
    g <- pe()
    d <- m$master |> filter(Gene == g) |>
      transmute(Hallmark,
        `Obs β (rev)` = fmt_n(Rev_Obs_Beta), `Obs P (rev)` = fmt_p(Rev_Obs_P),
        `Obs β (fwd)` = fmt_n(Fwd_Obs_Beta), `Obs P (fwd)` = fmt_p(Fwd_Obs_P),
        `MR β (rev)` = fmt_n(Rev_MR_Beta), `MR P (rev)` = fmt_p(Rev_MR_P), `MR sig (rev)` = Rev_MR_Sig,
        `MR β (fwd)` = fmt_n(Fwd_MR_Beta), `MR P (fwd)` = fmt_p(Fwd_MR_P), `MR sig (fwd)` = Fwd_MR_Sig,
        STEP1 = fmt_n(Estimate_STEP1), STEP2 = fmt_n(Estimate_STEP2),
        `Causal class` = Causal_Class)
    datatable(d, rownames = FALSE, options = list(dom = "t", scrollX = TRUE))
  })

  output$pe_cad <- renderUI({
    g <- pe(); d <- m$cad |> filter(Gene == g)
    if (!nrow(d)) return(p(class = "section-note", "Not among trial-triangulated CAD-tested proteins."))
    tagList(
      p(strong("Odds ratio: "), sprintf("%.3f (95%% CI %.2f–%.2f)", d$OR[1], d$OR_lo[1], d$OR_hi[1])),
      p(strong("FDR: "), signif(d$FDR[1], 3), " · ", d$Direction[1]),
      p(class = "section-note", paste0("n = ", d$SampleSize[1])))
  })

  output$pe_drug <- renderDT({
    g <- pe(); d <- m$drugbank |> filter(gene_symbol == g) |>
      transmute(Drug = drug_name, Action = actions, `Needed` = action_needed, DrugBank = drugbank_id)
    if (!nrow(d)) d <- data.frame(Drug = "No approved-drug match", Action = "", Needed = "", DrugBank = "")
    datatable(d, rownames = FALSE, options = list(dom = "t"))
  })

  output$pe_gtex <- renderUI({
    g <- pe(); d <- m$gtex |> filter(Gene == g) |> distinct(Term) |> arrange(Term)
    if (!nrow(d)) return(p(class = "section-note", "No tissue-specific enrichment recorded."))
    p(paste(gsub("_", " ", d$Term), collapse = " · "))
  })

  # ---- Downloads ----
  dl <- function(name, data) downloadHandler(
    filename = function() paste0(name, ".csv"),
    content = function(f) write.csv(data, f, row.names = FALSE))
  output$dl_master <- dl("mapd_integrated_evidence", m$master)
  output$dl_cad    <- dl("mapd_incident_cad", m$cad)
  output$dl_drug   <- dl("mapd_drugbank", m$drugbank)
  output$dl_gtex   <- dl("mapd_gtex", m$gtex)
}

shinyApp(ui, server)
