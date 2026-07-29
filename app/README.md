# MAP-D — Metabolic Atlas of the Proteome in Diabetes

Interactive companion to the Nature Metabolism manuscript
*"Assessing Type 2 Diabetes and GLP-1 agonist response trajectories with a proteomic
atlas of disease progression"* (Tangirala, Isaac, … Tierney & Patel).

The app guides a reader through the paper's triangulation narrative (Fig 1 → 5):
observational atlas → semaglutide trials → **bidirectional MR genetics** →
evidence funnel → incident CAD → drug targeting, plus a per-protein evidence explorer.

## Files

| File | Purpose |
|---|---|
| `build_data.R` | Converts `SUPP_TABLE_1_05_20_2026.xlsx` (single source of truth) → `data/mapd.rds`. Asserts the manuscript's headline numbers. |
| `data/mapd.rds` | Pre-built data bundle the app reads (only ~0.9 MB). |
| `plots.R` | All figure builders (pure functions; render without the Shiny runtime). |
| `app.R` | The Shiny app (bslib navbar, 9 sections). |
| `deploy.R` | Deploys the 3 runtime files to shinyapps.io. |
| `app.R.backup` | The previous (7-phenotype) app, preserved. |

## Rebuild the data (only if the workbook changes)

```r
Rscript build_data.R      # writes data/mapd.rds, prints sanity checks
```

## Run locally

```r
Rscript -e "shiny::runApp('.', port=7654, launch.browser=TRUE)"
```

Then open http://127.0.0.1:7654 . Packages: shiny, bslib, plotly, DT, dplyr,
tidyr, stringr, ggplot2.

## Deploy (public app)

```r
Rscript deploy.R          # -> https://btierneyshiny.shinyapps.io/mapd-visualizer/
```

Only `app.R`, `plots.R`, and `data/mapd.rds` are bundled — the large `data_packet/`
and the raw workbook are **not** uploaded.

## Data provenance

Every figure is computed from the manuscript supplementary workbook
(`SUPP_TABLE_1_05_20_2026.xlsx`): integrated forward/reverse observational +
Mendelian-randomization evidence, STEP 1/2 semaglutide effects, incident-CAD
associations, DrugBank matches, GTEx tissue enrichment, and LDSC genetic
correlations. Scope is the manuscript's three cardiometabolic hallmarks — BMI,
HbA1c, TG/HDL — across normoglycemia, prediabetes, and incident T2D.
