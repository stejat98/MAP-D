# MAP-D: Metabolic Atlas of Progression in Diabetes  

This repository contains the code used to build **MAP-D**, the **Metabolic Atlas of Progression in Diabetes**.  
MAP-D integrates UK Biobank proteomic and clinical data with STEP 1/2 semaglutide intervention trials to chart stage-resolved protein associations across metabolic disease progression. It identifies treatment-responsive vs. treatment-intransigent proteins, builds predictive models, and maps therapeutic opportunities using DrugBank.  

---

## Repository structure

```
mapd/
├─ README.md
├─ LICENSE
├─ renv.lock                     # R environment snapshot
├─ config/
│  └─ config.yaml                # File paths & parameters
├─ data/
│  ├─ raw/                       # Unprocessed inputs (UKB, Olink, DrugBank XML)
│  ├─ interim/                   # Intermediate merged/clean datasets
│  └─ processed/                 # Analysis-ready data
├─ results/
│  ├─ tables/                    # CSV/XLSX outputs (main + supplementary)
│  └─ figures/                   # Figures (main + supplementary)
├─ scripts/                      # Analysis modules
│  ├─ 00_setup_env.R
│  ├─ 01_ingest_ukb.R
│  ├─ 02_ingest_olink.R
│  ├─ 03_merge_clean.R
│  ├─ 04_pwas_by_group.R
│  ├─ 05_step_validation.R
│  ├─ 06_pca_kmeans.R
│  ├─ 07_prediction_lasso.R
│  ├─ 08_quadrant_mapping.R
│  ├─ 09_enrichment.R
│  ├─ 10_outcomes_incident.R
│  ├─ 11_drugbank_parse_xml.R
│  ├─ 12_drug_prioritization.R
│  ├─ 13_figures_main.R
│  ├─ 14_tables_main.R
│  └─ utils_common.R
├─ pipeline/
│  ├─ Makefile                   # Workflow orchestration
└─ tests/
   └─ test_utils.R
```

---

## Quick start

### 1. Restore environment
```r
install.packages("renv")
renv::init()
```

### 2. Configure
Edit `config/config.yaml` with your paths and parameters:
```yaml
paths:
  ukb_root: "/path/to/ukb/"
  olink_file: "data/raw/olink_data.txt"
  drugbank_xml: "data/raw/drugbank_full_database.xml"

analysis:
  phenotypes: ["BMI","HDL","LDL","TRIG_HDL_RATIO","DBP","SBP","HbA1c"]
  subgroups: ["Normoglycemic","Prediabetes","T2D"]
  train_frac: 0.70
  seed: 123
  lasso:
    nfolds: 10
    alpha: 1
```

### 3. Run pipeline
```bash
make all
```

---

## Scripts overview

- **04_pwas_by_group.R** → PWAS stratified by subgroup × phenotype.  
- **05_step_validation.R** → Map UKB vs STEP 1/2 proteomic changes.  
- **07_prediction_lasso.R** → Predictive modeling (demographic vs proteomic vs combined).  
- **08_quadrant_mapping.R** → Identify reversion vs intransigent proteins.  
- **10_outcomes_incident.R** → Incident CAD, CKD, NAFLD association testing.  
- **11_drugbank_parse_xml.R** → Parse DrugBank XML for drug–protein mapping.  
- **12_drug_prioritization.R** → Directionally matched drug prioritization.  

---

## Outputs

- **Main Tables**: `results/tables/`  
- **Main Figures**: `results/figures/`  
- **Drug–protein matches**: `results/tables/directionally_matched_drug_hits.csv`  
- **Cohort summary (Table 1)**: `results/tables/table1_summary.csv`  

---

## Data access

- **UK Biobank** and **Olink** raw data cannot be redistributed.  
- Place files under `data/raw/` as specified in `config.yaml`.  
- **DrugBank XML** requires a license and is not included.  

---

## Reproducibility

- Environment captured in `renv.lock`.  
- Deterministic random seeds.  
- Configurable through `config/config.yaml`.  
- Orchestration supported via **Make** or **{targets}**.  
