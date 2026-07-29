# Split-Sample GWAS & Mendelian Randomization: Walkthrough

## 1. Study Design Overview

This pipeline tests **causal relationships between plasma proteins and cardiometabolic hallmarks** using a split-sample two-sample Mendelian Randomization (MR) design.

**Hallmark traits in scope:** BMI, HbA1c, TRIG_HDL_RATIO

**Core idea:** The UK Biobank (UKB) sample is split into two non-overlapping subsets:
- **Proteomics sample** (~50K individuals with Olink data) -- used for protein GWAS
- **Held-out sample** (~300K+ remaining individuals) -- used for hallmark trait GWAS

This avoids sample overlap between exposure and outcome GWAS, satisfying a key assumption of two-sample MR. Only **cis-pQTLs** are used as instruments to minimize horizontal pleiotropy risk.

### Directions of MR

| Direction | Exposure (instruments from) | Outcome (effects looked up in) | Purpose |
|-----------|----------------------------|-------------------------------|---------|
| **Forward** | UKB-PPP cis-pQTLs (protein instruments) | UKB hallmark GWAS (held-out sample) | Protein -> Hallmark causal effect |
| **Reverse** | UKB hallmark GWAS (held-out sample instruments) | DECODE Iceland pQTLs | Hallmark -> Protein causal effect |
| **Forward (external replication)** | DECODE cis-pQTLs | UKB hallmark GWAS | Independent replication with zero UKB instrument overlap |

### Population

The GWAS and MR analyses here use the **full cohort only**: ~75-78K individuals for the
hallmark traits (held-out sample) and ~50K for the proteins (Olink sample).

---

## 2. Pipeline Root

```
/n/groups/patel/sivateja/regenie_pipeline/
```

---

## 3. Step-by-Step: Running the Split-Sample GWAS

### 3A. Generate Input Files (EID splitting happens here)

**Script:** `scripts/generate_inputs.R`

```bash
Rscript /n/groups/patel/sivateja/regenie_pipeline/scripts/generate_inputs.R
```

This script:
1. Loads the main UKB phenotype/covariate dataframe:
   `/n/groups/patel/sivateja/UKB/PEWAS_results/data_plus_GLP_complications_glycemic_status_HbA1c_adjusted.RDS`
2. Loads the Olink EID list:
   `/n/groups/patel/sivateja/olink_eids_for_proteins_gwas.RDS`
3. Loads validated proteins:
   `/n/groups/patel/sivateja/UKB/merged_validated_proteins_w_EntrezGeneSymbol_w_protein_var_code_UKB.csv`
4. **Splits** the sample:
   - Proteins: `eid %in% olink_eids_for_proteins_gwas`
   - Hallmarks: `eid %notin% olink_eids_for_proteins_gwas`
5. Writes `pheno.txt` and `covar.txt` to:

```
inputs/
└── full/
    ├── proteins_only/{PROTEIN_NAME}/pheno.txt, covar.txt
    └── hallmarks_heldout/{BMI,HbA1c,TRIG_HDL_RATIO}/pheno.txt, covar.txt
```

### 3B. Run Preflight Checks

```bash
Rscript /n/groups/patel/sivateja/regenie_pipeline/scripts/preflight_checks.R
```

Validates that all files exist, sample sizes look right, and EIDs don't overlap between the two samples.

### 3C. Submit GWAS Jobs (REGENIE)

**Main REGENIE runner:** `scripts/run_regenie_array.sh`
- Runs REGENIE Step 1 (null model on non-imputed genotypes) + Step 2 (association on imputed genotypes)
- Genotype data: `/n/no_backup2/patel/uk_biobank/ukb_genetics/22881/`
- Configuration: `--bsize 1000` (Step 1), `--bsize 400` (Step 2), 8 threads, 50G RAM

**Submit hallmark GWAS (BMI, HbA1c, TRIG_HDL_RATIO in held-out sample):**
```bash
bash /n/groups/patel/sivateja/regenie_pipeline/scripts/submit_hallmarks_full.sh
```

**Submit protein GWAS (Olink sample):**
```bash
bash /n/groups/patel/sivateja/regenie_pipeline/scripts/submit_proteins_full.sh
```

### 3D. Combine Per-Chromosome Results

```bash
bash /n/groups/patel/sivateja/regenie_pipeline/scripts/combine_regenie_results.sh
```

Concatenates chr1-22 into single gzipped files per phenotype/stratum.

---

## 4. GWAS Summary Statistics Locations

### Hallmark trait GWAS (held-out sample)

Located under `results/GWAS/{PHENOTYPE}/`:

| Trait | File |
|-------|------|
| BMI | `results/GWAS/BMI/BMI_full_all_chr.regenie.gz` |
| HbA1c | `results/GWAS/HbA1c/HbA1c_full_all_chr.regenie.gz` |
| TRIG_HDL_RATIO | `results/GWAS/TRIG_HDL_RATIO/TRIG_HDL_RATIO_full_all_chr.regenie.gz` |

**REGENIE format columns:** `CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N TEST BETA SE CHISQ LOG10P EXTRA`

### Protein cis-pQTL instruments (exposure for forward MR)

From the UKB-PPP publication (Sun et al. 2023, Nature):
- **File:** `41586_2023_6592_MOESM3_ESM.xlsx`, sheet "ST10"
- Contains genome-wide significant pQTLs with `cis_trans` column labelling each variant as cis or trans
- **For MR, we restrict to `cis_trans == "cis"` entries only** (variants near the encoding gene)

### DECODE Iceland pQTLs (external, no UKB overlap)

```
/n/groups/patel/IGLOO/DECODE/pQTL/final_somascan_smp/
```
Files named: `Proteomics_SMP_PC0_{seqid}_{GENE}_{ProteinName}_10032022.txt.gz`

These are full GWAS summary statistics from the Icelandic DECODE cohort (~35K), used as the outcome for reverse MR (zero sample overlap with UKB).

---

## 5. Forward MR: Proteins -> Hallmarks (cis-pQTL instruments)

### What it tests
Causal effect of protein levels on hallmark traits (e.g., does LEP causally increase BMI?).

### Instruments
**Cis-pQTLs** from UKB-PPP ST10 -- genome-wide significant SNPs near the protein-encoding gene. Using cis-only instruments reduces the risk of horizontal pleiotropy since the variant acts specifically through the protein.

### Outcome
UKB hallmark GWAS from the **held-out** (non-Olink) sample.

### Key scripts

**Cis-only forward MR (chunked, recommended):**

```bash
# Run cis-pQTL forward MR in 8 parallel chunks
bash /n/groups/patel/sivateja/regenie_pipeline/scripts/submit_cis_trans_chunked.sh
```

This calls `scripts/twosampleMR_cis_trans_chunked.R` with `qtl_type = "cis"`, which explicitly filters ST10 to `cis_trans == "cis"` before running MR.

**Per-hallmark scripts (all pQTLs -- used for initial exploration):**

| Script | Hallmark |
|--------|----------|
| `scripts/twosampleMR_proteins_bmi.R` | BMI |
| `scripts/twosampleMR_proteins_hba1c.R` | HbA1c |
| `scripts/twosampleMR_proteins_trig_hdl_ratio_full.R` | TRIG_HDL_RATIO |

### How the forward MR works (e.g., `twosampleMR_cis_trans_chunked.R` with `cis`)
1. Loads validated proteins from `merged_step1_2_ukb_T2D_coloc.csv`
2. Loads UKB-PPP pQTLs from ST10 and **filters to `cis_trans == "cis"`**
3. Maps validated proteins to their cis-pQTL instruments (by gene symbol via `UKBPPP_ProteinID`)
4. Formats exposure data for `TwoSampleMR::format_data()`
5. Loads the **held-out** hallmark GWAS (e.g., `results/GWAS/BMI/BMI_full_all_chr.regenie.gz`)
6. Harmonizes alleles (`TwoSampleMR::harmonise_data()`)
7. Runs MR methods: Wald ratio (1 SNP), IVW (2+), MR-Egger + Weighted Median (3+)
8. Outputs results CSVs

### Forward MR results

```
results/twosampleMR/
├── cis_trans_chunked/
│   ├── MR_cis_chunk1_of_8.csv  ... MR_cis_chunk8_of_8.csv    # cis-only forward MR
│   └── (trans chunks also exist but not used for main analysis)
├── MR_proteins_BMI_full.csv
├── MR_proteins_HbA1c_full.csv
├── MR_proteins_TRIG_HDL_RATIO_full.csv
├── harmonized_data_proteins_{PHENO}_full.csv
├── heterogeneity_proteins_{PHENO}_full.csv
├── pleiotropy_proteins_{PHENO}_full.csv
└── leaveoneout_proteins_{PHENO}_full.csv
```

---

## 6. Reverse MR: Hallmarks -> Proteins

### What it tests
Causal effect of hallmark traits on protein levels (e.g., does higher BMI causally change LEP levels?). This completes the bidirectional picture.

### Instruments
Genome-wide significant SNPs (P < 5e-8) from UKB hallmark GWAS (held-out sample), distance-pruned at 500kb to ensure independence.

### Outcome
**DECODE Iceland SomaScan pQTLs** -- a completely independent cohort (Icelandic population), ensuring zero sample overlap with UKB.

### Key scripts

| Script | Description |
|--------|-------------|
| `scripts/reverse_mr_bmi_to_all_proteins.R` | BMI -> all proteins via DECODE |
| `scripts/reverse_mr_phenotype_to_all_proteins.R` | Generic, takes phenotype as argument |

**Running reverse MR for each hallmark:**
```bash
Rscript scripts/reverse_mr_bmi_to_all_proteins.R
Rscript scripts/reverse_mr_phenotype_to_all_proteins.R HbA1c
Rscript scripts/reverse_mr_phenotype_to_all_proteins.R TRIG_HDL_RATIO
```

### How it works
1. Loads UKB hallmark GWAS (e.g., `results/GWAS/BMI/BMI_full_all_chr.regenie.gz`)
2. Extracts GW-significant SNPs (LOG10P > 7.3, i.e. P < 5e-8), keeps rsID-named variants
3. Distance-prunes at 500kb (keeping the most significant per window)
4. For **each protein**, loads its DECODE pQTL file from `/n/groups/patel/IGLOO/DECODE/pQTL/final_somascan_smp/`
5. Looks up the hallmark instruments' effects on protein levels in DECODE
6. Harmonizes alleles and runs MR (Wald ratio / IVW / MR-Egger)
7. Reports F-statistics for instrument strength

### Reverse MR results

```
results/twosampleMR/
├── bidirectional_bmi_proteins/
│   ├── bidirectional_MR_BMI_all_proteins.csv
│   ├── bidirectional_MR_BMI_proteins_summary.csv
│   └── protein_decode_mapping.csv          # maps protein names to DECODE files
├── bidirectional_hba1c_proteins/
│   └── bidirectional_MR_HbA1c_all_proteins.csv
└── bidirectional_trig_hdl_ratio_proteins/
    └── bidirectional_MR_TRIG_HDL_RATIO_all_proteins.csv
```

---

## 7. Forward MR with DECODE cis-pQTLs (External Replication)

### What it tests
Forward MR using **DECODE cis-pQTL instruments** (completely independent of UKB) -> UKB hallmark GWAS. This is a fully independent replication of the forward MR that uses zero UKB data for instrument selection.

### Key script
`scripts/bidirectional_mr_decode_trig_hdl.R`

### How it works
1. Loads UKB-PPP ST10 annotations and filters to `cis_trans == "cis"` to define cis regions per protein
2. For each protein, loads the DECODE full GWAS and extracts GW-significant SNPs within the cis window
3. Distance-prunes instruments (500kb)
4. Looks up instrument effects in UKB hallmark GWAS (held-out sample)
5. Runs MR (Wald ratio / IVW / MR-Egger depending on instrument count)

### Results
```
results/twosampleMR/bidirectional_decode_trig_hdl/
└── forward_MR_DECODE_cis_TRIG_HDL.csv
```

---

## 8. Summary of Key Data Flow

```
                    ┌─────────────────────────────┐
                    │   UK Biobank (~500K total)   │
                    └─────────┬───────────────────┘
                              │
              ┌───────────────┴────────────────┐
              │                                │
    ┌─────────▼──────────┐        ┌────────────▼───────────┐
    │  Olink Proteomics  │        │  Held-out (non-Olink)  │
    │    ~50K people     │        │     ~300K+ people      │
    │                    │        │                        │
    │  Protein GWAS      │        │  Hallmark GWAS         │
    │  (REGENIE)         │        │  (BMI, HbA1c,          │
    │                    │        │   TRIG_HDL_RATIO)      │
    └─────────┬──────────┘        └────────────┬───────────┘
              │                                │
              │  UKB-PPP cis-pQTLs (ST10)      │  GW-sig instruments
              │  = protein instruments         │  = hallmark instruments
              ▼                                ▼
    ┌──────────────────┐          ┌──────────────────────────┐
    │  FORWARD MR      │          │  REVERSE MR              │
    │  Protein -> BMI  │          │  BMI -> Protein          │
    │  (cis-pQTL instr │          │  (UKB hallmark instr ->  │
    │   -> held-out    │          │   DECODE pQTLs)          │
    │   hallmark GWAS) │          │                          │
    └──────────────────┘          └──────────────────────────┘
                                            │
                                            │ Uses external cohort:
                                            ▼
                              ┌──────────────────────────┐
                              │  DECODE Iceland pQTLs    │
                              │  (SomaScan, N~35K)       │
                              │  No UKB overlap          │
                              └──────────────────────────┘

    EXTERNAL REPLICATION (Forward):
    DECODE cis-pQTLs (instruments) -> UKB held-out hallmark GWAS (outcome)
    (Zero UKB data used for instrument selection)
```

---

## 9. Quick Reference: Running Everything End-to-End

```bash
cd /n/groups/patel/sivateja/regenie_pipeline

# --- GWAS ---

# 1. Generate split-sample input files
Rscript scripts/generate_inputs.R

# 2. Preflight checks
Rscript scripts/preflight_checks.R

# 3. Submit hallmark GWAS (held-out sample: BMI, HbA1c, TRIG_HDL_RATIO)
bash scripts/submit_hallmarks_full.sh

# 4. Submit protein GWAS (Olink sample)
bash scripts/submit_proteins_full.sh

# 5. Combine per-chromosome results
bash scripts/combine_regenie_results.sh

# --- FORWARD MR (cis-pQTL instruments -> held-out hallmark GWAS) ---

# 6a. Cis-only forward MR (chunked across all hallmarks)
bash scripts/submit_cis_trans_chunked.sh

# 6b. Or per-hallmark scripts
Rscript scripts/twosampleMR_proteins_bmi.R                     # BMI
Rscript scripts/twosampleMR_proteins_hba1c.R                   # HbA1c
Rscript scripts/twosampleMR_proteins_trig_hdl_ratio_full.R     # TRIG_HDL_RATIO

# --- REVERSE MR (hallmark instruments -> DECODE protein outcomes) ---

# 7. Reverse MR
Rscript scripts/reverse_mr_bmi_to_all_proteins.R
Rscript scripts/reverse_mr_phenotype_to_all_proteins.R HbA1c
Rscript scripts/reverse_mr_phenotype_to_all_proteins.R TRIG_HDL_RATIO

# --- EXTERNAL REPLICATION (DECODE cis-pQTLs -> UKB hallmark GWAS) ---

# 8. Forward MR with DECODE cis-pQTL instruments
Rscript scripts/bidirectional_mr_decode_trig_hdl.R

# --- CONCORDANCE ---

# 9. Validation / concordance figures
Rscript scripts/validation_concordance_all_phenotypes.R
```

---

## 10. Key External Data Dependencies

| Resource | Path | Description |
|----------|------|-------------|
| UKB-PPP pQTL supplement | `41586_2023_6592_MOESM3_ESM.xlsx` (sheet ST10) | Sun et al. 2023 pQTL catalog; filter `cis_trans == "cis"` |
| DECODE Iceland pQTLs | `/n/groups/patel/IGLOO/DECODE/pQTL/final_somascan_smp/` | Full GWAS per protein (SomaScan, ~35K Icelanders) |
| UKB genotypes (non-imputed) | `/n/no_backup2/patel/uk_biobank/ukb_genetics/22881/` | REGENIE Step 1 |
| UKB genotypes (imputed) | Same dir, `ukb_imp_chr*_v3.bgen` | REGENIE Step 2 |
| Protein-DECODE mapping | `results/twosampleMR/bidirectional_bmi_proteins/protein_decode_mapping.csv` | Maps validated proteins to DECODE filenames |
| Validated proteins | `/n/groups/patel/sivateja/UKB/merged_validated_proteins_w_EntrezGeneSymbol_w_protein_var_code_UKB.csv` | Proteins passing Step 1 + Step 2 validation |

---

## 11. R Package Requirements

- `TwoSampleMR` -- two-sample MR analysis framework
- `data.table` -- fast file I/O
- `readxl` -- reading UKB-PPP Excel supplement
- `dplyr`, `ggplot2`, `ggrepel`, `patchwork` -- data wrangling and figures

Library path: `/n/groups/patel/sivateja/.R/library`

---

## 12. GitHub-Ready Staging

Clean copies of scripts are staged at:
```
MAP_D_github_staging_OPTION_A/
├── gwas_regenie/scripts/     # REGENIE GWAS scripts
└── mr_twosample/scripts/     # All MR scripts (forward, reverse, DECODE replication)
```
