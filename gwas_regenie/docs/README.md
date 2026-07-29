# REGENIE GWAS Pipeline

This pipeline runs GWAS using REGENIE for validated plasma proteins and the cardiometabolic hallmark traits in the full cohort, with strict sample splitting for the split-sample MR design.

## Overview

- **Proteins GWAS**: Run in proteomics-only sample (`eid %in% olink_eids_for_proteins_gwas`)
- **Hallmarks GWAS**: Run in held-out sample (`eid %notin% olink_eids_for_proteins_gwas`)
- **Cohort**: Full population only

## Directory Structure

```
regenie_pipeline/
├── inputs/
│   └── full/
│       ├── proteins_only/      # Protein phenotypes (proteomics sample)
│       └── hallmarks_heldout/  # Hallmark traits (held-out sample)
├── scripts/
│   ├── generate_inputs.R        # Generate pheno.txt and covar.txt files
│   ├── run_regenie_array.sh     # REGENIE step 1 + step 2 (array job)
│   ├── submit_regenie.sh        # Submit array jobs
│   └── preflight_checks.R       # Validation and sanity checks
├── slurm/                        # SLURM log files
└── results/
    ├── step1/                    # Step 1 outputs (null models)
    └── step2/                    # Step 2 outputs (association results)
```

## Input Files

1. **Main phenotype/covariate dataframe**:
   `/path/to/project/UKB/PEWAS_results/data_plus_GLP_complications_glycemic_status_HbA1c_adjusted.RDS`

2. **Proteomics EID list**:
   `/path/to/project/olink_eids_for_proteins_gwas.RDS`

3. **Validated proteins list**:
   `/path/to/project/UKB/merged_validated_proteins_2.csv`

4. **Genotype data**:
   - Step 1: `/path/to/shared_data/UKB/gwas/ukb_nonimputed_snps`
   - Step 2: `/path/to/shared_data/UKB/gwas/UKBallchr`

## Hallmark Traits

- BMI
- HbA1c
- TRIG_HDL_RATIO (triglyceride / HDL ratio)

## Usage

### Step 1: Preflight Checks

```bash
Rscript /path/to/project/regenie_pipeline/scripts/preflight_checks.R
```

This validates input files, checks data structure, and reports sample sizes.

### Step 2: Generate Input Files

```bash
Rscript /path/to/project/regenie_pipeline/scripts/generate_inputs.R
```

This script:
- Loads the main dataframe, proteomics EIDs, and validated proteins
- Creates stratum-specific datasets
- Applies EID split rules (proteins in proteomics sample, hallmarks in held-out)
- Writes `pheno.txt` and `covar.txt` for each phenotype

### Step 3: Submit GWAS Jobs

```bash
/path/to/project/regenie_pipeline/scripts/submit_regenie.sh
```

This script:
- Generates a phenotype list from generated input files
- Submits SLURM array jobs for all phenotypes
- Each job runs REGENIE step 1 (null model) and step 2 (association testing)

## REGENIE Configuration

- **Step 1**: Uses `ukb_nonimputed_snps` with `--bsize 1000`, `--lowmem`
- **Step 2**: Uses `UKBallchr` with `--bsize 400`
- **Threads**: 8
- **Memory**: 50G
- **Time**: 7 days (long partition)

## Output Files

### Step 1 Outputs
- `regenie_step1_{phenotype}_pred.list`: Prediction file for step 2
- `regenie_step1_{phenotype}.log`: Log file

### Step 2 Outputs
- `regenie_step2_{phenotype}.regenie.gz`: Association results (gzipped)
- `regenie_step2_{phenotype}.log`: Log file

## Sample Split Enforcement

The pipeline strictly enforces the EID split:

- **Proteins**: Only individuals where `eid %in% olink_eids_for_proteins_gwas`
- **Hallmarks**: Only individuals where `eid %notin% olink_eids_for_proteins_gwas`

This ensures no overlap between protein and hallmark GWAS samples, supporting split-sample MR design.

## Monitoring Jobs

```bash
# Check job status
squeue -u $USER

# View logs
tail -f /path/to/project/regenie_pipeline/slurm/regenie_*.out
```

## Notes

- All work is performed within `/path/to/project/`
- Reference implementation (read-only): `/path/to/collab_data/UK_Biobank/BScripts/GWAS/REGENIE`
- The pipeline handles continuous traits by default; binary traits require `--bt` flag modification
