# REGENIE GWAS Pipeline - Status Update

## Overview
We're running genome-wide association studies (GWAS) using REGENIE v4.1.2 for **hallmark traits** and **proteins** across three glycemic strata (full cohort, prediabetes, diabetes). This uses a two-step approach: Step 1 fits a null model to account for population structure, and Step 2 performs association testing across all chromosomes.

## Analytical Approach
- **Method**: REGENIE (fast, scalable GWAS method)
- **Genotype data**: UK Biobank imputed data (22 chromosomes, ~1.2M variants/chr)
- **Sample split**: 
  - Hallmarks analyzed in held-out sample (n≈418K)
  - Proteins analyzed in proteomics-only sample (n≈55K)
- **Covariates**: Age, sex, genetic principal components (44 covariates)
- **Strata**: Full cohort, prediabetes subset, diabetes subset

## Phenotypes Analyzed
**Hallmarks (7 traits):** BMI, HbA1c, HDL, LDL, TRIG_HDL_RATIO, systolic_BP, diastolic_BP  
**Proteins (405 proteins):** Validated proteins from Olink panel

## Current Status

### ✅ Completed
- Fixed IID matching issues between phenotype and genotype files
- Updated pipeline to use correct genotype files (`/n/no_backup2/patel/uk_biobank/ukb_genetics/22881/`)
- Verified code works (test run successfully processed 349K individuals)
- All input files generated and validated

### 🚀 Running Now
**Hallmarks (21 jobs):**
- 7 phenotypes × 3 strata
- Estimated completion: 1-2.5 days per job
- Status: Queued and running

**Proteins (1,215 jobs):**
- 405 proteins × 3 strata  
- Optimized for fairshare (26 batches, sequential execution)
- Estimated completion: 4-8 hours per job (much faster due to smaller sample size)
- Status: Queued with dependencies (batched to manage cluster resources)

### 📊 Progress Tracking
- Monitor jobs: `squeue -u st320 | grep regenie`
- Logs: `/n/groups/patel/sivateja/regenie_pipeline/slurm/`
- Outputs: `/n/scratch/users/s/st320/regenie/`

## Expected Timeline
- **Hallmarks**: ~1-2 weeks (21 jobs, running in parallel)
- **Proteins**: ~2-3 weeks (1,215 jobs, batched sequentially for fairshare)
- **Total completion**: ~3-4 weeks for all analyses

## Outputs
Each analysis produces:
- Step 1: Null model predictions (used for Step 2)
- Step 2: Association results per chromosome (`.regenie.gz` files)
- Results include: variant ID, position, allele info, effect sizes, p-values, standard errors

## Next Steps
1. Monitor job completion
2. Combine chromosome-level results
3. Quality control and filtering
4. Prepare results for downstream analysis

---

*Last updated: Dec 22, 2025*
