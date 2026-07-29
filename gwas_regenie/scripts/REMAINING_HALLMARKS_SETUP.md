# Remaining Hallmarks Setup - Chromosome 15 Fix Applied

## Overview
This document describes the setup for running remaining hallmarks (HDL, LDL, TRIG_HDL_RATIO, systolic_BP, diastolic_BP) using the chromosome-specific configuration that worked successfully for BMI and HbA1c.

## Chromosome 15 Fix (Prediabetes Stratum)

### Issue Identified
- **Chromosome 15** had NaN errors at **block 49** specifically for **prediabetes stratum**
- Error: "Chi Square parameter was -nan" at block 49/6920
- Affected: BMI and HbA1c in prediabetes stratum
- **Diabetes and Full strata completed successfully** - issue is prediabetes-specific

### Solution Applied
The fix is **already built into** `run_regenie_array_per_chr.sh`:
- Automatically excludes block 49 variants (variants 19,201-19,600) for chr15 prediabetes
- Check: `if [ "$TARGET_CHR" = "15" ] && [ "$STRATUM" = "prediabetes" ]`
- Excludes ~400 variants from block 49 to prevent NaN errors

### Affected Jobs for Remaining Hallmarks
The following prediabetes chr15 jobs will automatically use the block 49 exclusion:
- Task 15: HDL chr15
- Task 37: LDL chr15  
- Task 59: TRIG_HDL_RATIO chr15
- Task 81: systolic_BP chr15
- Task 103: diastolic_BP chr15

**No manual intervention needed** - the run script handles this automatically.

## Files Created

### Phenotype Lists (110 lines each = 5 phenotypes × 22 chromosomes)
- `pheno_list_diabetes_per_chr_remaining.txt`
- `pheno_list_prediabetes_per_chr_remaining.txt` (includes chr15 with automatic fix)
- `pheno_list_full_per_chr_remaining.txt`

### Submit Scripts
- `submit_diabetes_per_chr_remaining.sh` - short partition, 6 hours
- `submit_prediabetes_per_chr_remaining.sh` - short partition, 6 hours (chr15 auto-fixed)
- `submit_full_per_chr_remaining.sh` - medium partition, 2 days
- `submit_remaining_hallmarks_all_strata.sh` - master script to submit all

## Configuration Details

### Format
- `stratum|hallmarks_heldout|phenotype|chromosome`
- Example: `prediabetes|hallmarks_heldout|HDL|15`

### Partitions & Time Limits
- **Diabetes/Prediabetes**: `short` partition, 6 hours (same as BMI/HbA1c)
- **Full**: `medium` partition, 2 days (same as BMI/HbA1c)

### Run Script
- Uses: `run_regenie_array_per_chr.sh`
- Automatically applies chr15 block 49 exclusion for prediabetes
- No modifications needed

## Usage

```bash
# Submit all strata (recommended)
cd /n/groups/patel/sivateja/regenie_pipeline/scripts
bash submit_remaining_hallmarks_all_strata.sh

# Or submit individually
bash submit_diabetes_per_chr_remaining.sh
bash submit_prediabetes_per_chr_remaining.sh  # chr15 auto-fixed
bash submit_full_per_chr_remaining.sh
```

## Monitoring

```bash
# Check job status
squeue -u st320

# Check prediabetes chr15 jobs (should show block 49 exclusion in logs)
squeue -u st320 -p short | grep regenie_prediabetes_perchr_remaining

# View logs for chr15 prediabetes jobs
tail -f slurm/regenie_prediabetes_perchr_remaining_*_15.out
tail -f slurm/regenie_prediabetes_perchr_remaining_*_37.out
tail -f slurm/regenie_prediabetes_perchr_remaining_*_59.out
tail -f slurm/regenie_prediabetes_perchr_remaining_*_81.out
tail -f slurm/regenie_prediabetes_perchr_remaining_*_103.out
```

## Expected Behavior

### For Prediabetes Chr15 Jobs
When running, you should see in the logs:
```
[PHENO_NAME] Creating exclude file for block 49 variants (fixing NaN error)...
[PHENO_NAME] Note: MAC filtering (5, 10, 100) all failed - excluding block 49 variants
[PHENO_NAME] Excluding 400 variants from block 49 (variants 19201-19600)
```

This confirms the fix is being applied automatically.

## Total Jobs
- **Diabetes**: 110 jobs (5 phenotypes × 22 chromosomes)
- **Prediabetes**: 110 jobs (5 phenotypes × 22 chromosomes) - chr15 auto-fixed
- **Full**: 110 jobs (5 phenotypes × 22 chromosomes)
- **Total**: 330 jobs

## Notes
- The chr15 block 49 exclusion is **automatic** - no manual intervention needed
- Only affects prediabetes stratum (diabetes and full completed successfully)
- The exclusion removes ~400 variants (0.014% of chr15) to prevent NaN errors
- This is the same fix that successfully resolved BMI and HbA1c chr15 issues
