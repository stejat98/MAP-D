# Step-by-Step Plan to Run REGENIE GWAS

## Overview
This plan walks you through running GWAS for validated proteins and hallmark traits using REGENIE, with proper EID splitting for split-sample MR design.

---

## Step 1: Preflight Checks ✅

**Purpose**: Validate all input files and data structure before running.

```bash
cd /n/groups/patel/sivateja
Rscript regenie_pipeline/scripts/preflight_checks.R
```

**What to check**:
- ✓ All input files exist
- ✓ Sample sizes look reasonable (~50k proteomics, ~300k held-out)
- ✓ Validated proteins and hallmarks are available
- ✓ GlycemicStatus distribution looks correct

**If errors**: Fix data issues before proceeding.

---

## Step 2: Generate Input Files 📝

**Purpose**: Create `pheno.txt` and `covar.txt` files for each phenotype with proper EID splitting.

```bash
cd /n/groups/patel/sivateja
Rscript regenie_pipeline/scripts/generate_inputs.R
```

**What this does**:
- Loads main RDS, proteomics EIDs, and validated proteins
- Creates datasets for each stratum (full, prediabetes, diabetes)
- Applies EID split:
  - **Proteins** → `eid %in% olink_eids_for_proteins_gwas`
  - **Hallmarks** → `eid %notin% olink_eids_for_proteins_gwas`
- Writes `pheno.txt` and `covar.txt` for each phenotype

**Expected output**:
- Files in `regenie_pipeline/inputs/{stratum}/{sample_type}/{phenotype}/`
- Each directory contains `pheno.txt` and `covar.txt`

**Time**: ~10-30 minutes depending on data size

---

## Step 3: Test Run 🧪

**Purpose**: Verify the pipeline works correctly on a small subset before full run.

```bash
cd /n/groups/patel/sivateja
bash regenie_pipeline/scripts/test_regenie.sh
```

**What this does**:
- Tests one protein and one hallmark
- Uses first 1000 individuals from each sample
- Runs REGENIE step 1 and step 2 with smaller block sizes
- Validates outputs are created correctly

**Expected output**:
- Test results in `regenie_pipeline/test_results/`
- Step 1 predictions and Step 2 association results

**What to verify**:
- ✓ No errors in logs
- ✓ Step 1 creates `.pred.list` file
- ✓ Step 2 creates `.regenie.gz` file
- ✓ Output files have reasonable size

**If test fails**: Debug issues before running full pipeline.

---

## Step 4: Submit Full GWAS Jobs 🚀

**Purpose**: Run REGENIE GWAS for all phenotypes across all strata.

### Option A: Interactive submission (with confirmation)

```bash
cd /n/groups/patel/sivateja
bash regenie_pipeline/scripts/submit_regenie.sh
```

This will:
1. Generate phenotype list from input files
2. Show you how many jobs will be submitted
3. Ask for confirmation
4. Submit SLURM array jobs

### Option B: Non-interactive submission

```bash
cd /n/groups/patel/sivateja
bash regenie_pipeline/scripts/submit_regenie.sh --yes
```

**What gets submitted**:
- One array job with tasks for each phenotype
- Each task runs REGENIE step 1 + step 2
- Jobs run on `long` partition (7 days max)

**Expected number of jobs**:
- ~(number of validated proteins + 7 hallmarks) × 3 strata × 2 sample types
- Example: ~100 proteins × 3 × 2 = ~600 tasks

---

## Step 5: Monitor Jobs 📊

**Check job status**:
```bash
squeue -u $USER
```

**View logs**:
```bash
# View latest log
tail -f /n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_*.out

# View specific log
tail -f /n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_<JOBID>_<TASKID>.out
```

**Check for errors**:
```bash
grep -i error /n/groups/patel/sivateja/regenie_pipeline/slurm/regenie_*.out
```

---

## Step 6: Verify Results ✅

**Check outputs**:
```bash
# Count completed jobs
ls -1 /n/groups/patel/sivateja/regenie_pipeline/results/step2/*/*/*/*.regenie.gz | wc -l

# Check one result file
zcat /n/groups/patel/sivateja/regenie_pipeline/results/step2/full/proteins_only/<protein_name>/regenie_step2_<protein_name>.regenie.gz | head
```

**Expected outputs**:
- Step 1: `regenie_step1_{phenotype}_pred.list` files
- Step 2: `regenie_step2_{phenotype}.regenie.gz` files (association results)

---

## Troubleshooting

### Issue: Input files not found
**Solution**: Run `generate_inputs.R` first

### Issue: REGENIE fails with memory error
**Solution**: Already using `--lowmem` flag. If still fails, may need to reduce `--bsize`

### Issue: Job times out
**Solution**: Jobs are set to 7 days. If still timing out, may need to split into smaller batches

### Issue: No phenotypes found
**Solution**: Check that `generate_inputs.R` completed successfully and created input files

---

## Expected Timeline

- **Preflight checks**: ~5 minutes
- **Generate inputs**: ~10-30 minutes
- **Test run**: ~30-60 minutes (depends on system)
- **Full GWAS**: ~1-3 days per phenotype (depends on sample size and system load)

---

## Quick Reference

```bash
# 1. Preflight
Rscript regenie_pipeline/scripts/preflight_checks.R

# 2. Generate inputs
Rscript regenie_pipeline/scripts/generate_inputs.R

# 3. Test
bash regenie_pipeline/scripts/test_regenie.sh

# 4. Submit
bash regenie_pipeline/scripts/submit_regenie.sh

# 5. Monitor
squeue -u $USER
tail -f regenie_pipeline/slurm/regenie_*.out
```
