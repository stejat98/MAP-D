# Step-by-Step Test Run Guide

## Prerequisites

1. **Wait for `generate_inputs.R` to complete**
   ```bash
   # Check if job is still running
   squeue -j 24020857
   
   # Check if input files were created
   ls /path/to/project/regenie_pipeline/inputs/full/proteins_only/
   ls /path/to/project/regenie_pipeline/inputs/full/hallmarks_heldout/
   ```

   You should see directories for at least one protein and one hallmark.

---

## Step 1: Check Input Files Exist

```bash
cd /path/to/project

# Check proteins
ls regenie_pipeline/inputs/full/proteins_only/ | head -5

# Check hallmarks  
ls regenie_pipeline/inputs/full/hallmarks_heldout/ | head -5

# Verify a protein has pheno.txt and covar.txt
ls regenie_pipeline/inputs/full/proteins_only/*/pheno.txt | head -1
ls regenie_pipeline/inputs/full/proteins_only/*/covar.txt | head -1
```

**Expected**: You should see at least one protein directory and all 7 hallmark directories (BMI, HbA1c, HDL, LDL, TRIG_HDL_RATIO, systolic_BP, diastolic_BP).

---

## Step 2: Submit Test Job (SLURM)

The test script runs REGENIE, so it needs to run on a compute node. Create a SLURM script:

```bash
cd /path/to/project
cat > regenie_pipeline/scripts/run_test_regenie.sh << 'EOF'
#!/bin/bash
#SBATCH -c 8
#SBATCH -t 2:00:00
#SBATCH --mem=50G
#SBATCH -p short
#SBATCH -o /path/to/project/regenie_pipeline/slurm/test_regenie_%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_email@example.com

cd /path/to/project
bash regenie_pipeline/scripts/test_regenie.sh
EOF

chmod +x regenie_pipeline/scripts/run_test_regenie.sh
```

---

## Step 3: Submit the Test Job

```bash
sbatch regenie_pipeline/scripts/run_test_regenie.sh
```

**Note the job ID** (e.g., `24020858`)

---

## Step 4: Monitor the Test Job

```bash
# Check job status
squeue -j <JOB_ID>

# Watch the log file (replace <JOB_ID> with actual job ID)
tail -f /path/to/project/regenie_pipeline/slurm/test_regenie_<JOB_ID>.out
```

**What to watch for**:
- ✓ "Test phenotype 1 (protein): [protein_name]"
- ✓ "Test phenotype 2 (hallmark): [hallmark_name]"
- ✓ "Running REGENIE step 1..."
- ✓ "Step 1 complete"
- ✓ "Running REGENIE step 2..."
- ✓ "Step 2 complete"

---

## Step 5: Verify Test Results

Once the job completes, check the outputs:

```bash
cd /path/to/project

# Check test output directory structure
ls -R regenie_pipeline/test_results/

# Check Step 1 outputs (prediction files)
ls regenie_pipeline/test_results/step1/proteins_only/*/
ls regenie_pipeline/test_results/step1/hallmarks_heldout/*/

# Check Step 2 outputs (association results)
ls regenie_pipeline/test_results/step2/proteins_only/*/
ls regenie_pipeline/test_results/step2/hallmarks_heldout/*/

# Verify Step 2 output files exist and have content
zcat regenie_pipeline/test_results/step2/proteins_only/*/regenie_step2_*.regenie.gz | head -5
zcat regenie_pipeline/test_results/step2/hallmarks_heldout/*/regenie_step2_*.regenie.gz | head -5
```

**Expected outputs**:
- Step 1: `regenie_step1_*_pred.list` files
- Step 2: `regenie_step2_*.regenie.gz` files (gzipped association results)

---

## Step 6: Check for Errors

```bash
# Check log files for errors
grep -i error regenie_pipeline/slurm/test_regenie_*.out
grep -i "failed\|killed" regenie_pipeline/slurm/test_regenie_*.out

# Check REGENIE logs
grep -i error regenie_pipeline/test_results/step1/*/*/*.log
grep -i error regenie_pipeline/test_results/step2/*/*/*.log
```

**If no errors found**: Test passed! ✅

---

## What the Test Does

1. **Finds first available protein and hallmark** from generated input files
2. **Creates subsets**: Uses first 1000 individuals from each sample
3. **Runs REGENIE Step 1**: Null model on chromosome 22 only (faster)
4. **Runs REGENIE Step 2**: Association testing on chromosome 22 only
5. **Validates outputs**: Checks that prediction files and association results are created

**Test time**: ~30-60 minutes (much faster than full run because it uses only chr22 and 1000 individuals)

---

## Troubleshooting

### Error: "No proteins found"
- **Solution**: Wait for `generate_inputs.R` to complete, or check if it failed

### Error: "REGENIE command not found"
- **Solution**: Make sure conda environment path is correct in `run_regenie_test.sh`

### Error: "File not found" for pgen files
- **Solution**: Verify genotype file paths exist:
  ```bash
  ls /path/to/shared_data/UKB/gwas/ukb_nonimputed_snps.pgen*
  ls /path/to/shared_data/UKB/gwas/UKBallchr.pgen*
  ```

### Test runs but produces no output
- **Solution**: Check REGENIE logs in the test_results directories for detailed error messages

---

## Next Steps After Successful Test

Once the test passes:
1. ✅ Pipeline is working correctly
2. ✅ Input files are formatted correctly
3. ✅ REGENIE can run successfully

Then proceed to submit full GWAS jobs:
```bash
bash regenie_pipeline/scripts/submit_regenie.sh
```
