# Monitor Generate Inputs Job

## Job ID: 24020857

### Quick Status Check
```bash
# Check if job is running/completed
squeue -j 24020857

# Or use sacct for more details
sacct -j 24020857 --format=JobID,JobName,State,ExitCode,Elapsed,MaxRSS,NodeList
```

### Watch Log File (once job starts)
```bash
# Watch the SLURM output log
tail -f /n/groups/patel/sivateja/regenie_pipeline/slurm/generate_inputs_24020857.out

# Or watch the main log file
tail -f /n/groups/patel/sivateja/regenie_pipeline/generate_inputs.log
```

### Check Progress
```bash
# Count how many input files have been created
find /n/groups/patel/sivateja/regenie_pipeline/inputs -name "pheno.txt" | wc -l

# Check specific strata
ls /n/groups/patel/sivateja/regenie_pipeline/inputs/full/proteins_only/ | wc -l
ls /n/groups/patel/sivateja/regenie_pipeline/inputs/full/hallmarks_heldout/
```

### Expected Outputs
- **Proteins**: Should create directories for each available protein (likely 0 if names don't match)
- **Hallmarks**: Should create 7 directories (BMI, HbA1c, HDL, LDL, TRIG_HDL_RATIO, systolic_BP, diastolic_BP)
- **Strata**: full, prediabetes, diabetes

### When Job Completes
1. Check exit code: `sacct -j 24020857 --format=ExitCode`
2. Verify files: `find regenie_pipeline/inputs -name "pheno.txt" | wc -l`
3. Check for errors: `grep -i error regenie_pipeline/slurm/generate_inputs_24020857.out`

### Next Steps After Completion
Once job completes successfully:
1. Run test: `sbatch regenie_pipeline/scripts/run_test_regenie.sh`
2. Then submit full GWAS: `bash regenie_pipeline/scripts/submit_regenie.sh`
