# REGENIE Test GWAS Resume Guide

## Quick Start

To resume the regenie test GWAS pipeline, simply run:

```bash
bash /path/to/project/regenie_pipeline/scripts/resume_regenie.sh
```

Or to skip confirmation prompt:

```bash
bash /path/to/project/regenie_pipeline/scripts/resume_regenie.sh --yes
```

## What the Resume Script Does

1. **Checks completed jobs**: Scans the output directory (`/path/to/scratch/regenie_test/`) to identify which phenotypes have completed Step 2 (final GWAS output)

2. **Identifies remaining work**: Creates a list of phenotypes that are:
   - Not started
   - Partially completed (Step 1 done but Step 2 missing)
   - Failed (empty or missing output files)

3. **Submits array jobs**: Submits SLURM array jobs for all remaining phenotypes

4. **Smart resuming**: The `run_regenie_array.sh` script automatically skips Step 1 if it's already completed, saving time and compute resources

## Output Structure

Outputs follow the scratch directory convention:
```
/path/to/scratch/regenie_test/
  <PHENOTYPE>/
    <STRATUM>/
      step1/          # Step 1 outputs (null model)
      step2/          # Step 2 outputs (association testing)
```

## Monitoring Jobs

Check job status:
```bash
squeue -u $USER
```

Check logs:
```bash
ls -lh /path/to/project/regenie_pipeline/slurm/
```

## Conda Environment Setup

All scripts automatically set the required conda environment variables:
- `HOME=/path/to/project`
- `CONDA_NO_PLUGINS=true`
- `CONDA_PKGS_DIRS=/path/to/project/.conda/pkgs`
- `CONDA_ENVS_PATH=/path/to/project/.conda/envs`

And activate the regenie environment:
- `/path/to/project/.conda/envs/regenie_env`

## Files Created

- `pheno_list.txt`: Full list of all phenotypes to process
- `pheno_list_remaining.txt`: List of phenotypes that still need processing
- `pheno_list_remaining_<JOB_ID>.txt`: Temporary list used by array jobs

## Troubleshooting

If a job fails, check the log files:
```bash
tail -f /path/to/project/regenie_pipeline/slurm/regenie_<JOB_ID>_<ARRAY_ID>.out
tail -f /path/to/project/regenie_pipeline/slurm/regenie_<JOB_ID>_<ARRAY_ID>.err
```

To manually check completion status:
```bash
# Check if a specific phenotype is complete
ls -lh /path/to/scratch/regenie_test/<PHENOTYPE>/<STRATUM>/step2/regenie_step2_<PHENOTYPE>.regenie.gz
```
