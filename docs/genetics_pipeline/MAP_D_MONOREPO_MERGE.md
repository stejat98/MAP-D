# Merging this staging folder into [MAP-D](https://github.com/stejat98/MAP-D)

This directory was generated from `regenie_pipeline` as **Option A**: add three top-level folders to the existing MAP-D repo.

**Regenerate this staging folder:** from `regenie_pipeline`, run `module load gcc/14.2.0 R/4.4.2 && Rscript build_github_staging.R`.

## Layout after merge

Place contents at the **root** of `MAP-D` (alongside existing `src/`, `main_analysis.R`, etc.):

| Staged folder | Target on GitHub |
|---------------|------------------|
| `gwas_regenie/` | `gwas_regenie/` |
| `mr_twosample/` | `mr_twosample/` |
| `figures_tables/` | `figures_tables/` |

## Steps

1. Clone MAP-D: `git clone https://github.com/stejat98/MAP-D.git && cd MAP-D`
2. Copy (or merge) each staged folder into the repo root.
3. Append `SUGGESTED_ROOT_GITIGNORE_APPEND.txt` to `.gitignore` (edit paths if needed).
4. Refactor hardcoded cluster paths in R/Python scripts to `config.R` or environment variables before pushing.
5. Add `renv` / `sessionInfo()` and document REGENIE + TwoSampleMR versions in README.

## Not included (by design)

- `results/`, `inputs/`, `test_results/`, `filtered_snps/`, `py_pkgs/`
- UK Biobank individual-level data
- SLURM `.out` / `.err` logs
- Analyses outside the BMI / HbA1c / TRIG_HDL_RATIO scope
- LDSC tooling: use the official [LDSC](https://github.com/bulik/ldsc) repository with `regenie_to_ldsc_munge_input.py`

