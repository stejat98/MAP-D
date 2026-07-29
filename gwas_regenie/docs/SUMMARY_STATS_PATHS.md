# REGENIE GWAS Summary Statistics - BMI and HbA1c

## Combined Summary Statistics Files (Recommended)

The per-chromosome results have been combined into single files following UK Biobank GWAS best practices (Neale lab style). These are the **recommended files to use** for downstream analyses.

### Combined Files (All Chromosomes)

**Location:** `/n/groups/patel/sivateja/regenie_pipeline/results/`

| Phenotype | Stratum | File Path | Variants | Size | Sample Size (N) |
|-----------|---------|-----------|----------|------|----------------|
| BMI | diabetes | `BMI_diabetes_all_chr.regenie.gz` | 20,532,508 | 581M | ~2,009 |
| BMI | prediabetes | `BMI_prediabetes_all_chr.regenie.gz` | 32,032,589 | 927M | ~8,311 |
| BMI | full | `BMI_full_all_chr.regenie.gz` | 48,835,193 | 1.5G | ~78,466 |
| HbA1c | diabetes | `HbA1c_diabetes_all_chr.regenie.gz` | 20,282,727 | 573M | ~1,935 |
| HbA1c | prediabetes | `HbA1c_prediabetes_all_chr.regenie.gz` | 32,070,846 | 933M | ~8,336 |
| HbA1c | full | `HbA1c_full_all_chr.regenie.gz` | 48,517,826 | 1.5G | ~74,986 |

**File Format:**
- **Format:** REGENIE output (space-separated, gzipped)
- **Columns:** `CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N TEST BETA SE CHISQ LOG10P EXTRA`
- **Sorted by:** CHROM (numeric), GENPOS (numeric)
- **Compression:** gzip

**Usage:**
```bash
# View first 20 lines
zcat /n/groups/patel/sivateja/regenie_pipeline/results/BMI_diabetes_all_chr.regenie.gz | head -20

# Search for specific variant
zcat /n/groups/patel/sivateja/regenie_pipeline/results/BMI_diabetes_all_chr.regenie.gz | grep -w 'rs12345'

# Count variants
zcat /n/groups/patel/sivateja/regenie_pipeline/results/BMI_diabetes_all_chr.regenie.gz | tail -n +2 | wc -l
```

---

## Per-Chromosome Files (Raw Output)

If you need access to individual chromosome files, they are located in the scratch directory:

**Location:** `/n/scratch/users/s/st320/regenie/`

### BMI - Diabetes
- **Step 1 (null model):** `/n/scratch/users/s/st320/regenie/BMI/diabetes/step1/`
- **Step 2 (associations):** `/n/scratch/users/s/st320/regenie/BMI/diabetes/step2/`
  - Files: `regenie_step2_BMI_chr{1..22}_BMI.regenie` (22 files, ~23-141M each)

### BMI - Prediabetes
- **Step 1 (null model):** `/n/scratch/users/s/st320/regenie/BMI/prediabetes/step1/`
- **Step 2 (associations):** `/n/scratch/users/s/st320/regenie/BMI/prediabetes/step2/`
  - Files: `regenie_step2_BMI_chr{1..22}_BMI.regenie` (22 files, ~37-224M each)

### HbA1c - Diabetes
- **Step 1 (null model):** `/n/scratch/users/s/st320/regenie/HbA1c/diabetes/step1/`
- **Step 2 (associations):** `/n/scratch/users/s/st320/regenie/HbA1c/diabetes/step2/`
  - Files: `regenie_step2_HbA1c_chr{1..22}_HbA1c.regenie` (22 files)

### HbA1c - Prediabetes
- **Step 1 (null model):** `/n/scratch/users/s/st320/regenie/HbA1c/prediabetes/step1/`
- **Step 2 (associations):** `/n/scratch/users/s/st320/regenie/HbA1c/prediabetes/step2/`
  - Files: `regenie_step2_HbA1c_chr{1..22}_HbA1c.regenie` (22 files)

### HbA1c - Full
- **Step 1 (null model):** `/n/scratch/users/s/st320/regenie/HbA1c/full/step1/`
- **Step 2 (associations):** `/n/scratch/users/s/st320/regenie/HbA1c/full/step2/`
  - Files: `regenie_step2_HbA1c_chr{1..22}_HbA1c.regenie` (22 files)

### BMI - Full
- **Step 1 (null model):** `/n/scratch/users/s/st320/regenie/BMI/full/step1/`
- **Step 2 (associations):** `/n/scratch/users/s/st320/regenie/BMI/full/step2/`
  - Files: `regenie_step2_BMI_chr{1..22}_BMI.regenie` (22 files)

**Note:** Scratch directory files may be purged after a period of inactivity. The combined files in `/n/groups/patel/sivateja/regenie_pipeline/results/` are the permanent location.

---

## Processing Steps

### What Was Done

1. **GWAS Analysis:** REGENIE Step 1 (null model) and Step 2 (association testing) were run per chromosome
2. **Combination:** Per-chromosome results were combined into single files:
   - Header from chromosome 1
   - Data from all 22 chromosomes (skipping headers)
   - Sorted by chromosome and genomic position
   - Compressed with gzip

### Script Used

The combination was performed using:
```bash
/n/groups/patel/sivateja/regenie_pipeline/scripts/combine_regenie_results.sh
```

To re-run combination (if needed):
```bash
# Combine all BMI and HbA1c results
/n/groups/patel/sivateja/regenie_pipeline/scripts/combine_regenie_results.sh

# Combine specific phenotype/stratum
/n/groups/patel/sivateja/regenie_pipeline/scripts/combine_regenie_results.sh BMI diabetes
```

---

## Best Practices (UK Biobank / Neale Lab Style)

The combined files follow standard UK Biobank GWAS practices:

1. **Single file per analysis:** All chromosomes combined into one file
2. **Sorted:** By chromosome (numeric) and genomic position (numeric)
3. **Compressed:** gzip compression for storage efficiency
4. **Standard format:** REGENIE output format with standard columns
5. **Permanent location:** Stored in `/n/groups/patel/sivateja/regenie_pipeline/results/` (not scratch)

This matches the approach used by:
- Neale Lab UK Biobank GWAS (https://www.nealelab.is/uk-biobank)
- Standard REGENIE output conventions
- UK Biobank GWAS best practices

---

## Column Descriptions

| Column | Description |
|--------|-------------|
| CHROM | Chromosome number |
| GENPOS | Genomic position (bp) |
| ID | Variant ID (rsID or chr:pos:ref:alt) |
| ALLELE0 | Reference allele |
| ALLELE1 | Alternate allele |
| A1FREQ | Frequency of alternate allele |
| N | Sample size |
| TEST | Test type (ADD for additive) |
| BETA | Effect size (beta coefficient) |
| SE | Standard error of beta |
| CHISQ | Chi-square statistic |
| LOG10P | -log10(p-value) |
| EXTRA | Additional information (usually NA) |

---

## Verification & Safety Checks

The combination script includes multiple safety checks to ensure no mix-ups:

1. **File naming verification:** Each file must contain the correct phenotype name (BMI or HbA1c) in the filename
2. **Chromosome verification:** Each chromosome file is checked to ensure it contains data for the correct chromosome
3. **Sample size verification:** Sample sizes are consistent within each stratum:
   - **Full stratum:** ~75-78K individuals (largest)
   - **Prediabetes stratum:** ~8-8.3K individuals (medium)
   - **Diabetes stratum:** ~1.9-2K individuals (smallest)
4. **Chromosome range:** All files span chromosomes 1-22 (autosomes only)
5. **Sorting:** Files are sorted by chromosome (numeric) and genomic position (numeric)

**Verified sample sizes (from first variant in each file):**
- BMI_diabetes: N=2,009
- BMI_prediabetes: N=8,311
- BMI_full: N=78,466
- HbA1c_diabetes: N=1,935
- HbA1c_prediabetes: N=8,336
- HbA1c_full: N=74,986

## Notes

- **Sample sizes:** 
  - Full stratum: ~75-78K individuals
  - Prediabetes stratum: ~8-8.3K individuals
  - Diabetes stratum: ~1.9-2K individuals
- **All 22 chromosomes:** Complete results for autosomes (chr1-22)
- **MAC filtering:** Variants with MAC < 5 were excluded (REGENIE default)
- **Format:** Space-separated, gzipped text files
- **Compatibility:** Standard REGENIE format, compatible with most GWAS tools (PLINK, FUMA, LocusZoom, etc.)
- **File integrity:** All combined files verified to start at CHROM=1 and end at CHROM=22
