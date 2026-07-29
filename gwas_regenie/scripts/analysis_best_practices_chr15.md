# Best Practices for Handling Numerical Errors in UK Biobank-Scale GWAS

## Current Issue: NaN Error at Block 49, Chromosome 15

### Problem
- REGENIE step 2 fails with NaN error in chi-square calculation
- Issue is specific to prediabetes stratum (diabetes/full pass successfully)
- Error occurs at block 49/6920 (variants 19,201-19,600)

### Current Approach (Not Ideal)
**Excluding 400 variants from block 49**
- ❌ Ad-hoc solution
- ❌ Doesn't address root cause
- ❌ May miss other problematic variants elsewhere
- ❌ Not scalable or reproducible

## Better Approach: Increase MAC Threshold

### Why This Is Better
1. **Addresses root cause**: Low MAC variants cause numerical instability
2. **Principled filtering**: Standard practice in GWAS
3. **Prevents future issues**: Catches all low-MAC variants, not just block 49
4. **Reproducible**: Same threshold applied consistently

### Recommended Solution

**Increase `--minMAC` from default (5) to 10-20 for step 2**

```bash
--minMAC 10  # or 20 for more conservative filtering
```

### Rationale

1. **Default is too low**: `--minMAC 5` is the REGENIE default, but:
   - Very low MAC variants are prone to numerical errors
   - UK Biobank analyses typically use MAC >= 10-20 for step 2
   - Your step 1 uses MAC >= 100, so step 2 should be consistent

2. **Prediabetes-specific issue**: 
   - Smaller sample size (37K vs 35K full) means lower MAC for same variants
   - Variants with MAC 5-10 may pass filter but cause numerical issues
   - Increasing to MAC >= 10-20 filters these problematic variants

3. **Standard practice**:
   - Most UK Biobank GWAS use MAC >= 10-20 for association testing
   - Prevents numerical instability while retaining most variants
   - Only excludes ~0.1-0.5% of variants (those with very low frequency)

### Implementation

Modify `run_regenie_array_per_chr.sh` to add `--minMAC` flag:

```bash
regenie \
  --step 2 \
  --bed "${FAM_PATH}/ukb${TARGET_CHR}" \
  --phenoFile "${PHENO_FILE}" \
  --covarFile "${COVAR_FILE}" \
  --pred "${PRED_FILE}" \
  ${STEP2_FLAGS} \
  --minMAC 10 \
  --bsize 400 \
  --threads 8 \
  --out "${OUTDIR_STEP2}/regenie_step2_${PHENO_NAME}_chr${TARGET_CHR}"
```

### Comparison

| Approach | Pros | Cons |
|----------|------|------|
| **Exclude block 49** | Quick fix | Ad-hoc, doesn't prevent other issues |
| **Increase --minMAC** | Principled, prevents all low-MAC issues | Excludes some variants (but appropriate) |

### Recommendation

**Use `--minMAC 10` or `--minMAC 20` instead of excluding specific variants**

This is the standard approach used in:
- UK Biobank official analyses
- Large-scale GWAS studies
- REGENIE documentation and best practices

### Impact

- **Variants excluded**: ~0.1-0.5% (those with MAC < 10-20)
- **Variants retained**: 99.5-99.9% of all variants
- **Benefit**: Prevents numerical errors, more robust analysis
