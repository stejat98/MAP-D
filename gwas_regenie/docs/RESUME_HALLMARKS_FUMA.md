# Resume: Generate FUMA p0.1 files for 5 Hallmarks (Full Stratum)

## Status (Jan 15, 2026)
- Per-chromosome REGENIE results exist for all 5 hallmarks
- Location: `/n/scratch/users/s/st320/regenie/{PHENO}/full/step2/`
- All have 22 chromosome files (complete)

## Phenotypes to Process
1. HDL
2. LDL
3. systolic_BP
4. diastolic_BP
5. TRIG_HDL_RATIO

---

## RUN THESE THREE STEPS (with progress output)

```bash
cd /n/groups/patel/sivateja/regenie_pipeline
```

### Step 1: Concatenate per-chromosome files (~5-10 min)
```bash
./scripts/step1_concatenate_hallmarks.sh
```

### Step 2: Convert to FUMA format (~10-15 min)
```bash
./scripts/step2_convert_to_fuma.sh
```

### Step 3: Filter to P < 0.1 (~5 min)
```bash
./scripts/step3_filter_p0.1.sh
```

---

## Expected Output Files (in results/)
For each phenotype:
- `{PHENO}_full_all_chr.regenie.gz` - Combined REGENIE output
- `{PHENO}_full_all_chr.fuma.gz` - FUMA format  
- `{PHENO}_full_all_chr.p0.1.fuma.gz` - P < 0.1 filtered for FUMA upload

## Notes
- Each script shows progress: [1/5], [2/5], etc.
- Each script shows per-chromosome progress
- BMI and HbA1c already have FUMA files (done previously)
- This adds 5 more phenotypes × 3 files each = 15 new files
