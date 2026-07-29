# Chromosome 15 Block 49 Investigation Summary

## Error Details

**Failed Jobs:**
- Task 15: BMI, Chromosome 15 (Job 25084946_15)
- Task 37: HbA1c, Chromosome 15 (Job 25084946_37)

**Error Message:**
```
block [49/6920] : ERROR: Throw location unknown (consider using BOOST_THROW_EXCEPTION)
Dynamic exception type: boost::wrapexcept<std::domain_error>
std::exception::what: Error in function boost::math::cdf(const chi_squared_distribution<double>&, double): Chi Square parameter was -nan, but must be > 0 !
```

## Block 49 Details

- **Block number:** 49 out of 6,920 total blocks
- **Block size:** 400 variants per block
- **Variant range:** Variants 19,201 to 19,600
- **Genomic region:** Chromosome 15, positions ~20,858,089 to ~20,867,998
- **Total variants in chr15:** 2,767,971

## Key Finding

**Diabetes and Full strata successfully completed block 49:**
- ✓ diabetes BMI chr15: PASSED block 49
- ✓ full BMI chr15: PASSED block 49  
- ✓ diabetes HbA1c chr15: PASSED block 49
- ✓ full HbA1c chr15: PASSED block 49

This indicates the issue is **stratum-specific to prediabetes**, not a general data quality problem.

## Possible Causes

1. **Smaller sample size in prediabetes:** May lead to numerical instability
2. **Different phenotype distribution:** Could affect test statistic calculations
3. **Covariate patterns:** Different covariate distributions in prediabetes samples
4. **Variant-sample interaction:** Specific variants may have issues only in prediabetes samples

## Variants in Block 49

Block 49 contains variants in the region ~20.8-20.9 Mb on chromosome 15, including:
- Standard SNPs (rs IDs)
- Indels (e.g., 15:20858209_TA_T, 15:20858252_AT_A)

## Action Plan

1. **Resubmit failed jobs** - May succeed if transient issue
2. **Monitor closely** - If fails again, investigate prediabetes-specific factors
3. **If persistent:** Consider:
   - Checking variant MAF/MAC in prediabetes samples specifically
   - Verifying phenotype/covariate distributions
   - Using alternative statistical approaches for problematic variants

## Resubmission

Use script: `resubmit_prediabetes_failed.sh`
- Tasks 15, 37: Resubmit with standard parameters
- Task 24: Resubmit with 12-hour time limit (was timeout)
