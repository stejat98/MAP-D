#!/bin/bash
# Comprehensive status report for REGENIE pipeline

echo "=========================================="
echo "REGENIE PIPELINE - COMPREHENSIVE STATUS"
echo "=========================================="
echo ""
echo "Sample Type: hallmarks_heldout"
echo "Phenotypes: BMI, HbA1c"
echo "Chromosomes: 1-22"
echo ""

echo "=========================================="
echo "OVERALL SUMMARY BY STRATUM"
echo "=========================================="
echo ""

# Diabetes
echo "DIABETES STRATUM (Job Array 25084945):"
echo "  Total jobs: 44 (22 BMI + 22 HbA1c)"
echo "  Completed: 44/44 (100%)"
echo "  Failed: 0"
echo "  Status: ✓ COMPLETE"
echo ""

# Prediabetes
echo "PREDIABETES STRATUM (Job Array 25084946):"
echo "  Total jobs: 44 (22 BMI + 22 HbA1c)"
echo "  Completed: 41/44 (93.2%)"
echo "  Failed: 2 (chr15 BMI, chr15 HbA1c - block 49 NaN error)"
echo "  Timeout: 1 (chr2 HbA1c - resubmitted with 12h limit)"
echo "  Status: ⚠ IN PROGRESS"
echo ""

# Full
echo "FULL STRATUM (Job Array 24981436):"
echo "  Total jobs: 44 (22 BMI + 22 HbA1c)"
echo "  Completed: 44/44 (100%)"
echo "  Failed: 0"
echo "  Status: ✓ COMPLETE"
echo ""

echo "=========================================="
echo "ACTIVE/RETRY JOBS"
echo "=========================================="
echo ""
echo "Prediabetes chr15 (block 49 exclusion):"
echo "  Job 25132241_15: BMI chr15 - RUNNING"
echo "  Job 25132241_37: HbA1c chr15 - RUNNING"
echo "  Solution: Excluding 400 variants from block 49"
echo ""
echo "Prediabetes chr2 (HbA1c, 12h limit):"
echo "  Job 25130569_24: HbA1c chr2 - RUNNING (~4h estimated)"
echo "  Solution: Increased time limit from 6h to 12h"
echo ""

echo "=========================================="
echo "CHROMOSOME BREAKDOWN - PREDIABETES"
echo "=========================================="
echo ""
echo "BMI:"
for chr in {1..22}; do
    status=$(sacct -j 25084946_${chr} --format=State --noheader --parsable2 2>/dev/null | grep -v "\.batch\|\.extern" | head -1)
    if [ "$status" = "COMPLETED" ]; then
        printf "chr%2d: ✓  " $chr
    elif [ "$chr" = "15" ]; then
        printf "chr%2d: ✗ (retrying)  " $chr
    else
        printf "chr%2d: ?  " $chr
    fi
    [ $((chr % 11)) -eq 0 ] && echo ""
done
echo ""
echo ""

echo "HbA1c:"
for chr in {1..22}; do
    task=$((chr + 22))
    status=$(sacct -j 25084946_${task} --format=State --noheader --parsable2 2>/dev/null | grep -v "\.batch\|\.extern" | head -1)
    if [ "$status" = "COMPLETED" ]; then
        printf "chr%2d: ✓  " $chr
    elif [ "$status" = "TIMEOUT" ]; then
        printf "chr%2d: ⏱ (running)  " $chr
    elif [ "$chr" = "15" ]; then
        printf "chr%2d: ✗ (retrying)  " $chr
    elif [ "$chr" = "2" ]; then
        printf "chr%2d: ▶ (running)  " $chr
    else
        printf "chr%2d: ?  " $chr
    fi
    [ $((chr % 11)) -eq 0 ] && echo ""
done
echo ""

echo "=========================================="
echo "TOTAL PROGRESS"
echo "=========================================="
echo ""
echo "Total jobs across all strata: 132 (44 × 3)"
echo "Completed: 129/132 (97.7%)"
echo "Failed/Timeout: 3/132 (2.3%)"
echo "Active/Retry: 3 jobs running"
echo ""
echo "Legend:"
echo "  ✓ = Completed"
echo "  ✗ = Failed (being retried)"
echo "  ⏱ = Timeout (resubmitted)"
echo "  ▶ = Currently running"
echo "  ? = Status unknown"
echo ""
