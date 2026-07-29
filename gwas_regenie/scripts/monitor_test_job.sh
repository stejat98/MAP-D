#!/bin/bash
# Monitor test job and alert when complete
# Usage: bash monitor_test_job.sh [JOB_ID] [check_interval_seconds]

JOB_ID="${1:-24315831}"
CHECK_INTERVAL="${2:-300}"  # Default: check every 5 minutes
EMAIL="sivateja_tangirala@hms.harvard.edu"

echo "=========================================="
echo "Monitoring Test Job: $JOB_ID"
echo "Check interval: $CHECK_INTERVAL seconds"
echo "=========================================="
echo ""

while true; do
    # Check if job is still running
    JOB_STATUS=$(squeue -j $JOB_ID -h -o %T 2>/dev/null)
    
    if [ -z "$JOB_STATUS" ]; then
        # Job finished (not in queue)
        echo "=========================================="
        echo "Job $JOB_ID has finished!"
        echo "=========================================="
        echo ""
        
        # Check if Step 1 completed
        if [ -f /n/scratch/users/s/st320/regenie_test/BMI/diabetes/step1/regenie_step1_BMI_pred.list ]; then
            echo "✓ Step 1 COMPLETE"
            ls -lh /n/scratch/users/s/st320/regenie_test/BMI/diabetes/step1/regenie_step1_BMI_pred.list
        else
            echo "✗ Step 1 prediction file not found - may have failed"
        fi
        
        # Check if Step 2 completed
        if [ -f /n/scratch/users/s/st320/regenie_test/BMI/diabetes/step2/regenie_step2_BMI.regenie.gz ]; then
            echo ""
            echo "✓✓✓ TEST JOB COMPLETE - Both steps finished!"
            ls -lh /n/scratch/users/s/st320/regenie_test/BMI/diabetes/step2/regenie_step2_BMI.regenie.gz
            echo ""
            echo "Output file size: $(du -h /n/scratch/users/s/st320/regenie_test/BMI/diabetes/step2/regenie_step2_BMI.regenie.gz | cut -f1)"
            echo "Number of variants: $(zcat /n/scratch/users/s/st320/regenie_test/BMI/diabetes/step2/regenie_step2_BMI.regenie.gz 2>/dev/null | wc -l)"
            echo ""
            echo "=========================================="
            echo "SUCCESS! Test job completed successfully."
            echo "Ready to submit full pipeline."
            echo "=========================================="
            
            # Send email notification
            echo "Test job $JOB_ID completed successfully. Ready to submit full pipeline." | \
                mail -s "REGENIE Test Job Complete - Ready for Full Pipeline" "$EMAIL" 2>/dev/null || true
            
            exit 0
        else
            echo ""
            echo "✗ Step 2 output file not found - checking logs..."
            echo ""
            echo "Recent output:"
            tail -50 /n/groups/patel/sivateja/regenie_pipeline/slurm/test_bmi_diabetes_${JOB_ID}.out 2>/dev/null | tail -20
            echo ""
            echo "Errors (if any):"
            cat /n/groups/patel/sivateja/regenie_pipeline/slurm/test_bmi_diabetes_${JOB_ID}.err 2>/dev/null || echo "No errors file"
            echo ""
            echo "=========================================="
            echo "WARNING: Test job may have failed"
            echo "Check logs before submitting full pipeline"
            echo "=========================================="
            
            # Send email notification
            echo "Test job $JOB_ID finished but Step 2 output not found. Check logs before submitting full pipeline." | \
                mail -s "REGENIE Test Job Finished - Check Required" "$EMAIL" 2>/dev/null || true
            
            exit 1
        fi
    else
        # Job still running
        RUNTIME=$(squeue -j $JOB_ID -h -o %M 2>/dev/null)
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Job $JOB_ID still running (Runtime: $RUNTIME)"
        
        # Show progress if log exists
        if [ -f /n/scratch/users/s/st320/regenie_test/BMI/diabetes/step1/regenie_step1_BMI.log ]; then
            CURRENT_BLOCK=$(grep "block \[" /n/scratch/users/s/st320/regenie_test/BMI/diabetes/step1/regenie_step1_BMI.log | tail -1 | grep -o "block \[[0-9]*\]" | grep -o "[0-9]*")
            if [ ! -z "$CURRENT_BLOCK" ]; then
                PROGRESS=$(echo "scale=1; $CURRENT_BLOCK * 100 / 454" | bc 2>/dev/null || echo "?")
                echo "  Progress: Block $CURRENT_BLOCK of 454 (~${PROGRESS}%)"
            fi
        fi
        
        echo "  Next check in $CHECK_INTERVAL seconds..."
        echo ""
    fi
    
    sleep $CHECK_INTERVAL
done
