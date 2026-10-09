#!/bin/bash
set -euo pipefail

# Example script to run evolvdx/evolvdx-methylseq on AWS Batch
# Prerequisites:
#   - AWS CLI configured with appropriate credentials
#   - nextflow-env conda environment activated: micromamba activate nextflow-env
#   - AWS Batch compute environment and job queue set up

# Parameters
PARAMS_FILE="conf/aws/aws.params.json"
PROFILE="aws_500gb"  # Use 'aws' for standard queue, 'aws_500gb' for high-memory runs
WORK_DIR="s3://your-bucket/work/nextflow/$(date +%Y%m%d_%H%M%S)"

echo "Starting evolvdx/evolvdx-methylseq run on AWS Batch"
echo "Profile: $PROFILE"
echo "Params file: $PARAMS_FILE"
echo "Work directory: $WORK_DIR"

nextflow run evolvdx/evolvdx-methylseq \
    -profile "$PROFILE" \
    -params-file "$PARAMS_FILE" \
    -work-dir "$WORK_DIR" \
    -resume

echo "Run complete! Results available at $(jq -r '.outdir' "$PARAMS_FILE")"
