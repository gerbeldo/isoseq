I#!/bin/bash

set -euo pipefail

# Check if all required arguments are provided
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <sample_name> <mapped_bam_file> <genome_reference_file> <s3_destination>"
    echo ""
    echo "Example: $0 sample1 mapped.bam genome.fa s3://my-bucket/results"
    exit 1
fi

# sample name used as output prefix
SAMPLE_NAME="$1"
MAPPED_BAM_FILE="$2"
GENOME_REFERENCE="$3"
S3_DESTINATION="$4"

OUTDIR="longcallr_output"
SYNC_TIMESTAMP=$(date -u +%Y%m%dT%H%M%SZ)
S3_TARGET="${S3_DESTINATION%/}/${SAMPLE_NAME}_${SYNC_TIMESTAMP}"

# Create output directory if it doesn't exist
mkdir -p "${OUTDIR}"

echo "==================================================================="
echo "Running longcallR for sample: ${SAMPLE_NAME}"
echo "==================================================================="
echo "Input BAM file: $MAPPED_BAM_FILE"
echo "Genome reference file: $GENOME_REFERENCE"
echo "Output directory: $OUTDIR"
echo "S3 sync target: ${S3_TARGET}"
echo "==================================================================="

longcallr \
    -b "${MAPPED_BAM_FILE}" \
    -f "${GENOME_REFERENCE}" \
    -p hifi-isoseq \
    -t $(nproc)  \
    -o "${OUTDIR}/${SAMPLE_NAME}"

echo "Syncing results to S3: ${S3_TARGET}"
aws s3 sync "${OUTDIR}" "${S3_TARGET}" --quiet

echo "Script finished."
