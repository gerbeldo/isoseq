#!/bin/bash
set -euo pipefail


# Usage: ./run-minimap2.sh <sample_name> <fasta_path> <reference_genome> <annotation_gtf> <s3_bucket>

# Check if correct number of arguments provided
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <sample_name> <fasta_path> <reference_genome> <annotation_gtf> <s3_bucket>"
    echo ""
    echo "Arguments:"
    echo "  sample_name       : Name of the sample (used for output naming)"
    echo "  fasta_path        : Path to input FASTA file (can be concatenated)"
    echo "  reference_genome  : Path to reference genome FASTA"
    echo "  annotation_gtf    : Path to annotation file (e.g., GENCODE v48 GTF) - for reference only"
    echo "  s3_bucket         : S3 bucket path (e.g., s3://my-bucket/complete-genomics-results)"
    exit 1
fi

# Assign input arguments
SAMPLE_NAME=$1
FASTA_PATH=$2
REFERENCE=$3
ANNOTATION=$4
S3_BUCKET=$5

# Generate timestamp for S3 folder
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
S3_OUTPUT_PATH="${S3_BUCKET}/${SAMPLE_NAME}_${TIMESTAMP}"

# Output directory and files
OUT_DIR="alignment_${SAMPLE_NAME}"
mkdir -p "${OUT_DIR}"

BAM_OUT="${OUT_DIR}/${SAMPLE_NAME}_aligned.bam"
UNMAPPED_BAM="${OUT_DIR}/${SAMPLE_NAME}_unmapped.bam"
LOG_FILE="${OUT_DIR}/${SAMPLE_NAME}_alignment.log"

echo "Starting alignment for sample: ${SAMPLE_NAME}" | tee -a "${LOG_FILE}"
echo "Input FASTA: ${FASTA_PATH}" | tee -a "${LOG_FILE}"
echo "Reference: ${REFERENCE}" | tee -a "${LOG_FILE}"
echo "Annotation: ${ANNOTATION}" | tee -a "${LOG_FILE}"
echo "Output directory: ${OUT_DIR}" | tee -a "${LOG_FILE}"
echo "S3 output path: ${S3_OUTPUT_PATH}" | tee -a "${LOG_FILE}"

# Run minimap2 alignment
# Parameters explained:
# -ax splice:hq    : Preset for high-quality long-read spliced alignment (for Complete Genomics data)
# -t 8             : Number of threads (adjust based on your system)
# -G 500k          : Maximum intron length (500kb; default is 200k)
# -p 0.5           : Minimum secondary-to-primary score ratio (default 0.8; 0.5 is more permissive)
# -N 10            : Retain at most 10 secondary alignments (default is 5)
# --cs=long        : Include long-format cs tag for detailed alignment information
# -a               : Output in SAM format (includes unmapped reads)
# -Y               : Use soft clipping for supplementary alignments
# --eqx            : Write =/X CIGAR operators for match/mismatch distinction
# -u f             : Only consider forward transcript strand (use if strand-specific)

echo "Running minimap2 alignment..." | tee -a "${LOG_FILE}"

minimap2 \
    -ax splice:hq \
    -t 8 \
    -G 500k \
    -p 0.5 \
    -N 10 \
    --cs=long \
    --eqx \
    -Y \
    -a \
    "${REFERENCE}" \
    "${FASTA_PATH}" \
    2>> "${LOG_FILE}" | \
samtools view -bh - | \
samtools sort -@ 4 -m 4G -o "${BAM_OUT}" -

# Index the BAM file
echo "Indexing BAM file..." | tee -a "${LOG_FILE}"
samtools index "${BAM_OUT}"

# Extract unmapped reads to separate BAM file
echo "Extracting unmapped reads..." | tee -a "${LOG_FILE}"
samtools view -b -f 4 "${BAM_OUT}" > "${UNMAPPED_BAM}"
samtools index "${UNMAPPED_BAM}"

# Generate alignment statistics
echo "Generating alignment statistics..." | tee -a "${LOG_FILE}"
samtools flagstat "${BAM_OUT}" > "${OUT_DIR}/${SAMPLE_NAME}_flagstat.txt"
samtools stats "${BAM_OUT}" > "${OUT_DIR}/${SAMPLE_NAME}_stats.txt"

# Count unmapped reads
UNMAPPED_COUNT=$(samtools view -c -f 4 "${BAM_OUT}")
TOTAL_COUNT=$(samtools view -c "${BAM_OUT}")
echo "Total reads: ${TOTAL_COUNT}" | tee -a "${LOG_FILE}"
echo "Unmapped reads: ${UNMAPPED_COUNT}" | tee -a "${LOG_FILE}"

echo "Alignment complete!" | tee -a "${LOG_FILE}"
echo "Output BAM: ${BAM_OUT}" | tee -a "${LOG_FILE}"
echo "Unmapped BAM: ${UNMAPPED_BAM}" | tee -a "${LOG_FILE}"

# Sync results to S3
echo "Syncing results to S3: ${S3_OUTPUT_PATH}" | tee -a "${LOG_FILE}"
aws s3 sync --quiet "${OUT_DIR}" "${S3_OUTPUT_PATH}"

if [ $? -eq 0 ]; then
    echo "S3 sync completed successfully!" | tee -a "${LOG_FILE}"
else
    echo "WARNING: S3 sync failed. Check AWS credentials and bucket permissions." | tee -a "${LOG_FILE}"
    exit 1
fi

echo "Pipeline complete. All results synced to: ${S3_OUTPUT_PATH}" | tee -a "${LOG_FILE}"
