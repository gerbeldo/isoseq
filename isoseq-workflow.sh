#!/bin/bash

################################################################################
# PacBio Iso-Seq Bulk RNA-seq Analysis Pipeline (Single Sample)
# For processing individual samples in parallel
# Input: Full-length (FL/FLNC) reads BAM generated after primer removal/refine
#
# Usage: ./isoseq_pipeline.sh <sample_name> <fl_bam> <s3_destination>
# Example: ./isoseq_pipeline.sh disease_sample disease.refined.fl.bam s3://my-bucket/analysis
################################################################################

set -e  # Exit on error
set -u  # Exit on undefined variable
set -o pipefail  # Exit on pipe failure

################################################################################
# PARSE COMMAND LINE ARGUMENTS
################################################################################

if [ $# -ne 3 ]; then
    echo "ERROR: Incorrect number of arguments"
    echo ""
    echo "Usage: $0 <sample_name> <fl_bam> <s3_destination>"
    echo ""
    echo "Arguments:"
    echo "  sample_name      Name for this sample (e.g., 'disease', 'control')"
    echo "  fl_bam           Path to input full-length (FL/FLNC) reads BAM file"
    echo "  s3_destination   S3 bucket + prefix for results (e.g., s3://bucket/path)"
    echo ""
    echo "Example:"
    echo "  $0 disease disease_sample.refined.fl.bam s3://my-bucket/analysis"
    echo "  $0 control control_sample.refined.fl.bam s3://my-bucket/analysis"
    echo ""
    exit 1
fi

SAMPLE_NAME="$1"
FL_BAM="$2"
S3_DESTINATION="$3"
SYNC_TIMESTAMP=$(date -u +%Y%m%dT%H%M%SZ)
S3_TARGET="${S3_DESTINATION%/}/${SAMPLE_NAME}_${SYNC_TIMESTAMP}"

################################################################################
# CONFIGURATION - MODIFY THESE VARIABLES AS NEEDED
################################################################################

# Reference assets directory and filenames
GENOME_DIR="/home/ubuntu/genomes"
REFERENCE_GENOME_BASENAME="GRCh38.p14.genome.fa"
ANNOTATION_GTF_BASENAME="gencode.v48.primary_assembly.annotation.gtf"

# Reference genome files (provide both gzipped and uncompressed paths)
REFERENCE_GENOME_GZ="${GENOME_DIR}/${REFERENCE_GENOME_BASENAME}.gz"
REFERENCE_GENOME_FA="${GENOME_DIR}/${REFERENCE_GENOME_BASENAME}"

# sorted references
ANNOTATION_GTF_SORTED="${GENOME_DIR}/gencode.v48.primary_assembly.annotation.sorted.gtf"
REFERENCE_GENOME_FAI="${REFERENCE_GENOME_FA}.fai"

# Annotation files (provide both gzipped and uncompressed paths)
ANNOTATION_GTF_GZ="${GENOME_DIR}/${ANNOTATION_GTF_BASENAME}.gz"
ANNOTATION_GTF="${GENOME_DIR}/${ANNOTATION_GTF_BASENAME}"

# Output directory (will create subdirectory for this sample)
OUTDIR="isoseq_output/${SAMPLE_NAME}"
mkdir -p "${OUTDIR}"
OUTDIR="$(cd "${OUTDIR}" && pwd)"

# Thread/CPU settings
THREADS=16

# Prerequisite:
#   Run ./prepare-reference.sh "${ANNOTATION_GTF}" "${REFERENCE_GENOME_FA}" once
#   before executing this workflow so pigeon indexes are available.

# Log file
LOGFILE="${OUTDIR}/${SAMPLE_NAME}_pipeline.log"

################################################################################
# VALIDATE INPUTS
################################################################################

echo "==================================================================="
echo "PacBio Iso-Seq Pipeline - Sample: ${SAMPLE_NAME}"
echo "==================================================================="
echo "Started at: $(date)"
echo ""

# Check if input file exists
if [ ! -f "${FL_BAM}" ]; then
    echo "ERROR: Input full-length BAM file not found: ${FL_BAM}"
    exit 1
fi

if [[ "${S3_DESTINATION}" != s3://* ]]; then
    echo "ERROR: S3 destination must start with s3:// (got: ${S3_DESTINATION})"
    exit 1
fi

# Check if reference files exist
if [ ! -f "${REFERENCE_GENOME_GZ}" ]; then
    echo "ERROR: gzipped reference genome not found: ${REFERENCE_GENOME_GZ}"
    exit 1
fi

if [[ "${REFERENCE_GENOME_GZ}" != *.gz ]]; then
    echo "ERROR: REFERENCE_GENOME_GZ must point to a .gz file: ${REFERENCE_GENOME_GZ}"
    exit 1
fi

if [ ! -f "${REFERENCE_GENOME_FA}" ]; then
    echo "ERROR: uncompressed reference genome not found: ${REFERENCE_GENOME_FA}"
    exit 1
fi

if [[ "${REFERENCE_GENOME_FA}" == *.gz ]]; then
    echo "ERROR: REFERENCE_GENOME_FA must point to an uncompressed FASTA: ${REFERENCE_GENOME_FA}"
    exit 1
fi

if [ ! -f "${ANNOTATION_GTF_GZ}" ]; then
    echo "ERROR: gzipped annotation GTF not found: ${ANNOTATION_GTF_GZ}"
    exit 1
fi

if [[ "${ANNOTATION_GTF_GZ}" != *.gz ]]; then
    echo "ERROR: ANNOTATION_GTF_GZ must point to a .gz file: ${ANNOTATION_GTF_GZ}"
    exit 1
fi

if [ ! -f "${ANNOTATION_GTF}" ]; then
    echo "ERROR: uncompressed annotation GTF not found: ${ANNOTATION_GTF}"
    exit 1
fi

if [[ "${ANNOTATION_GTF}" == *.gz ]]; then
    echo "ERROR: ANNOTATION_GTF must point to an uncompressed GTF: ${ANNOTATION_GTF}"
    exit 1
fi

################################################################################
# STEP 0: Reference path setup
################################################################################

echo "Using reference assets:"
echo "  Reference directory:       ${GENOME_DIR}"
echo "  Annotation (gzipped):      ${ANNOTATION_GTF_GZ}"
echo "  Annotation (uncompressed): ${ANNOTATION_GTF}"
echo "  Reference (gzipped):       ${REFERENCE_GENOME_GZ}"
echo "  Reference (uncompressed):  ${REFERENCE_GENOME_FA}"
echo ""

################################################################################
# STEP 0: Setup and Prerequisites
################################################################################

echo "==================================================================="
echo "Step 0: Setting up directories and checking prerequisites"
echo "==================================================================="

mkdir -p ${OUTDIR}/{03_cluster,04_mapping,05_collapse,06_classification}

TMP_BASE="${OUTDIR}/tmp"
mkdir -p "${TMP_BASE}"
case "${TMP_BASE}" in
    */) TMPDIR="${TMP_BASE}" ;;
    *)  TMPDIR="${TMP_BASE}/" ;;
esac
export TMPDIR
export TMP="${TMPDIR}"
export TEMP="${TMPDIR}"

# Redirect all output to log file while still showing on screen
exec > >(tee -a ${LOGFILE})
exec 2>&1

echo "Sample name: ${SAMPLE_NAME}"
echo "Input file: ${FL_BAM}"
echo "Output directory: ${OUTDIR}"
echo "Reference genome (pbmm2): ${REFERENCE_GENOME_GZ}"
echo "Reference genome (pigeon): ${REFERENCE_GENOME_FA}"
echo "Annotation GTF (gzipped): ${ANNOTATION_GTF_GZ}"
echo "Annotation GTF (pigeon): ${ANNOTATION_GTF}"
echo "S3 destination prefix: ${S3_DESTINATION}"
echo "S3 sync target: ${S3_TARGET}"
echo "Threads: ${THREADS}"
echo "Temporary directory: ${TMPDIR}"
echo ""

# Check if required tools are installed
command -v isoseq >/dev/null 2>&1 || { echo "ERROR: isoseq not found. Install via: conda install isoseq"; exit 1; }
command -v pbmm2 >/dev/null 2>&1 || { echo "ERROR: pbmm2 not found. Install via: conda install pbmm2"; exit 1; }
command -v pigeon >/dev/null 2>&1 || { echo "ERROR: pigeon not found. Install via: conda install pbpigeon"; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools not found. Install via: conda install samtools"; exit 1; }
command -v aws >/dev/null 2>&1 || { echo "ERROR: aws CLI not found. Install via: https://docs.aws.amazon.com/cli/"; exit 1; }

echo "All required tools found."
echo ""

FLNC_BAM="${FL_BAM}"

sync_to_s3() {
    local phase="$1"
    echo "-------------------------------------------------------------------"
    echo "Syncing output directory to S3 (${phase})"
    echo "Destination: ${S3_TARGET}"
    aws s3 sync "${OUTDIR}" "${S3_TARGET}"
    echo "S3 sync completed at: $(date)"
    echo ""
}


################################################################################
# STEP 1: Clustering - Isoform-level clustering (isoseq cluster2)
################################################################################

echo "==================================================================="
echo "Step 1: Clustering isoforms with cluster2"
echo "==================================================================="
echo "Started at: $(date)"

isoseq cluster2 ${FLNC_BAM} \
    ${OUTDIR}/03_cluster/${SAMPLE_NAME}.clustered.bam \
    --num-threads ${THREADS} \
    --log-level DEBUG

CLUSTERED_BAM="${OUTDIR}/03_cluster/${SAMPLE_NAME}.clustered.bam"

if [ ! -f "${CLUSTERED_BAM}" ]; then
    echo "ERROR: Cluster2 did not produce expected output: ${CLUSTERED_BAM}"
    exit 1
fi

echo "Clustering completed at: $(date)"
sync_to_s3 "Step 1 - clustering"
echo ""

################################################################################
# STEP 2: Mapping - Align to reference genome (pbmm2)
################################################################################

echo "==================================================================="
echo "Step 2: Mapping transcripts to reference genome with pbmm2"
echo "==================================================================="
echo "Started at: $(date)"

MAPPED_BAM="${OUTDIR}/04_mapping/${SAMPLE_NAME}.mapped.bam"
mkdir -p "$(dirname "${MAPPED_BAM}")"

pbmm2 align ${REFERENCE_GENOME_GZ} ${CLUSTERED_BAM} \
    ${MAPPED_BAM} \
    --preset ISOSEQ \
    --sort \
    -j ${THREADS} \
    --log-level DEBUG

if [ ! -f "${MAPPED_BAM}" ]; then
    echo "ERROR: pbmm2 did not produce expected output: ${MAPPED_BAM}"
    exit 1
fi

# Index BAM file
echo "Indexing BAM file..."
samtools index ${MAPPED_BAM}

echo "Mapping completed at: $(date)"
sync_to_s3 "Step 2 - mapping"
echo ""

################################################################################
# STEP 3: Collapse - Remove redundant transcripts (isoseq collapse)
################################################################################

echo "==================================================================="
echo "Step 3: Collapsing redundant transcripts"
echo "==================================================================="
echo "Started at: $(date)"

isoseq collapse ${MAPPED_BAM} ${FLNC_BAM} \
    ${OUTDIR}/05_collapse/${SAMPLE_NAME}.collapsed.gff \
    --do-not-collapse-extra-5exons \
    --log-level DEBUG

COLLAPSED_GFF="${OUTDIR}/05_collapse/${SAMPLE_NAME}.collapsed.gff"

if [ ! -f "${COLLAPSED_GFF}" ]; then
    echo "ERROR: Collapse did not produce expected output: ${COLLAPSED_GFF}"
    exit 1
fi

echo "Collapse completed at: $(date)"
sync_to_s3 "Step 3 - collapse"
echo ""

################################################################################
# STEP 4: Verify prepared reference files
################################################################################

echo "==================================================================="
echo "Step 4: Verifying prepared reference files for pigeon"
echo "==================================================================="

if [ ! -f "${ANNOTATION_GTF_SORTED}" ] || [ ! -f "${REFERENCE_GENOME_FAI}" ]; then
    echo "ERROR: Prepared reference files not found."
    echo "Run ./prepare-reference.sh \"${ANNOTATION_GTF}\" \"${REFERENCE_GENOME_FA}\" before executing this workflow."
    exit 1
fi

echo "Reference files verified."
sync_to_s3 "Step 4 - reference verification"
echo ""

################################################################################
# STEP 5: Sort transcript GFF files
################################################################################

echo "==================================================================="
echo "Step 5: Sorting transcript GFF file"
echo "==================================================================="
echo "Started at: $(date)"

pigeon prepare ${COLLAPSED_GFF}

SORTED_GFF="${OUTDIR}/05_collapse/${SAMPLE_NAME}.collapsed.sorted.gff"

if [ ! -f "${SORTED_GFF}" ]; then
    echo "ERROR: Pigeon prepare did not produce expected output: ${SORTED_GFF}"
    exit 1
fi

echo "GFF sorting completed at: $(date)"
sync_to_s3 "Step 5 - GFF sorting"
echo ""

################################################################################
# STEP 6: Classification - Classify isoforms (pigeon classify)
################################################################################

echo "==================================================================="
echo "Step 6: Classifying transcripts with pigeon"
echo "==================================================================="
echo "Started at: $(date)"

pigeon classify ${SORTED_GFF} ${ANNOTATION_GTF_SORTED} ${REFERENCE_GENOME_FA} \
    --flnc ${OUTDIR}/05_collapse/${SAMPLE_NAME}.collapsed.flnc_count.txt \
    --log-level DEBUG \
    --out-dir ${OUTDIR}/06_classification

CLASSIFICATION_FILE="${OUTDIR}/06_classification/${SAMPLE_NAME}_classification.txt"

if [ ! -f "${CLASSIFICATION_FILE}" ]; then
    echo "ERROR: Pigeon classify did not produce expected output: ${CLASSIFICATION_FILE}"
    exit 1
fi

echo "Classification completed at: $(date)"
sync_to_s3 "Step 6 - classification"
echo ""

################################################################################
# STEP 7: Filter - Remove low-quality transcripts (pigeon filter)
################################################################################

echo "==================================================================="
echo "Step 7: Filtering transcripts with pigeon"
echo "==================================================================="
echo "Started at: $(date)"

pigeon filter \
    ${CLASSIFICATION_FILE} \
    --isoforms ${SORTED_GFF} \
    --log-level DEBUG


FILTERED_CLASSIFICATION="${OUTDIR}/06_classification/${SAMPLE_NAME}_classification.filtered_lite_classification.txt"

if [ ! -f "${FILTERED_CLASSIFICATION}" ]; then
    echo "WARNING: Pigeon filter did not produce expected output: ${FILTERED_CLASSIFICATION}"
    echo "This may be normal if no transcripts were filtered."
fi

echo "Filtering completed at: $(date)"
sync_to_s3 "Step 7 - filtering"
echo ""

################################################################################
# STEP 8: Generate summary report
################################################################################

echo "==================================================================="
echo "Step 8: Generating summary report"
echo "==================================================================="

SUMMARY="${OUTDIR}/${SAMPLE_NAME}_analysis_summary.txt"

cat > ${SUMMARY} << SUMMARY_EOF
================================================================================
PacBio Iso-Seq Analysis Summary
================================================================================

Sample Name: ${SAMPLE_NAME}
Input File: ${FL_BAM}
Analysis Date: $(date)

================================================================================
Pipeline Statistics
================================================================================

SUMMARY_EOF

# Extract statistics from log files
echo "Processing Statistics:" >> ${SUMMARY}
echo "" >> ${SUMMARY}

# Input full-length reads
if [ -f "${FLNC_BAM}" ]; then
    echo "1. Input FL/FLNC reads:" >> ${SUMMARY}
    FLNC_COUNT=$(samtools view -c ${FLNC_BAM} 2>/dev/null || echo "N/A")
    echo "   Full-length reads available for clustering: ${FLNC_COUNT}" >> ${SUMMARY}
    echo "" >> ${SUMMARY}
fi

# Clustered transcripts
if [ -f "${CLUSTERED_BAM}" ]; then
    echo "2. Clustering:" >> ${SUMMARY}
    CLUSTERED_COUNT=$(samtools view -c ${CLUSTERED_BAM} 2>/dev/null || echo "N/A")
    echo "   Clustered transcripts: ${CLUSTERED_COUNT}" >> ${SUMMARY}
    echo "" >> ${SUMMARY}
fi

# Mapped transcripts
if [ -f "${MAPPED_BAM}" ]; then
    echo "3. Mapping:" >> ${SUMMARY}
    MAPPED_COUNT=$(samtools view -c -F 4 ${MAPPED_BAM} 2>/dev/null || echo "N/A")
    UNMAPPED_COUNT=$(samtools view -c -f 4 ${MAPPED_BAM} 2>/dev/null || echo "N/A")
    echo "   Mapped transcripts: ${MAPPED_COUNT}" >> ${SUMMARY}
    echo "   Unmapped transcripts: ${UNMAPPED_COUNT}" >> ${SUMMARY}
    echo "" >> ${SUMMARY}
fi

# Collapsed transcripts
if [ -f "${COLLAPSED_GFF}" ]; then
    echo "4. Collapse:" >> ${SUMMARY}
    COLLAPSED_COUNT=$(grep -v "^#" ${COLLAPSED_GFF} | wc -l)
    echo "   Collapsed transcripts (unique isoforms): ${COLLAPSED_COUNT}" >> ${SUMMARY}
    echo "" >> ${SUMMARY}
fi

# Classified transcripts
if [ -f "${CLASSIFICATION_FILE}" ]; then
    echo "5. Classification:" >> ${SUMMARY}
    TOTAL_TRANSCRIPTS=$(tail -n +2 ${CLASSIFICATION_FILE} | wc -l)
    echo "   Total classified transcripts: ${TOTAL_TRANSCRIPTS}" >> ${SUMMARY}

    # Breakdown by category
    echo "" >> ${SUMMARY}
    echo "   Breakdown by structural category:" >> ${SUMMARY}
    tail -n +2 ${CLASSIFICATION_FILE} | cut -f6 | sort | uniq -c | sort -rn | awk '{printf "     %-30s: %s\n", $2, $1}' >> ${SUMMARY}
    echo "" >> ${SUMMARY}
fi

# Filtered transcripts
if [ -f "${FILTERED_CLASSIFICATION}" ]; then
    echo "6. Filtering:" >> ${SUMMARY}
    FILTERED_COUNT=$(tail -n +2 ${FILTERED_CLASSIFICATION} | wc -l)
    echo "   High-quality filtered transcripts: ${FILTERED_COUNT}" >> ${SUMMARY}
    echo "" >> ${SUMMARY}
fi

cat >> ${SUMMARY} << SUMMARY_EOF

================================================================================
Output Files
================================================================================

Key output files for downstream analysis:

1. Classification results:
   ${CLASSIFICATION_FILE}

2. Filtered classification (high-quality only):
   ${FILTERED_CLASSIFICATION}

3. Transcript structures (GFF):
   ${COLLAPSED_GFF}
   ${SORTED_GFF}

4. Mapped reads (BAM):
   ${MAPPED_BAM}
   ${MAPPED_BAM}.bai

5. Full analysis log:
   ${LOGFILE}

================================================================================
Pipeline completed successfully!
Finished at: $(date)
================================================================================
SUMMARY_EOF

# Display summary
cat ${SUMMARY}
sync_to_s3 "Step 8 - summary"

################################################################################
# COMPLETION
################################################################################

echo ""
echo "==================================================================="
echo "Pipeline completed successfully for sample: ${SAMPLE_NAME}"
echo "==================================================================="
echo ""
echo "Total runtime: $SECONDS seconds"
echo ""
echo "Summary saved to: ${SUMMARY}"
echo "Full log saved to: ${LOGFILE}"
echo "S3 sync destination: ${S3_TARGET}"
echo ""
