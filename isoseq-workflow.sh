#!/bin/bash

################################################################################
# PacBio Iso-Seq Bulk RNA-seq Analysis Pipeline (Single Sample)
# For processing individual samples in parallel
# Input: HiFi reads from PacBio Revio
#
# Usage: ./isoseq_pipeline.sh <sample_name> <hifi_reads.bam>
# Example: ./isoseq_pipeline.sh disease_sample disease.hifi_reads.bam
################################################################################

set -e  # Exit on error
set -u  # Exit on undefined variable
set -o pipefail  # Exit on pipe failure

################################################################################
# PARSE COMMAND LINE ARGUMENTS
################################################################################

if [ $# -ne 2 ]; then
    echo "ERROR: Incorrect number of arguments"
    echo ""
    echo "Usage: $0 <sample_name> <hifi_reads_bam>"
    echo ""
    echo "Arguments:"
    echo "  sample_name      Name for this sample (e.g., 'disease', 'control')"
    echo "  hifi_reads_bam   Path to input HiFi reads BAM file"
    echo ""
    echo "Example:"
    echo "  $0 disease disease_sample.hifi_reads.bam"
    echo "  $0 control control_sample.hifi_reads.bam"
    echo ""
    exit 1
fi

SAMPLE_NAME="$1"
HIFI_READS="$2"

################################################################################
# CONFIGURATION - MODIFY THESE VARIABLES AS NEEDED
################################################################################

# Reference assets directory and filenames
GENOME_DIR="/home/ubuntu/genomes"
REFERENCE_GENOME_BASENAME="GRCh38.p14.genome.fa"
ANNOTATION_GTF_BASENAME="gencode.v49.primary_assembly.annotation.gtf"

# Reference genome files (provide both gzipped and uncompressed paths)
REFERENCE_GENOME_GZ="${GENOME_DIR}/${REFERENCE_GENOME_BASENAME}.gz"
REFERENCE_GENOME_FA="${GENOME_DIR}/${REFERENCE_GENOME_BASENAME}"

# sorted references
ANNOTATION_GTF_SORTED="${GENOME_DIR}/gencode.v49.primary_assembly.annotation.gtf"
REFERENCE_GENOME_FAI="${REFERENCE_GENOME_FA}.fai"

# Annotation files (provide both gzipped and uncompressed paths)
ANNOTATION_GTF_GZ="${GENOME_DIR}/${ANNOTATION_GTF_BASENAME}.gz"
ANNOTATION_GTF="${GENOME_DIR}/${ANNOTATION_GTF_BASENAME}"

# Primer file for Iso-Seq (standard primers)
PRIMERS="/home/ubuntu/data/ismb_workshop/primers.fasta"

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
if [ ! -f "${HIFI_READS}" ]; then
    echo "ERROR: Input HiFi reads file not found: ${HIFI_READS}"
    exit 1
fi

# Check if reference files exist
if [ ! -f "${REFERENCE_GENOME_GZ}" ]; then
    echo "ERROR: gzipped reference genome not found: ${REFERENCE_GENOME_GZ}"
    echo "Download from: https://www.gencodegenes.org/human/"
    exit 1
fi

if [[ "${REFERENCE_GENOME_GZ}" != *.gz ]]; then
    echo "ERROR: REFERENCE_GENOME_GZ must point to a .gz file: ${REFERENCE_GENOME_GZ}"
    exit 1
fi

if [ ! -f "${REFERENCE_GENOME_FA}" ]; then
    echo "ERROR: uncompressed reference genome not found: ${REFERENCE_GENOME_FA}"
    echo "Download from: https://www.gencodegenes.org/human/"
    exit 1
fi

if [[ "${REFERENCE_GENOME_FA}" == *.gz ]]; then
    echo "ERROR: REFERENCE_GENOME_FA must point to an uncompressed FASTA: ${REFERENCE_GENOME_FA}"
    exit 1
fi

if [ ! -f "${ANNOTATION_GTF_GZ}" ]; then
    echo "ERROR: gzipped annotation GTF not found: ${ANNOTATION_GTF_GZ}"
    echo "Download from: https://www.gencodegenes.org/human/"
    exit 1
fi

if [[ "${ANNOTATION_GTF_GZ}" != *.gz ]]; then
    echo "ERROR: ANNOTATION_GTF_GZ must point to a .gz file: ${ANNOTATION_GTF_GZ}"
    exit 1
fi

if [ ! -f "${ANNOTATION_GTF}" ]; then
    echo "ERROR: uncompressed annotation GTF not found: ${ANNOTATION_GTF}"
    echo "Download from: https://www.gencodegenes.org/human/"
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

mkdir -p ${OUTDIR}/{01_primers,02_refine,03_cluster,04_mapping,05_collapse,06_classification}

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
echo "Input file: ${HIFI_READS}"
echo "Output directory: ${OUTDIR}"
echo "Reference genome (pbmm2): ${REFERENCE_GENOME_GZ}"
echo "Reference genome (pigeon): ${REFERENCE_GENOME_FA}"
echo "Annotation GTF (gzipped): ${ANNOTATION_GTF_GZ}"
echo "Annotation GTF (pigeon): ${ANNOTATION_GTF}"
echo "Threads: ${THREADS}"
echo "Temporary directory: ${TMPDIR}"
echo ""

# Check if required tools are installed
command -v lima >/dev/null 2>&1 || { echo "ERROR: lima not found. Install via: conda install lima"; exit 1; }
command -v isoseq >/dev/null 2>&1 || { echo "ERROR: isoseq not found. Install via: conda install isoseq"; exit 1; }
command -v pbmm2 >/dev/null 2>&1 || { echo "ERROR: pbmm2 not found. Install via: conda install pbmm2"; exit 1; }
command -v pigeon >/dev/null 2>&1 || { echo "ERROR: pigeon not found. Install via: conda install pbpigeon"; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools not found. Install via: conda install samtools"; exit 1; }

echo "All required tools found."
echo ""

# Create standard primer file if it doesn't exist
if [ ! -f "${PRIMERS}" ]; then
    echo "Creating standard Iso-Seq primer file: ${PRIMERS}"
    cat > ${PRIMERS} << 'EOF'
>primer_5p
AAGCAGTGGTATCAACGCAGAGTACATGGGG
>primer_3p
GTACTCTGCGTTGATACCACTGCTT
EOF
    echo "Primer file created."
    echo ""
fi

################################################################################
# STEP 1: Primer Removal (lima)
################################################################################

echo "==================================================================="
echo "Step 1: Removing primers with lima"
echo "==================================================================="
echo "Started at: $(date)"

lima ${HIFI_READS} ${PRIMERS} \
    ${OUTDIR}/01_primers/${SAMPLE_NAME}.fl.bam \
    --isoseq \
    --dump-clips \
    --peek-guess \
    --num-threads ${THREADS} \
    --log-level INFO \
    --log-file ${OUTDIR}/01_primers/${SAMPLE_NAME}.lima.log

# Get the output BAM file (lima creates primer_5p--primer_3p.bam)
FL_BAM="${OUTDIR}/01_primers/${SAMPLE_NAME}.fl.primer_5p--primer_3p.bam"

if [ ! -f "${FL_BAM}" ]; then
    echo "ERROR: Lima did not produce expected output: ${FL_BAM}"
    exit 1
fi

echo "Primer removal completed at: $(date)"
echo ""

################################################################################
# STEP 2: Refine - Remove polyA tails and concatemers (isoseq refine)
################################################################################

echo "==================================================================="
echo "Step 2: Refining reads - removing polyA tails and concatemers"
echo "==================================================================="
echo "Started at: $(date)"

isoseq refine ${FL_BAM} ${PRIMERS} \
    ${OUTDIR}/02_refine/${SAMPLE_NAME}.flnc.bam \
    --num-threads ${THREADS} \
    --log-level INFO \
    --log-file ${OUTDIR}/02_refine/${SAMPLE_NAME}.refine.log

FLNC_BAM="${OUTDIR}/02_refine/${SAMPLE_NAME}.flnc.bam"

if [ ! -f "${FLNC_BAM}" ]; then
    echo "ERROR: Refine did not produce expected output: ${FLNC_BAM}"
    exit 1
fi

echo "Refine completed at: $(date)"
echo ""

################################################################################
# STEP 3: Clustering - Isoform-level clustering (isoseq cluster2)
################################################################################

echo "==================================================================="
echo "Step 3: Clustering isoforms with cluster2"
echo "==================================================================="
echo "Started at: $(date)"

isoseq cluster2 ${FLNC_BAM} \
    ${OUTDIR}/03_cluster/${SAMPLE_NAME}.clustered.bam \
    --num-threads ${THREADS} \
    --log-level INFO

CLUSTERED_BAM="${OUTDIR}/03_cluster/${SAMPLE_NAME}.clustered.bam"

if [ ! -f "${CLUSTERED_BAM}" ]; then
    echo "ERROR: Cluster2 did not produce expected output: ${CLUSTERED_BAM}"
    exit 1
fi

echo "Clustering completed at: $(date)"
echo ""

################################################################################
# STEP 4: Mapping - Align to reference genome (pbmm2)
################################################################################

echo "==================================================================="
echo "Step 4: Mapping transcripts to reference genome with pbmm2"
echo "==================================================================="
echo "Started at: $(date)"

MAPPED_BAM="${OUTDIR}/04_mapping/${SAMPLE_NAME}.mapped.bam"
mkdir -p "$(dirname "${MAPPED_BAM}")"

pbmm2 align ${REFERENCE_GENOME_GZ} ${CLUSTERED_BAM} \
    ${MAPPED_BAM} \
    --preset ISOSEQ \
    --sort \
    -j ${THREADS} \
    --log-level INFO

if [ ! -f "${MAPPED_BAM}" ]; then
    echo "ERROR: pbmm2 did not produce expected output: ${MAPPED_BAM}"
    exit 1
fi

# Index BAM file
echo "Indexing BAM file..."
samtools index ${MAPPED_BAM}

echo "Mapping completed at: $(date)"
echo ""

################################################################################
# STEP 5: Collapse - Remove redundant transcripts (isoseq collapse)
################################################################################

echo "==================================================================="
echo "Step 5: Collapsing redundant transcripts"
echo "==================================================================="
echo "Started at: $(date)"

isoseq collapse ${MAPPED_BAM} ${FLNC_BAM} \
    ${OUTDIR}/05_collapse/${SAMPLE_NAME}.collapsed.gff \
    --do-not-collapse-extra-5exons \
    --log-level INFO

COLLAPSED_GFF="${OUTDIR}/05_collapse/${SAMPLE_NAME}.collapsed.gff"

if [ ! -f "${COLLAPSED_GFF}" ]; then
    echo "ERROR: Collapse did not produce expected output: ${COLLAPSED_GFF}"
    exit 1
fi

echo "Collapse completed at: $(date)"
echo ""

################################################################################
# STEP 6: Verify prepared reference files
################################################################################

echo "==================================================================="
echo "Step 6: Verifying prepared reference files for pigeon"
echo "==================================================================="

if [ ! -f "${ANNOTATION_GTF_SORTED}" ] || [ ! -f "${REFERENCE_GENOME_FAI}" ]; then
    echo "ERROR: Prepared reference files not found."
    echo "Run ./prepare-reference.sh \"${ANNOTATION_GTF}\" \"${REFERENCE_GENOME_FA}\" before executing this workflow."
    exit 1
fi

echo "Reference files verified."
echo ""

################################################################################
# STEP 7: Sort transcript GFF files
################################################################################

echo "==================================================================="
echo "Step 7: Sorting transcript GFF file"
echo "==================================================================="
echo "Started at: $(date)"

pigeon prepare ${COLLAPSED_GFF}

SORTED_GFF="${OUTDIR}/05_collapse/${SAMPLE_NAME}.collapsed.sorted.gff"

if [ ! -f "${SORTED_GFF}" ]; then
    echo "ERROR: Pigeon prepare did not produce expected output: ${SORTED_GFF}"
    exit 1
fi

echo "GFF sorting completed at: $(date)"
echo ""

################################################################################
# STEP 8: Classification - Classify isoforms (pigeon classify)
################################################################################

echo "==================================================================="
echo "Step 8: Classifying transcripts with pigeon"
echo "==================================================================="
echo "Started at: $(date)"

pigeon classify ${SORTED_GFF} ${ANNOTATION_GTF_SORTED} ${REFERENCE_GENOME_FA} \
    --flnc ${OUTDIR}/05_collapse/${SAMPLE_NAME}.collapsed.flnc_count.txt \
    --log-level INFO \
    --out-dir ${OUTDIR}/06_classification/${SAMPLE_NAME}

CLASSIFICATION_FILE="${OUTDIR}/06_classification/${SAMPLE_NAME}_classification.txt"

if [ ! -f "${CLASSIFICATION_FILE}" ]; then
    echo "ERROR: Pigeon classify did not produce expected output: ${CLASSIFICATION_FILE}"
    exit 1
fi

echo "Classification completed at: $(date)"
echo ""

################################################################################
# STEP 9: Filter - Remove low-quality transcripts (pigeon filter)
################################################################################

echo "==================================================================="
echo "Step 9: Filtering transcripts with pigeon"
echo "==================================================================="
echo "Started at: $(date)"

pigeon filter \
    ${CLASSIFICATION_FILE} \
    --isoforms ${SORTED_GFF}

FILTERED_CLASSIFICATION="${OUTDIR}/06_classification/${SAMPLE_NAME}_classification.filtered_lite_classification.txt"

if [ ! -f "${FILTERED_CLASSIFICATION}" ]; then
    echo "WARNING: Pigeon filter did not produce expected output: ${FILTERED_CLASSIFICATION}"
    echo "This may be normal if no transcripts were filtered."
fi

echo "Filtering completed at: $(date)"
echo ""

################################################################################
# STEP 10: Generate summary report
################################################################################

echo "==================================================================="
echo "Step 10: Generating summary report"
echo "==================================================================="

SUMMARY="${OUTDIR}/${SAMPLE_NAME}_analysis_summary.txt"

cat > ${SUMMARY} << SUMMARY_EOF
================================================================================
PacBio Iso-Seq Analysis Summary
================================================================================

Sample Name: ${SAMPLE_NAME}
Input File: ${HIFI_READS}
Analysis Date: $(date)

================================================================================
Pipeline Statistics
================================================================================

SUMMARY_EOF

# Extract statistics from log files
echo "Processing Statistics:" >> ${SUMMARY}
echo "" >> ${SUMMARY}

# Full-length reads from lima
if [ -f "${OUTDIR}/01_primers/${SAMPLE_NAME}.lima.log" ]; then
    echo "1. Primer Removal (lima):" >> ${SUMMARY}
    FL_READS=$(grep -c "^>" ${FL_BAM} 2>/dev/null || samtools view -c ${FL_BAM} 2>/dev/null || echo "N/A")
    echo "   Full-length reads after primer removal: ${FL_READS}" >> ${SUMMARY}
    echo "" >> ${SUMMARY}
fi

# FLNC reads from refine
if [ -f "${FLNC_BAM}" ]; then
    echo "2. Refine (FLNC reads):" >> ${SUMMARY}
    FLNC_COUNT=$(samtools view -c ${FLNC_BAM} 2>/dev/null || echo "N/A")
    echo "   Full-length non-chimeric (FLNC) reads: ${FLNC_COUNT}" >> ${SUMMARY}
    echo "" >> ${SUMMARY}
fi

# Clustered transcripts
if [ -f "${CLUSTERED_BAM}" ]; then
    echo "3. Clustering:" >> ${SUMMARY}
    CLUSTERED_COUNT=$(samtools view -c ${CLUSTERED_BAM} 2>/dev/null || echo "N/A")
    echo "   Clustered transcripts: ${CLUSTERED_COUNT}" >> ${SUMMARY}
    echo "" >> ${SUMMARY}
fi

# Mapped transcripts
if [ -f "${MAPPED_BAM}" ]; then
    echo "4. Mapping:" >> ${SUMMARY}
    MAPPED_COUNT=$(samtools view -c -F 4 ${MAPPED_BAM} 2>/dev/null || echo "N/A")
    UNMAPPED_COUNT=$(samtools view -c -f 4 ${MAPPED_BAM} 2>/dev/null || echo "N/A")
    echo "   Mapped transcripts: ${MAPPED_COUNT}" >> ${SUMMARY}
    echo "   Unmapped transcripts: ${UNMAPPED_COUNT}" >> ${SUMMARY}
    echo "" >> ${SUMMARY}
fi

# Collapsed transcripts
if [ -f "${COLLAPSED_GFF}" ]; then
    echo "5. Collapse:" >> ${SUMMARY}
    COLLAPSED_COUNT=$(grep -v "^#" ${COLLAPSED_GFF} | wc -l)
    echo "   Collapsed transcripts (unique isoforms): ${COLLAPSED_COUNT}" >> ${SUMMARY}
    echo "" >> ${SUMMARY}
fi

# Classified transcripts
if [ -f "${CLASSIFICATION_FILE}" ]; then
    echo "6. Classification:" >> ${SUMMARY}
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
    echo "7. Filtering:" >> ${SUMMARY}
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
Next Steps for Your Analysis
================================================================================

For your compound heterozygote study:

1. Extract your gene of interest from the classification file
2. Compare FL (full-length) read counts between disease and control samples
3. Examine isoform structures in IGV using the GFF and BAM files
4. Look for evidence of nonsense-mediated decay (reduced transcript levels)
5. For allele-specific expression, examine SNPs in the BAM file to phase variants

Recommended commands:

  # Extract your gene (replace GENE_NAME with your gene symbol)
  grep "GENE_NAME" ${CLASSIFICATION_FILE} > my_gene_classification.txt

  # View in IGV
  # Load: ${MAPPED_BAM} and ${COLLAPSED_GFF}
  # Navigate to your gene coordinates

================================================================================
Pipeline completed successfully!
Finished at: $(date)
================================================================================
SUMMARY_EOF

# Display summary
cat ${SUMMARY}

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
echo ""
