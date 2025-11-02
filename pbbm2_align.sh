#!/bin/bash

set -e  # Exit on error
set -u  # Exit on undefined variable
set -o pipefail  # Exit on pipe failure

SAMPLE_NAME="$1"
FL_BAM="$2"
S3_DESTINATION="$3"
SYNC_TIMESTAMP=$(date -u +%Y%m%dT%H%M%SZ)
S3_TARGET="${S3_DESTINATION%/}/${SAMPLE_NAME}_${SYNC_TIMESTAMP}"

REFERENCE_GENOME_GZ="genomes/GRCh38.p14.genome.fa.gz"

THREADS=32

OUTDIR="pbmm2_output"
mkdir -p "${OUTDIR}"
OUTDIR="$(cd "${OUTDIR}" && pwd)"

MAPPED_BAM="${OUTDIR}/${SAMPLE_NAME}.mapped.bam"
mkdir -p "$(dirname "${MAPPED_BAM}")"

echo "Running pbmm2:"
echo "  reference: ${REFERENCE_GENOME_GZ}"
echo "  input BAM: ${FL_BAM}"
echo "  output BAM: ${MAPPED_BAM}"
echo "  threads: ${THREADS}"

pbmm2 align "${REFERENCE_GENOME_GZ}" "${FL_BAM}" \
    "${MAPPED_BAM}" \
    --preset ISOSEQ \
    --sort \
    -j "${THREADS}" \
    --log-level INFO

# Index BAM file
echo "Indexing BAM file..."
samtools index -@ "${THREADS}" "${MAPPED_BAM}"

echo "syncing to s3: ${S3_TARGET}"
aws s3 sync "${OUTDIR}" "${S3_TARGET}" --quiet
