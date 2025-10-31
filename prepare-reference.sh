#!/bin/bash

################################################################################
# Prepare reference assets for the Iso-Seq workflow.
# Sorts the annotation GTF and ensures the reference FASTA is indexed so that
# downstream pigeon commands can run without additional setup.
#
# Usage: ./prepare-reference.sh <annotation_gtf> <reference_fasta>
################################################################################

set -euo pipefail

if [ $# -ne 2 ]; then
    echo "Usage: $0 <annotation_gtf> <reference_fasta>"
    exit 1
fi

ANNOTATION_GTF="$1"
REFERENCE_FASTA="$2"

if [ ! -f "${ANNOTATION_GTF}" ]; then
    echo "ERROR: Annotation GTF not found: ${ANNOTATION_GTF}"
    exit 1
fi

if [[ "${ANNOTATION_GTF}" == *.gz ]]; then
    echo "ERROR: Annotation GTF must be uncompressed for pigeon: ${ANNOTATION_GTF}"
    exit 1
fi

if [ ! -f "${REFERENCE_FASTA}" ]; then
    echo "ERROR: Reference FASTA not found: ${REFERENCE_FASTA}"
    exit 1
fi

if [[ "${REFERENCE_FASTA}" == *.gz ]]; then
    echo "ERROR: Reference FASTA must be uncompressed for pigeon: ${REFERENCE_FASTA}"
    exit 1
fi

command -v pigeon >/dev/null 2>&1 || { echo "ERROR: pigeon not found in PATH."; exit 1; }

echo "Preparing annotation and reference for pigeon..."
echo "  Annotation: ${ANNOTATION_GTF}"
echo "  Reference:  ${REFERENCE_FASTA}"

pigeon prepare "${ANNOTATION_GTF}" "${REFERENCE_FASTA}"

ANNOTATION_SORTED="${ANNOTATION_GTF}.sorted.gtf"
REFERENCE_FAI="${REFERENCE_FASTA}.fai"

if [ ! -f "${ANNOTATION_SORTED}" ]; then
    echo "WARNING: Expected sorted annotation not found at ${ANNOTATION_SORTED}."
else
    echo "Sorted annotation generated: ${ANNOTATION_SORTED}"
fi

if [ ! -f "${REFERENCE_FAI}" ]; then
    echo "WARNING: FASTA index not found at ${REFERENCE_FAI}. pigeon prepare should generate this."
else
    echo "FASTA index generated: ${REFERENCE_FAI}"
fi

echo "Reference preparation complete."
