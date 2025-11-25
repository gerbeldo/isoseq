#!/bin/bash
set -euo pipefail

# usage: ./run_whatshap.sh <input.vcf> <sample.cram> <reference.fasta> <region>
# example: ./run_whatshap.sh variants.vcf sample.cram hg19.fasta "chr17:43044295-43170245"

if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <input.vcf> <sample.cram> <reference.fasta> <region>"
    echo "Example: $0 variants.vcf sample.cram hg19.fasta 'chr17:43044295-43170245'"
    exit 1
fi

VCF_GZ=$1
INPUT_CRAM=$2
REFERENCE=$3
REGION=$4
THREADS=$(nproc)

# Extract gene name from region for output naming (e.g., chr17:100-200 -> chr17_100_200)
GENE_NAME=$(echo "$REGION" | sed 's/[:-]/_/g')
OUTPUT_PREFIX="${GENE_NAME}_phased"

echo "=== WhatsHap Phasing Pipeline ==="
echo "Input VCF: $VCF_GZ"
echo "Input CRAM: $INPUT_CRAM"
echo "Reference: $REFERENCE"
echo "Region: $REGION"
echo "Threads: $THREADS"
echo ""


# extract relevant region to speed stuff up
echo "Extracting region ${REGION}..."
tabix -h "$VCF_GZ" "$REGION" > "${OUTPUT_PREFIX}_region.vcf"

# phase variats
echo "Phasing variants with WhatsHap..."
whatshap phase \
    -o "${OUTPUT_PREFIX}.vcf" \
    --reference "$REFERENCE" \
    "${OUTPUT_PREFIX}_region.vcf" \
    "$INPUT_CRAM"

echo "Compressing phased VCF..."
bgzip -c "${OUTPUT_PREFIX}.vcf" > "${OUTPUT_PREFIX}.vcf.gz"
tabix -p vcf "${OUTPUT_PREFIX}.vcf.gz"

# stats
echo "output phasing stats..."
whatshap stats "${OUTPUT_PREFIX}.vcf.gz" > "${OUTPUT_PREFIX}_stats.txt"

echo "add haplotype tag to visualize in IGV..."
whatshap haplotag \
    -o "${OUTPUT_PREFIX}_haplotagged.bam" \
    --reference "$REFERENCE" \
    --output-threads "$THREADS" \
    --regions "$REGION" \
    "${OUTPUT_PREFIX}.vcf.gz" \
    "$INPUT_CRAM"

# index
echo "Indexing haplotagged BAM..."
samtools index -@ "$THREADS" "${OUTPUT_PREFIX}_haplotagged.bam"


# print
echo ""
echo "=== Phasing Statistics ==="
cat "${OUTPUT_PREFIX}_stats.txt"

echo ""
echo "=== Complete! ==="
echo "Output files:"
echo "  - Phased VCF: ${OUTPUT_PREFIX}.vcf.gz"
echo "  - Statistics: ${OUTPUT_PREFIX}_stats.txt"
echo "  - Haplotagged BAM: ${OUTPUT_PREFIX}_haplotagged.bam"
echo "  - Haplotagged BAM index: ${OUTPUT_PREFIX}_haplotagged.bam.bai"
echo ""
echo "To inspect variants:"
echo "  zcat ${OUTPUT_PREFIX}.vcf.gz | grep -v '^##'"
