#!/bin/bash

ml bcftools/1.9
ml htslib/1.9

# Input directory
BASE_DIR="/Node4/Research/Project1/PGI-Marmoset-2022-12/Workspace/jhcha/machine_learning_marmoset/Journal_analysis/analysis/vcf_genotyping"
VCF_DIR="/Node4/Research/Project1/PGI-Marmoset-2022-12/Workspace/cjp81/mpileup_call_gvcf/results/gvcf2individual_vcf"

for SAMPLE_DIR in "$BASE_DIR"/*_vcf_parsing_output; do
    SAMPLE_NAME=$(basename "$SAMPLE_DIR" | sed 's/_vcf_parsing_output//')
    FINAL_DIR="$SAMPLE_DIR/final"
    ANNOTATION_FILE="$FINAL_DIR/annotation.txt"

    # Check the annotation.txt file
    if [[ ! -f "$ANNOTATION_FILE" ]]; then
        echo "Skipping $SAMPLE_NAME - annotation.txt not found."
        continue
    fi

    echo "Processing $SAMPLE_NAME ..."

    # 1. bgzip + tabix
    bgzip -c "$ANNOTATION_FILE" > "$ANNOTATION_FILE.gz"
    tabix -s 1 -b 2 -e 2 "$ANNOTATION_FILE.gz"

    # Extract GT from existing VCF
    VCF_FILE="$VCF_DIR/${SAMPLE_NAME}.vcf.gz"
    GT_ORIGINAL_TSV="$FINAL_DIR/${SAMPLE_NAME}_GT_original.tsv"

    if [[ ! -f "$VCF_FILE" ]]; then
        echo "Skipping $SAMPLE_NAME - VCF file not found."
        continue
    fi

    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t[%GT]\n' "$VCF_FILE" > "$GT_ORIGINAL_TSV"

    # GT update (merge with awk)
    GT_UPDATED_TSV="$FINAL_DIR/${SAMPLE_NAME}_GT_updated.tsv"
    awk 'NR==FNR{a[$1,$2]=$3; next} ($1,$2) in a {$10=a[$1,$2]":"$11} 1' "$ANNOTATION_FILE" "$GT_ORIGINAL_TSV" > "$GT_UPDATED_TSV"

    # Remove colons (:) and clean up tabs
    awk 'BEGIN { OFS = "\t" }
    {
        if ($NF ~ /:$/) {
            sub(/:$/, "", $NF)
        }
        print
    }' "$GT_UPDATED_TSV" > "${GT_UPDATED_TSV}.tmp" && mv "${GT_UPDATED_TSV}.tmp" "$GT_UPDATED_TSV"

    # (Optional) Reflect in VCF
    VCF_UPDATED="$FINAL_DIR/${SAMPLE_NAME}_GT_updated.vcf.gz"
    bcftools reheader -s "$GT_UPDATED_TSV" "$VCF_FILE" -o "$VCF_UPDATED"

    echo "Completed: $SAMPLE_NAME"
done

echo " All samples processed successfully!"

