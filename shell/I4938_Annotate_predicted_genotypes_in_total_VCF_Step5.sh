#!/bin/bash

ml bcftools/1.9
ml htslib/1.9

# 상위 디렉터리 경로
BASE_DIR="/Node4/Research/Project1/PGI-Marmoset-2022-12/Workspace/jhcha/machine_learning_marmoset/Journal_analysis/analysis/vcf_genotyping"
VCF_DIR="/Node4/Research/Project1/PGI-Marmoset-2022-12/Workspace/cjp81/mpileup_call_gvcf/results/gvcf2individual_vcf"

# 모든 샘플 폴더 반복 처리
for SAMPLE_DIR in "$BASE_DIR"/*_vcf_parsing_output; do
    # 샘플 이름 추출
    SAMPLE_NAME=$(basename "$SAMPLE_DIR" | sed 's/_vcf_parsing_output//')
    FINAL_DIR="$SAMPLE_DIR/final"
    ANNOTATION_FILE="$FINAL_DIR/annotation.txt"

    # annotation.txt 파일 확인
    if [[ ! -f "$ANNOTATION_FILE" ]]; then
        echo "Skipping $SAMPLE_NAME - annotation.txt not found."
        continue
    fi

    echo "Processing $SAMPLE_NAME ..."

    # 1. bgzip + tabix
    bgzip -c "$ANNOTATION_FILE" > "$ANNOTATION_FILE.gz"
    tabix -s 1 -b 2 -e 2 "$ANNOTATION_FILE.gz"

    # 2. 기존 VCF에서 GT 추출
    VCF_FILE="$VCF_DIR/${SAMPLE_NAME}.vcf.gz"
    GT_ORIGINAL_TSV="$FINAL_DIR/${SAMPLE_NAME}_GT_original.tsv"

    if [[ ! -f "$VCF_FILE" ]]; then
        echo "Skipping $SAMPLE_NAME - VCF file not found."
        continue
    fi

    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t[%GT]\n' "$VCF_FILE" > "$GT_ORIGINAL_TSV"

    # 3. GT 업데이트 (awk)
    GT_UPDATED_TSV="$FINAL_DIR/${SAMPLE_NAME}_GT_updated.tsv"
    awk 'NR==FNR{a[$1,$2]=$3; next} ($1,$2) in a {$10=a[$1,$2]":"$11} 1' "$ANNOTATION_FILE" "$GT_ORIGINAL_TSV" > "$GT_UPDATED_TSV"

    # 4. 업데이트된 GT 정보를 반영한 새로운 VCF 생성
    VCF_UPDATED="$FINAL_DIR/${SAMPLE_NAME}_GT_updated.vcf.gz"
    bcftools reheader -s "$GT_UPDATED_TSV" "$VCF_FILE" -o "$VCF_UPDATED"

    echo "Completed: $SAMPLE_NAME"
done

echo "All samples processed successfully!"

