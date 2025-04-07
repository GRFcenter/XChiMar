
#!/bin/bash

# 상위 디렉토리 경로
BASE_DIR="/Node4/Research/Project1/PGI-Marmoset-2022-12/Workspace/jhcha/machine_learning_marmoset/Journal_analysis/analysis/vcf_genotyping"

# _vcf_parsing_output으로 끝나는 디렉토리만 찾아서 반복
for dir in "$BASE_DIR"/*_vcf_parsing_output; do
  # 디렉토리인지 확인
  if [ -d "$dir" ]; then
    final_dir="$dir/final"
    if [ ! -d "$final_dir" ]; then
      mkdir -p "$final_dir"
      echo "Created: $final_dir"
    else
      echo "Already exists: $final_dir"
    fi
  fi
done

