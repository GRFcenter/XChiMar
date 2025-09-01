
#!/bin/bash

# Input directory
BASE_DIR="/Node4/Research/Project1/PGI-Marmoset-2022-12/Workspace/jhcha/machine_learning_marmoset/Journal_analysis/analysis/vcf_genotyping"

for dir in "$BASE_DIR"/*_vcf_parsing_output; do
  # Directory identify
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

