#!/bin/bash

# =====================================================
# XChiMar: Chimerism-aware VCF Genotyping Pipeline
# With config.yaml for portability across environments
# =====================================================

# --------- 0. Read config.yaml ------------
CONFIG_FILE=${1:-"config.yaml"}

if [ ! -f "$CONFIG_FILE" ]; then
  echo "Error: config.yaml not found."
  exit 1
fi

# Parse config.yaml to extract paths
eval $(python3 <<EOF
import yaml
with open("$CONFIG_FILE") as f:
    cfg = yaml.safe_load(f)
print(f'INPUT_DIR="{cfg.get("input_dir")}"')
print(f'OUTPUT_BASE_DIR="{cfg.get("output_dir")}"')
print(f'MODEL_PATH="{cfg.get("model_path")}"')
EOF
)

echo "Config Loaded:"
echo "   Input VCF Dir     : $INPUT_DIR"
echo "   Output Dir        : $OUTPUT_BASE_DIR"
echo "   Model Path        : $MODEL_PATH"

ml bcftools/1.9
ml htslib/1.9

# --------- 1. VCF Parsing ------------
echo "Parsing VCF files from: $INPUT_DIR"

# Create output base directory if it doesn't exist
mkdir -p "$OUTPUT_BASE_DIR"

# Loop through all VCF files in the directory
for INPUT_FILE in "$INPUT_DIR"/*.vcf.gz; do
	    # Extract sample name from file (remove path and extension)
	        SAMPLE_NAME=$(basename "$INPUT_FILE" .vcf.gz)
		    
		    # Define output directory and file
		        OUTPUT_DIR="$OUTPUT_BASE_DIR/${SAMPLE_NAME}_vcf_parsing_output"
			    OUTPUT_FILE="$OUTPUT_DIR/${SAMPLE_NAME}_output.txt"
			        
			        # Create output directory if it doesn't exist
				    mkdir -p "$OUTPUT_DIR"

				        # Process the file
					    echo "Processing $INPUT_FILE..."
					        {
							      # Add column names to the output file
							            echo -e "#CHROM\tPOS\tREF\tALT\tQUAL\tGT\tDP\tAD"
								          # Extract data using bcftools query
									        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t[%GT]\t[%DP]\t[%AD]\n' "$INPUT_FILE"
										    } > "$OUTPUT_FILE" 2>> error.log

										        # Check if the output file is empty
											    if [ ! -s "$OUTPUT_FILE" ]; then
												            echo "Failed to process $INPUT_FILE. Check error.log for details."
													            rm -f "$OUTPUT_FILE"
														        else
																        echo "Processed $INPUT_FILE -> $OUTPUT_FILE"
																	    fi
																    done

# --------- 2. Make final directories ------------
echo "Creating 'final' directories.."
for dir in "$OUTPUT_BASE_DIR"/*_vcf_parsing_output; do
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

# --------- 3. Run R Prediction Pipeline ------------
echo "Running R script for genotyping prediction..."

export R_OUTPUT_BASE_DIR="$OUTPUT_BASE_DIR"
export R_MODEL_PATH="$MODEL_PATH"

LD_PRELOAD=/BiO/Access/home/jhcha/local/gcc-13.3.0/lib64/libstdc++.so.6 Rscript - "$R_OUTPUT_BASE_DIR" "$R_MODEL_PATH" <<'EOF'

# Load necessary libraries
library(data.table)
library(xgboost)
 
# get command line args (args[1]: input_dir, args[2]: model_path)
args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
model_path <- args[2]

cat("DEBUG: input_dir =", input_dir, "\n")
cat("DEBUG: model_path =", model_path, "\n")


# Find all *_output.txt files recursively under the input directory
file_list <- list.files(input_dir, pattern = "_output.txt$", full.names = TRUE, recursive = TRUE)

cat("DEBUG: Found", length(file_list), "files\n")

#  
final_model <- readRDS(model_path)

#######################################################################################
# Define function to calculate Abratio
#######################################################################################

calculate_abratio_ref <- function(row) {
  ad_values <- as.numeric(unlist(strsplit(row["AD"], ",")))
  ref <- ifelse(length(ad_values) > 0, ad_values[1], 0)
  dp <- as.numeric(row["DP"])
  if (dp > 0) {
    return(ref / dp)
  } else {
    return(NA)
  }
}

################################################################## 
# Run genotyping automation pipeline after processing all parsing files
################################################################## 

for (file_path in file_list){
  # Load file
  data <- fread(file_path, header = TRUE, sep = "\t")
  data <- as.data.frame(data)
  data$predicted_Y <- as.character(NA)

  # Calculate Abratio alt
  data$"abratio" <- 1 - apply(data, 1, calculate_abratio_ref)

  # Remove GT values equal to 1/2
  sub_data <- data[!data$GT == "1/2", ]

  # Set Y value
  sub_data$Y <- NA
  sub_data[which(sub_data$abratio >= 0.8), "Y"] <- "alt_alt"
  sub_data[which(sub_data$abratio > 0.2 & sub_data$abratio < 0.8), "Y"] <- "alt_ref"
  sub_data[which(sub_data$abratio <= 0.2), "Y"] <- "ref_ref"
  
  # Select required columns
  res_data <- sub_data[, c("#CHROM", "POS", "abratio", "Y")]

  # Rename columns
  colnames(res_data) <- c("CHROM", "POS", "abratio", "Y")

  # Create test feature matrix and labels
  x_test <- as.matrix(sub_data[, "abratio", drop = FALSE])
  y_test <- sub_data$Y

  # Make predictions on test data
  test_predictions <- predict(final_model, x_test)   
  test_probabilities <- predict(final_model, x_test, type = "prob") 

  # Generate confusion matrix
  test_conf_matrix <- table(test_predictions, y_test)

  # Remove "ref_ref" row -> VCF files don't have ref_ref unlike gVCF
  # test_conf_matrix <- test_conf_matrix[rownames(test_conf_matrix) != "ref_ref", ]

  # Accuracy
  test_accuracy <- sum(diag(test_conf_matrix)) / sum(test_conf_matrix)
  
  # Balanced Accuracy  
  test_balanced_accuracy <- mean(diag(prop.table(test_conf_matrix, 1)))
  
  # Precision, Recall, F1-score 
  test_precision <- diag(test_conf_matrix) / rowSums(test_conf_matrix)
  test_recall <- diag(test_conf_matrix) / colSums(test_conf_matrix)
  test_f1_score <- 2 * (test_precision * test_recall) / (test_precision + test_recall)

  ####### Genotyping ########
  res_data$"predicted_Y"<- test_predictions 
  
  res_data$id <- paste(res_data$"CHROM", res_data$POS, sep = "_")
  
  data$id <- paste(data$"#CHROM", data$POS, sep = "_")  

  # Convert predicted_Y to character before assignment
  data$predicted_Y <- as.character(data$predicted_Y) 
  res_data$predicted_Y <- as.character(res_data$predicted_Y)

  # Update predicted_Y values based on ID match
  data[match(res_data$id, data$id), "predicted_Y"] <- res_data$predicted_Y
  
  data[which(data$GT == "1/2"), "predicted_Y"] <- "alt1_alt2"
  data$predicted_Y <- as.factor(data$predicted_Y)  

  # Convert predicted_Y to VCF GT format: 0/1, 1/1, 1/2
  data$predicted_Y <- as.character(data$predicted_Y)
  data[which(data$"predicted_Y"=="alt_alt"), "predicted_Y"] <- "1/1"
  data[which(data$"predicted_Y"=="alt_ref"), "predicted_Y"] <- "0/1"
  data[which(data$"predicted_Y"=="alt1_alt2"), "predicted_Y"] <- "1/2"
  data$predicted_Y <- as.factor(data$predicted_Y) 

  # Below operations annotate the predicted GT values into the actual VCF file through R.
  
  # Set path to final directory
  final_dir <- file.path(dirname(file_path), "final")
  
  
  if (!dir.exists(final_dir)) {
    dir.create(final_dir, recursive = TRUE)
  }

  # Set file paths
  annotation_file <- file.path(final_dir, "annotation.txt")
  metrics_file <- file.path(final_dir, "metrics.txt")

  # Save annotation file
  fwrite(data[, c("#CHROM", "POS", "predicted_Y")], annotation_file, sep = "\t", col.names = FALSE)
  
  # Save metrics file
  writeLines(c(
    paste("Test Data Accuracy:", test_accuracy),
    paste("Balanced Accuracy:", test_balanced_accuracy),
    paste("Precision (per class):", paste(test_precision, collapse = ", ")),
    paste("Recall (per class):", paste(test_recall, collapse = ", ")),
    paste("F1-score (per class):", paste(test_f1_score, collapse = ", "))
  ), metrics_file)

  cat("Saved to:", final_dir, "\n")
}
EOF

# --------- 4. Reflect into VCF ------------
echo "Updating VCF files..."
#!/bin/bash

ml bcftools/1.9
ml htslib/1.9

# Input directory

# Reuse paths read from config.yaml
BASE_DIR="$OUTPUT_BASE_DIR"
VCF_DIR="$INPUT_DIR"

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
