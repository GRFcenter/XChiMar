# input and output directories
INPUT_DIR="/Node4/Research/Project1/PGI-Marmoset-2022-12/Workspace/cjp81/mpileup_call_gvcf/results/gvcf2individual_vcf"
OUTPUT_BASE_DIR="/Node4/Research/Project1/PGI-Marmoset-2022-12/Workspace/jhcha/machine_learning_marmoset/Journal_analysis/analysis/vcf_genotyping"

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

																    echo "All processing complete."

