# Genotyping ML Pipeline (XChiMar)

The **XChiMar pipeline** consists of five steps to extract features from VCF files, train the XGBoost model, and apply genotype correction to gVCF-derived genotypes.

This pipeline was developed for the GigaScience 2025 submission and supports reproducible machine learning-based genotyping correction in chimeric marmosets.

To deploy this pipeline, place all five scripts and the relevant input files in the `genotyping_ml_pipeline/` directory as shown in the structure below. Each script can be executed independently to follow the full workflow.

---

## ⚙️ Step-by-Step Workflow

### ✅ Step 1: Parse Genotype Features from gVCF
📄 **Script**: [`process_total_Vcf_parsing_Step1.sh`](./process_total_Vcf_parsing_Step1.sh)

Extracts variant features from a standard VCF file generated using `GATK` or parsed with `bcftools query`.

**Output**: `{sample_file}.txt` file per sample, containing the following tab-delimited columns:
```
CHROM	POS	REF	ALT	QUAL	GT	DP	AD
```
where:
- `GT`: genotype
- `DP`: total depth
- `AD`: allele depth

These parsed files serve as machine learning test inputs for downstream genotype correction.

🛠️ **Command:**
```bash
bash process_total_Vcf_parsing_Step1.sh
```

---

### ✅ Step 2: Prepare Output Directory Structure *(Optional)*
📄 **Script**: [`Final_directory_create_code_Step2.sh`](./Final_directory_create_code_Step2.sh)

Creates directory paths to organize future outputs:
- Corrected VCFs
- Performance metrics
- Logs and model outputs

🛠️ **Command:**
```bash
bash Final_directory_create_code_Step2.sh
```

---

### ✅ Step 3: Train XGBoost Model with Nested CV
📄 **Script**: [`Xgboost_nested_cv_code_Step3.R`](./Xgboost_nested_cv_code_Step3.R)

- Uses parsed `.txt` file from I4938 (hair follicle sample)
- Performs 5-fold nested cross-validation
- Evaluates on both internal and independent test data (e.g., 1722300M)

🛠️ **Command:**
```bash
Rscript Xgboost_nested_cv_code_Step3.R
```

📦 **Inside the script (example):**
```r
library(xgboost)
data <- read.table("I4938_parsed.txt", header=TRUE)
model <- xgboost(data = as.matrix(data[, -1]), label = data[,1], nrounds = 100)
```

---

### ✅ Step 4: Predict Genotypes with Trained Model
📄 **Script**: [`VCF_genotyping_Pipeline_code_Step4.R`](./VCF_genotyping_Pipeline_code_Step4.R)

- Calculates allele balance (AB) ratios: `ALT / (REF + ALT)`
- Applies the trained model from Step 3
- Outputs results in:
  - `annotation.txt`: Predicted genotypes
  - `metrics.txt`: Performance metrics (precision, recall, F1-score, etc.)

🛠️ **Command:**
```bash
Rscript VCF_genotyping_Pipeline_code_Step4.R
```

---

### ✅ Step 5: Annotate and Finalize Corrected VCF
📄 **Script**: [`I4938_Annotate_predicted_genotypes_in_total_VCF_Step5.sh`](./I4938_Annotate_predicted_genotypes_in_total_VCF_Step5.sh)

- Updates `GT` field in original gVCF using `annotation.txt`
- Cleans GT formatting by removing colons
- Compresses and indexes final `.vcf.gz` file using `bgzip` and `tabix`

🛠️ **Command:**
```bash
bash I4938_Annotate_predicted_genotypes_in_total_VCF_Step5.sh
```

---

## 📁 File Structure

To replicate this pipeline, arrange your directory as follows:

```plaintext
genotyping_ml_pipeline/
├── process_total_Vcf_parsing_Step1.sh
├── Final_directory_create_code_Step2.sh
├── Xgboost_nested_cv_code_Step3.R
├── VCF_genotyping_Pipeline_code_Step4.R
├── I4938_Annotate_predicted_genotypes_in_total_VCF_Step5.sh
├── I4938_parsed.txt
├── annotation.txt
├── metrics.txt
└── README.md
```

---

## 🔒 License

This software and training pipeline are provided for **academic and non-commercial use only**.  
Commercial use is strictly prohibited without written permission.

Please contact **[grf@pgi.re.kr]** for licensing inquiries.

