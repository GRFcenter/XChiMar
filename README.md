<p align="center">
  <img src="./XChiMar_logo3.png" alt="XChiMar Logo" width="300"/>
</p>

## XChiMar: Machine Learning-Based Genotype Correction in Chimeric Marmosets

**XChiMar** is a machine learning-based genotyping correction tool specifically designed to address chimerism-induced errors in next-generation sequencing (NGS) data from _Callithrix_ marmosets. It leverages an XGBoost classifier trained on allele balance metrics to improve variant calling accuracy in chimeric, non-blood tissues such as hair follicles, tail, and ear.

## Introduction
In chimeric marmosets, especially in non-blood tissues, standard variant callers tend to misclassify variant genotypes due to DNA contamination from the co-twin’s genome. We trained a robust model using nested cross-validation on 5.4 million gold-standard SNVs and demonstrated that our correction reduces genotyping errors significantly while preserving true heterozygosity.

## Features
- **Optimized for chimeric organisms** : Corrects systematic genotyping biases arising from twin-derived DNA contamination.
- **Machine learning powered** : Trained with high-confidence genotype labels using XGBoost.
- **Input-agnostic** : Works with standard VCF files derived from tools like ‘GATK’ or `bcftools mpileup`.
- **Interpretable predictions ** : Provides insights into false positives and false negatives caused by heterozygous overcalls.
- Validated across 56 chimeric common marmosets and applied to public datasets across _Callithrix jacchus_, _C. geoffroyi_, and _C. kuhlii_. 

## Installation
Linux shell environment
- bash  
- bcftools >= 1.9
- htslib >= 1.9
- awk
- R >= 4.3.3
- Phython3

R environment
- install.packages("data.table")
- install.packages("xgboost")
- install.packages("caret")
- install.packages("pROC")
- install.packages("PRROC")
  
## Usage
- The run_XChiMar.sh code loads the pre-trained xgboost model (XChiMar); xgboost_final_model.rds using hair data from Marmoset I4938 
  samples, and can be used through a simple shell command. 
- It can be run with the following command: bash run_XChiMar.sh config.yml
- At this time, the data and directory paths suitable for the user environment are specified in the config.yml file, considering user 
  convenience and enabling flexible use.
  
## Example
- We provide example_data and code in https://github.com/GRFcenter/XChiMar/.
- Using the data ("./GRFcenter/XChiMar/example_data/input_data/1612151M.vcf.gz"), you can obtain the genotype prediction 
  performance using the XChiMar model and the updated genotype vcf file as a result.
- If you set the path as follows in config.yml, you can easily use it according to your environment as follows.
__________________________________________________________________________________________________
- # config.yaml
# Just change the path to suit your environment.

# 1. Input directory containing VCF files
input_dir: "/your/path/to/input/vcf_dir"
# ex) input_dir: "C:/Users/user/OneDrive/Giga_science/XChiMar/example_data/input_data"

# 2. Output directory to save results (subfolders for each sample are automatically created)
output_dir: "/your/path/to/output/results"
# ex) output_dir: "C:/Users/user/OneDrive/Giga_science/XChiMar/example_data/output_data"

# 3. Trained XGBoost model RDS file path
model_path: "/your/path/to/xgboost_final_model.rds"
# ex) model_path: "C:/Users/user/OneDrive/Giga_science/XChiMar/code/xgboost_final_model.rds"
__________________________________________________________________________________________________
- The actual XChiMar can be easily run as follows: bash run_XChiMar.sh config.yml
  
## LICENSE
This project is licensed under the [MIT License](LICENSE) with a patent notice. See the LICENSE file for details.
For commercial use or questions regarding the pending patent, please contact [grf@pgi.re.kr](mailto:grf@pgi.re.kr).

## Citation
 1. XGBoost: A Scalable Tree Boosting System, https://doi.org/10.1145/2939672.2939785
    
## Contact
pgi@grf.re.kr

## Version

