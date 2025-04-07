############################################################################### 
# 1. 파싱 한 vcf 파일들 R에서 불러와서 전처리 및 xgboost 로 전처리 하기
############################################################################### 

# 필요한 라이브러리 로드
library(data.table)

# VCF 파싱 출력 파일이 있는 최상위 디렉토리
input_dir <- "/Node4/Research/Project1/PGI-Marmoset-2022-12/Workspace/jhcha/machine_learning_marmoset/Journal_analysis/analysis/vcf_genotyping"

# 모든 하위 디렉토리를 포함하여 *_output.txt 파일 찾기
file_list <- list.files(input_dir, pattern = "_output.txt$", full.names = TRUE, recursive = TRUE)


#######################################################################################
# 2. Abratio 계산 함수 정의
#######################################################################################
  
calculate_abratio <- function(row) {
  # AD 값 처리
  ad_values <- as.numeric(unlist(strsplit(row["AD"], ",")))
  
  # Reference 값 (첫 번째 AD 값)
  ref <- ifelse(length(ad_values) > 0, ad_values[1], 0)
  
  # DP (깊이) 값
  dp <- as.numeric(row["DP"])
  
  # Abratio 계산
  if (dp > 0) {
    return(ref / dp)
  } else {
    return(NA) # DP가 0일 경우 NA 반환
  }
}



################################################################## 
# 3. 모든 파일 파싱 파일들 위 코드 실행 후, Genotyping 자동화 pipeline 코드
################################################################## 

for (file_path in file_list){
  # 파일 불러오기
  data <- fread(file_path, header = TRUE, sep = "\t")
  data <- as.data.frame(data)
  data$predicted_Y <- as.character(NA)


  # Abratio 계산
  data$"abratio" <- apply(data, 1, calculate_abratio)

  # GT 값에서 1/2 제외
  sub_data <- data[!data$GT == "1/2", ]


  # Y 값 설정
  sub_data$Y <- NA
  sub_data[which(sub_data$GT == "0/1"), "Y"] <- "alt_ref"
  sub_data[which(sub_data$GT == "1/1"), "Y"] <- "alt_alt"

  # 필요한 컬럼 선택
  res_data <- sub_data[, c("#CHROM", "POS", "abratio", "Y")]


  # 컬럼명 변경
  colnames(res_data) <- c("CHROM", "POS", "abratio", "Y")
 

  # 테스트 데이터 feature matrix 및 labels 생성
  x_test <- as.matrix(res_data[, "abratio", drop = FALSE])
  y_test <- res_data$Y


  # 테스트 데이터 예측
  test_predictions <- predict(final_model, x_test)  # 클래스 예측값
  test_probabilities <- predict(final_model, x_test, type = "prob") 
	

	# 혼동 행렬 생성

	test_conf_matrix <- table(test_predictions, y_test)

	# "ref_ref" 행 제거 -> vcf 파일은 ref_ref 데이터가 gvcf 와 달리 없기 때문에
	test_conf_matrix <- test_conf_matrix[rownames(test_conf_matrix) != "ref_ref", ]

	# Accuracy 계산
	test_accuracy <- sum(diag(test_conf_matrix)) / sum(test_conf_matrix)
	# cat("\nTest Data Accuracy:", test_accuracy, "\n")
	 
	# Balanced Accuracy 계산
	test_balanced_accuracy <- mean(diag(prop.table(test_conf_matrix, 1)))
	# cat("Balanced Accuracy:", test_balanced_accuracy, "\n")
	 	
	# Precision, Recall, F1-score 계산
	test_precision <- diag(test_conf_matrix) / rowSums(test_conf_matrix)
	test_recall <- diag(test_conf_matrix) / colSums(test_conf_matrix)
	test_f1_score <- 2 * (test_precision * test_recall) / (test_precision + test_recall)

	# cat("\nPrecision (per class):", test_precision, "\n")
 	# cat("Recall (per class):", test_recall, "\n")
 	# cat("F1-score (per class):", test_f1_score, "\n")

 
	####### Genotyping ########
	res_data$"predicted_Y"<- test_predictions 
	 
	
	res_data$id <- paste(res_data$"CHROM", res_data$POS, sep = "_")
	
	data$id <- paste(data$"#CHROM", data$POS, sep = "_")  

	# predicted_Y를 character 타입으로 변환 후 할당
	data$predicted_Y <- as.character(data$predicted_Y) 
	res_data$predicted_Y <- as.character(res_data$predicted_Y)
 

	# ID 매칭하여 predicted_Y 값 업데이트
	data[match(res_data$id, data$id), "predicted_Y"] <- res_data$predicted_Y
  
	data[which(data$GT == "1/2"), "predicted_Y"] <- "alt1_alt2"
	data$predicted_Y <- as.factor(data$predicted_Y)  # 다시 Factor로 변환 (필요 시)

	# 다시 predictd_Y 값을 0/1, 1/1, 1/2 로 변환
	data$predicted_Y <- as.character(data$predicted_Y)
	data[which(data$"predicted_Y"=="alt_alt"), "predicted_Y"] <- "1/1"
	data[which(data$"predicted_Y"=="alt_ref"), "predicted_Y"] <- "0/1"
	data[which(data$"predicted_Y"=="alt1_alt2"), "predicted_Y"] <- "1/2"
	data$predicted_Y <- as.factor(data$predicted_Y) 

	# 위에서 R 로 작업한 data로 실제 vcf 파일 GT Prediction 값을 annotation 을 하기 위해서 아래 작업함.	
 

 	 # final 디렉터리 경로 설정
  	final_dir <- file.path(dirname(file_path), "final")

  	# 파일 경로 지정
  	annotation_file <- file.path(final_dir, "annotation.txt")
	metrics_file <- file.path(final_dir, "metrics.txt")

  	# 파일 저장
  	fwrite(data[, c("#CHROM", "POS", "predicted_Y")], annotation_file, sep = "\t", col.names = FALSE)

 	 writeLines(c(
    		paste("Test Data Accuracy:", test_accuracy),
    		paste("Balanced Accuracy:", test_balanced_accuracy),
    		paste("Precision (per class):", paste(test_precision, collapse = ", ")),
    		paste("Recall (per class):", paste(test_recall, collapse = ", ")),
    		paste("F1-score (per class):", paste(test_f1_score, collapse = ", "))
  	), metrics_file)

 	cat("Saved to:", final_dir, "\n")
}
 