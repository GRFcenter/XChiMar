##############################################################################################
## 0. Bulid xgboost best model using 5-fold nested cross validation.      ##  2025.04.07   ##
##############################################################################################
 
R # R ver. 4.3.3 - 실행
 
###########################################################################
# 1. Training data (Reference 파일 gold standard 정보 불러오기)
###########################################################################

# I4938 개체 (Hair 정보, blood 가 섞이지 않아 키메리즘이 없는 gold standard csv  파일로 모델 학습하기) 

# 필요한 라이브러리 로드
library(data.table)

# 파일 경로 설정
file_path <- "/Node4/Research/Project1/PGI-Marmoset-2022-12/Workspace/jekim/chimerism/4938-1/I4938_GoldStandard_mod.csv"

# CSV 파일 불러오기
gold_standard_data <- fread(file_path, header = TRUE, sep = ",")

# 데이터 확인
head(gold_standard_data)

# 데이터 변환
transformed_data <- gold_standard_data[, .(CHROM = Chrom, POS = Pos, abratio= AB_Ratio , Y = True_Genotype)]

table(transformed_data$Y)
# HET HOM_ALT HOM_REF 
# 2528227 1665828 1264597


# True_Genotype 값을 annotation 형식에 맞게 변환
transformed_data[Y == "HET", Y := "alt_ref"]
transformed_data[Y == "HOM_ALT", Y := "alt_alt"]
transformed_data[Y == "HOM_REF", Y := "ref_ref"]

# 결과 확인 
transformed_data<- as.data.frame(transformed_data)
head(transformed_data)
#         CHROM   POS abratio 	Y
# 1 NC_025586.1  7872  0.0588     alt_alt
# 2 NC_048383.1  3043  0.2500     alt_ref
# 3 NC_048383.1  3073  0.4286     alt_ref
# 4 NC_048383.1  3088  0.2308     alt_alt
# 5 NC_048383.1  3117  0.2069     alt_ref
# 6 NC_048383.1  3132  0.9111     ref_ref
  
######################################################
# 2. xgboost 모델 구축 (5-fold Nested cross validation method training)
######################################################

# 필요한 패키지 로드
library(caret)
library(xgboost)

set.seed(123)

# Step 1: 데이터 분할 (Training 80%, Test 20% - 클래스 균형 유지)
train_index <- createDataPartition(transformed_data$Y, p = 0.8, list = FALSE) 
train_data <- transformed_data[train_index, ]
test_data <- transformed_data[-train_index, ]

# Step 2: Outer 5-Fold CV (Train Data에서 5-Fold 생성)
outer_folds <- createFolds(train_data$Y, k = 5, returnTrain = TRUE) 

# 결과 저장 리스트
nested_results <- list()

# Outer Loop (5-Fold)
for (i in 1:length(outer_folds)) {
  
  cat("\n===== Outer Fold", i, "=====\n")
  
  # Outer Fold: Train/Validation Set 분할
  train_index_inner <- outer_folds[[i]]
  inner_train_data <- train_data[train_index_inner, ]
  validation_data <- train_data[-train_index_inner, ]
  
  # X, Y 분리 (`abratio` 사용)
  x_inner_train <- as.matrix(inner_train_data[, "abratio", drop = FALSE])
  y_inner_train <- factor(inner_train_data$Y)  
  x_validation <- as.matrix(validation_data[, "abratio", drop = FALSE])
  y_validation <- factor(validation_data$Y)  

  # Step 3: Inner 5-Fold CV for Hyperparameter Tuning (`caret` 사용)
  train_control <- trainControl(
    method = "cv", 
    number = 5,
    verboseIter = FALSE,
    savePredictions = "final",
    classProbs = TRUE  # 🚀 확률 기반 예측을 위해 classProbs 사용
  )

  # XGBoost 모델 파라미터 설정
  xgb_grid <- expand.grid(
    nrounds = 100,         
    eta = 0.1,             
    max_depth = 6,         
    gamma = 0,             
    colsample_bytree = 0.8, 
    min_child_weight = 1,  
    subsample = 0.8        
  )

  # Inner 5-Fold CV를 통해 최적 하이퍼파라미터 찾기 (`caret` 사용)
  xgb_tuned <- train(
    x_inner_train, y_inner_train,
    method = "xgbTree",
    trControl = train_control,
    tuneGrid = xgb_grid,
    metric = "Accuracy"
  )

  # Step 4: 예측값을 확률 기반으로 예측
  pred_prob <- predict(xgb_tuned, x_validation, type = "prob")

  # 가장 높은 확률을 가진 클래스로 변환
  predictions <- colnames(pred_prob)[max.col(pred_prob, ties.method = "first")]

  # Factor 변환하여 실제 y와 비교 가능하게 맞추기
  predictions <- factor(predictions, levels = levels(y_validation))

  # 정확도 계산 (NA 값 제외)
  accuracy <- mean(predictions == y_validation, na.rm = TRUE)  

  cat("Validation Accuracy (Fold", i, "):", accuracy, "\n")

  # 결과 저장
  nested_results[[i]] <- list(
    best_model = xgb_tuned,  # `finalModel`이 아닌 `xgb_tuned` 그대로 저장
    accuracy = accuracy
  )
}

# Step 5: 최적 모델로 Test Set 평가
best_outer_fold <- which.max(sapply(nested_results, function(x) x$accuracy))
final_model <- nested_results[[best_outer_fold]]$best_model

# Test 데이터 변환
test_data$Y <- factor(test_data$Y)  
y_test <- test_data$Y  

# 최적 모델을 사용하여 **Test Set(20%)** 예측 수행
test_pred_prob <- predict(final_model, as.matrix(test_data[, "abratio", drop = FALSE]), type = "prob")

# 가장 높은 확률을 가진 클래스로 변환
test_predictions <- colnames(test_pred_prob)[max.col(test_pred_prob, ties.method = "first")]

# Factor 변환하여 실제 y와 비교 가능하게 맞추기
test_predictions <- factor(test_predictions, levels = levels(y_test))

# Step 6. 최종 테스트 정확도 계산
test_accuracy <- mean(test_predictions == y_test, na.rm = TRUE)

cat("\n===== Final Test Set Accuracy:", test_accuracy, "=====\n")
# ===== Final Test Set Accuracy: 0.9969599  =====


###############################################
# 3. Test Set(20%)에 대한 추가 성능 지표 계산
###############################################

# Confusion Matrix (Test Set 기준)
test_conf_matrix <- table(Predicted = test_predictions, Actual = y_test)
cat("\n=== Confusion Matrix (Test Data) ===\n")
print(test_conf_matrix)
#    Actual
Predicted alt_alt alt_ref ref_ref
  alt_alt  331086    1107       0
  alt_ref    2079  504405       0
  ref_ref       0     133  252919


# Accuracy 계산
test_accuracy <- sum(diag(test_conf_matrix)) / sum(test_conf_matrix)
cat("\nTest Data Accuracy:", test_accuracy, "\n")
# Test Data Accuracy: 0.9969599 


# Balanced Accuracy 계산
test_balanced_accuracy <- mean(diag(prop.table(test_conf_matrix, 1)))
cat("Balanced Accuracy:", test_balanced_accuracy, "\n")
# Balanced Accuracy: 0.9973457 


# Precision, Recall, F1-score 계산
test_precision <- diag(test_conf_matrix) / rowSums(test_conf_matrix)
test_recall <- diag(test_conf_matrix) / colSums(test_conf_matrix)
test_f1_score <- 2 * (test_precision * test_recall) / (test_precision + test_recall)

cat("\nPrecision (per class):", test_precision, "\n")
# Precision (per class): 0.9966676 0.9958952 0.9994744 

cat("Recall (per class):", test_recall, "\n")
# Recall (per class): 0.9937598 0.9975477 1 

cat("F1-score (per class):", test_f1_score, "\n")
# F1-score (per class): 0.9952116 0.9967208 0.9997371 

# Error Rate 계산
test_error_rate <- 1 - test_accuracy
cat("Error Rate:", test_error_rate, "\n")
# Error Rate: 0.003040132 


# Multiclass AUC 계산 (Test Set 기준)
library(pROC)
test_roc <- multiclass.roc(response = y_test, predictor = as.matrix(test_pred_prob))
cat("Multiclass AUC:", test_roc$auc, "\n")
# Multiclass AUC: 0.9996271 


# Area Under Precision-Recall Curve (AUPRC, Test Set 기준)
library(PRROC)
test_auprc <- sapply(1:ncol(test_pred_prob), function(i) {
  PR <- pr.curve(
    scores.class0 = test_pred_prob[, i], 
    weights.class0 = (y_test == levels(y_test)[i]), 
    curve = TRUE
  )
  PR$auc.integral
})
cat("AUPRC (per class):", test_auprc, "\n")
# AUPRC (per class):  0.9990208 0.9987269 0.9999979


#############################################################################
# 4. Test data (1722300M) 으로 외부 데이터에 대해 평가하기
#############################################################################
 
execution_time <- system.time({

 # 데이터 로드
file_path <- "/Node4/Research/Project1/PGI-Marmoset-2022-12/Workspace/jhcha/machine_learning_marmoset/1722300M_gvcf_parsing_output/merged_chr_final_dat/combined_fin.TEST.dat.txt"

  test_data_300M <- as.data.frame(fread(file_path))
  colnames(test_data_300M) <- c("CHROM", "POS", "abratio", "Y")

  test_data_300M$id <- paste(test_data_300M$CHROM, test_data_300M$POS, sep = "_")
  test_data_300M <- test_data_300M[, c("id", "abratio", "Y")]
  test_data_300M$Y <- as.factor(test_data_300M$Y)
  test_data_300M[which(test_data_300M$abratio <= 0.2),"Y"] <- "alt_alt"

  # 예측 준비
  x_test <- as.matrix(test_data_300M[, "abratio", drop = FALSE])
  y_test <- test_data_300M$Y

  # 예측 수행
  test_predictions <- predict(final_model, x_test)
  test_probabilities <- predict(final_model, x_test, type = "prob")

  # 평가 지표
  test_conf_matrix <- table(Predicted = test_predictions, Actual = test_data_300M$Y)
  cat("\n=== Confusion Matrix (Test Data) ===\n")
  print(test_conf_matrix)

  test_accuracy <- sum(diag(test_conf_matrix)) / sum(test_conf_matrix)
  cat("\nTest Data Accuracy:", test_accuracy, "\n")

  test_balanced_accuracy <- mean(diag(prop.table(test_conf_matrix, 1)))
  cat("Balanced Accuracy:", test_balanced_accuracy, "\n")

  test_precision <- diag(test_conf_matrix) / rowSums(test_conf_matrix)
  test_recall <- diag(test_conf_matrix) / colSums(test_conf_matrix)
  test_f1_score <- 2 * (test_precision * test_recall) / (test_precision + test_recall)

  cat("\nPrecision (per class):", test_precision, "\n")
  cat("Recall (per class):", test_recall, "\n")
  cat("F1-score (per class):", test_f1_score, "\n")

  test_error_rate <- 1 - test_accuracy
  cat("Error Rate:", test_error_rate, "\n")

  # AUC
  library(pROC)
  test_roc <- multiclass.roc(response = test_data_300M$Y, predictor = as.matrix(test_probabilities))
  cat("Multiclass AUC:", test_roc$auc, "\n")

  # AUPRC
  library(PRROC)
  test_auprc <- sapply(1:ncol(test_probabilities), function(i) {
    PR <- pr.curve(
      scores.class0 = test_probabilities[, i],
      weights.class0 = (test_data_300M$Y == levels(test_data_300M$Y)[i]),
      curve = TRUE
    )
    PR$auc.integral
  })
  cat("AUPRC (per class):", test_auprc, "\n")

}) 

# === Confusion Matrix (Test Data) ===
#         Actual
# Predicted  alt_alt  alt_ref  ref_ref
#  alt_alt    70114        0        0
#  alt_ref    11872   295456    31801
#  ref_ref      742    61141 10481953

# Test Data Accuracy: 0.9903629 
# Balanced Accuracy: 0.955117 

# Precision (per class): 1 0.8712201 0.9941309 
# Recall (per class): 0.8475244 0.8285431 0.9969753 
# F1-score (per class): 0.9174703 0.8493459 0.9955511 
# Error Rate: 0.009637108 
# Multiclass AUC: 0.9844706 
# AUPRC (per class): 0.9619854 0.8980108 0.9996453 


# 실행 시간 출력
cat("\n[총 실행 시간 (초)]\n")
print(execution_time)
#   user  system elapsed 
# 243.191  15.266  99.613


###################################################################
# 5. Ab-ratio 기준 - "Total region", "0~0.2이하", "0.8이상~1.0", 영역 genotype 비율 보기
###################################################################

#############################################
## 5-1. I4938 Training dataset                                      ##
#############################################

#################################
# Train data in total region SNP   
#################################

table(transformed_data$Y)    # gatk 기반 클래스 별 genotyping 비율
# alt_alt alt_ref ref_ref 
# 1665828 2528227 1264597 


predictions <- predict(final_model, transformed_data)
transformed_data$"predicted_Y" <- predictions 

table(transformed_data$"predicted_Y")
# alt_alt alt_ref ref_ref 
# 1660900 2532410 1265342 
 

#################################
# Train data <= 0.2 ratio site SNP #
#################################

data_0.2 <- transformed_data[which(transformed_data$abratio <= 0.2),] 
 
table(data_0.2$Y)
# alt_alt	 alt_ref   ref_ref 
# 1657574    9717      0

table(data_0.2$predicted_Y)
# alt_alt	 alt_ref   ref_ref 
# 1660900    6391       0 


#################################
# Train data >= 0.8 ratio site SNP # 
#################################

data_0.8 <-  transformed_data[which(transformed_data$abratio >= 0.8),] 
table(data_0.8$Y)
# alt_alt  alt_ref   ref_ref 
#     41    8818 1264597 

table(data_0.8$predicted_Y)
# alt_alt  alt_ref  ref_ref 
#  0        8114  1265342

 
#############################################
## 5-2. 1722300M Test dataset                                     ##
#############################################
 
#################################
# Test data in total region SNP   
#################################


table(test_data_300M$Y)
#     alt_alt   alt_ref   ref_ref 
#    82728   356597 10513754 

test_predictions <- predict(final_model, test_data_300M)
test_data_300M$"predicted_Y" <- test_predictions
  
table(test_data_300M$"predicted_Y")
# alt_alt	  alt_ref   ref_ref 
#   70114   339129  10543836 


#################################
# Test data <= 0.2 ratio site SNP #
#################################

test_data_0.2 <- test_data_300M[which(test_data_300M$abratio <= 0.2),] 

table(test_data_0.2$Y)
# alt_alt   alt_ref ref_ref 
#  73489       0       0 

table(test_data_0.2$predicted_Y)
# alt_alt    alt_ref     ref_ref 
#  70114    3375       0 


#################################
# Test data >= 0.8 ratio site SNP   #
#################################

test_data_0.8 <- test_data_300M[which(test_data_300M$abratio >= 0.8),] 

table(test_data_0.8$Y)
# alt_alt  alt_ref    ref_ref 
#   778    69869 10484125 
 
table(test_data_0.8$predicted_Y)
# alt_alt  alt_ref  ref_ref 
#       0   10936 10543836 
