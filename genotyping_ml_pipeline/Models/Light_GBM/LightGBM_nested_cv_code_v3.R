#########################################################################################################
# 0.  Machine Learning # Build Marmoset Chimerism Classification Model       ##    LightGBM model
#########################################################################################################

R # R ver. 4.3.3 - execution
 
###########################################################################
# 1. Training Data (Load Reference Gold Standard File)
###########################################################################

# I4938 individual (Hair sample only, no blood contamination, chimerism-free) used to train gold standard model

# Load required libraries
# install.packages("lightgbm", repos = "https://cran.r-project.org")
library(data.table)
library(lightgbm)
library(caret)
library(pROC)


# Set file path

transformed_data <- read.table(
  "/Node4/Research/Project1/PGI-Marmoset-2022-12/Workspace/jhcha/machine_learning_marmoset/Journal_analysis/workspace/I4938_data/final/final_input_data.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)


#  Convert True_Genotype to annotation-compatible format
transformed_data$Y[transformed_data$Y == "HET"] <- "alt_ref"
transformed_data$Y[transformed_data$Y == "HOM_ALT"] <- "alt_alt"
transformed_data$Y[transformed_data$Y == "HOM_REF"] <- "ref_ref"


transformed_data$Y <- as.factor(transformed_data$Y)

# Check results

head(transformed_data)
         CHROM   POS abratio 	Y
# 1 NC_025586.1  7872  0.0588     alt_alt
# 2 NC_048383.1  3043  0.2500     alt_ref
# 3 NC_048383.1  3073  0.4286     alt_ref
# 4 NC_048383.1  3088  0.2308     alt_alt
# 5 NC_048383.1  3117  0.2069     alt_ref
# 6 NC_048383.1  3132  0.9111     ref_ref


# Define class levels (keep consistent with num_class)
class_labels <- levels(transformed_data$Y)
num_classes <- length(class_labels)

 
######################################################
# 2.  Build LightGBM Model (5-fold Nested Cross Validation)
######################################################

set.seed(123)
#  Split data

train_index <- createDataPartition(transformed_data$Y, p = 0.8, list = FALSE)
train_data <- transformed_data[train_index, ]
test_data <- transformed_data[-train_index, ]
dim(test_data)

#  Set up Outer 5-Fold CV
outer_folds <- createFolds(train_data$Y, k = 5, returnTrain = TRUE)
nested_results <- list()

#  **Outer Loop (5-Fold)**
for (i in 1:length(outer_folds)) {
  
  cat("\n===== Outer Fold", i, "=====\n")

  #  Split Train/Validation Set
  train_index_inner <- outer_folds[[i]]
  inner_train_data <- train_data[train_index_inner, ]
  validation_data <- train_data[-train_index_inner, ]

  # Create LightGBM Datasets
  dtrain <- lgb.Dataset(data = as.matrix(inner_train_data[, "abratio", drop = FALSE]), 
                        label = as.numeric(inner_train_data$Y) - 1)
  dvalid <- lgb.Dataset(data = as.matrix(validation_data[, "abratio", drop = FALSE]), 
                        label = as.numeric(validation_data$Y) - 1)

  # Set LightGBM parameters (consistent num_class)
  params <- list(
    objective = "multiclass",
    metric = "multi_logloss",
    num_class = num_classes, 
    boosting = "gbdt",   #  Use GBDT instead of DART for early stopping
    learning_rate = 0.1,
    num_leaves = 31,     #  More deeper model (2 -> 31)
    max_depth = -1,      #  Auto-depth (-1)
    min_data_in_leaf = 20,
    lambda_l1 = 0.1,
    lambda_l2 = 0.1
  )

  #  Train LightGBM model
  lgb_model <- lgb.train(
    params = params,
    data = dtrain,
    nrounds = 100,
    valids = list(validation = dvalid),
    verbose = 1
  )

  # Evaluate on Validation Set
  pred_prob <- predict(lgb_model, as.matrix(validation_data[, "abratio", drop = FALSE]))
  
  # Convert prediction to class with highest probability
  pred_labels <- apply(matrix(pred_prob, ncol = num_classes, byrow = TRUE), 1, which.max) - 1
  true_labels <- as.numeric(validation_data$Y) - 1

  # Calculate accuracy
  accuracy <- mean(pred_labels == true_labels, na.rm = TRUE)
  cat("Validation Accuracy (Fold", i, "):", accuracy, "\n")

  # Save results
  nested_results[[i]] <- list(
    best_model = lgb_model,
    accuracy = accuracy
  )
}
 

#  Select best model (with highest validation accuracy)
best_outer_fold <- which.max(sapply(nested_results, function(x) x$accuracy))
final_model <- nested_results[[best_outer_fold]]$best_model

#  **Evaluate final model on 20% Test Set**
test_pred_prob <- predict(final_model, as.matrix(test_data[, "abratio", drop = FALSE]))

#  Convert probabilities to class labels
test_pred_labels <- apply(matrix(test_pred_prob, ncol = num_classes, byrow = TRUE), 1, which.max) - 1
test_true_labels <- as.numeric(test_data$Y) - 1

#  Confusion Matrix 
test_conf_matrix <- table(Predicted = test_pred_labels, Actual = test_true_labels)
print(test_conf_matrix)
#   Actual
# Predicted      0      1      2
#        0 178240 270533 135239
#        1  80000 122396  60748
#        2  72502 111675  56309


#  Calculate final test accuracy
test_accuracy <- sum(diag(test_conf_matrix)) / sum(test_conf_matrix)
cat("\n===== Final Test Set Accuracy:", test_accuracy, "=====\n")
# ===== Final Test Set Accuracy: 0.3281824 =====


###############################################
# 3. Additional Performance Metrics for Test Set(20%) 
###############################################


#  Accuracy 
test_accuracy <- sum(diag(test_conf_matrix)) / sum(test_conf_matrix)
cat("\nTest Data Accuracy:", test_accuracy, "\n")
# Test Data Accuracy: 0.3281824 


#  Balanced Accuracy  
test_balanced_accuracy <- mean(diag(prop.table(test_conf_matrix, 1)))
cat("Balanced Accuracy:", test_balanced_accuracy, "\n")
# Balanced Accuracy: 0.3348251

#   Precision, Recall, F1-score  
test_precision <- diag(test_conf_matrix) / rowSums(test_conf_matrix)
test_recall <- diag(test_conf_matrix) / colSums(test_conf_matrix)
test_f1_score <- 2 * (test_precision * test_recall) / (test_precision + test_recall)

cat("\nPrecision (per class):", test_precision, "\n")
# Precision (per class): 0.3051992 0.4651294 0.2341467 

cat("Recall (per class):", test_recall, "\n")
# Recall (per class): 0.5389095 0.2425585 0.2231863 

cat("F1-score (per class):", test_f1_score, "\n")
# F1-score (per class): 0.3897004 0.3188442 0.2285351 

# Error Rate  
test_error_rate <- 1 - test_accuracy
cat("Error Rate:", test_error_rate, "\n")
# Error Rate: 0.6718176 

# Multiclass AUC 
library(pROC)

class_labels <- levels(test_data$Y)  # Real Y class level 
test_pred_prob <- as.data.frame(matrix(test_pred_prob, ncol = length(class_labels), byrow = TRUE))
colnames(test_pred_prob) <- class_labels   

test_data$Y <- factor(test_data$Y, levels = class_labels)

#   Multiclass AUC  
test_roc <- multiclass.roc(response = test_data$Y, predictor = test_pred_prob)
cat("Multiclass AUC:", test_roc$auc, "\n")
# Multiclass AUC: 0.5003747 


#  Area Under Precision-Recall Curve 
library(PRROC)
 
class_labels <- levels(test_data$Y)  
test_pred_prob <- as.data.frame(matrix(unlist(test_pred_prob), ncol = length(class_labels), byrow = TRUE))
colnames(test_pred_prob) <- class_labels   

#  AUPRC (Precision-Recall Curve)
test_auprc <- sapply(class_labels, function(class_name) {
  #  Transform 1/0 vetor
  binary_labels <- as.numeric(test_data$Y == class_name)
  
 
  valid_idx <- !is.na(binary_labels) & !is.na(unlist(test_pred_prob[[class_name]]))
  if (sum(valid_idx) == 0) return(0) 

  # PR Curve  
  PR <- pr.curve(
    scores.class0 = as.numeric(unlist(test_pred_prob[[class_name]])[valid_idx]),   
    weights.class0 = binary_labels[valid_idx],  
    curve = TRUE
  )

  if (is.nan(PR$auc.integral)) {
    return(0)
  } else {
    return(PR$auc.integral)
  }
})

#  AUPRC  
cat("AUPRC (per class):", test_auprc, "\n")
AUPRC (per class):  0.3005513 0.4694609 0.2300914 
  

#############################################################################
# 4. Evaluate Model on External Test Data (1722300M)
#############################################################################

library(data.table)
library(pROC)
library(PRROC)

execution_time <- system.time({
  
file_path <- "/Node4/Research/Project1/PGI-Marmoset-2022-12/Workspace/jhcha/machine_learning_marmoset/1722300M_gvcf_parsing_output/merged_chr_final_dat/combined_fin.dat.txt"
  test_data_300M <- as.data.frame(fread(file_path))
  colnames(test_data_300M) <- c("CHROM", "POS", "abratio", "Y")

  test_data_300M$id <- paste(test_data_300M$CHROM, test_data_300M$POS, sep = "_")
  test_data_300M <- test_data_300M[, c("id", "abratio", "Y")]
  test_data_300M$Y <- as.factor(test_data_300M$Y)
  test_data_300M[which(test_data_300M$abratio <= 0.2), "Y"] <- "alt_alt"

  
  x_test <- as.matrix(test_data_300M[, "abratio", drop = FALSE])
  y_test <- test_data_300M$Y
  num_classes <- length(levels(y_test))

 
  test_probabilities <- predict(final_model, x_test)

  test_pred_labels <- apply(matrix(test_probabilities, ncol = num_classes, byrow = TRUE), 1, which.max)
  true_labels <- as.numeric(y_test)

  # 3. Confusion Matrix
  test_conf_matrix <- table(Predicted = test_pred_labels, Actual = true_labels)
  cat("\n=== Confusion Matrix (Test Data) ===\n")
  print(test_conf_matrix)

  #  4. Accuracy, Balanced Accuracy
  test_accuracy <- sum(diag(test_conf_matrix)) / sum(test_conf_matrix)
  cat("\nTest Data Accuracy:", test_accuracy, "\n")

  test_balanced_accuracy <- mean(diag(prop.table(test_conf_matrix, 1)))
  cat("Balanced Accuracy:", test_balanced_accuracy, "\n")

  # 5. Precision, Recall, F1-score
  test_precision <- diag(test_conf_matrix) / rowSums(test_conf_matrix)
  test_recall <- diag(test_conf_matrix) / colSums(test_conf_matrix)
  test_f1_score <- 2 * (test_precision * test_recall) / (test_precision + test_recall)

  cat("\nPrecision (per class):", test_precision, "\n")
  cat("Recall (per class):", test_recall, "\n")
  cat("F1-score (per class):", test_f1_score, "\n")

  #  6. Error Rate
  test_error_rate <- 1 - test_accuracy
  cat("Error Rate:", test_error_rate, "\n")

  # 7. Multiclass AUC
  class_labels <- levels(y_test)
  test_prob_df <- as.data.frame(matrix(test_probabilities, ncol = num_classes, byrow = TRUE))
  colnames(test_prob_df) <- class_labels
  test_data_300M$Y <- factor(y_test, levels = class_labels)

  test_roc <- multiclass.roc(response = test_data_300M$Y, predictor = test_prob_df)
  cat("Multiclass AUC:", test_roc$auc, "\n")

  # 8. AUPRC (Precision-Recall Curve)
  test_auprc <- sapply(class_labels, function(class_name) {
    binary_labels <- as.numeric(test_data_300M$Y == class_name)
    valid_idx <- !is.na(binary_labels) & !is.na(test_prob_df[[class_name]])
    if (sum(valid_idx) == 0) return(0)
    PR <- pr.curve(
      scores.class0 = as.numeric(test_prob_df[[class_name]][valid_idx]),
      weights.class0 = binary_labels[valid_idx],
      curve = TRUE
    )
    if (is.nan(PR$auc.integral)) return(0) else return(PR$auc.integral)
  })

  cat("AUPRC (per class):", test_auprc, "\n")

})

 

# === Confusion Matrix (Test Data) ===
#          Actual
# Predicted        1        2        3
#        1    79210   341010 10057691
#        2     1986     8581   252345
#        3     1532     7006   203718

# Test Data Accuracy: 0.02661434
# Balanced Accuracy: 0.3333243

# Precision (per class): 0.007559713 0.0326383 0.959775
# Recall (per class): 0.9574751 0.02406358 0.01937633
# F1-score (per class): 0.01500099 0.02770258 0.03798579
# Error Rate: 0.9733857
# Multiclass AUC: 0.4998842
# AUPRC (per class): 0.008068789 0.03248985 0.9594391


# Print execution time
cat("\n[Total Execution Time (seconds)]\n")
print(execution_time)
#  user  system elapsed
# 209.996   8.852 122.20

 