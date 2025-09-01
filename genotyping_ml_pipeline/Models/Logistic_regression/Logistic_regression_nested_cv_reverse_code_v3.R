#########################################################################################################
# 0.  Machine Learning # Build Marmoset Chimerism Classification Model     ##    Logistic Regression Model
#########################################################################################################
										 
R # R ver. 4.3.3 - execution
 
###########################################################################
# 1. Training data (Load reference gold standard file)
###########################################################################

# Use I4938 individual (Hair sample, not mixed with blood → chimerism-free) to train gold standard model

# Load required library
library(data.table)

transformed_data <- read.table(
  "/Node4/Research/Project1/PGI-Marmoset-2022-12/Workspace/jhcha/machine_learning_marmoset/Journal_analysis/workspace/I4938_data/final/final_input_data.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

transformed_data$Y[transformed_data$Y == "HET"] <- "alt_ref"
transformed_data$Y[transformed_data$Y == "HOM_ALT"] <- "alt_alt"
transformed_data$Y[transformed_data$Y == "HOM_REF"] <- "ref_ref"

transformed_data$abratio <- 1 - transformed_data$abratio 

##################################################################
## 2. Build Logistic Regression Model (5-fold Nested Cross Validation)
##################################################################
 
# Load required packages
library(caret)


set.seed(123)

# Step 1: Split data (Train 80%, Test 20% - maintain class balance)
train_index <- createDataPartition(transformed_data$Y, p = 0.8, list = FALSE) 
train_data <- transformed_data[train_index, ]
test_data <- transformed_data[-train_index, ]

# Step 2: Outer 5-Fold CV (generate 5-folds from training data)
outer_folds <- createFolds(train_data$Y, k = 5, returnTrain = TRUE) 

# Store results
nested_results <- list()

# Outer Loop (5-Fold)
for (i in 1:length(outer_folds)) {
  
  cat("\n===== Outer Fold", i, "=====\n")
  
  # Split train/validation set for outer fold
  train_index_inner <- outer_folds[[i]]
  inner_train_data <- train_data[train_index_inner, ]
  validation_data <- train_data[-train_index_inner, ]
  
  # Separate X, Y using abratio
  x_inner_train <- inner_train_data[, "abratio", drop = FALSE]
  y_inner_train <- factor(inner_train_data$Y)  
  x_validation <- validation_data[, "abratio", drop = FALSE]
  y_validation <- factor(validation_data$Y)  

  # Step 3: Inner 5-Fold CV (using caret)
  train_control <- trainControl(
    method = "cv", 
    number = 5,
    verboseIter = FALSE,
    savePredictions = "final",
    classProbs = TRUE  # required for probability-based prediction
  )

  # Train logistic regression model
  logistic_model <- train(
    Y ~ abratio,       
    data = inner_train_data,  
    method = "multinom",  
    trControl = train_control,  
    trace = FALSE        
  )

  # Evaluate on validation set
  pred_prob <- predict(logistic_model, newdata = x_validation, type = "prob")

  # Convert to class with highest probability
  predictions <- colnames(pred_prob)[max.col(pred_prob, ties.method = "first")]

  # Convert to factor to match actual labels
  predictions <- factor(predictions, levels = levels(y_validation))

  # Calculate accuracy (excluding NA)
  accuracy <- mean(predictions == y_validation, na.rm = TRUE)  

  cat("Validation Accuracy (Fold", i, "):", accuracy, "\n")

  # Store results
  nested_results[[i]] <- list(
    best_model = logistic_model,   
    accuracy = accuracy
  )
}

# Step 4: Evaluate best model on Test Set
best_outer_fold <- which.max(sapply(nested_results, function(x) x$accuracy))
final_model <- nested_results[[best_outer_fold]]$best_model

# Transform test data
test_data$Y <- factor(test_data$Y)  
y_test <- test_data$Y  

#  Step 5:  Predict on test set (Training data set 20%) using best model
test_pred_prob <- predict(final_model, newdata = test_data, type = "prob")

#  Convert to class with highest probability
test_predictions <- colnames(test_pred_prob)[max.col(test_pred_prob, ties.method = "first")]


test_predictions <- factor(test_predictions, levels = levels(y_test))

#  Final test accuracy
test_accuracy <- mean(test_predictions == y_test, na.rm = TRUE)

cat("\n===== Final Test Set Accuracy:", test_accuracy, "=====\n")
#  0.9701427

###############################################
# 3. Calculate Additional Performance Metrics on Test Set(20%)
###############################################

#  Confusion Matrix 
test_conf_matrix <- table(Predicted = test_predictions, Actual = y_test)
cat("\n=== Confusion Matrix (Test Data) ===\n")
print(test_conf_matrix)
#    Actual
# Predicted alt_alt alt_ref ref_ref
#  alt_alt  328254    3180     132
#  alt_ref    2225  486936   12186
#  ref_ref     263   14488  239978



#  Accuracy  
test_accuracy <- sum(diag(test_conf_matrix)) / sum(test_conf_matrix)
cat("\nTest Data Accuracy:", test_accuracy, "\n")
# Test Data Accuracy: 0.9701427 



#  Balanced Accuracy 
test_balanced_accuracy <- mean(diag(prop.table(test_conf_matrix, 1)))
cat("Balanced Accuracy:", test_balanced_accuracy, "\n")
# Balanced Accuracy: 0.967786 


#  Precision, Recall, F1-score  
test_precision <- diag(test_conf_matrix) / rowSums(test_conf_matrix)
test_recall <- diag(test_conf_matrix) / colSums(test_conf_matrix)
test_f1_score <- 2 * (test_precision * test_recall) / (test_precision + test_recall)

cat("\nPrecision (per class):", test_precision, "\n")
# Precision (per class):  0.990011 0.9712554 0.9420914 

cat("Recall (per class):", test_recall, "\n")
# Recall (per class): 0.9924775 0.9649864 0.9511764 


cat("F1-score (per class):", test_f1_score, "\n")
# F1-score (per class): 0.9912427 0.9681108 0.9466121 

#  Error Rate  
test_error_rate <- 1 - test_accuracy
cat("Error Rate:", test_error_rate, "\n")
# Error Rate: 0.02985725 


#   Multiclass AUC  (Test Set)
library(pROC)
test_roc <- multiclass.roc(response = y_test, predictor = as.matrix(test_pred_prob))
cat("Multiclass AUC:", test_roc$auc, "\n")
# Multiclass AUC: 0.9937326



#   Area Under Precision-Recall Curve (AUPRC, Test Set)
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
# AUPRC (per class):  0.9943172 0.9890813 0.9641187


#############################################################################
# 4. External Evaluation using Test Data (1722300M)
#############################################################################
 
 
execution_time <- system.time({

file_path <- "/Node4/Research/Project1/PGI-Marmoset-2022-12/Workspace/jhcha/machine_learning_marmoset/1722300M_gvcf_parsing_output/merged_chr_final_dat/combined_fin.dat.txt"

  test_data_300M <- as.data.frame(fread(file_path))
  colnames(test_data_300M) <- c("CHROM", "POS", "abratio", "Y")

  test_data_300M$id <- paste(test_data_300M$CHROM, test_data_300M$POS, sep = "_")
  test_data_300M <- test_data_300M[, c("id", "abratio", "Y")]
  test_data_300M$Y <- as.factor(test_data_300M$Y)
  test_data_300M[which(test_data_300M$abratio <= 0.2),"Y"] <- "alt_alt"
  test_data_300M[which(test_data_300M$abratio > 0.2 & test_data_300M$abratio < 0.8), "Y"] <- "alt_ref"
  test_data_300M[which(test_data_300M$abratio >= 0.8),"Y"] <- "ref_ref"
    
  test_data_300M$abratio <- 1- test_data_300M$abratio

   
  x_test <- as.matrix(test_data_300M[, "abratio", drop = FALSE])
  y_test <- test_data_300M$Y

  
  test_predictions <- predict(final_model, x_test)
  test_probabilities <- predict(final_model, x_test, type = "prob")

  # Confusion Matrix
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



 
cat("\n[Total Execution Time (seconds)]\n")
print(execution_time)
#  user  system elapsed 
# 129.262  11.737 130.119 


#  === Confusion Matrix (Test Data) ===
#         Actual
# Predicted  alt_alt  alt_ref  ref_ref
#   alt_alt    70586        0        0
#   alt_ref    11403   295934    31801
#   ref_ref      739    60663 10481953

# Test Data Accuracy: 0.9904496 
# Balanced Accuracy: 0.9555942 

# Precision (per class): 1 0.8726064 0.9941762 
# Recall (per class): 0.8532299 0.8298836 0.9969753 
# F1-score (per class): 0.9208031 0.850709 0.9955738 
# Error Rate: 0.009550374 
# Multiclass AUC: 0.9754847 
# AUPRC (per class): 0.9683265 0.891182 0.9999465 

