################################################################################### 
# 0.  Machine Learning # Build Chimerism Classification Model in Marmoset      ##    DNN
################################################################################### 

R # R ver. 4.3.3  - execution

###########################################################################
# 1.  Training data (Load gold standard reference file)
###########################################################################
  
# Use I4938 individual (hair sample, not mixed with blood → chimerism-free) to train gold standard model

#  Load required libraries
library(nnet)
library(data.table)
library(pROC)
library(PRROC)
library(caret)

transformed_data <- read.table(
  "/Node4/Research/Project1/PGI-Marmoset-2022-12/Workspace/jhcha/machine_learning_marmoset/Journal_analysis/workspace/I4938_data/final/final_input_data.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# Convert True_Genotype to annotation-compatible format
transformed_data$Y[transformed_data$Y == "HET"] <- "alt_ref"
transformed_data$Y[transformed_data$Y == "HOM_ALT"] <- "alt_alt"
transformed_data$Y[transformed_data$Y == "HOM_REF"] <- "ref_ref"
  
transformed_data$Y <- as.factor(transformed_data$Y)  #   Convert Y to factor


###############################################################
# 2. Build and train MLP model (5-fold Nested cross validation method)
###############################################################

# Split data (Train 80% - Test 20%)
set.seed(123)
train_index <- createDataPartition(transformed_data$Y, p = 0.8, list = FALSE)
train_data <- transformed_data[train_index, ]
test_data <- transformed_data[-train_index, ]

#  Outer 5-Fold CV  
set.seed(123)
outer_folds <- createFolds(train_data$Y, k = 5, returnTrain = TRUE)
nested_results <- list()

#  **Outer Loop (5-Fold)**
for (i in 1:length(outer_folds)) {
  
  cat("\n===== Outer Fold", i, "=====\n")
  
  #   Outer Fold: split Train/Validation set
  train_index_inner <- outer_folds[[i]]
  inner_train_data <- train_data[train_index_inner, ]
  validation_data <- train_data[-train_index_inner, ]

  # Set up Inner 5-Fold CV
  train_control <- trainControl(
    method = "cv", 
    number = 5,
    verboseIter = FALSE,
    savePredictions = "final",
    classProbs = TRUE
  )

  # Train MLP model (using nnet)
  set.seed(123)
  nnet_model <- nnet(
    Y ~ abratio,
    data = inner_train_data,
    size = 5,         #  number of neurons in hidden layer
    decay = 0.01,     # L2 regularization
    maxit = 200       # max iterations
  )

  #  Evaluate on Validation Set
  pred_labels <- predict(nnet_model, validation_data, type = "class")
  accuracy <- mean(pred_labels == validation_data$Y, na.rm = TRUE)
  
  cat("Validation Accuracy (Fold", i, "):", accuracy, "\n")

  # Save result
  nested_results[[i]] <- list(
    best_model = nnet_model,
    accuracy = accuracy
  )
}

#   Select best model (with highest accuracy)
best_outer_fold <- which.max(sapply(nested_results, function(x) x$accuracy))
final_model <- nested_results[[best_outer_fold]]$best_model

#  Evaluate final model on Test Set (20%)
test_pred_labels <- predict(final_model, test_data, type = "class")

cat("\n=== Confusion Matrix (Test Data) ===\n")
test_conf_matrix <- table(Predicted = test_pred_labels, Actual = test_data$Y)
print(test_conf_matrix)
# === Confusion Matrix (Test Data) ===
# Actual
# Predicted alt_alt alt_ref ref_ref
  alt_alt    327670    2211     117
  alt_ref    2809  488293   12785
  ref_ref     263   14100  239394



#  Final test accuracy calculation
test_accuracy <- sum(diag(test_conf_matrix)) / sum(test_conf_matrix)
cat("\n===== Final Test Set Accuracy:", test_accuracy, "=====\n")
# ===== Final Test Set Accuracy: 0.9703165 =====

 

###############################################
# 3. Additional performance metrics on Test Set (20%)
###############################################

 
## Balanced Accuracy  ## 
test_balanced_accuracy <- mean(diag(prop.table(test_conf_matrix, 1)))
cat("Balanced Accuracy:", test_balanced_accuracy, "\n")
# Balanced Accuracy: 0.9684655 



## Precision, Recall, F1-score  ##
 
test_precision <- diag(test_conf_matrix) / rowSums(test_conf_matrix)
test_recall <- diag(test_conf_matrix) / colSums(test_conf_matrix)
test_f1_score <- 2 * (test_precision * test_recall) / (test_precision + test_recall)

cat("\nPrecision (per class):", test_precision, "\n")
# Precision (per class): 0.9929454 0.9690526 0.9433986 

cat("Recall (per class):", test_recall, "\n")
# 0.9907118 0.9676756 0.9488617 

cat("F1-score (per class):", test_f1_score, "\n")
# F1-score (per class): 0.9918273 0.9683636 0.9461222 


##  Error Rate  ##
test_error_rate <- 1 - test_accuracy
cat("Error Rate:", test_error_rate, "\n")
# Error Rate: 0.02968348 


## Multiclass AUC  (One-vs-Rest) ##
library(pROC)

# Predict probabilities (nnet outputs probabilities with type="raw")
test_pred_prob <- predict(final_model, test_data, type = "raw")

# Convert to data.frame
test_pred_prob <- as.data.frame(test_pred_prob)
colnames(test_pred_prob) <- levels(test_data$Y)   

##  Multiclass AUC   ##
test_roc <- multiclass.roc(response = test_data$Y, predictor = test_pred_prob)
cat("Multiclass AUC (Test Data):", test_roc$auc, "\n")
# Multiclass AUC (Test Data): 0.9937702


## AUPRC (Precision-Recall Curve)  ##
library(PRROC)

## AUPRC (Precision-Recall Curve)
test_auprc <- sapply(colnames(test_pred_prob), function(class_name) {
  binary_labels <- as.numeric(test_data$Y == class_name)

  #  NA value remove
  valid_idx <- !is.na(binary_labels) & !is.na(test_pred_prob[[class_name]])
  if (sum(valid_idx) == 0) return(0)   
  
  #   PR Curve  
  PR <- pr.curve(
    scores.class0 = as.numeric(test_pred_prob[[class_name]][valid_idx]),  
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
# AUPRC (per class):   0.9943029 0.9890921 0.9642094 


#############################################################################
# 4. External evaluation with test data (1722300M)
#############################################################################

library(data.table)
library(pROC)
library(PRROC)

execution_time <- system.time({

  # Load test data
  file_path <- "/Node4/Research/Project1/PGI-Marmoset-2022-12/Workspace/jhcha/machine_learning_marmoset/1722300M_gvcf_parsing_output/merged_chr_final_dat/combined_fin.dat.txt"
  test_data <- as.data.frame(fread(file_path))
  colnames(test_data) <- c("CHROM", "POS", "abratio", "Y")

  test_data$id <- paste(test_data$CHROM, test_data$POS, sep = "_")
  test_data <- test_data[, c("id", "abratio", "Y")]
  test_data$Y <- as.factor(test_data$Y)
  test_data[which(test_data$abratio <= 0.2), "Y"] <- "alt_alt"

  # Predict
  test_pred_labels <- predict(final_model, test_data, type = "class")
  test_pred_prob <- predict(final_model, test_data, type = "raw")

  # Confusion matrix
  test_conf_matrix <- table(Predicted = test_pred_labels, Actual = test_data$Y)
  cat("\n=== Confusion Matrix (Test Data) ===\n")
  print(test_conf_matrix)

  # Accuracy
  test_accuracy <- sum(diag(test_conf_matrix)) / sum(test_conf_matrix)
  cat("\nTest Data Accuracy:", test_accuracy, "\n")

  # Balanced Accuracy
  test_balanced_accuracy <- mean(diag(prop.table(test_conf_matrix, 1)))
  cat("Balanced Accuracy:", test_balanced_accuracy, "\n")

  # Precision, Recall, F1
  test_precision <- diag(test_conf_matrix) / rowSums(test_conf_matrix)
  test_recall <- diag(test_conf_matrix) / colSums(test_conf_matrix)
  test_f1_score <- 2 * (test_precision * test_recall) / (test_precision + test_recall)

  cat("\nPrecision (per class):", test_precision, "\n")
  cat("Recall (per class):", test_recall, "\n")
  cat("F1-score (per class):", test_f1_score, "\n")

  # Error Rate
  test_error_rate <- 1 - test_accuracy
  cat("Error Rate:", test_error_rate, "\n")

  # AUC
  test_pred_prob <- as.data.frame(test_pred_prob)
  colnames(test_pred_prob) <- levels(test_data$Y)
  test_data$Y <- factor(test_data$Y, levels = levels(test_data$Y))

  test_roc <- multiclass.roc(response = test_data$Y, predictor = test_pred_prob)
  cat("Multiclass AUC:", test_roc$auc, "\n")

  # AUPRC
  test_auprc <- sapply(colnames(test_pred_prob), function(class_name) {
    binary_labels <- as.numeric(test_data$Y == class_name)
    valid_idx <- !is.na(binary_labels) & !is.na(test_pred_prob[[class_name]])
    if (sum(valid_idx) == 0) return(0)
    PR <- pr.curve(
      scores.class0 = test_pred_prob[[class_name]][valid_idx],
      weights.class0 = binary_labels[valid_idx],
      curve = TRUE
    )
    if (is.nan(PR$auc.integral)) return(0) else return(PR$auc.integral)
  })

  cat("AUPRC (per class):", test_auprc, "\n")

})


# === Confusion Matrix (Test Data) ===
#          Actual
# Predicted  alt_alt  alt_ref  ref_ref
#   alt_alt    70031        0        0
#   alt_ref    11955   295481    31801
#   ref_ref      742    61116 10481953

# Test Data Accuracy: 0.9903576
# Balanced Accuracy: 0.9550499

# Precision (per class): 1 0.8710164 0.9941332
# Recall (per class): 0.8465211 0.8286133 0.9969753
# F1-score (per class): 0.9168821 0.8492859 0.9955522
# Error Rate: 0.009642403
# Multiclass AUC: 0.9854124
# AUPRC (per class): 0.9623818 0.8935943 0.9998566
 

cat("\n[Total Execution Time (seconds)]\n")
print(execution_time)
# user  system elapsed
# 216.366  23.560 229.914


