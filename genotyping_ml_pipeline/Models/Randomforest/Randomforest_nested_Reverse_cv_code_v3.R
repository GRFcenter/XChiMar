###################################################################################
# 0. Machine Learning # Build Chimerism Classification Model in Marmoset     ##    Random forest reverse model
###################################################################################
R # R ver. 4.3.3 - Execution 				
						
LD_PRELOAD=/BiO/Access/home/jhcha/local/gcc-13.3.0/lib64/libstdc++.so.6 R

###########################################################################
# 1. Training data (Load gold standard reference file)
###########################################################################

# Use I4938 individual (Hair sample only, no blood contamination - gold standard csv) to train the model

# install.packages("caret")
# install.packages("randomForest")
# install.packages("e1071")
# install.packages("data.table")

#  Load required libraries
library(data.table)
library(caret)
library(randomForest)
library(e1071)
 
transformed_data <- read.table(
  "/Node4/Research/Project1/PGI-Marmoset-2022-12/Workspace/jhcha/machine_learning_marmoset/Journal_analysis/workspace/I4938_data/final/final_input_data.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

transformed_data$Y[transformed_data$Y == "HET"] <- "alt_ref"
transformed_data$Y[transformed_data$Y == "HOM_ALT"] <- "alt_alt"
transformed_data$Y[transformed_data$Y == "HOM_REF"] <- "ref_ref"

head(transformed_data)
#        CHROM  POS    abratio       Y
#1 NC_048383.1 3043 0.50000000 alt_ref
# 2 NC_048383.1 3073 0.60000000 alt_ref
# 3 NC_048383.1 3088 0.18181818 alt_alt
# 4 NC_048383.1 3117 0.05263158 alt_ref
# 5 NC_048383.1 3132 0.89090909 ref_ref
# 6 NC_048383.1 3134 0.77777778 alt_ref

transformed_data$abratio <- 1 - transformed_data$abratio 

###############################################################
## 2. Build Random forest model (5-fold Nested cross validation method training)
###############################################################
 
library(caret)
library(ranger)

set.seed(123)

# Step 1:  Split data (Training 80%, Test 20% - maintain class balance)
train_index <- createDataPartition(transformed_data$Y, p = 0.8, list = FALSE) 
train_data <- transformed_data[train_index, ]
test_data <- transformed_data[-train_index, ]

# Step 2: Outer 5-Fold CV (generate 5 folds from training data)
outer_folds <- createFolds(train_data$Y, k = 5, returnTrain = TRUE) 


nested_results <- list()

# Outer Loop (5-Fold)
for (i in 1:length(outer_folds)) {
  
  cat("\n===== Outer Fold", i, "=====\n")
  
  # Split  outer fold into Train/Validation sets
  train_index_inner <- outer_folds[[i]]
  inner_train_data <- train_data[train_index_inner, ]
  validation_data <- train_data[-train_index_inner, ]
  
  # Split X and Y 
  x_inner_train <- as.matrix(inner_train_data[, "abratio", drop = FALSE])
  y_inner_train <- factor(inner_train_data$Y)  
  x_validation <- as.matrix(validation_data[, "abratio", drop = FALSE])
  y_validation <- factor(validation_data$Y)  


  # Step 3: Inner 5-Fold CV using caret
  train_control <- trainControl(
    method = "cv", 
    number = 5,
    verboseIter = FALSE,
    savePredictions = "none",
    classProbs = TRUE
  )

  # Random Forest model training (using ranger)*
  rf_model <- train(
    Y ~ abratio,  
    data = inner_train_data,  
    method = "ranger",  #  `randomForest`  `ranger`  
    trControl = train_control,
    tuneGrid = data.frame(mtry = 1, splitrule = "gini", min.node.size = 10)  # `ranger`  
  )


  #  Step 4. Evaluate model on validation set
  predictions <- predict(rf_model, x_validation)

  accuracy <- mean(predictions == y_validation)
  
  cat("Validation Accuracy (Fold", i, "):", accuracy, "\n")

  # Save result
  nested_results[[i]] <- list(
    best_model = rf_model,
    accuracy = accuracy
  )
}

# Step 5: Evaluate final model on test set
best_outer_fold <- which.max(sapply(nested_results, function(x) x$accuracy))
final_model <- nested_results[[best_outer_fold]]$best_model

 
test_data$Y <- factor(test_data$Y)  
y_test <- test_data$Y  

# Predict on **Test Set(20%)**  
test_predictions <- predict(final_model, test_data)

# Calculate final test accuracy
test_accuracy <- mean(test_predictions == y_test)

cat("\n===== Final Test Set Accuracy:", test_accuracy, "=====\n")
# ===== Final Test Set Accuracy: 0.97079 =====


###############################################
# 3. Calculate additional performance metrics on Test Set (20%)
###############################################

#  Confusion Matrix  ##

confusionMatrix(final_model)
# Cross-Validated (5 fold) Confusion Matrix 
# (entries are percentual average cell counts across resamples)
 
#       Reference
# Prediction alt_alt alt_ref ref_ref
#   alt_alt    30.1     0.2     0.0
#   alt_ref     0.3    44.9     1.1
#   ref_ref     0.0     1.3    22.1

# Accuracy (average) : 0.9706

 
## Confusion Matrix   ##
test_conf_matrix <- table(Predicted = test_predictions, Actual = y_test)
cat("\n=== Confusion Matrix (Test Data) ===\n")
print(test_conf_matrix)
#		   Actual
# Predicted alt_alt alt_ref ref_ref
#  alt_alt  327533    2065     114
#  alt_ref    2949  488211   12054
#  ref_ref     260   14328  240128



##  Accuracy   ##
test_accuracy <- sum(diag(test_conf_matrix)) / sum(test_conf_matrix)
cat("\nTest Data Accuracy:", test_accuracy, "\n")
# Test Data Accuracy: 0.97079


## Balanced Accuracy   ##
test_balanced_accuracy <- mean(diag(prop.table(test_conf_matrix, 1)))
cat("Balanced Accuracy:", test_balanced_accuracy, "\n")
# Balanced Accuracy: 0.9687684

##  Precision, Recall, F1-score   ##
test_precision <- diag(test_conf_matrix) / rowSums(test_conf_matrix)
test_recall <- diag(test_conf_matrix) / colSums(test_conf_matrix)
test_f1_score <- 2 * (test_precision * test_recall) / (test_precision + test_recall)

cat("\nPrecision (per class):", test_precision, "\n")
# Precision (per class): 0.9933912 0.9701856 0.9427284 

cat("Recall (per class):", test_recall, "\n")
# Recall (per class):  0.9902976 0.9675131 0.9517709

cat("F1-score (per class):", test_f1_score, "\n")
#  F1-score (per class): 0.991842 0.9688475 0.9472281 

##  Error Rate   ##
test_error_rate <- 1 - test_accuracy
cat("Error Rate:", test_error_rate, "\n")
# Error Rate: 0.02920998  
  

## Multiclass AUC  (Test Set) ##

library(pROC)

class_labels <- levels(y_test)

y_test <- factor(y_test, levels = class_labels)

test_pred_prob <- predict(final_model, test_data, type = "prob")

test_pred_prob <- as.data.frame(test_pred_prob)
colnames(test_pred_prob) <- class_labels   

# Multiclass AUC 
test_roc <- multiclass.roc(response = y_test, predictor = as.matrix(test_pred_prob))
cat("Multiclass AUC:", test_roc$auc, "\n")
# Multiclass AUC:  0.9939828


## Area Under Precision-Recall Curve (AUPRC) 
library(PRROC)

test_auprc <- sapply(class_labels, function(class_name) {
  binary_labels <- as.numeric(y_test == class_name)  
   
  PR <- pr.curve(
    scores.class0 = test_pred_prob[[class_name]],  
    weights.class0 = binary_labels, 
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
# AUPRC (per class):  0.994458 0.9900069 0.9651267


#############################################################################
# 4. External evaluation using test data (1722300M)
#############################################################################
 
execution_time <- system.time({

  # Load data

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

  # Prepare for prediction
  x_test <- as.matrix(test_data_300M[, "abratio", drop = FALSE])
  y_test <- test_data_300M$Y

  # Make predictions
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
 
#  
# === Confusion Matrix (Test Data) ===
#         Actual
# Predicted  alt_alt  alt_ref  ref_ref
#  alt_alt    68191       24        0
#  alt_ref     5284   323769    31545
#  ref_ref       14     1025 10523227

# Test Data Accuracy: 0.9965405
# Balanced Accuracy: 0.9658054

# Precision (per class): 0.9996482 0.8978669 0.9999013
# Recall (per class): 0.9279076 0.9967705 0.9970113
# F1-score (per class): 0.9624428 0.9447372 0.9984542
# Error Rate: 0.003459484
# Multiclass AUC: 0.9992599
# AUPRC (per class): 0.9997897 0.9975316 0.9999985

# execution_time 
#   user    system   elapsed 
#11244.554   268.297  5827.319

 