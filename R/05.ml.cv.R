# --- ML Cross-Validation Wrapper ---
#' Run k-fold Cross-Validation for ML genomic prediction
#'
#' Performs k-fold cross-validation for machine learning genomic prediction models.
#'
#' @param pheno Data frame with phenotype data
#' @param geno Numeric matrix of genotype (rownames = individual IDs)
#' @param trait_name Column name of trait in pheno
#' @param id_col_name Column name of sample IDs in pheno
#' @param model_type One of "RF", "XGB", "MLP", "CNN", or "PLS"
#' @param k Number of CV folds (default 5)
#' @param seed Random seed for reproducibility
#' @param ... Additional arguments passed to train_ml_model
#'
#' @return List with:
#' \item{cv_gebv}{Data frame with predicted GEBVs for all individuals}
#' \item{cv_correlation}{Pearson correlation between observed and predicted values}
#' \item{cv_mse}{Mean squared error}
#'
#' @importFrom stats cor
#' @export
run_ml_cv <- function(pheno, geno, trait_name, id_col_name,
                      model_type = "RF", k = 5, seed = 123, ...) {
  
  set.seed(seed)
  
  # Align genotype rows to phenotype IDs
  geno <- geno[match(as.character(pheno[[id_col_name]]), rownames(geno)), , drop = FALSE]
  stopifnot(all(pheno[[id_col_name]] == rownames(geno)))
  
  # Initialize storage
  y_true <- pheno[[trait_name]]
  y_pred <- rep(NA, length(y_true))
  
  # Create folds (only on non-NA observations)
  non_na_idx <- which(!is.na(y_true))
  if (length(non_na_idx) == 0) {
    stop("No valid (non-NA) trait values found.")
  }
  folds <- sample(rep(1:k, length.out = length(non_na_idx)))
  
  for (i in 1:k) {
    message(sprintf("CV fold %d/%d", i, k))
    
    test_idx <- non_na_idx[which(folds == i)]
    train_idx <- setdiff(non_na_idx, test_idx)
    
    pheno_train <- pheno[train_idx, , drop = FALSE]
    geno_train <- geno[train_idx, , drop = FALSE]
    
    # Train on training set
    res_train <- run_ml_model(
      pheno = pheno_train,
      geno = geno_train,
      trait_name = trait_name,
      id_col_name = id_col_name,
      model_type = model_type,
      ...
    )
    
    # Predict on test set using trained model
    geno_test <- geno[test_idx, , drop = FALSE]
    y_pred[test_idx] <- predict_ml_model(
      res_train$model_fit,
      geno_test,
      model_type = model_type
    )
  }
  
  # Compute CV metrics
  cv_correlation <- stats::cor(y_true, y_pred, use = "complete.obs")
  cv_mse <- mean((y_true - y_pred)^2, na.rm = TRUE)
  
  cv_gebv <- data.frame(ID = pheno[[id_col_name]], gebv = y_pred)
  
  return(list(
    cv_gebv = cv_gebv,
    cv_correlation = cv_correlation,
    cv_mse = cv_mse
  ))
}