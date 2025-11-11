#' Cross-Validate Bayesian Genomic Prediction using BGLR
#'
#' Performs k-fold cross-validation for Bayesian genomic prediction models.
#'
#' @param geno Numeric genotype matrix with rownames as individual IDs.
#' @param pheno Data frame containing phenotype data with individual IDs.
#' @param trait_name Character string naming the trait column in pheno.
#' @param id_col_name Character string naming the ID column in pheno.
#' @param model_type Bayesian model type. Options: "BayesA", "BayesB", "BayesC", "BL". Default "BayesA".
#' @param n_iter Number of MCMC iterations. Default 5000.
#' @param burn_in Number of burn-in iterations. Default 1000.
#' @param k_folds Number of cross-validation folds. Default 5.
#' @param seed Random seed for reproducibility. Default 123.
#'
#' @return A list containing:
#' \item{y_true}{Observed trait values.}
#' \item{y_pred}{Predicted trait values from cross-validation.}
#' \item{cv_correlation}{Pearson correlation between observed and predicted values.}
#' \item{cv_mse}{Mean squared error.}
#'
#' @importFrom BGLR BGLR
#' @importFrom stats cor
#' @export
run_bglr_cv <- function(geno, pheno, trait_name, id_col_name,
                        model_type = "BayesA",
                        n_iter = 5000, burn_in = 1000,
                        k_folds = 5, seed = 123) {
  
  set.seed(seed)
  
  # Align phenotype and genotype
  pheno <- pheno[match(rownames(geno), pheno[[id_col_name]]), ]
  y_true <- pheno[[trait_name]]
  n <- length(y_true)
  
  # Create folds
  folds <- sample(rep(1:k_folds, length.out = n))
  
  # Store predictions
  y_pred <- rep(NA, n)
  
  for(f in 1:k_folds) {
    message("Running fold ", f, " of ", k_folds)
    
    train_idx <- which(folds != f)
    test_idx  <- which(folds == f)
    
    # Training data
    y_train <- y_true[train_idx]
    geno_train <- geno[train_idx, ]
    
    # Fit BGLR
    ETA <- list(list(X = scale(geno_train), model = model_type))
    fit <- BGLR::BGLR(y = y_train, ETA = ETA,
                       nIter = n_iter, burnIn = burn_in,
                       verbose = FALSE)
    
    # Test data: scale using training mean/sd
    geno_test <- scale(geno[test_idx, ],
                       center = colMeans(geno_train),
                       scale  = apply(geno_train, 2, sd))
    
    # Predict
    y_pred[test_idx] <- geno_test %*% fit$ETA[[1]]$b
  }
  
  # Compute CV metrics
  cv_cor <- stats::cor(y_true, y_pred)
  cv_mse <- mean((y_true - y_pred)^2)
  
  list(
    y_true = y_true,
    y_pred = y_pred,
    cv_correlation = cv_cor,
    cv_mse = cv_mse
  )
}