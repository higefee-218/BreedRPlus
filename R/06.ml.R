# --- 1. Train a ML model ---
#' Train a machine learning model for genomic prediction
#'
#' Internal function to train various ML models. Not exported.
#'
#' @param X_train Training genotype matrix
#' @param y_train Training trait values
#' @param model_type Model type: "RF", "XGB", "MLP", "CNN", or "PLS"
#' @param ... Additional arguments for model-specific parameters
#' @return Trained model object
#' @importFrom ranger ranger
#' @importFrom xgboost xgb.DMatrix xgb.train
#' @import keras
#' @importFrom pls plsr
#' @importFrom stats predict
#' @importFrom magrittr %>%
#' @noRd
train_ml_model <- function(X_train, y_train, model_type = "RF", ...) {
  model_type <- toupper(model_type)
  
  if (model_type == "RF") {
    n_trees <- list(...)$n_trees
    if (is.null(n_trees)) n_trees <- 500
    
    train_data <- data.frame(y = y_train, X_train)
    model_fit <- ranger::ranger(
      dependent.variable.name = "y",
      data = train_data,
      num.trees = n_trees,
      importance = 'none'
    )
    
  } else if (model_type == "XGB") {
    n_rounds <- list(...)$n_rounds
    if (is.null(n_rounds)) n_rounds <- 100
    
    dtrain <- xgboost::xgb.DMatrix(data = X_train, label = y_train)
    params <- list(objective = "reg:squarederror", eta = 0.1, max_depth = 6)
    
    model_fit <- xgboost::xgb.train(
      params = params,
      data = dtrain,
      nrounds = n_rounds,
      verbose = 0
    )
    
  } else if (model_type == "MLP") {
    epochs <- list(...)$epochs
    if (is.null(epochs)) epochs <- 50
    
    # Scale X_train
    X_scaled <- scale(X_train)
    scale_params <- attr(X_scaled, "scaled:center")
    scale_sd <- attr(X_scaled, "scaled:scale")
    
    keras::k_clear_session()
    n_features <- ncol(X_scaled)
    model_fit <- keras_model_sequential() %>%
      layer_dense(units = 64, activation = 'relu', input_shape = c(n_features)) %>%
      layer_dropout(rate = 0.3) %>%
      layer_dense(units = 32, activation = 'relu') %>%
      layer_dense(units = 1)
    
    model_fit %>% compile(loss = 'mse', optimizer = optimizer_rmsprop())
    model_fit %>% fit(X_scaled, y_train, epochs = epochs, batch_size = 32,
                      validation_split = 0.1, verbose = 0)
    
    model_fit <- list(model = model_fit, center = scale_params, scale = scale_sd)
    
  } else if (model_type == "CNN") {
    epochs <- list(...)$epochs
    if (is.null(epochs)) epochs <- 50
    
    # Scale X_train
    X_scaled <- scale(X_train)
    scale_params <- attr(X_scaled, "scaled:center")
    scale_sd <- attr(X_scaled, "scaled:scale")
    
    keras::k_clear_session()
    n_features <- ncol(X_scaled)
    X_array <- array_reshape(X_scaled, c(nrow(X_scaled), n_features, 1))
    
    model_fit_nn <- keras_model_sequential() %>%
      layer_conv_1d(filters = 32, kernel_size = 3, activation = 'relu',
                    input_shape = c(n_features, 1)) %>%
      layer_max_pooling_1d(pool_size = 2) %>%
      layer_flatten() %>%
      layer_dense(units = 50, activation = 'relu') %>%
      layer_dense(units = 1)
    
    model_fit_nn %>% compile(loss = 'mse', optimizer = optimizer_rmsprop())
    model_fit_nn %>% fit(X_array, y_train, epochs = epochs, batch_size = 32,
                         validation_split = 0.1, verbose = 0)
    
    model_fit <- list(model = model_fit_nn, center = scale_params, scale = scale_sd)
    
  } else if (model_type == "PLS") {
    n_comp <- list(...)$n_comp
    if (is.null(n_comp)) n_comp <- min(50, ncol(X_train))
    
    model_fit <- pls::plsr(y_train ~ X_train, ncomp = n_comp, validation = "none")
    
  } else {
    stop("model_type must be one of 'RF','XGB','MLP','CNN','PLS'.")
  }
  
  return(model_fit)
}

# --- 2. Predict using a trained ML model ---
#' Predict using a trained ML model
#'
#' Internal function to make predictions from trained ML models. Not exported.
#'
#' @param model_fit Trained model object
#' @param X_new New genotype matrix for prediction
#' @param model_type Model type: "RF", "XGB", "MLP", "CNN", or "PLS"
#' @return Vector of predictions
#' @importFrom xgboost xgb.DMatrix
#' @importFrom stats predict
#' @import keras
#' @noRd
predict_ml_model <- function(model_fit, X_new, model_type = "RF") {
  model_type <- toupper(model_type)
  
  if (model_type %in% c("RF", "XGB", "PLS")) {
    if (model_type == "RF") {
      pred <- stats::predict(model_fit, data.frame(X_new))$predictions
    } else if (model_type == "XGB") {
      dpredict <- xgboost::xgb.DMatrix(data = X_new)
      pred <- stats::predict(model_fit, dpredict)
    } else if (model_type == "PLS") {
      n_comp <- model_fit$ncomp
      pred <- stats::predict(model_fit, newdata = X_new, ncomp = n_comp)[,,1]
    }
    
  } else if (model_type == "MLP") {
    X_scaled <- scale(X_new, center = model_fit$center, scale = model_fit$scale)
    pred <- as.vector(stats::predict(model_fit$model, X_scaled))
    
  } else if (model_type == "CNN") {
    X_scaled <- scale(X_new, center = model_fit$center, scale = model_fit$scale)
    n_features <- ncol(X_scaled)
    X_array <- array_reshape(X_scaled, c(nrow(X_scaled), n_features, 1))
    pred <- as.vector(stats::predict(model_fit$model, X_array))
  }
  
  return(pred)
}

# --- 3. Run ML model and return GEBVs ---
#' Run Machine Learning Model for Genomic Prediction
#'
#' Trains a machine learning model and returns genomic estimated breeding values (GEBVs).
#'
#' @param pheno Data frame with phenotype data
#' @param geno Numeric genotype matrix (rownames = individual IDs)
#' @param trait_name Column name of trait in pheno
#' @param id_col_name Column name of sample IDs in pheno
#' @param model_type One of "RF", "XGB", "MLP", "CNN", or "PLS"
#' @param ... Additional arguments passed to train_ml_model
#'
#' @return List containing:
#' \item{gebv}{Data frame with IDs and predicted GEBVs}
#' \item{model_fit}{The trained model object}
#'
#' @importFrom stats cor
#' @export
run_ml_model <- function(pheno, geno, trait_name, id_col_name, model_type = "RF", ...) {
  
  # Align genotype rows to phenotype IDs
  geno <- geno[match(as.character(pheno[[id_col_name]]), rownames(geno)), , drop = FALSE]
  stopifnot(all(pheno[[id_col_name]] == rownames(geno)))
  
  # Extract trait values (remove NAs for training)
  y <- pheno[[trait_name]]
  valid_idx <- which(!is.na(y))
  
  if (length(valid_idx) == 0) {
    stop("No valid (non-NA) trait values found.")
  }
  
  # Train model on valid observations
  X_train <- geno[valid_idx, , drop = FALSE]
  y_train <- y[valid_idx]
  
  model_fit <- train_ml_model(X_train, y_train, model_type = model_type, ...)
  
  # Predict for all individuals
  y_pred <- predict_ml_model(model_fit, geno, model_type = model_type)
  
  gebv_df <- data.frame(ID = pheno[[id_col_name]], gebv = y_pred)
  
  return(list(gebv = gebv_df, model_fit = model_fit))
}