# gblup_utils.R
# Utilities for GBLUP: fit, predict for all, and k-fold CV
#
# Requires: bigsnpr, rrBLUP, stats
#
# Functions:
# - run_gblup(...)               : fit GBLUP and optionally predict for all individuals
# - cv_gblup(...)                : k-fold CV using GBLUP; returns per-fold metrics and predictions

#' Run GBLUP and optionally predict for all individuals (phenotyped + unphenotyped)
#'
#' @param qc_results List with elements:
#'   snp_obj (bigSNP), pheno (data.frame aligned one-row-per-genotype), optional keep_ind (row indices)
#' @param trait_name Trait column name (character)
#' @param id_col_name ID column name (character)
#' @param fixed_effects Formula for fixed effects (default ~1)
#' @param drop_missing If TRUE (default) drops individuals with NA phenotype for fitting.
#' @param predict_all If TRUE (default FALSE) compute GEBVs for all individuals (including those with missing phenotypes).
#' @param tol_diag small ridge added to diagonal for numerical stability (default 1e-6)
#' @return list with elements: gebv (data.frame ID, gebv for requested animals), variances, model_fit, G (full G over keep_ind)
#' @importFrom bigsnpr bed bed_scaleBinom bed_tcrossprodSelf
#' @importFrom rrBLUP mixed.solve
#' @importFrom stats model.matrix
#' @export
run_gblup <- function(qc_results,
                      trait_name,
                      id_col_name,
                      fixed_effects = ~1,
                      drop_missing = TRUE,
                      predict_all = FALSE,
                      tol_diag = 1e-6) {
  if (is.null(qc_results)) stop("qc_results must be provided.")
  if (is.null(qc_results$snp_obj)) stop("qc_results$snp_obj not found.")
  if (is.null(qc_results$pheno)) stop("qc_results$pheno not found.")
  pheno <- qc_results$pheno

  if (!trait_name %in% colnames(pheno)) stop("trait_name not found in pheno.")
  if (!id_col_name %in% colnames(pheno)) stop("id_col_name not found in pheno.")

  # get keep_ind (rows in genotype order). default to all individuals
  n_total <- nrow(qc_results$snp_obj$fam)
  keep_ind <- qc_results$keep_ind
  if (is.null(keep_ind)) keep_ind <- seq_len(n_total)

  # Subset phenotype to keep_ind rows (we assume pheno is aligned to genotype order)
  ph_sub <- pheno[keep_ind, , drop = FALSE]
  # ensure numeric trait
  ph_sub[[trait_name]] <- as.numeric(ph_sub[[trait_name]])
  y_full <- ph_sub[[trait_name]]
  ids_full <- ph_sub[[id_col_name]]

  # Build X design matrix for all keep_ind rows
  X_full <- stats::model.matrix(fixed_effects, data = ph_sub)

  # Count missing
  na_count <- sum(is.na(y_full))

  # Compute G (full over keep_ind)
  bedfile_path <- qc_results$snp_obj$bedfile
  if (is.null(bedfile_path)) stop("snp_obj$bedfile is missing. Re-run load_plink() that sets snp_obj$bedfile.")
  obj.bed <- bigsnpr::bed(bedfile_path)

  message("Computing G matrix (this may take some seconds)...")
  G_full <- bigsnpr::bed_tcrossprodSelf(obj.bed, ind.row = keep_ind, fun.scaling = bigsnpr::bed_scaleBinom)
  diag(G_full) <- diag(G_full) + tol_diag
  message("G computed. Dimensions: ", paste(dim(G_full), collapse = " x "))

  # Decide training set
  if (na_count > 0) {
    if (drop_missing) {
      obs_idx <- which(!is.na(y_full))
      message(sprintf("Dropping %d missing phenotypes; fitting on %d observed individuals.", na_count, length(obs_idx)))
    } else {
      stop("Missing phenotypes present. Set drop_missing = TRUE to fit only on observed individuals.")
    }
  } else {
    obs_idx <- seq_along(y_full)
  }

  # Observed data subsets
  y_obs <- y_full[obs_idx]
  X_obs <- X_full[obs_idx, , drop = FALSE]
  ids_obs <- ids_full[obs_idx]

  # Subset G to observed individuals for fitting
  G_obs <- G_full[obs_idx, obs_idx, drop = FALSE]

  # Fit mixed model using rrBLUP::mixed.solve
  message("Fitting mixed model (rrBLUP::mixed.solve)...")
  fit <- rrBLUP::mixed.solve(y = y_obs, X = X_obs, K = G_obs)
  message("Model fitted.")

  # Extract variance components and heritability
  sigma_g <- as.numeric(fit$Vu)
  sigma_e <- as.numeric(fit$Ve)
  h2 <- sigma_g / (sigma_g + sigma_e)
  variances <- list(sigma_g = sigma_g, sigma_e = sigma_e, h2 = h2)

  # GEBVs for observed individuals (from fit$u)
  u_obs <- as.vector(fit$u)
  gebv_obs_df <- data.frame(ID = ids_obs, gebv = u_obs, stringsAsFactors = FALSE)
  names(gebv_obs_df)[1] <- id_col_name

  out <- list(gebv_obs = gebv_obs_df,
              variances = variances,
              model_fit = fit,
              G = G_full,
              keep_ind = keep_ind,
              ids_full = ids_full)

  # If predict_all requested, compute GEBVs for all individuals (observed + unobserved)
  if (predict_all) {
    message("Computing GEBVs for all individuals (observed + unobserved) ...")

    # Solve system: u_all = G_all_obs %*% inv(G_obs + (Ve/Vu) I) %*% (y_obs - X_obs beta)
    ratio <- sigma_e / sigma_g
    # regularizer matrix is ratio * I
    # Build matrix to invert: (G_obs + ratio * I). For stability, add tiny tol_diag if needed.
    A <- G_obs + ratio * diag(nrow(G_obs))
    # compute right-hand side vector r = y_obs - X_obs %*% beta
    beta_hat <- as.numeric(fit$beta)
    rvec <- as.numeric(y_obs - X_obs %*% beta_hat)

    # Precompute invA %*% rvec (solve instead of invert)
    invA_r <- solve(A, rvec)

    # G_all_obs: rows = all individuals, cols = observed individuals
    G_all_obs <- G_full[, obs_idx, drop = FALSE]    # dimension: n_all x n_obs

    u_all <- as.vector(G_all_obs %*% invA_r)

    # Build dataframe for all individuals (aligned to ids_full)
    gebv_all_df <- data.frame(ID = ids_full, gebv = u_all, stringsAsFactors = FALSE)
    names(gebv_all_df)[1] <- id_col_name

    out$gebv_all <- gebv_all_df
    message("GEBVs for all individuals computed.")
  }

  return(out)
}


#' Cross-validate GBLUP (k-fold)
#'
#' @param qc_results same as run_gblup
#' @param trait_name trait column name
#' @param id_col_name id column name
#' @param k number of folds (default 5)
#' @param fixed_effects fixed effects formula
#' @param seed random seed for fold assignment
#' @param stratify_by optional column name in qc_results$pheno to stratify folds (e.g., family) (NULL default)
#' @param drop_missing if TRUE, individuals with NA trait are removed before fold assignment
#' @param return_fold_preds if TRUE include predictions for each test set row in output
#' @return list containing fold-level metrics, overall metrics, and optionally fold-level predictions
#' @importFrom bigsnpr bed bed_scaleBinom bed_tcrossprodSelf
#' @importFrom rrBLUP mixed.solve
#' @importFrom stats model.matrix cor lm coef
#' @export
cv_gblup <- function(qc_results,
                     trait_name,
                     id_col_name,
                     k = 5,
                     fixed_effects = ~1,
                     seed = 2025,
                     stratify_by = NULL,
                     drop_missing = TRUE,
                     return_fold_preds = TRUE) {
  if (is.null(qc_results)) stop("qc_results must be provided.")
  if (is.null(qc_results$snp_obj)) stop("qc_results$snp_obj not found.")
  if (is.null(qc_results$pheno)) stop("qc_results$pheno not found.")
  pheno <- qc_results$pheno

  if (!trait_name %in% colnames(pheno)) stop("trait_name not found in pheno.")
  if (!id_col_name %in% colnames(pheno)) stop("id_col_name not found in pheno.")

  # keep_ind default
  n_total <- nrow(qc_results$snp_obj$fam)
  keep_ind <- qc_results$keep_ind
  if (is.null(keep_ind)) keep_ind <- seq_len(n_total)

  ph_sub <- pheno[keep_ind, , drop = FALSE]
  ph_sub[[trait_name]] <- as.numeric(ph_sub[[trait_name]])
  ids_all <- ph_sub[[id_col_name]]
  y_all <- ph_sub[[trait_name]]

  # remove rows with NA if requested
  if (drop_missing) {
    valid_idx <- which(!is.na(y_all))
    if (length(valid_idx) < length(y_all)) {
      message(sprintf("Removing %d rows with missing trait before CV.", length(y_all) - length(valid_idx)))
    }
  } else {
    valid_idx <- seq_along(y_all)
    if (any(is.na(y_all))) stop("Missing phenotypes present and drop_missing = FALSE.")
  }

  # Prepare folds
  set.seed(seed)
  n_valid <- length(valid_idx)
  fold_ids <- rep_len(seq_len(k), n_valid)
  fold_ids <- sample(fold_ids, n_valid, replace = FALSE)

  # If stratification requested and available, do stratified assignment
  if (!is.null(stratify_by) && stratify_by %in% colnames(ph_sub)) {
    # simple stratification by splitting within each strata
    strata <- ph_sub[[stratify_by]][valid_idx]
    fold_ids <- integer(n_valid)
    unique_strata <- unique(strata)
    for (s in unique_strata) {
      idx_s <- which(strata == s)
      n_s <- length(idx_s)
      tmp <- rep_len(seq_len(k), n_s)
      fold_ids[idx_s] <- sample(tmp, n_s)
    }
  }

  # Compute full G once (over keep_ind)
  bedfile_path <- qc_results$snp_obj$bedfile
  if (is.null(bedfile_path)) stop("snp_obj$bedfile missing; re-run load_plink() that sets it.")
  obj.bed <- bigsnpr::bed(bedfile_path)
  G_full <- bigsnpr::bed_tcrossprodSelf(obj.bed, ind.row = keep_ind, fun.scaling = bigsnpr::bed_scaleBinom)
  diag(G_full) <- diag(G_full) + 1e-6

  # Precompute X design matrix for all rows (will subset per fold)
  X_full <- stats::model.matrix(fixed_effects, data = ph_sub)

  # Storage for fold metrics and predictions
  fold_results <- vector("list", k)
  all_preds <- list()

  for (fold in seq_len(k)) {
    # indices (in valid_idx space) assigned to this fold are test
    test_positions <- which(fold_ids == fold)
    if (length(test_positions) == 0) {
      warning("Fold ", fold, " has zero test samples; skipping.")
      next
    }
    test_idx_global <- valid_idx[test_positions]   # indices relative to ph_sub / G_full
    train_idx_global <- setdiff(valid_idx, test_idx_global)

    # Build training subsets
    y_train <- y_all[train_idx_global]
    X_train <- X_full[train_idx_global, , drop = FALSE]
    ids_train <- ids_all[train_idx_global]

    # G matrices
    G_train <- G_full[train_idx_global, train_idx_global, drop = FALSE]
    G_test_train <- G_full[test_idx_global, train_idx_global, drop = FALSE]

    # Fit model on training set
    fit_tr <- rrBLUP::mixed.solve(y = y_train, X = X_train, K = G_train)

    # Extract variances and beta
    Vu <- as.numeric(fit_tr$Vu)
    Ve <- as.numeric(fit_tr$Ve)
    beta_hat <- as.numeric(fit_tr$beta)
    # compute residual vector r = y_train - X_train %*% beta
    r_train <- as.numeric(y_train - X_train %*% beta_hat)

    # Solve A = (G_train + (Ve/Vu) I)
    ratio <- Ve / Vu
    A <- G_train + ratio * diag(nrow(G_train))
    # compute invA * r_train
    invA_r <- solve(A, r_train)

    # Predict u for test samples
    u_test <- as.vector(G_test_train %*% invA_r)

    # Predicted phenotype = X_test %*% beta_hat + u_test
    X_test <- X_full[test_idx_global, , drop = FALSE]
    yhat_test <- as.vector(X_test %*% beta_hat + u_test)
    yobs_test <- y_all[test_idx_global]
    ids_test <- ids_all[test_idx_global]

    # Metrics
    cor_val <- if (length(unique(yobs_test)) > 1) cor(yhat_test, yobs_test) else NA_real_
    rmse_val <- sqrt(mean((yhat_test - yobs_test)^2))
    # bias: regression slope of observed ~ predicted
    lm_slope <- tryCatch({
      coef(lm(yobs_test ~ yhat_test))[2]
    }, error = function(e) NA_real_)

    fold_results[[fold]] <- data.frame(
      fold = fold,
      n_train = length(y_train),
      n_test = length(yobs_test),
      cor = cor_val,
      rmse = rmse_val,
      bias_slope = lm_slope,
      Vu = Vu,
      Ve = Ve,
      h2_est = Vu / (Vu + Ve),
      stringsAsFactors = FALSE
    )

    if (return_fold_preds) {
      all_preds[[fold]] <- data.frame(ID = ids_test, yobs = yobs_test, yhat = yhat_test, fold = fold, stringsAsFactors = FALSE)
    }
  } # end fold loop

  # Combine fold results
  fold_df <- do.call(rbind, fold_results)
  overall <- data.frame(
    mean_cor = mean(fold_df$cor, na.rm = TRUE),
    mean_rmse = mean(fold_df$rmse, na.rm = TRUE),
    mean_bias_slope = mean(fold_df$bias_slope, na.rm = TRUE),
    mean_h2 = mean(fold_df$h2_est, na.rm = TRUE)
  )

  out <- list(fold_results = fold_df, overall = overall)
  if (return_fold_preds) out$predictions <- do.call(rbind, all_preds)

  return(out)
}