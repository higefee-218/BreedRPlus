#' Run Bayesian Genomic Prediction using BGLR
#'
#' Fits a Bayesian regression model for genomic prediction using the BGLR package.
#' Supports various Bayesian models including BayesA, BayesB, BayesC, and BL.
#'
#' @param obj A bigSNP object containing genotype data (must have $genotypes and $fam).
#' @param pheno Data frame containing phenotype data with individual IDs.
#' @param trait_name Character string naming the trait column in pheno.
#' @param id_col_name Character string naming the ID column in pheno.
#' @param model_type Bayesian model type. Options: "BayesA", "BayesB", "BayesC", "BL" (Bayesian Lasso). Default "BayesA".
#' @param n_iter Number of MCMC iterations. Default 5000.
#' @param burn_in Number of burn-in iterations. Default 1000.
#'
#' @return A list containing:
#' \item{gebv}{Data frame with individual IDs and genomic estimated breeding values.}
#' \item{snp_effects}{Data frame with SNP effects.}
#' \item{residual_var}{Residual variance estimate.}
#' \item{model_fit}{The full BGLR model fit object.}
#'
#' @importFrom BGLR BGLR
#' @export
run_bglr <- function(obj, pheno, trait_name, id_col_name,
                     model_type = "BayesA",
                     n_iter = 5000, burn_in = 1000) {

  # Extract numeric genotype matrix from FBM.code256
  geno <- obj$genotypes[, ]  # forces numeric
  rownames(geno) <- as.character(obj$fam$sample.ID)
  colnames(geno) <- paste0("SNP", seq_len(ncol(geno)))

  # Align phenotype to genotype IDs
  idx <- match(rownames(geno), pheno[[id_col_name]])
  pheno_aligned <- pheno[idx, ]
  if (!all(rownames(geno) == pheno_aligned[[id_col_name]])) {
    stop("Genotype IDs and phenotype IDs do not match.")
  }

  # Trait vector
  y <- pheno_aligned[[trait_name]]

  # BGLR model
  ETA <- list(list(X = scale(geno), model = model_type))
  fit <- BGLR::BGLR(y = y, ETA = ETA,
                     nIter = n_iter, burnIn = burn_in, verbose = TRUE)

  # Extract SNP effects safely
  snp_effects <- fit$ETA[[1]]$b
  if (is.null(snp_effects) || length(snp_effects) != ncol(geno)) {
    snp_effects <- rep(0, ncol(geno))
  }
  snp_effects_df <- data.frame(SNP = colnames(geno), Effect = snp_effects)

  # Extract GEBVs
  gebv <- fit$yHat
  gebv_df <- data.frame(ID = pheno_aligned[[id_col_name]], GEBV = gebv)

  return(list(gebv = gebv_df,
              snp_effects = snp_effects_df,
              residual_var = fit$varE,
              model_fit = fit))
}