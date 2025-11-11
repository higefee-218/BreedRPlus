# Example Workflow for breedRplus Package
# This script demonstrates how to use the package with the example dataset

library(breedRplus)
library(readxl)  # For reading Excel files

# ============================================================================
# Step 1: Load the example data
# ============================================================================

# Get paths to example data files
bed_file <- system.file("extdata", "500ind.30K.bed", package = "breedRplus")
bim_file <- system.file("extdata", "500ind.30K.bim", package = "breedRplus")
fam_file <- system.file("extdata", "500ind.30K.fam", package = "breedRplus")
pheno_file <- system.file("extdata", "500_ind_pheno.xlsx", package = "breedRplus")

# Load phenotype data
pheno <- read_excel(pheno_file)

# Check phenotype structure
head(pheno)
str(pheno)

# ============================================================================
# Step 2: Load PLINK files and align phenotypes
# ============================================================================

# Load PLINK files into bigSNP format
# Note: You'll need to identify the ID column name in your phenotype file
# Replace "ID" with the actual column name containing individual IDs

data <- load_plink(
  bed_file = bed_file,
  bim_file = bim_file,
  fam_file = fam_file,
  pheno = pheno,
  id_col_name = "ID",  # Update this to match your phenotype file
  impute_method = "mode",
  backingfile = tempfile(fileext = ".rds")
)

# Check the loaded data
str(data)
head(data$pheno)

# ============================================================================
# Step 3: Run GBLUP
# ============================================================================

# Identify the trait column name (update this to match your data)
trait_name <- "Trait"  # Update this to your actual trait column name

# Fit GBLUP model
gblup_result <- run_gblup(
  qc_results = data,
  trait_name = trait_name,
  id_col_name = "ID",  # Update if different
  fixed_effects = ~1,
  drop_missing = TRUE,
  predict_all = TRUE
)

# View results
print(gblup_result$variances)
head(gblup_result$gebv_obs)
head(gblup_result$gebv_all)

# Cross-validation
cv_gblup_result <- cv_gblup(
  qc_results = data,
  trait_name = trait_name,
  id_col_name = "ID",
  k = 5,
  seed = 2025
)

print(cv_gblup_result$overall)
print(cv_gblup_result$fold_results)

# ============================================================================
# Step 4: Run Bayesian Methods
# ============================================================================

# Extract genotype matrix for Bayesian methods
geno_matrix <- data$snp_obj$genotypes[, ]
rownames(geno_matrix) <- as.character(data$snp_obj$fam$sample.ID)

# Fit Bayesian model
bayes_result <- run_bglr(
  obj = data$snp_obj,
  pheno = data$pheno,
  trait_name = trait_name,
  id_col_name = "ID",
  model_type = "BayesA",
  n_iter = 5000,
  burn_in = 1000
)

head(bayes_result$gebv)

# Cross-validation
bayes_cv <- run_bglr_cv(
  geno = geno_matrix,
  pheno = data$pheno,
  trait_name = trait_name,
  id_col_name = "ID",
  model_type = "BayesA",
  n_iter = 5000,
  burn_in = 1000,
  k_folds = 5,
  seed = 123
)

print(bayes_cv$cv_correlation)
print(bayes_cv$cv_mse)

# ============================================================================
# Step 5: Run Machine Learning Models
# ============================================================================

# Random Forest
ml_rf <- run_ml_model(
  pheno = data$pheno,
  geno = geno_matrix,
  trait_name = trait_name,
  id_col_name = "ID",
  model_type = "RF",
  n_trees = 500
)

head(ml_rf$gebv)

# Cross-validation
ml_cv_rf <- run_ml_cv(
  pheno = data$pheno,
  geno = geno_matrix,
  trait_name = trait_name,
  id_col_name = "ID",
  model_type = "RF",
  k = 5,
  seed = 123
)

print(ml_cv_rf$cv_correlation)
print(ml_cv_rf$cv_mse)

# XGBoost
ml_cv_xgb <- run_ml_cv(
  pheno = data$pheno,
  geno = geno_matrix,
  trait_name = trait_name,
  id_col_name = "ID",
  model_type = "XGB",
  k = 5,
  seed = 123,
  n_rounds = 100
)

print(ml_cv_xgb$cv_correlation)

# ============================================================================
# Notes:
# ============================================================================
# 1. Update "ID" and "Trait" column names to match your actual data
# 2. Adjust model parameters (n_iter, burn_in, k, etc.) as needed
# 3. For large datasets, consider using predict_all = FALSE in run_gblup()
# 4. Machine learning models may require more computational resources

