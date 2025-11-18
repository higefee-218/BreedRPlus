# breedRplus

Genomic Prediction Tools for Animal and Plant Breeding

## Description

`breedRplus` is a comprehensive R package for genomic prediction in breeding programs. It provides tools for:

- **GBLUP**: Genomic Best Linear Unbiased Prediction with cross-validation
  - Implements VanRaden's Method 1 for G-matrix computation
  - Uses REML via `rrBLUP` for variance component estimation
  - Efficient computation using `bigsnpr` package
- **Bayesian Methods**: BGLR-based Bayesian regression (BayesA, BayesB, BayesC, BL)
- **Machine Learning**: Random Forest, XGBoost, Neural Networks (MLP, CNN), and PLS regression
  - Uses `ranger` for Random Forest, `xgboost` for gradient boosting
  - Uses `keras` for neural networks (MLP, CNN)
  - Uses `pls` for Partial Least Squares regression
- **Data Preparation**: PLINK file loading and genotype imputation

## Installation

### From Source

```r
# Install devtools if not already installed
install.packages("devtools")

# Install breedRplus
devtools::install("path/to/breedRplus")
```

### Building the Package

```r
# Build the package
devtools::build()

# Or build and check
devtools::check()
```

## Example Data

The package includes example data files in `inst/extdata/`:
- `500ind.30K.bed`, `500ind.30K.bim`, `500ind.30K.fam` - PLINK format genotype files
- `500_ind_pheno.xlsx` - Phenotype data

See `inst/example_workflow.R` for a complete example using this dataset.

## Quick Start

### 1. Load PLINK Data

```r
library(breedRplus)

# Load PLINK files
data <- load_plink(
  bed_file = "genotypes.bed",
  bim_file = "genotypes.bim",
  fam_file = "genotypes.fam",
  pheno = phenotype_df,
  id_col_name = "ID"
)
```

### 2. Run GBLUP

```r
# Fit GBLUP model
gblup_result <- run_gblup(
  qc_results = data,
  trait_name = "yield",
  id_col_name = "ID",
  predict_all = TRUE
)

# Access variance components and heritability
print(gblup_result$variances)
# $sigma_g: Genetic variance
# $sigma_e: Residual variance  
# $h2: Heritability

# View predicted breeding values
head(gblup_result$gebv_all)

# Cross-validation for accuracy assessment
cv_result <- cv_gblup(
  qc_results = data,
  trait_name = "yield",
  id_col_name = "ID",
  k = 5,
  return_fold_preds = TRUE
)

# Access accuracy metrics
print(cv_result$overall)
# mean_cor: Mean correlation (accuracy metric)
# mean_rmse: Mean root mean squared error
# mean_bias_slope: Mean bias
# mean_h2: Mean heritability

# View per-fold results
print(cv_result$fold_results)
```

### 3. Run Bayesian Methods

```r
# Extract genotype matrix for Bayesian CV (required for run_bglr_cv)
geno_matrix <- data$snp_obj$genotypes[, ]
rownames(geno_matrix) <- as.character(data$snp_obj$fam$sample.ID)

# Fit Bayesian model
bayes_result <- run_bglr(
  obj = data$snp_obj,
  pheno = data$pheno,
  trait_name = "yield",
  id_col_name = "ID",
  model_type = "BayesA",
  n_iter = 5000,
  burn_in = 1000
)

# Cross-validation
bayes_cv <- run_bglr_cv(
  geno = geno_matrix,
  pheno = data$pheno,
  trait_name = "yield",
  id_col_name = "ID",
  model_type = "BayesA",
  n_iter = 5000,
  burn_in = 1000,
  k_folds = 5,
  seed = 123
)

# View accuracy metrics
print(bayes_cv$cv_correlation)
print(bayes_cv$cv_mse)

# Plot predictions vs observed (save to file)
png("bayesian_cv_predictions.png", width = 800, height = 600, res = 100)
par(mar = c(4, 4, 3, 2))
plot(bayes_cv$y_true, bayes_cv$y_pred,
     xlab = "Observed", ylab = "Predicted",
     main = "Bayesian Model Cross-Validation Predictions",
     pch = 19, col = "blue", cex = 0.7)
abline(0, 1, col = "red", lty = 2, lwd = 2)
cor_val <- cor(bayes_cv$y_true, bayes_cv$y_pred, use = "complete.obs")
text(x = min(bayes_cv$y_true, na.rm = TRUE),
     y = max(bayes_cv$y_pred, na.rm = TRUE),
     labels = paste("r =", round(cor_val, 3)),
     pos = 4, col = "darkgreen", font = 2, cex = 1.2)
dev.off()
```

### 4. Run Machine Learning Models

**Note: We recommend using the Python version for ML methods** (especially neural networks) due to:
- Easier TensorFlow installation
- Better performance for neural networks
- Native implementation without R→Python overhead
- More active development and features

**See `ML_PYTHON_GUIDE.md` for a complete Python ML guide, or `README_PYTHON.md` and `breedrplus-py/` for the full Python package documentation.**

**Python Example (Recommended):**

```python
from breedrplus import load_plink, run_ml_model, run_ml_cv
import pandas as pd

# Load data
data = load_plink(
    bed_file="genotypes.bed",
    bim_file="genotypes.bim",
    fam_file="genotypes.fam",
    pheno=phenotype_df,
    id_col_name="ID"
)

# Extract genotype matrix
geno_matrix = data['snp_obj']['genotypes']

# Train ML model (Random Forest example)
ml_result = run_ml_model(
    pheno=data['pheno'],
    geno=geno_matrix,
    trait_name="yield",
    id_col_name="ID",
    model_type="RF",  # Options: "RF", "XGB", "MLP", "CNN", "PLS"
    n_trees=500       # RF-specific parameter
)

# View predictions
print(ml_result['gebv'].head())

# Cross-validation
ml_cv = run_ml_cv(
    pheno=data['pheno'],
    geno=geno_matrix,
    trait_name="yield",
    id_col_name="ID",
    model_type="RF",
    k=5,
    seed=123,
    n_trees=500
)

# View accuracy metrics
print(f"Correlation: {ml_cv['cv_correlation']:.4f}")
print(f"MSE: {ml_cv['cv_mse']:.4f}")
```

**R Example (Also Available):**

```r
# Extract genotype matrix (required for ML functions)
geno_matrix <- data$snp_obj$genotypes[, ]
rownames(geno_matrix) <- as.character(data$snp_obj$fam$sample.ID)

# Train ML model (Random Forest example)
ml_result <- run_ml_model(
  pheno = data$pheno,
  geno = geno_matrix,
  trait_name = "yield",
  id_col_name = "ID",
  model_type = "RF",  # Options: "RF", "XGB", "MLP", "CNN", "PLS"
  n_trees = 500       # RF-specific parameter
)

# View predictions
head(ml_result$gebv)

# Cross-validation
ml_cv <- run_ml_cv(
  pheno = data$pheno,
  geno = geno_matrix,
  trait_name = "yield",
  id_col_name = "ID",
  model_type = "RF",
  k = 5,
  seed = 123,
  n_trees = 500
)

# View accuracy metrics
print(ml_cv$cv_correlation)
print(ml_cv$cv_mse)
```

## Main Functions

- `load_plink()`: Load PLINK files into bigSNP format
- `run_gblup()`: Fit GBLUP model and predict GEBVs
- `cv_gblup()`: Cross-validate GBLUP model
- `run_bglr()`: Fit Bayesian genomic prediction model
- `run_bglr_cv()`: Cross-validate Bayesian model
- `run_ml_model()`: Train machine learning model
- `run_ml_cv()`: Cross-validate machine learning model

## Dependencies

- `bigsnpr`: For handling large-scale genotype data
- `rrBLUP`: For GBLUP implementation
- `BGLR`: For Bayesian methods
- `ranger`: For Random Forest
- `xgboost`: For gradient boosting
- `keras`: For neural networks
- `pls`: For partial least squares regression
- `magrittr`: For pipe operator

## License

GPL-3

## Example Workflow

A complete example workflow is available in `inst/example_workflow.R`. To run it:

```r
# After installing the package
source(system.file("example_workflow.R", package = "breedRplus"))
```

Or access the example data files:

```r
# Get paths to example data
bed_file <- system.file("extdata", "500ind.30K.bed", package = "breedRplus")
bim_file <- system.file("extdata", "500ind.30K.bim", package = "breedRplus")
fam_file <- system.file("extdata", "500ind.30K.fam", package = "breedRplus")
pheno_file <- system.file("extdata", "500_ind_pheno.xlsx", package = "breedRplus")
```

## Documentation

- **User Manual (English)**: See `USER_MANUAL.md` for comprehensive documentation with examples
- **用户手册 (中文)**: 查看 `USER_MANUAL_CN.md` 获取完整的中文文档和示例
- **Function Help**: `?function_name` (e.g., `?run_gblup`)
- **Example Workflow**: `inst/example_workflow.R`

## Author

Fei Ge, thanks to the LMMs! 

