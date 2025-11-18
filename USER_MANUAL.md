# breedRplus User Manual

**Genomic Prediction Tools for Animal and Plant Breeding**

Version 0.0.0.9000

---

## Table of Contents

1. [Introduction](#introduction)
2. [Installation](#installation)
3. [Quick Start Guide](#quick-start-guide)
4. [Data Preparation](#data-preparation)
5. [GBLUP Methods](#gblup-methods)
6. [Bayesian Methods](#bayesian-methods)
7. [Machine Learning Methods](#machine-learning-methods)
8. [Complete Workflow Examples](#complete-workflow-examples)
9. [Troubleshooting](#troubleshooting)
10. [Best Practices](#best-practices)
11. [References](#references)

---

## Introduction

`breedRplus` is a comprehensive R package designed for genomic prediction in plant and animal breeding programs. It provides a unified interface for multiple genomic prediction methods:

- **GBLUP** (Genomic Best Linear Unbiased Prediction)
- **Bayesian Methods** (BayesA, BayesB, BayesC, Bayesian Lasso)
- **Machine Learning** (Random Forest, XGBoost, Neural Networks, PLS)

### Key Features

- Efficient handling of large-scale genotype data using `bigsnpr`
- Support for PLINK format genotype files
- Cross-validation utilities for model evaluation
- Prediction for both phenotyped and unphenotyped individuals
- Comprehensive documentation and examples

---

## Installation

### Prerequisites

Before installing `breedRplus`, ensure you have R (≥ 4.0) installed. The package requires several dependencies:

- `bigsnpr` - For handling large genotype data
- `rrBLUP` - For GBLUP implementation
- `BGLR` - For Bayesian methods
- `ranger` - For Random Forest
- `xgboost` - For gradient boosting
- `keras` - For neural networks
- `pls` - For partial least squares

### Installing from GitHub

```r
# Install devtools if not already installed
if (!require("devtools")) install.packages("devtools")

# Install breedRplus
devtools::install_github("higefee-218/breedRplus")

# Load the package
library(breedRplus)
```

### Installing from Source

```r
# Download and extract the package
# Then install:
devtools::install("path/to/breedRplus")
```

---

## Quick Start Guide

### Example: Complete Workflow in 5 Steps

```r
library(breedRplus)
library(readxl)  # For reading Excel files

# Step 1: Load your data
data <- load_plink(
  bed_file = "genotypes.bed",
  bim_file = "genotypes.bim",
  fam_file = "genotypes.fam",
  pheno = read_excel("phenotypes.xlsx"),
  id_col_name = "ID"
)

# Step 2: Run GBLUP
gblup_result <- run_gblup(
  qc_results = data,
  trait_name = "Yield",
  id_col_name = "ID",
  predict_all = TRUE
)

# Step 3: View results
print(gblup_result$variances)  # Variance components
head(gblup_result$gebv_all)    # Predicted breeding values

# Step 4: Cross-validation
cv_result <- cv_gblup(
  qc_results = data,
  trait_name = "Yield",
  id_col_name = "ID",
  k = 5
)

# Step 5: Evaluate performance
print(cv_result$overall)  # Overall CV metrics
```

---

## Data Preparation

### Loading PLINK Files

The `load_plink()` function is your entry point for loading genotype data.

#### Function: `load_plink()`

**Purpose**: Load PLINK format genotype files (.bed, .bim, .fam) into R and optionally align phenotype data.

**Syntax**:
```r
load_plink(
  bed_file,
  bim_file,
  fam_file,
  pheno = NULL,
  id_col_name = NULL,
  impute_method = "mode",
  backingfile = tempfile()
)
```

**Parameters**:

| Parameter | Type | Description |
|-----------|------|-------------|
| `bed_file` | character | Path to the .bed file (binary genotype data) |
| `bim_file` | character | Path to the .bim file (variant information) |
| `fam_file` | character | Path to the .fam file (sample information) |
| `pheno` | data.frame | Optional phenotype data frame |
| `id_col_name` | character | Name of ID column in pheno (required if pheno provided) |
| `impute_method` | character | Imputation method: "mode", "mean0", "mean2", "random" (default: "mode") |
| `backingfile` | character | Path for bigSNP backing file (default: temporary file) |

**Returns**: A list containing:
- `snp_obj`: The bigSNP object with genotype data
- `pheno`: Phenotype data aligned to genotype order (if provided)

**Example 1: Load Genotypes Only**

```r
# Load PLINK files without phenotypes
data <- load_plink(
  bed_file = "genotypes.bed",
  bim_file = "genotypes.bim",
  fam_file = "genotypes.fam"
)

# Check the structure
str(data$snp_obj)
head(data$snp_obj$fam)  # Sample information
```

**Example 2: Load Genotypes with Phenotypes**

```r
# Load phenotype data
pheno_df <- read.csv("phenotypes.csv")

# Load PLINK files and align phenotypes
data <- load_plink(
  bed_file = "genotypes.bed",
  bim_file = "genotypes.bim",
  fam_file = "genotypes.fam",
  pheno = pheno_df,
  id_col_name = "AnimalID"
)

# Verify alignment
head(data$pheno)
```

**Example 3: Using Example Data**

```r
# Access example data included with the package
bed_file <- system.file("extdata", "500ind.30K.bed", package = "breedRplus")
bim_file <- system.file("extdata", "500ind.30K.bim", package = "breedRplus")
fam_file <- system.file("extdata", "500ind.30K.fam", package = "breedRplus")
pheno_file <- system.file("extdata", "500_ind_pheno.xlsx", package = "breedRplus")

# Load data
library(readxl)
pheno <- read_excel(pheno_file)

data <- load_plink(
  bed_file = bed_file,
  bim_file = bim_file,
  fam_file = fam_file,
  pheno = pheno,
  id_col_name = "ID"  # Update to match your data
)
```

**Important Notes**:

- The function automatically imputes missing genotypes using the specified method
- Phenotype data is automatically aligned to match genotype sample order
- Only individuals present in both genotype and phenotype data are retained
- The backing file stores genotype data efficiently for large datasets

---

## GBLUP Methods

### Overview

GBLUP (Genomic Best Linear Unbiased Prediction) is a mixed model approach that uses a genomic relationship matrix (G-matrix) to predict breeding values.

### GBLUP Methodology

The `breedRplus` package implements the standard GBLUP method as described by VanRaden (2008). The implementation uses the following components:

#### 1. **Genomic Relationship Matrix (G-matrix) Computation**

The genomic relationship matrix is computed using **VanRaden's Method 1** (VanRaden, 2008):

**Formula**: G = Z Z' / m

Where:
- **Z** is the standardized genotype matrix: Z = (X - 2p) / √[2p(1-p)]
  - X: Genotype matrix (coded as 0, 1, 2 for AA, Aa, aa)
  - p: Allele frequencies (computed from the data)
- **m** is the number of markers (SNPs)

**Implementation Details**:
- Uses `bigsnpr::bed_tcrossprodSelf()` with `bed_scaleBinom` scaling function
- This efficiently computes the G-matrix directly from PLINK .bed files
- A small value (default: 1e-6) is added to the diagonal for numerical stability
- The G-matrix represents the genomic similarity between all pairs of individuals

#### 2. **Mixed Model Framework**

The GBLUP model follows the standard mixed model:

**y = Xβ + Zu + e**

Where:
- **y**: Vector of phenotypic observations
- **X**: Design matrix for fixed effects
- **β**: Vector of fixed effect coefficients
- **Z**: Incidence matrix linking observations to genetic values (identity matrix in GBLUP)
- **u**: Vector of random genetic effects (GEBVs) ~ N(0, Gσ²g)
- **e**: Vector of residuals ~ N(0, Iσ²e)

**Variance Components**:
- **σ²g**: Genetic variance (variance of breeding values)
- **σ²e**: Residual variance (environmental + error variance)
- **Heritability**: h² = σ²g / (σ²g + σ²e)

#### 3. **Model Fitting**

The mixed model is solved using **REML (Restricted Maximum Likelihood)** via `rrBLUP::mixed.solve()`:
- Estimates variance components (σ²g, σ²e) using REML
- Solves Henderson's mixed model equations to obtain:
  - Fixed effect estimates (β)
  - Best Linear Unbiased Predictions (BLUP) of breeding values (u)

#### 4. **Prediction for Unphenotyped Individuals**

When `predict_all = TRUE`, GEBVs are predicted for all individuals (including those without phenotypes) using:

**u_all = G_all,obs × (G_obs + λI)⁻¹ × (y_obs - X_obsβ)**

Where:
- **G_all,obs**: Genomic relationships between all individuals and observed individuals
- **G_obs**: Genomic relationship matrix for observed individuals
- **λ**: Ratio of residual to genetic variance (σ²e / σ²g)
- **y_obs**: Observed phenotypes
- **X_obsβ**: Fixed effects for observed individuals

#### 5. **Cross-Validation**

K-fold cross-validation is performed by:
1. Randomly dividing individuals into k folds
2. For each fold:
   - Training set: All individuals except the fold
   - Test set: Individuals in the fold
   - Fit GBLUP on training set
   - Predict test set individuals using genomic relationships
3. Compute accuracy metrics (correlation, RMSE, bias) for each fold
4. Average metrics across all folds

**Key Features**:
- **Stratified CV**: Optional stratification by groups (e.g., families) ensures each fold has representatives from all groups
- **Unbiased Accuracy**: CV provides unbiased estimates of prediction accuracy

#### Software Implementation

- **G-matrix computation**: `bigsnpr` package (Privé et al., 2018)
  - Efficient handling of large-scale genotype data
  - Direct computation from PLINK .bed files
- **Mixed model solving**: `rrBLUP` package (Endelman, 2011)
  - REML estimation of variance components
  - BLUP computation for breeding values

#### References

- **VanRaden, P. M.** (2008). Efficient methods to compute genomic predictions. *Journal of dairy science*, 91(11), 4414-4423.
- **Endelman, J. B.** (2011). Ridge regression and other kernels for genomic selection with R package rrBLUP. *The Plant Genome*, 4(3), 250-255.
- **Privé, F., et al.** (2018). Efficient analysis of large-scale genome-wide data with two R packages: bigstatsr and bigsnpr. *Bioinformatics*, 34(16), 2781-2787.

### Function: `run_gblup()`

**Purpose**: Fit a GBLUP model and predict genomic estimated breeding values (GEBVs).

**Syntax**:

```r
run_gblup(
  qc_results,
  trait_name,
  id_col_name,
  fixed_effects = ~1,
  drop_missing = TRUE,
  predict_all = FALSE,
  tol_diag = 1e-6
)
```

**Parameters**:

| Parameter | Type | Description |
|-----------|------|-------------|
| `qc_results` | list | Output from `load_plink()` |
| `trait_name` | character | Name of trait column in pheno |
| `id_col_name` | character | Name of ID column in pheno |
| `fixed_effects` | formula | Fixed effects formula (default: `~1` for intercept only) |
| `drop_missing` | logical | Drop individuals with missing phenotypes? (default: TRUE) |
| `predict_all` | logical | Predict for all individuals including unphenotyped? (default: FALSE) |
| `tol_diag` | numeric | Small value added to G-matrix diagonal for stability (default: 1e-6) |

**Returns**: A list containing:
- `gebv_obs`: Data frame with GEBVs for observed individuals (ID, gebv)
- `gebv_all`: Data frame with GEBVs for all individuals (if `predict_all=TRUE`)
- `variances`: List with variance components:
  - `sigma_g`: Genetic variance (σ²g)
  - `sigma_e`: Residual/Environmental variance (σ²e)
  - `h2`: Heritability (h² = σ²g / (σ²g + σ²e))
- `model_fit`: The full mixed model fit object from `rrBLUP::mixed.solve()`
- `G`: The genomic relationship matrix
- `keep_ind`: Indices of individuals used
- `ids_full`: Full list of individual IDs

**Example 1: Basic GBLUP**

```r
# Run GBLUP with intercept only
gblup_result <- run_gblup(
  qc_results = data,
  trait_name = "Yield",
  id_col_name = "ID"
)

# View variance components and heritability
print(gblup_result$variances)
# $sigma_g
# [1] 0.4523
# $sigma_e
# [1] 0.6789
# $h2
# [1] 0.3998

# Access individual components
gblup_result$variances$sigma_g  # Genetic variance
gblup_result$variances$sigma_e  # Residual variance
gblup_result$variances$h2        # Heritability (proportion of variance due to genetics)

# View GEBVs for observed individuals
head(gblup_result$gebv_obs)
```

**Example 2: GBLUP with Fixed Effects**

```r
# Run GBLUP with fixed effects (e.g., batch, sex)
gblup_result <- run_gblup(
  qc_results = data,
  trait_name = "Yield",
  id_col_name = "ID",
  fixed_effects = ~ Batch + Sex
)

# View results
print(gblup_result$variances)
```

**Example 3: Predict for All Individuals**

```r
# Predict GEBVs for both phenotyped and unphenotyped individuals
gblup_result <- run_gblup(
  qc_results = data,
  trait_name = "Yield",
  id_col_name = "ID",
  predict_all = TRUE  # Important!
)

# GEBVs for observed individuals
head(gblup_result$gebv_obs)

# GEBVs for ALL individuals (including unphenotyped)
head(gblup_result$gebv_all)
```

### Function: `cv_gblup()`

**Purpose**: Perform k-fold cross-validation for GBLUP models.

**Syntax**:
```r
cv_gblup(
  qc_results,
  trait_name,
  id_col_name,
  k = 5,
  fixed_effects = ~1,
  seed = 2025,
  stratify_by = NULL,
  drop_missing = TRUE,
  return_fold_preds = TRUE
)
```

**Parameters**:

| Parameter | Type | Description |
|-----------|------|-------------|
| `qc_results` | list | Output from `load_plink()` |
| `trait_name` | character | Name of trait column |
| `id_col_name` | character | Name of ID column |
| `k` | integer | Number of CV folds (default: 5) |
| `fixed_effects` | formula | Fixed effects formula |
| `seed` | integer | Random seed for reproducibility |
| `stratify_by` | character | Optional column name for stratified CV (e.g., "Family") |
| `drop_missing` | logical | Remove individuals with missing trait? |
| `return_fold_preds` | logical | Include fold-level predictions? |

**Returns**: A list containing:
- `fold_results`: Data frame with metrics for each fold:
  - `fold`: Fold number
  - `n_train`: Number of training samples
  - `n_test`: Number of test samples
  - `cor`: Correlation between predicted and observed (accuracy metric)
  - `rmse`: Root mean squared error
  - `bias_slope`: Bias slope (regression of observed on predicted)
  - `Vu`: Genetic variance estimate
  - `Ve`: Residual variance estimate
  - `h2_est`: Heritability estimate
- `overall`: Data frame with overall CV metrics (averaged across folds):
  - `mean_cor`: Mean correlation (primary accuracy metric)
  - `mean_rmse`: Mean root mean squared error
  - `mean_bias_slope`: Mean bias slope
  - `mean_h2`: Mean heritability estimate
- `predictions`: Data frame with individual predictions (if `return_fold_preds=TRUE`):
  - `ID`: Individual ID
  - `yobs`: Observed trait value
  - `yhat`: Predicted trait value
  - `fold`: Fold number

**Example: 5-Fold Cross-Validation**

```r
# Run 5-fold cross-validation
cv_result <- cv_gblup(
  qc_results = data,
  trait_name = "Yield",
  id_col_name = "ID",
  k = 5,
  seed = 2025,
  return_fold_preds = TRUE
)

# View overall accuracy metrics
print(cv_result$overall)
#   mean_cor mean_rmse mean_bias_slope mean_h2
# 1   0.6523    1.2345           0.98   0.3998

# Access individual accuracy metrics
accuracy <- cv_result$overall$mean_cor      # Correlation (accuracy)
rmse <- cv_result$overall$mean_rmse         # Root mean squared error
bias <- cv_result$overall$mean_bias_slope   # Bias slope
h2_cv <- cv_result$overall$mean_h2          # Mean heritability

cat("GBLUP Accuracy (Correlation):", round(accuracy, 4), "\n")
cat("RMSE:", round(rmse, 4), "\n")
cat("Bias Slope:", round(bias, 4), "\n")
cat("Heritability:", round(h2_cv, 4), "\n")

# View fold-specific results
print(cv_result$fold_results)
#   fold n_train n_test      cor     rmse bias_slope        Vu        Ve     h2_est
# 1    1     400    100   0.6456   1.2456      0.975  0.4523  0.6789   0.3998
# 2    2     400    100   0.6589   1.2234      0.985  0.4489  0.6821   0.3971
# 3    3     400    100   0.6512   1.2345      0.978  0.4501  0.6798   0.3985
# 4    4     400    100   0.6498   1.2412      0.982  0.4512  0.6805   0.3989
# 5    5     400    100   0.6556   1.2278      0.980  0.4498  0.6812   0.3978

# View individual predictions
head(cv_result$predictions)
#      ID     yobs      yhat fold
# 1 IND001   12.34   11.89    1
# 2 IND002   15.67   15.23    1
# ...
```

**Interpreting Accuracy Metrics:**

- **`mean_cor` (Correlation)**: Primary accuracy metric
  - Range: -1 to 1, higher is better
  - Typical values: 0.3-0.9 for genomic prediction
  - Represents the correlation between predicted and observed values
  - This is the standard accuracy metric in genomic selection

- **`mean_rmse` (Root Mean Squared Error)**: Prediction error
  - Lower is better
  - Units match your trait units
  - Measures average prediction error magnitude

- **`mean_bias_slope`**: Bias in predictions
  - Ideal value: 1.0 (no bias)
  - < 1.0: Over-prediction (predictions too high)
  - > 1.0: Under-prediction (predictions too low)

- **`mean_h2`**: Mean heritability estimate across folds
  - Range: 0 to 1
  - Proportion of phenotypic variance due to genetics

**Example: Stratified Cross-Validation**

```r
# Stratify by family to ensure each fold has representatives from all families
cv_result <- cv_gblup(
  qc_results = data,
  trait_name = "Yield",
  id_col_name = "ID",
  k = 5,
  stratify_by = "Family",  # Column name in pheno
  seed = 2025
)

# Access accuracy
print(cv_result$overall$mean_cor)
```

### Assessing GBLUP Accuracy and Variance Components

#### Accessing Variance Components from `run_gblup()`

After running `run_gblup()`, you can access variance components and heritability:

```r
# Run GBLUP
gblup_result <- run_gblup(
  qc_results = data,
  trait_name = "Yield",
  id_col_name = "ID",
  predict_all = TRUE
)

# Method 1: Print all variance components
print(gblup_result$variances)
# $sigma_g
# [1] 0.4523
# $sigma_e
# [1] 0.6789
# $h2
# [1] 0.3998

# Method 2: Access individual components
sigma_g <- gblup_result$variances$sigma_g  # Genetic variance
sigma_e <- gblup_result$variances$sigma_e  # Residual variance
h2 <- gblup_result$variances$h2            # Heritability

# Method 3: Format nicely
cat("Genetic Variance (σ²g):", round(sigma_g, 4), "\n")
cat("Residual Variance (σ²e):", round(sigma_e, 4), "\n")
cat("Heritability (h²):", round(h2, 4), "\n")
```

**Understanding Variance Components:**

- **`sigma_g` (σ²g)**: Genetic variance - variance due to genetic differences among individuals
- **`sigma_e` (σ²e)**: Residual variance - variance due to environmental effects and measurement error
- **`h2` (h²)**: Heritability - proportion of phenotypic variance due to genetics
  - h² = σ²g / (σ²g + σ²e)
  - Range: 0 to 1
  - h² = 0: No genetic contribution
  - h² = 1: All variance is genetic
  - Typical values: 0.2-0.8 for many traits

#### Assessing Accuracy from Cross-Validation

**Cross-validation provides unbiased accuracy estimates.** Always use `cv_gblup()` to assess model accuracy:

```r
# Run cross-validation
cv_result <- cv_gblup(
  qc_results = data,
  trait_name = "Yield",
  id_col_name = "ID",
  k = 5,
  seed = 2025,
  return_fold_preds = TRUE
)

# Access overall accuracy
accuracy <- cv_result$overall$mean_cor
cat("GBLUP Accuracy (Correlation):", round(accuracy, 4), "\n")

# Complete accuracy summary
cat("=== GBLUP Accuracy Summary ===\n")
cat("Correlation (Accuracy):", round(cv_result$overall$mean_cor, 4), "\n")
cat("RMSE:", round(cv_result$overall$mean_rmse, 4), "\n")
cat("Bias Slope:", round(cv_result$overall$mean_bias_slope, 4), "\n")
cat("Heritability:", round(cv_result$overall$mean_h2, 4), "\n")
```

#### Complete Workflow: Accuracy and Variance Components

```r
# ============================================================================
# Step 1: Run Cross-Validation (for unbiased accuracy)
# ============================================================================
cv_result <- cv_gblup(
  qc_results = data,
  trait_name = "Yield",
  id_col_name = "ID",
  k = 5,
  seed = 2025,
  return_fold_preds = TRUE
)

# Extract accuracy metrics
accuracy_cor <- cv_result$overall$mean_cor
accuracy_rmse <- cv_result$overall$mean_rmse
accuracy_bias <- cv_result$overall$mean_bias_slope

# Print accuracy results
cat("=== GBLUP Accuracy (from Cross-Validation) ===\n")
cat("Correlation (Accuracy):", round(accuracy_cor, 4), "\n")
cat("RMSE:", round(accuracy_rmse, 4), "\n")
cat("Bias Slope:", round(accuracy_bias, 4), "\n")
cat("Heritability:", round(cv_result$overall$mean_h2, 4), "\n")

# ============================================================================
# Step 2: Run Full GBLUP (for predictions on all individuals)
# ============================================================================
gblup_result <- run_gblup(
  qc_results = data,
  trait_name = "Yield",
  id_col_name = "ID",
  fixed_effects = ~1,
  drop_missing = TRUE,
  predict_all = TRUE
)

# Get variance components
cat("\n=== Variance Components ===\n")
cat("Genetic Variance (σ²g):", round(gblup_result$variances$sigma_g, 4), "\n")
cat("Residual Variance (σ²e):", round(gblup_result$variances$sigma_e, 4), "\n")
cat("Heritability (h²):", round(gblup_result$variances$h2, 4), "\n")

# ============================================================================
# Step 3: Visualize CV Results (Optional)
# ============================================================================
# Plot predicted vs observed from CV
if (!is.null(cv_result$predictions)) {
  plot(cv_result$predictions$yhat, cv_result$predictions$yobs,
       xlab = "Predicted", ylab = "Observed",
       main = "GBLUP Cross-Validation Accuracy")
  abline(0, 1, col = "red", lty = 2)
}
```

**Important Notes:**

1. **Accuracy from CV is unbiased**: The correlation from `cv_gblup()` is the standard accuracy metric in genomic prediction
2. **Variance components from full model**: Use `run_gblup()` to get variance components for all data
3. **Higher correlation = better accuracy**: Typical accuracy ranges from 0.3-0.9 depending on trait heritability and sample size
4. **Heritability affects accuracy**: Higher heritability generally leads to higher prediction accuracy

---

## Bayesian Methods

### Overview

Bayesian methods use Markov Chain Monte Carlo (MCMC) to estimate SNP effects and predict breeding values. `breedRplus` supports several Bayesian models via the BGLR package.

### Function: `run_bglr()`

**Purpose**: Fit a Bayesian genomic prediction model.

**Syntax**:
```r
run_bglr(
  obj,
  pheno,
  trait_name,
  id_col_name,
  model_type = "BayesA",
  n_iter = 5000,
  burn_in = 1000
)
```

**Parameters**:

| Parameter | Type | Description |
|-----------|------|-------------|
| `obj` | bigSNP | bigSNP object from `load_plink()` |
| `pheno` | data.frame | Phenotype data frame |
| `trait_name` | character | Name of trait column |
| `id_col_name` | character | Name of ID column |
| `model_type` | character | Model type: "BayesA", "BayesB", "BayesC", or "BL" (Bayesian Lasso) |
| `n_iter` | integer | Number of MCMC iterations (default: 5000) |
| `burn_in` | integer | Number of burn-in iterations (default: 1000) |

**Returns**: A list containing:
- `gebv`: Data frame with GEBVs
- `snp_effects`: Data frame with SNP effects
- `residual_var`: Residual variance estimate
- `model_fit`: Full BGLR model object

**Model Types**:

- **BayesA**: Assumes t-distributed SNP effects (good for traits with many small effects)
- **BayesB**: Assumes many SNPs have zero effect (sparse model)
- **BayesC**: Similar to BayesB but with different prior
- **BL** (Bayesian Lasso): Uses L1 regularization (similar to LASSO)

**Example 1: BayesA Model**

```r
# Fit BayesA model
bayes_result <- run_bglr(
  obj = data$snp_obj,
  pheno = data$pheno,
  trait_name = "Yield",
  id_col_name = "ID",
  model_type = "BayesA",
  n_iter = 5000,
  burn_in = 1000
)

# View GEBVs
head(bayes_result$gebv)

# View top SNP effects
head(bayes_result$snp_effects[order(abs(bayes_result$snp_effects$Effect), 
                                     decreasing = TRUE), ])
```

**Example 2: Bayesian Lasso**

```r
# Fit Bayesian Lasso model
bayes_result <- run_bglr(
  obj = data$snp_obj,
  pheno = data$pheno,
  trait_name = "Yield",
  id_col_name = "ID",
  model_type = "BL",  # Bayesian Lasso
  n_iter = 10000,     # More iterations for BL
  burn_in = 2000
)
```

### Function: `run_bglr_cv()`

**Purpose**: Perform k-fold cross-validation for Bayesian models.

**Syntax**:
```r
run_bglr_cv(
  geno,
  pheno,
  trait_name,
  id_col_name,
  model_type = "BayesA",
  n_iter = 5000,
  burn_in = 1000,
  k_folds = 5,
  seed = 123
)
```

**Parameters**: Similar to `run_bglr()`, plus:
- `geno`: Numeric genotype matrix (not bigSNP object)
- `k_folds`: Number of CV folds
- `seed`: Random seed

**Returns**: A list with:
- `y_true`: Observed values
- `y_pred`: Predicted values
- `cv_correlation`: Pearson correlation
- `cv_mse`: Mean squared error

**Example: Bayesian Cross-Validation**

```r
# Extract genotype matrix (required for run_bglr_cv)
geno_matrix <- data$snp_obj$genotypes[, ]
rownames(geno_matrix) <- as.character(data$snp_obj$fam$sample.ID)

# Run cross-validation
bayes_cv <- run_bglr_cv(
  geno = geno_matrix,
  pheno = data$pheno,
  trait_name = "Yield",
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
```

**Plotting Cross-Validation Results:**

To visualize the relationship between observed and predicted values:

```r
# Save plot to PNG file (recommended to avoid plot window size issues)
png("bayesian_cv_predictions.png", width = 800, height = 600, res = 100)

# Set plot margins
par(mar = c(4, 4, 3, 2))

# Create scatter plot
plot(bayes_cv$y_true, bayes_cv$y_pred,
     xlab = "Observed", ylab = "Predicted",
     main = "Bayesian Model Cross-Validation Predictions",
     pch = 19, col = "blue", cex = 0.7)

# Add 1:1 line (perfect prediction line)
abline(0, 1, col = "red", lty = 2, lwd = 2)

# Add correlation text
cor_val <- cor(bayes_cv$y_true, bayes_cv$y_pred, use = "complete.obs")
text(x = min(bayes_cv$y_true, na.rm = TRUE),
     y = max(bayes_cv$y_pred, na.rm = TRUE),
     labels = paste("r =", round(cor_val, 3)),
     pos = 4, col = "darkgreen", font = 2, cex = 1.2)

# Close the PNG device
dev.off()

cat("Plot saved as 'bayesian_cv_predictions.png'\n")
```

**Alternative: Plot directly to screen** (if plot window is large enough):

```r
# Close any existing plots
graphics.off()

# Set margins
par(mar = c(4, 4, 3, 2))

# Plot
plot(bayes_cv$y_true, bayes_cv$y_pred,
     xlab = "Observed", ylab = "Predicted",
     main = "Bayesian Model Cross-Validation Predictions",
     pch = 19, col = "blue")
abline(0, 1, col = "red", lty = 2, lwd = 2)
```

---

## Machine Learning Methods

### Overview

`breedRplus` supports several machine learning algorithms for genomic prediction:
- **RF**: Random Forest
- **XGB**: XGBoost (Gradient Boosting)
- **MLP**: Multi-Layer Perceptron (Neural Network)
- **CNN**: 1D Convolutional Neural Network
- **PLS**: Partial Least Squares Regression

### **Recommendation: Use Python for ML Methods**

**We strongly recommend using the Python version of `breedrplus` for machine learning methods**, especially for neural networks (MLP, CNN). Here's why:

#### **Advantages of Python for ML:**

1. **Easier Installation**: 
   - TensorFlow installation is straightforward: `pip install tensorflow`
   - No R→Python bridge complications
   - No deprecated package issues (R's `keras` is deprecated in favor of `keras3`)

2. **Better Performance**:
   - Native TensorFlow/Keras implementation (no R wrapper overhead)
   - More efficient memory usage for large datasets
   - Better GPU support for neural networks

3. **More Features & Updates**:
   - Python ML libraries are more actively developed
   - Access to latest TensorFlow/Keras features
   - Better documentation and community support

4. **Production Ready**:
   - Python is standard for ML in production
   - Better integration with ML pipelines
   - More deployment options

#### **When to Use R for ML:**

- You're already working entirely in R and prefer consistency
- You only need RF, XGB, or PLS (no TensorFlow dependency)
- Small to medium datasets where performance difference is negligible

#### **Hybrid Approach (Recommended):**

- **Use R for**: GBLUP, Bayesian methods, RF, XGB, PLS
- **Use Python for**: MLP, CNN (neural networks)

This gives you the best of both worlds: statistical methods in R, ML methods in Python.

### ML Methodology and Packages

The package uses the following R packages for machine learning implementations:

#### **1. Random Forest (RF)**
- **Package**: `ranger` - Fast implementation of Random Forest
- **Method**: Ensemble of decision trees with bagging
- **Default Parameters**: 500 trees
- **Use Case**: Good baseline, handles non-linear relationships well

#### **2. XGBoost (XGB)**
- **Package**: `xgboost` - Extreme Gradient Boosting
- **Method**: Gradient boosting with tree-based learners
- **Default Parameters**: 100 rounds, learning rate (eta) = 0.1, max_depth = 6
- **Use Case**: Often achieves high accuracy, handles complex patterns

#### **3. Multi-Layer Perceptron (MLP)**
- **Package**: `keras` - Deep learning framework
- **Method**: Feedforward neural network with 3 layers (64 → 32 → 1 units)
- **Architecture**: 
  - Input layer: All SNP features
  - Hidden layer 1: 64 units with ReLU activation + dropout (0.3)
  - Hidden layer 2: 32 units with ReLU activation
  - Output layer: 1 unit (continuous prediction)
- **Default Parameters**: 50 epochs, batch_size = 32, validation_split = 0.1
- **Use Case**: Captures complex non-linear relationships

#### **4. Convolutional Neural Network (CNN)**
- **Package**: `keras` - Deep learning framework
- **Method**: 1D convolutional layers for sequence-like genotype data
- **Architecture**:
  - Conv1D layer: 32 filters, kernel_size = 3
  - Max pooling: pool_size = 2
  - Flatten layer
  - Dense layer: 50 units with ReLU
  - Output layer: 1 unit
- **Default Parameters**: 50 epochs, batch_size = 32
- **Use Case**: Can capture local patterns in genotype sequences

#### **5. Partial Least Squares (PLS)**
- **Package**: `pls` - Partial Least Squares Regression
- **Method**: Dimensionality reduction and regression
- **Default Parameters**: Up to 50 components (or number of features if smaller)
- **Use Case**: Handles high-dimensional data, reduces overfitting

**Supporting Packages**:
- `magrittr`: Pipe operator (`%>%`) for keras model building
- `stats`: Prediction and correlation functions

### **Python Implementation (Recommended)**

The Python version uses native implementations:
- **RF**: `scikit-learn` RandomForestRegressor
- **XGB**: Native `xgboost` package
- **MLP/CNN**: Native `tensorflow/keras`
- **PLS**: `scikit-learn` PLSRegression

**Installation**:
```bash
pip install numpy pandas scipy scikit-learn xgboost tensorflow bed-reader
```

**Python Usage Example**:
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

# Train Random Forest model
rf_result = run_ml_model(
    pheno=data['pheno'],
    geno=geno_matrix,
    trait_name="Yield",
    id_col_name="ID",
    model_type="RF",
    n_trees=500
)

# View predictions
print(rf_result['gebv'].head())

# Cross-validation
rf_cv = run_ml_cv(
    pheno=data['pheno'],
    geno=geno_matrix,
    trait_name="Yield",
    id_col_name="ID",
    model_type="RF",
    k=5,
    seed=123,
    n_trees=500
)

# View accuracy
print(f"Correlation: {rf_cv['cv_correlation']:.4f}")
print(f"MSE: {rf_cv['cv_mse']:.4f}")

# Train Neural Network (MLP) - Much easier in Python!
mlp_result = run_ml_model(
    pheno=data['pheno'],
    geno=geno_matrix,
    trait_name="Yield",
    id_col_name="ID",
    model_type="MLP",
    epochs=50
)

# Cross-validation for MLP
mlp_cv = run_ml_cv(
    pheno=data['pheno'],
    geno=geno_matrix,
    trait_name="Yield",
    id_col_name="ID",
    model_type="MLP",
    k=5,
    seed=123,
    epochs=50
)
```

**See `ML_PYTHON_GUIDE.md` for a quick-start Python ML guide, or `README_PYTHON.md` and `breedrplus-py/` directory for complete Python documentation.**

### Function: `run_ml_model()` (R Version - Also Available)

**Note**: While R implementation is available, we recommend using Python for ML methods (see above).

**Purpose**: Train a machine learning model and predict GEBVs.

**Syntax**:
```r
run_ml_model(
  pheno,
  geno,
  trait_name,
  id_col_name,
  model_type = "RF",
  ...
)
```

**Parameters**:

| Parameter | Type | Description |
|-----------|------|-------------|
| `pheno` | data.frame | Phenotype data |
| `geno` | matrix | Numeric genotype matrix (rownames = IDs) |
| `trait_name` | character | Name of trait column |
| `id_col_name` | character | Name of ID column |
| `model_type` | character | Model type: "RF", "XGB", "MLP", "CNN", or "PLS" |
| `...` | | Model-specific parameters (see below) |

**Model-Specific Parameters**:

- **RF**: `n_trees` (default: 500) - Number of trees
- **XGB**: `n_rounds` (default: 100) - Number of boosting rounds
- **MLP/CNN**: `epochs` (default: 50) - Training epochs
- **PLS**: `n_comp` (default: min(50, ncol(geno))) - Number of components

**Returns**: A list with:
- `gebv`: Data frame with predicted GEBVs
- `model_fit`: Trained model object

**Example 1: Random Forest**

```r
# Extract genotype matrix
geno_matrix <- data$snp_obj$genotypes[, ]
rownames(geno_matrix) <- as.character(data$snp_obj$fam$sample.ID)

# Train Random Forest model
rf_result <- run_ml_model(
  pheno = data$pheno,
  geno = geno_matrix,
  trait_name = "Yield",
  id_col_name = "ID",
  model_type = "RF",
  n_trees = 500
)

# View predictions
head(rf_result$gebv)
```

**Example 2: XGBoost**

```r
# Train XGBoost model
xgb_result <- run_ml_model(
  pheno = data$pheno,
  geno = geno_matrix,
  trait_name = "Yield",
  id_col_name = "ID",
  model_type = "XGB",
  n_rounds = 100
)
```

**Example 3: Neural Network (MLP)**

```r
# Train Multi-Layer Perceptron
mlp_result <- run_ml_model(
  pheno = data$pheno,
  geno = geno_matrix,
  trait_name = "Yield",
  id_col_name = "ID",
  model_type = "MLP",
  epochs = 50
)
```

**Example 4: Convolutional Neural Network**

```r
# Train 1D-CNN
cnn_result <- run_ml_model(
  pheno = data$pheno,
  geno = geno_matrix,
  trait_name = "Yield",
  id_col_name = "ID",
  model_type = "CNN",
  epochs = 50
)
```

**Example 5: Partial Least Squares**

```r
# Train PLS model
pls_result <- run_ml_model(
  pheno = data$pheno,
  geno = geno_matrix,
  trait_name = "Yield",
  id_col_name = "ID",
  model_type = "PLS",
  n_comp = 50
)
```

### Function: `run_ml_cv()` (R Version - Also Available)

**Note**: While R implementation is available, we recommend using Python for ML methods (see above).

**Purpose**: Perform k-fold cross-validation for machine learning models.

**Syntax**:
```r
run_ml_cv(
  pheno,
  geno,
  trait_name,
  id_col_name,
  model_type = "RF",
  k = 5,
  seed = 123,
  ...
)
```

**Parameters**: Similar to `run_ml_model()`, plus:
- `k`: Number of CV folds
- `seed`: Random seed

**Returns**: A list with:
- `cv_gebv`: Data frame with CV predictions
- `cv_correlation`: Pearson correlation
- `cv_mse`: Mean squared error

**Example: Compare Multiple ML Models**

```r
# Extract genotype matrix
geno_matrix <- data$snp_obj$genotypes[, ]
rownames(geno_matrix) <- as.character(data$snp_obj$fam$sample.ID)

# Compare models
models <- c("RF", "XGB", "MLP", "PLS")
results <- list()

for (model in models) {
  cat("Running", model, "...\n")
  cv_result <- run_ml_cv(
    pheno = data$pheno,
    geno = geno_matrix,
    trait_name = "Yield",
    id_col_name = "ID",
    model_type = model,
    k = 5,
    seed = 123
  )
  results[[model]] <- cv_result
}

# Compare results
comparison <- data.frame(
  Model = models,
  Correlation = sapply(results, function(x) x$cv_correlation),
  MSE = sapply(results, function(x) x$cv_mse)
)

print(comparison)
```

**Plotting ML Cross-Validation Results:**

To visualize the relationship between observed and predicted values:

```r
# Run cross-validation
ml_cv <- run_ml_cv(
  pheno = data$pheno,
  geno = geno_matrix,
  trait_name = "Yield",
  id_col_name = "ID",
  model_type = "RF",
  k = 5,
  seed = 123,
  n_trees = 500
)

# Extract observed values
y_true <- data$pheno[[trait_name]][match(ml_cv$cv_gebv$ID, data$pheno[[id_col_name]])]
y_pred <- ml_cv$cv_gebv$gebv

# Save plot to PNG file
png("ml_cv_predictions.png", width = 800, height = 600, res = 100)

# Set plot margins
par(mar = c(4, 4, 3, 2))

# Create scatter plot
plot(y_true, y_pred,
     xlab = "Observed", ylab = "Predicted",
     main = "ML Model Cross-Validation Predictions (Random Forest)",
     pch = 19, col = "blue", cex = 0.7)

# Add 1:1 line (perfect prediction line)
abline(0, 1, col = "red", lty = 2, lwd = 2)

# Add correlation text
cor_val <- cor(y_true, y_pred, use = "complete.obs")
text(x = min(y_true, na.rm = TRUE),
     y = max(y_pred, na.rm = TRUE),
     labels = paste("r =", round(cor_val, 3)),
     pos = 4, col = "darkgreen", font = 2, cex = 1.2)

# Close the PNG device
dev.off()

cat("Plot saved as 'ml_cv_predictions.png'\n")
cat("Correlation:", round(ml_cv$cv_correlation, 4), "\n")
cat("MSE:", round(ml_cv$cv_mse, 4), "\n")
```

**Alternative: Plot directly to screen** (if plot window is large enough):

```r
# Close any existing plots
graphics.off()

# Set margins
par(mar = c(4, 4, 3, 2))

# Extract values
y_true <- data$pheno[[trait_name]][match(ml_cv$cv_gebv$ID, data$pheno[[id_col_name]])]
y_pred <- ml_cv$cv_gebv$gebv

# Plot
plot(y_true, y_pred,
     xlab = "Observed", ylab = "Predicted",
     main = "ML Model Cross-Validation Predictions",
     pch = 19, col = "blue")
abline(0, 1, col = "red", lty = 2, lwd = 2)
```

### **Important Notes for R Users:**

1. **TensorFlow Installation**: If using MLP or CNN in R, you'll need to install TensorFlow:
   ```r
   # Install keras3 (recommended) or keras (deprecated)
   install.packages("keras3")
   library(keras3)
   install_keras()  # Installs TensorFlow
   ```

2. **Performance**: Neural networks (MLP, CNN) will be slower in R due to R→Python bridge overhead

3. **Compatibility**: The `keras` package is deprecated; use `keras3` if available

4. **Recommendation**: For neural networks, strongly consider using the Python version instead

---

## Complete Workflow Examples

### Example 1: Full GBLUP Workflow

```r
library(breedRplus)
library(readxl)

# Step 1: Load data
pheno <- read_excel("phenotypes.xlsx")
data <- load_plink(
  bed_file = "genotypes.bed",
  bim_file = "genotypes.bim",
  fam_file = "genotypes.fam",
  pheno = pheno,
  id_col_name = "AnimalID"
)

# Step 2: Run GBLUP with fixed effects
gblup_result <- run_gblup(
  qc_results = data,
  trait_name = "MilkYield",
  id_col_name = "AnimalID",
  fixed_effects = ~ Herd + Parity,
  predict_all = TRUE
)

# Step 3: Cross-validation
cv_result <- cv_gblup(
  qc_results = data,
  trait_name = "MilkYield",
  id_col_name = "AnimalID",
  fixed_effects = ~ Herd + Parity,
  k = 10,
  seed = 2025
)

# Step 4: Evaluate and save results
print("Variance Components:")
print(gblup_result$variances)

print("Cross-Validation Results:")
print(cv_result$overall)

# Save GEBVs
write.csv(gblup_result$gebv_all, "GEBVs_all_animals.csv", row.names = FALSE)
```

### Example 2: Comparing Multiple Methods

```r
library(breedRplus)

# Load data
data <- load_plink(
  bed_file = "genotypes.bed",
  bim_file = "genotypes.bim",
  fam_file = "genotypes.fam",
  pheno = read.csv("phenotypes.csv"),
  id_col_name = "ID"
)

# Extract genotype matrix for ML methods
geno_matrix <- data$snp_obj$genotypes[, ]
rownames(geno_matrix) <- as.character(data$snp_obj$fam$sample.ID)

# Compare methods
methods <- list()

# 1. GBLUP
cat("Running GBLUP...\n")
methods$GBLUP <- cv_gblup(
  qc_results = data,
  trait_name = "Yield",
  id_col_name = "ID",
  k = 5
)

# 2. BayesA
cat("Running BayesA...\n")
methods$BayesA <- run_bglr_cv(
  geno = geno_matrix,
  pheno = data$pheno,
  trait_name = "Yield",
  id_col_name = "ID",
  model_type = "BayesA",
  k_folds = 5
)

# 3. Random Forest
cat("Running Random Forest...\n")
methods$RF <- run_ml_cv(
  pheno = data$pheno,
  geno = geno_matrix,
  trait_name = "Yield",
  id_col_name = "ID",
  model_type = "RF",
  k = 5
)

# 4. XGBoost
cat("Running XGBoost...\n")
methods$XGB <- run_ml_cv(
  pheno = data$pheno,
  geno = geno_matrix,
  trait_name = "Yield",
  id_col_name = "ID",
  model_type = "XGB",
  k = 5
)

# Compare results
comparison <- data.frame(
  Method = names(methods),
  Correlation = c(
    methods$GBLUP$overall$mean_cor,
    methods$BayesA$cv_correlation,
    methods$RF$cv_correlation,
    methods$XGB$cv_correlation
  ),
  MSE = c(
    methods$GBLUP$overall$mean_rmse^2,
    methods$BayesA$cv_mse,
    methods$RF$cv_mse,
    methods$XGB$cv_mse
  )
)

print(comparison)
```

### Example 3: Prediction for Selection Candidates

```r
library(breedRplus)

# Load data (includes both training and selection candidates)
data <- load_plink(
  bed_file = "all_genotypes.bed",
  bim_file = "all_genotypes.bim",
  fam_file = "all_genotypes.fam",
  pheno = read.csv("training_phenotypes.csv"),  # Only training animals have phenotypes
  id_col_name = "ID"
)

# Run GBLUP and predict for ALL individuals
gblup_result <- run_gblup(
  qc_results = data,
  trait_name = "Yield",
  id_col_name = "ID",
  predict_all = TRUE  # Predict for unphenotyped individuals
)

# Identify selection candidates (those without phenotypes)
training_ids <- data$pheno$ID
all_ids <- data$snp_obj$fam$sample.ID
candidate_ids <- setdiff(all_ids, training_ids)

# Extract GEBVs for selection candidates
candidate_gebvs <- gblup_result$gebv_all[
  gblup_result$gebv_all$ID %in% candidate_ids,
]

# Rank candidates
candidate_gebvs <- candidate_gebvs[order(candidate_gebvs$gebv, decreasing = TRUE), ]

# Select top 10%
top_10_percent <- head(candidate_gebvs, n = round(0.1 * nrow(candidate_gebvs)))

print("Top 10% Selection Candidates:")
print(top_10_percent)
```

---

## Troubleshooting

### Common Issues and Solutions

#### Issue 1: "File not found" Error

**Problem**: `load_plink()` cannot find PLINK files.

**Solution**:
```r
# Check file paths
file.exists("genotypes.bed")  # Should return TRUE

# Use absolute paths if needed
bed_file <- "/full/path/to/genotypes.bed"
```

#### Issue 2: ID Mismatch Between Genotype and Phenotype

**Problem**: Phenotype IDs don't match genotype sample IDs.

**Solution**:
```r
# Check IDs
geno_ids <- data$snp_obj$fam$sample.ID
pheno_ids <- data$pheno$ID

# Find mismatches
setdiff(geno_ids, pheno_ids)
setdiff(pheno_ids, geno_ids)

# Ensure ID column name is correct
data <- load_plink(..., id_col_name = "correct_column_name")
```

#### Issue 3: Memory Issues with Large Datasets

**Problem**: Running out of memory with large genotype files.

**Solutions**:
- Use `bigsnpr` efficiently (already implemented)
- Consider subsetting SNPs before analysis
- Increase available RAM
- Use a computing cluster for very large datasets

#### Issue 4: Cross-Validation Takes Too Long

**Problem**: CV is slow, especially for Bayesian methods.

**Solutions**:
```r
# Reduce number of folds
cv_result <- cv_gblup(..., k = 3)  # Instead of k = 10

# Reduce iterations for Bayesian methods
bayes_cv <- run_bglr_cv(..., n_iter = 2000, burn_in = 500)

# Use faster methods for initial screening
rf_cv <- run_ml_cv(..., model_type = "RF")  # Fast
```

#### Issue 5: Keras/Neural Network Errors

**Problem**: Errors when running MLP or CNN models.

**Solutions**:
- Ensure TensorFlow is installed: `install.packages("tensorflow")` and `tensorflow::install_tensorflow()`
- Check that Python is available: `reticulate::py_available()`
- Try reducing number of epochs or features

#### Issue 6: Low Prediction Accuracy

**Problem**: CV correlations are low (< 0.3).

**Possible Causes**:
- Small training population size
- Low heritability trait
- Poor genotype-phenotype alignment
- Need for quality control (remove low-quality SNPs)

**Solutions**:
- Increase training population size
- Check data quality
- Try different models
- Consider including fixed effects

---

## Best Practices

### 1. Data Quality Control

- **Before analysis**: Check for missing data, outliers, and data quality
- **Genotype QC**: Consider filtering SNPs by minor allele frequency (MAF), call rate, etc.
- **Phenotype QC**: Remove obvious outliers, check for data entry errors

### 2. Model Selection

- **Start with GBLUP**: It's fast and provides a good baseline
- **Try multiple methods**: Different methods work better for different traits
- **Use cross-validation**: Always validate your models

### 3. Cross-Validation Strategy

- **Use appropriate k**: 5-10 folds is usually sufficient
- **Stratify when possible**: If you have family structure, use `stratify_by`
- **Set random seed**: For reproducibility
- **Report both fold-level and overall metrics**

### 4. Computational Efficiency

- **For large datasets**: Start with faster methods (GBLUP, RF)
- **For Bayesian methods**: Start with fewer iterations to test, then increase
- **Parallel processing**: Consider parallelizing CV folds (future feature)

### 5. Interpretation of Results

- **Heritability (h²)**: Should be biologically reasonable for your trait
- **CV correlation**: Higher is better, but depends on trait heritability
- **Bias**: Slope close to 1 indicates unbiased predictions
- **Compare methods**: Use multiple metrics (correlation, MSE, bias)

### 6. Reporting

Always report:
- Sample size (training and validation)
- Number of SNPs
- Model parameters (iterations, folds, etc.)
- Variance components (for GBLUP)
- Cross-validation metrics
- Prediction accuracy

---

## References

### Key Papers

- **GBLUP**: VanRaden, P. M. (2008). Efficient methods to compute genomic predictions. *Journal of dairy science*, 91(11), 4414-4423.

- **Bayesian Methods**: de los Campos, G., et al. (2013). Whole-genome regression and prediction methods applied to plant and animal breeding. *Genetics*, 193(2), 327-345.

- **Machine Learning**: González-Recio, O., et al. (2014). Machine learning methods and predictive ability metrics for genome-wide prediction of complex traits. *Livestock Science*, 166, 217-231.

### Package Dependencies

- `bigsnpr`: Privé, F., et al. (2018). Efficient analysis of large-scale genome-wide data with two R packages: bigstatsr and bigsnpr. *Bioinformatics*, 34(16), 2781-2787.

- `BGLR`: Pérez, P., & de los Campos, G. (2014). Genome-wide regression and prediction with the BGLR statistical package. *Genetics*, 198(2), 483-495.

---

## Getting Help

### Documentation

- Package help: `?function_name` (e.g., `?run_gblup`)
- Vignettes: `browseVignettes("breedRplus")`
- Example workflow: `system.file("example_workflow.R", package = "breedRplus")`

### Support

- GitHub Issues: https://github.com/higefee-218/BreedRPlus.git
- Email: higefee@gmail.com

### Contributing

Contributions are welcome! Please see the contributing guidelines in the repository.

---

## Appendix: Function Quick Reference

| Function | Purpose | Key Parameters |
|----------|---------|----------------|
| `load_plink()` | Load PLINK files | `bed_file`, `bim_file`, `fam_file`, `pheno` |
| `run_gblup()` | Fit GBLUP model | `qc_results`, `trait_name`, `predict_all` |
| `cv_gblup()` | GBLUP cross-validation | `qc_results`, `trait_name`, `k` |
| `run_bglr()` | Fit Bayesian model | `obj`, `pheno`, `model_type` |
| `run_bglr_cv()` | Bayesian cross-validation | `geno`, `pheno`, `model_type`, `k_folds` |
| `run_ml_model()` | Train ML model | `pheno`, `geno`, `model_type` |
| `run_ml_cv()` | ML cross-validation | `pheno`, `geno`, `model_type`, `k` |

---

**Last Updated**: 2025

**Package Version**: 0.0.0.9000

