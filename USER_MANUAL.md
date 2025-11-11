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
- `gebv_obs`: Data frame with GEBVs for observed individuals
- `gebv_all`: Data frame with GEBVs for all individuals (if `predict_all=TRUE`)
- `variances`: List with `sigma_g` (genetic variance), `sigma_e` (residual variance), `h2` (heritability)
- `model_fit`: The full mixed model fit object
- `G`: The genomic relationship matrix

**Example 1: Basic GBLUP**

```r
# Run GBLUP with intercept only
gblup_result <- run_gblup(
  qc_results = data,
  trait_name = "Yield",
  id_col_name = "ID"
)

# View variance components
print(gblup_result$variances)
# $sigma_g
# [1] 0.45
# $sigma_e
# [1] 0.55
# $h2
# [1] 0.45

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
- `fold_results`: Data frame with metrics for each fold
- `overall`: Data frame with overall CV metrics
- `predictions`: Data frame with predictions for each test set (if `return_fold_preds=TRUE`)

**Example: 5-Fold Cross-Validation**

```r
# Run 5-fold cross-validation
cv_result <- cv_gblup(
  qc_results = data,
  trait_name = "Yield",
  id_col_name = "ID",
  k = 5,
  seed = 2025
)

# View overall metrics
print(cv_result$overall)
#   mean_cor mean_rmse mean_bias_slope mean_h2
# 1    0.65      2.34           0.98    0.45

# View fold-specific results
print(cv_result$fold_results)
#   fold n_train n_test   cor   rmse bias_slope    Vu    Ve    h2_est
# 1    1     400    100  0.64  2.45       0.97  0.44  0.56     0.44
# 2    2     400    100  0.66  2.23       0.99  0.46  0.54     0.46
# ...

# View predictions
head(cv_result$predictions)
```

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
```

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
# Extract genotype matrix
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

# View results
print(bayes_cv$cv_correlation)
print(bayes_cv$cv_mse)

# Plot predictions vs observed
plot(bayes_cv$y_true, bayes_cv$y_pred,
     xlab = "Observed", ylab = "Predicted",
     main = "Bayesian Model Predictions")
abline(0, 1, col = "red")
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

### Function: `run_ml_model()`

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

### Function: `run_ml_cv()`

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

- GitHub Issues: [Your repository URL]
- Email: [Your email]

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

