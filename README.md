# breedRplus

Genomic Prediction Tools for Animal and Plant Breeding

## Description

`breedRplus` is a comprehensive R package for genomic prediction in breeding programs. It provides tools for:

- **GBLUP**: Genomic Best Linear Unbiased Prediction with cross-validation
- **Bayesian Methods**: BGLR-based Bayesian regression (BayesA, BayesB, BayesC, BL)
- **Machine Learning**: Random Forest, XGBoost, Neural Networks (MLP, CNN), and PLS regression
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

# Cross-validation
cv_result <- cv_gblup(
  qc_results = data,
  trait_name = "yield",
  id_col_name = "ID",
  k = 5
)
```

### 3. Run Bayesian Methods

```r
# Fit Bayesian model
bayes_result <- run_bglr(
  obj = data$snp_obj,
  pheno = data$pheno,
  trait_name = "yield",
  id_col_name = "ID",
  model_type = "BayesA"
)

# Cross-validation
bayes_cv <- run_bglr_cv(
  geno = genotype_matrix,
  pheno = data$pheno,
  trait_name = "yield",
  id_col_name = "ID",
  model_type = "BayesA",
  k_folds = 5
)
```

### 4. Run Machine Learning Models

```r
# Train ML model
ml_result <- run_ml_model(
  pheno = data$pheno,
  geno = genotype_matrix,
  trait_name = "yield",
  id_col_name = "ID",
  model_type = "RF"  # or "XGB", "MLP", "CNN", "PLS"
)

# Cross-validation
ml_cv <- run_ml_cv(
  pheno = data$pheno,
  geno = genotype_matrix,
  trait_name = "yield",
  id_col_name = "ID",
  model_type = "RF",
  k = 5
)
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

[Fei Ge]

