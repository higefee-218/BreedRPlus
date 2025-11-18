# Machine Learning with Python - Quick Guide

## Why Use Python for ML?

We **strongly recommend using Python** for machine learning methods in `breedrplus`, especially for neural networks (MLP, CNN). Here's why:

- ✅ **Easier Installation**: `pip install tensorflow` (no R→Python bridge issues)
- ✅ **Better Performance**: Native TensorFlow/Keras (no wrapper overhead)
- ✅ **More Features**: Access to latest ML libraries and updates
- ✅ **Production Ready**: Standard for ML deployment

## Installation

### Step 1: Install Python Dependencies

```bash
pip install numpy pandas scipy scikit-learn xgboost tensorflow bed-reader
```

Or install the full package:

```bash
cd breedrplus-py
pip install -e .
```

### Step 2: Verify Installation

```python
import tensorflow as tf
print(f"TensorFlow version: {tf.__version__}")

# Test that it works
tf.constant("Hello, TensorFlow!")
```

## Quick Start

### Complete ML Workflow in Python

```python
from breedrplus import load_plink, run_ml_model, run_ml_cv
import pandas as pd
import numpy as np

# ============================================================================
# Step 1: Load Data
# ============================================================================
data = load_plink(
    bed_file="genotypes.bed",
    bim_file="genotypes.bim",
    fam_file="genotypes.fam",
    pheno=phenotype_df,  # pandas DataFrame
    id_col_name="ID"
)

# Extract genotype matrix
geno_matrix = data['snp_obj']['genotypes']

# ============================================================================
# Step 2: Train ML Models
# ============================================================================

# Random Forest
rf_result = run_ml_model(
    pheno=data['pheno'],
    geno=geno_matrix,
    trait_name="Yield",
    id_col_name="ID",
    model_type="RF",
    n_trees=500
)

print("RF Predictions:")
print(rf_result['gebv'].head())

# XGBoost
xgb_result = run_ml_model(
    pheno=data['pheno'],
    geno=geno_matrix,
    trait_name="Yield",
    id_col_name="ID",
    model_type="XGB",
    n_rounds=100
)

# Neural Network (MLP) - Much easier in Python!
mlp_result = run_ml_model(
    pheno=data['pheno'],
    geno=geno_matrix,
    trait_name="Yield",
    id_col_name="ID",
    model_type="MLP",
    epochs=50
)

# Convolutional Neural Network (CNN)
cnn_result = run_ml_model(
    pheno=data['pheno'],
    geno=geno_matrix,
    trait_name="Yield",
    id_col_name="ID",
    model_type="CNN",
    epochs=50
)

# Partial Least Squares
pls_result = run_ml_model(
    pheno=data['pheno'],
    geno=geno_matrix,
    trait_name="Yield",
    id_col_name="ID",
    model_type="PLS",
    n_comp=50
)

# ============================================================================
# Step 3: Cross-Validation
# ============================================================================

# Compare all models
models = ["RF", "XGB", "MLP", "CNN", "PLS"]
results = {}

for model in models:
    print(f"\nRunning {model} cross-validation...")
    cv_result = run_ml_cv(
        pheno=data['pheno'],
        geno=geno_matrix,
        trait_name="Yield",
        id_col_name="ID",
        model_type=model,
        k=5,
        seed=123
    )
    results[model] = cv_result
    print(f"{model} - Correlation: {cv_result['cv_correlation']:.4f}, MSE: {cv_result['cv_mse']:.4f}")

# Compare results
comparison = pd.DataFrame({
    'Model': models,
    'Correlation': [results[m]['cv_correlation'] for m in models],
    'MSE': [results[m]['cv_mse'] for m in models]
})
print("\nModel Comparison:")
print(comparison)

# ============================================================================
# Step 4: Visualize Results
# ============================================================================

import matplotlib.pyplot as plt

# Plot predictions vs observed for best model
best_model = comparison.loc[comparison['Correlation'].idxmax(), 'Model']
best_cv = results[best_model]

# Extract observed values
y_true = data['pheno']['Yield'].values
y_pred = best_cv['cv_gebv']['gebv'].values

# Create plot
plt.figure(figsize=(8, 6))
plt.scatter(y_true, y_pred, alpha=0.6, color='blue')
plt.plot([y_true.min(), y_true.max()], [y_true.min(), y_true.max()], 
         'r--', lw=2, label='1:1 line')
plt.xlabel('Observed', fontsize=12)
plt.ylabel('Predicted', fontsize=12)
plt.title(f'{best_model} Cross-Validation Predictions', fontsize=14)
plt.legend()

# Add correlation text
cor_val = best_cv['cv_correlation']
plt.text(0.05, 0.95, f'r = {cor_val:.3f}', 
         transform=plt.gca().transAxes,
         fontsize=12, verticalalignment='top',
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

plt.tight_layout()
plt.savefig('ml_cv_predictions.png', dpi=100)
plt.close()

print(f"\nPlot saved as 'ml_cv_predictions.png'")
```

## Model-Specific Parameters

### Random Forest (RF)
```python
run_ml_model(
    ...,
    model_type="RF",
    n_trees=500  # Number of trees (default: 500)
)
```

### XGBoost (XGB)
```python
run_ml_model(
    ...,
    model_type="XGB",
    n_rounds=100  # Number of boosting rounds (default: 100)
)
```

### Multi-Layer Perceptron (MLP)
```python
run_ml_model(
    ...,
    model_type="MLP",
    epochs=50  # Training epochs (default: 50)
)
```

### Convolutional Neural Network (CNN)
```python
run_ml_model(
    ...,
    model_type="CNN",
    epochs=50  # Training epochs (default: 50)
)
```

### Partial Least Squares (PLS)
```python
run_ml_model(
    ...,
    model_type="PLS",
    n_comp=50  # Number of components (default: min(50, n_features))
)
```

## Cross-Validation Example

```python
# 5-fold cross-validation for Random Forest
rf_cv = run_ml_cv(
    pheno=data['pheno'],
    geno=geno_matrix,
    trait_name="Yield",
    id_col_name="ID",
    model_type="RF",
    k=5,           # Number of folds
    seed=123,      # Random seed
    n_trees=500    # Model-specific parameter
)

# Access results
print(f"Correlation: {rf_cv['cv_correlation']:.4f}")
print(f"MSE: {rf_cv['cv_mse']:.4f}")
print(f"\nPredictions:")
print(rf_cv['cv_gebv'].head())
```

## Comparing Multiple Models

```python
# Compare all models
models = ["RF", "XGB", "MLP", "CNN", "PLS"]
results = {}

for model in models:
    print(f"Running {model}...")
    cv_result = run_ml_cv(
        pheno=data['pheno'],
        geno=geno_matrix,
        trait_name="Yield",
        id_col_name="ID",
        model_type=model,
        k=5,
        seed=123
    )
    results[model] = cv_result

# Create comparison table
comparison = pd.DataFrame({
    'Model': models,
    'Correlation': [results[m]['cv_correlation'] for m in models],
    'MSE': [results[m]['cv_mse'] for m in models]
}).sort_values('Correlation', ascending=False)

print("\nModel Performance Comparison:")
print(comparison)
```

## Troubleshooting

### TensorFlow Not Found
```bash
pip install tensorflow
```

### Import Errors
Make sure you're in the correct directory and have installed the package:
```bash
cd breedrplus-py
pip install -e .
```

### Memory Issues
For very large datasets, consider:
- Using fewer epochs for neural networks
- Reducing number of trees for RF
- Using a subset of features

## Integration with R Workflow

You can use Python for ML while keeping other analyses in R:

1. **In R**: Run GBLUP, Bayesian methods
2. **Export data**: Save genotype matrix and phenotypes
3. **In Python**: Run ML models
4. **Compare results**: Combine outputs from both

## See Also

- `README_PYTHON.md` - Complete Python package documentation
- `breedrplus-py/` - Python package source code
- `example_workflow.py` - Complete example workflow

