# breedRplus 用户手册

**动物和植物育种的基因组预测工具**

版本 0.0.0.9000

---

## 目录

1. [简介](#简介)
2. [安装](#安装)
3. [快速入门指南](#快速入门指南)
4. [数据准备](#数据准备)
5. [GBLUP方法](#gblup方法)
6. [贝叶斯方法](#贝叶斯方法)
7. [机器学习方法](#机器学习方法)
8. [完整工作流程示例](#完整工作流程示例)
9. [故障排除](#故障排除)
10. [最佳实践](#最佳实践)
11. [参考文献](#参考文献)

---

## 简介

`breedRplus` 是一个专为动物和植物育种项目设计的综合R包，用于基因组预测。它提供了多种基因组预测方法的统一接口：

- **GBLUP**（基因组最佳线性无偏预测）
- **贝叶斯方法**（BayesA、BayesB、BayesC、贝叶斯Lasso）
- **机器学习**（随机森林、XGBoost、神经网络、PLS）

### 主要特性

- 使用 `bigsnpr` 高效处理大规模基因型数据
- 支持PLINK格式的基因型文件
- 用于模型评估的交叉验证工具
- 可预测已表型和未表型个体
- 全面的文档和示例

---

## 安装

### 系统要求

在安装 `breedRplus` 之前，请确保已安装 R（≥ 4.0）。该包需要以下依赖项：

- `bigsnpr` - 用于处理大型基因型数据
- `rrBLUP` - 用于GBLUP实现
- `BGLR` - 用于贝叶斯方法
- `ranger` - 用于随机森林
- `xgboost` - 用于梯度提升
- `keras` - 用于神经网络
- `pls` - 用于偏最小二乘回归

### 从GitHub安装

```r
# 如果尚未安装devtools，请先安装
if (!require("devtools")) install.packages("devtools")

# 安装breedRplus
devtools::install_github("higefee-218/breedRplus")

# 加载包
library(breedRplus)
```

### 从源代码安装

```r
# 下载并解压包
# 然后安装：
devtools::install("path/to/breedRplus")
```

---

## 快速入门指南

### 示例：5步完整工作流程

```r
library(breedRplus)
library(readxl)  # 用于读取Excel文件

# 步骤1：加载数据
data <- load_plink(
  bed_file = "genotypes.bed",
  bim_file = "genotypes.bim",
  fam_file = "genotypes.fam",
  pheno = read_excel("phenotypes.xlsx"),
  id_col_name = "ID"
)

# 步骤2：运行GBLUP
gblup_result <- run_gblup(
  qc_results = data,
  trait_name = "Yield",
  id_col_name = "ID",
  predict_all = TRUE
)

# 步骤3：查看结果
print(gblup_result$variances)  # 方差组分
head(gblup_result$gebv_all)    # 预测的育种值

# 步骤4：交叉验证
cv_result <- cv_gblup(
  qc_results = data,
  trait_name = "Yield",
  id_col_name = "ID",
  k = 5
)

# 步骤5：评估性能
print(cv_result$overall)  # 整体CV指标
```

---

## 数据准备

### 加载PLINK文件

`load_plink()` 函数是加载基因型数据的入口点。

#### 函数：`load_plink()`

**用途**：将PLINK格式的基因型文件（.bed、.bim、.fam）加载到R中，并可选择性地对齐表型数据。

**语法**：
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

**参数**：

| 参数 | 类型 | 描述 |
|------|------|------|
| `bed_file` | character | .bed文件路径（二进制基因型数据） |
| `bim_file` | character | .bim文件路径（变异信息） |
| `fam_file` | character | .fam文件路径（样本信息） |
| `pheno` | data.frame | 可选的表型数据框 |
| `id_col_name` | character | pheno中ID列的名称（如果提供pheno则必需） |
| `impute_method` | character | 插补方法："mode"、"mean0"、"mean2"、"random"（默认："mode"） |
| `backingfile` | character | bigSNP备份文件路径（默认：临时文件） |

**返回值**：包含以下内容的列表：
- `snp_obj`：包含基因型数据的bigSNP对象
- `pheno`：与基因型顺序对齐的表型数据（如果提供）

**示例1：仅加载基因型**

```r
# 加载PLINK文件，不包含表型
data <- load_plink(
  bed_file = "genotypes.bed",
  bim_file = "genotypes.bim",
  fam_file = "genotypes.fam"
)

# 检查结构
str(data$snp_obj)
head(data$snp_obj$fam)  # 样本信息
```

**示例2：加载基因型和表型**

```r
# 加载表型数据
pheno_df <- read.csv("phenotypes.csv")

# 加载PLINK文件并对齐表型
data <- load_plink(
  bed_file = "genotypes.bed",
  bim_file = "genotypes.bim",
  fam_file = "genotypes.fam",
  pheno = pheno_df,
  id_col_name = "AnimalID"
)

# 验证对齐
head(data$pheno)
```

**示例3：使用示例数据**

```r
# 访问包中包含的示例数据
bed_file <- system.file("extdata", "500ind.30K.bed", package = "breedRplus")
bim_file <- system.file("extdata", "500ind.30K.bim", package = "breedRplus")
fam_file <- system.file("extdata", "500ind.30K.fam", package = "breedRplus")
pheno_file <- system.file("extdata", "500_ind_pheno.xlsx", package = "breedRplus")

# 加载数据
library(readxl)
pheno <- read_excel(pheno_file)

data <- load_plink(
  bed_file = bed_file,
  bim_file = bim_file,
  fam_file = fam_file,
  pheno = pheno,
  id_col_name = "ID"  # 更新以匹配您的数据
)
```

**重要提示**：

- 函数会自动使用指定方法插补缺失基因型
- 表型数据会自动对齐以匹配基因型样本顺序
- 仅保留同时存在于基因型和表型数据中的个体
- 备份文件高效存储大型数据集的基因型数据

---

## GBLUP方法

### 概述

GBLUP（基因组最佳线性无偏预测）是一种混合模型方法，使用基因组关系矩阵（G矩阵）来预测育种值。

### GBLUP方法论

`breedRplus` 包实现了VanRaden (2008)描述的标准GBLUP方法。实现使用以下组件：

#### 1. **基因组关系矩阵（G矩阵）计算**

使用**VanRaden方法1**（VanRaden, 2008）计算基因组关系矩阵：

**公式**：G = Z Z' / m

其中：
- **Z** 是标准化基因型矩阵：Z = (X - 2p) / √[2p(1-p)]
  - X：基因型矩阵（编码为0、1、2，分别代表AA、Aa、aa）
  - p：等位基因频率（从数据计算）
- **m** 是标记（SNP）数量

**实现细节**：
- 使用 `bigsnpr::bed_tcrossprodSelf()` 和 `bed_scaleBinom` 缩放函数
- 直接从PLINK .bed文件高效计算G矩阵
- 在对角线上添加小值（默认：1e-6）以提高数值稳定性
- G矩阵表示所有个体对之间的基因组相似性

#### 2. **混合模型框架**

GBLUP模型遵循标准混合模型：

**y = Xβ + Zu + e**

其中：
- **y**：表型观测值向量
- **X**：固定效应的设计矩阵
- **β**：固定效应系数向量
- **Z**：将观测值链接到遗传值的关联矩阵（在GBLUP中为单位矩阵）
- **u**：随机遗传效应向量（GEBV）~ N(0, Gσ²g)
- **e**：残差向量 ~ N(0, Iσ²e)

**方差组分**：
- **σ²g**：遗传方差（育种值的方差）
- **σ²e**：残差方差（环境+误差方差）
- **遗传力**：h² = σ²g / (σ²g + σ²e)

#### 3. **模型拟合**

使用**REML（限制性最大似然）**通过 `rrBLUP::mixed.solve()` 求解混合模型：
- 使用REML估计方差组分（σ²g, σ²e）
- 求解Henderson混合模型方程以获得：
  - 固定效应估计（β）
  - 育种值的最佳线性无偏预测（BLUP）（u）

#### 4. **未表型个体的预测**

当 `predict_all = TRUE` 时，使用以下公式预测所有个体（包括没有表型的个体）的GEBV：

**u_all = G_all,obs × (G_obs + λI)⁻¹ × (y_obs - X_obsβ)**

其中：
- **G_all,obs**：所有个体与已观测个体之间的基因组关系
- **G_obs**：已观测个体的基因组关系矩阵
- **λ**：残差与遗传方差的比率（σ²e / σ²g）
- **y_obs**：观测表型
- **X_obsβ**：已观测个体的固定效应

#### 5. **交叉验证**

K折交叉验证通过以下步骤进行：
1. 将个体随机分为k折
2. 对于每一折：
   - 训练集：除该折外的所有个体
   - 测试集：该折中的个体
   - 在训练集上拟合GBLUP
   - 使用基因组关系预测测试集个体
3. 计算每折的准确性指标（相关性、RMSE、偏差）
4. 对所有折的指标取平均

**关键特性**：
- **分层CV**：可选的按组（例如，家系）分层，确保每折都有来自所有组的代表
- **无偏准确性**：CV提供预测准确性的无偏估计

#### 软件实现

- **G矩阵计算**：`bigsnpr` 包（Privé et al., 2018）
  - 高效处理大规模基因型数据
  - 直接从PLINK .bed文件计算
- **混合模型求解**：`rrBLUP` 包（Endelman, 2011）
  - 方差组分的REML估计
  - 育种值的BLUP计算

#### 参考文献

- **VanRaden, P. M.** (2008). Efficient methods to compute genomic predictions. *Journal of dairy science*, 91(11), 4414-4423.
- **Endelman, J. B.** (2011). Ridge regression and other kernels for genomic selection with R package rrBLUP. *The Plant Genome*, 4(3), 250-255.
- **Privé, F., et al.** (2018). Efficient analysis of large-scale genome-wide data with two R packages: bigstatsr and bigsnpr. *Bioinformatics*, 34(16), 2781-2787.

### 函数：`run_gblup()`

**用途**：拟合GBLUP模型并预测基因组估计育种值（GEBV）。

**语法**：
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

**参数**：

| 参数 | 类型 | 描述 |
|------|------|------|
| `qc_results` | list | `load_plink()` 的输出 |
| `trait_name` | character | pheno中性状列的名称 |
| `id_col_name` | character | pheno中ID列的名称 |
| `fixed_effects` | formula | 固定效应公式（默认：`~1` 仅截距） |
| `drop_missing` | logical | 删除缺失表型的个体？（默认：TRUE） |
| `predict_all` | logical | 预测所有个体包括未表型个体？（默认：FALSE） |
| `tol_diag` | numeric | 添加到G矩阵对角线的较小值以提高稳定性（默认：1e-6） |

**返回值**：包含以下内容的列表：
- `gebv_obs`：已观测个体的GEBV数据框
- `gebv_all`：所有个体的GEBV数据框（如果 `predict_all=TRUE`）
- `variances`：包含 `sigma_g`（遗传方差）、`sigma_e`（残差方差）、`h2`（遗传力）的列表
- `model_fit`：完整的混合模型拟合对象
- `G`：基因组关系矩阵

**示例1：基本GBLUP**

```r
# 仅使用截距运行GBLUP
gblup_result <- run_gblup(
  qc_results = data,
  trait_name = "Yield",
  id_col_name = "ID"
)

# 查看方差组分
print(gblup_result$variances)
# $sigma_g
# [1] 0.45
# $sigma_e
# [1] 0.55
# $h2
# [1] 0.45

# 查看已观测个体的GEBV
head(gblup_result$gebv_obs)
```

**示例2：带固定效应的GBLUP**

```r
# 运行带固定效应的GBLUP（例如批次、性别）
gblup_result <- run_gblup(
  qc_results = data,
  trait_name = "Yield",
  id_col_name = "ID",
  fixed_effects = ~ Batch + Sex
)

# 查看结果
print(gblup_result$variances)
```

**示例3：预测所有个体**

```r
# 预测已表型和未表型个体的GEBV
gblup_result <- run_gblup(
  qc_results = data,
  trait_name = "Yield",
  id_col_name = "ID",
  predict_all = TRUE  # 重要！
)

# 已观测个体的GEBV
head(gblup_result$gebv_obs)

# 所有个体的GEBV（包括未表型个体）
head(gblup_result$gebv_all)
```

### 函数：`cv_gblup()`

**用途**：对GBLUP模型进行k折交叉验证。

**语法**：
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

**参数**：

| 参数 | 类型 | 描述 |
|------|------|------|
| `qc_results` | list | `load_plink()` 的输出 |
| `trait_name` | character | 性状列名称 |
| `id_col_name` | character | ID列名称 |
| `k` | integer | CV折数（默认：5） |
| `fixed_effects` | formula | 固定效应公式 |
| `seed` | integer | 随机种子以确保可重复性 |
| `stratify_by` | character | 用于分层CV的可选列名（例如"Family"） |
| `drop_missing` | logical | 删除缺失性状的个体？ |
| `return_fold_preds` | logical | 包含折级预测？ |

**返回值**：包含以下内容的列表：
- `fold_results`：每折指标的数据框
- `overall`：整体CV指标的数据框
- `predictions`：每个测试集预测的数据框（如果 `return_fold_preds=TRUE`）

**示例：5折交叉验证**

```r
# 运行5折交叉验证
cv_result <- cv_gblup(
  qc_results = data,
  trait_name = "Yield",
  id_col_name = "ID",
  k = 5,
  seed = 2025
)

# 查看整体指标
print(cv_result$overall)
#   mean_cor mean_rmse mean_bias_slope mean_h2
# 1    0.65      2.34           0.98    0.45

# 查看折特定结果
print(cv_result$fold_results)
#   fold n_train n_test   cor   rmse bias_slope    Vu    Ve    h2_est
# 1    1     400    100  0.64  2.45       0.97  0.44  0.56     0.44
# 2    2     400    100  0.66  2.23       0.99  0.46  0.54     0.46
# ...

# 查看预测
head(cv_result$predictions)
```

**示例：分层交叉验证**

```r
# 按家系分层，确保每折都有来自所有家系的代表
cv_result <- cv_gblup(
  qc_results = data,
  trait_name = "Yield",
  id_col_name = "ID",
  k = 5,
  stratify_by = "Family",  # pheno中的列名
  seed = 2025
)
```

### **评估GBLUP准确性和方差组分**

#### 访问方差组分和遗传力

```r
# 运行GBLUP
gblup_result <- run_gblup(
  qc_results = data,
  trait_name = "Yield",
  id_col_name = "ID"
)

# 访问方差组分
variances <- gblup_result$variances
sigma_g <- variances$sigma_g
sigma_e <- variances$sigma_e
h2 <- variances$h2

# 打印结果
cat("遗传方差 (σ²g):", round(sigma_g, 4), "\n")
cat("残差方差 (σ²e):", round(sigma_e, 4), "\n")
cat("遗传力 (h²):", round(h2, 4), "\n")
```

**理解方差组分：**

- **`sigma_g` (σ²g)**：遗传方差 - 由于个体间遗传差异导致的方差
- **`sigma_e` (σ²e)**：残差方差 - 由于环境效应和测量误差导致的方差
- **`h2` (h²)**：遗传力 - 由于遗传导致的表型方差比例
  - h² = σ²g / (σ²g + σ²e)
  - 范围：0到1
  - h² = 0：无遗传贡献
  - h² = 1：所有方差都是遗传的
  - 典型值：许多性状为0.2-0.8

#### 从交叉验证评估准确性

**交叉验证提供无偏的准确性估计。**始终使用 `cv_gblup()` 评估模型准确性：

```r
# 运行交叉验证
cv_result <- cv_gblup(
  qc_results = data,
  trait_name = "Yield",
  id_col_name = "ID",
  k = 5,
  seed = 2025,
  return_fold_preds = TRUE
)

# 访问整体准确性
accuracy <- cv_result$overall$mean_cor
cat("GBLUP准确性（相关性）:", round(accuracy, 4), "\n")

# 完整准确性摘要
cat("=== GBLUP准确性摘要 ===\n")
cat("相关性（准确性）:", round(cv_result$overall$mean_cor, 4), "\n")
cat("RMSE:", round(cv_result$overall$mean_rmse, 4), "\n")
cat("偏差斜率:", round(cv_result$overall$mean_bias_slope, 4), "\n")
cat("遗传力:", round(cv_result$overall$mean_h2, 4), "\n")
```

#### 完整工作流程：准确性和方差组分

```r
# ============================================================================
# 步骤1：运行交叉验证（用于无偏准确性）
# ============================================================================
cv_result <- cv_gblup(
  qc_results = data,
  trait_name = "Yield",
  id_col_name = "ID",
  k = 5,
  seed = 2025,
  return_fold_preds = TRUE
)

# ============================================================================
# 步骤2：运行完整GBLUP（用于方差组分和所有个体的预测）
# ============================================================================
gblup_result <- run_gblup(
  qc_results = data,
  trait_name = "Yield",
  id_col_name = "ID",
  predict_all = TRUE
)

# ============================================================================
# 步骤3：报告结果
# ============================================================================
cat("=== GBLUP结果摘要 ===\n\n")

cat("--- 准确性（来自交叉验证）---\n")
cat("相关性:", round(cv_result$overall$mean_cor, 4), "\n")
cat("RMSE:", round(cv_result$overall$mean_rmse, 4), "\n")
cat("偏差斜率:", round(cv_result$overall$mean_bias_slope, 4), "\n\n")

cat("--- 方差组分（来自完整模型）---\n")
cat("遗传方差 (σ²g):", round(gblup_result$variances$sigma_g, 4), "\n")
cat("残差方差 (σ²e):", round(gblup_result$variances$sigma_e, 4), "\n")
cat("遗传力 (h²):", round(gblup_result$variances$h2, 4), "\n\n")

cat("--- 预测 ---\n")
cat("已观测个体数:", nrow(gblup_result$gebv_obs), "\n")
cat("所有个体数:", nrow(gblup_result$gebv_all), "\n")
```

---

## 贝叶斯方法

### 概述

贝叶斯方法使用马尔可夫链蒙特卡洛（MCMC）来估计SNP效应并预测育种值。`breedRplus` 通过BGLR包支持多种贝叶斯模型。

### 函数：`run_bglr()`

**用途**：拟合贝叶斯基因组预测模型。

**语法**：
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

**参数**：

| 参数 | 类型 | 描述 |
|------|------|------|
| `obj` | bigSNP | 来自 `load_plink()` 的bigSNP对象 |
| `pheno` | data.frame | 表型数据框 |
| `trait_name` | character | 性状列名称 |
| `id_col_name` | character | ID列名称 |
| `model_type` | character | 模型类型："BayesA"、"BayesB"、"BayesC"或"BL"（贝叶斯Lasso） |
| `n_iter` | integer | MCMC迭代次数（默认：5000） |
| `burn_in` | integer | 预烧迭代次数（默认：1000） |

**返回值**：包含以下内容的列表：
- `gebv`：包含GEBV的数据框
- `snp_effects`：包含SNP效应的数据框
- `residual_var`：残差方差估计
- `model_fit`：完整的BGLR模型对象

**模型类型**：

- **BayesA**：假设SNP效应服从t分布（适用于具有许多小效应的性状）
- **BayesB**：假设许多SNP效应为零（稀疏模型）
- **BayesC**：与BayesB类似但先验不同
- **BL**（贝叶斯Lasso）：使用L1正则化（类似于LASSO）

**示例1：BayesA模型**

```r
# 拟合BayesA模型
bayes_result <- run_bglr(
  obj = data$snp_obj,
  pheno = data$pheno,
  trait_name = "Yield",
  id_col_name = "ID",
  model_type = "BayesA",
  n_iter = 5000,
  burn_in = 1000
)

# 查看GEBV
head(bayes_result$gebv)

# 查看前几个SNP效应
head(bayes_result$snp_effects[order(abs(bayes_result$snp_effects$Effect), 
                                     decreasing = TRUE), ])
```

**示例2：贝叶斯Lasso**

```r
# 拟合贝叶斯Lasso模型
bayes_result <- run_bglr(
  obj = data$snp_obj,
  pheno = data$pheno,
  trait_name = "Yield",
  id_col_name = "ID",
  model_type = "BL",  # 贝叶斯Lasso
  n_iter = 10000,     # BL需要更多迭代
  burn_in = 2000
)
```

### 函数：`run_bglr_cv()`

**用途**：对贝叶斯模型进行k折交叉验证。

**语法**：
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

**参数**：与 `run_bglr()` 类似，另外：
- `geno`：数值型基因型矩阵（不是bigSNP对象）
- `k_folds`：CV折数
- `seed`：随机种子

**返回值**：包含以下内容的列表：
- `y_true`：观测值
- `y_pred`：预测值
- `cv_correlation`：Pearson相关系数
- `cv_mse`：均方误差

**示例：贝叶斯交叉验证**

```r
# 提取基因型矩阵
geno_matrix <- data$snp_obj$genotypes[, ]
rownames(geno_matrix) <- as.character(data$snp_obj$fam$sample.ID)

# 运行交叉验证
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

# 查看结果
print(bayes_cv$cv_correlation)
print(bayes_cv$cv_mse)
```

**绘制交叉验证结果：**

为了可视化观测值与预测值之间的关系：

```r
# 保存图为PNG文件（推荐，避免图形窗口大小问题）
png("bayesian_cv_predictions.png", width = 800, height = 600, res = 100)

# 设置图形边距
par(mar = c(4, 4, 3, 2))

# 创建散点图
plot(bayes_cv$y_true, bayes_cv$y_pred,
     xlab = "观测值", ylab = "预测值",
     main = "贝叶斯模型交叉验证预测",
     pch = 19, col = "blue", cex = 0.7)

# 添加1:1线（完美预测线）
abline(0, 1, col = "red", lty = 2, lwd = 2)

# 添加相关性文本
cor_val <- cor(bayes_cv$y_true, bayes_cv$y_pred, use = "complete.obs")
text(x = min(bayes_cv$y_true, na.rm = TRUE),
     y = max(bayes_cv$y_pred, na.rm = TRUE),
     labels = paste("r =", round(cor_val, 3)),
     pos = 4, col = "darkgreen", font = 2, cex = 1.2)

# 关闭PNG设备
dev.off()

cat("图已保存为 'bayesian_cv_predictions.png'\n")
```

**替代方法：直接绘制到屏幕**（如果图形窗口足够大）：

```r
# 关闭任何现有图形
graphics.off()

# 设置边距
par(mar = c(4, 4, 3, 2))

# 绘图
plot(bayes_cv$y_true, bayes_cv$y_pred,
     xlab = "观测值", ylab = "预测值",
     main = "贝叶斯模型交叉验证预测",
     pch = 19, col = "blue")
abline(0, 1, col = "red", lty = 2, lwd = 2)
```

---

## 机器学习方法

### 概述

`breedRplus` 支持多种机器学习算法进行基因组预测：
- **RF**：随机森林
- **XGB**：XGBoost（梯度提升）
- **MLP**：多层感知器（神经网络）
- **CNN**：一维卷积神经网络
- **PLS**：偏最小二乘回归

### **推荐：使用Python进行机器学习方法**

**我们强烈建议使用Python版本的 `breedrplus` 进行机器学习方法**，特别是神经网络（MLP、CNN）。原因如下：

#### **Python用于ML的优势：**

1. **安装更简单**：
   - TensorFlow安装简单：`pip install tensorflow`
   - 没有R→Python桥接的复杂性
   - 没有已弃用包的问题（R的 `keras` 已弃用，推荐使用 `keras3`）

2. **性能更好**：
   - 原生TensorFlow/Keras实现（无R包装器开销）
   - 对大型数据集的内存使用更高效
   - 对神经网络的GPU支持更好

3. **更多功能和更新**：
   - Python ML库开发更活跃
   - 可访问最新的TensorFlow/Keras功能
   - 文档和社区支持更好

4. **生产就绪**：
   - Python是ML生产环境的标准
   - 与ML管道集成更好
   - 更多部署选项

#### **何时使用R进行ML：**

- 您完全在R环境中工作，希望保持一致性
- 您只需要RF、XGB或PLS（无需TensorFlow依赖）
- 中小型数据集，性能差异可忽略

#### **混合方法（推荐）：**

- **使用R进行**：GBLUP、贝叶斯方法、RF、XGB、PLS
- **使用Python进行**：MLP、CNN（神经网络）

这样您可以两全其美：R中的统计方法，Python中的ML方法。

### ML方法和包

该包使用以下R包进行机器学习实现：

#### **1. 随机森林（RF）**
- **包**：`ranger` - 快速随机森林实现
- **方法**：带装袋的决策树集成
- **默认参数**：500棵树
- **使用场景**：良好的基线，能很好地处理非线性关系

#### **2. XGBoost（XGB）**
- **包**：`xgboost` - 极端梯度提升
- **方法**：基于树的梯度提升学习器
- **默认参数**：100轮，学习率（eta）= 0.1，max_depth = 6
- **使用场景**：通常达到高精度，处理复杂模式

#### **3. 多层感知器（MLP）**
- **包**：`keras` - 深度学习框架
- **方法**：3层前馈神经网络（64 → 32 → 1 单元）
- **架构**：
  - 输入层：所有SNP特征
  - 隐藏层1：64个单元，ReLU激活 + dropout（0.3）
  - 隐藏层2：32个单元，ReLU激活
  - 输出层：1个单元（连续预测）
- **默认参数**：50个epochs，batch_size = 32，validation_split = 0.1
- **使用场景**：捕获复杂的非线性关系

#### **4. 卷积神经网络（CNN）**
- **包**：`keras` - 深度学习框架
- **方法**：用于类序列基因型数据的1D卷积层
- **架构**：
  - Conv1D层：32个滤波器，kernel_size = 3
  - 最大池化：pool_size = 2
  - 展平层
  - 密集层：50个单元，ReLU
  - 输出层：1个单元
- **默认参数**：50个epochs，batch_size = 32
- **使用场景**：可捕获基因型序列中的局部模式

#### **5. 偏最小二乘（PLS）**
- **包**：`pls` - 偏最小二乘回归
- **方法**：降维和回归
- **默认参数**：最多50个成分（或特征数，取较小值）
- **使用场景**：处理高维数据，减少过拟合

**支持包**：
- `magrittr`：用于keras模型构建的管道操作符（`%>%`）
- `stats`：预测和相关性函数

### **Python实现（推荐）**

Python版本使用原生实现：
- **RF**：`scikit-learn` RandomForestRegressor
- **XGB**：原生 `xgboost` 包
- **MLP/CNN**：原生 `tensorflow/keras`
- **PLS**：`scikit-learn` PLSRegression

**安装**：
```bash
pip install numpy pandas scipy scikit-learn xgboost tensorflow bed-reader
```

**Python使用示例**：
```python
from breedrplus import load_plink, run_ml_model, run_ml_cv
import pandas as pd

# 加载数据
data = load_plink(
    bed_file="genotypes.bed",
    bim_file="genotypes.bim",
    fam_file="genotypes.fam",
    pheno=phenotype_df,
    id_col_name="ID"
)

# 提取基因型矩阵
geno_matrix = data['snp_obj']['genotypes']

# 训练随机森林模型
rf_result = run_ml_model(
    pheno=data['pheno'],
    geno=geno_matrix,
    trait_name="Yield",
    id_col_name="ID",
    model_type="RF",
    n_trees=500
)

# 查看预测
print(rf_result['gebv'].head())

# 交叉验证
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

# 查看准确性
print(f"Correlation: {rf_cv['cv_correlation']:.4f}")
print(f"MSE: {rf_cv['cv_mse']:.4f}")

# 训练神经网络（MLP）- 在Python中更容易！
mlp_result = run_ml_model(
    pheno=data['pheno'],
    geno=geno_matrix,
    trait_name="Yield",
    id_col_name="ID",
    model_type="MLP",
    epochs=50
)

# MLP交叉验证
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

**请参阅 `ML_PYTHON_GUIDE.md` 获取Python ML快速入门指南，或 `README_PYTHON.md` 和 `breedrplus-py/` 目录获取完整的Python文档。**

### 函数：`run_ml_model()`（R版本 - 也可用）

**注意**：虽然R实现可用，但我们建议使用Python进行ML方法（见上文）。

**用途**：训练机器学习模型并预测GEBV。

**语法**：
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

**参数**：

| 参数 | 类型 | 描述 |
|------|------|------|
| `pheno` | data.frame | 表型数据 |
| `geno` | matrix | 数值型基因型矩阵（行名=ID） |
| `trait_name` | character | 性状列名称 |
| `id_col_name` | character | ID列名称 |
| `model_type` | character | 模型类型："RF"、"XGB"、"MLP"、"CNN"或"PLS" |
| `...` | | 模型特定参数（见下文） |

**模型特定参数**：

- **RF**：`n_trees`（默认：500）- 树的数量
- **XGB**：`n_rounds`（默认：100）- 提升轮数
- **MLP/CNN**：`epochs`（默认：50）- 训练轮数
- **PLS**：`n_comp`（默认：min(50, ncol(geno))）- 成分数

**返回值**：包含以下内容的列表：
- `gebv`：包含预测GEBV的数据框
- `model_fit`：训练好的模型对象

**示例1：随机森林**

```r
# 提取基因型矩阵
geno_matrix <- data$snp_obj$genotypes[, ]
rownames(geno_matrix) <- as.character(data$snp_obj$fam$sample.ID)

# 训练随机森林模型
rf_result <- run_ml_model(
  pheno = data$pheno,
  geno = geno_matrix,
  trait_name = "Yield",
  id_col_name = "ID",
  model_type = "RF",
  n_trees = 500
)

# 查看预测
head(rf_result$gebv)
```

**示例2：XGBoost**

```r
# 训练XGBoost模型
xgb_result <- run_ml_model(
  pheno = data$pheno,
  geno = geno_matrix,
  trait_name = "Yield",
  id_col_name = "ID",
  model_type = "XGB",
  n_rounds = 100
)
```

**示例3：神经网络（MLP）**

```r
# 训练多层感知器
mlp_result <- run_ml_model(
  pheno = data$pheno,
  geno = geno_matrix,
  trait_name = "Yield",
  id_col_name = "ID",
  model_type = "MLP",
  epochs = 50
)
```

**示例4：卷积神经网络**

```r
# 训练一维CNN
cnn_result <- run_ml_model(
  pheno = data$pheno,
  geno = geno_matrix,
  trait_name = "Yield",
  id_col_name = "ID",
  model_type = "CNN",
  epochs = 50
)
```

**示例5：偏最小二乘**

```r
# 训练PLS模型
pls_result <- run_ml_model(
  pheno = data$pheno,
  geno = geno_matrix,
  trait_name = "Yield",
  id_col_name = "ID",
  model_type = "PLS",
  n_comp = 50
)
```

### 函数：`run_ml_cv()`（R版本 - 也可用）

**注意**：虽然R实现可用，但我们建议使用Python进行ML方法（见上文）。

**用途**：对机器学习模型进行k折交叉验证。

**语法**：
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

**参数**：与 `run_ml_model()` 类似，另外：
- `k`：CV折数
- `seed`：随机种子

**返回值**：包含以下内容的列表：
- `cv_gebv`：包含CV预测的数据框
- `cv_correlation`：Pearson相关系数
- `cv_mse`：均方误差

**示例：比较多种ML模型**

```r
# 提取基因型矩阵
geno_matrix <- data$snp_obj$genotypes[, ]
rownames(geno_matrix) <- as.character(data$snp_obj$fam$sample.ID)

# 比较模型
models <- c("RF", "XGB", "MLP", "PLS")
results <- list()

for (model in models) {
  cat("运行", model, "...\n")
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

# 比较结果
comparison <- data.frame(
  Model = models,
  Correlation = sapply(results, function(x) x$cv_correlation),
  MSE = sapply(results, function(x) x$cv_mse)
)

print(comparison)
```

### **绘制交叉验证结果**

**方法1：保存为PNG文件（推荐）**

```r
# 提取观测值
y_true <- data$pheno[[trait_name]][match(ml_cv$cv_gebv$ID, data$pheno[[id_col_name]])]
y_pred <- ml_cv$cv_gebv$gebv

# 保存图为PNG文件
png("ml_cv_predictions.png", width = 800, height = 600, res = 100)

# 设置图形边距
par(mar = c(4, 4, 3, 2))

# 创建散点图
plot(y_true, y_pred,
     xlab = "观测值", ylab = "预测值",
     main = "ML模型交叉验证预测（随机森林）",
     pch = 19, col = "blue", cex = 0.7)

# 添加1:1线（完美预测线）
abline(0, 1, col = "red", lty = 2, lwd = 2)

# 添加相关性文本
cor_val <- cor(y_true, y_pred, use = "complete.obs")
text(x = min(y_true, na.rm = TRUE),
     y = max(y_pred, na.rm = TRUE),
     labels = paste("r =", round(cor_val, 3)),
     pos = 4, col = "darkgreen", font = 2, cex = 1.2)

# 关闭PNG设备
dev.off()

cat("图已保存为 'ml_cv_predictions.png'\n")
cat("相关性:", round(ml_cv$cv_correlation, 4), "\n")
cat("MSE:", round(ml_cv$cv_mse, 4), "\n")
```

**方法2：直接绘制到屏幕**（如果图形窗口足够大）：

```r
# 关闭任何现有图形
graphics.off()

# 设置边距
par(mar = c(4, 4, 3, 2))

# 提取值
y_true <- data$pheno[[trait_name]][match(ml_cv$cv_gebv$ID, data$pheno[[id_col_name]])]
y_pred <- ml_cv$cv_gebv$gebv

# 绘图
plot(y_true, y_pred,
     xlab = "观测值", ylab = "预测值",
     main = "ML模型交叉验证预测",
     pch = 19, col = "blue")
abline(0, 1, col = "red", lty = 2, lwd = 2)
```

### **R用户的重要提示：**

1. **TensorFlow安装**：如果在R中使用MLP或CNN，您需要安装TensorFlow：
   ```r
   # 安装keras3（推荐）或keras（已弃用）
   install.packages("keras3")
   library(keras3)
   install_keras()  # 安装TensorFlow
   ```

2. **性能**：由于R→Python桥接开销，神经网络（MLP、CNN）在R中会更慢

3. **兼容性**：`keras` 包已弃用；如果可用，请使用 `keras3`

4. **推荐**：对于神经网络，强烈建议使用Python版本

---

## 完整工作流程示例

### 示例1：完整GBLUP工作流程

```r
library(breedRplus)
library(readxl)

# 步骤1：加载数据
pheno <- read_excel("phenotypes.xlsx")
data <- load_plink(
  bed_file = "genotypes.bed",
  bim_file = "genotypes.bim",
  fam_file = "genotypes.fam",
  pheno = pheno,
  id_col_name = "AnimalID"
)

# 步骤2：运行带固定效应的GBLUP
gblup_result <- run_gblup(
  qc_results = data,
  trait_name = "MilkYield",
  id_col_name = "AnimalID",
  fixed_effects = ~ Herd + Parity,
  predict_all = TRUE
)

# 步骤3：交叉验证
cv_result <- cv_gblup(
  qc_results = data,
  trait_name = "MilkYield",
  id_col_name = "AnimalID",
  fixed_effects = ~ Herd + Parity,
  k = 10,
  seed = 2025
)

# 步骤4：评估并保存结果
print("方差组分：")
print(gblup_result$variances)

print("交叉验证结果：")
print(cv_result$overall)

# 保存GEBV
write.csv(gblup_result$gebv_all, "GEBVs_all_animals.csv", row.names = FALSE)
```

### 示例2：比较多种方法

```r
library(breedRplus)

# 加载数据
data <- load_plink(
  bed_file = "genotypes.bed",
  bim_file = "genotypes.bim",
  fam_file = "genotypes.fam",
  pheno = read.csv("phenotypes.csv"),
  id_col_name = "ID"
)

# 提取基因型矩阵用于ML方法
geno_matrix <- data$snp_obj$genotypes[, ]
rownames(geno_matrix) <- as.character(data$snp_obj$fam$sample.ID)

# 比较方法
methods <- list()

# 1. GBLUP
cat("运行GBLUP...\n")
methods$GBLUP <- cv_gblup(
  qc_results = data,
  trait_name = "Yield",
  id_col_name = "ID",
  k = 5
)

# 2. BayesA
cat("运行BayesA...\n")
methods$BayesA <- run_bglr_cv(
  geno = geno_matrix,
  pheno = data$pheno,
  trait_name = "Yield",
  id_col_name = "ID",
  model_type = "BayesA",
  k_folds = 5
)

# 3. 随机森林
cat("运行随机森林...\n")
methods$RF <- run_ml_cv(
  pheno = data$pheno,
  geno = geno_matrix,
  trait_name = "Yield",
  id_col_name = "ID",
  model_type = "RF",
  k = 5
)

# 4. XGBoost
cat("运行XGBoost...\n")
methods$XGB <- run_ml_cv(
  pheno = data$pheno,
  geno = geno_matrix,
  trait_name = "Yield",
  id_col_name = "ID",
  model_type = "XGB",
  k = 5
)

# 比较结果
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

### 示例3：选择候选个体的预测

```r
library(breedRplus)

# 加载数据（包括训练和选择候选个体）
data <- load_plink(
  bed_file = "all_genotypes.bed",
  bim_file = "all_genotypes.bim",
  fam_file = "all_genotypes.fam",
  pheno = read.csv("training_phenotypes.csv"),  # 仅训练动物有表型
  id_col_name = "ID"
)

# 运行GBLUP并预测所有个体
gblup_result <- run_gblup(
  qc_results = data,
  trait_name = "Yield",
  id_col_name = "ID",
  predict_all = TRUE  # 预测未表型个体
)

# 识别选择候选个体（没有表型的个体）
training_ids <- data$pheno$ID
all_ids <- data$snp_obj$fam$sample.ID
candidate_ids <- setdiff(all_ids, training_ids)

# 提取选择候选个体的GEBV
candidate_gebvs <- gblup_result$gebv_all[
  gblup_result$gebv_all$ID %in% candidate_ids,
]

# 排序候选个体
candidate_gebvs <- candidate_gebvs[order(candidate_gebvs$gebv, decreasing = TRUE), ]

# 选择前10%
top_10_percent <- head(candidate_gebvs, n = round(0.1 * nrow(candidate_gebvs)))

print("前10%选择候选个体：")
print(top_10_percent)
```

---

## 故障排除

### 常见问题和解决方案

#### 问题1："文件未找到"错误

**问题**：`load_plink()` 无法找到PLINK文件。

**解决方案**：
```r
# 检查文件路径
file.exists("genotypes.bed")  # 应返回TRUE

# 如需要，使用绝对路径
bed_file <- "/完整路径/to/genotypes.bed"
```

#### 问题2：基因型和表型之间的ID不匹配

**问题**：表型ID与基因型样本ID不匹配。

**解决方案**：
```r
# 检查ID
geno_ids <- data$snp_obj$fam$sample.ID
pheno_ids <- data$pheno$ID

# 查找不匹配
setdiff(geno_ids, pheno_ids)
setdiff(pheno_ids, geno_ids)

# 确保ID列名正确
data <- load_plink(..., id_col_name = "正确的列名")
```

#### 问题3：大型数据集的内存问题

**问题**：处理大型基因型文件时内存不足。

**解决方案**：
- 高效使用 `bigsnpr`（已实现）
- 考虑在分析前对SNP进行子集化
- 增加可用RAM
- 对于非常大的数据集，使用计算集群

#### 问题4：交叉验证耗时过长

**问题**：CV很慢，特别是对于贝叶斯方法。

**解决方案**：
```r
# 减少折数
cv_result <- cv_gblup(..., k = 3)  # 而不是 k = 10

# 减少贝叶斯方法的迭代次数
bayes_cv <- run_bglr_cv(..., n_iter = 2000, burn_in = 500)

# 使用更快的方法进行初步筛选
rf_cv <- run_ml_cv(..., model_type = "RF")  # 快速
```

#### 问题5：Keras/神经网络错误

**问题**：运行MLP或CNN模型时出错。

**解决方案**：
- 确保已安装TensorFlow：`install.packages("tensorflow")` 和 `tensorflow::install_tensorflow()`
- 检查Python是否可用：`reticulate::py_available()`
- 尝试减少轮数或特征数

#### 问题6：预测准确性低

**问题**：CV相关系数低（< 0.3）。

**可能原因**：
- 训练群体规模小
- 低遗传力性状
- 基因型-表型对齐差
- 需要质量控制（删除低质量SNP）

**解决方案**：
- 增加训练群体规模
- 检查数据质量
- 尝试不同模型
- 考虑包含固定效应

---

## 最佳实践

### 1. 数据质量控制

- **分析前**：检查缺失数据、异常值和数据质量
- **基因型QC**：考虑按次要等位基因频率（MAF）、检出率等过滤SNP
- **表型QC**：删除明显异常值，检查数据录入错误

### 2. 模型选择

- **从GBLUP开始**：它快速且提供良好的基线
- **尝试多种方法**：不同方法对不同性状效果更好
- **使用交叉验证**：始终验证您的模型

### 3. 交叉验证策略

- **使用适当的k**：5-10折通常足够
- **尽可能分层**：如果您有家系结构，使用 `stratify_by`
- **设置随机种子**：为了可重复性
- **报告折级和整体指标**

### 4. 计算效率

- **对于大型数据集**：从更快的方法开始（GBLUP、RF）
- **对于贝叶斯方法**：从较少迭代开始测试，然后增加
- **并行处理**：考虑并行化CV折（未来功能）

### 5. 结果解释

- **遗传力（h²）**：应该对您的性状在生物学上合理
- **CV相关系数**：越高越好，但取决于性状遗传力
- **偏差**：斜率接近1表示预测无偏
- **比较方法**：使用多个指标（相关系数、MSE、偏差）

### 6. 报告

始终报告：
- 样本量（训练和验证）
- SNP数量
- 模型参数（迭代次数、折数等）
- 方差组分（对于GBLUP）
- 交叉验证指标
- 预测准确性

---

## 参考文献

### 关键论文

- **GBLUP**：VanRaden, P. M. (2008). Efficient methods to compute genomic predictions. *Journal of dairy science*, 91(11), 4414-4423.

- **贝叶斯方法**：de los Campos, G., et al. (2013). Whole-genome regression and prediction methods applied to plant and animal breeding. *Genetics*, 193(2), 327-345.

- **机器学习**：González-Recio, O., et al. (2014). Machine learning methods and predictive ability metrics for genome-wide prediction of complex traits. *Livestock Science*, 166, 217-231.

### 包依赖

- `bigsnpr`：Privé, F., et al. (2018). Efficient analysis of large-scale genome-wide data with two R packages: bigstatsr and bigsnpr. *Bioinformatics*, 34(16), 2781-2787.

- `BGLR`：Pérez, P., & de los Campos, G. (2014). Genome-wide regression and prediction with the BGLR statistical package. *Genetics*, 198(2), 483-495.

---

## 获取帮助

### 文档

- 包帮助：`?function_name`（例如，`?run_gblup`）
- 小插图：`browseVignettes("breedRplus")`
- 示例工作流程：`system.file("example_workflow.R", package = "breedRplus")`

### 支持

- GitHub Issues
- 电子邮件：[higefee@gmail.com]

### 贡献

欢迎贡献！请查看仓库中的贡献指南。

---

## 附录：函数快速参考

| 函数 | 用途 | 关键参数 |
|------|------|----------|
| `load_plink()` | 加载PLINK文件 | `bed_file`, `bim_file`, `fam_file`, `pheno` |
| `run_gblup()` | 拟合GBLUP模型 | `qc_results`, `trait_name`, `predict_all` |
| `cv_gblup()` | GBLUP交叉验证 | `qc_results`, `trait_name`, `k` |
| `run_bglr()` | 拟合贝叶斯模型 | `obj`, `pheno`, `model_type` |
| `run_bglr_cv()` | 贝叶斯交叉验证 | `geno`, `pheno`, `model_type`, `k_folds` |
| `run_ml_model()` | 训练ML模型 | `pheno`, `geno`, `model_type` |
| `run_ml_cv()` | ML交叉验证 | `pheno`, `geno`, `model_type`, `k` |

---

**最后更新**：2025

**包版本**：0.0.0.9000

