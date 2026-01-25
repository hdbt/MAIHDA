# MAIHDA: Multilevel Analysis of Individual Heterogeneity and Discriminatory Accuracy
[![R-CMD-check](https://github.com/hdbt/MAIHDA/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/hdbt/MAIHDA/actions/workflows/R-CMD-check.yaml)
[![R](https://img.shields.io/badge/R-%3E%3D3.5.0-blue)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Codecov test coverage](https://codecov.io/gh/hdbt/MAIHDA/branch/main/graph/badge.svg)](https://app.codecov.io/gh/hdbt/MAIHDA?branch=main)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

## Overview

The MAIHDA package provides a comprehensive toolkit for conducting Multilevel Analysis of Individual Heterogeneity and Discriminatory Accuracy (MAIHDA). This approach is particularly valuable for examining intersectional inequalities in health and social outcomes by considering the joint effects of multiple social categories (e.g., gender, race, socioeconomic status).

## Key Features

- **Create Intersectional Strata**: Automatically generate strata from multiple categorical variables
- **Flexible Model Fitting**: Support for both lme4 (frequentist) and brms (Bayesian) engines
- **Comprehensive Summaries**: Variance partition coefficients (VPC/ICC) and stratum-specific estimates
- **Multiple Prediction Types**: Individual-level and stratum-level predictions
- **Rich Visualizations**: Caterpillar plots, VPC visualizations, and observed vs. shrunken estimates
- **Model Comparison**: Compare models with bootstrap confidence intervals for VPC/ICC
- **Proportional Change in Variance (PVC)**: Quantify how much between-stratum variance is explained by additional predictors

## Installation

You can install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("hdbt/MAIHDA")
```

## Quick Start

```r
library(MAIHDA)

# 1. Create intersectional strata
data <- data.frame(
  gender = sample(c("Male", "Female"), 1000, replace = TRUE),
  race = sample(c("White", "Black", "Hispanic"), 1000, replace = TRUE),
  age = rnorm(1000, 50, 10),
  health_outcome = rnorm(1000, 100, 15)
)

strata_result <- make_strata(data, vars = c("gender", "race"))

# 2. Fit MAIHDA model
model <- fit_maihda(
  health_outcome ~ age + (1 | stratum),
  data = strata_result$data,
  engine = "lme4"
)

# 3. Summarize with variance partition
summary_result <- summary_maihda(model, bootstrap = TRUE, n_boot = 500)
print(summary_result)

# 4. Make predictions
pred_strata <- predict_maihda(model, type = "strata")
pred_individual <- predict_maihda(model, type = "individual")

# 5. Visualize results
plot_maihda(model, type = "caterpillar")
plot_maihda(model, type = "vpc")
plot_maihda(model, type = "obs_vs_shrunken")
```

## Main Functions

### `make_strata()`
Creates intersectional strata from multiple categorical variables with optional minimum count filtering.

### `fit_maihda()`
Fits multilevel models using either lme4 (default) or brms engine. Supports various families including gaussian, binomial, and poisson.

### `summary_maihda()`
Provides comprehensive model summaries including:
- Variance Partition Coefficient (VPC/ICC)
- Variance components decomposition
- Stratum-specific random effects
- Optional bootstrap confidence intervals

### `predict_maihda()`
Makes predictions at two levels:
- **Individual**: Full predictions including random effects
- **Strata**: Stratum-specific random effects with uncertainty

### `plot_maihda()`
Creates various visualizations:
- **Caterpillar plots**: Displays stratum random effects with confidence intervals
- **VPC plots**: Visualizes variance partitioning
- **Observed vs. Shrunken**: Shows shrinkage of stratum estimates

### `compare_maihda()`
Compares VPC/ICC across multiple models with optional bootstrap confidence intervals.

### `calculate_pvc()`
Calculates the proportional change in between-stratum variance (PVC) between two models. This measures how much of the between-stratum variance from a baseline model is explained (or changed) by adding additional predictors in a second model:
- Formula: PVC = (Var_model1 - Var_model2) / Var_model1
- Supports bootstrap confidence intervals
- Works with both lme4 and brms engines

## Example: Intersectional Health Inequalities

```r
# Create strata from gender and race
strata_result <- make_strata(health_data, vars = c("gender", "race", "education"))

# Fit model adjusting for age
model <- fit_maihda(
  health_outcome ~ age + (1 | stratum),
  data = strata_result$data
)

# Get variance partition coefficient
summary <- summary_maihda(model, bootstrap = TRUE, n_boot = 1000)

# VPC of 0.15 means 15% of variance is between strata
# This indicates substantial intersectional inequality

# Visualize which strata have higher/lower outcomes
plot_maihda(model, type = "caterpillar")
```

## Using brms for Bayesian Inference

```r
# Requires brms package
model_brms <- fit_maihda(
  health_outcome ~ age + (1 | stratum),
  data = strata_result$data,
  engine = "brms",
  chains = 4,
  iter = 2000
)

summary_brms <- summary_maihda(model_brms)
```

## Model Comparison with Bootstrap

```r
# Fit competing models
model1 <- fit_maihda(outcome ~ age + (1 | stratum), data = data1)
model2 <- fit_maihda(outcome ~ age + gender + (1 | stratum), data = data2)

# Compare with bootstrap CI
comparison <- compare_maihda(
  model1, model2,
  model_names = c("Base", "With Gender"),
  bootstrap = TRUE,
  n_boot = 1000
)

# Visualize comparison
plot_comparison(comparison)
```

## Calculating Proportional Change in Variance (PVC)

```r
# Fit baseline model
model1 <- fit_maihda(outcome ~ age + (1 | stratum), data = data)

# Fit model with additional predictor
model2 <- fit_maihda(outcome ~ age + gender + (1 | stratum), data = data)

# Calculate PVC without bootstrap
pvc_result <- calculate_pvc(model1, model2)
print(pvc_result)

# Calculate PVC with bootstrap confidence intervals
pvc_boot <- calculate_pvc(model1, model2, bootstrap = TRUE, n_boot = 1000)
print(pvc_boot)

# Interpretation: A PVC of 0.25 means that model2 explains 25% of the 
# between-stratum variance that was present in model1
```

## Documentation

For detailed documentation and examples, see the package vignette:

```r
vignette("introduction", package = "MAIHDA")
```

## Dependencies

**Required:**
- R (>= 3.5.0)
- lme4 (>= 1.1-27)
- ggplot2 (>= 3.3.0)
- dplyr (>= 1.0.0)
- tidyr (>= 1.1.0)
- stats, methods, tibble, rlang

**Optional:**
- brms (>= 2.15.0) - for Bayesian models
- boot (>= 1.3-20) - for bootstrap confidence intervals

## References

- Evans, C. R., Williams, D. R., Onnela, J. P., & Subramanian, S. V. (2018). A multilevel approach to modeling health inequalities at the intersection of multiple social identities. *Social Science & Medicine*, 203, 64-73.

- Merlo, J. (2018). Multilevel analysis of individual heterogeneity and discriminatory accuracy (MAIHDA) within an intersectional framework. *Social Science & Medicine*, 203, 74-80.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use this package in your research, please cite:

```
  Bulut (2025). MAIHDA: Multilevel Analysis of Individual Heterogeneity and Discriminatory Accuracy
R package version 0.1.3
  https://github.com/hdbt/MAIHDA

A BibTeX entry for LaTeX users is:

@Manual{Bulut2025MAIHDA,
  title  = {MAIHDA: Multilevel Analysis of Individual Heterogeneity and Discriminatory Accuracy},
  author = {Hamid Bulut},
  year   = {2025},
  note   = {R package version 0.1.3},
  url    = {https://github.com/hdbt/MAIHDA}
}

```
