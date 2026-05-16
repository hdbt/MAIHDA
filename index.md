# MAIHDA: Multilevel Analysis of Individual Heterogeneity and Discriminatory Accuracy

[![R-CMD-check](https://github.com/hdbt/MAIHDA/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/hdbt/MAIHDA/actions/workflows/R-CMD-check.yaml)
[![R](https://img.shields.io/badge/R-%3E%3D4.1.0-blue)](https://www.r-project.org/)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Codecov test
coverage](https://codecov.io/gh/hdbt/MAIHDA/branch/main/graph/badge.svg)](https://app.codecov.io/gh/hdbt/MAIHDA?branch=main)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

## Overview

The MAIHDA package provides a comprehensive toolkit for conducting
Multilevel Analysis of Individual Heterogeneity and Discriminatory
Accuracy (MAIHDA). This approach is particularly valuable for examining
intersectional inequalities in health and social outcomes by considering
the joint effects of multiple social categories (e.g., gender, race,
socioeconomic status).

## Key Features

- **Create Intersectional Strata**: Automatically generate strata from
  multiple categorical variables
- **[Interactive Dashboard](https://hdbt.shinyapps.io/shiny/)**

: A fully-featured Shiny application
([`run_maihda_app()`](https://hdbt.github.io/MAIHDA/reference/run_maihda_app.md))
for no-code exploratory data analysis and model fitting - **Model
Fitting**: Support for both lme4 and brms (Bayesian) engines -
**Summaries & Decompositions**: Variance partition coefficients
(VPC/ICC), stratum-specific estimates, and stepwise Proportional Change
in Variance (PCV) - **Multiple Prediction Types**: Individual-level and
stratum-level predictions - **Visualizations**: Predicted stratum
values, VPC visualizations, risk/effect diagnostics, and observed
vs. shrunken estimates - **Model Comparison**: Compare models with
robust bootstrap confidence intervals for VPC/ICC - **Proportional
Change in Variance (PVC)**: Quantify how much between-stratum variance
is explained by additional predictors

## Installation

You can install the development version from GitHub:

``` r

# install.packages("devtools")
devtools::install_github("hdbt/MAIHDA")
```

or the current stable version from CRAN:

``` r

install.packages("MAIHDA")
```

## Quick Start

``` r

library(MAIHDA)
data("maihda_health_data")

# 1. Fit a MAIHDA model
model <- fit_maihda(
  BMI ~ Age + Poverty + (1 | Gender:Race:Education),
  data = maihda_health_data
)

# 2. Summarize the model (Variance Partition Coefficient, stratum estimates)
summary(model)

# 3. Visualize results (Predicted outputs, VPC, Risk vs. Effect, etc.)
plot(model)
```

## Main Functions

### `make_strata()`

Creates intersectional strata from multiple categorical variables with
optional minimum count filtering.

### `fit_maihda()`

Fits multilevel models using either lme4 (default) or brms engine.
Supports various families including gaussian, binomial, and poisson.

### `summary()`

Provides comprehensive model summaries including:

- Variance Partition Coefficient (VPC/ICC)
- Variance components decomposition
- Stratum-specific random effects
- Optional bootstrap confidence intervals

### `predict_maihda()`

Makes predictions at two levels:

- **Individual**: Full predictions including random effects
- **Strata**: Stratum-specific random effects with uncertainty

### `plot()`

Creates various visualizations:

- **Predicted**: Predicted stratum values with confidence intervals
- **VPC plots**: Visualizes variance partitioning
- **Observed vs. Shrunken**: Shows shrinkage of stratum estimates
- **Risk vs Effect**: Visualizes baseline strata risk against
  intersectional effect size
- **Effect Decomp**: Additive versus specific interaction decompositions

### `maihda_ternary_plot()`

Generates an advanced ternary plot showing the relative contribution of
intercept, additive effects, and unique intersectional interaction
effects to the overall outcome for each stratum.

### `plot_prediction_deviation_panels()`

Creates an advanced, publication-ready two-panel dashboard for
visualizing predicted values and identifying extreme/deviant cases in
individual predictions.

### `compare_maihda()`

Compares VPC/ICC across multiple models with optional bootstrap
confidence intervals.

### `calculate_pvc()`

Calculates the proportional change in between-stratum variance (PVC)
between two models. This measures how much of the between-stratum
variance from a baseline model is explained (or changed) by adding
additional predictors in a second model:- Formula: PVC = (Var_model1 -
Var_model2) / Var_model1 - Works with both lme4 and brms engines -
Supports bootstrap confidence intervals for lme4 models

### `stepwise_pcv()`

Evaluates multiple sequential models by iteratively adding covariates
step-by-step to quantify precisely which variables explain the
structural inequalities.

### `run_maihda_app()`

Launches a locally-hosted, interactive Shiny Dashboard that exposes the
core functionalities for data modeling, visualization, and summarization
visually.

## Example: Intersectional Health Inequalities

``` r

# Fit model adjusting for age, automatically creating strata from gender, race, and education
model <- fit_maihda(
  BMI ~ Age + (1 | Gender:Race:Education),
  data = maihda_health_data
)

# Get variance partition coefficient
summary <- summary(model, bootstrap = TRUE, n_boot = 1000)

# VPC of 0.15 means 15% of variance is between strata
# This indicates substantial intersectional inequality

# Visualize which strata have higher/lower outcomes using new advanced plots
plot(model, type = "predicted")
plot(model, type = "risk_vs_effect")

# Map out where the specific intersectional variance is emerging
plot(model, type = "ternary")

# Or run it with no type to see them all!
# plot(model)
```

## Using brms for Bayesian Inference

``` r

# Requires brms package
model_brms <- fit_maihda(
  BMI ~ Age + (1 | Gender:Race:Education),
  data = maihda_health_data,
  engine = "brms",
  chains = 4,
  iter = 2000
)

summary_brms <- summary(model_brms)
```

## Model Comparison with Bootstrap

``` r

# Fit competing models
model1 <- fit_maihda(outcome ~ age + (1 | gender:race), data = data1)
model2 <- fit_maihda(outcome ~ age + gender + (1 | gender:race), data = data2)

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

``` r

# Fit baseline model
model1 <- fit_maihda(outcome ~ age + (1 | gender:race), data = data)

# Fit model with additional predictor
model2 <- fit_maihda(outcome ~ age + gender + (1 | gender:race), data = data)

# Calculate PVC without bootstrap
pvc_result <- calculate_pvc(model1, model2)
print(pvc_result)

# Calculate PVC with bootstrap confidence intervals
pvc_boot <- calculate_pvc(model1, model2, bootstrap = TRUE, n_boot = 1000)
print(pvc_boot)

# Interpretation: A PVC of 0.25 means that model2 explains 25% of the
# between-stratum variance that was present in model1
```

## Stepwise Proportional Change in Variance (PCV)

``` r

# Incrementally track the explained variables interactively
stepwise_results <- stepwise_pcv(
  data = data,
  outcome = "health_outcome",
  vars = c("gender", "race", "age")
)
print(stepwise_results)
```

## Interactive Shiny App

You can access a live, cloud-hosted version of the MAIHDA interactive
dashboard directly in your browser without installing R:
**<https://hdbt.shinyapps.io/shiny/>**

Alternatively, you can run all analyses described above in the browser
locally by typing:

``` r

library(MAIHDA)
run_maihda_app()
```

## Documentation

For detailed documentation and examples, see the package vignette:

``` r

vignette("introduction", package = "MAIHDA")
```

## Dependencies

**Required:**

- R (\>= 4.1.0)
- lme4 (\>= 1.1-27)
- ggplot2 (\>= 3.4.0)
- dplyr (\>= 1.0.0)
- tidyr (\>= 1.1.0)
- stats, methods, tibble, rlang

**Optional:**

- brms (\>= 2.15.0) - for Bayesian models
- boot (\>= 1.3-20) - for bootstrap confidence intervals

## References

- Evans, C. R., Williams, D. R., Onnela, J. P., & Subramanian, S. V.
  (2018). A multilevel approach to modeling health inequalities at the
  intersection of multiple social identities. *Social Science &
  Medicine*, 203, 64-73.

- Merlo, J. (2018). Multilevel analysis of individual heterogeneity and
  discriminatory accuracy (MAIHDA) within an intersectional framework.
  *Social Science & Medicine*, 203, 74-80.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the
[LICENSE](https://hdbt.github.io/MAIHDA/LICENSE) file for details.

## Citation

If you use this package in your research, please cite:

``` text
Bulut (2025). *MAIHDA: Multilevel Analysis of Individual Heterogeneity and Discriminatory Accuracy.* R package version 0.1.9, https://github.com/hdbt/MAIHDA. doi: 10.32614/CRAN.package.MAIHDA
```

A BibTeX entry for LaTeX users is:

``` bibtex
@Manual{Bulut2025MAIHDA,
  title  = {MAIHDA: Multilevel Analysis of Individual Heterogeneity and Discriminatory Accuracy},
  author = {Hamid Bulut},
  year   = {2025},
  note   = {R package version 0.1.9},
  url    = {https://github.com/hdbt/MAIHDA},
  doi    = {10.32614/CRAN.package.MAIHDA}
}
```
