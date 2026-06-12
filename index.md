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

- **One-call Workflow**:
  [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md) fits
  the null *and* adjusted models, summarises the VPC/ICC and the PCV
  decomposition (additive vs. intersectional), and (optionally) compares
  across a higher-level group in a single call
- **Create Intersectional Strata**: Automatically generate strata from
  multiple categorical variables
- **[Interactive Dashboard](https://hdbt.shinyapps.io/shiny/)**: A
  fully-featured Shiny application
  ([`run_maihda_app()`](https://hdbt.github.io/MAIHDA/reference/run_maihda_app.md))
  for no-code exploratory data analysis and model fitting
- **Model Fitting**: Support for both lme4 and brms (Bayesian) engines
- **Design-Weighted MAIHDA**: Survey/sampling weights via
  `sampling_weights` – weighted pseudo-maximum-likelihood through the
  `WeMix` engine (or likelihood weights under brms) for complex-survey
  data such as NHANES or PISA, with design-consistent fixed-effect
  standard errors and a design-weighted VPC, PCV, stratum summaries and
  AUC
- **Summaries & Decompositions**: Variance partition coefficients
  (VPC/ICC), stratum-specific estimates, and stepwise Proportional
  Change in Variance (PCV)
- **Multiple Prediction Types**: Individual-level and stratum-level
  predictions
- **Visualizations**: Predicted stratum values, VPC visualizations,
  mean-prediction vs. stratum-effect diagnostics, and observed
  vs. shrunken estimates
- **Model Comparison**: Compare models with parametric-bootstrap
  confidence intervals for VPC/ICC
- **Group Comparison**:
  [`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
  contrasts intersectional inequality (VPC/ICC) across levels of a
  higher-level variable such as country or region
- **Contextual Cross-Classified MAIHDA**: `context = "school"` crosses
  the intersectional strata with a place/institution level (school,
  hospital, region) in one model and partitions the unexplained variance
  into between-stratum vs. between-context vs. residual – the
  cross-classified MAIHDA of the literature
- **Proportional Change in Variance (PCV)**: Quantify the proportional
  change in between-stratum variance when additional predictors are
  added (it is variance “explained” only when the models are nested;
  otherwise it is a model-dependent change), or read the
  additive/interaction split off a single model with
  `decomposition = "crossed-dimensions"`

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

> **Note:** the CRAN release can lag behind this repository. Some
> features documented below
> (e.g. [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md)
> and
> [`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md))
> were added after the current CRAN version, so install the development
> version from GitHub if you need them. \## Quick Start

``` r

library(MAIHDA)
data("maihda_health_data")

# Everything in one call: null + adjusted fit, VPC/ICC summary, and PCV decomposition.
# Write the strata variables as additive fixed effects -- the lme4-native, fully-
# specified ADJUSTED model. maihda() uses that as the adjusted model and derives the
# NULL by dropping those main effects. (Omit them and maihda() adds them for you,
# printing a message; the result is identical.)
analysis <- maihda(BMI ~ Age + Gender + Race + Education + (1 | Gender:Race:Education),
                   data = maihda_health_data)
analysis                       # VPC/ICC (null) and PCV (null -> adjusted)
analysis$formula               # null:     BMI ~ Age + (1 | stratum)
analysis$adjusted_formula      # adjusted: BMI ~ Age + Gender + Race + Education + (1 | stratum)
summary(analysis)              # variance components
analysis$pcv                   # proportional change in between-stratum variance

# plot() routes each view to the right model: the VPC/shrinkage views use the null
# model; the additive-vs-intersectional views use the adjusted model.
plot(analysis, type = "vpc")
plot(analysis, type = "effect_decomp")
```

Prefer the individual building blocks? They are all still available:

``` r

# 1. Fit a MAIHDA model
model <- fit_maihda(
  BMI ~ Age + Poverty + (1 | Gender:Race:Education),
  data = maihda_health_data
)

# 2. Summarize the model (Variance Partition Coefficient, stratum estimates)
summary(model)

# 3. Visualize results (Predicted outputs, VPC, Mean Prediction vs. Stratum Effect, etc.)
plot(model)
```

## Main Functions

### `maihda()`

A single high-level entry point that runs the standard two-model MAIHDA
workflow: it fits the **null** model (covariates only) and the
**adjusted** model (plus the dimensions’ additive main effects – write
them in the formula, or let
[`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md) add them
with a message), summarises the null-model VPC/ICC, and reports the
**PCV** (the additive share of the intersectional inequality). When a
`group` is supplied it also runs this decomposition within each group.
Alternatively, `decomposition = "crossed-dimensions"` reads the
additive/interaction split off a *single* model that enters each
dimension’s main effect as a random intercept. Returns one
`maihda_analysis` object with
[`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html), and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods
([`plot()`](https://rdrr.io/r/graphics/plot.default.html) routes the
VPC/shrinkage views to the null model and the additive-vs-intersectional
views to the adjusted model). It is intrinsically a decomposition and
has no single-model mode – use
[`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
for a single fit.

### `make_strata()`

Creates intersectional strata from multiple categorical variables with
optional minimum count filtering.

### `fit_maihda()`

Fits multilevel models using the lme4 (default), brms, or WeMix
(design-weighted, via `sampling_weights`) engine. Supports various
families including gaussian, binomial, and poisson.

### Contextual cross-classified MAIHDA (`context =`)

The MAIHDA literature’s *cross-classified* design crosses individuals’
intersectional strata with a higher-level **context** – hospitals
(patient survival), schools (student achievement), neighbourhoods. Pass
`context = "school"` to
[`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
or [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md) to
fit `outcome ~ covars + (1 | stratum) + (1 | school)` in one model; the
summary then splits the unexplained variance into **between-stratum**
(intersectional, net of the context), **between-context** (the general
contextual effect), and **residual**, and `plot(type = "context_vpc")`
visualises the partition.

``` r

data(maihda_country_data)
# Strata (gender x SES) crossed with country in ONE model -- contrast with
# group = "country", which instead fits an independent model per country.
a <- maihda(math ~ 1 + (1 | gender:ses), data = maihda_country_data,
            context = "country")
a$summary$context$vpc_context_total  # the country share of unexplained variance
plot(a, type = "context_vpc")
```

A context with few levels (like these 6 countries) weakly identifies its
variance; prefer many-level contexts or `engine = "brms"` for serious
use.

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
- **Mean Prediction vs. Stratum Effect**: Plots each stratum’s mean
  predicted outcome against its stratum random effect (the direction of
  “worse”/“better” depends on the outcome, so it is not framed as risk)
- **Effect Decomp**: Additive versus specific interaction decompositions

### `maihda_ternary_plot()`

Generates a ternary diagnostic plot. For each stratum it normalizes
three magnitudes to sum to 1: the additive signal (how far the
fixed-effect-only prediction sits from the grand mean), the
intersection-specific signal (the magnitude of the stratum random
effect), and the uncertainty in that estimate. It is a relative-signal
diagnostic, not a formal variance decomposition.

### `plot_prediction_deviation_panels()`

Creates an advanced, publication-ready two-panel dashboard for
visualizing predicted values and highlighting the most notable cases or
strata. What counts as notable depends on the model type — the largest
deviation from the mean prediction (Gaussian/Poisson), the largest
absolute deviance residual (binomial), or the most surprising
observation (ordinal) — and the labelled points are not
regression-diagnostic outliers.

### `compare_maihda()`

Compares VPC/ICC across multiple models with optional bootstrap
confidence intervals.

### `compare_maihda_groups()`

Compares intersectional inequality (VPC/ICC and between-/within-stratum
variance) across the levels of a higher-level grouping variable such as
country, region, or survey wave, fitting a stratified MAIHDA model per
group. Visualize with `plot(result, type = "vpc")`. The bundled
`maihda_country_data` (OECD PISA 2018; gender × socioeconomic-status
strata across six countries) is built to demonstrate this.

### `calculate_pvc()`

Calculates the proportional change in between-stratum variance (PCV)
between two models. It is the share of the baseline model’s
between-stratum variance explained by the second model only when the
second nests the first (adding predictors on the same outcome, sample
and strata); otherwise it is a model-dependent change in variance.

- Formula: PCV = (Var_model1 - Var_model2) / Var_model1
- Works with both lme4 and brms engines
- Supports bootstrap confidence intervals for lme4 models

### `stepwise_pcv()`

Evaluates multiple sequential models by iteratively adding covariates
step-by-step. Each step’s PCV is the change in between-stratum variance
contributed by a predictor given the variables already in the model, so
it is order-dependent rather than an order-invariant “unique”
contribution.

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

# Get the variance partition coefficient (VPC/ICC)
summary <- summary(model, bootstrap = TRUE, n_boot = 1000)

# The VPC is the share of the *unexplained* (between-stratum + residual) variance
# that lies between strata -- here conditional on Age, and on the latent scale for
# non-Gaussian models (e.g. logistic uses pi^2/3 as the level-1 variance). Read
# from the null model `outcome ~ 1 + (1 | stratum)` it is the total between-stratum
# share. NOTE: the between-stratum variation reflects the combined additive +
# interaction differences across strata; it represents the *pure* intersectional
# (interaction) effect only once the additive main effects of the strata variables
# are added to the model.

# To get that pure intersectional part, fit the ADJUSTED model -- the same model PLUS
# the strata variables as additive fixed effects -- and read the PCV (the additive
# share of the intersectional inequality). The one-call maihda() does this for you;
# here it is spelled out with the building blocks.
adjusted <- fit_maihda(
  BMI ~ Age + Gender + Race + Education + (1 | Gender:Race:Education),
  data = maihda_health_data
)
calculate_pvc(model, adjusted)  # null -> adjusted: proportional change in variance

# Visualize which strata have higher/lower outcomes using new advanced plots
plot(model, type = "predicted")
plot(model, type = "risk_vs_effect")

# Map out where the intersection-specific signal is concentrated
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

## Design-Weighted MAIHDA (Survey Data)

For complex-survey data (NHANES, PISA, …), pass the sampling-weight
column via `sampling_weights`. Survey weights are **not** lme4
`weights=` (those are precision weights), so the fit routes through
[`WeMix::mix()`](https://american-institutes-for-research.github.io/WeMix/reference/mix.html)
– weighted pseudo-maximum-likelihood (Rabe-Hesketh & Skrondal 2006) –
and `engine = "lme4"` with sampling weights is an error rather than a
silent misfit.

``` r

# One call: design-weighted null + adjusted models and PCV. The engine switches
# to "wemix" automatically (install WeMix from CRAN).
analysis <- maihda(
  outcome ~ age + (1 | gender:race:education),
  data = survey_data,
  sampling_weights = "person_weight"
)
analysis            # design-weighted VPC/ICC and PCV
summary(analysis)   # design-consistent (sandwich) SEs for the fixed effects

# Works across the toolkit: stepwise PCV, group comparison, prediction, plots,
# and the design-weighted AUC for binary outcomes.
stepwise_pcv(strata_data, "outcome", c("gender", "race", "education"),
             sampling_weights = "person_weight")

# Bayesian alternative: weights enter as likelihood weights (pseudo-posterior --
# point estimates are design-consistent; credible intervals are not design-based).
fit_maihda(outcome ~ age + (1 | gender:race:education), data = survey_data,
           engine = "brms", sampling_weights = "person_weight")
```

Limitations (explicit, not silent): the wemix engine covers the
canonical `gaussian(identity)` / `binomial(logit)` MAIHDA with a single
`(1 | stratum)` intercept; crossed random effects (`context =`,
`decomposition = "crossed-dimensions"`) and bootstrap intervals require
lme4/brms. A design-based interval would need replicate weights, which
is a possible future extension.

## Model Comparison with Bootstrap

``` r

# Compare nested models on the SAME data and strata (null vs covariate-adjusted).
# VPCs are only comparable when models share an outcome, family/link, analytic
# sample (the same rows), and strata; compare_maihda() warns otherwise.
strata <- make_strata(maihda_sim_data, vars = c("gender", "race"))
null_model <- fit_maihda(health_outcome ~ 1 + (1 | stratum), data = strata$data)
adj_model  <- fit_maihda(health_outcome ~ age + (1 | stratum), data = strata$data)

# Compare with bootstrap CI
comparison <- compare_maihda(
  null_model, adj_model,
  model_names = c("Null", "Adjusted"),
  bootstrap = TRUE,
  n_boot = 1000
)

# Visualize comparison
plot(comparison)
```

## Calculating Proportional Change in Variance (PCV)

``` r

# Fit baseline model
model1 <- fit_maihda(outcome ~ age + (1 | gender:race), data = data)

# Fit model with additional predictor
model2 <- fit_maihda(outcome ~ age + gender + (1 | gender:race), data = data)

# Calculate PCV without bootstrap
pvc_result <- calculate_pvc(model1, model2)
print(pvc_result)

# Calculate PCV with bootstrap confidence intervals
pvc_boot <- calculate_pvc(model1, model2, bootstrap = TRUE, n_boot = 1000)
print(pvc_boot)

# Interpretation: A PCV of 0.25 means the between-stratum variance is 25% lower in
# model2 than in model1. This is "variance explained" by the added predictors only
# when model2 nests model1 (same outcome, sample and strata); otherwise read it as
# a model-dependent change in between-stratum variance.
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
- reformulas, patchwork, ggrepel, tidyselect, stats, tibble, rlang

**Optional:**

- brms (\>= 2.15.0) - for Bayesian models
- WeMix (\>= 4.0.0) - for design-weighted (survey) models via
  `sampling_weights`
- ggtern - for ternary diagrams
- shiny, bslib, DT, plotly, shinyjs, shinycssloaders, future, promises -
  for the interactive dashboard
  ([`run_maihda_app()`](https://hdbt.github.io/MAIHDA/reference/run_maihda_app.md))
- haven - for uploading SPSS (.sav) / Stata (.dta) files in the
  dashboard

Bootstrap confidence intervals use a parametric bootstrap via
`lme4::simulate()` /
[`lme4::refit()`](https://rdrr.io/pkg/lme4/man/refit.html); no external
bootstrap package is required.

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
Bulut (2025). *MAIHDA: Multilevel Analysis of Individual Heterogeneity and Discriminatory Accuracy.* R package version 0.1.11, https://github.com/hdbt/MAIHDA. doi: 10.32614/CRAN.package.MAIHDA
```

A BibTeX entry for LaTeX users is:

``` bibtex
@Manual{Bulut2025MAIHDA,
  title  = {MAIHDA: Multilevel Analysis of Individual Heterogeneity and Discriminatory Accuracy},
  author = {Hamid Bulut},
  year   = {2025},
  note   = {R package version 0.1.11},
  url    = {https://github.com/hdbt/MAIHDA},
  doi    = {10.32614/CRAN.package.MAIHDA}
}
```
