# Introduction to MAIHDA

## Introduction

The **MAIHDA** package provides specialized tools for conducting
Multilevel Analysis of Individual Heterogeneity and Discriminatory
Accuracy. This modern epidemiological approach is highly effective for
investigating intersectional health inequalities and understanding how
joint social categories (e.g., Race x Gender x Education) influence
individual outcomes.

By utilizing multilevel mixed-effects models (via `lme4` or `brms`),
MAIHDA allows researchers to: 1. Automatically construct intersectional
strata. 2. Estimate between-stratum variance and Variance Partition
Coefficients (VPC). 3. Evaluate the Proportional Change in Variance
(PCV) to understand how much inequalities are driven by additive main
effects versus unique intersectional effects. 4. Launch an interactive
Shiny Dashboard for code-free analysis.

## Installation

``` r

# Released version (CRAN):
install.packages("MAIHDA")

# Development version (GitHub):
# install.packages("remotes")
# remotes::install_github("hdbt/MAIHDA")
```

## Quick Start: the `maihda()` one-call workflow

If you just want the standard analysis,
[`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md) runs the
whole pipeline in a single call – it fits the model, summarises the
VPC/ICC and variance components, and (optionally) compares
intersectional inequality across a higher-level grouping variable. It
returns one object with [`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html), and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods.

``` r

library(MAIHDA)
data("maihda_health_data")

# Fit + variance summary in one call
analysis <- maihda(BMI ~ Age + (1 | Gender:Race:Education), data = maihda_health_data)
analysis                      # concise report (VPC/ICC, #strata)
summary(analysis)             # full variance components
plot(analysis, type = "vpc")  # any plot.maihda_model type works here

# Add a higher-level grouping variable to also compare across its levels.
# maihda_country_data (PISA 2018) has a real country grouping:
data("maihda_country_data")
by_country <- maihda(math ~ 1 + (1 | gender:ses), data = maihda_country_data,
                     group = "country")
by_country
plot(by_country, type = "group_vpc")
```

For a fuller walkthrough of the group comparison workflow, see the
dedicated [group comparison
vignette](https://hdbt.github.io/MAIHDA/articles/group_comparison.md).

The sections below walk through the same steps individually with
[`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
and friends, which you can use directly when you need finer control.

## Real-World Example Analysis

The package includes a pedagogical subset of the National Health and
Nutrition Examination Survey (`maihda_health_data`). We will use this to
examine how Body Mass Index (BMI) varies across intersectional
demographic groups.

> **Note.** The bundled `maihda_health_data` and `maihda_country_data`
> are for teaching only. They are non-random subsamples, and the package
> ignores the surveys’ weights and complex sampling design (and uses a
> single PISA plausible value), so the results below are **not**
> survey-representative and should not be read as substantive population
> findings.

### Step 1 & 2: Create Intersectional Strata and Fit a Null MAIHDA Model

Use
[`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
to combine multiple social categories directly in the random effect
formula. Providing the variables in the random effect combined with `:`
allows the function to automatically build the intersectional strata on
the fly and fit the multilevel model.

``` r

library(MAIHDA)

# Load the built-in NHANES dataset
data("maihda_health_data")

# PVC compares variance across models, so both models must use the same
# analytic sample. Keep complete cases for all variables used below.
health_complete <- maihda_health_data[complete.cases(
  maihda_health_data[, c("BMI", "Age", "Gender", "Race", "Education", "Poverty")]
), ]

# Fit the initial Null model with auto-generated strata
model_null <- fit_maihda(
  BMI ~ 1 + (1 | Gender:Race:Education),
  data = health_complete,
  engine = "lme4"
)

# Summarize the variance components (VPC)
summary_null <- summary(model_null)
print(summary_null)
```

**Interpretation:** The resulting Variance Partition Coefficient (VPC or
ICC) tells us what percentage of the total variance in BMI in the
population lies *between* the intersectional social groups, rather than
just *within* them.

### Step 3: Evaluate Proportional Change in Variance (PCV)

To understand if these intersectional inequalities are simply the sum of
their parts (additive), we evaluate how much variance is explained by
adding main-effects to the model.

If the variance drops substantially (high PCV), much of the
between-stratum inequality is explained by the additive main effects. If
it remains or even *increases* (negative PCV), little is explained by
the main effects. Interpret this cautiously: the PCV is a
model-dependent variance change, and a low or negative value does not by
itself prove a “true” intersectional interaction – it can also reflect
suppression, rescaling (on the latent scale for non-Gaussian models),
sample composition, or estimation uncertainty.

``` r

# Fit an adjusted model
model_adj <- fit_maihda(
  BMI ~ Age + Gender + Race + Education + Poverty + (1 | Gender:Race:Education),
  data = health_complete
)

# Calculate PCV with Parametric Bootstrap Confidence Intervals
pcv_result <- calculate_pvc(model_null, model_adj, bootstrap = TRUE, n_boot = 500)
print(pcv_result)
```

### Step 4: Stepwise PCV Decomposition

Often, researchers want to know exactly *which* variable explained the
variance. Use the
[`stepwise_pcv()`](https://hdbt.github.io/MAIHDA/reference/stepwise_pcv.md)
function to add covariates one-by-one and track the variance
dynamically.

``` r

# Run a stepwise variance decomposition using the prepared data with strata
stepwise_results <- stepwise_pcv(
  data = model_null$original_data,
  outcome = "BMI",
  vars = c("Age", "Gender", "Race", "Education", "Poverty")
)

print(stepwise_results)
```

Negative step PCVs in this table highlight “unmasking” or suppression
effects: adding a variable caused the intersectional groups to push
further apart mathematically, revealing hidden structural inequalities.

### Step 5: Visualizations

The package provides multiple pre-configured, advanced visualization
options for checking your model estimates natively mirroring the Shiny
application logic:

``` r

# Predicted stratum values with 95% CIs
plot(model_adj, type = "predicted")

# Variance partition (VPC) visualization
plot(model_adj, type = "vpc")

# Mean predicted outcome against the stratum random effect
plot(model_adj, type = "risk_vs_effect")

# Additive versus Intersectional Effect decomposition
plot(model_adj, type = "effect_decomp")

# Ternary diagnostic: additive vs intersection-specific signal vs uncertainty
plot(model_adj, type = "ternary")

# Individual Prediction Deviance Dashboard
plot(model_adj, type = "prediction_deviation")
```

### Step 6: Compare Intersectional Inequality Across Groups

When the data spans several higher-level contexts – countries, regions,
survey waves – you often want to know whether the *degree* of
intersectional inequality differs across them.
[`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
fits a separate intercept-only MAIHDA model within each level of a
grouping variable and reports the VPC/ICC and variance components side
by side.

By default the intersectional strata are defined once on the full data
(`shared_strata = TRUE`), so a stratum denotes the same combination in
every group and the VPCs are directly comparable. This is a *stratified*
analysis (one model per group), not a cross-classified model.

The bundled `maihda_country_data` (OECD PISA 2018) is designed for
exactly this: it asks how much the joint **gender x
socioeconomic-status** inequality in mathematics achievement differs
across six countries.

``` r

data("maihda_country_data")

group_cmp <- compare_maihda_groups(
  math ~ 1 + (1 | gender:ses),
  data  = maihda_country_data,
  group = "country"
)
group_cmp

# VPC by country (add bootstrap = TRUE for per-group confidence intervals)
plot(group_cmp, type = "vpc")

# Between- vs within-stratum variance share by country
plot(group_cmp, type = "components")
```

Groups smaller than `min_group_n`, or with fewer than two populated
strata, are skipped with a warning rather than producing an unstable
VPC; a group with no between-stratum variation reports a VPC of 0
instead of erroring.

## Interactive Shiny App

The MAIHDA package ships with a fully-featured, interactive Shiny
Dashboard.

You can upload your own data (CSV, SPSS `.sav`, Stata `.dta`),
dynamically select variables, and compute Stepwise PCV tables and
prediction plots.

``` r

# Launch the interactive interface
run_maihda_app()
```

## References

- Evans, C. R., Williams, D. R., Onnela, J. P., & Subramanian, S. V.
  (2018). A multilevel approach to modeling health inequalities at the
  intersection of multiple social identities. *Social Science &
  Medicine*, 203, 64-73.

- Merlo, J. (2018). Multilevel analysis of individual heterogeneity and
  discriminatory accuracy (MAIHDA) within an intersectional framework.
  *Social Science & Medicine*, 203, 74-80.
