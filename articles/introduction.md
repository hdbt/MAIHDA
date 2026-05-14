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

You can install the development version of MAIHDA from GitHub:

``` r

install.packages("MAIHDA")
# Or for the latest development version:
# install.packages("remotes")
# remotes::install_github("hamidbulut/MAIHDA")
```

## Real-World Example Analysis

The package includes a pedagogical subset of the National Health and
Nutrition Examination Survey (`maihda_health_data`). We will use this to
examine how Body Mass Index (BMI) varies across intersectional
demographic groups.

### Step 1 & 2: Create Intersectional Strata and Fit a Null MAIHDA Model

Use
[`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
to combine multiple social categories directly in the random effect
formula. Providing the variables in the random effect combined with `+`
allows the function to automatically build the intersectional strata on
the fly and fit the multilevel model.

``` r

library(MAIHDA)

# Load the built-in NHANES dataset
data("maihda_health_data")

# Fit the initial Null model with auto-generated strata
model_null <- fit_maihda(
  BMI ~ 1 + (1 | Gender:Race:Education),
  data = maihda_health_data,
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

If the variance drops significantly (High PCV), the inequalities are
largely explained by the additive characteristics. If the variance
remains or even *increases* (Negative PCV), it signifies strong, unique
intersectional interactions that cannot be explained away by simple main
effects.

``` r

# Fit an adjusted model
model_adj <- fit_maihda(
  BMI ~ Age + Gender + Race + Education + Poverty + (1 | Gender:Race:Education),
  data = maihda_health_data
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

# Run a stepwise variance decomposition (using the auto-generated strata from model_null)
stepwise_results <- stepwise_pcv(
  data = model_null$data,
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

# Caterpillar plot of predictions for stratum random effects (with 95% CIs)
plot(model_adj, type = "predicted")

# Variance partition (VPC) visualization
plot(model_adj, type = "vpc")

# Bivariate risk against stratum-level intersectional effect
plot(model_adj, type = "risk_vs_effect")

# Additive versus Intersectional Effect decomposition
plot(model_adj, type = "effect_decomp")

# Ternary Plot of Variances
plot(model_adj, type = "ternary")

# Individual Prediction Deviance Dashboard
plot(model_adj, type = "prediction_deviation")
```

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
