# MAIHDA: Multilevel Analysis of Individual Heterogeneity and Discriminatory Accuracy

[![CRAN
status](https://www.r-pkg.org/badges/version/MAIHDA)](https://CRAN.R-project.org/package=MAIHDA)
[![CRAN
downloads](https://cranlogs.r-pkg.org/badges/MAIHDA)](https://CRAN.R-project.org/package=MAIHDA)
[![R-CMD-check](https://github.com/hdbt/MAIHDA/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/hdbt/MAIHDA/actions/workflows/R-CMD-check.yaml)
[![R](https://img.shields.io/badge/R-%3E%3D4.1.0-blue)](https://www.r-project.org/)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Codecov test
coverage](https://codecov.io/gh/hdbt/MAIHDA/branch/main/graph/badge.svg)](https://app.codecov.io/gh/hdbt/MAIHDA?branch=main)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-blue.svg)](https://lifecycle.r-lib.org/articles/figures/lifecycle-stable.svg)

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
- **[Interactive Dashboard](https://hdbt.shinyapps.io/shiny/)**: A
  fully-featured Shiny application
  ([`run_maihda_app()`](https://hdbt.github.io/MAIHDA/reference/run_maihda_app.md))
  for no-code exploratory data analysis and model fitting
- **Model Fitting**: Multiple engines — `lme4` (frequentist), `brms`
  (Bayesian), `WeMix` (design-weighted pseudo-likelihood), and `ordinal`
  ([`ordinal::clmm`](https://rdrr.io/pkg/ordinal/man/clmm.html), for
  cumulative/ordered-factor outcomes)
- **Design-Weighted MAIHDA**: Survey/sampling weights via
  `sampling_weights` fit a weighted pseudo-likelihood (WeMix for lme4, a
  pseudo-posterior for brms) targeting a population-weighted estimand.
  The WeMix path gives design-consistent (sandwich) SEs for the **fixed
  effects**; both paths give population-weighted point estimates for the
  VPC, PCV, stratum summaries and AUC. The intersectional strata are
  treated as exhaustively sampled population cells (level-2 weights =
  1), and design-based (replicate-weight / linearised) variances for the
  variance components are not computed — so this is weighting to the
  population, not full complex-survey support for arbitrary
  clustered/stratified/replicate-weight designs
- **Summaries & Decompositions**: Variance partition coefficients
  (VPC/ICC), stratum-specific estimates, and stepwise Proportional
  Change in Variance (PCV)
- **Multiple Prediction Types**: Individual-level and stratum-level
  predictions
- **Visualizations**: Predicted stratum values, VPC visualizations,
  mean-prediction vs. stratum-effect diagnostics, and observed
  vs. shrunken estimates
- **Group Comparison**:
  [`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
  contrasts intersectional inequality (VPC/ICC) across levels of a
  higher-level variable such as country or region
- **Contextual Cross-Classified MAIHDA**: `context = "school"` crosses
  the intersectional strata with a place/institution level (school,
  hospital, region)

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

## Quick Start: 30 Seconds to MAIHDA

``` r

library(MAIHDA)
data("maihda_health_data")

# Everything in one call: null + adjusted fit, VPC/ICC summary, and PCV decomposition.
analysis <- maihda(BMI ~ Age + Gender + Race + Education + (1 | Gender:Race:Education),
                   data = maihda_health_data)
analysis
summary(analysis)

plot(analysis, type = "vpc")
plot(analysis, type = "effect_decomp")
plot(analysis) # plots all plot types
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

### `maihda_table()`

Assembles the two standard MAIHDA write-up deliverables from a fitted
[`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md) analysis
(or a single
[`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
model) in one call: (a) a **model-results table** contrasting the null
(Model 1) and adjusted (Model 2) fits — intercept, between-stratum
variance and SD, VPC/ICC, the PCV, and (for a binary outcome) the AUC
and Median Odds Ratio — and (b) a **ranked-strata table** ordering every
stratum by its predicted outcome, with conditional intervals and the
stratum random effect. The `$models` data frame is numeric and
export-ready ([`write.csv()`](https://rdrr.io/r/utils/write.table.html)
/ [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html));
[`print()`](https://rdrr.io/r/base/print.html) renders the familiar
“estimate \[low, high\]” layout plus the top/bottom strata. It adapts to
every fit type (crossed-dimensions shares, contextual
`Context share (VPC)`, ordinal thresholds) and engine
(lme4/brms/WeMix/ordinal).

``` r

analysis <- maihda(BMI ~ Age + Gender + Race + (1 | Gender:Race), data = maihda_health_data)
tab <- maihda_table(analysis)
tab          # printed: Model 1 vs Model 2 table + highest/lowest strata
tab$models   # numeric, export-ready results table
tab$strata   # all strata ranked by predicted BMI
```

### `make_strata()`

Creates intersectional strata from multiple categorical variables with
optional minimum count filtering.

### `fit_maihda()`

Fits multilevel models using the lme4 (default), brms, WeMix
(design-weighted, via `sampling_weights`), or ordinal
([`ordinal::clmm`](https://rdrr.io/pkg/ordinal/man/clmm.html), selected
automatically for ordered-factor outcomes) engine. Supports various
families including gaussian, binomial, poisson, negbinomial
(overdispersed counts; theta estimated via
[`lme4::glmer.nb()`](https://rdrr.io/pkg/lme4/man/glmer.nb.html) or
brms’s `shape` parameter, with the VPC using the latent-scale level-1
variance of Nakagawa, Johnson & Schielzeth 2017), and ordinal/cumulative
(proportional-odds models for ordered outcomes; latent-scale VPC with
pi^2/3 logit / 1 probit level-1 variance, response-scale predictions as
expected category scores).

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
confidence intervals, and (by default, `ic = TRUE`) appends relative-fit
information criteria — AIC/BIC for the likelihood engines, WAIC/LOOIC
for brms — for comparing model structures.

### `maihda_ic()`

Reports the relative-fit information criteria for one or more models (or
a [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md)
analysis, expanded into its null and adjusted models) to help choose
between model *structures* — a question the VPC/PCV do not address.
AIC/BIC for the likelihood engines (lme4, ordinal) and the Bayesian
WAIC/LOOIC for brms, with a `delta` column from the best model. REML
`lmer` fits are refitted with ML so AIC/BIC are comparable across models
with different fixed effects (the null-vs-adjusted case).

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
- Point estimate works across all fitting engines (`lme4`, `brms`,
  `WeMix`, `ordinal`)
- Supports bootstrap confidence intervals for lme4 models

### `stepwise_pcv()`

Evaluates multiple sequential models by iteratively adding covariates
step-by-step. Each step’s PCV is the change in between-stratum variance
contributed by a predictor given the variables already in the model, so
it is order-dependent rather than an order-invariant “unique”
contribution. For a binary outcome it also reports the
discriminatory-accuracy trajectory (`AUC`, the step/total change in AUC,
and `MOR`) alongside the PCV, so the strata’s discriminatory accuracy
can be tracked as covariates are added.

### Longitudinal (growth-curve) MAIHDA

Supply `id` and `time` to
[`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)/[`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md)
for a 3-level growth-curve MAIHDA (occasions within individuals within
strata, with random intercept *and* slope on time at both levels), the
life-course extension of Bell, Evans, Holman & Leckie (2024). The
between-stratum VPC is then time-varying, and
`maihda(decomposition = "longitudinal")` reports the PCV separately for
the baseline (intercept) and the slope variance — the
additive-vs-multiplicative split of the intersectional *trajectory*
inequality.

``` r

data(maihda_long_data)   # 600 persons x 5 waves; gender x ethnicity x education strata

# Time-varying VPC from a 3-level growth model
m <- fit_maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
                data = maihda_long_data, id = "id", time = "wave")
summary(m)                         # baseline VPC + the VPC trajectory over waves
plot(m, type = "vpc_trajectory")   # VPC(t) curve
plot(m, type = "trajectories")     # predicted per-stratum mean trajectories

# Additive-vs-multiplicative PCV (null vs adjusted growth model)
a <- maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
            data = maihda_long_data, id = "id", time = "wave",
            decomposition = "longitudinal")
a$pcv                              # PCV_intercept (baseline) and PCV_slope (trajectory)
plot(a, type = "pcv_trajectory")
```

### `run_maihda_app()`

Launches a locally-hosted, interactive Shiny Dashboard that exposes the
core functionalities for data modeling, visualization, and summarization
visually.

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
# point estimates target the population-weighted estimand, not full design-based
# inference; credible intervals are not design-based).
fit_maihda(outcome ~ age + (1 | gender:race:education), data = survey_data,
           engine = "brms", sampling_weights = "person_weight")
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

## Tidy output (`broom`)

[`tidy()`](https://generics.r-lib.org/reference/tidy.html) and
[`glance()`](https://generics.r-lib.org/reference/glance.html) methods
turn a fitted model or analysis into tidy data for tables (`gt`,
`flextable`) and `ggplot2`, reading the quantities
[`summary()`](https://rdrr.io/r/base/summary.html) already computes (no
new statistics, no refit):

``` r

analysis <- maihda(BMI ~ Age + Gender + Race + Education + (1 | Gender:Race:Education),
                   data = maihda_health_data)

# One-row headline: VPC/ICC, PCV, AUC/MOR (binary), additive/interaction shares, n.
# This row has no equivalent in broom.mixed/easystats -- PCV needs the null+adjusted pair.
glance(analysis)

# Per-estimate tidy frames (broom's term/estimate/std.error/conf.low/conf.high shape):
tidy(analysis)                          # stratum (intersection) effects, with labels
tidy(analysis, component = "variance")  # variance components
tidy(analysis, component = "fixed", which = "adjusted")

# A caterpillar plot in two lines:
library(ggplot2)
tidy(analysis) |>
  ggplot(aes(reorder(label, estimate), estimate, ymin = conf.low, ymax = conf.high)) +
  geom_pointrange() + coord_flip()
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
Bulut (2025). *MAIHDA: Multilevel Analysis of Individual Heterogeneity and Discriminatory Accuracy.* R package version 0.1.12, https://github.com/hdbt/MAIHDA. doi: 10.32614/CRAN.package.MAIHDA
```

A BibTeX entry for LaTeX users is:

``` bibtex
@Manual{Bulut2025MAIHDA,
  title  = {MAIHDA: Multilevel Analysis of Individual Heterogeneity and Discriminatory Accuracy},
  author = {Hamid Bulut},
  year   = {2025},
  note   = {R package version 0.1.12},
  url    = {https://github.com/hdbt/MAIHDA},
  doi    = {10.32614/CRAN.package.MAIHDA}
}
```
