# Information criteria for MAIHDA models

Reports the relative-fit information criteria for one or more MAIHDA
models, to help choose between model *structures* (different covariate
sets, strata definitions, or families) – a question the VPC/ICC and PCV
do not address. The criteria reported depend on the engine: **AIC** and
**BIC** for the likelihood engines (`lme4`, and
[`ordinal::clmm`](https://rdrr.io/pkg/ordinal/man/clmm.html)), and the
Bayesian **WAIC** and **LOOIC** (leave-one-out information criterion)
for `brms`. Lower is better for all four.

## Usage

``` r
maihda_ic(..., model_names = NULL)
```

## Arguments

- ...:

  One or more `maihda_model` objects (from
  [`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md))
  or `maihda_analysis` objects (from
  [`maihda`](https://hdbt.github.io/MAIHDA/reference/maihda.md)). A
  `maihda_analysis` contributes its null model and, when present, its
  adjusted model as separate rows.

- model_names:

  Optional character vector of names, one per `...` argument. A
  `maihda_analysis` argument's null/adjusted rows are suffixed from its
  name.

## Value

A `data.frame` of class `maihda_ic` with one row per model and the
columns that apply: `model`, `n` (analytic sample size), `estimator`,
`df` (number of parameters; likelihood engines), `logLik`, `AIC`, `BIC`
(likelihood engines), `WAIC`, `LOOIC` (brms), and – when more than one
model is supplied – `delta` (the difference from the best model on the
primary criterion: AIC for the likelihood engines, LOOIC for brms).
Columns that are entirely `NA` across the supplied models are dropped.

## Details

**REML vs ML.** `lmer` fits Gaussian models by REML by default, and a
REML log-likelihood (hence its AIC/BIC) is *not* comparable across
models with different fixed effects – exactly the canonical MAIHDA
null-vs-adjusted comparison. When more than one model is supplied,
`maihda_ic()` therefore refits any REML `lmer` model with maximum
likelihood ([`refitML`](https://rdrr.io/pkg/lme4/man/refitML.html))
before computing AIC/BIC, matching the behaviour of
[`anova()`](https://rdrr.io/r/stats/anova.html) on `lme4` models; the
`estimator` column records when this happened. For a single model the
criterion is reported as fitted (the `estimator` column then reads
`"REML"`).

**Comparability.** Like the VPC, information criteria are only
comparable across models fitted to the *same* analytic sample (same rows
and outcome) with the *same* weights – prior (precision) weights and
sampling (design) weights each change which likelihood, or
pseudo-likelihood, is being maximised, so the criteria of a weighted and
an unweighted fit of the identical model are not on a common scale.
AIC/BIC additionally require the same response distribution – they are
not comparable across families (e.g. a Gaussian vs a Poisson fit), nor
between the likelihood engines and `brms` (AIC/BIC vs WAIC/LOOIC are
different scales). When the supplied models differ in any of these
respects `maihda_ic()` warns and omits the `delta` column, still
reporting each model's own criteria;
[`compare_maihda`](https://hdbt.github.io/MAIHDA/reference/compare_maihda.md)
warns on the same grounds.

**Predictive target of the Bayesian criteria.**
[`brms::waic()`](https://mc-stan.org/loo/reference/waic.html) and
[`brms::loo()`](https://mc-stan.org/loo/reference/loo.html) are computed
from pointwise log-likelihoods *conditional on the fitted random
effects*, so WAIC/LOOIC assess prediction of new observations *within
the strata (and persons or contexts) already represented in the data* –
not performance for a new, unseen intersectional stratum. Choosing
between strata definitions on LOOIC therefore compares conditional
predictive fit; generalisation to new strata is a leave-one-group-out
cross-validation question (e.g.
[`brms::kfold()`](https://mc-stan.org/loo/reference/kfold-generic.html)
with `group = "stratum"`), which this package does not wrap. AIC/BIC for
the likelihood engines are instead computed from the *marginal*
likelihood (random effects integrated out) – a further reason the
likelihood and Bayesian criteria are never comparable with each other.

**Design-weighted fits.** For the `wemix` (design-weighted) engine the
criteria are reported as `NA`: a pseudo-likelihood with sampling weights
does not define a standard AIC/BIC. A `brms` fit with `sampling_weights`
is treated the same way: its sampling weights enter as likelihood
weights, giving a pseudo-posterior whose weighted pointwise
log-likelihoods are not log predictive densities, so WAIC/LOOIC are
likewise reported as `NA` (the `estimator` column reads
`"Bayesian (weighted pseudo-posterior)"`).

## See also

[`compare_maihda`](https://hdbt.github.io/MAIHDA/reference/compare_maihda.md),
which reports these criteria alongside the VPC/ICC, and
[`calculate_pcv`](https://hdbt.github.io/MAIHDA/reference/calculate_pcv.md)
for the variance decomposition.

## Examples

``` r
# \donttest{
strata <- make_strata(maihda_sim_data, vars = c("gender", "race"))
null_model <- fit_maihda(health_outcome ~ 1 + (1 | stratum), data = strata$data)
adj_model  <- fit_maihda(health_outcome ~ age + (1 | stratum), data = strata$data)

# AIC/BIC for two nested structures (REML lmer fits are ML-refitted first)
maihda_ic(null_model, adj_model, model_names = c("Null", "Adjusted"))
#> MAIHDA Information Criteria
#> ===========================
#> 
#>     model   n            estimator df logLik  AIC  BIC delta
#>      Null 500 ML (refit from REML)  3  -1918 3843 3855 41.68
#>  Adjusted 500 ML (refit from REML)  4  -1897 3801 3818  0.00
#> 
#> delta = difference from the best model on AIC (lower is better).
#> REML lmer fit(s) were refitted with ML so AIC/BIC are comparable across different fixed effects.
#> Information criteria are only comparable across models fitted to the same analytic sample with the same weights (and, for AIC/BIC, the same family).
#> 

# Or straight from a one-call maihda() analysis (null + adjusted rows)
a <- maihda(health_outcome ~ age + gender + race + (1 | gender:race),
            data = maihda_sim_data)
maihda_ic(a)
#> MAIHDA Information Criteria
#> ===========================
#> 
#>              model   n            estimator df logLik  AIC  BIC delta
#>      Model1 (Null) 500 ML (refit from REML)  4  -1897 3801 3818 11.36
#>  Model1 (Adjusted) 500 ML (refit from REML)  8  -1887 3790 3823  0.00
#> 
#> delta = difference from the best model on AIC (lower is better).
#> REML lmer fit(s) were refitted with ML so AIC/BIC are comparable across different fixed effects.
#> Information criteria are only comparable across models fitted to the same analytic sample with the same weights (and, for AIC/BIC, the same family).
#> 
# }
```
