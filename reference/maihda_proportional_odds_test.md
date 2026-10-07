# Parametric-bootstrap proportional-odds test for a cumulative MAIHDA fit

Tests the proportional-odds (parallel-lines) assumption of a fitted
cumulative `clmm` MAIHDA model by calibrating the omnibus
nominal-effects likelihood-ratio statistic against its own null
distribution under the fitted model.

## Usage

``` r
maihda_proportional_odds_test(object, n_sim = 199, seed = NULL)
```

## Arguments

- object:

  A `maihda_model` fitted with `engine = "ordinal"`.

- n_sim:

  Number of parametric-bootstrap replicates (default 199).

- seed:

  Optional integer seed, for a reproducible bootstrap.

## Value

An object of class `maihda_po_test`: a list with `lrt`, `df`, `n_terms`,
`p_value` (the bootstrap p-value), `p_chisq`, `n_sim` (replicates that
produced a usable statistic), `n_failed`, and `null_lrt` (the simulated
null statistics). `p_value` is the one to report and the only one
printed. `p_chisq` is the chi-squared p-value of the same statistic; it
is not calibrated for a mixed model (see Details) and is returned for
comparison only.

## Details

The statistic is the omnibus nominal-effects likelihood-ratio statistic:
the fixed-effect part of the model is refitted twice with
[`ordinal::clm()`](https://rdrr.io/pkg/ordinal/man/clm.html) – once with
all covariates proportional, once with every covariate entering as a
threshold-specific (nominal) effect – and twice the log-likelihood
difference is taken.

Its usual chi-squared reference
([`ordinal::nominal_test()`](https://rdrr.io/pkg/ordinal/man/nominal.test.html))
is not valid here. After the stratum random intercept is integrated out,
the marginal cumulative-logit slopes generally differ across thresholds
once the stratum variance is non-zero. The chi-squared reference also
treats the observations as independent, whereas observations that share
a stratum share its random effect. A correctly specified model is
therefore rejected too often, and more often as the sample grows.

The null distribution is instead simulated under the fitted `clmm`: each
replicate redraws the stratum random effects from \\N(0, \tau^2)\\ at
the fitted variance, forms the conditional category probabilities from
the fitted thresholds and location predictor, redraws the response, and
recomputes the statistic. The p-value is \\(1 + \\\\T_b \ge T\_{obs}\\)
/ (1 + B)\\ over the \\B\\ replicates that refitted. `n_sim` sets only
its Monte Carlo resolution: the smallest attainable value is \\1 / (B +
1)\\.

This is a parametric-bootstrap approximation, not an exact test. The
replicates are drawn at the fitted stratum variance, which is estimated
from only as many units as there are strata, so the test is least
reliable when the strata are few, whatever `n_sim` is.

The test is expensive – every replicate refits two `clm()` models – and
is not run at fit time.

## See also

[`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md),
[`maihda_cumulative`](https://hdbt.github.io/MAIHDA/reference/maihda_cumulative.md)

## Examples

``` r
# \donttest{
strata <- make_strata(maihda_sim_data, vars = c("gender", "race"))
d <- strata$data
d$y <- factor(cut(d$health_outcome, 3), labels = 1:3, ordered = TRUE)
m <- fit_maihda(y ~ age + (1 | stratum), data = d, family = "ordinal")
#> fit_maihda(): ordinal (cumulative) family; using engine = "ordinal" (ordinal::clmm). Set 'engine' explicitly to silence this message or to choose engine = "brms".
maihda_proportional_odds_test(m, n_sim = 99, seed = 1)
#> Proportional-odds test (parametric bootstrap under the fitted clmm)
#> 
#>   Nominal-effects LRT : 0.011 on 1 df over 1 covariate(s)
#>   Bootstrap p-value   : 0.9500  (99 replicates)
# }
```
