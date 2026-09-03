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
null statistics). `p_value` is the only p-value the print method shows.
`p_chisq` is the uncalibrated chi-squared p-value that the removed
automatic screen used; it is retained on the object for comparison but
deliberately not printed, and it is not evidence against the fitted
model.

## Details

The statistic is the ordinary omnibus nominal-effects LRT: the
fixed-effect part of the model is refitted twice with
[`ordinal::clm()`](https://rdrr.io/pkg/ordinal/man/clm.html) – once with
all covariates proportional, once with every covariate entering as a
threshold-specific (nominal) effect – and twice the log-likelihood
difference is taken. Because `clm()` has no random effect, that
statistic is computed on the *marginal* fit.

Referring it to a chi-squared distribution, as an ordinary
[`ordinal::nominal_test()`](https://rdrr.io/pkg/ordinal/man/nominal.test.html)
would, is not valid here. A conditional cumulative model with a normal
random intercept does not in general remain an ordinary
proportional-odds model after the random intercept is marginalised away:
the implied marginal cumulative-logit slopes differ across thresholds
for any non-zero stratum variance. The chi-squared null is therefore
false under the correctly specified model, and its rejection rate grows
without bound in the sample size – at a stratum VPC near 7 percent, a
plausible MAIHDA value, a correctly specified model is rejected about a
quarter of the time at `n = 96,000`. Stratum heterogeneity and genuine
non-proportional odds are not separable by the fixed-only statistic
alone.

The one exception is an exactly symmetric threshold configuration, where
the marginal slopes coincide and the fixed-only statistic is valid: for
symmetric \\u\\ the map \\\eta \mapsto\\ logit
\\E\[\mathrm{plogis}(\eta - u)\]\\ is odd, so its derivative is even and
thresholds placed symmetrically about the location share a slope. Three
categories cut at \\-c\\ and \\+c\\ is the case that arises in practice.
Asymmetric thresholds are markedly worse rather than better, so this
exception narrows the problem without softening it.

This function removes that confounding by simulating the null
distribution *under the fitted `clmm` itself*: each replicate redraws
the stratum random effects from \\N(0, \tau^2)\\ at the fitted variance,
forms the conditional category probabilities from the fitted thresholds
and location predictor, redraws the ordinal response, and recomputes the
same fixed-only statistic. The reported p-value is \\(1 + \\\\T_b \ge
T\_{obs}\\) / (1 + B)\\, so proportional-odds data generated from the
fitted model rejects at the nominal rate by construction.

`ordinal` supplies no
[`simulate()`](https://rdrr.io/r/stats/simulate.html) method for `clmm`,
so the simulation is built directly from the fitted thresholds, location
coefficients and random-effect variance.

The test is opt-in because it is expensive: every replicate refits two
`clm()` models, so the cost is roughly `n_sim` times the cost of the
fixed-only refit and grows with the sample size. It is not run
automatically at fit time.

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
