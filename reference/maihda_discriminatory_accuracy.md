# Discriminatory accuracy of a binary MAIHDA model

Bundles the individual-level discriminatory-accuracy summaries for a
binomial MAIHDA model: the AUC / C-statistic (how well the model's
predicted probabilities separate cases from non-cases) and the Median
Odds Ratio. Applied to a strata-only (null) model, the AUC is the
discriminatory accuracy of the intersectional strata themselves –
Merlo's central quantity; comparing it with an adjusted model shows
whether individual covariates beyond stratum membership sharpen
classification. The AUC is computed for any binomial link; the Median
Odds Ratio is reported only for the logit link and is `NA` otherwise
(e.g. for a probit fit), since the MOR is an odds-ratio-scale quantity.

**The AUC is apparent (in-sample).** It scores the same observations
used to estimate the model – the fixed effects, the variance components,
and the shrunken stratum BLUPs – so it is an *apparent* (resubstitution)
AUC and is optimistically biased, the more so with small or sparse
strata. It is reported as the conventional descriptive MAIHDA
discriminatory accuracy (Merlo 2018), **not** as a cross-validated
estimate of out-of-sample predictive discrimination; the returned object
carries `apparent = TRUE` and
[`print()`](https://rdrr.io/r/base/print.html) says so. For genuine
predictive accuracy, validate with an out-of-fold (group-aware) scheme
or an optimism correction.

Aggregated-binomial fits are supported in all three forms R accepts – an
lme4 `cbind(success, failure)` response, an lme4 response in \[0, 1\]
whose trial counts are supplied as `weights =` (see
[`glm`](https://rdrr.io/r/stats/glm.html)), and a brms `y | trials(n)`
response: the AUC is the count-weighted C-statistic over the implied
individual-level 0/1 data, and `n_case` / `n_control` are the total
successes / failures. All three spellings of the same data give the same
AUC. This includes individual 0/1 records collapsed to frequency cells
(every row all-success or all-failure), whose weights are trial counts
like any others; see `binomial_weights` to override that reading.

**Scope.** When the model carries random effects *beyond* the
intersectional partition – a contextual `(1 | school)` from
`fit_maihda(context = )` or an explicit extra grouping such as
`(1 | site)` – the headline `auc` is the *intersectional-scope*
concordance: it **excludes** those other random effects but keeps the
fixed effects plus the stratum random effect (and, for a
crossed-dimensions fit, the additive dimension effects). The concordance
of the full model including the other random effects is reported
separately as `auc_full` (`auc_scope = "intersectional"`). **Caveat –
this is not strata-only discrimination.** The intersectional-scope score
retains the *entire* fixed-effects predictor, so when the model is
adjusted for individual-level covariates (e.g. `age`, `income` that vary
*within* strata) those covariates enter this AUC too; it is the
concordance of the adjusted fixed effects plus the intersectional random
effect(s), not of the strata alone, and it then matches the
between-stratum MOR's scope only when the fixed part is intercept-only.
For a strata-only discriminatory accuracy, score the null (strata-only)
model. For the canonical single-`(1 | stratum)` model with no other
random effects, the full and intersectional scopes coincide and only
`auc` is reported (`auc_scope = "model"`, `auc_full` absent).

## Usage

``` r
maihda_discriminatory_accuracy(
  model,
  binomial_weights = c("auto", "trials", "analytic")
)
```

## Arguments

- model:

  A `maihda_model` from
  [`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
  fitted with a `binomial` family – including an aggregated response (an
  lme4 `cbind(success, failure)` or a brms `y | trials(n)`) – or the
  `bernoulli` family a binary 0/1 outcome is fit with under
  `engine = "brms"`.

- binomial_weights:

  How to read non-unit `weights =` on an **lme4** binomial fit when
  computing the AUC. `"auto"` (default) reads integral weights as
  **trial counts**, which is what R documents a binomial prior weight to
  be ([`glm`](https://rdrr.io/r/stats/glm.html): "For a binomial GLM
  prior weights are used to give the number of trials"), and leaves
  non-integral weights – which cannot be counts – on the
  observation-level path, flagged `precision_weights_ignored`.
  `"trials"` forces the trial-count reading even for non-integral
  weights (the fractional case/control mass is carried as it stands,
  with a warning). `"analytic"` forces the ordinary observation-level
  concordance, ignoring the weights and setting
  `precision_weights_ignored = TRUE`; it is an error for a response
  carrying values strictly between 0 and 1, which has no
  observation-level case/control reading. Has no effect on a
  `cbind(success, failure)` response (its denominator is structural), on
  a brms `y | trials(n)` fit, or on a design-weighted
  (`sampling_weights`) fit.

## Value

An object of class `maihda_da`: a list with `auc`, `auc_scope`,
`auc_full`, `mor`, `n_case`, `n_control`, `family`, `link`, `engine`,
`weighted`, `weight_type`, `precision_weights_ignored` and `apparent`
(always `TRUE` – the AUC is in-sample; see Description). `mor` is `NA`
for a non-logit binomial link, where the AUC is still reported. For an
aggregated-binomial fit `n_case` / `n_control` are the total successes /
failures. `weighted` is `TRUE` only when the AUC is a genuinely weighted
(population-mass) concordance, with `weight_type` `"sampling"` for a
design-weighted fit (each observation contributes its sampling weight as
case/control mass, estimating the population discriminatory accuracy);
`n_case` / `n_control` stay unweighted observation counts. `weight_type`
is `NULL` for an unweighted AUC. An **aggregated-binomial** fit is
reported `weighted = FALSE` with `weight_type = NULL`: its
count-weighted AUC equals the ordinary individual-level concordance over
the implied 0/1 data (the trial counts are real observations, not
sampling weights), so it is not a design-weighted population quantity.
Aggregation is recognised from an lme4 `cbind(success, failure)` matrix
response, from non-unit integral lme4 `weights=` on a response in \[0,
1\] (the trial counts, see [`glm`](https://rdrr.io/r/stats/glm.html)),
or from a brms `y | trials(n)` term. The weights-based form covers a
response of exactly 0/1 – individual records collapsed to frequency
cells – because for the binomial family a prior weight *is* a trial
count: the weighted log-likelihood of one 0/1 row is exactly that of `w`
trials sharing its outcome, so the fit equals the row-expanded one and
its AUC should too. Pass `binomial_weights` to override. When a response
times its trial counts is *not* a whole number of successes – a
malformed binomial, which `glmer` itself warns about as "non-integer
\#successes" – the fractional case/control mass is carried into the AUC
as it stands, with a warning, rather than rounded into whole
observations that were never collected; `n_case` / `n_control` are then
fractional. **Weights not read as trial counts** – non-integral
`weights=`, or any weights under `binomial_weights = "analytic"` – are
ignored by the AUC, which is then the ordinary observation-level
concordance (`weighted = FALSE`), with
`precision_weights_ignored = TRUE` flagging that such weights were
present.

## References

Merlo, J. (2018). Multilevel analysis of individual heterogeneity and
discriminatory accuracy (MAIHDA) within an intersectional framework.
*Social Science & Medicine*, 203, 74-80.

## See also

[`maihda_auc`](https://hdbt.github.io/MAIHDA/reference/maihda_auc.md),
[`maihda_mor`](https://hdbt.github.io/MAIHDA/reference/maihda_mor.md)

## Examples

``` r
# \donttest{
# Obese (Yes/No) by intersectional strata of Gender x Race
strata <- make_strata(maihda_health_data, vars = c("Gender", "Race"))
d <- maihda_health_data
d$stratum <- strata$data$stratum
m <- fit_maihda(Obese ~ (1 | stratum), data = d, family = "binomial")
#> Binary outcome 'Obese' recoded to 0/1: 'No' = 0 (reference), 'Yes' = 1 (modeled event). Set the factor levels (or supply a 0/1 outcome) to control which level is the event.
maihda_discriminatory_accuracy(m)
#> Discriminatory accuracy (binomial MAIHDA)
#>   AUC (C-statistic): 0.571
#>   Median Odds Ratio: 1.482
#>   Cases / controls:  1077 / 1923
#>   (AUC is apparent / in-sample: scored on the same rows used to fit the
#>   model, so it is optimistically biased -- more so with sparse strata. It
#>   is a descriptive measure, not cross-validated out-of-sample discrimination.)
#> 
# }
```
