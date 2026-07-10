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

Aggregated-binomial fits are supported on both engines that fit them –
an lme4 `cbind(success, failure)` response and a brms `y | trials(n)`
response: the AUC is the count-weighted C-statistic over the implied
individual-level 0/1 data, and `n_case` / `n_control` are the total
successes / failures.

**Scope.** When the model carries random effects *beyond* the
intersectional partition – a contextual `(1 | school)` from
`fit_maihda(context = )` or an explicit extra grouping such as
`(1 | site)` – the headline `auc` is the concordance of the
*intersectional strata* (fixed effects plus the stratum random effect
and, for a crossed-dimensions fit, the additive dimension effects),
matching the scope of the MOR; the concordance of the full model
including the other random effects is reported separately as `auc_full`
(`auc_scope = "strata"`). For the canonical single-`(1 | stratum)` model
the two coincide and only `auc` is reported (`auc_scope = "model"`,
`auc_full` absent).

## Usage

``` r
maihda_discriminatory_accuracy(model)
```

## Arguments

- model:

  A `maihda_model` from
  [`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
  fitted with a `binomial` family – including an aggregated response (an
  lme4 `cbind(success, failure)` or a brms `y | trials(n)`) – or the
  `bernoulli` family a binary 0/1 outcome is fit with under
  `engine = "brms"`.

## Value

An object of class `maihda_da`: a list with `auc`, `auc_scope`,
`auc_full`, `mor`, `n_case`, `n_control`, `family`, `link`, `engine`,
`weighted` and `weight_type`. `mor` is `NA` for a non-logit binomial
link, where the AUC is still reported. For an aggregated-binomial fit
`n_case` / `n_control` are the total successes / failures. `weighted` is
`TRUE` when the AUC is a weighted concordance – `weight_type`
`"sampling"` for a design-weighted fit (each observation contributes its
sampling weight as case/control mass, estimating the population
discriminatory accuracy) or `"precision"` for an lme4 fit with
non-integer precision weights; `n_case` / `n_control` stay unweighted
observation counts either way. `weight_type` is `NULL` for an unweighted
AUC.

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
# }
```
