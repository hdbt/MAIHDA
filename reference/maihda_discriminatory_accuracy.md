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

## Usage

``` r
maihda_discriminatory_accuracy(model)
```

## Arguments

- model:

  A `maihda_model` from
  [`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
  fitted with a `binomial` family (lme4 engine).

## Value

An object of class `maihda_da`: a list with `auc`, `mor`, `n_case`,
`n_control`, `family`, `link` and `engine`. `mor` is `NA` for a
non-logit binomial link, where the AUC is still reported.

## References

Merlo, J. (2018). Multilevel analysis of individual heterogeneity and
discriminatory accuracy (MAIHDA) within an intersectional framework.
*Social Science & Medicine*, 203, 74-80.

## See also

[`maihda_auc`](https://hdbt.github.io/MAIHDA/reference/maihda_auc.md),
[`maihda_mor`](https://hdbt.github.io/MAIHDA/reference/maihda_mor.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Obese (Yes/No) by intersectional strata of Gender x Race
strata <- make_strata(maihda_health_data, vars = c("Gender", "Race"))
d <- maihda_health_data
d$stratum <- strata$data$stratum
m <- fit_maihda(Obese ~ (1 | stratum), data = d, family = "binomial")
maihda_discriminatory_accuracy(m)
} # }
```
