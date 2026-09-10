# Posterior-mean cumulative thresholds of a brms cumulative fit

The brms analogue of `clmm`'s `object$model$alpha`: the ordered cut
points \\\alpha_k\\ of a
[`brms::cumulative()`](https://paulbuerkner.com/brms/reference/brmsfamily.html)
fit, read as posterior means from
[`brms::fixef()`](https://rdrr.io/pkg/nlme/man/fixed.effects.html),
where they appear as `Intercept[1]`, `Intercept[2]`, ... (any location
predictors are separate rows and are dropped here). Returned in
threshold order so they pair with the
`brms::posterior_linpred(re_formula = NA)` location – which excludes the
thresholds, exactly the latent \\\eta\\ that
[`maihda_ordinal_category_probs`](https://hdbt.github.io/MAIHDA/reference/maihda_ordinal_category_probs.md)
expects (\\P(Y \le k) = g^{-1}(\alpha_k - \eta)\\).

## Usage

``` r
maihda_brms_ordinal_thresholds(model)
```

## Arguments

- model:

  A fitted `brmsfit` from
  [`brms::cumulative()`](https://paulbuerkner.com/brms/reference/brmsfamily.html).

## Value

A numeric vector of thresholds (length \\K-1\\ for \\K\\ categories).

## Details

Unlike `clmm`'s `$alpha`, this needs no expansion under a structured
`threshold`. brms keeps the constrained parameterisation in the model
block only – an equidistant fit samples `first_Intercept` and `delta` –
and its generated quantities write the full `b_Intercept[1..K-1]`
vector, which is what `fixef()` reports. Verified on a fitted
equidistant model: `fixef()` carries all \\K-1\\ `Intercept[k]` rows
(the raw `delta` is not among them) and the cut points come back equally
spaced. Contrast
[`maihda_clmm_cutpoints`](https://hdbt.github.io/MAIHDA/reference/maihda_clmm_cutpoints.md),
which must do that expansion itself.
