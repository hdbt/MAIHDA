# Fit a cumulative MAIHDA model via ordinal::clmm

Internal engine call for `fit_maihda(engine = "ordinal")`. Builds the
analytic sample (complete cases on the model variables) so the stored
`data` matches the rows clmm fits, then calls
[`ordinal::clmm()`](https://rdrr.io/pkg/ordinal/man/clmm.html) with
`Hess = TRUE` (needed for the threshold standard errors). The analytic
frame is passed by NAME (bound in a private environment) so the call
clmm stores stays one line: ordinal's `print`/`summary` methods deparse
`call$data`, and embedding the frame there made printing an ordinal fit
dump the whole data set. Any formula
[`offset()`](https://rdrr.io/r/stats/offset.html) term is moved ahead of
the `(1 | stratum)` term first, because clmm reads the offset by
position in the bar-free variables list (see
`maihda_offset_before_bars`).

## Usage

``` r
maihda_fit_clmm(formula, data, family, dot_vals)
```

## Arguments

- formula:

  The resolved model formula (with `(1 | stratum)`).

- data:

  The data (after strata creation / response preparation).

- family:

  The cumulative family marker (link "logit" or "probit").

- dot_vals:

  Named list of evaluated `...` arguments forwarded to
  [`ordinal::clmm()`](https://rdrr.io/pkg/ordinal/man/clmm.html) (e.g.
  `nAGQ`, `control`).

## Value

A list with `model` (the `clmm` fit), `data` (the analytic data frame
actually fitted) and `coding` (the fit's factor levels and contrast
matrices, taken while the options that shaped the design are in force;
`NULL` if they could not be read).

## Details

A Laplace fit that ends with every stratum effect at exactly zero, or
with a between-stratum SD below 1% of the latent residual SD, is
compared with the boundary fit: `clmm()`'s general routine
(`clmm.control(useMatrix = TRUE)`) started from the first fit's
thresholds and coefficients and an SD of zero. The boundary fit is the
one returned unless its log-likelihood is below the first fit's by more
than 1e-9 of its size.

The model is fitted on its fixed-effects design with every column
centred and scaled, the boundary fit included, and returned as the
`clmm()` fit of the formula at that maximum. A covariate far from zero
or on a large or tiny scale can leave `clmm()`'s optimiser at its
starting values or short of the maximum, and `clmm()` takes the gradient
and Hessian by differences with a fixed step, too coarse for the
coefficient of such a covariate. The returned fit carries the gradient
and Hessian of the rescaled fit, mapped exactly, so its standard errors
are those of that fit. Under `threshold = "symmetric2"` the thresholds
cannot take up a shift, and the columns are scaled only. A fit given a
`start` is made as the formula stands, since the start is in the
formula's own terms.
