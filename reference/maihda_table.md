# Canonical MAIHDA results table and ranked-strata table

Assembles the two standard MAIHDA write-up deliverables from a fitted
analysis in one call: (a) a **model-results table** contrasting the null
and adjusted models (intercept, between-stratum variance and SD,
VPC/ICC, the PCV, and – for a binary outcome – the AUC and Median Odds
Ratio), and (b) a **ranked-strata table** ordering the intersectional
strata by their predicted outcome, so the best- and worst-off strata can
be read directly. It introduces no new estimator: the model-results
table reuses the quantities from
[`summary()`](https://rdrr.io/r/base/summary.html) (calling
[`summary()`](https://rdrr.io/r/base/summary.html) itself for a bare
[`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
model), and the ranked-strata table reuses the same stratum predictions
as `plot(type = "predicted")`, so the table agrees exactly with
[`summary()`](https://rdrr.io/r/base/summary.html) and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Usage

``` r
maihda_table(
  x,
  n_strata = 10L,
  scale = c("response", "link"),
  which = c("null", "adjusted"),
  digits = 3,
  ...
)
```

## Arguments

- x:

  A `maihda_analysis` from
  [`maihda`](https://hdbt.github.io/MAIHDA/reference/maihda.md) (the
  usual input), or a single `maihda_model` from
  [`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md).

- n_strata:

  Number of strata to show at each end (top and bottom) in the printed
  ranked-strata table. The returned `$strata` always holds all strata.
  Default 10.

- scale:

  Scale for the predicted stratum values: `"response"` (default) or
  `"link"`. For a cumulative (ordinal) model the response scale is the
  expected category score.

- which:

  For a two-model analysis, which model's predictions to rank the strata
  by: `"null"` (default) or `"adjusted"`. Ignored for a
  crossed-dimensions analysis or a single model.

- digits:

  Number of decimal places for the
  [`print()`](https://rdrr.io/r/base/print.html) method. Default 3.

- ...:

  For a `maihda_model` input, additional arguments passed to
  [`summary.maihda_model`](https://hdbt.github.io/MAIHDA/reference/summary.maihda_model.md)
  (e.g. `bootstrap = TRUE`); ignored for a `maihda_analysis` input,
  whose summaries are already computed.

## Value

An object of class `maihda_table`: a list with

- models:

  a data frame of the model-results table (statistics in rows; one
  estimate column per model, each with `*_lower`/`*_upper` interval
  columns)

- strata:

  a data frame of all strata ranked by predicted outcome, with `rank`,
  `stratum`, `label`, `n`, the predicted value and its conditional
  interval, and the stratum random effect and its interval; `NULL` if
  the stratum predictions could not be computed

- model_keys, model_labels:

  the estimate-column keys and their display labels

- family, engine, mode, scale, ranked_by, n_obs, n_strata_total,
  context_vars:

  metadata used by [`print()`](https://rdrr.io/r/base/print.html)

## Details

The model-results table is mostly numeric and export-ready (e.g.
`write.csv(maihda_table(a)$models, ...)` or pass it to
[`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html)):
statistics are rows, models are columns, and each estimate has
accompanying `*_lower`/`*_upper` columns that hold the
confidence/credible interval when one is available (the VPC bootstrap or
posterior interval, and the bootstrap PCV interval) and `NA` otherwise.
The intercept and the variance/SD rows are point estimates. The
[`print()`](https://rdrr.io/r/base/print.html) method renders the same
table in the familiar “estimate \[low, high\]” layout.

For a `"crossed-dimensions"` analysis (one model, no null/adjusted pair)
the results table has a single estimate column and gains “Additive
share” / “Interaction share” rows instead of the PCV. For a contextual
cross-classified analysis (`maihda(context = )`) it gains a “Context
share (VPC)” row. A bare
[`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
model is also accepted and yields a single-model table (no PCV).

The ranked-strata table ranks every stratum by its model-predicted
outcome (on the `scale` requested), using the same stratum predictions
as `plot(type = "predicted")`: the predicted value carries the
conditional (random-effect) interval, and the stratum random effect
(BLUP) is reported alongside it. By default the ranking uses the
**null** model – the headline intersectional inequality (which strata
fare best/worst overall); set `which = "adjusted"` to rank by the
adjusted model instead. The full ranked table is returned in `$strata`;
[`print()`](https://rdrr.io/r/base/print.html) shows the top and bottom
`n_strata`.

## See also

[`maihda`](https://hdbt.github.io/MAIHDA/reference/maihda.md),
[`summary.maihda_model`](https://hdbt.github.io/MAIHDA/reference/summary.maihda_model.md),
[`calculate_pvc`](https://hdbt.github.io/MAIHDA/reference/calculate_pvc.md),
[`maihda_discriminatory_accuracy`](https://hdbt.github.io/MAIHDA/reference/maihda_discriminatory_accuracy.md).

## Examples

``` r
# \donttest{
data(maihda_health_data)
a <- maihda(BMI ~ Age + Gender + Race + (1 | Gender:Race), data = maihda_health_data)

tab <- maihda_table(a)
tab                 # printed: model-results table + top/bottom strata
#> MAIHDA Results Table
#> ====================
#> 
#> Engine: lme4 | Family: gaussian | Mode: two-model
#> Observations: 3000 | Strata: 10
#> 
#> Model results:
#>                 Statistic Null (Model 1) Adjusted (Model 2)
#>                 Intercept         28.034             30.103
#>  Between-stratum variance          2.738              1.379
#>        Between-stratum SD          1.655              1.174
#>                   VPC/ICC          0.058              0.030
#>    PCV (null -> adjusted)                             0.821
#> 
#> Strata ranked by predicted value (null model):
#>  Rank           Stratum    N               Predicted              Stratum RE
#>     1    female × Black  182 32.123 [31.198, 33.047]    3.285 [2.361, 4.210]
#>     2  female × Mexican   99 29.557 [28.344, 30.771]   0.808 [-0.405, 2.021]
#>     3    male × Mexican  143 29.291 [28.259, 30.323]   0.487 [-0.545, 1.519]
#>     4      male × White  990 29.178 [28.768, 29.589]   0.272 [-0.138, 0.682]
#>     5   male × Hispanic   75 28.974 [27.610, 30.338]   0.192 [-1.172, 1.556]
#>     6      male × Black  154 28.913 [27.915, 29.911]   0.043 [-0.955, 1.041]
#>     7 female × Hispanic   91 28.623 [27.365, 29.881]  -0.127 [-1.385, 1.131]
#>     8    female × White 1044 28.424 [28.024, 28.824] -0.501 [-0.901, -0.101]
#>     9      male × Other  111 26.767 [25.612, 27.921] -2.015 [-3.170, -0.861]
#>    10    female × Other  111 26.357 [25.202, 27.511] -2.444 [-3.598, -1.289]
#>   Predicted intervals are conditional (random-effect) only; Stratum RE is the stratum BLUP.
#> 
#> Estimates are point values unless a [low, high] interval is shown (VPC/PCV).
tab$models          # the numeric, export-ready results table
#>                  statistic        null null_lower null_upper    adjusted
#> 1                Intercept 28.03406530         NA         NA 30.10251452
#> 2 Between-stratum variance  2.73828213         NA         NA  1.37921627
#> 3       Between-stratum SD  1.65477555         NA         NA  1.17440039
#> 4                  VPC/ICC  0.05845177         NA         NA  0.03032209
#> 5   PCV (null -> adjusted)          NA         NA         NA  0.82112349
#>   adjusted_lower adjusted_upper
#> 1             NA             NA
#> 2             NA             NA
#> 3             NA             NA
#> 4             NA             NA
#> 5             NA             NA
tab$strata          # all strata ranked by predicted BMI
#>    rank stratum             label    n predicted predicted_lower
#> 1     1       9    female × Black  182  32.12263        31.19779
#> 2     2       4  female × Mexican   99  29.55745        28.34416
#> 3     3       8    male × Mexican  143  29.29069        28.25871
#> 4     4       5      male × White  990  29.17842        28.76803
#> 5     5       1   male × Hispanic   75  28.97383        27.61007
#> 6     6       2      male × Black  154  28.91298        27.91492
#> 7     7       6 female × Hispanic   91  28.62295        27.36517
#> 8     8       3    female × White 1044  28.42387        28.02407
#> 9     9       7      male × Other  111  26.76676        25.61216
#> 10   10      10    female × Other  111  26.35652        25.20192
#>    predicted_upper random_effect   re_lower   re_upper
#> 1         33.04747    3.28538677  2.3605470  4.2102265
#> 2         30.77074    0.80780868 -0.4054814  2.0210987
#> 3         30.32267    0.48681069 -0.5451697  1.5187911
#> 4         29.58881    0.27209916 -0.1382894  0.6824878
#> 5         30.33760    0.19198013 -1.1717834  1.5557437
#> 6         29.91103    0.04292594 -0.9551305  1.0409824
#> 7         29.88074   -0.12717821 -1.3849632  1.1306068
#> 8         28.82367   -0.50097391 -0.9007737 -0.1011741
#> 9         27.92136   -2.01515829 -3.1697567 -0.8605599
#> 10        27.51112   -2.44370096 -3.5982994 -1.2891025

# write.csv(tab$models, "results.csv", row.names = FALSE)
# }
```
