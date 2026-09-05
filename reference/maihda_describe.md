# Describe the MAIHDA sample before fitting

Builds the standard “Table 1” descriptives of a MAIHDA analytic sample
*before* any model is fitted: the total and complete-case sample sizes,
the observed vs. expected intersectional strata (and which are empty or
small), the distribution of each stratum-defining dimension,
family-aware per-stratum outcome summaries, missing-data accounting, and
– for a contextual cross-classified design – the context units. It is
the pre-model counterpart of
[`maihda_table`](https://hdbt.github.io/MAIHDA/reference/maihda_table.md).

The description comes from the *same machinery* as the model: the strata
are built by
[`make_strata`](https://hdbt.github.io/MAIHDA/reference/make_strata.md)
via the same formula-shorthand resolution as
[`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md),
the analytic sample reproduces the engines' own row selection (missing
outcome, covariates, stratum dimensions, context, and – for a weighted
design – missing/non-positive sampling weights), and the outcome family
uses the same auto-detection as the fitters. Stratum IDs, labels, and
counts therefore match a subsequent
[`make_strata()`](https://hdbt.github.io/MAIHDA/reference/make_strata.md)
/
[`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
/ [`maihda`](https://hdbt.github.io/MAIHDA/reference/maihda.md) on the
same formula and data exactly.

## Usage

``` r
maihda_describe(
  x,
  data = NULL,
  context = NULL,
  family = NULL,
  sampling_weights = NULL,
  flag_stratum_n = 20,
  include_empty_strata = TRUE,
  autobin = TRUE,
  digits = 3,
  weights = NULL
)
```

## Arguments

- x:

  A model formula (with the intersectional shorthand
  `outcome ~ covars + (1 | var1:var2)` or `... + (1 | stratum)` after
  [`make_strata`](https://hdbt.github.io/MAIHDA/reference/make_strata.md)),
  or a fitted `maihda_model` / `maihda_analysis` to describe post hoc.

- data:

  A data frame with the model variables; required when `x` is a formula,
  and must be omitted for a fitted-model input.

- context:

  Optional character vector naming higher-level context column(s),
  exactly as in
  [`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md);
  adds the per-unit `$context` table and the contextual identification
  checks. Must be omitted for a fitted-model input (taken from the fit).

- family:

  `NULL` (default) auto-detects the family the same way the fitters do –
  a binary outcome becomes `"binomial"`, an ordered factor with 3+
  levels becomes the cumulative (ordinal) model, anything else is
  described as `"gaussian"`. Otherwise any family specification
  [`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
  accepts (e.g. `"poisson"`, `"negbinomial"`, a family object). Must be
  omitted for a fitted-model input.

- sampling_weights:

  Optional single character string naming a numeric column of `data`
  with individual sampling (design) weights, as in
  [`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md);
  adds weighted counts and outcome summaries. Must be omitted for a
  fitted-model input.

- flag_stratum_n:

  Strata with `n` at or below this threshold are *flagged* as small
  (column `small` of `$strata` and a `$warnings` entry) – never dropped.
  This deliberately differs from
  [`make_strata()`](https://hdbt.github.io/MAIHDA/reference/make_strata.md)'s
  `min_n`, which drops. Default 20.

- include_empty_strata:

  Logical; when `TRUE` (default) the zero-count combinations of the
  observed dimension levels are appended to `$strata` with
  `stratum = NA` and `n = 0`.

- autobin:

  Passed to
  [`make_strata`](https://hdbt.github.io/MAIHDA/reference/make_strata.md)
  when the formula shorthand builds the strata, so the description
  matches a fit with the same `autobin` setting. Ignored for a
  fitted-model input (the strata are already built). Default `TRUE`.

- digits:

  Decimal places used by the
  [`print()`](https://rdrr.io/r/base/print.html) method. Default 3.

- weights:

  Optional precision weights, given as a bare column name of `data`
  exactly as in
  [`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md),
  so the same call describes and fits the same sample. Rows whose weight
  is missing, zero, negative or non-finite are excluded from the
  analytic sample, as the engines exclude them. For a **binomial** model
  the weights also supply the denominator of R's second
  aggregated-binomial idiom (see
  [`?glm`](https://rdrr.io/r/stats/glm.html)), in which the trial counts
  ride in `weights` rather than in a `cbind(successes, failures)`
  response; the outcome is then summarised as events out of trials
  instead of one trial per row. That reading applies when the weights
  are non-unit whole numbers – the same rule
  [`maihda_discriminatory_accuracy`](https://hdbt.github.io/MAIHDA/reference/maihda_discriminatory_accuracy.md)
  uses, so the two never report different sample sizes for one model –
  and it covers a proportion response and a 0/1 frequency-cell response
  alike. Non-integral weights are not trial counts and are left alone.
  The `binomial_weights` overrides available on the AUC have no
  counterpart here: they select an estimand for the concordance, not a
  way to count a sample. Mutually exclusive with `sampling_weights`,
  which are design weights and mean something different. Must be omitted
  for a fitted-model input (the fit already carries its weights). Placed
  last in the argument list for backward compatibility; supply it by
  name.

## Value

An object of class `maihda_describe`: a list of export-ready data frames
(pass to [`write.csv()`](https://rdrr.io/r/utils/write.table.html) or
[`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html)) plus
metadata:

- overview:

  one row: `n_total`, `n_missing_outcome`, `pct_missing_outcome`,
  `n_rows_missing_dimensions` (rows outside every stratum),
  `n_analytic`, `n_strata_observed`, `n_strata_expected`,
  `n_empty_strata`; plus `n_invalid_weights`/`sum_weights` for a
  weighted description and
  `n_individuals`/`median_occasions_per_individual` for a fitted
  longitudinal model

- dimensions:

  counts and percentages for each level of each MAIHDA dimension (`NULL`
  when the dimensions are unknown, e.g. a hand-built `stratum` column)

- strata:

  one row per stratum: `stratum`, `label`, the dimension values, `n`,
  `n_analytic`, `n_missing_outcome`, `pct_missing_outcome`, the
  family-aware outcome summary columns, and the `small`/`empty` flags

- outcome_overall:

  the same family-aware outcome summary over the whole sample

- outcome_levels:

  the outcome's category distribution (binomial / ordinal outcomes;
  `NULL` otherwise)

- context:

  one row per context unit: `context`, `level`, `n`, `n_analytic`, `pct`
  (`NULL` without a context)

- missingness:

  per variable (outcome, each dimension, each context, the weight
  column): `n_missing` and `pct_missing`

- warnings:

  data-quality flags as a data frame of `check` / `message` rows (zero
  rows when clean): auto-binned or ID-like or linear-numeric dimensions,
  empty/small strata, rows lost to missing dimensions, high or
  concentrated outcome missingness, and weakly identified contexts

- observations:

  a slim per-row frame (stratum and the outcome numerator/denominator on
  the summary scale) kept so
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) is
  self-contained

- outcome, family, family_detected, outcome_summary, event_level,
  strata_vars, context_vars, sampling_weights, source, engine,
  longitudinal, response_recoding, call, ...:

  metadata used by [`print()`](https://rdrr.io/r/base/print.html)

## Details

**Family-aware outcome summaries.** The per-stratum (and overall)
outcome summary depends on the resolved family, and the applied summary
type is recorded in `$outcome_summary` – a 0/1 outcome is never
described with a Gaussian mean/SD as if continuous:

- gaussian (and other continuous families): mean, SD, median, min, max;

- binomial (including an aggregated outcome written either as
  `cbind(success, failure)` or, for `brms`, as `success | trials(n)`):
  event count, trials, and the observed proportion. Trials are the
  binomial denominator, so a row with a missing or non-positive trial
  count has no observed outcome;

- poisson / negative binomial: the same numeric summaries, read as the
  observed count mean (rate);

- cumulative (ordinal): the mean and (lower) median category score
  (categories scored 1..K in level order, the same convention as the
  package's plots), with the full category distribution in
  `$outcome_levels`.

**Intersectional sparsity.** `$overview` contrasts `n_strata_observed`
with `n_strata_expected` (the full Cartesian product of the observed
dimension levels) and `n_empty_strata`; with
`include_empty_strata = TRUE` the empty combinations are enumerated as
rows of `$strata` (with `stratum = NA`, so real stratum IDs stay
identical to a subsequent fit). This sparsity is exactly what motivates
the Bayesian engine for small-cell designs (`engine = "brms"`).

**Missing data.** `$overview` separates the rows lost to a missing
stratum dimension
([`make_strata()`](https://hdbt.github.io/MAIHDA/reference/make_strata.md)
routes them to the NA stratum) from rows with a missing outcome, and
reports the resulting complete-case `n_analytic` – computed with the
same row mask the engines apply, so it equals
`nrow(fit_maihda(...)$data)`.

**Weights.** With `sampling_weights`, weighted counts and weighted
outcome means/proportions are reported alongside the unweighted ones,
and rows with a missing or non-positive weight are excluded from the
analytic sample (as the weighted engines do). As elsewhere in the
package, these are population-weighted *point* summaries – not a complex
survey design; no design-based variances are computed.

**Units.** `maihda_describe()` treats rows as the unit of description.
For longitudinal (long-format) data a row is a measurement occasion, not
a person; describing a fitted longitudinal model
(`fit_maihda(id =, time =)`) additionally reports `n_individuals` and
the median occasions per individual in `$overview`.

**Fitted-model input.** A `maihda_model` (or the `maihda_analysis`
bundle from
[`maihda`](https://hdbt.github.io/MAIHDA/reference/maihda.md)) can be
described post hoc; the formula, data, family, context, and sampling
weights are taken from the fit, so the description covers the exact
analytic sample the model used. For the `wemix`/`ordinal` engines the
stored data already *is* the pre-filtered analytic sample, so total and
analytic counts coincide there.

## See also

[`maihda_table`](https://hdbt.github.io/MAIHDA/reference/maihda_table.md)
for the post-model tables,
[`make_strata`](https://hdbt.github.io/MAIHDA/reference/make_strata.md),
[`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md),
[`plot`](https://hdbt.github.io/MAIHDA/reference/plot.maihda_describe.md)
for the descriptive plots.

## Examples

``` r
data(maihda_health_data)
desc <- maihda_describe(BMI ~ Age + (1 | Gender:Race:Education),
                        data = maihda_health_data, flag_stratum_n = 20)
desc
#> MAIHDA Sample Description
#> =========================
#> 
#> Outcome: BMI | Family: gaussian | Summary: mean/SD (gaussian)
#> Dimensions: Gender, Race, Education
#> 
#> Sample:
#>   Rows (total):                    3,000
#>   Missing outcome:                 0 (0.0%)
#>   Outside strata (missing dims):   0
#>   Analytic sample (complete case): 3,000
#> 
#> Observed outcome (3,000 non-missing): mean 28.883 (SD 6.743), median 27.865, range 15.800 to 81.250
#> 
#> Strata (50 observed / 50 expected, 0 empty):
#>   Cell sizes: min 1, median 25.5, max 349; 21 strata at/below n = 20 (flagged)
#>   Smallest strata:
#>     female × Black × 8th Grade             n = 1
#>     female × Mexican × College Grad        n = 5
#>     female × Other × 9 - 11th Grade        n = 7
#>     male × Black × 8th Grade               n = 8
#>     female × Other × High School           n = 9
#>   (full table in $strata)
#> 
#> Dimensions:
#>   Gender: female 1,527 (50.9%), male 1,473 (49.1%)
#>   Race: Black 336 (11.2%), Hispanic 166 (5.5%), Mexican 242 (8.1%), White
#>     2,034 (67.8%), Other 222 (7.4%)
#>   Education: 8th Grade 196 (6.5%), 9 - 11th Grade 366 (12.2%), High School
#>     632 (21.1%), Some College 960 (32.0%), College Grad 846 (28.2%)
#> 
#> Warnings:
#>   ! 21 strata have n <= 20; smallest: 'female × Black × 8th Grade' (n=1),
#>     'female × Mexican × College Grad' (n=5), 'female × Other × 9 - 11th
#>     Grade' (n=7). Partial pooling shrinks small-stratum estimates toward
#>     the mean; expect wide intervals for these cells.
desc$overview
#>   n_total n_missing_outcome pct_missing_outcome n_rows_missing_dimensions
#> 1    3000                 0                   0                         0
#>   n_analytic n_strata_observed n_strata_expected n_empty_strata
#> 1       3000                50                50              0
head(desc$strata)
#>   stratum                          label Gender     Race    Education   n
#> 1       1 male × Hispanic × Some College   male Hispanic Some College  25
#> 2       2    male × Black × College Grad   male    Black College Grad  20
#> 3       3  female × White × College Grad female    White College Grad 335
#> 4       4    male × Hispanic × 8th Grade   male Hispanic    8th Grade  10
#> 5       5   female × Mexican × 8th Grade female  Mexican    8th Grade  33
#> 6       6    male × White × College Grad   male    White College Grad 325
#>   n_analytic n_missing_outcome pct_missing_outcome outcome_mean outcome_sd
#> 1         25                 0                   0     28.42960   4.591365
#> 2         20                 0                   0     30.27050   5.229145
#> 3        335                 0                   0     26.99206   6.344460
#> 4         10                 0                   0     31.30700   5.657263
#> 5         33                 0                   0     31.58515   5.430962
#> 6        325                 0                   0     28.18175   4.953313
#>   outcome_median outcome_min outcome_max small empty
#> 1          28.01       19.50       35.65 FALSE FALSE
#> 2          28.45       25.10       45.80  TRUE FALSE
#> 3          25.60       17.50       48.91 FALSE FALSE
#> 4          30.87       24.80       42.40  TRUE FALSE
#> 5          30.70       21.07       43.70 FALSE FALSE
#> 6          27.40       18.55       44.42 FALSE FALSE

# Binary outcome: auto-detected binomial, proportions instead of means
desc_bin <- maihda_describe(Obese ~ (1 | Gender:Race:Education),
                            data = maihda_health_data)
desc_bin$outcome_levels
#>   level    n  pct
#> 1    No 1923 64.1
#> 2   Yes 1077 35.9

# Contextual cross-classified design (country as context)
data(maihda_country_data)
desc_ctx <- maihda_describe(math ~ escs + (1 | gender:ses),
                            data = maihda_country_data, context = "country")
desc_ctx$context
#>   context          level   n n_analytic      pct
#> 1 country        Finland 600        600 16.66667
#> 2 country        Germany 600        600 16.66667
#> 3 country United Kingdom 600        600 16.66667
#> 4 country          Italy 600        600 16.66667
#> 5 country          Japan 600        600 16.66667
#> 6 country         Mexico 600        600 16.66667

# Export-ready tables
# write.csv(desc$strata, "table1_strata.csv", row.names = FALSE)
```
