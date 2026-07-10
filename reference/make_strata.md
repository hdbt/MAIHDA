# Create Strata from Multiple Variables

This function creates strata (intersectional categories) from multiple
categorical variables in a dataset.

## Usage

``` r
make_strata(data, vars, sep = " × ", min_n = 1, autobin = TRUE)
```

## Arguments

- data:

  A data frame containing the variables to create strata from.

- vars:

  Character vector of variable names to use for creating strata.

- sep:

  Separator to use between variable values when creating stratum labels.
  Default is " \u00d7 " (a mathematical multiplication sign).

- min_n:

  Minimum number of observations required for a stratum to be included.
  Strata with fewer observations will be coded as NA. Default is 1.

- autobin:

  Logical indicating whether to automatically bin numeric grouping
  variables with more than 10 unique values into 3 categories
  (tertiles). Default is TRUE. When this happens a
  [`message()`](https://rdrr.io/r/base/message.html) is emitted, because
  the resulting strata are data-dependent (tertile cut-points depend on
  the sample) and a continuous variable placed in the grouping term is
  usually unintended. Set `autobin = FALSE` to disable, or bin the
  variable yourself for explicit, reproducible cut-points.

## Value

An object of class `maihda_strata`: a list with elements

- data:

  The original data frame with an added 'stratum' column. The
  strata_info is also attached as an attribute for use by fit_maihda()

- strata_info:

  A data frame with information about each stratum including counts and
  the combination of variable values

- vars:

  The stratum-defining variable names, as supplied.

- sep:

  The label separator used.

- min_n:

  The minimum stratum size applied.

- autobin_info:

  A named list with the `breaks` and `labels` of each auto-binned
  numeric variable (empty when nothing was binned).

## Details

If any of the specified variables has a missing value (NA) for a given
observation, that observation will be assigned to the NA stratum
(stratum = NA), rather than creating a stratum that includes the missing
value.

The strata_info data frame is also attached as an attribute to the data,
which allows fit_maihda() to automatically capture stratum labels for
use in plots and summaries.

When `autobin` discretises a numeric grouping variable `v`, the
adjusted-model and prediction machinery later add an internal factor
column named `.maihda_dim_<v>`; the `.maihda_dim_` prefix is therefore
reserved. `make_strata()` errors if `data` already holds the
`.maihda_dim_<v>` column for a variable it is about to auto-bin, so an
existing user column is never silently overwritten (rename it, or pass
`autobin = FALSE`).

A numeric grouping variable with 10 or fewer unique values (or any
numeric when `autobin = FALSE`) is kept as its raw values: distinct
values still define distinct strata, but the variable is *not* binned.
If such a variable is then used to build an adjusted MAIHDA
([`maihda`](https://hdbt.github.io/MAIHDA/reference/maihda.md)), it
enters the adjusted model as a linear fixed effect rather than
categorical main effects, and
[`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md) warns
(when it has three or more distinct values). Wrap category codes in
[`factor()`](https://rdrr.io/r/base/factor.html) before creating strata
if you want them treated categorically.

## Examples

``` r
# Create strata from gender and race variables
result <- make_strata(maihda_sim_data, vars = c("gender", "race"))
print(result$strata_info)
#>   stratum             label   n gender     race
#> 1       1      Male × White 152   Male    White
#> 2       2    Female × White 150 Female    White
#> 3       3   Male × Hispanic  52   Male Hispanic
#> 4       4 Female × Hispanic  32 Female Hispanic
#> 5       5    Female × Black  50 Female    Black
#> 6       6      Male × Black  44   Male    Black
#> 7       7      Male × Asian  13   Male    Asian
#> 8       8    Female × Asian   7 Female    Asian
```
