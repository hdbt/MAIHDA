# Create Strata from Multiple Variables

This function creates strata (intersectional categories) from multiple
categorical variables in a dataset.

## Usage

``` r
make_strata(data, vars, sep = "_", min_n = 1, autobin = TRUE)
```

## Arguments

- data:

  A data frame containing the variables to create strata from.

- vars:

  Character vector of variable names to use for creating strata.

- sep:

  Separator to use between variable values when creating stratum labels.
  Default is "\_".

- min_n:

  Minimum number of observations required for a stratum to be included.
  Strata with fewer observations will be coded as NA. Default is 1.

- autobin:

  Logical indicating whether to automatically bin numeric grouping
  variables with more than 10 unique values into 3 categories
  (tertiles). Default is TRUE.

## Value

A list with two elements:

- data:

  The original data frame with an added 'stratum' column. The
  strata_info is also attached as an attribute for use by fit_maihda()

- strata_info:

  A data frame with information about each stratum including counts and
  the combination of variable values

## Details

If any of the specified variables has a missing value (NA) for a given
observation, that observation will be assigned to the NA stratum
(stratum = NA), rather than creating a stratum that includes the missing
value.

The strata_info data frame is also attached as an attribute to the data,
which allows fit_maihda() to automatically capture stratum labels for
use in plots and summaries.

## Examples

``` r
# Create strata from gender and race variables
result <- make_strata(maihda_sim_data, vars = c("gender", "race"))
print(result$strata_info)
#>   stratum           label   n gender     race
#> 1       1    Female_Asian   7 Female    Asian
#> 2       2    Female_Black  50 Female    Black
#> 3       3 Female_Hispanic  32 Female Hispanic
#> 4       4    Female_White 150 Female    White
#> 5       5      Male_Asian  13   Male    Asian
#> 6       6      Male_Black  44   Male    Black
#> 7       7   Male_Hispanic  52   Male Hispanic
#> 8       8      Male_White 152   Male    White
```
