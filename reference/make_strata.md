# Create Strata from Multiple Variables

This function creates strata (intersectional categories) from multiple
categorical variables in a dataset.

## Usage

``` r
make_strata(data, vars, sep = "_", min_n = 1)
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

## Value

A list with two elements:

- data:

  The original data frame with an added 'stratum' column

- strata_info:

  A data frame with information about each stratum including counts and
  the combination of variable values

## Examples

``` r
if (FALSE) { # \dontrun{
# Create strata from gender and race variables
data <- data.frame(
  gender = c("M", "F", "M", "F"),
  race = c("White", "Black", "White", "Black"),
  outcome = c(1, 2, 3, 4)
)
result <- make_strata(data, vars = c("gender", "race"))
print(result$strata_info)
} # }
```
