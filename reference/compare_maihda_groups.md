# Compare MAIHDA Metrics Across Levels of a Grouping Variable

Fits a separate intercept-only MAIHDA model within each level of a
higher-level grouping variable (for example country, region, or survey
wave) and reports how the variance partition coefficient (VPC/ICC) and
the between-/within-stratum variance components differ across those
groups.

## Usage

``` r
compare_maihda_groups(
  formula,
  data,
  group,
  engine = "lme4",
  family = "gaussian",
  shared_strata = TRUE,
  min_group_n = 30,
  bootstrap = FALSE,
  n_boot = 1000,
  conf_level = 0.95,
  autobin = TRUE,
  ...
)
```

## Arguments

- formula:

  A model formula. Either the shorthand intersectional form
  `outcome ~ covars + (1 | var1:var2)` (strata are built automatically)
  or `outcome ~ covars + (1 | stratum)` when `data` already contains a
  `stratum` column from
  [`make_strata`](https://hdbt.github.io/MAIHDA/reference/make_strata.md).

- data:

  A data frame containing the variables in `formula` and the grouping
  variable.

- group:

  Character string naming the grouping variable in `data` (e.g.
  `"country"`). A separate model is fitted for each non-missing level.

- engine:

  Modeling engine, "lme4" (default) or "brms".

- family:

  Model family. Default "gaussian". As in
  [`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md),
  a binary outcome is auto-detected once on the full data and switched
  to "binomial" (with a warning) so every group uses the same family.

- shared_strata:

  Logical. When TRUE (default) intersectional strata are defined once on
  the full data so that a stratum denotes the same combination in every
  group and VPCs are directly comparable; strata absent from a given
  group are simply unused there. When FALSE, strata are rebuilt
  independently within each group (stratum identities are then not
  comparable across groups).

- min_group_n:

  Minimum number of rows a group must have to be modelled. Smaller
  groups are skipped with a warning. Default 30.

- bootstrap:

  Logical; compute per-group parametric-bootstrap VPC confidence
  intervals. lme4 engine only. Default FALSE.

- n_boot:

  Number of bootstrap samples when `bootstrap = TRUE`. Default 1000.

- conf_level:

  Confidence level for bootstrap intervals. Default 0.95.

- autobin:

  Logical passed to
  [`make_strata`](https://hdbt.github.io/MAIHDA/reference/make_strata.md)
  controlling tertile binning of numeric grouping variables. Default
  TRUE.

- ...:

  Additional arguments passed to
  [`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
  (and on to `lmer`/`glmer`).

## Value

A `data.frame` of class `maihda_group_comparison` with one row per group
and columns `group`, `n`, `n_strata`, `vpc`, `var_between`, `var_other`,
`var_residual`, `status` (and `ci_lower`/`ci_upper` when
`bootstrap = TRUE`). `var_other` is the variance of any additional
random effects and is 0 for the canonical single-stratum model. Groups
that were skipped or failed have `NA` metrics and an explanatory
`status`.

## Details

This answers the question "is intersectional inequality larger in some
groups than others?" by estimating one VPC per group. It is a stratified
analysis: each group is modelled independently. It is *not* a
cross-classified model and does not adjust the strata for the grouping
variable.

Robustness: a group with fewer than `min_group_n` rows is always skipped
with a warning. A group with fewer than two populated strata is also
skipped (VPC is undefined with a single stratum) when the stratum
membership is known before fitting – that is, when
`shared_strata = TRUE` or `data` already carries a `stratum` column.
Under `shared_strata = FALSE` strata are rebuilt inside each group, so a
degenerate single-stratum group is instead reported with a "fit failed"
status rather than a pre-fit skip. A singular fit yields a VPC of 0
rather than an error (unlike
[`calculate_pvc`](https://hdbt.github.io/MAIHDA/reference/calculate_pvc.md)).
A hard fit failure in one group records `NA` and a status note without
aborting the whole comparison.

## See also

[`compare_maihda`](https://hdbt.github.io/MAIHDA/reference/compare_maihda.md)
for comparing different models on the same data;
[`plot_group_comparison`](https://hdbt.github.io/MAIHDA/reference/plot_group_comparison.md)
for visualising the result.

## Examples

``` r
# \donttest{
data(maihda_health_data)
# Using Education as a stand-in higher-level grouping variable
cmp <- compare_maihda_groups(
  BMI ~ Age + (1 | Gender:Race),
  data = maihda_health_data,
  group = "Education"
)
print(cmp)
#> MAIHDA Group Comparison
#> =======================
#> 
#> Group variable: Education 
#> Engine: lme4  | Family: gaussian  | Strata: shared/global 
#> 
#>           group   n n_strata     vpc var_between var_other var_residual status
#>       8th Grade 196       10 0.11236       3.819         0        30.17     ok
#>  9 - 11th Grade 366       10 0.03875       2.004         0        49.72     ok
#>    College Grad 846       10 0.09503       3.308         0        31.50     ok
#>     High School 632       10 0.05250       2.666         0        48.11     ok
#>    Some College 960       10 0.03748       1.978         0        50.80     ok
plot_group_comparison(cmp, type = "vpc")

# }
```
