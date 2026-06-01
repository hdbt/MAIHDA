# Run a Complete MAIHDA Analysis

A single high-level entry point that runs the standard MAIHDA workflow
and returns one bundled object: it fits the multilevel model, summarises
the variance partition (VPC/ICC) and components, and – when a
higher-level grouping variable is supplied – also compares
intersectional inequality across that variable's levels.

## Usage

``` r
maihda(
  formula,
  data,
  group = NULL,
  engine = "lme4",
  family = "gaussian",
  autobin = TRUE,
  shared_strata = TRUE,
  min_group_n = 30,
  bootstrap = FALSE,
  n_boot = 1000,
  conf_level = 0.95,
  ...
)
```

## Arguments

- formula:

  A model formula, using either the intersectional shorthand
  `outcome ~ covars + (1 | var1:var2)` or `... + (1 | stratum)` when
  `data` already has a `stratum` column from
  [`make_strata`](https://hdbt.github.io/MAIHDA/reference/make_strata.md).

- data:

  A data frame with the model variables (and the `group` variable, if
  used).

- group:

  Optional character string naming a higher-level grouping variable
  (e.g. `"country"`). When supplied,
  [`compare_maihda_groups`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
  is run and attached to the result.

- engine:

  Modeling engine, "lme4" (default) or "brms".

- family:

  Model family. Default "gaussian". As in
  [`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md),
  a binary outcome is auto-detected when `family` is left at the
  default, and the same resolved family is then used for the group
  comparison so all models agree.

- autobin:

  Logical passed to
  [`make_strata`](https://hdbt.github.io/MAIHDA/reference/make_strata.md);
  tertile-bins numeric grouping variables. Default TRUE.

- shared_strata:

  Logical, forwarded to
  [`compare_maihda_groups`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
  when `group` is supplied: build strata once on the full data so VPCs
  are comparable across groups (TRUE, default) or rebuild them within
  each group.

- min_group_n:

  Minimum group size for the per-group comparison, forwarded to
  [`compare_maihda_groups`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md).
  Default 30.

- bootstrap:

  Logical; compute parametric-bootstrap VPC confidence intervals (lme4
  only) for both the overall summary and the per-group comparison.
  Default FALSE.

- n_boot:

  Number of bootstrap samples when `bootstrap = TRUE`.

- conf_level:

  Confidence level for bootstrap intervals. Default 0.95.

- ...:

  Additional arguments passed to
  [`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
  (and on to `lmer`/`glmer`).

## Value

An object of class `maihda_analysis`: a list with

- model:

  the fitted `maihda_model` (see
  [`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md))

- summary:

  the `maihda_summary` (VPC/ICC, variance components, stratum estimates)

- groups:

  a `maihda_group_comparison` when `group` is supplied, otherwise `NULL`

- formula, group_var, call:

  bookkeeping for printing

## Details

This is a convenience wrapper around
[`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md),
[`summary.maihda_model`](https://hdbt.github.io/MAIHDA/reference/summary.maihda_model.md)
and
[`compare_maihda_groups`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md);
it always returns the same `maihda_analysis` structure (the `groups`
slot is simply `NULL` when `group` is not given), so downstream code
never has to branch on the return type.

## See also

[`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
for the single-model fitter,
[`compare_maihda_groups`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
for the group comparison, and
[`summary.maihda_model`](https://hdbt.github.io/MAIHDA/reference/summary.maihda_model.md)
for the variance summary.

## Examples

``` r
# \donttest{
data(maihda_health_data)

# One call: fit + VPC summary
a <- maihda(BMI ~ Age + (1 | Gender:Race), data = maihda_health_data)
a
#> MAIHDA Analysis
#> ===============
#> 
#> Formula: BMI ~ Age + (1 | stratum) 
#> Engine: lme4 | Family: gaussian
#> VPC/ICC: 0.0585
#> Strata: 10
#> 
#> Use summary() for variance components and plot(type = ...) for figures.
plot(a, type = "vpc")


# Add a higher-level grouping variable to also compare across its levels.
# maihda_country_data has a real country grouping (PISA achievement data):
data(maihda_country_data)
a2 <- maihda(math ~ 1 + (1 | gender:ses), data = maihda_country_data,
             group = "country")
a2
#> MAIHDA Analysis
#> ===============
#> 
#> Formula: math ~ (1 | stratum) 
#> Engine: lme4 | Family: gaussian
#> VPC/ICC: 0.1493
#> Strata: 6
#> 
#> Group comparison by 'country':
#> MAIHDA Group Comparison
#> =======================
#> 
#> Group variable: country 
#> Engine: lme4  | Family: gaussian  | Strata: shared/global 
#> 
#>           group   n n_strata     vpc var_between var_other var_residual status
#>         Finland 600        6 0.10994       785.8         0         6361     ok
#>         Germany 600        6 0.14448      1271.6         0         7529     ok
#>           Italy 600        6 0.11890      1065.3         0         7895     ok
#>           Japan 600        6 0.13344      1032.3         0         6704     ok
#>          Mexico 600        6 0.13649       771.5         0         4881     ok
#>  United Kingdom 600        6 0.06011       470.5         0         7357     ok
#> 
#> Use summary() for variance components and plot(type = ...) for figures.
plot(a2, type = "group_vpc")

# }
```
