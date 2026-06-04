# Comparing Intersectional Inequality Across Groups

## Introduction

MAIHDA is often used to estimate how much outcome variation lies between
intersectional strata in a single population. The group comparison
workflow asks a related contextual question: is intersectional
inequality larger in some higher-level groups than in others?

Examples include comparing gender-by-class inequality across countries,
race-by-gender inequality across regions, or intersectional inequality
across survey waves. In the MAIHDA package, this is handled by
[`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md) with a
`group` argument or, when you want only the group comparison table, by
[`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md).

## Example data: countries, gender, and socioeconomic status

The package includes `maihda_country_data`, a teaching dataset based on
a balanced subset of OECD PISA 2018 data. Each row is a student, the
outcome is a mathematics score, and the intersectional strata are formed
by `gender` and `ses`.

The higher-level grouping variable is `country`, which lets us ask
whether the VPC/ICC for gender-by-SES inequality differs across
countries.

``` r

library(MAIHDA)
data("maihda_country_data")

country_counts <- as.data.frame(table(maihda_country_data$country))
names(country_counts) <- c("country", "n")
country_counts
#>          country   n
#> 1        Finland 600
#> 2        Germany 600
#> 3 United Kingdom 600
#> 4          Italy 600
#> 5          Japan 600
#> 6         Mexico 600

table(maihda_country_data$gender, maihda_country_data$ses)
#>         
#>          Low Medium High
#>   female 619    571  588
#>   male   581    629  612
```

## One-call workflow

Use [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md) when
you want the usual model summary and the group comparison in one object.
The formula below defines six intersectional strata (`gender:ses`) and
asks for a separate VPC/ICC estimate within each country.

``` r

analysis <- maihda(
  math ~ 1 + (1 | gender:ses),
  data = maihda_country_data,
  group = "country"
)

analysis
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
```

The printed report includes the overall VPC/ICC first, then the
country-level comparison. The group comparison is also stored in
`analysis$groups`, so it can be inspected or saved as a regular data
frame.

``` r

group_results <- as.data.frame(analysis$groups)
group_results[order(group_results$vpc, decreasing = TRUE),
              c("group", "n", "n_strata", "vpc", "var_between",
                "var_residual", "status")]
#>            group   n n_strata        vpc var_between var_residual status
#> 2        Germany 600        6 0.14448319   1271.5810     7529.311     ok
#> 5         Mexico 600        6 0.13649162    771.5419     4881.127     ok
#> 4          Japan 600        6 0.13344144   1032.3346     6703.902     ok
#> 3          Italy 600        6 0.11889899   1065.3186     7894.544     ok
#> 1        Finland 600        6 0.10994297    785.7610     6361.226     ok
#> 6 United Kingdom 600        6 0.06011138    470.5471     7357.372     ok
```

In this example all countries have 600 students and all six
gender-by-SES strata are populated. Higher VPC/ICC values indicate that
a larger share of mathematics score variation lies between the
intersectional strata within that country.

## Visualizing the comparison

The group VPC plot orders countries by the estimated intersectional
VPC/ICC. If the comparison was run with `bootstrap = TRUE`, the same
plot will include per-country confidence intervals.

``` r

plot(analysis, type = "group_vpc")
```

![](group_comparison_files/figure-html/group-vpc-plot-1.png)

The variance-components view shows the same comparison as shares of the
total model variance. The orange segment is the between-stratum
component, which is the numerator of the VPC/ICC.

``` r

plot(analysis, type = "group_components")
```

![](group_comparison_files/figure-html/group-components-plot-1.png)

## Direct group comparison

If you only need the group comparison table and plots, call
[`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
directly. This fits the same stratified per-country models used by
`maihda(group = "country")`.

``` r

group_cmp <- compare_maihda_groups(
  math ~ 1 + (1 | gender:ses),
  data = maihda_country_data,
  group = "country"
)

group_cmp
plot(group_cmp, type = "vpc")
plot(group_cmp, type = "components")
```

By default, `shared_strata = TRUE`. That means the `gender:ses` strata
are defined once on the pooled data, then reused in every country. A
“female:Low” stratum therefore refers to the same type of stratum in
each country, which is what makes the VPC/ICC estimates directly
comparable.

Set `shared_strata = FALSE` only when you intentionally want each group
to rebuild its own strata. In that case, stratum identities are no
longer shared across groups, so interpretation should focus on
within-group structure rather than direct stratum-to-stratum comparison.

## Adding bootstrap intervals

For publication-oriented summaries, you can request per-group parametric
bootstrap confidence intervals with the `lme4` engine. This can take
noticeably longer because each group requires repeated model refits, so
the example is not run when the vignette is built.

``` r

group_cmp_boot <- compare_maihda_groups(
  math ~ 1 + (1 | gender:ses),
  data = maihda_country_data,
  group = "country",
  bootstrap = TRUE,
  n_boot = 500
)

plot(group_cmp_boot, type = "vpc")
```
