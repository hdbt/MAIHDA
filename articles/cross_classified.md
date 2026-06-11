# Cross-classified MAIHDA: additive vs. interaction in one model

## Two ways to split intersectional inequality

Intersectional MAIHDA asks how much of the variation in an outcome lies
*between* intersectional strata, and how much of that between-stratum
variation is **additive** (explained by the constituent dimensions’ main
effects) versus a genuinely **intersectional interaction** (over and
above the additive parts).

The package offers two estimators for that split, selected with the
`decomposition` argument of
[`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md):

- **`"two-model"` (default)** – the classic approach. Fit a *null* model
  and an *adjusted* model that adds the dimensions’ additive main
  effects as **fixed** effects, and read the additive share from the
  **PCV** (the proportional change in between-stratum variance). See
  [`vignette("introduction", package = "MAIHDA")`](https://hdbt.github.io/MAIHDA/articles/introduction.md).

- **`"cross-classified"`** – a **single** model that enters each
  dimension’s additive main effect as a **random** intercept:

  ``` math
  y_i = \beta_0 + \mathbf{x}_i\boldsymbol\beta
    + u^{(1)}_{d_1[i]} + \dots + u^{(K)}_{d_K[i]} + u^{(\text{stratum})}_{s[i]} + e_i.
  ```

  Each dimension’s random-effect variance is its additive contribution,
  and the full intersection (`stratum`) random-effect variance is the
  interaction *beyond* additive. The additive and interaction shares of
  the total between-strata variance are read directly from this one fit.

This is exactly the cross-classified random-effects design used in the
published cross-classified MAIHDA studies. It is fit with ordinary
multilevel software – `lme4` for a frequentist fit or `brms` for a
Bayesian one – so it does **not** require any special toolchain.

## Running a cross-classified analysis

``` r

library(MAIHDA)
data("maihda_health_data")

cc <- maihda(
  BMI ~ Age + (1 | Gender:Race:Education),
  data = maihda_health_data,
  decomposition = "cross-classified"
)
#> boundary (singular) fit: see help('isSingular')
cc
#> MAIHDA Analysis
#> ===============
#> 
#> Decomposition:   cross-classified (single model)
#> Formula:         BMI ~ Age + (1 | Gender) + (1 | Race) + (1 | Education) + (1 |      stratum)
#> Engine: lme4 | Family: gaussian
#> Fit diagnostics:
#>   Singular fit: at least one variance component is estimated at (or near) zero.
#>     The between-stratum variance and any VPC/PCV derived from it may be unreliable.
#>   Convergence warnings reported by lme4:
#>     - boundary (singular) fit: see help('isSingular')
#> 
#> 
#> VPC/ICC: 0.0707
#> 
#> Additive vs. Intersectional Decomposition (cross-classified):
#>   Additive (sum of dimension main effects) variance: 2.1996
#>   Intersectional interaction variance:               1.1024
#>   Total between-strata variance:                     3.3020
#>   Additive share of between-strata variance:    66.6%
#>   Interaction share of between-strata variance: 33.4%
#>   Per-dimension additive variance:
#>     Gender: 0.0000
#>     Race: 1.8589
#>     Education: 0.3407
#>   Note: the additive share is the cross-classified analogue of the PCV but a
#>   different estimator; interpret the interaction share cautiously.
#> 
#> Strata: 50
#> 
#> Use summary() for variance components and plot(type = ...) for figures.
```

[`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md) builds
the cross-classified model for you from the intersectional shorthand: it
reads the dimensions (`Gender`, `Race`, `Education`) from the grouping
term, adds one additive random intercept per dimension plus the
intersection random intercept, and fits the single model:

``` r

cc$formula
#> BMI ~ Age + (1 | Gender) + (1 | Race) + (1 | Education) + (1 | 
#>     stratum)
#> <environment: 0x55fb0c6773e8>
```

The partition is on `cc$decomposition` (and printed above):

``` r

cc$decomposition$additive_var        # sum of the dimension random-effect variances
#> [1] 2.199613
cc$decomposition$interaction_var     # the intersection random-effect variance
#> [1] 1.10236
cc$decomposition$additive_share      # additive share of the between-strata variance
#> [1] 0.666151
cc$decomposition$interaction_share   # the complement: the interaction share
#> [1] 0.333849
cc$decomposition$per_dim             # additive variance per dimension
#>       Gender         Race    Education 
#> 6.632123e-08 1.858914e+00 3.406986e-01
```

[`summary()`](https://rdrr.io/r/base/summary.html) shows the full
variance-components table (one row per dimension, the interaction, and
the residual) alongside the decomposition:

``` r

summary(cc$model)
#> MAIHDA Model Summary
#> ====================
#> 
#> Fit diagnostics:
#>   Singular fit: at least one variance component is estimated at (or near) zero.
#>     The between-stratum variance and any VPC/PCV derived from it may be unreliable.
#>   Convergence warnings reported by lme4:
#>     - boundary (singular) fit: see help('isSingular')
#> 
#> 
#> Variance Partition Coefficient (VPC/ICC):
#>   Estimate: 0.0707
#> 
#> Variance Components:
#>                   component  variance        sd proportion
#>            Additive: Gender 6.632e-08 0.0002575  1.420e-09
#>              Additive: Race 1.859e+00 1.3634201  3.981e-02
#>         Additive: Education 3.407e-01 0.5836939  7.297e-03
#>  Intersectional interaction 1.102e+00 1.0499335  2.361e-02
#>   Within-stratum (residual) 4.339e+01 6.5871467  9.293e-01
#>                       Total 4.669e+01 6.8331892  1.000e+00
#> 
#> Additive vs. Intersectional Decomposition (cross-classified):
#>   Additive (sum of dimension main effects) variance: 2.1996
#>   Intersectional interaction variance:               1.1024
#>   Total between-strata variance:                     3.3020
#>   Additive share of between-strata variance:    66.6%
#>   Interaction share of between-strata variance: 33.4%
#>   Per-dimension additive variance:
#>     Gender: 0.0000
#>     Race: 1.8589
#>     Education: 0.3407
#>   Note: the additive share is the cross-classified analogue of the PCV but a
#>   different estimator; interpret the interaction share cautiously.
#> 
#> Fixed Effects:
#>         term estimate
#>  (Intercept) 28.18811
#>          Age  0.01508
#> 
#> Stratum Estimates (first 10):
#>  stratum stratum_id                           label random_effect     se
#>        1          1  male × Hispanic × Some College      -0.12058 0.8578
#>        2          2     male × Black × College Grad       0.16241 0.8783
#>        3          3   female × White × College Grad      -1.13921 0.5571
#>        4          4     male × Hispanic × 8th Grade       0.46491 0.9478
#>        5          5    female × Mexican × 8th Grade       0.98670 0.8241
#>        6          6     male × White × College Grad      -0.13330 0.5588
#>        7          7    female × White × High School      -0.50071 0.5816
#>        8          8     male × White × Some College       1.08841 0.5517
#>        9          9 female × White × 9 - 11th Grade       0.55145 0.6631
#>       10         10 female × Hispanic × High School      -0.06806 0.8681
#>  lower_95 upper_95
#>  -1.80195  1.56079
#>  -1.55907  1.88389
#>  -2.23108 -0.04735
#>  -1.39286  2.32269
#>  -0.62859  2.60200
#>  -1.22856  0.96195
#>  -1.64070  0.63928
#>   0.00712  2.16969
#>  -0.74820  1.85110
#>  -1.76951  1.63340
#>   ... and 40 more strata
```

### Figures

The variance-partition and decomposition figures are
cross-classified-aware: the VPC plot shows one additive slice per
dimension plus the interaction and residual, and the deviation
decomposition splits each stratum’s deviation into its additive
(dimension random effects) and interaction (intersection random effect)
parts.

``` r

plot(cc, type = "vpc")            # per-dimension additive + interaction + residual
plot(cc, type = "effect_decomp")  # additive vs. interaction, per stratum
plot(cc, type = "ternary")        # additive / interaction / uncertainty per stratum
```

## Comparing across a higher-level group

Pass a `group` to decompose within each level of a higher-level variable
– here, how the additive-vs-interaction split differs across countries
in the PISA data:

``` r

data("maihda_country_data")
cc_grp <- maihda(
  math ~ 1 + (1 | gender:ses),
  data = maihda_country_data,
  group = "country",
  decomposition = "cross-classified"
)
plot(cc_grp, type = "group_additive_share")  # additive share by country
plot(cc_grp, type = "group_components")       # additive / interaction / residual
```

## A Bayesian fit

Set `engine = "brms"` for a Bayesian cross-classified fit; the additive
and interaction shares then carry posterior **credible** intervals (no
bootstrap needed). This is the recommended engine when dimensions have
few levels (see the caveats below).

``` r

cc_b <- maihda(
  BMI ~ Age + (1 | Gender:Race:Education),
  data = maihda_health_data,
  engine = "brms",
  decomposition = "cross-classified"
)
cc_b$decomposition$additive_share_ci
```

## Two important caveats

1.  **The additive share is not the PCV.** The two-model PCV and the
    cross-classified additive share target the same idea – how much of
    the between-strata variance is additive – but with **different
    estimators** (fixed main effects across two models vs. a single
    model’s random main-effect variances, which are partially pooled).
    They will be close but not numerically identical; do not treat one
    as a check on the other.

2.  **Few-level dimensions are poorly identified.** A dimension’s
    additive variance is estimated from its handful of levels. A binary
    dimension (e.g. a two-level sex variable) contributes a variance
    estimated from just two groups, which `lme4` will often pin at the
    boundary (a *singular fit*), inflating the additive share and
    shrinking the interaction toward zero. Watch for the singular-fit
    note in the output, and prefer `engine = "brms"`, whose
    weakly-informative priors regularise these variances, when
    dimensions are few-levelled.

## Genuinely cross-classified *data*

The same machinery also serves a genuinely cross-classified data
structure – where units belong to two non-nested higher-level groupings
(for example strata crossed with region, or an age–period–cohort
design). Build the model with
[`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
using a `stratum` column from
[`make_strata()`](https://hdbt.github.io/MAIHDA/reference/make_strata.md)
plus the extra grouping factor, e.g.
`outcome ~ covars + (1 | stratum) + (1 | region)`, and
[`summary()`](https://rdrr.io/r/base/summary.html) reports each grouping
factor’s variance component. The `decomposition = "cross-classified"`
mode of [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md)
is the *additive/interaction* special case of this, where the extra
grouping factors are the strata’s own dimensions.
