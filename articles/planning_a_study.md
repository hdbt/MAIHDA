# Planning a MAIHDA analysis

## Before you fit

through those design decisions, with small runnable checks you can do on
your own data first to evaluate dimensions, the number of strata, and
the analytic sample

``` r

library(MAIHDA)
data("maihda_health_data")
```

## Is MAIHDA the right tool?

MAIHDA is for questions of the form *“how much of the variation in an
outcome lies between people’s intersectional social positions, and how
much of that is more than the sum of its parts?”* It is well suited
when:

- you have several **categorical** social dimensions (gender,
  race/ethnicity, education, class, …) whose **joint** categories define
  the strata;
- the outcome is measured at the **individual** level;
- you have enough individuals to populate the cells (see below).

## The central tradeoff: more dimensions means emptier cells

Strata are the **cross-product** of the dimensions, so cell counts fall
off fast as you add dimensions.
[`make_strata()`](https://hdbt.github.io/MAIHDA/reference/make_strata.md)
builds the strata and returns a `strata_info` table of counts you can
inspect *before* modelling:

``` r

s2 <- make_strata(maihda_health_data, vars = c("Gender", "Race"))
nrow(s2$strata_info)                       # number of strata
#> [1] 10
summary(s2$strata_info$n)              # cell-size distribution
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>      75     102     127     300     175    1044
```

Add education and the same sample splits into many more, smaller cells:

``` r

s3 <- make_strata(maihda_health_data, vars = c("Gender", "Race", "Education"))
nrow(s3$strata_info)
#> [1] 50
summary(s3$strata_info$n)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>    1.00   13.25   25.50   60.00   45.50  349.00
sum(s3$strata_info$n < 10)             # how many strata have < 10 people
#> [1] 5
```

Each extra dimension multiplies the number of strata and divides the
people among them. Small cells are not fatal, (partial pooling shrinkage
is exactly what protects MAIHDA against noisy small strata) but they
have consequences (next section). A useful rule: choose the fewest
dimensions that answer your question, and look at the cell-size
distribution before committing.

## What sparse cells do: singular fits

When cells get very small the maximum-likelihood (`lme4`) estimate of
the between-stratum variance can collapse to the boundary ( a singular
fit) and report a VPC of (near) zero with no uncertainty. The package
records this and surfaces it in a “Fit diagnostics” note rather than
letting it pass silently:

``` r

over <- fit_maihda(
  BMI ~ 1 + (1 | Gender:Race:Education),
  data = maihda_health_data[1:60, ]       # deliberately too few people per stratum
)
#> boundary (singular) fit: see help('isSingular')
over
#> MAIHDA Model
#> ============
#> 
#> Engine: lme4 
#> Family: gaussian 
#> Formula: BMI ~ (1 | stratum) 
#> 
#> Fit diagnostics:
#>   Singular fit: at least one variance component is estimated at (or near) zero.
#>     The between-stratum variance and any VPC/PCV derived from it may be unreliable.
#>   Convergence warnings reported by lme4:
#>     - boundary (singular) fit: see help('isSingular')
#> 
#> 
#> Underlying model:
#> Linear mixed model fit by REML ['lmerMod']
#> Formula: BMI ~ (1 | stratum)
#>    Data: data
#> REML criterion at convergence: 386.8857
#> Random effects:
#>  Groups   Name        Std.Dev.
#>  stratum  (Intercept) 0.000   
#>  Residual             6.203   
#> Number of obs: 60, groups:  stratum, 24
#> Fixed Effects:
#> (Intercept)  
#>        28.8  
#> optimizer (nloptwrap) convergence code: 0 (OK) ; 0 optimizer warnings; 1 lme4 warnings
```

If you see a singular-fit note, do not read the VPC as a clean zero. The
solution is to collapse dimensions or categories (fewer, larger cells),
or to use `engine = "brms"`, whose weakly-informative priors regularise
the variance off the boundary and return an honest interval, the subject
of the [Bayesian sparse
vignette](https://hdbt.github.io/MAIHDA/articles/bayesian_sparse_maihda.md).

## Continuous variables and the analytic sample

- **Keep continuous variables out of the strata.** A continuous variable
  in the grouping term gives one stratum per value.
  [`make_strata()`](https://hdbt.github.io/MAIHDA/reference/make_strata.md)
  will auto-bin a numeric dimension into tertiles (with a
  [`message()`](https://rdrr.io/r/base/message.html)), but a continuous
  *covariate* belongs in the fixed part of the formula, not the strata.

## What the summaries can and cannot tell you

| Quantity | Answers | Does *not* answer |
|----|----|----|
| **VPC/ICC** | share of variance between strata | the *amount* of between-stratum variation (a share can rise just because the residual fell) |
| **PCV** | additive share of the between-stratum variance | a causal decomposition; a negative PCV is not proof of hidden inequality |
| **Discriminatory accuracy (AUC/MOR)** | how well strata predict the *individual* outcome | how large the *group* differences are (a high VPC can go with modest AUC) |

## Which engine, which design?

- **`lme4` (default)** – fast frequentist fits for adequately-sized
  cells.
- **`brms`** – Bayesian; preferred when cells are sparse or dimensions
  have few levels (regularising priors, honest intervals).

For extensions beyond the cross-sectional case, see the [crossed random
effects](https://hdbt.github.io/MAIHDA/articles/cross_classified.md)
(dimensions/contexts) and
[longitudinal](https://hdbt.github.io/MAIHDA/articles/longitudinal.md)
vignettes.

## A suggested learning path

1.  [Introduction to
    MAIHDA](https://hdbt.github.io/MAIHDA/articles/introduction.md) –
    the end-to-end workflow.
2.  [Interpreting MAIHDA plots and
    diagnostics](https://hdbt.github.io/MAIHDA/articles/interpreting_plots.md).
3.  [Finding the interactions that
    matter](https://hdbt.github.io/MAIHDA/articles/finding_interactions.md).
4.  [Reporting MAIHDA
    results](https://hdbt.github.io/MAIHDA/articles/reporting_results.md)
    – tidy output and tables.
5.  Specialised designs: [binary
    outcomes](https://hdbt.github.io/MAIHDA/articles/binary_outcomes.md),
    [group
    comparison](https://hdbt.github.io/MAIHDA/articles/group_comparison.md),
    [survey
    weights](https://hdbt.github.io/MAIHDA/articles/design_weighted.md),
    and [Bayesian /
    sparse](https://hdbt.github.io/MAIHDA/articles/bayesian_sparse_maihda.md).

## References

- Evans, C. R., Leckie, G., Subramanian, S. V., Bell, A., & Merlo, J.
  (2024). A tutorial for conducting intersectional multilevel analysis
  of individual heterogeneity and discriminatory accuracy (MAIHDA).
  *SSM - Population Health*, 26,
  101664. 
