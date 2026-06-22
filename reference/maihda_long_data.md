# Simulated Longitudinal Data for MAIHDA

A simulated long-format (one row per person-occasion) panel for
demonstrating **longitudinal / growth-curve MAIHDA**
(`fit_maihda(id =, time =)` and
`maihda(decomposition = "longitudinal")`). 600 individuals are each
measured over five waves on a continuous wellbeing score, within 12
intersectional strata defined by gender x ethnicity x education. The
between-stratum *trajectory* differences are constructed to be mostly
additive (each dimension's main effect on the baseline level and on the
rate of change) with one genuine multiplicative interaction, so the
longitudinal PCV – a high but sub-1 `PCV_slope` – is demonstrable.

## Usage

``` r
maihda_long_data
```

## Format

A data frame with 3000 rows (600 persons x 5 waves) and 8 variables:

- id:

  Person identifier (level 2); repeated across waves.

- wave:

  Measurement occasion, 0 to 4 (the numeric time variable).

- gender:

  Gender (`Women`/`Men`); a stratum dimension.

- ethnicity:

  Ethnicity (`EthA`/`EthB`/`EthC`); a stratum dimension.

- education:

  Education (`Low`/`High`); a stratum dimension.

- age:

  Baseline age in years, a time-invariant covariate.

- wellbeing:

  Continuous wellbeing outcome (the growth-curve response).

- low_wellbeing:

  Binary companion outcome (1 = bottom 40% of wellbeing), for the
  logistic longitudinal path.

## Source

Simulated for the purpose of the MAIHDA package. The growth structure
follows the longitudinal MAIHDA of Bell, Evans, Holman & Leckie (2024)
[doi:10.1016/j.socscimed.2024.116955](https://doi.org/10.1016/j.socscimed.2024.116955)
.

## Examples

``` r
data(maihda_long_data)
# \donttest{
# Time-varying VPC from a 3-level growth model:
m <- fit_maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
                data = maihda_long_data, id = "id", time = "wave")
#> Warning: Model failed to converge with max|grad| = 0.00213664 (tol = 0.002, component 1)
#>   See ?lme4::convergence and ?lme4::troubleshooting.
summary(m)
#> MAIHDA Model Summary
#> ====================
#> 
#> Fit diagnostics:
#>   Convergence warnings reported by lme4:
#>     - Model failed to converge with max|grad| = 0.00213664 (tol = 0.002, component 1)
#>   See ?lme4::convergence and ?lme4::troubleshooting.
#> 
#> 
#> Variance Partition Coefficient (VPC/ICC) at baseline (wave = 0):
#>   Estimate: 0.3986
#> 
#> Variance Components:
#>                                            component  variance     sd
#>                Between-stratum: intercept (time = 0)  0.407878 0.6387
#>                        Between-stratum: slope (wave)  0.042860 0.2070
#>          Between-stratum: intercept-slope covariance  0.088498     NA
#>        Between-individual (id): intercept (time = 0)  0.255609 0.5056
#>                Between-individual (id): slope (wave)  0.019360 0.1391
#>  Between-individual (id): intercept-slope covariance -0.001107     NA
#>                                    Within (residual)  0.359816 0.5998
#> 
#> Time-varying VPC/ICC (between-stratum share over wave):
#>   range 0.3986 to 0.6629 across wave in [0, 4].
#>   The between-stratum variance is a function of time (random intercept +
#>   slope), so the VPC varies; it depends on where time is zeroed. See
#>   plot(type = "vpc_trajectory") for the full curve.
#> 
#> Fixed Effects:
#>         term  estimate
#>  (Intercept)  4.986466
#>         wave -0.006306
#> 
#> Stratum baseline (intercept) deviations (first 10):
#>  stratum stratum_id               label random_effect      se lower_95 upper_95
#>        1          1    Men × EthB × Low      -0.63796 0.08953 -0.81344  -0.4625
#>        2          2  Women × EthB × Low      -0.06459 0.08595 -0.23305   0.1039
#>        3          3   Men × EthA × High       0.86072 0.07783  0.70819   1.0133
#>        4          4   Men × EthB × High       0.09121 0.10260 -0.10988   0.2923
#>        5          5 Women × EthA × High       0.93602 0.08804  0.76345   1.1086
#>        6          6 Women × EthC × High       0.35851 0.13545  0.09302   0.6240
#>        7          7   Men × EthC × High      -0.11460 0.13287 -0.37503   0.1458
#>        8          8    Men × EthA × Low      -0.32848 0.07050 -0.46667  -0.1903
#>        9          9  Women × EthA × Low       0.15992 0.07326  0.01633   0.3035
#>       10         10    Men × EthC × Low      -0.82710 0.13287 -1.08753  -0.5667
#>   ... and 2 more strata
#>   (random slope not shown; use predict(type = "strata") for the per-stratum intercept and slope, or plot(type = "trajectories")).

# Additive-vs-multiplicative PCV (null vs adjusted growth model):
a <- maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
            data = maihda_long_data, id = "id", time = "wave",
            decomposition = "longitudinal")
#> Warning: Model failed to converge with max|grad| = 0.00213664 (tol = 0.002, component 1)
#>   See ?lme4::convergence and ?lme4::troubleshooting.
a$pcv
#> Longitudinal PCV (additive vs. multiplicative intersectionality)
#> ================================================================
#> 
#> Baseline (wave = 0) variance: 0.4079 (null) -> 0.0332 (adjusted)
#>   PCV_intercept: 91.9% of the baseline between-stratum inequality is additive.
#> Slope (wave) variance:          0.0429 (null) -> 0.0058 (adjusted)
#>   PCV_slope:     86.6% of the *trajectory* between-stratum inequality is additive
#>                  (the remainder is the multiplicative/interaction part).
#> 
#> The PCV is the share of the null model's between-stratum (trajectory) variance
#> explained by the dimensions' additive main effects and their time interactions;
#> a high PCV_slope means trajectory inequalities are 'mostly additive'.
# }
```
