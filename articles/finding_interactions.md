# Finding interaction patterns

## Overview

A common follow-up question is which strata depart most clearly from the
additive expectation.

For this exploratory diagnostic,
[`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md) can
compute the adjusted-model stratum random effects, intervals, and flags.
When you use `interactions = "BH"`, the flags use the Benjamini-Hochberg
adjustment. The diagnostic is stored in the analysis object and can be
reused by the effect-decomposition and predicted-strata plots.

## Run a standard analysis and choose the multiplicity rule

The default is `adjust = "none"`. This reports each stratum-level
interval and flag without a multiple-testing correction.

If the goal is to scan all strata and highlight a smaller set for
follow-up, use an adjustment such as Benjamini-Hochberg:

``` r

library(MAIHDA)
model_bh <- maihda(
  BMI ~ Age + Gender + Race + Education + (1 | Gender:Race:Education),
  data = maihda_health_data,
  interactions = "BH" # Benjamini-Hochberg adjustment
)
```

The printed output reports how many strata were flagged and which
adjustment rule was used. The full table is stored in
`model_bh$interactions`.

``` r

model_bh$interactions
#> Strata with credibly non-zero intersectional interaction
#> ========================================================
#> 
#> 1 of 50 strata flagged (95% interval; BH-adjusted p-values).
#> Model: adjusted (two-model); interaction on the link (latent) scale.
#> 
#>  stratum                       label   n interaction     se  lower upper
#>        8 male × White × Some College 328       1.359 0.3448 0.6836 2.035
#>    p_value p_adjusted flagged direction
#>  8.056e-05   0.004028    TRUE     above
#> 
#> Interaction BLUPs are shrunken (partially pooled) estimates; treat flags as
#>   exploratory. See ?maihda_interactions.
```

Each row is one stratum. The main columns are:

- `interaction`: the adjusted-model stratum random effect, on the model
  scale.
- `lower` and `upper`: the interval for that random effect.
- `direction`: whether the stratum is above or below the additive
  expectation.
- `flagged`: whether the stratum passes the selected screening rule.

For frequentist fits, the table also includes the conditional standard
error, p-value, and adjusted p-value when a correction is requested.

## Highlight flagged strata

The plotting methods can reuse the stored diagnostic. Because `model_bh`
was fitted with `interactions = "BH"`, `highlight_interactions = TRUE`
uses the Benjamini-Hochberg flags. In the effect-decomposition plot, the
labels also follow that same flagged set.

``` r

plot(model_bh, type = "effect_decomp", highlight_interactions = TRUE)
```

![](finding_interactions_files/figure-html/plot-decomp-1.png)

The same flags can be reused in the predicted-strata view.

``` r

plot(model_bh, type = "predicted", highlight_interactions = TRUE)
```

![](finding_interactions_files/figure-html/plot-predicted-1.png)

If the analysis was fitted without a stored interaction diagnostic, pass
the adjustment rule directly:

``` r

plot(model, type = "effect_decomp", highlight_interactions = "BH")
```

## See also

- [Introduction to
  MAIHDA](https://hdbt.github.io/MAIHDA/articles/introduction.md), the
  main two-model workflow.
- [Interpreting MAIHDA plots and
  diagnostics](https://hdbt.github.io/MAIHDA/articles/interpreting_plots.md),
  how to read the effect-decomposition and predicted-strata plots.
- [Bayesian MAIHDA for sparse
  intersections](https://hdbt.github.io/MAIHDA/articles/bayesian_sparse_maihda.md),
  sparse cells and uncertainty in variance components.

## References

- Evans, C. R., Williams, D. R., Onnela, J. P., & Subramanian, S. V.
  (2018). A multilevel approach to modeling health inequalities at the
  intersection of multiple social identities. *Social Science &
  Medicine*, 203, 64-73.

- Merlo, J. (2018). Multilevel analysis of individual heterogeneity and
  discriminatory accuracy (MAIHDA) within an intersectional framework.
  *Social Science & Medicine*, 203, 74-80.
