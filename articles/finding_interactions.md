# Finding the interactions that matter

## Where is the intersectionality?

The headline VPC tells you *how much* of the variation lies between
intersectional strata; the PCV tells you how much of that is the
*additive* sum of the dimensions’ main effects. Neither tells you
**which strata** are genuinely intersectional – whose outcome departs
from what the additive parts predict. That departure is the stratum
random effect (BLUP) of the **adjusted** model, and
[`maihda_interactions()`](https://hdbt.github.io/MAIHDA/reference/maihda_interactions.md)
packages it into a flag of which strata are credibly non-zero.

``` r

library(MAIHDA)
data("maihda_health_data")

a <- maihda(
  BMI ~ Age + Gender + Race + Education + (1 | Gender:Race:Education),
  data = maihda_health_data
)
```

## `maihda_interactions()`

Pass the [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md)
analysis and the function uses the **adjusted** model automatically
(that is the only model on which the stratum effect isolates the *pure*
interaction). It returns one row per stratum – the interaction BLUP, its
interval, and whether it is flagged – sorted flagged-first, then by
magnitude:

``` r

ints <- maihda_interactions(a)
ints
#> Strata with credibly non-zero intersectional interaction
#> ========================================================
#> 
#> 4 of 50 strata flagged (95% interval; no multiplicity correction).
#> Model: adjusted (two-model); interaction on the link (latent) scale.
#> 
#>  stratum                         label   n interaction     se   lower   upper
#>       26  female × Black × High School  46      1.5953 0.7230  0.1783  3.0122
#>        8   male × White × Some College 328      1.3593 0.3448  0.6836  2.0350
#>        7  female × White × High School 232     -1.0650 0.4016 -1.8521 -0.2779
#>        3 female × White × College Grad 335     -0.9967 0.3415 -1.6660 -0.3274
#>    p_value flagged direction
#>  2.734e-02    TRUE     above
#>  8.056e-05    TRUE     above
#>  8.001e-03    TRUE     below
#>  3.515e-03    TRUE     below
#> 
#> Flagging many strata inflates false positives; for a screening error-rate
#>   story use adjust = "BH" (FDR). Interaction BLUPs are shrunken estimates,
#>   so correction is optional -- see ?maihda_interactions.
```

The number of flagged strata is on the object’s attributes:

``` r

sprintf("%d of %d strata flagged", attr(ints, "n_flagged"), attr(ints, "n_strata"))
#> [1] "4 of 50 strata flagged"
```

Each row carries the `interaction` (BLUP, on the model’s link scale),
its `lower`/`upper` interval, the stratum size `n`, the `direction`
(above/below the additive expectation), and – for the frequentist
engines – the conditional `se` and `p_value`.

### Screening many strata: FDR

With many strata you may want an explicit error-rate control. The
default is `adjust = "none"` (see the caveat below on why), but for a
screen across many strata the **false-discovery rate** (`adjust = "BH"`)
matches the goal – discovery – better than family-wise methods:

``` r

maihda_interactions(a, adjust = "BH")
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

## Highlighting the flagged strata on the plots

Pass `highlight_interactions = TRUE` to
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) and the flagged
strata are ringed and starred on the BLUP-based views. On a
[`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md) analysis
the flags are computed once from the adjusted model and reused:

``` r

plot(a, type = "effect_decomp", highlight_interactions = TRUE)
```

![](finding_interactions_files/figure-html/plot-decomp-1.png)

``` r

plot(a, type = "predicted", highlight_interactions = TRUE)
```

![](finding_interactions_files/figure-html/plot-predicted-1.png)

You can also pass a precomputed `maihda_interactions` object (to reuse a
specific `conf_level`/`adjust`) instead of `TRUE`.

## How to read it – and what not to conclude

- **It must come from the adjusted model.** On a null model the stratum
  effect is the *total* between-stratum deviation (additive +
  interaction), so passing a bare null model is flagged with a warning.
  Passing the
  [`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md)
  analysis handles this for you.
- **It is on the link (latent) scale.** For a logistic model the
  interaction is a log-odds departure; the additive/interaction split is
  only exact on that scale.
- **Treat the flags as exploratory, not confirmatory.** The BLUP is
  already a shrunken (partially pooled) estimate, which protects against
  spurious extremes, so the relevant risk is sign/magnitude error rather
  than family-wise Type I error. That is *why* `adjust = "none"` is the
  default; corrections applied to shrunken estimates are conservative.
  Use `"BH"` when you need a discovery-rate story, and read the flagged
  strata as hypotheses to probe, not confirmed interactions.
- **A singular fit voids the flags.** When the between-stratum variance
  is pinned at the boundary the BLUP SEs are unreliable; the object
  records `attr(., "singular")` and “not flagged” is then not evidence
  of “no interaction” (see the [Bayesian sparse
  vignette](https://hdbt.github.io/MAIHDA/articles/bayesian_sparse_maihda.md)).

## Bayesian evidence (brms)

For a `brms` fit the function uses the **exact posterior tail** rather
than a normal approximation: a credible interval and the probability of
direction `pd = P(BLUP > 0)`. `adjust` is inert (the Bayesian answer is
multiplicity-free).

``` r

a_brms <- maihda(BMI ~ Age + Gender + Race + Education + (1 | Gender:Race:Education),
                 data = maihda_health_data, engine = "brms")
maihda_interactions(a_brms)   # gains a `pd` column; no p-values
```

## See also

- [Reporting MAIHDA
  results](https://hdbt.github.io/MAIHDA/articles/reporting_results.md)
  – tidy output and publication tables.
- [Interpreting MAIHDA plots and
  diagnostics](https://hdbt.github.io/MAIHDA/articles/interpreting_plots.md)
  – the effect-decomposition and predicted views in depth.
- [Bayesian MAIHDA for sparse
  intersections](https://hdbt.github.io/MAIHDA/articles/bayesian_sparse_maihda.md)
  – when ML cannot see a real interaction.

## References

- Evans, C. R., Williams, D. R., Onnela, J. P., & Subramanian, S. V.
  (2018). A multilevel approach to modeling health inequalities at the
  intersection of multiple social identities. *Social Science &
  Medicine*, 203, 64-73.
- Gelman, A., Hill, J., & Yajima, M. (2012). Why we (usually) don’t have
  to worry about multiple comparisons. *Journal of Research on
  Educational Effectiveness*, 5(2), 189-211.
