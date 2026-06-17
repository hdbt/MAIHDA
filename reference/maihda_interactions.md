# Flag strata with credibly non-zero intersectional interaction

Reports, for each intersectional stratum, the **interaction** component
of its outcome – the stratum random effect (BLUP) of an *adjusted*
MAIHDA model, i.e. how far the stratum departs from the additive
main-effects prediction of its defining dimensions – and **flags** the
strata whose interaction is credibly different from zero. This is the
heart of "where is there genuine intersectionality": a flagged stratum
is one whose joint identity produces an outcome the additive parts do
not.

## Usage

``` r
maihda_interactions(object, conf_level = 0.95, adjust = "none", ...)
```

## Arguments

- object:

  A `maihda_analysis` from
  [`maihda`](https://hdbt.github.io/MAIHDA/reference/maihda.md)
  (preferred – its adjusted / crossed-dimensions model is used
  automatically) or a `maihda_model` from
  [`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
  (which should be the *adjusted* model; a null model is accepted but
  warned about).

- conf_level:

  Confidence / credible level for the interval and the flag. Default
  0.95.

- adjust:

  Multiple-comparison adjustment for the per-stratum p-values
  (frequentist engines only): `"none"` (default) or any method accepted
  by [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html) (e.g. `"BH"`,
  `"holm"`, `"bonferroni"`). Ignored for `brms` (with a message), which
  uses the posterior tail directly.

- ...:

  Currently unused.

## Value

An object of class `maihda_interactions` (a data frame), one row per
stratum, sorted flagged-first then by `abs(interaction)`. Columns common
to every engine: `stratum`, `label`, `n` (stratum size), `interaction`
(the BLUP), `lower`/`upper` (the interval), `flagged` (logical), and
`direction` (`"above"`/`"below"` the additive expectation). Frequentist
fits add `se` and `p_value` (and `p_adjusted` when `adjust != "none"`);
`brms` adds `pd` (probability of direction). Attributes record
`conf_level`, `adjust`, `engine`, `model_type`, `n_strata`, `n_flagged`,
`scale` and `singular`.

## Details

**It must be read off the adjusted model.** Only when the dimensions'
additive main effects are in the model (the *adjusted* model of the
two-model decomposition, or the crossed-dimensions model) does the
stratum random effect isolate the *pure interaction*. On a null model
the stratum random effect is the total between-stratum deviation
(additive + interaction), so passing one is flagged with a warning.
Passing a [`maihda`](https://hdbt.github.io/MAIHDA/reference/maihda.md)
result uses the right model automatically.

**Frequentist vs. Bayesian evidence.** For the frequentist engines
(`lme4`, `wemix`, `ordinal`) the flag comes from the BLUP's conditional
standard error: a Wald interval at `conf_level` and a two-sided p-value,
with an optional multiplicity correction (`adjust`). For `brms` the full
posterior is already available, so the *exact* posterior tail is used –
a credible interval at `conf_level` and the probability of direction
`pd = P(BLUP > 0)` – and `adjust` is not applied (the Bayesian answer is
multiplicity-free).

**Multiplicity: partial pooling and a correction are different things,
and the experts disagree.**

- *Shrinkage (magnitude/sign).* The stratum BLUP is partially pooled, so
  extreme values are regularised toward the grand mean, attenuating
  exaggerated-magnitude and wrong-sign (Type M/S) error (Gelman & Carlin
  2014). Gelman, Hill & Yajima (2012) argue this shrinkage *usually
  substitutes* for a classical multiple-comparisons correction (the
  problem can "disappear entirely" in the hierarchical model); on that
  view the flag/no-flag step itself is what to avoid – the null of an
  *exactly* zero interaction is rarely the question (McShane, Gelman et
  al. 2019) – so report the estimate and its interval.

- *Whether to correct.* If you do want an error-rate screen, whether a
  correction is warranted depends on the *inferential structure* of the
  claim – the joint hypothesis, not the number of strata (Rubin 2021).
  Each stratum as its own pre-specified hypothesis ("does *this* stratum
  interact?") is *individual* testing and needs none – **only** if you
  do not also read the flags collectively. Once the question is "is
  there an interaction *somewhere*?" – which an automated all-strata
  scan effectively is – it is *disjunction* testing and a correction
  applies.

`adjust = "none"` is the default because the table is formally a set of
individual hypotheses; **for the common exploratory scan of all strata,
prefer `adjust = "BH"`**. Choosing FDR over family-wise
(`"bonferroni"`/`"holm"`) matches a screening goal (the expected
*proportion* of false discoveries) – this is the package's choice, not a
recommendation of Rubin (2021), who raises FDR only to distinguish it
from the family-wise rate. The flag itself is a Wald test on a shrunken
BLUP whose conditional SE treats the variance components as known, so it
(and any `adjust` on it) is an explicit, approximate *screen*, not a
procedure inheriting an exact guarantee from the model. Lead with the
interval (and, for `brms`, the probability of direction); the
substantive question is often not whether an interaction differs from
zero but whether it exceeds a smallest interaction of interest (an
equivalence/SESOI reading; Lakens, Scheel & Isager 2018), read from the
interval.

The interaction is reported on the model's link (latent) scale – a
log-odds deviation for a logistic model, etc. – because the
additive/interaction split is only exact there.

## References

Evans, C. R., Williams, D. R., Onnela, J. P., & Subramanian, S. V.
(2018). A multilevel approach to modeling health inequalities at the
intersection of multiple social identities. *Social Science & Medicine*,
203, 64-73.

Merlo, J. (2018). Multilevel analysis of individual heterogeneity and
discriminatory accuracy (MAIHDA) within an intersectional framework.
*Social Science & Medicine*, 203, 74-80.

Gelman, A., Hill, J., & Yajima, M. (2012). Why we (usually) don't have
to worry about multiple comparisons. *Journal of Research on Educational
Effectiveness*, 5(2), 189-211.

Gelman, A., & Carlin, J. (2014). Beyond power calculations: assessing
Type S (sign) and Type M (magnitude) errors. *Perspectives on
Psychological Science*, 9(6), 641-651.

Rubin, M. (2021). When to adjust alpha during multiple testing: a
consideration of disjunction, conjunction, and individual testing.
*Synthese*, 199(3-4), 10969-11000.
[doi:10.1007/s11229-021-03276-4](https://doi.org/10.1007/s11229-021-03276-4)

McShane, B. B., Gal, D., Gelman, A., Robert, C., & Tackett, J. L.
(2019). Abandon statistical significance. *The American Statistician*,
73(sup1), 235-245.

Lakens, D., Scheel, A. M., & Isager, P. M. (2018). Equivalence testing
for psychological research: a tutorial. *Advances in Methods and
Practices in Psychological Science*, 1(2), 259-269.

## See also

[`maihda`](https://hdbt.github.io/MAIHDA/reference/maihda.md),
[`calculate_pvc`](https://hdbt.github.io/MAIHDA/reference/calculate_pvc.md),
[`summary.maihda_model`](https://hdbt.github.io/MAIHDA/reference/summary.maihda_model.md);
and `plot(..., highlight_interactions = TRUE)` to mark the flagged
strata on the effect-decomposition / predicted / shrinkage plots.

## Examples

``` r
# \donttest{
data(maihda_health_data)
a <- maihda(BMI ~ Age + Gender + Race + (1 | Gender:Race),
            data = maihda_health_data)
maihda_interactions(a)                 # which strata interact (95%, no correction)
#> Strata with credibly non-zero intersectional interaction
#> ========================================================
#> 
#> 4 of 10 strata flagged (95% interval; no multiplicity correction).
#> Model: adjusted (two-model); interaction on the link (latent) scale.
#> 
#>  stratum          label    n interaction     se   lower   upper  p_value
#>        2   male × Black  154     -1.2902 0.4870 -2.2446 -0.3357 0.008067
#>        9 female × Black  182      1.2902 0.4540  0.4003  2.1800 0.004487
#>        3 female × White 1044     -0.6003 0.2025 -0.9972 -0.2035 0.003025
#>        5   male × White  990      0.6003 0.2077  0.1932  1.0075 0.003855
#>  flagged direction
#>     TRUE     below
#>     TRUE     above
#>     TRUE     below
#>     TRUE     above
#> 
#> Flagging many strata inflates false positives; for a screening error-rate
#>   story use adjust = "BH" (FDR). Interaction BLUPs are shrunken estimates,
#>   so correction is optional -- see ?maihda_interactions.
maihda_interactions(a, adjust = "BH")  # FDR-controlled screen
#> Strata with credibly non-zero intersectional interaction
#> ========================================================
#> 
#> 4 of 10 strata flagged (95% interval; BH-adjusted p-values).
#> Model: adjusted (two-model); interaction on the link (latent) scale.
#> 
#>  stratum          label    n interaction     se   lower   upper  p_value
#>        2   male × Black  154     -1.2902 0.4870 -2.2446 -0.3357 0.008067
#>        9 female × Black  182      1.2902 0.4540  0.4003  2.1800 0.004487
#>        3 female × White 1044     -0.6003 0.2025 -0.9972 -0.2035 0.003025
#>        5   male × White  990      0.6003 0.2077  0.1932  1.0075 0.003855
#>  p_adjusted flagged direction
#>     0.02017    TRUE     below
#>     0.01496    TRUE     above
#>     0.01496    TRUE     below
#>     0.01496    TRUE     above
#> 
#> Interaction BLUPs are shrunken (partially pooled) estimates; treat flags as
#>   exploratory. See ?maihda_interactions.
# }
```
