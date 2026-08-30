# Flag strata with credibly non-zero intersectional interaction

Reports, for each intersectional stratum, the **interaction** component
of its outcome – the stratum random effect (BLUP) of an *adjusted*
MAIHDA model, i.e. how far the stratum departs from the additive
main-effects prediction of its defining dimensions – and **flags** the
strata whose interaction is credibly different from zero. This is the
heart of "where is there genuine intersectionality": a flagged stratum
is one whose outcome departs credibly from what the additive parts
predict – a descriptive statement about the stratum's outcome, not a
causal claim about identity.

## Usage

``` r
maihda_interactions(
  object,
  conf_level = 0.95,
  adjust = "BH",
  rope = NULL,
  scale = c("link", "response"),
  ...
)
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
  (frequentist engines only): `"BH"` (default; false-discovery rate) or
  any method accepted by
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html), including
  `"none"` for the uncorrected, per-stratum individual-testing view.
  Ignored for `brms` (which uses the posterior tail directly; a message
  is shown only if you set it explicitly).

- rope:

  Optional equivalence region (a "smallest interaction of interest") for
  an "is the interaction *negligible*?" reading (Schuirmann 1987;
  Kruschke 2018), read on the requested `scale` – log-odds under the
  default, probability points under `scale = "response"` for a logistic
  fit. `NULL` (default) gives only the usual zero-centred flag. A single
  positive number `d` means the symmetric region `c(-d, d)`; or supply
  `c(lower, upper)`. When set, the result gains a `decision` column
  classifying each stratum from its `conf_level` interval relative to
  the region: `"relevant"` (interval entirely outside it),
  `"negligible"` (entirely inside it), or `"inconclusive"` (straddling a
  bound).

- scale:

  Scale the interaction is reported on. `"link"` (default) gives the
  stratum BLUP \\u_j\\ on the model's link scale – for a Gaussian fit
  the outcome's own units, for a logistic fit a log-odds departure.
  `"response"` instead gives Evans et al.'s (2024, sec. 2.5.1) \\\pi^B_j
  = \pi_j - \pi^A_j\\: the stratum's total predicted outcome minus the
  outcome implied by its additive main effects alone, so a logistic fit
  reports a difference in *probability*, a count fit a difference in
  expected count, and a cumulative (ordinal) fit a difference in
  expected category score – the printed column is named for whichever it
  is. The two scales flag the same strata (see the section below); an
  identity-link fit returns the same numbers either way. For a
  `decomposition = "crossed-dimensions"` model the additive baseline is
  the dimension random effects, matching that decomposition.

- ...:

  Currently unused.

## Value

An object of class `maihda_interactions` (a data frame), one row per
stratum, sorted flagged-first then by `abs(interaction)`. Columns common
to every engine: `stratum`, `label`, `n` (stratum size), `interaction`
(the BLUP, under the default `scale`), `lower`/`upper` (the interval),
`flagged` (logical), and `direction` (`"above"`/`"below"` the additive
expectation). Under `scale = "response"` the
`interaction`/`lower`/`upper` columns hold \\\pi^B_j\\ and its interval,
and `se` is dropped: the conditional standard error is a link-scale
quantity, and the response-scale interval is deliberately not symmetric
about the estimate, so no single SE would reproduce it. Frequentist fits
add `se` (link scale only) and `p_value` (and `p_adjusted` when
`adjust != "none"`). `p_value` is a *conditional* screening statistic –
a Wald tail on the shrunken BLUP's conditional SE, with the variance
components treated as known – **not** a calibrated frequentist p-value;
it is *conservative* (stochastically large under a true null), so the
BH-adjusted flag under-flags truly-null strata rather than delivering an
exact error-rate guarantee (see Details). `brms` instead adds `pd`
(probability of direction, `max(P(>0), P(<0))` in `[0.5, 1]`). When
`rope` is set, a `decision` column
(`"relevant"`/`"negligible"`/`"inconclusive"`) is added. Attributes
record `conf_level`, `adjust`, `rope`, `engine`, `model_type`,
`n_strata`, `n_flagged`, `scale`, `link` (the model's link name, so a
caller can tell a log-odds departure from one in the outcome's own
units), `response_kind` (what the model's response scale is –
`"probability"`, `"count"`, `"score"` or `"response"` – which the link
alone does not settle, since a cumulative fit is logit-linked but scores
categories), `singular`, and – on a non-identity link under the default
`scale` – `response_interaction`, the outcome-scale estimates keyed by
stratum. That last attribute exists so
[`print()`](https://rdrr.io/r/base/print.html) can show the
outcome-scale size beside the link-scale one; the columns are the same
either way, and `scale = "response"` is how to obtain those numbers,
with their interval, as data.

## Details

**It must be read off the adjusted model.** Only when the dimensions'
additive main effects are in the model (the *adjusted* model of the
two-model decomposition, or the crossed-dimensions model) does the
stratum random effect isolate the *pure interaction*. On a null model
the stratum random effect is the total between-stratum deviation
(additive + interaction), so passing one is flagged with a warning. The
opposite mis-specification is flagged too: a model that adds a *fixed*
interaction among the dimensions (e.g. `var1 * var2`) absorbs the
intersectional effect into fixed cell means, so the stratum random
effect is no longer the pure interaction. Passing a
[`maihda`](https://hdbt.github.io/MAIHDA/reference/maihda.md) result
uses the right model automatically.

**Frequentist vs. Bayesian evidence.** For the frequentist engines
(`lme4`, `wemix`, `ordinal`) the flag comes from the BLUP's conditional
standard error: a Wald interval at `conf_level` and a two-sided p-value,
with an optional multiplicity correction (`adjust`). For `brms` the full
posterior is already available, so the *exact* posterior tail is used –
a credible interval at `conf_level` and the probability of direction
`pd = max(P(BLUP > 0), P(BLUP < 0))` (in `[0.5, 1]`; the sign is in
`direction`) – and `adjust` is not applied: the posterior is already
partially pooled, the hierarchical-shrinkage answer to multiplicity
(Gelman, Hill & Yajima 2012). That is a per-stratum answer; the marginal
intervals carry no formal joint error-rate guarantee for the collective
claim "an interaction exists *somewhere*" (the disjunction reading below
applies to Bayesian flags too).

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
  claim. Each stratum as its own pre-specified hypothesis ("does *this*
  stratum interact?") is *individual* testing and needs none – **only**
  if you do not also read the flags collectively. Once the question is
  "is there an interaction *somewhere*?" – which an automated all-strata
  scan effectively is – it is *disjunction* testing and a correction
  applies.

**`adjust = "BH"` is the default**: fitting and flagging every stratum
in one call is the disjunction/screening case, where controlling the
expected *proportion* of false discoveries (FDR) is the appropriate
goal. Pass `adjust = "none"` only when each stratum is a genuine,
pre-specified individual hypothesis. The FDR choice (over family-wise
`"bonferroni"`/`"holm"`) is this package's, matching that screening
goal. The flag itself is a Wald test on a shrunken BLUP whose
conditional SE treats the variance components as known, so it (and any
`adjust` on it) is an explicit, *conservative* screen – strict, not
liberal. Partial pooling deflates a truly-null stratum's BLUP more than
its conditional SE, so the null z-statistic is sub-normal (variance
about the shrinkage fraction, below 1) and the screen under-flags rather
than over-flags, degenerating to no flags at the singular/zero-variance
boundary. Lead with the interval (and, for `brms`, the probability of
direction); the substantive question is often not whether an interaction
differs from zero but whether it exceeds a smallest interaction of
interest (an equivalence reading; Schuirmann 1987; Kruschke 2018), read
from the interval.

The interaction is reported on the model's link (latent) scale – a
log-odds deviation for a logistic model, etc. – because the
additive/interaction split is only exact there. On a non-identity link
that makes it a *multiplicative* departure, which is a different
quantity from a shift in the outcome itself; see the section below
before reporting it.

## What the interaction means on a non-identity link

For a Gaussian (identity-link) fit the stratum BLUP \\u_j\\ *is* the
part of the stratum's mean outcome attributable to interaction, in the
outcome's own units. For a logistic fit it is not: \\u_j\\ is a
deviation in **log-odds**, so what this function reports is the
*multiplicative*, odds-scale interaction. (What follows is written for
the logistic case, the one the literature works in. The same logic holds
on any non-identity link, in its own units: a Poisson fit's \\u_j\\ is a
log-rate departure and `scale = "response"` gives it as a difference in
expected counts, and a cumulative (ordinal) fit's is a latent shift
reported as a difference in expected category score. The printed column
is named for whichever applies.) Evans et al. (2024, section 2.5.1) put
it directly – in a logistic MAIHDA one can "no longer directly interpret
\\u_j\\ as the change in mean outcomes (i.e., shift in probabilities)
attributable to interaction effects".

**Three quantities, of which this function reports two.** Keeping them
apart is the whole of the difficulty:

1.  **The multiplicative interaction, \\u_j\\** – the departure from
    additivity *in log-odds*, i.e. whether a dimension multiplies the
    odds by the same factor at every level of the others. This is what
    the default `scale = "link"` reports and what `flagged` is about.

2.  **The same departure in outcome units, \\\pi^B_j\\** – what
    `scale = "response"` returns. Evans et al. (2024) define it as
    \$\$\pi^B_j = \pi_j - \pi^A_j, \qquad \pi_j =
    \mathrm{logit}^{-1}(x_j'\beta + u_j), \qquad \pi^A_j =
    \mathrm{logit}^{-1}(x_j'\beta),\$\$ the gap between a stratum's
    total predicted probability and the probability implied by the
    additive main effects alone, and rank-plot \\\pi^B_j\\ where the
    linear case plots \\u_j\\. This is quantity 1 *re-expressed*, not a
    second finding: it is zero exactly when \\u_j\\ is zero and it flags
    the same strata, so it changes the units you report, not what you
    may conclude.

3.  **The additive (risk-difference) interaction** – whether a dimension
    adds the same number of percentage points of risk at every level of
    the others. **Neither of the above reports this**, and it is
    generally non-zero even where \\u_j\\ is exactly zero, because the
    logistic curve is steeper in the middle than in the tails. With a
    \\-2\\ baseline and \\+0.7\\ for each of two dimensions and no
    interaction at all (\\u_j = 0\\ throughout), the second dimension
    still adds 9.5 percentage points of risk at one level of the first
    and 14.0 at the other. That excess is a property of the link, not a
    finding.

So "no strata flagged" supports "no credible **multiplicative**
interaction", not "no interaction": additivity in log-odds does not
carry over to probabilities, so quantity 3 is generally non-zero
regardless, and it is often the one a policy audience cares about.
(*Generally*, not always. The risk-difference interaction is positive
below the curve's midpoint and negative above it, so it passes through
zero: with two dimensions of effect \\a\\ and \\b\\ it vanishes exactly
at an intercept of \\-(a + b)/2\\, where the two comparisons straddle
the midpoint symmetrically. That is a single configuration, not the
general case, and the flags say nothing about which one you are in
either way.)

**Both scales flag the same strata.** Writing \\g(u)\\ for the map from
a stratum's BLUP to its \\\pi^B_j\\, \\g\\ is strictly increasing with
\\g(0) = 0\\, so estimate, interval endpoints and zero all carry across
together. `flagged`, `direction`, `p_value` / `p_adjusted` and `pd` are
therefore identical under either `scale`, and the response-scale
interval is the *exact* image of the link-scale one rather than a
delta-method approximation – inheriting its conditionality (fixed
effects and variance components held at their point estimates), which is
why no simulation step is needed. What does change is the size and the
ranking: the same log-odds departure is worth more probability near
\\\pi = 0.5\\ than in the tail, so strata are ordered by the quantity
actually reported.

**If quantity 3 is what you need**, no argument here will give it to you
– fit a linear probability model instead:
`fit_maihda(..., family = "gaussian")` on a 0/1 outcome, which
[`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md)
signposts when it auto-detects a binary outcome. There the BLUP *is* the
risk-difference interaction by construction, at the usual costs
(predictions outside `[0, 1]`, heteroskedastic residuals).

## References

Evans, C. R., Williams, D. R., Onnela, J. P., & Subramanian, S. V.
(2018). A multilevel approach to modeling health inequalities at the
intersection of multiple social identities. *Social Science & Medicine*,
203, 64-73.

Merlo, J. (2018). Multilevel analysis of individual heterogeneity and
discriminatory accuracy (MAIHDA) within an intersectional framework.
*Social Science & Medicine*, 203, 74-80.

Evans, C. R., Leckie, G., Subramanian, S. V., Bell, A., & Merlo, J.
(2024). A tutorial for conducting intersectional multilevel analysis of
individual heterogeneity and discriminatory accuracy (MAIHDA). *SSM -
Population Health*, 26, 101664.
[doi:10.1016/j.ssmph.2024.101664](https://doi.org/10.1016/j.ssmph.2024.101664)

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

Schuirmann, D. J. (1987). A comparison of the two one-sided tests
procedure and the power approach for assessing the equivalence of
average bioavailability. *Journal of Pharmacokinetics and
Biopharmaceutics*, 15(6), 657-680.

Kruschke, J. K. (2018). Rejecting or accepting parameter values in
Bayesian estimation. *Advances in Methods and Practices in Psychological
Science*, 1(2), 270-280.

## See also

[`maihda`](https://hdbt.github.io/MAIHDA/reference/maihda.md),
[`calculate_pcv`](https://hdbt.github.io/MAIHDA/reference/calculate_pcv.md),
[`summary.maihda_model`](https://hdbt.github.io/MAIHDA/reference/summary.maihda_model.md);
and `plot(..., highlight_interactions = TRUE)` to mark the flagged
strata on the effect-decomposition / predicted / shrinkage plots.

## Examples

``` r
# \donttest{
data(maihda_health_data)
a <- maihda(BMI ~ Age + Gender + Race + (1 | Gender:Race),
            data = maihda_health_data)
maihda_interactions(a)                  # FDR-screened (default adjust = "BH")
#> ── Intersectional interactions ─────────────────────────────────────────────────
#> 4 of 10 strata flagged (95% interval; BH-adjusted conservative p-values).
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
#> 
maihda_interactions(a, adjust = "none") # uncorrected per-stratum individual view
#> ── Intersectional interactions ─────────────────────────────────────────────────
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
#> 
maihda_interactions(a, rope = 0.1)      # equivalence: |interaction| within 0.1?
#> ── Intersectional interactions ─────────────────────────────────────────────────
#> 4 of 10 strata flagged (95% interval; BH-adjusted conservative p-values).
#> Model: adjusted (two-model); interaction on the link (latent) scale.
#> Equivalence vs ROPE [-0.1, 0.1]: 4 relevant | 0 negligible | 6 inconclusive.
#> 
#>  stratum          label    n interaction     se   lower   upper  p_value
#>        2   male × Black  154     -1.2902 0.4870 -2.2446 -0.3357 0.008067
#>        9 female × Black  182      1.2902 0.4540  0.4003  2.1800 0.004487
#>        3 female × White 1044     -0.6003 0.2025 -0.9972 -0.2035 0.003025
#>        5   male × White  990      0.6003 0.2077  0.1932  1.0075 0.003855
#>  p_adjusted flagged direction decision
#>     0.02017    TRUE     below relevant
#>     0.01496    TRUE     above relevant
#>     0.01496    TRUE     below relevant
#>     0.01496    TRUE     above relevant
#> 
#> Interaction BLUPs are shrunken (partially pooled) estimates; treat flags as
#>   exploratory. See ?maihda_interactions.
#> 

# On a logistic fit the BLUP is a log-odds departure; scale = "response" reports
# Evans et al.'s (2024) pi_j - pi^A_j, the same interaction in probability points.
b <- maihda(Obese ~ Gender + Race + (1 | Gender:Race),
            data = maihda_health_data, family = "binomial")
#> Binary outcome 'Obese' recoded to 0/1: 'No' = 0 (reference), 'Yes' = 1 (modeled event). Set the factor levels (or supply a 0/1 outcome) to control which level is the event.
maihda_interactions(b, scale = "response")
#> ── Intersectional interactions ─────────────────────────────────────────────────
#> 0 of 10 strata flagged (95% interval; BH-adjusted conservative p-values).
#> Model: adjusted (two-model); interaction on the response (outcome) scale.
#> A probability difference (pi_j - pi^A_j, Evans et al. 2024): the
#>   interaction carried onto the outcome scale; flags match the
#>   log-odds scale.
#> 
#> No strata show interaction credibly different from zero at this level.
#> 
#> Interaction BLUPs are shrunken (partially pooled) estimates; treat flags as
#>   exploratory. See ?maihda_interactions.
#> 
# }
```
