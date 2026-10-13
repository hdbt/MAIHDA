# Longitudinal MAIHDA proportional change in variance (PCV)

Compares the stratum-level random-effect covariance block of the null
growth model with that of the adjusted model (null + dimension main
effects + their `dim:time` interactions). Reports the PCV in the
baseline variance and in the instantaneous-slope variance at baseline –
the additive-vs-multiplicative split of the intersectional trajectory
inequality (Bell, Evans, Holman & Leckie 2024) – and the time-specific
PCV over the supplied times. Both are evaluated at the observed baseline
time (`ref_time`), so they are invariant to how the time axis is coded
(for linear growth the slope variance is the same at every time, so this
reduces to the slope-variance cell).

## Usage

``` r
maihda_longitudinal_pcv(
  null_model,
  adjusted_model,
  times = NULL,
  estimation = c("fitted", "ML")
)
```

## Arguments

- null_model, adjusted_model:

  Longitudinal `maihda_model`s from a
  `maihda(decomposition = "longitudinal")` pair.

- times:

  Optional numeric times for the time-specific PCV; defaults to the null
  model's reporting grid.

- estimation:

  Variance-estimation basis, `"fitted"` (default) or `"ML"`; see
  [`calculate_pcv`](https://hdbt.github.io/MAIHDA/reference/calculate_pcv.md).

## Value

An object of class `maihda_long_pcv`. `variance_summary` names what
`var_baseline_*`, `var_slope_*` and the variances in `pcv_t` are
(`"estimate"` for lme4, `"posterior median"` for brms), and
`Sigma_summary` what the two stored blocks are (`"estimate"` or
`"posterior mean"`).

## Details

As in
[`calculate_pcv`](https://hdbt.github.io/MAIHDA/reference/calculate_pcv.md),
the `estimation` argument selects the variance-estimation basis. With
`estimation = "fitted"` (the default) each fit's own REML covariance
block is used, matching the single-model summaries. With
`estimation = "ML"`, REML `lmer` growth fits are refitted with maximum
likelihood ([`refitML`](https://rdrr.io/pkg/lme4/man/refitML.html))
before the comparison, which overstates both PCVs when strata are few.
The stored models (and the single-model summaries computed from them,
e.g. the time-varying VPC) always keep their REML fit; `ml_refit` on the
result records whether an ML refit fully applied, and `estimation_used`
records the basis ACTUALLY used – `"fitted"`, `"ML"`, or `"mixed"` when
an ML refit was requested but a model kept its REML fit (a boundary skip
or a failed refit), so a mixed REML/ML comparison is not mistaken for a
clean one. glmer (GLMM) and brms fits are already on the ML / posterior
scale, so the choice does not affect them.

For a `brms` pair every variance compared is the posterior median of its
per-draw value: the between-stratum variance \\a(t)' \Sigma a(t)\\ at
the baseline and at each reporting time, and the slope variance. The
stored blocks `Sigma_stratum_null` and `Sigma_stratum_adjusted` are the
posterior means [`summary()`](https://rdrr.io/r/base/summary.html)
tabulates, so those variances are not quadratic forms of them. The
posterior of a between-stratum variance is right-skewed when strata are
few, and more so in the adjusted model, so the PCV of a `brms` pair is
understated on average with few strata – by more when the additive share
is moderate, when the strata are small, or, for the slope, when the
strata differ little in slope next to the precision of a stratum's slope
– and can be negative where the true value is positive. The `lme4`
(REML) decomposition is less biased then; neither is precise with so few
strata.
