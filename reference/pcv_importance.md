# Order-Invariant PCV Attribution Across Predictors (Shapley / Dominance)

Apportions the between-stratum variance reduction (the PCV) among a set
of predictors. Where
[`stepwise_pcv`](https://hdbt.github.io/MAIHDA/reference/stepwise_pcv.md)
adds the variables one at a time – so each variable's contribution
depends on its entry order – `pcv_importance()` treats the PCV as a
value function over variable subsets and attributes it *fairly*: the
flagship `"shapley"` method averages each variable's marginal PCV over
every possible entry order, `"dominance"` adds Budescu's pairwise
dominance detail (its general dominance weights coincide with the
Shapley values), and `"sequential"` reproduces the order-dependent path
for continuity. All methods satisfy the **efficiency** identity: the
contributions sum exactly to the full model's total PCV.

Two attribution targets are useful in a MAIHDA. Passing the *stratum
dimensions* (e.g. `c("gender", "race", "education")`) splits the
additive share – the canonical null-to-adjusted PCV – fairly across the
dimensions: which social dimension drives the additive between-stratum
inequality. Passing *individual-level covariates* gives an
order-invariant counterpart to
[`stepwise_pcv()`](https://hdbt.github.io/MAIHDA/reference/stepwise_pcv.md):
how much each covariate explains of the between-stratum variance
(subject to the latent-scale caveat below for non-Gaussian families).

## Usage

``` r
pcv_importance(
  data,
  outcome,
  vars,
  method = c("shapley", "sequential", "dominance"),
  approx = c("exact", "montecarlo"),
  n_perm = 2000,
  engine = "lme4",
  family = "gaussian",
  context = NULL,
  sampling_weights = NULL,
  bootstrap = FALSE,
  n_boot = 1000,
  conf_level = 0.95
)
```

## Arguments

- data:

  Data frame with observations. Ensure
  [`make_strata`](https://hdbt.github.io/MAIHDA/reference/make_strata.md)
  was run first so the `stratum` variable exists.

- outcome:

  Character string; the dependent variable.

- vars:

  Character vector of predictors (stratum dimensions and/or covariates)
  among which the PCV is apportioned. Order does not affect the
  `"shapley"` and `"dominance"` results; it defines the path for
  `"sequential"`.

- method:

  Attribution method: `"shapley"` (default; order-invariant Shapley
  values), `"sequential"` (the order-dependent one-at-a-time path, kept
  for continuity with
  [`stepwise_pcv`](https://hdbt.github.io/MAIHDA/reference/stepwise_pcv.md)),
  or `"dominance"` (Budescu dominance analysis: general dominance –
  equal to the Shapley values – plus the conditional and complete
  dominance detail).

- approx:

  For `method = "shapley"` only: `"exact"` fits all \\2^k - 1\\
  non-empty variable subsets (plus the null); `"montecarlo"` samples
  `n_perm` random entry orders instead. When `approx` is not supplied it
  defaults to `"exact"` for up to 10 variables and switches to
  `"montecarlo"` (with a message) beyond. `"sequential"` needs only the
  \\k\\ path models and ignores `approx`; `"dominance"` requires every
  subset, so it errors if `approx = "montecarlo"` is requested.

- n_perm:

  Number of random permutations for `approx = "montecarlo"`.
  Default 2000. Distinct subset fits are cached, so the number of *model
  fits* is bounded by `min(2^k, n_perm * (k - 1) + k + 1)`; set a seed
  ([`set.seed()`](https://rdrr.io/r/base/Random.html)) for reproducible
  sampling.

- engine:

  Modeling engine ("lme4", "brms", "wemix", or "ordinal"); default
  "lme4". Resolved exactly as in
  [`stepwise_pcv`](https://hdbt.github.io/MAIHDA/reference/stepwise_pcv.md):
  switches to "wemix" when `sampling_weights` is supplied and to
  "ordinal" for an ordinal family or ordered-factor outcome. Exact
  attribution with `engine = "brms"` refits a Stan model for *every*
  subset and is strongly discouraged beyond a handful of variables (see
  Details).

- family:

  Error distribution and link function. Default "gaussian"; a binary or
  ordered-factor outcome is auto-detected when `family` is left
  unspecified, mirroring
  [`stepwise_pcv`](https://hdbt.github.io/MAIHDA/reference/stepwise_pcv.md).

- context:

  Optional higher-level context column(s) (e.g. `"school"`), forwarded
  to every subset fit so each model carries the crossed `(1 | context)`
  intercept alongside the stratum effect; the attribution is then of the
  between-stratum PCV **net of** the context, exactly as in
  [`stepwise_pcv`](https://hdbt.github.io/MAIHDA/reference/stepwise_pcv.md).
  `lme4`/`brms` engines only.

- sampling_weights:

  Optional name of a sampling-weight column for design-weighted fits;
  see
  [`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md).
  The weight column joins the complete-case filter so every subset fit
  uses the same analytic sample.

- bootstrap:

  Logical; compute parametric-bootstrap confidence intervals for each
  contribution by refitting *every* subset model on responses simulated
  from the full model. **lme4 engine and exact attribution only**; the
  cost is `n_boot` times the number of subset models (up to
  `n_boot * 2^k` refits). Default FALSE.

- n_boot:

  Number of bootstrap draws if `bootstrap = TRUE`. Default 1000.

- conf_level:

  Confidence level for bootstrap intervals. Default 0.95.

## Value

An object of class `maihda_pcv_importance`: a list with

- importance:

  Data frame with one row per variable (in the order of `vars`):
  `Contribution` (the variable's share of the total PCV, on the PCV
  scale), `Share` (`Contribution / total_pcv`; `NA` when the total PCV
  is zero), plus `MC_SE` (Monte-Carlo standard errors,
  `approx = "montecarlo"` only) and `CI_lower`/`CI_upper`
  (`bootstrap = TRUE` only).

- total_pcv:

  The full-model PCV, \\(V_0 - V\_{full}) / V_0\\; the contributions sum
  to this value (efficiency).

- null_variance, full_variance:

  Between-stratum variance of the null and the full (all-`vars`) model.

- subsets:

  Data frame of every subset model fit: the variables in the subset, its
  size, between-stratum variance, and PCV \\v(S)\\.

- conditional, complete_dominance:

  `method = "dominance"` only: the conditional dominance matrix
  (variables x adjustment-set size) and the pairwise complete-dominance
  matrix.

- method, approx, n_perm, n_fits, n_obs, engine, family, context,
  bootstrap, conf_level, n_boot_ok:

  Metadata; `n_fits` counts the distinct models fit (including the
  null).

## Details

**Value function and efficiency.** Write \\V_0\\ for the null model's
between-stratum variance and \\V(S)\\ for the between-stratum variance
after adding the variable subset \\S\\ as fixed effects (each REML
`lmer` fit is refit with ML first, exactly as in
[`calculate_pcv`](https://hdbt.github.io/MAIHDA/reference/calculate_pcv.md)
and
[`stepwise_pcv`](https://hdbt.github.io/MAIHDA/reference/stepwise_pcv.md),
so all attributions live on the same scale as the rest of the package).
The value function is \\v(S) = (V_0 - V(S)) / V_0\\ – the total PCV of
the model that adds \\S\\ – and the Shapley contribution of variable
\\i\\ averages its marginal PCV \\v(S \cup \\i\\) - v(S)\\ over all
subsets with the usual Shapley weights. Because every permutation's
marginals telescope to \\v(N)\\, the contributions of *every* method
here (including the Monte-Carlo approximation, for any `n_perm`) sum
exactly to the full-model total PCV. This is the multilevel-PCV analogue
of the LMG / Shorrocks-Shapley decomposition of \\R^2\\ (Groemping 2006;
Shorrocks 2013).

**Sequential method vs.
[`stepwise_pcv()`](https://hdbt.github.io/MAIHDA/reference/stepwise_pcv.md).**
The `"sequential"` contributions are the *increments in total PCV* along
the entry order, \\v(\\x_1..x_i\\) - v(\\x_1..x\_{i-1}\\)\\ – i.e.
`diff(c(0, Total_PCV))` of the
[`stepwise_pcv`](https://hdbt.github.io/MAIHDA/reference/stepwise_pcv.md)
table, so they sum to the same total. They are *not* the `Step_PCV`
column, which normalises each step by the *previous* step's variance
rather than by \\V_0\\.

**Exact vs. Monte-Carlo cost.** The exact attribution fits \\2^k - 1\\
subset models plus the null, with every fit cached and reused across all
marginal differences: 255 fits at \\k = 8\\, 1023 at \\k = 10\\. That is
feasible for `lme4` but quickly infeasible beyond, and *each* of those
fits is a separate Stan run under `engine = "brms"` – exact attribution
on brms is therefore strongly discouraged except for very small \\k\\.
The Monte-Carlo route samples `n_perm` entry orders (unbiased for the
Shapley values, since a uniform random permutation reproduces the
Shapley weights), reports a per-variable Monte-Carlo standard error, and
warns when the largest `MC_SE` exceeds 0.01 on the PCV scale (increase
`n_perm`).

**Latent-scale families and rescaling.** For binomial/ordinal (and, in
attenuated form, count) families, adding a predictor that varies
*within* strata rescales the latent metric, so its marginal PCV mixes
explained variance with rescaling (see the latent-scale note in
[`calculate_pcv`](https://hdbt.github.io/MAIHDA/reference/calculate_pcv.md));
Shapley values of the PCV inherit this. Attribution among the stratum
*dimensions* – constant within each stratum – is largely unaffected.

**Suppression and negative contributions.** The PCV can be negative, so
a contribution (and its `Share`) can be negative or exceed 100%.
Efficiency still holds; contributions are reported *signed* and no
non-negative normalisation is applied. A negative Shapley contribution
flags a suppressor-style variable whose inclusion tends to *raise* the
between-stratum variance.

**Bootstrap.** With `bootstrap = TRUE` the whole attribution is
bootstrapped: responses are simulated from the full model and every
subset model is refit per draw
([`lme4::refit`](https://rdrr.io/pkg/lme4/man/refit.html)), giving
percentile intervals per contribution – `n_boot` times `n_fits` refits,
so gate it by cost. Available for `engine = "lme4"` with exact
attribution (any `method`); the Monte-Carlo approximation already
carries its own sampling error and is not combined with the bootstrap.
As in
[`calculate_pcv`](https://hdbt.github.io/MAIHDA/reference/calculate_pcv.md),
`refit()` holds a `glmer.nb` dispersion parameter fixed at its original
estimate.

## Reproducibility

`approx = "montecarlo"` draws random permutations with the session RNG;
call [`set.seed()`](https://rdrr.io/r/base/Random.html) first for
reproducible results. The exact methods are deterministic.

## References

Budescu, D. V. (1993). Dominance analysis: a new approach to the problem
of relative importance of predictors in multiple regression.
*Psychological Bulletin*, 114(3), 542-551.

Groemping, U. (2006). Relative importance for linear regression in R:
the package relaimpo. *Journal of Statistical Software*, 17(1), 1-27.

Shorrocks, A. F. (2013). Decomposition procedures for distributional
analysis: a unified framework based on the Shapley value. *Journal of
Economic Inequality*, 11, 99-126.

## See also

[`stepwise_pcv`](https://hdbt.github.io/MAIHDA/reference/stepwise_pcv.md)
for the sequential path with the discriminatory-accuracy trajectory,
[`calculate_pcv`](https://hdbt.github.io/MAIHDA/reference/calculate_pcv.md)
for the two-model PCV and the latent-scale caveat, and
[`maihda`](https://hdbt.github.io/MAIHDA/reference/maihda.md) for the
canonical null-vs-adjusted decomposition.

## Examples

``` r
# \donttest{
strata <- make_strata(maihda_sim_data, c("gender", "race"))

# Order-invariant split of the PCV across the two dimensions and age
imp <- pcv_importance(strata$data, "health_outcome", c("gender", "race", "age"))
print(imp)
#> PCV Attribution Across Predictors
#> =================================
#> 
#> Method:  Shapley values (exact)
#> Outcome: health_outcome   Engine: lme4 (gaussian(identity))
#> Analytic sample: 500 observations; 8 models fit (incl. null).
#> 
#> Between-stratum variance: null 23.241439 -> full model 0.000000
#> Total PCV (null -> all variables): 1.0000
#> 
#>  Variable Contribution  Share
#>    gender       0.0143  +1.4%
#>      race       0.9814 +98.1%
#>       age       0.0043  +0.4%
#>     Total       1.0000 100.0%
#> 
#> Contributions are shares of the null model's between-stratum variance
#> explained (the PCV scale) and sum to the full-model Total PCV (efficiency).
#> 
plot(imp)


# Dominance analysis: adds conditional / complete dominance detail
pcv_importance(strata$data, "health_outcome", c("gender", "race", "age"),
               method = "dominance")
#> PCV Attribution Across Predictors
#> =================================
#> 
#> Method:  Dominance analysis (general dominance = Shapley)
#> Outcome: health_outcome   Engine: lme4 (gaussian(identity))
#> Analytic sample: 500 observations; 8 models fit (incl. null).
#> 
#> Between-stratum variance: null 23.241439 -> full model 0.000000
#> Total PCV (null -> all variables): 1.0000
#> 
#>  Variable Contribution  Share
#>    gender       0.0143  +1.4%
#>      race       0.9814 +98.1%
#>       age       0.0043  +0.4%
#>     Total       1.0000 100.0%
#> 
#> Contributions are shares of the null model's between-stratum variance
#> explained (the PCV scale) and sum to the full-model Total PCV (efficiency).
#> 
#> Conditional dominance (average marginal PCV by adjustment-set size):
#>        size_0  size_1 size_2
#> gender 0.0146  0.0072 0.0210
#> race   0.9927  0.9744 0.9773
#> age    0.0157 -0.0028 0.0000
#> Complete dominance: race > gender, race > age

# Monte-Carlo approximation (set a seed for reproducibility)
set.seed(42)
pcv_importance(strata$data, "health_outcome", c("gender", "race", "age"),
               approx = "montecarlo", n_perm = 500)
#> PCV Attribution Across Predictors
#> =================================
#> 
#> Method:  Shapley values (Monte-Carlo, 500 permutations)
#> Outcome: health_outcome   Engine: lme4 (gaussian(identity))
#> Analytic sample: 500 observations; 8 models fit (incl. null).
#> 
#> Between-stratum variance: null 23.241439 -> full model 0.000000
#> Total PCV (null -> all variables): 1.0000
#> 
#>  Variable Contribution  Share  MC_SE
#>    gender       0.0141  +1.4% 0.0003
#>      race       0.9820 +98.2% 0.0005
#>       age       0.0039  +0.4% 0.0005
#>     Total       1.0000 100.0%       
#> 
#> Contributions are shares of the null model's between-stratum variance
#> explained (the PCV scale) and sum to the full-model Total PCV (efficiency).
#> 
# }
```
