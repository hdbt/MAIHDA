# Response-scale VPC for a binomial or count MAIHDA model

Computes the variance partition coefficient on the scale of the observed
response: the probability scale for a binomial MAIHDA model, the count
scale for a Poisson or negative-binomial one.

For a **binomial** model it uses the simulation method of Goldstein,
Browne & Rasbash (2002). Stratum random effects \\u \sim N(0,
\sigma^2_u)\\ are simulated and converted to predicted probabilities \\p
= g^{-1}(\eta + u)\\ (with \\\eta\\ the fixed-part linear predictor);
the VPC is then the between-stratum variance of \\p\\ as a share of the
total (between + the binomial within-stratum variance
\\\overline{p(1-p)}\\). Unlike the latent-scale VPC (fixed level-1
variance \\\pi^2/3\\ for the logit), the response-scale VPC depends on
the overall outcome prevalence, so report it as a complement to – not a
replacement for – the latent-scale value.

For a **Poisson or negative-binomial** model with the log link it
evaluates the exact expressions of Leckie, Browne, Goldstein, Merlo &
Austin (2020). With \\\lambda = \exp(\eta + \sigma^2_u / 2)\\ the
marginal expected count, the between-stratum variance of the counts is
\\\lambda^2 (\exp(\sigma^2_u) - 1)\\ and their within-stratum variance
is \\\lambda + \alpha \lambda^2 \exp(\sigma^2_u)\\, where \\\alpha = 1 /
\theta\\ is the negative-binomial overdispersion (0 for a Poisson
model); the VPC is the between-stratum variance as a share of their sum.
It is the share of the variance of the counts themselves that lies
between strata, and it involves no approximation of the level-1
variance. The latent-scale VPC of
[`summary.maihda_model`](https://hdbt.github.io/MAIHDA/reference/summary.maihda_model.md)
refers to the log scale instead and depends on the `count_approximation`
of
[`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md),
so the two answer different questions. Report the count-scale VPC when
the partition should describe the observed counts, or when the expected
count is low and the latent-scale approximations disagree.

## Usage

``` r
maihda_vpc_response(model, n_sim = 10000, seed = NULL)
```

## Arguments

- model:

  A binomial, Poisson or negative-binomial `maihda_model` (lme4 engine)
  from
  [`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md).

- n_sim:

  Number of Monte-Carlo draws of the stratum random effect (\>= 100) for
  a binomial model. Default 10000.

- seed:

  Optional integer seed for the binomial simulation.

## Value

An object of class `maihda_vpc_response`: a list with `estimate`,
`scale = "response"`, `method` (`"simulation"` for a binomial model,
`"closed-form"` for a count model), `var_between` (the latent-scale
between-stratum variance) and `var_other` (the summed latent-scale
variance of any non-stratum random intercepts, 0 when there are none).

For a binomial model it also holds `n_sim`, `lp_fixed` (the mean
fixed-part linear predictor), `inner_method` (`"gauss-hermite"` when
non-stratum effects are integrated out, else `"none"`) and `inner_nodes`
(the number of Gauss-Hermite quadrature nodes used for that inner
integral, `NA` when there are no non-stratum effects).

For a count model it also holds `family`, the `median`, `quartiles` and
`range` of the observation-level VPCs (all equal to `estimate` when
\\\lambda\\ is the same for every observation), `lambda` (the mean
marginal expected count), `lambda_range` and `theta` (the
negative-binomial dispersion, `Inf` for a Poisson model).

## Details

**Binomial models.** The fixed part \\\eta\\ is collapsed to a single
value – the mean linear predictor \\\bar\eta\\ over the analytic sample
– before the random effect is simulated around it. The result is
therefore a VPC *evaluated at the mean covariate profile* (a
conditional-at-mean estimate), not one marginalised over the empirical
covariate distribution. For the canonical strata-only (null) model
\\\eta\\ is constant (the intercept), so the two coincide and the value
is exact. For an *adjusted* model (one with covariates) they can differ,
because the inverse link is nonlinear and \\g^{-1}(\bar\eta) \neq
\overline{g^{-1}(\eta)}\\: read the response-scale VPC from the null
model, or interpret an adjusted value as conditional on the average
covariate profile rather than as a covariate-averaged (marginal) VPC.

The method is binomial-link agnostic: it maps the simulated stratum
effects through whichever inverse link the model uses (logit, probit,
cloglog, ...), so a non-logit binomial fit is computed on its own scale
rather than rejected. Only the family is required to be binomial.

When the model carries random intercepts *beyond* the stratum (a
contextual `(1 | school)` or an explicit `(1 | site)`), the simulation
integrates over them: the reported estimate is the *stratum share*
\\Var(E\[p \mid u\_{stratum}\])\\ of the total response-scale variance,
where the total includes the variation the other random effects induce
in \\p\\ plus the binomial within-variance. Simulating the stratum
effect alone would overstate the stratum share. The inner integral over
the combined non-stratum effect is evaluated by *deterministic
Gauss-Hermite quadrature*, not an inner Monte-Carlo sample, so it
contributes no simulation error of its own and increasing `n_sim`
converges the whole estimate (the outer stratum draw is then the only
Monte-Carlo dimension). The number of quadrature nodes used is reported
in `inner_nodes`.

**Count models.** \\\lambda\\ is computed for each observation of the
analytic sample from its fixed-part linear predictor, any offset
included, so with covariates or an offset the VPC differs between
observations. The `estimate` is the mean of the observation-level VPCs;
their median, quartiles and range are returned with it, as are the mean
and range of \\\lambda\\. Prior `weights` of the fit enter these
summaries as frequencies. A null model without an offset has a single
\\\lambda\\ and a single VPC. `n_sim` and `seed` are not used.

The expressions are those of a model whose only random effect is the
stratum random intercept. A count model with a `context` or any other
additional random effect, a crossed-dimensions fit, a longitudinal fit,
and a link other than the log are refused.

## References

Goldstein, H., Browne, W., & Rasbash, J. (2002). Partitioning variation
in multilevel models. *Understanding Statistics*, 1(4), 223-231.

Leckie, G., Browne, W. J., Goldstein, H., Merlo, J., & Austin, P. C.
(2020). Partitioning variation in multilevel models for count data.
*Psychological Methods*, 25(6), 787-801.
[doi:10.1037/met0000265](https://doi.org/10.1037/met0000265)

## See also

[`maihda_discriminatory_accuracy`](https://hdbt.github.io/MAIHDA/reference/maihda_discriminatory_accuracy.md),
[`summary.maihda_model`](https://hdbt.github.io/MAIHDA/reference/summary.maihda_model.md)

## Examples

``` r
# \donttest{
strata <- make_strata(maihda_health_data, vars = c("Gender", "Race"))
d <- maihda_health_data
d$stratum <- strata$data$stratum
m <- fit_maihda(Obese ~ (1 | stratum), data = d, family = "binomial")
#> Binary outcome 'Obese' recoded to 0/1: 'No' = 0 (reference), 'Yes' = 1 (modeled event). Set the factor levels (or supply a 0/1 outcome) to control which level is the event.
maihda_vpc_response(m, seed = 1)
#> Response-scale VPC (simulation method)
#>   VPC: 0.0370
#>   10000 simulated stratum effects; between-stratum variance 0.1700 (latent scale).
#> 

# A count outcome: the VPC on the count scale, in closed form
set.seed(1)
u <- rnorm(max(d$stratum), sd = 0.4)
d$visits <- rpois(nrow(d), exp(0.5 + u[d$stratum]))
m_count <- fit_maihda(visits ~ (1 | stratum), data = d, family = "poisson")
maihda_vpc_response(m_count)
#> Count-scale VPC (closed form, Leckie et al. 2020)
#>   VPC: 0.1463
#>   Marginal expected count 1.8; between-stratum variance 0.0911 (log scale).
#> 
# }
```
