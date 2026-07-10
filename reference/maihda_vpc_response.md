# Response-scale VPC for a binomial MAIHDA model

Computes the variance partition coefficient on the response
(probability) scale for a binomial MAIHDA model, using the simulation
method of Goldstein, Browne & Rasbash (2002). Stratum random effects \\u
\sim N(0, \sigma^2_u)\\ are simulated and converted to predicted
probabilities \\p = g^{-1}(\eta + u)\\ (with \\\eta\\ the fixed-part
linear predictor); the VPC is then the between-stratum variance of \\p\\
as a share of the total (between + the binomial within-stratum variance
\\\overline{p(1-p)}\\).

Unlike the latent-scale VPC (fixed level-1 variance \\\pi^2/3\\ for the
logit), the response-scale VPC depends on the overall outcome
prevalence, so report it as a complement to – not a replacement for –
the latent-scale value.

## Usage

``` r
maihda_vpc_response(model, n_sim = 10000, seed = NULL)
```

## Arguments

- model:

  A binomial `maihda_model` (lme4 engine) from
  [`fit_maihda`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md).

- n_sim:

  Number of Monte-Carlo draws of the stratum random effect (\>= 100).
  Default 10000.

- seed:

  Optional integer seed for reproducibility.

## Value

An object of class `maihda_vpc_response`: a list with `estimate`,
`scale = "response"`, `method = "simulation"`, `n_sim`, `var_between`
(the latent-scale between-stratum variance), `var_other` (the summed
latent-scale variance of any non-stratum random intercepts, 0 when there
are none) and `lp_fixed` (the mean fixed-part linear predictor).

## Details

The fixed part \\\eta\\ is collapsed to a single value – the mean linear
predictor \\\bar\eta\\ over the analytic sample – before the random
effect is simulated around it. The result is therefore a VPC *evaluated
at the mean covariate profile* (a conditional-at-mean estimate), not one
marginalised over the empirical covariate distribution. For the
canonical strata-only (null) model \\\eta\\ is constant (the intercept),
so the two coincide and the value is exact. For an *adjusted* model (one
with covariates) they can differ, because the inverse link is nonlinear
and \\g^{-1}(\bar\eta) \neq \overline{g^{-1}(\eta)}\\: read the
response-scale VPC from the null model, or interpret an adjusted value
as conditional on the average covariate profile rather than as a
covariate-averaged (marginal) VPC.

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
effect alone would overstate the stratum share.

## References

Goldstein, H., Browne, W., & Rasbash, J. (2002). Partitioning variation
in multilevel models. *Understanding Statistics*, 1(4), 223-231.

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
# }
```
