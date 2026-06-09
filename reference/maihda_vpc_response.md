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
(the latent-scale between-stratum variance) and `lp_fixed` (the mean
fixed-part linear predictor).

## References

Goldstein, H., Browne, W., & Rasbash, J. (2002). Partitioning variation
in multilevel models. *Understanding Statistics*, 1(4), 223-231.

## See also

[`maihda_discriminatory_accuracy`](https://hdbt.github.io/MAIHDA/reference/maihda_discriminatory_accuracy.md),
[`summary.maihda_model`](https://hdbt.github.io/MAIHDA/reference/summary.maihda_model.md)

## Examples

``` r
if (FALSE) { # \dontrun{
strata <- make_strata(maihda_health_data, vars = c("Gender", "Race"))
d <- maihda_health_data
d$stratum <- strata$data$stratum
m <- fit_maihda(Obese ~ (1 | stratum), data = d, family = "binomial")
maihda_vpc_response(m, seed = 1)
} # }
```
