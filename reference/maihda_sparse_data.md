# Sparse Intersectional Data for Bayesian MAIHDA

A simulated cross-sectional dataset in which many intersectional strata
hold only a handful of individuals, built with a known intersectional
interaction. With strata this small the split of the between-stratum
variance into additive and interaction parts is weakly identified. An
lme4 fit of the crossed-dimensions decomposition can be singular; a brms
fit of it (`engine = "brms"`) gives the split a posterior median and a
credible interval, both of which depend on the prior on the standard
deviations.

## Usage

``` r
maihda_sparse_data
```

## Format

A data frame with 240 rows and 6 variables:

- gender:

  Strata dimension (Women/Men).

- ethnicity:

  Strata dimension (White/Black/Asian).

- education:

  Strata dimension (Low/High).

- age_group:

  Strata dimension (Young/Mid/Older).

- y:

  A continuous (Gaussian) outcome with residual SD 1. The 36 stratum
  effects have a variance of 0.31 (a VPC of 0.24), of which 41% is the
  intersectional interaction.

- event:

  A binary outcome (No/Yes), ~46% "Yes". On the latent scale its stratum
  effects have a variance of 1.32 (a VPC of 0.29), again 41%
  interaction.

`attr(maihda_sparse_data, "truth")` holds, for each outcome, the
intercept and the 36 stratum effects with their additive, interaction
and between-stratum variances, the interaction share and the VPC.

## Source

Simulated; see `data-raw/maihda_sparse_data.R`.

## Details

Four dimensions form 36 intersectional strata with deliberately skewed
sizes (median 6 individuals, 12 of 36 cells below 5, two singletons).
The stratum effects were drawn with the interaction aimed at 40% of the
between-stratum variance; the 36 effects drawn carry 41% on both
outcomes (on the latent scale for `event`), by the mean-square split of
the section Realised strata of
[`maihda`](https://hdbt.github.io/MAIHDA/reference/maihda.md). The
interaction share a crossed-dimensions fit reports is a ratio of
variance parameters, and the four dimensions' variances, each estimated
from two or three levels, add up to more than the spread of the main
effects among the 36 strata: fitted to these stratum effects with 1,000
individuals in every stratum, lme4 puts the share at about 28% on `y`
and 35% on `event`, the gap between the two being that of REML and
maximum likelihood (the two sets of effects differ only in scale). In
[`vignette("bayesian_sparse_maihda")`](https://hdbt.github.io/MAIHDA/articles/bayesian_sparse_maihda.md)
both lme4 fits are singular through a dimension's variance, not the
interaction's, and put the share at about 23% on `y` and 3% on `event`.
The brms posterior medians are below 41% as well; under the priors the
vignette sets on the standard deviations, their credible intervals run
from near zero to about 41% or more.

## Note

A purely illustrative dataset. The dimension labels are arbitrary and
the interaction is constructed, not estimated from any real population –
its only purpose is to make the sparse-cell behaviour of the ML and
Bayesian estimators visible against a known answer.

## Examples

``` r
data(maihda_sparse_data)
attr(maihda_sparse_data, "truth")$gaussian$interaction_share  # 0.41
#> [1] 0.4056439

# lme4: the interaction share as a point estimate, from a fit that is singular here:
# m_lme4 <- maihda(y ~ 1 + (1 | gender:ethnicity:education:age_group),
#                  data = maihda_sparse_data, decomposition = "crossed-dimensions")
#
# brms: a posterior median and a credible interval, both resting on the prior:
# m_brms <- maihda(y ~ 1 + (1 | gender:ethnicity:education:age_group),
#                  data = maihda_sparse_data, decomposition = "crossed-dimensions",
#                  engine = "brms",
#                  prior = brms::set_prior("normal(0, 0.5)", class = "sd"))
# See vignette("bayesian_sparse_maihda").
```
