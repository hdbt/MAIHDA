#' Sparse Intersectional Data for Bayesian MAIHDA
#'
#' A simulated cross-sectional dataset in which many intersectional strata hold
#' only a handful of individuals, built with a known intersectional interaction.
#' With strata this small the split of the between-stratum variance into additive
#' and interaction parts is weakly identified. An lme4 fit of the crossed-dimensions
#' decomposition can be singular; a brms fit of it (\code{engine = "brms"}) gives
#' the split a posterior median and a credible interval, both of which depend on
#' the prior on the standard deviations.
#'
#' Four dimensions form 36 intersectional strata with deliberately skewed sizes
#' (median 6 individuals, 12 of 36 cells below 5, two singletons). The stratum
#' effects were drawn with the interaction aimed at 40\% of the between-stratum
#' variance; the 36 effects drawn carry 41\% on both outcomes (on the latent scale
#' for \code{event}), by the mean-square split of the section Realised strata of
#' \code{\link{maihda}}. The interaction share a crossed-dimensions fit reports is a
#' ratio of variance parameters, and the four dimensions' variances, each estimated
#' from two or three levels, add up to more than the spread of the main effects
#' among the 36 strata: fitted to these stratum effects with 1,000 individuals in
#' every stratum, lme4 puts the share at about 28\% on \code{y} and 35\% on
#' \code{event}, the gap between the two being that of REML and maximum likelihood
#' (the two sets of effects differ only in scale). In
#' \code{vignette("bayesian_sparse_maihda")} both lme4 fits are singular through a
#' dimension's variance, not the interaction's, and put the share at about 23\% on
#' \code{y} and 3\% on \code{event}. The brms posterior medians are below 41\% as
#' well; under the priors the vignette sets on the standard deviations, their
#' credible intervals run from near zero to about 41\% or more.
#'
#' @format A data frame with 240 rows and 6 variables:
#' \describe{
#'   \item{gender}{Strata dimension (Women/Men).}
#'   \item{ethnicity}{Strata dimension (White/Black/Asian).}
#'   \item{education}{Strata dimension (Low/High).}
#'   \item{age_group}{Strata dimension (Young/Mid/Older).}
#'   \item{y}{A continuous (Gaussian) outcome with residual SD 1. The 36 stratum
#'     effects have a variance of 0.31 (a VPC of 0.24), of which 41\% is the
#'     intersectional interaction.}
#'   \item{event}{A binary outcome (No/Yes), ~46\% "Yes". On the latent scale its
#'     stratum effects have a variance of 1.32 (a VPC of 0.29), again 41\%
#'     interaction.}
#' }
#' \code{attr(maihda_sparse_data, "truth")} holds, for each outcome, the intercept
#' and the 36 stratum effects with their additive, interaction and between-stratum
#' variances, the interaction share and the VPC.
#'
#' @source Simulated; see \code{data-raw/maihda_sparse_data.R}.
#'
#' @note A purely illustrative dataset. The dimension labels are arbitrary and the
#'   interaction is constructed, not estimated from any real population -- its only
#'   purpose is to make the sparse-cell behaviour of the ML and Bayesian estimators
#'   visible against a known answer.
#'
#' @examples
#' data(maihda_sparse_data)
#' attr(maihda_sparse_data, "truth")$gaussian$interaction_share  # 0.41
#'
#' # lme4: the interaction share as a point estimate, from a fit that is singular here:
#' # m_lme4 <- maihda(y ~ 1 + (1 | gender:ethnicity:education:age_group),
#' #                  data = maihda_sparse_data, decomposition = "crossed-dimensions")
#' #
#' # brms: a posterior median and a credible interval, both resting on the prior:
#' # m_brms <- maihda(y ~ 1 + (1 | gender:ethnicity:education:age_group),
#' #                  data = maihda_sparse_data, decomposition = "crossed-dimensions",
#' #                  engine = "brms",
#' #                  prior = brms::set_prior("normal(0, 0.5)", class = "sd"))
#' # See vignette("bayesian_sparse_maihda").
"maihda_sparse_data"
