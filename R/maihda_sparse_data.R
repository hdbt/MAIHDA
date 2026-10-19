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
#' effects were generated with the interaction at 40\% of the between-stratum
#' variance on both outcomes, on the latent scale for \code{event}. In
#' \code{vignette("bayesian_sparse_maihda")} neither engine recovers it. Both lme4
#' fits are singular through a dimension's variance, not
#' the interaction's, and put the interaction share at about 23\% on \code{y} and
#' 3\% on \code{event}. The brms posterior medians are below 40\% as well; under
#' the priors the vignette sets on the standard deviations, their credible
#' intervals run from near zero to past 40\%.
#'
#' @format A data frame with 240 rows and 6 variables:
#' \describe{
#'   \item{gender}{Strata dimension (Women/Men).}
#'   \item{ethnicity}{Strata dimension (White/Black/Asian).}
#'   \item{education}{Strata dimension (Low/High).}
#'   \item{age_group}{Strata dimension (Young/Mid/Older).}
#'   \item{y}{A continuous (Gaussian) outcome. True between-stratum VPC 0.26, of
#'     which 40\% is the intersectional interaction.}
#'   \item{event}{A binary outcome (No/Yes), ~46\% "Yes". Its latent-scale
#'     between-stratum VPC is 0.31, again 40\% interaction.}
#' }
#' The exact generative truth is also attached as
#' \code{attr(maihda_sparse_data, "truth")} (additive/interaction variances, shares,
#' and VPCs for each outcome).
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
#' attr(maihda_sparse_data, "truth")$gaussian$interaction_share  # 0.40
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
