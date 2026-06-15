#' Sparse Intersectional Data for Bayesian MAIHDA
#'
#' A simulated cross-sectional dataset built to showcase **Bayesian (brms) MAIHDA
#' for sparse intersections** -- the regime where many intersectional strata each
#' hold only a handful of individuals. There the maximum-likelihood (lme4) estimate
#' of the *interaction* between-stratum variance collapses to a singular fit with no
#' uncertainty, so the additive-vs-interaction split is both unstable and falsely
#' precise; weakly-informative priors (`engine = "brms"`) regularise the variance
#' off the boundary and return a calibrated credible interval.
#'
#' The data carry a **known, non-trivial interaction** so the vignette can claim
#' *recovery* rather than merely report numbers: 4 dimensions form 36 intersectional
#' strata with deliberately skewed sizes (median 6 individuals, 12 of 36 cells below
#' 5, two singletons), and the true interaction accounts for **40% of the
#' between-stratum variance** on both outcomes. On the binary outcome a genuine 40%
#' interaction is read by lme4 as roughly 3% -- a spurious "fully additive" result
#' that is purely a small-cell artifact.
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
#' # ML over-shrinks the interaction under sparse cells (a singular fit):
#' # m_lme4 <- maihda(y ~ 1 + (1 | gender:ethnicity:education:age_group),
#' #                  data = maihda_sparse_data, decomposition = "crossed-dimensions")
#' #
#' # Weakly-informative priors regularise it and report honest uncertainty:
#' # m_brms <- maihda(y ~ 1 + (1 | gender:ethnicity:education:age_group),
#' #                  data = maihda_sparse_data, decomposition = "crossed-dimensions",
#' #                  engine = "brms",
#' #                  prior = brms::set_prior("normal(0, 0.5)", class = "sd"))
#' # See vignette("bayesian_sparse_maihda").
"maihda_sparse_data"
