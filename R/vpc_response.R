# Response-scale (observed / probability-scale) VPC for binomial MAIHDA models.
#
# The package's default VPC for a logistic model is on the latent scale, with a
# fixed level-1 variance (pi^2/3 for the logit link). That is scale-free but not on
# an interpretable probability metric. This file adds the response-scale VPC via the
# simulation method of Goldstein, Browne & Rasbash (2002): simulate stratum random
# effects, map them through the inverse link to probabilities, and express the
# between-stratum variance of those probabilities as a share of the total. It is an
# interpretable complement to the latent-scale VPC, not a replacement (it depends on
# the overall prevalence).

#' Response-scale VPC for a binomial MAIHDA model
#'
#' @description
#' Computes the variance partition coefficient on the response (probability) scale
#' for a binomial MAIHDA model, using the simulation method of Goldstein, Browne &
#' Rasbash (2002). Stratum random effects \eqn{u \sim N(0, \sigma^2_u)} are
#' simulated and converted to predicted probabilities \eqn{p = g^{-1}(\eta + u)}
#' (with \eqn{\eta} the fixed-part linear predictor); the VPC is then the
#' between-stratum variance of \eqn{p} as a share of the total
#' (between + the binomial within-stratum variance \eqn{\overline{p(1-p)}}).
#'
#' Unlike the latent-scale VPC (fixed level-1 variance \eqn{\pi^2/3} for the logit),
#' the response-scale VPC depends on the overall outcome prevalence, so report it as
#' a complement to -- not a replacement for -- the latent-scale value.
#'
#' @details
#' The fixed part \eqn{\eta} is collapsed to a single value -- the mean linear
#' predictor \eqn{\bar\eta} over the analytic sample -- before the random effect is
#' simulated around it. The result is therefore a VPC \emph{evaluated at the mean
#' covariate profile} (a conditional-at-mean estimate), not one marginalised over the
#' empirical covariate distribution. For the canonical strata-only (null) model
#' \eqn{\eta} is constant (the intercept), so the two coincide and the value is
#' exact. For an \emph{adjusted} model (one with covariates) they can differ, because
#' the inverse link is nonlinear and \eqn{g^{-1}(\bar\eta) \neq \overline{g^{-1}(\eta)}}:
#' read the response-scale VPC from the null model, or interpret an adjusted value as
#' conditional on the average covariate profile rather than as a covariate-averaged
#' (marginal) VPC.
#'
#' @param model A binomial \code{maihda_model} (lme4 engine) from
#'   \code{\link{fit_maihda}}.
#' @param n_sim Number of Monte-Carlo draws of the stratum random effect (>= 100).
#'   Default 10000.
#' @param seed Optional integer seed for reproducibility.
#'
#' @return An object of class \code{maihda_vpc_response}: a list with
#'   \code{estimate}, \code{scale = "response"}, \code{method = "simulation"},
#'   \code{n_sim}, \code{var_between} (the latent-scale between-stratum variance) and
#'   \code{lp_fixed} (the mean fixed-part linear predictor).
#'
#' @references
#' Goldstein, H., Browne, W., & Rasbash, J. (2002). Partitioning variation in
#' multilevel models. \emph{Understanding Statistics}, 1(4), 223-231.
#'
#' @seealso \code{\link{maihda_discriminatory_accuracy}}, \code{\link{summary.maihda_model}}
#'
#' @examples
#' \dontrun{
#' strata <- make_strata(maihda_health_data, vars = c("Gender", "Race"))
#' d <- maihda_health_data
#' d$stratum <- strata$data$stratum
#' m <- fit_maihda(Obese ~ (1 | stratum), data = d, family = "binomial")
#' maihda_vpc_response(m, seed = 1)
#' }
#'
#' @export
maihda_vpc_response <- function(model, n_sim = 10000, seed = NULL) {
  if (!inherits(model, "maihda_model")) {
    stop("'model' must be a maihda_model object from fit_maihda().", call. = FALSE)
  }
  if (!identical(model$engine, "lme4")) {
    stop("maihda_vpc_response() is currently implemented for the lme4 engine.",
         call. = FALSE)
  }
  fam <- maihda_model_family_name(model)
  if (!identical(fam, "binomial")) {
    stop("The response-scale VPC is only defined for binomial (logistic) MAIHDA ",
         "models; this model uses family = '", fam, "'.", call. = FALSE)
  }
  if (!is.numeric(n_sim) || length(n_sim) != 1 || !is.finite(n_sim) || n_sim < 100) {
    stop("'n_sim' must be a single number >= 100.", call. = FALSE)
  }

  var_between <- tryCatch(extract_between_variance(model), error = function(e) NA_real_)
  if (!is.numeric(var_between) || length(var_between) != 1 ||
      !is.finite(var_between) || var_between < 0) {
    return(structure(
      list(estimate = NA_real_, scale = "response", method = "simulation",
           n_sim = n_sim, var_between = var_between, lp_fixed = NA_real_),
      class = "maihda_vpc_response"))
  }

  fitted_model <- model$model
  linkinv <- stats::family(fitted_model)$linkinv
  # Fixed-part linear predictor (random effects excluded), collapsed to its sample
  # mean. For the canonical null / strata-only model this is exactly the intercept,
  # so the VPC below is exact; with covariates this is the mean fixed linear
  # predictor, so the VPC is a conditional-at-mean estimate (evaluated at the average
  # covariate profile) rather than one integrated over the covariate distribution --
  # see the @details section of the function documentation.
  lp_fixed <- mean(stats::predict(fitted_model, re.form = NA, type = "link"), na.rm = TRUE)

  if (!is.null(seed)) {
    set.seed(seed)
  }
  u <- stats::rnorm(n_sim, mean = 0, sd = sqrt(var_between))
  p <- linkinv(lp_fixed + u)
  v_between <- stats::var(p)
  v_within <- mean(p * (1 - p))
  vpc <- v_between / (v_between + v_within)

  structure(
    list(estimate = vpc, scale = "response", method = "simulation",
         n_sim = n_sim, var_between = var_between, lp_fixed = lp_fixed),
    class = "maihda_vpc_response"
  )
}

#' @export
print.maihda_vpc_response <- function(x, ...) {
  cat("Response-scale VPC (simulation method)\n")
  cat(sprintf("  VPC: %s\n",
              if (is.finite(x$estimate)) sprintf("%.4f", x$estimate) else "NA"))
  cat(sprintf("  %d simulated stratum effects; between-stratum variance %.4f (latent scale).\n",
              x$n_sim, x$var_between))
  invisible(x)
}
