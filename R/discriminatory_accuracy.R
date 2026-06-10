# Discriminatory accuracy for binary MAIHDA models.
#
# The "DA" in MAIHDA is discriminatory accuracy: how well the intersectional
# strata separate individuals who do and do not have a binary outcome. The VPC
# summarises *variation* between strata; discriminatory accuracy summarises
# *prediction* at the individual level, and the two can diverge sharply (a high
# between-stratum VPC can still translate into only modest individual-level AUC).
# This file promotes the helpers previously sketched in the "binary_outcomes"
# vignette to first-class, tested, exported functions.

#' Area under the ROC curve (C-statistic), rank-based
#'
#' @description
#' Computes the AUC / C-statistic as the Mann-Whitney U statistic: the
#' probability that a randomly chosen case (\code{y == 1}) is assigned a higher
#' predicted value than a randomly chosen non-case (\code{y == 0}), with ties
#' counting as one half. This needs no external package. An AUC of 0.5 is chance;
#' 1 is perfect separation.
#'
#' @param prob Numeric vector of predicted probabilities (or any score where
#'   larger means more case-like).
#' @param y Observed binary outcome as 0/1 numeric or logical, the same length as
#'   \code{prob}.
#'
#' @return A single number in \code{[0, 1]}, or \code{NA_real_} if either class is
#'   absent.
#'
#' @references
#' Merlo, J., Wagner, P., Ghith, N., & Leckie, G. (2016). An original stepwise
#' multilevel logistic regression analysis of discriminatory accuracy: the case of
#' neighbourhoods and health. \emph{PLOS ONE}, 11(4), e0153778.
#'
#' @examples
#' maihda_auc(c(0.1, 0.4, 0.35, 0.8), c(0, 0, 1, 1))
#'
#' @export
maihda_auc <- function(prob, y) {
  if (!is.numeric(prob)) {
    stop("'prob' must be a numeric vector of predicted probabilities/scores.",
         call. = FALSE)
  }
  if (is.logical(y)) {
    y <- as.integer(y)
  }
  y <- suppressWarnings(as.numeric(y))
  if (length(prob) != length(y)) {
    stop("'prob' and 'y' must have the same length (", length(prob), " vs ",
         length(y), ").", call. = FALSE)
  }

  keep <- !(is.na(prob) | is.na(y))
  prob <- prob[keep]
  y <- y[keep]

  if (!all(y %in% c(0, 1))) {
    stop("'y' must be a binary 0/1 (or logical) outcome. For a factor outcome, ",
         "convert it to 0/1 first (e.g. as.numeric(factor) - 1).", call. = FALSE)
  }

  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) {
    return(NA_real_)
  }

  r <- rank(prob)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

#' Median Odds Ratio (MOR) for a logistic MAIHDA model
#'
#' @description
#' The Median Odds Ratio translates the between-stratum variance of a logistic
#' MAIHDA model onto the odds-ratio scale: the median relative change in the odds
#' of the outcome when comparing two individuals from randomly chosen strata
#' (higher- vs lower-risk). \code{MOR = exp(sqrt(2 * V_A) * qnorm(0.75))}, where
#' \code{V_A} is the between-stratum (latent, logit-scale) variance. An MOR of 1
#' indicates no between-stratum heterogeneity. The MOR is defined only for the
#' \strong{logit} link (it is the median \emph{odds} ratio); a non-logit binomial
#' fit such as \code{probit} is rejected, because its latent variance is on a
#' different scale and the \code{exp(...)} above would not be an odds ratio.
#'
#' @param model A \code{maihda_model} from \code{\link{fit_maihda}} fitted with a
#'   \code{binomial} family and a \strong{logit} link.
#'
#' @return A single number (the MOR, \eqn{\ge 1}), or \code{NA_real_} if the
#'   between-stratum variance is unavailable.
#'
#' @references
#' Larsen, K., & Merlo, J. (2005). Appropriate assessment of neighborhood effects
#' on individual health: integrating random and fixed effects in multilevel
#' logistic regression. \emph{American Journal of Epidemiology}, 161(1), 81-88.
#'
#' @seealso \code{\link{maihda_discriminatory_accuracy}}
#'
#' @export
maihda_mor <- function(model) {
  if (!inherits(model, "maihda_model")) {
    stop("'model' must be a maihda_model object from fit_maihda().", call. = FALSE)
  }
  fam <- maihda_model_family_name(model)
  if (!identical(fam, "binomial")) {
    stop("The Median Odds Ratio is only defined for binomial (logistic) MAIHDA ",
         "models; this model uses family = '", fam, "'.", call. = FALSE)
  }
  link <- maihda_model_link_name(model)
  if (!identical(link, "logit")) {
    stop("The Median Odds Ratio is defined only for the logit link -- it is the ",
         "median *odds* ratio, derived from the logistic latent variance -- but this ",
         "model uses a binomial '", link, "' link, whose latent variance is on a ",
         "different scale. Use maihda_auc() / maihda_discriminatory_accuracy() for ",
         "link-agnostic discriminatory accuracy.", call. = FALSE)
  }

  v_a <- tryCatch(extract_between_variance(model), error = function(e) NA_real_)
  if (!is.numeric(v_a) || length(v_a) != 1 || !is.finite(v_a) || v_a < 0) {
    return(NA_real_)
  }

  exp(sqrt(2 * v_a) * stats::qnorm(0.75))
}

#' Discriminatory accuracy of a binary MAIHDA model
#'
#' @description
#' Bundles the individual-level discriminatory-accuracy summaries for a binomial
#' MAIHDA model: the AUC / C-statistic (how well the model's predicted
#' probabilities separate cases from non-cases) and the Median Odds Ratio. Applied
#' to a strata-only (null) model, the AUC is the discriminatory accuracy of the
#' intersectional strata themselves -- Merlo's central quantity; comparing it with
#' an adjusted model shows whether individual covariates beyond stratum membership
#' sharpen classification. The AUC is computed for any binomial link; the Median
#' Odds Ratio is reported only for the logit link and is \code{NA} otherwise (e.g.
#' for a probit fit), since the MOR is an odds-ratio-scale quantity.
#'
#' @param model A \code{maihda_model} from \code{\link{fit_maihda}} fitted with a
#'   \code{binomial} family (lme4 engine).
#'
#' @return An object of class \code{maihda_da}: a list with \code{auc}, \code{mor},
#'   \code{n_case}, \code{n_control}, \code{family}, \code{link} and \code{engine}.
#'   \code{mor} is \code{NA} for a non-logit binomial link, where the AUC is still
#'   reported.
#'
#' @references
#' Merlo, J. (2018). Multilevel analysis of individual heterogeneity and
#' discriminatory accuracy (MAIHDA) within an intersectional framework.
#' \emph{Social Science & Medicine}, 203, 74-80.
#'
#' @seealso \code{\link{maihda_auc}}, \code{\link{maihda_mor}}
#'
#' @examples
#' \dontrun{
#' # Obese (Yes/No) by intersectional strata of Gender x Race
#' strata <- make_strata(maihda_health_data, vars = c("Gender", "Race"))
#' d <- maihda_health_data
#' d$stratum <- strata$data$stratum
#' m <- fit_maihda(Obese ~ (1 | stratum), data = d, family = "binomial")
#' maihda_discriminatory_accuracy(m)
#' }
#'
#' @export
maihda_discriminatory_accuracy <- function(model) {
  if (!inherits(model, "maihda_model")) {
    stop("'model' must be a maihda_model object from fit_maihda().", call. = FALSE)
  }
  fam <- maihda_model_family_name(model)
  if (!identical(fam, "binomial")) {
    stop("Discriminatory accuracy (AUC / MOR) is only defined for binomial ",
         "(logistic) MAIHDA models; this model uses family = '", fam, "'.",
         call. = FALSE)
  }

  prob <- predict_maihda(model, type = "individual", scale = "response")
  y <- maihda_da_observed_response(model)

  auc <- maihda_auc(prob, y)
  # The AUC is link-agnostic (rank-based on predicted probabilities), but the MOR is
  # defined only for the logit link. For other binomial links (e.g. probit) report
  # the AUC with mor = NA rather than an odds ratio that is off the model's scale.
  link <- maihda_model_link_name(model)
  mor <- if (identical(link, "logit")) maihda_mor(model) else NA_real_

  structure(
    list(
      auc = auc,
      mor = mor,
      n_case = sum(y == 1, na.rm = TRUE),
      n_control = sum(y == 0, na.rm = TRUE),
      family = fam,
      link = link,
      engine = model$engine
    ),
    class = "maihda_da"
  )
}

#' @export
print.maihda_da <- function(x, ...) {
  cat("Discriminatory accuracy (binomial MAIHDA)\n")
  cat(sprintf("  AUC (C-statistic): %s\n",
              if (is.finite(x$auc)) sprintf("%.3f", x$auc) else "NA"))
  mor_str <- if (is.finite(x$mor)) {
    sprintf("%.3f", x$mor)
  } else if (!is.null(x$link) && !identical(x$link, "logit")) {
    sprintf("NA (requires the logit link; model uses '%s')", x$link)
  } else {
    "NA"
  }
  cat(sprintf("  Median Odds Ratio: %s\n", mor_str))
  cat(sprintf("  Cases / controls:  %d / %d\n", x$n_case, x$n_control))
  invisible(x)
}

# ---- internal helpers -------------------------------------------------------

# Resolve the family name ("binomial"/"gaussian"/"poisson") of a maihda_model,
# tolerating either a stored family object/string or falling back to the fitted
# model's family.
maihda_model_family_name <- function(model) {
  fam <- model$family
  if (is.list(fam) && !is.null(fam$family)) {
    return(fam$family)
  }
  if (is.character(fam) && length(fam) == 1) {
    return(fam)
  }
  ff <- tryCatch(maihda_family(model$model), error = function(e) NULL)
  if (!is.null(ff) && !is.null(ff$family)) {
    return(ff$family)
  }
  NA_character_
}

# Resolve the link name ("logit"/"probit"/...) of a maihda_model, preferring the
# stored family object and falling back to the fitted model's family. Used to gate
# the Median Odds Ratio, which is defined only for the logit link.
maihda_model_link_name <- function(model) {
  fam <- model$family
  if (is.list(fam) && !is.null(fam$link)) {
    return(fam$link)
  }
  ff <- tryCatch(maihda_family(model$model), error = function(e) NULL)
  if (!is.null(ff) && !is.null(ff$link)) {
    return(ff$link)
  }
  NA_character_
}

# Observed 0/1 response aligned with predict_maihda()'s individual predictions.
# For lme4, getME(, "y") returns the numeric 0/1 response used in fitting (the
# approach used in the binary_outcomes vignette). For other engines, fall back to
# the response column of the model frame, coerced to 0/1.
maihda_da_observed_response <- function(model) {
  if (identical(model$engine, "lme4")) {
    return(as.numeric(lme4::getME(model$model, "y")))
  }

  resp <- all.vars(model$formula)[1]
  y <- model$data[[resp]]
  if (is.logical(y)) {
    return(as.integer(y))
  }
  if (is.factor(y)) {
    return(as.integer(y) - 1L)
  }
  as.numeric(y)
}
