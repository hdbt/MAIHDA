#' Calculate Proportional Change in Between-Stratum Variance (PCV)
#'
#' Calculates the proportional change in between-stratum variance (PCV) between
#' two MAIHDA models. The PCV measures how much the between-stratum variance
#' changes when moving from one model to another, and is calculated as:
#' PCV = (Var_model1 - Var_model2) / Var_model1.
#'
#' @param model1 A maihda_model object from \code{fit_maihda()}. This is the
#'   reference model (typically a simpler or baseline model).
#' @param model2 A maihda_model object from \code{fit_maihda()}. This is the
#'   comparison model (typically a more complex model with additional predictors).
#' @param bootstrap Logical indicating whether to compute bootstrap confidence
#'   intervals for the PCV. Default is FALSE. \strong{lme4 engine only}: the
#'   parametric bootstrap relies on lme4's \code{simulate()}/\code{refit()}, so
#'   for the brms, wemix, and ordinal engines the PCV is reported as a point
#'   estimate and \code{bootstrap = TRUE} is an error (see Details).
#' @param n_boot Number of bootstrap samples if bootstrap = TRUE. Default is 1000.
#' @param conf_level Confidence level for bootstrap intervals. Default is 0.95.
#'
#' @return A list containing:
#'   \item{pcv}{The estimated proportional change in variance}
#'   \item{pvc}{Deprecated duplicate of \code{pcv}, kept so code written against
#'     the historical \code{calculate_pvc()} spelling keeps working; it will be
#'     removed in a future release}
#'   \item{var_model1}{Between-stratum variance from model1}
#'   \item{var_model2}{Between-stratum variance from model2}
#'   \item{ci_lower}{Lower bound of confidence interval (if bootstrap = TRUE)}
#'   \item{ci_upper}{Upper bound of confidence interval (if bootstrap = TRUE)}
#'   \item{bootstrap}{Logical indicating if bootstrap was used}
#'
#' @details
#' The PCV is the proportional change in between-stratum variance when moving from
#' model1 to model2: a positive value means model2 has lower between-stratum
#' variance, a negative value means higher. It is the share of model1's
#' between-stratum variance \emph{explained} by model2 only in the canonical nested
#' case, where model2 adds fixed-effect predictors to model1 on the same outcome,
#' analytic sample and strata. The function does not require nesting, so for
#' non-nested models the PCV is simply a model-dependent difference in variance,
#' not an explained proportion.
#'
#' \strong{REML vs ML.} \code{lmer} fits Gaussian models by REML, whose
#' between-stratum variance estimate is \emph{not} comparable across models with
#' different fixed effects -- exactly the canonical null-vs-adjusted PCV, where the
#' adjusted model adds the dimensions' main effects. \code{calculate_pcv()} therefore
#' refits any REML \code{lmer} model with maximum likelihood
#' (\code{\link[lme4]{refitML}}) before reading the variances (and before the
#' parametric bootstrap, so the interval matches), matching \code{\link{maihda_ic}}
#' and \code{anova()} on \code{lme4} models. Using REML estimates here biases the PCV
#' (it overstates the residual between-stratum variance of the adjusted model). GLMM
#' fits (\code{glmer}) and the brms/wemix/ordinal engines are already on the
#' maximum-likelihood scale and are unaffected; single-model VPC/ICC summaries keep
#' their REML fit, since that comparison-free quantity is not subject to the pitfall.
#'
#' \strong{Latent-scale families and rescaling.} For families whose level-1
#' variance is a fixed latent-scale constant -- binomial/Bernoulli
#' (\eqn{\pi^2/3} logit, 1 probit) and the cumulative (ordinal) model -- the
#' linear predictor is identified only up to scale. Adding predictors that
#' explain \emph{within-stratum} (individual-level) variation cannot shrink that
#' fixed level-1 variance; the latent scale stretches instead, inflating the
#' coefficients and the between-stratum variance alike (Bauer 2009; Mood 2010).
#' Part of a null-vs-adjusted change in the between-stratum variance is then
#' rescaling rather than genuinely explained variance, so latent-scale PCVs tend
#' to be understated and can turn negative on this account alone. The canonical
#' MAIHDA adjusted model -- which adds the stratum dimensions' main effects,
#' constant \emph{within} each stratum -- is largely unaffected, but the caveat
#' is first-order whenever an added predictor varies within strata (an
#' individual-level covariate, as in the \code{\link{stepwise_pcv}} steps that
#' add one). The count families' level-1 variance is not a fixed constant, but
#' as with any non-identity link the same non-collapsibility logic applies in
#' attenuated form. Gaussian identity-link PCVs are not subject to this.
#'
#' When bootstrap = TRUE, the function uses a parametric bootstrap: it simulates
#' new responses from model2 and refits both models with \code{lme4::refit()} for
#' each simulated response to obtain confidence intervals for the PCV estimate.
#' For negative-binomial models (\code{glmer.nb}) \code{refit()} holds the
#' dispersion parameter theta fixed at its original estimate, so the interval is
#' conditional on the estimated theta.
#'
#' The bootstrap is available for the \code{lme4} engine only. For the other
#' engines the PCV is a \emph{point estimate}: a brms fit's posterior credible
#' interval (reported by \code{\link{summary.maihda_model}}) covers a single
#' fit's VPC/ICC, not the PCV, which compares two separately fitted models -- no
#' posterior interval for the PCV itself is computed -- and a design-based
#' (wemix) interval would require replicate weights.
#'
#' @examples
#' \donttest{
#' # Create strata and fit two models
#' strata_result <- make_strata(maihda_sim_data, c("gender", "race"))
#' model1 <- fit_maihda(health_outcome ~ age + (1 | stratum), data = strata_result$data)
#' model2 <- fit_maihda(health_outcome ~ age + gender + (1 | stratum), data = strata_result$data)
#'
#' # Calculate PCV without bootstrap
#' pcv_result <- calculate_pcv(model1, model2)
#' print(pcv_result$pcv)
#'
#' # Calculate PCV with bootstrap CI
#' # pcv_boot <- calculate_pcv(model1, model2, bootstrap = TRUE, n_boot = 500)
#' # print(pcv_boot)
#' }
#'
#' @seealso \code{\link{stepwise_pcv}} for the sequential (one-variable-at-a-time)
#'   PCV, and \code{\link{maihda}} which computes the canonical null-vs-adjusted
#'   PCV automatically.
#'
#' @references
#' Bauer, D. J. (2009). A note on comparing the estimates of models for
#' cluster-correlated or longitudinal data with binary or ordinal outcomes.
#' \emph{Psychometrika}, 74(1), 97-105.
#'
#' Mood, C. (2010). Logistic regression: why we cannot do what we think we can
#' do, and what we can do about it. \emph{European Sociological Review}, 26(1),
#' 67-82.
#' @export
#' @importFrom lme4 lmer glmer VarCorr
calculate_pcv <- function(model1, model2, bootstrap = FALSE,
                         n_boot = 1000, conf_level = 0.95) {
  # Input validation
  if (!inherits(model1, "maihda_model")) {
    stop("'model1' must be a maihda_model object from fit_maihda()")
  }

  if (!inherits(model2, "maihda_model")) {
    stop("'model2' must be a maihda_model object from fit_maihda()")
  }

  if (model1$engine != model2$engine) {
    stop("Both models must use the same engine.")
  }

  if (!is.logical(bootstrap) || length(bootstrap) != 1 || is.na(bootstrap)) {
    stop("'bootstrap' must be TRUE or FALSE.", call. = FALSE)
  }
  if (bootstrap) {
    bootstrap_args <- maihda_validate_bootstrap_args(n_boot, conf_level)
    n_boot <- bootstrap_args$n_boot
    conf_level <- bootstrap_args$conf_level
  }

  validate_pcv_models(model1, model2)

  # REML vs ML: lmer fits Gaussian models by REML, whose between-stratum variance is
  # NOT comparable across models with different fixed effects -- exactly the canonical
  # null-vs-adjusted PCV. Refit any REML lmer fit with ML before the comparison (and
  # before the parametric bootstrap below, which reuses these fits, so the interval
  # matches the point estimate), mirroring maihda_ic() and anova.merMod. GLMM fits
  # (glmer) and the brms/wemix/ordinal engines are already ML / unaffected.
  model1 <- maihda_pcv_refit_ml(model1)
  model2 <- maihda_pcv_refit_ml(model2)

  # Extract between-stratum variance from both models
  var1 <- extract_between_variance(model1)
  var2 <- extract_between_variance(model2)

  # Validate variances
  if (is.na(var1) || is.na(var2)) {
    stop("Unable to extract variance components from one or both models")
  }

  if (var1 <= 0) {
    stop("Between-stratum variance in model1 is zero or negative. PCV cannot be calculated. ",
         "This may indicate a singular fit or no between-stratum variation.")
  }

  # Calculate PCV
  pcv <- (var1 - var2) / var1

  # Create result object. The 'pvc' element and the trailing "pvc_result" class
  # duplicate 'pcv'/"pcv_result" for code written against the historical
  # calculate_pvc() spelling.
  result <- list(
    pcv = pcv,
    pvc = pcv,
    var_model1 = var1,
    var_model2 = var2,
    bootstrap = FALSE
  )

  # Bootstrap confidence intervals if requested
  if (bootstrap) {
    pcv_ci <- bootstrap_pcv(model1, model2, n_boot, conf_level)
    result$ci_lower <- pcv_ci[1]
    result$ci_upper <- pcv_ci[2]
    result$bootstrap <- TRUE
    result$conf_level <- conf_level
    result$n_boot_ok <- attr(pcv_ci, "n_ok")
    result$mc_se <- attr(pcv_ci, "mc_se")
  }

  class(result) <- c("pcv_result", "pvc_result")
  return(result)
}

#' Deprecated: use calculate_pcv()
#'
#' \code{calculate_pvc()} is the former name of \code{\link{calculate_pcv}}: the
#' statistic is the PCV (proportional change in variance), but the historical
#' function name transposed the acronym. \code{calculate_pvc()} now forwards to
#' \code{calculate_pcv()} with a deprecation warning and will be removed in a
#' future release.
#'
#' @inheritParams calculate_pcv
#' @return See \code{\link{calculate_pcv}}.
#' @export
calculate_pvc <- function(model1, model2, bootstrap = FALSE,
                          n_boot = 1000, conf_level = 0.95) {
  .Deprecated("calculate_pcv")
  calculate_pcv(model1, model2, bootstrap = bootstrap,
                n_boot = n_boot, conf_level = conf_level)
}

# REML lmer between-stratum variance estimates are not comparable across models with
# different fixed effects (the canonical null-vs-adjusted PCV), because the REML
# criterion conditions on the fixed-effects design. Refit a REML lmer fit with ML
# (lme4::refitML) before any cross-model variance comparison, matching maihda_ic() and
# anova.merMod. Non-REML fits (glmer / the GLMM families) and the brms/wemix/ordinal
# engines are returned unchanged. Longitudinal fits are returned unchanged TOO, but
# only because their null-vs-adjusted comparison applies the same rule itself:
# maihda_longitudinal_pcv() ML-refits its growth models via
# maihda_longitudinal_refit_ml() (longitudinal.R), whose boundary guard reads the
# whole stratum growth block rather than the single intercept cell checked here.
# Single-model VPC/ICC summaries deliberately keep their REML fit.
maihda_pcv_refit_ml <- function(model) {
  if (!inherits(model, "maihda_model") || !identical(model$engine, "lme4") ||
      !is.null(model$longitudinal_info)) {
    return(model)
  }
  is_reml <- tryCatch(isTRUE(lme4::isREML(model$model)), error = function(e) FALSE)
  if (!is_reml) return(model)
  # Skip the ML refit only when the STRATUM variance is itself on the boundary
  # (effectively zero): there its between-stratum variance is ~0 under either
  # criterion, re-optimising only adds optimiser instability, and the refit would
  # nudge an exact-zero variance off the boundary, masking the zero-variance guard in
  # calculate_pcv(). A fit can be globally singular because a NON-stratum random
  # effect (e.g. an extra grouping factor like (1 | site)) sits on the boundary while
  # the stratum variance is comfortably nonzero -- testing global isSingular() here
  # wrongly skipped those, leaving the REML-vs-ML discrepancy this corrects in place.
  # If refitML() itself fails, the tryCatch below keeps the original (REML) fit.
  if (maihda_stratum_at_boundary_lme4(model$model)) return(model)
  refit <- tryCatch(lme4::refitML(model$model), error = function(e) NULL)
  if (!is.null(refit)) model$model <- refit
  model
}

# TRUE when the stratum random-intercept variance sits on the lower boundary
# (effectively zero). Reproduces lme4::isSingular()'s relative tolerance for the
# stratum component alone: its scaled SD parameter theta = sd_stratum / sigma_resid
# is below `tol`, which is exactly how lme4 flags a single intercept term as singular.
# Used by maihda_pcv_refit_ml() to decide whether an ML refit is worthwhile -- it is
# not when the stratum is itself at zero, but it IS when only another grouping factor
# is on the boundary. Returns TRUE (skip the refit) if the stratum variance cannot be
# read at all; calculate_pcv()'s own zero-variance guard then surfaces the problem.
maihda_stratum_at_boundary_lme4 <- function(model, tol = 1e-4) {
  stratum_var <- tryCatch(maihda_stratum_variance_lme4(model),
                          error = function(e) NA_real_)
  if (is.na(stratum_var) || stratum_var <= 0) return(TRUE)
  sigma_resid <- tryCatch(stats::sigma(model), error = function(e) NA_real_)
  if (is.na(sigma_resid) || sigma_resid <= 0) {
    # No usable residual scale: fall back to an absolute near-zero test.
    return(stratum_var <= .Machine$double.eps)
  }
  sqrt(stratum_var) / sigma_resid < tol
}

#' Extract Between-Stratum Variance
#'
#' Internal function to extract between-stratum variance from a MAIHDA model.
#'
#' @param model A maihda_model object
#'
#' @return Numeric value of between-stratum variance
#' @keywords internal
#' @importFrom lme4 VarCorr
extract_between_variance <- function(model) {
  engine <- model$engine
  fitted_model <- model$model

  # A longitudinal (growth-curve) fit has a random intercept AND slope on the
  # stratum, so the between-stratum variance is a function of time -- there is no
  # single scalar to return. Route the user to the time-varying decomposition.
  if (!is.null(model$longitudinal_info)) {
    stop("This is a longitudinal MAIHDA: the between-stratum variance is ",
         "time-varying (random intercept + slope), so a single PCV/VPC scalar is ",
         "undefined. Use maihda(decomposition = \"longitudinal\") for the ",
         "additive-vs-multiplicative PCV (baseline and slope), or summary() for ",
         "the time-varying VPC.", call. = FALSE)
  }

  if (engine == "lme4") {
    maihda_validate_intercept_only_random_effects_lme4(
      fitted_model,
      context = "PCV calculations"
    )
    return(maihda_stratum_variance_lme4(fitted_model))

  } else if (engine == "wemix") {
    # The wemix engine only ever fits the canonical single intercept-only
    # (1 | stratum) structure (enforced at fit time), so no random-slope
    # validation is needed here.
    return(maihda_wemix_variances(model)$stratum)

  } else if (engine == "ordinal") {
    # Like wemix, the ordinal engine enforces the canonical single
    # intercept-only (1 | stratum) structure at fit time.
    return(maihda_clmm_variances(model)$stratum)

  } else if (engine == "brms") {
    if (!requireNamespace("brms", quietly = TRUE)) {
      stop("Package 'brms' is required to work with brms models. Please install it with: install.packages('brms')")
    }
    maihda_validate_intercept_only_random_effects_brms(
      brms::VarCorr(fitted_model),
      context = "PCV calculations"
    )
    return(maihda_stratum_variance_brms(fitted_model))

  } else {
    stop("Unsupported engine: ", engine,
         ". Supported engines are 'lme4', 'brms', 'wemix', and 'ordinal'.")
  }
}

validate_pcv_models <- function(model1, model2) {
  response1 <- paste(deparse(model1$formula[[2]]), collapse = "")
  response2 <- paste(deparse(model2$formula[[2]]), collapse = "")
  if (!identical(response1, response2)) {
    stop("PCV requires both models to use the same outcome. ",
         "Model 1 uses '", response1, "' and Model 2 uses '", response2, "'.",
         call. = FALSE)
  }

  # Canonical "family(link)" keys (maihda_model_family_key falls back to the
  # wrapper-recorded family for engines where stats::family() is undefined, and
  # normalises labels so e.g. two glmer.nb fits with different estimated thetas
  # -- reported as "Negative Binomial(<theta>)" -- still compare equal).
  fam_key1 <- maihda_model_family_key(model1)
  fam_key2 <- maihda_model_family_key(model2)
  if (!identical(fam_key1, fam_key2)) {
    stop("PCV requires both models to use the same model family and link. ",
         "Model 1 uses ", fam_key1, " and Model 2 uses ", fam_key2, ".",
         call. = FALSE)
  }

  # nobs()/model.frame() are undefined for WeMixResults; the wrapper helpers fall
  # back to the stored analytic frame (the fitted rows) so the n, row-identity and
  # response checks below still apply to those engines instead of degrading to a
  # silent pass.
  n1 <- maihda_wrapper_nobs(model1)
  n2 <- maihda_wrapper_nobs(model2)
  if (is.finite(n1) && is.finite(n2) && n1 != n2) {
    stop("PCV requires both models to use the same analytic sample. ",
         "Model 1 used ", n1, " observations and Model 2 used ", n2, ".",
         call. = FALSE)
  }

  rows1 <- maihda_wrapper_row_ids(model1)
  rows2 <- maihda_wrapper_row_ids(model2)
  if (!is.null(rows1) && !is.null(rows2) && !identical(rows1, rows2)) {
    stop("PCV requires both models to use the same analytic sample in the same row order.",
         call. = FALSE)
  }

  # Content fingerprint: catch unrelated datasets that share n and default 1:n
  # row names but hold different responses.
  fp1 <- maihda_wrapper_response_fingerprint(model1)
  fp2 <- maihda_wrapper_response_fingerprint(model2)
  if (!is.na(fp1) && !is.na(fp2) && !identical(fp1, fp2)) {
    stop("PCV requires both models to be fitted to the same analytic data; the ",
         "outcome values differ between the two models.", call. = FALSE)
  }

  # Prior weights change the variance estimates, so a PCV across fits that used
  # different weights would compare incomparable quantities. (An unweighted fit and
  # an explicit weights = rep(1, n) fit are treated as equal.)
  if (!identical(maihda_weight_fingerprint(model1$model),
                 maihda_weight_fingerprint(model2$model))) {
    stop("PCV requires both models to use the same prior weights; the two models ",
         "were fit with different weights.", call. = FALSE)
  }

  # Same idea for SAMPLING weights (design-weighted fits): a PCV across fits that
  # used different design weights -- or one weighted and one unweighted fit --
  # would compare incomparable variance estimates. The fingerprint covers the
  # column name and its values on the analytic rows.
  if (!identical(maihda_sampling_weight_fingerprint(model1),
                 maihda_sampling_weight_fingerprint(model2))) {
    stop("PCV requires both models to use the same sampling weights; the two ",
         "models were fit with different (or differently specified) ",
         "sampling weights.", call. = FALSE)
  }

  if (!"stratum" %in% names(model1$data) || !"stratum" %in% names(model2$data)) {
    stop("PCV requires both models to include a 'stratum' column in their analytic data.",
         call. = FALSE)
  }
  row_strata1 <- as.character(model1$data$stratum)
  row_strata2 <- as.character(model2$data$stratum)
  if (!identical(row_strata1, row_strata2)) {
    stop("PCV requires both models to assign each analytic row to the same stratum.",
         call. = FALSE)
  }

  strata1 <- unique(as.character(model1$data$stratum))
  strata2 <- unique(as.character(model2$data$stratum))
  strata1 <- sort(strata1[!is.na(strata1)])
  strata2 <- sort(strata2[!is.na(strata2)])
  if (!identical(strata1, strata2)) {
    stop("PCV requires both models to use the same stratum definitions.",
         call. = FALSE)
  }

  info1 <- model1$strata_info
  info2 <- model2$strata_info
  if (!is.null(info1) && !is.null(info2) &&
      all(c("stratum", "label") %in% names(info1)) &&
      all(c("stratum", "label") %in% names(info2))) {
    labels1 <- info1$label[order(as.character(info1$stratum))]
    labels2 <- info2$label[order(as.character(info2$stratum))]
    if (!identical(labels1, labels2)) {
      stop("PCV requires both models to use the same stratum labels.",
           call. = FALSE)
    }
  }

  invisible(TRUE)
}

#' Bootstrap PCV
#'
#' Internal function to compute bootstrap confidence intervals for PCV.
#'
#' @param model1 First maihda_model object
#' @param model2 Second maihda_model object
#' @param n_boot Number of bootstrap samples
#' @param conf_level Confidence level
#'
#' @return A vector with lower and upper confidence bounds
#' @keywords internal
#' @importFrom lme4 lmer glmer VarCorr
bootstrap_pcv <- function(model1, model2, n_boot, conf_level) {
  engine <- model1$engine
  if (engine != "lme4") {
    # Engine-specific, honest guidance. The PCV compares TWO separately fitted
    # models, so a brms fit's posterior credible interval -- which summary()
    # reports for a single fit's VPC/ICC -- does not carry over to the PCV;
    # calculate_pcv() reports the PCV as a point estimate for every non-lme4
    # engine. The actionable advice is therefore bootstrap = FALSE, not an
    # engine switch (the old message told brms users to "refit with brms").
    hint <- switch(engine,
      brms = paste0(
        "For brms models the PCV is reported as a point estimate (the ratio of ",
        "posterior-mean between-stratum variances of two separately fitted ",
        "models); no posterior interval for the PCV itself is computed -- ",
        "summary()'s credible interval covers a single fit's VPC/ICC, not the ",
        "PCV. Call calculate_pcv() with bootstrap = FALSE."),
      wemix = paste0(
        "A design-based PCV interval would require replicate weights (not ",
        "implemented); the design-weighted PCV is reported as a point ",
        "estimate. Call calculate_pcv() with bootstrap = FALSE."),
      ordinal = paste0(
        "ordinal::clmm has no simulate()/refit() machinery; the PCV is ",
        "reported as a point estimate. Call calculate_pcv() with ",
        "bootstrap = FALSE."),
      paste0("The PCV is reported as a point estimate for this engine; call ",
             "calculate_pcv() with bootstrap = FALSE.")
    )
    stop("Bootstrap PCV intervals are only available for lme4 models (the ",
         "parametric bootstrap relies on lme4's simulate()/refit()). ", hint,
         call. = FALSE)
  }

  # Initialise to NA so iterations whose refit() throws — and never reach the
  # assignment inside the tryCatch body — stay NA rather than the numeric() default of 0.
  # The error handler runs in its own scope and cannot write back to this vector,
  # so the initial value is what survives a failure.
  pcv_boot <- rep(NA_real_, n_boot)

  # Parametric Bootstrap: Simulate new responses from the adjusted model (model2)
  # This mathematically preserves the hierarchical structure (random effects)
  # and the fixed-effects distributions, unlike naive row-resampling.
  sim_data <- stats::simulate(model2$model, nsim = n_boot)

  for (i in 1:n_boot) {
    tryCatch({
      # Fast parametric refitting with the newly simulated response vector
      boot_model1 <- lme4::refit(model1$model, newresp = sim_data[[i]])
      boot_model2 <- lme4::refit(model2$model, newresp = sim_data[[i]])

      # Extract variances
      var1 <- maihda_stratum_variance_lme4(boot_model1)
      var2 <- maihda_stratum_variance_lme4(boot_model2)

      # Calculate PCV
      pcv_boot[i] <- if (is.finite(var1) && var1 > 0) (var1 - var2) / var1 else NA_real_
    }, error = function(e) NULL)
  }

  # Reduce to an interval, requiring a minimum number of successful refits and
  # warning on a high failure rate.
  ci <- maihda_bootstrap_ci(pcv_boot, n_boot, conf_level, "PCV")

  return(ci)
}

#' Print method for PCV results
#'
#' @param x A pcv_result object from \code{\link{calculate_pcv}}
#' @param ... Additional arguments
#' @return No return value, called for side effects.
#' @export
print.pcv_result <- function(x, ...) {
  # An object saved by a pre-rename version of the package carries only 'pvc'.
  pcv <- if (is.null(x$pcv)) x$pvc else x$pcv
  pal <- maihda_palette()
  cat(pal$bold("Proportional Change in Variance (PCV)"), "\n", sep = "")
  cat("=====================================\n\n")

  if (x$bootstrap) {
    conf_pct <- if (!is.null(x$conf_level)) x$conf_level * 100 else 95
    cat(sprintf("PCV: %s [%.4f, %.4f]\n",
                pal$accent(sprintf("%.4f", pcv)), x$ci_lower, x$ci_upper))
    cat(pal$muted(sprintf("(Bootstrap %.0f%% CI)\n", conf_pct)))
    if (!is.null(x$mc_se) && is.finite(x$mc_se)) {
      cat(pal$muted(sprintf("(%d successful bootstrap draws; Monte Carlo SE %.4f)\n",
                  as.integer(x$n_boot_ok), x$mc_se)))
    }
    cat("\n")
  } else {
    cat(sprintf("PCV: %s\n\n", pal$accent(sprintf("%.4f", pcv))))
  }

  cat(pal$bold("Between-stratum variance:"), "\n", sep = "")
  cat(sprintf("  Model 1: %.6f\n", x$var_model1))
  cat(sprintf("  Model 2: %.6f\n", x$var_model2))
  cat(sprintf("  Change:  %.6f (%.2f%%)\n",
              x$var_model1 - x$var_model2,
              pcv * 100))

  cat(pal$muted("\nInterpretation (PCV is the proportional change in between-stratum\nvariance between the models):\n"))
  if (pcv > 0) {
    cat(sprintf("  Between-stratum variance is %.1f%% lower in Model 2 than in Model 1.\n",
                pcv * 100))
  } else if (pcv < 0) {
    cat(sprintf("  Between-stratum variance is %.1f%% higher in Model 2 than in Model 1.\n",
                abs(pcv) * 100))
  } else {
    cat("  No change in between-stratum variance between models.\n")
  }

  invisible(x)
}

# Objects saved by a pre-rename version of the package carry only the "pvc_result"
# class, so the legacy S3 method must stay registered for them to print.
#' @rdname print.pcv_result
#' @export
print.pvc_result <- print.pcv_result

# Total between-context variance held in a stepwise model (the crossed contextual
# random intercept(s) added by fit_maihda(context = )). Reads the per-group intercept
# variances and sums the context component(s), so multiple context columns collapse to
# one total -- the share reported alongside the net-of-context stratum PCV. Reads from
# the same fitted object stepwise_pcv() uses for the stratum variance (the ML-refit fit
# for a REML lmer model), keeping the two variances on one scale. Returns NA for an
# engine without a recoverable per-group split (wemix/ordinal never reach here: they
# reject context up front).
maihda_stepwise_context_variance <- function(mod, context) {
  v <- switch(mod$engine,
    lme4 = tryCatch(maihda_random_variances_lme4(mod$model),
                    error = function(e) NULL),
    brms = tryCatch(maihda_random_variances_brms(mod$model),
                    error = function(e) NULL),
    NULL)
  if (is.null(v)) {
    return(NA_real_)
  }
  present <- intersect(context, names(v))
  if (length(present) == 0) {
    return(NA_real_)
  }
  sum(v[present], na.rm = TRUE)
}

#' Stepwise Proportional Change in Variance (PCV)
#'
#' @description
#' Estimates the proportional change in variance (PCV) sequentially by fitting
#' intermediate (partially-adjusted) models, adding each predictor one-by-one. The
#' step-specific PCV is the change in between-stratum variance contributed by a
#' predictor \emph{given the variables already in the model}. Because the steps are
#' sequential it is order-dependent: it reflects each variable's marginal,
#' model-dependent change, not an order-invariant \dQuote{unique} contribution.
#'
#' @param data Data frame with observations. Ensure `make_strata()` was run first
#'   so the `stratum` variable exists.
#' @param outcome Character string; the dependent variable.
#' @param vars Character vector; predictors (strata groupings & covariates) to
#'   add sequentially to the model.
#' @param engine Modeling engine ("lme4", "brms", "wemix", or "ordinal"). Default
#'   is "lme4"; switches to "wemix" automatically when \code{sampling_weights} is
#'   supplied, and to "ordinal" for an ordinal family or ordered-factor outcome.
#' @param family Error distribution and link function. Default is "gaussian".
#' @param context Optional character vector naming one or more higher-level
#'   \emph{context} columns in \code{data} (e.g. \code{"school"}, \code{"region"}),
#'   forwarded to every step's \code{\link{fit_maihda}} call so that each model
#'   carries a crossed contextual random intercept, \code{(1 | context)}, alongside
#'   the intersectional stratum effect -- the stepwise analogue of
#'   \code{maihda(context = )} (the literature's contextual cross-classified MAIHDA).
#'   The reported \code{Step_PCV} / \code{Total_PCV} are then the between-stratum PCV
#'   \strong{net of} the context (the context intercept is held in the null model and
#'   every step, so the stratum variance is isolated from it), and an extra
#'   \code{Context_Variance} column reports the between-context variance held at each
#'   step. Context columns join the complete-case filter so every step shares one
#'   analytic sample. A context variable may not be a stratum dimension, \code{"stratum"}
#'   itself, or one of \code{vars}. Supported by the \code{lme4} and \code{brms} engines
#'   only (\code{wemix} and \code{ordinal} fit no crossed random effects). See
#'   \code{\link{fit_maihda}} for the full contextual-model semantics. Default
#'   \code{NULL} (no context).
#' @param sampling_weights Optional name of a sampling-weight column for
#'   design-weighted stepwise fits; see \code{\link{fit_maihda}}. The weight
#'   column joins the complete-case filter so every step uses the same analytic
#'   sample.
#'
#' @return A data.frame (class \code{maihda_stepwise}) showing the sequential
#'   models, the between-stratum variance at each step, and both the step-specific
#'   and total PCV. For a \strong{binary} (binomial/Bernoulli) outcome it also carries
#'   the discriminatory-accuracy trajectory: \code{AUC} (the C-statistic of each
#'   step's model -- step 0 is the strata-only discriminatory accuracy),
#'   \code{Step_AUC} and \code{Total_AUC} (the \emph{absolute} change in AUC,
#'   delta-AUC, versus the previous step and versus the null), and \code{MOR} (the
#'   Median Odds Ratio, logit link only). These columns are absent for non-binary
#'   outcomes. When \code{context} is supplied, a \code{Context_Variance} column
#'   reports the between-context variance held at each step, and the
#'   discriminatory-accuracy trajectory is omitted even for a binary outcome -- the
#'   AUC would include the context random effect and so mismatch the net-of-context
#'   \code{Step_PCV} / \code{Total_PCV}, exactly as \code{\link{summary.maihda_model}}
#'   drops it for a contextual fit.
#'
#' @details
#' All models are fit on the complete cases for `outcome`, `stratum`, and all
#' variables in `vars` so that each sequential variance comparison uses the same
#' analytic sample.
#'
#' For a binary (or ordinal) outcome the sequential between-stratum variances
#' live on the latent scale, whose level-1 variance is fixed by the link: a step
#' that adds an \emph{individual-level} variable (one varying within strata)
#' rescales the latent metric itself, so its \code{Step_PCV} mixes explained
#' variance with rescaling and can be understated or negative on that account --
#' see the \dQuote{Latent-scale families and rescaling} note in
#' \code{\link{calculate_pcv}}. Steps that add a stratum-constant dimension are
#' largely unaffected.
#'
#' For a binary outcome the table additionally tracks discriminatory accuracy
#' (Merlo et al. 2016): \code{AUC} is each model's C-statistic and \code{Step_AUC} /
#' \code{Total_AUC} are its \emph{absolute} change (delta-AUC), in contrast to the
#' \emph{proportional} \code{Step_PCV} / \code{Total_PCV}. The \code{MOR} is reported
#' for the logit link (\code{NA} otherwise) and is a monotone transform of the
#' between-stratum variance already in \code{Variance}. For a design-weighted fit
#' (\code{sampling_weights}) the AUC is the design-weighted (population) C-statistic.
#' Reuses \code{\link{maihda_discriminatory_accuracy}} on each step's fitted model, so
#' no additional models are fit. Note that adding a \emph{stratum-defining} dimension
#' (one already encoded by the strata) typically leaves the AUC essentially unchanged:
#' it re-partitions the between-stratum variance (so the PCV and MOR move) but not the
#' per-stratum predicted ranking the rank-based AUC depends on. The AUC trajectory is
#' therefore most informative for individual-level covariates that vary \emph{within}
#' strata.
#'
#' @references
#' Merlo, J., Wagner, P., Ghith, N., & Leckie, G. (2016). An original stepwise
#' multilevel logistic regression analysis of discriminatory accuracy: the case of
#' neighbourhoods and health. \emph{PLOS ONE}, 11(4), e0153778.
#'
#' @examples
#' \donttest{
#' strata_result <- make_strata(maihda_sim_data, c("gender", "race"))
#' stepwise_pcv(strata_result$data, "health_outcome", c("gender", "race", "age"))
#'
#' # Contextual cross-classified stepwise PCV: strata crossed with a higher-level
#' # context (country). Step_PCV / Total_PCV are then net of the country intercept,
#' # and a Context_Variance column reports the between-country variance per step.
#' cc <- make_strata(maihda_country_data, c("gender", "ses"))
#' stepwise_pcv(cc$data, "math", c("gender", "ses"), context = "country")
#' }
#'
#' @export
stepwise_pcv <- function(data, outcome, vars, engine = "lme4", family = "gaussian",
                         context = NULL, sampling_weights = NULL) {

  if (!"stratum" %in% names(data)) {
    stop("Variable 'stratum' not found in data. Please run make_strata() first.")
  }
  # Validate the higher-level context column(s) up front, mirroring fit_maihda()
  # (which each step calls). The clash-with-stratum-dimension and
  # clash-with-fixed-part checks fire per-step inside fit_maihda(); here we only
  # normalise and confirm existence so the context columns can join the shared
  # complete-case filter below (and every step then holds the (1 | context) effect).
  context <- maihda_validate_context(context, data)
  # Sampling weights select the design-weighted engine, mirroring fit_maihda();
  # the weight column joins the complete-case filter below so all steps share one
  # analytic sample.
  if (!is.null(sampling_weights)) {
    sampling_weights <- maihda_validate_sampling_weights(sampling_weights, data)
    if (missing(engine)) {
      engine <- "wemix"
      message("stepwise_pcv(): 'sampling_weights' supplied; using engine = \"wemix\" ",
              "(design-weighted pseudo-maximum-likelihood via WeMix).")
    } else if (identical(engine, "lme4")) {
      stop("Sampling weights are not supported by engine = \"lme4\" (lme4's ",
           "weights are precision weights, not sampling weights). Use ",
           "engine = \"wemix\" or \"brms\".", call. = FALSE)
    }
  }
  required_vars <- unique(c(outcome, "stratum", vars, context, sampling_weights))
  missing_vars <- setdiff(required_vars, names(data))
  if (length(missing_vars) > 0) {
    stop("Variables not found in data: ", paste(missing_vars, collapse = ", "))
  }

  strata_info <- attr(data, "strata_info")
  strata_vars <- attr(data, "strata_vars")
  strata_sep <- attr(data, "strata_sep")
  strata_autobin_info <- attr(data, "strata_autobin_info")

  complete_idx <- stats::complete.cases(data[, required_vars, drop = FALSE])
  if (!any(complete_idx)) {
    stop("No complete cases remain after filtering outcome, stratum, and stepwise variables.")
  }
  if (!all(complete_idx)) {
    data <- data[complete_idx, , drop = FALSE]
    attr(data, "strata_info") <- strata_info
    attr(data, "strata_vars") <- strata_vars
    attr(data, "strata_sep") <- strata_sep
    attr(data, "strata_autobin_info") <- strata_autobin_info
  }

  # Reconstruct auto-binned stratum dimensions so each step enters the SAME tertile
  # factor that defines the strata (matching the adjusted model and core maihda()),
  # rather than a raw linear term for an auto-binned numeric dimension. The displayed
  # "Added_Variable" keeps the original variable name; only the model term changes.
  model_terms <- vars
  if (length(vars) > 0 && !is.null(strata_autobin_info) &&
      length(strata_autobin_info) > 0 && any(vars %in% names(strata_autobin_info))) {
    adj <- maihda_adjusted_terms(vars, strata_autobin_info, data)
    data <- adj$data
    attr(data, "strata_info") <- strata_info
    attr(data, "strata_vars") <- strata_vars
    attr(data, "strata_sep") <- strata_sep
    attr(data, "strata_autobin_info") <- strata_autobin_info
    model_terms <- adj$terms
  }

  # Auto-detect a binary outcome when family is left at the default, mirroring
  # fit_maihda()/maihda(). Otherwise a binary outcome would silently be fit on the
  # Gaussian (linear) scale for numeric 0/1, or error for a factor. An ordered
  # factor likewise selects the cumulative (ordinal) model.
  if (missing(family) && maihda_is_binary_vector(data[[outcome]])) {
    warning("The outcome variable appears to be binary. Using family = 'binomial' ",
            "for the stepwise PCV. Specify family = 'gaussian' explicitly for a ",
            "linear probability model.", call. = FALSE)
    family <- "binomial"
  } else if (missing(family) && is.ordered(data[[outcome]]) &&
             nlevels(droplevels(data[[outcome]])) >= 3) {
    warning("The outcome variable is an ordered factor. Using the cumulative ",
            "(ordinal) model, family = 'ordinal', for the stepwise PCV.",
            call. = FALSE)
    family <- "ordinal"
  }

  # Ordinal family <-> engine handshake, mirroring fit_maihda(): the per-step
  # fits receive 'engine' explicitly, so fit_maihda()'s own auto-switch could
  # never fire through them.
  if (maihda_family_is_ordinal(
    if (is.function(family)) tryCatch(family(), error = function(e) NULL) else family
  )) {
    if (missing(engine) && is.null(sampling_weights)) {
      engine <- "ordinal"
      message("stepwise_pcv(): ordinal (cumulative) family; using engine = ",
              "\"ordinal\" (ordinal::clmm).")
    } else if (identical(engine, "lme4")) {
      stop("lme4 cannot fit a cumulative (ordinal) model. Use engine = ",
           "\"ordinal\" (ordinal::clmm, the default for this family) or ",
           "engine = \"brms\" (brms::cumulative).", call. = FALSE)
    }
  }

  # Context is a crossed random effect, which the design-weighted (wemix) and
  # cumulative (ordinal) engines cannot fit -- mirroring fit_maihda()'s guards
  # (which each step would otherwise hit at the null fit). Checked on the RESOLVED
  # engine, so a valid context + sampling_weights + engine = "brms" still passes.
  if (!is.null(context) && engine %in% c("wemix", "ordinal")) {
    stop("stepwise_pcv(): engine = \"", engine, "\" does not support 'context' ",
         "(it fits no crossed random effects). Use engine = \"lme4\" or ",
         "\"brms\" for a contextual cross-classified stepwise PCV.", call. = FALSE)
  }

  results <- data.frame(
    Step = integer(length(vars) + 1),
    Model = character(length(vars) + 1),
    Added_Variable = character(length(vars) + 1),
    Variance = numeric(length(vars) + 1),
    Step_PCV = numeric(length(vars) + 1),
    Total_PCV = numeric(length(vars) + 1),
    stringsAsFactors = FALSE
  )

  # Model 0: Null Model. Each step compares the stratum variance across models that
  # differ in fixed effects, so refit any REML lmer fit with ML first (see
  # calculate_pcv()); a no-op for glmer/wemix/ordinal binary fits used for the DA
  # trajectory below.
  null_fmla <- maihda_formula_with_stratum(outcome)
  null_mod <- maihda_pcv_refit_ml(fit_maihda(null_fmla, data, engine = engine,
                         family = family, context = context,
                         sampling_weights = sampling_weights))
  null_var <- extract_between_variance(null_mod)

  # Discriminatory-accuracy trajectory (binary outcomes only): read the AUC and MOR
  # off each step's already-fitted model, so no extra fits are needed. Whether it
  # applies is decided once from the null fit's resolved family; the vectors stay NA
  # (and the columns are dropped) for any other family, so the non-binomial table is
  # byte-for-byte unchanged. maihda_discriminatory_accuracy() already handles
  # aggregated-binomial, design-weighted (wemix) AUC and the logit-only MOR gate, so a
  # design-weighted stepwise yields the design-weighted AUC here too.
  #
  # Suppressed for a CONTEXTUAL fit (context = ), mirroring summary.maihda_model():
  # the AUC is built from full predictions that INCLUDE the context random effect, so
  # it would not match the net-of-context (stratum-only) Step_PCV/Total_PCV columns
  # here. A contextual binary stepwise therefore reports the PCV trajectory only.
  da_applies <- isTRUE(maihda_model_family_name(null_mod) %in%
                       c("binomial", "bernoulli")) && is.null(context)
  auc_vec <- rep(NA_real_, length(vars) + 1)
  mor_vec <- rep(NA_real_, length(vars) + 1)
  read_da <- function(m) {
    da <- tryCatch(maihda_discriminatory_accuracy(m), error = function(e) NULL)
    if (is.null(da)) c(NA_real_, NA_real_) else c(da$auc, da$mor)
  }
  if (da_applies) {
    da0 <- read_da(null_mod)
    auc_vec[1] <- da0[1]
    mor_vec[1] <- da0[2]
  }

  # Between-context variance trajectory (only when context = is supplied): the
  # crossed contextual random intercept held in every model, read from the same
  # (ML-refit, for lme4) fit used for the stratum variance. Surfaced as the
  # Context_Variance column so the context share is visible alongside the
  # net-of-context stratum PCV. NA (and the column dropped) otherwise, so the
  # non-contextual table is byte-for-byte unchanged.
  ctx_var_vec <- rep(NA_real_, length(vars) + 1)
  if (!is.null(context)) {
    ctx_var_vec[1] <- maihda_stepwise_context_variance(null_mod, context)
  }

  results[1, ] <- list(
    Step = 0,
    Model = "Null Model",
    Added_Variable = "None (Intercept only)",
    Variance = null_var,
    Step_PCV = 0,
    Total_PCV = 0
  )

  prev_var <- null_var

  # Sequentially add variables (using the reconstructed model term for any auto-binned
  # dimension, while reporting the original variable name in the table).
  current_terms <- character(0)

  for (i in seq_along(vars)) {
    current_terms <- c(current_terms, model_terms[i])

    fmla <- maihda_formula_with_stratum(outcome, current_terms)
    mod <- maihda_pcv_refit_ml(fit_maihda(fmla, data, engine = engine,
                      family = family, context = context,
                      sampling_weights = sampling_weights))

    curr_var <- extract_between_variance(mod)

    if (da_applies) {
      dai <- read_da(mod)
      auc_vec[i + 1] <- dai[1]
      mor_vec[i + 1] <- dai[2]
    }

    if (!is.null(context)) {
      ctx_var_vec[i + 1] <- maihda_stepwise_context_variance(mod, context)
    }

    step_pcv <- if (prev_var > 0) (prev_var - curr_var) / prev_var else NA
    total_pcv <- if (null_var > 0) (null_var - curr_var) / null_var else NA

    results[i + 1, ] <- list(
      Step = i,
      Model = sprintf("Model %d", i),
      Added_Variable = vars[i],
      Variance = curr_var,
      Step_PCV = step_pcv,
      Total_PCV = total_pcv
    )

    prev_var <- curr_var
  }

  # Attach the discriminatory-accuracy trajectory as extra columns only for a binary
  # outcome (so the gaussian/poisson/ordinal table is unchanged). Step_AUC / Total_AUC
  # are ABSOLUTE changes in AUC (delta-AUC) -- versus the previous step and versus the
  # null -- unlike the PROPORTIONAL Step_PCV / Total_PCV; the null row anchors both at 0.
  if (da_applies) {
    results$AUC <- auc_vec
    results$Step_AUC <- c(0, diff(auc_vec))
    results$Total_AUC <- auc_vec - auc_vec[1]
    results$MOR <- mor_vec
  }

  # Attach the between-context variance held in every model as an extra column, only
  # when context = was supplied (so the non-contextual table is unchanged). Placed
  # next to the net-of-context stratum Variance so the partition reads together.
  # A contextual fit never carries the AUC columns above (da_applies is FALSE when
  # context is set), so the column set here is exactly the PCV trajectory plus this.
  if (!is.null(context)) {
    results$Context_Variance <- ctx_var_vec
    base_cols <- setdiff(names(results), "Context_Variance")
    results <- results[append(base_cols, "Context_Variance",
                              after = match("Variance", base_cols))]
  }

  class(results) <- c("maihda_stepwise", "data.frame")
  return(results)
}

#' Print a stepwise MAIHDA table
#'
#' @param x A \code{maihda_stepwise} object from \code{\link{stepwise_pcv}}.
#' @param ... Additional arguments (not used).
#' @return Invisibly, \code{x}.
#' @export
print.maihda_stepwise <- function(x, ...) {
  print(as.data.frame(x), row.names = FALSE, digits = 4)
  if (all(c("AUC", "Step_AUC", "Total_AUC") %in% names(x))) {
    cat(maihda_palette()$muted(paste0(
        "\nStep_PCV / Total_PCV are proportional changes in between-stratum variance;\n",
        "Step_AUC / Total_AUC are absolute changes in AUC (delta-AUC). MOR is the\n",
        "median odds ratio (logit link only).\n")))
  }
  if ("Context_Variance" %in% names(x)) {
    cat(maihda_palette()$muted(paste0(
        "\nStep_PCV / Total_PCV are the between-stratum PCV NET OF the context random\n",
        "intercept held in every model; Context_Variance is that (summed)\n",
        "between-context variance at each step.\n")))
  }
  invisible(x)
}
