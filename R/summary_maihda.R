#' Add Stratum Labels to Estimates
#'
#' Internal helper function to merge stratum labels into stratum estimates.
#'
#' @param stratum_estimates Data frame with stratum estimates
#' @param strata_info Data frame with stratum information including labels
#' @return Data frame with labels merged in
#' @keywords internal
add_stratum_labels <- function(stratum_estimates, strata_info) {
  if (is.null(strata_info) || !"stratum" %in% names(strata_info) || !"label" %in% names(strata_info)) {
    return(stratum_estimates)
  }

  idx <- match(as.character(stratum_estimates$stratum), as.character(strata_info$stratum))
  stratum_estimates$label <- strata_info$label[idx]

  col_order <- c("stratum", "stratum_id", "label", "random_effect", "se", "lower_95", "upper_95")
  stratum_estimates <- stratum_estimates[, col_order[col_order %in% names(stratum_estimates)]]

  return(stratum_estimates)
}

#' Summarize MAIHDA Model
#'
#' Provides a summary of a MAIHDA model including variance partition coefficients
#' (VPC/ICC) and stratum-specific estimates.
#'
#' @section Interpreting the VPC/ICC: The VPC is the between-stratum variance
#'   divided by the total \emph{unexplained} variance (between-stratum + residual);
#'   it is a conditional/residual ICC that excludes variance captured by the fixed
#'   effects, so for models with covariates it is conditional on them. It is most
#'   commonly read from the null model \code{outcome ~ 1 + (1 | stratum)}, where it
#'   is the total between-stratum share. For non-Gaussian families the level-1
#'   (residual) variance uses a latent/distributional approximation (e.g.
#'   \eqn{\pi^2/3} for logistic), so the VPC is on that latent scale. The stratum
#'   random effects represent the total between-stratum deviation; they equal the
#'   \emph{pure} intersectional (interaction) component only when the additive main
#'   effects of the strata variables are included in the model.
#'
#' @param object A maihda_model object from \code{fit_maihda()}.
#' @param bootstrap Logical indicating whether to compute parametric bootstrap
#'   confidence intervals for VPC/ICC. Default is FALSE. Supported for lme4
#'   models only; \code{brms} models always return a posterior credible interval
#'   (see Details), so \code{bootstrap = TRUE} is rejected for them.
#' @param n_boot Number of bootstrap samples if bootstrap = TRUE. Default is 1000.
#' @param conf_level Confidence level for the VPC/ICC interval -- the lme4
#'   bootstrap CI or the brms posterior credible interval. Default is 0.95.
#' @param ... Additional arguments (not currently used).
#'
#' @return A maihda_summary object containing:
#'   \item{vpc}{Variance Partition Coefficient (ICC); for lme4 with
#'     \code{bootstrap = TRUE} and for all brms models this includes
#'     \code{ci_lower}/\code{ci_upper}/\code{conf_level}}
#'   \item{variance_components}{Data frame of variance components}
#'   \item{stratum_estimates}{Data frame of stratum-specific random effects with labels if available}
#'   \item{fixed_effects}{Fixed effects estimates}
#'   \item{model_summary}{Original model summary}
#'
#' @note
#' For \code{lme4} models a VPC/ICC interval is obtained from a parametric
#' bootstrap (\code{bootstrap = TRUE}). For \code{brms} models the VPC/ICC is
#' summarised directly from the posterior draws: the reported estimate is the
#' posterior median of the per-draw VPC (\eqn{E[\sigma^2]}-based, not the biased
#' \eqn{E[\sigma]^2}) and the interval is a central credible interval at
#' \code{conf_level} (default 95\%), so no \code{bootstrap} argument is needed.
#' The variance-components table reports the posterior-mean variance components,
#' so the stratum proportion shown there may differ slightly from the headline
#' VPC because the median of a ratio is not the ratio of means. For non-Gaussian
#' \code{brms} families the level-1 (residual) variance uses the usual
#' latent-scale approximation; for \code{poisson(log)} it is evaluated at the
#' posterior-mean fitted values rather than per draw to avoid an expensive
#' \eqn{ndraws \times nobs} computation.
#'
#' @examples
#' \donttest{
#' strata_result <- make_strata(maihda_sim_data, vars = c("gender", "race"))
#' model <- fit_maihda(health_outcome ~ age + (1 | stratum), data = strata_result$data)
#' summary_result <- summary(model)
#'
#' # With bootstrap CI
#' # summary_boot <- summary(model, bootstrap = TRUE, n_boot = 50)
#' }
#'
#' @export
#' @importFrom lme4 VarCorr fixef ranef
#' @importFrom stats vcov confint
summary.maihda_model <- function(object, bootstrap = FALSE, n_boot = 1000,
                          conf_level = 0.95, ...) {
  if (!inherits(object, "maihda_model")) {
    stop("'object' must be a maihda_model object from fit_maihda()")
  }

  if (!is.logical(bootstrap) || length(bootstrap) != 1 || is.na(bootstrap)) {
    stop("'bootstrap' must be TRUE or FALSE.", call. = FALSE)
  }
  if (bootstrap) {
    bootstrap_args <- maihda_validate_bootstrap_args(n_boot, conf_level)
    n_boot <- bootstrap_args$n_boot
    conf_level <- bootstrap_args$conf_level
  }

  engine <- object$engine
  model <- object$model

  # Extract variance components and calculate VPC
  if (engine == "lme4") {
    # Extract variance components
    vc <- lme4::VarCorr(model)
    var_random <- maihda_stratum_variance_lme4(model)
    var_total_random <- maihda_total_random_variance_lme4(model)
    var_other_random <- max(0, var_total_random - var_random)
    var_residual <- maihda_residual_variance_lme4(model, vc)

    # Calculate VPC (ICC)
    vpc <- var_random / (var_random + var_other_random + var_residual)

    # Create variance components data frame
    variance_components <- maihda_variance_components_table(
      var_random, var_other_random, var_residual
    )

    # Bootstrap confidence intervals for VPC if requested
    if (bootstrap) {
      vpc_ci <- bootstrap_vpc(model, object$data, object$formula, n_boot, conf_level)
      vpc_result <- list(
        estimate = vpc,
        ci_lower = vpc_ci[1],
        ci_upper = vpc_ci[2],
        conf_level = conf_level,
        bootstrap = TRUE,
        method = "bootstrap"
      )
    } else {
      vpc_result <- list(
        estimate = vpc,
        bootstrap = FALSE
      )
    }

    # Extract fixed effects
    fixed_effects <- data.frame(
      term = names(lme4::fixef(model)),
      estimate = lme4::fixef(model),
      row.names = NULL
    )

    stratum_estimates <- maihda_stratum_ranef_lme4(model)
    stratum_estimates <- add_stratum_labels(stratum_estimates, object$strata_info)

    # Get model summary
    model_summary <- summary(model)

  } else if (engine == "brms") {
    if (bootstrap) {
      stop("Bootstrap VPC confidence intervals are only supported for lme4 models. ",
           "brms summaries already return a posterior credible interval for the ",
           "VPC/ICC, so 'bootstrap = TRUE' is not needed.",
           call. = FALSE)
    }

    # Verify brms is available
    if (!requireNamespace("brms", quietly = TRUE)) {
      stop("Package 'brms' is required to summarize brms models. Please install it with: install.packages('brms')")
    }

    conf_level <- maihda_validate_conf_level(conf_level)

    # Summarise the VPC/ICC from posterior draws (E[sd^2], with a credible
    # interval) rather than from the posterior summary SDs (E[sd]^2, no interval).
    vpc_draws <- maihda_vpc_draws_brms(model, conf_level = conf_level)

    # Components table reports the posterior-mean variance of each component.
    variance_components <- maihda_variance_components_table(
      vpc_draws$var_stratum, vpc_draws$var_other_random, vpc_draws$var_residual
    )

    vpc_result <- list(
      estimate = vpc_draws$vpc$estimate,
      ci_lower = vpc_draws$vpc$ci_lower,
      ci_upper = vpc_draws$vpc$ci_upper,
      conf_level = conf_level,
      bootstrap = FALSE,
      method = "posterior"
    )

    # Extract fixed effects
    fixed_effects <- brms::fixef(model)

    stratum_estimates <- maihda_stratum_ranef_brms(model)
    stratum_estimates <- add_stratum_labels(stratum_estimates, object$strata_info)

    model_summary <- summary(model)
  }

  # Create summary object
  result <- structure(
    list(
      vpc = vpc_result,
      variance_components = variance_components,
      stratum_estimates = stratum_estimates,
      fixed_effects = fixed_effects,
      model_summary = model_summary,
      engine = engine
    ),
    class = "maihda_summary"
  )

  return(result)
}

#' Bootstrap VPC/ICC
#'
#' Internal function to compute bootstrap confidence intervals for VPC.
#'
#' @param model An lme4 model object
#' @param data The data used to fit the model
#' @param formula The model formula
#' @param n_boot Number of bootstrap samples
#' @param conf_level Confidence level
#'
#' @return A vector with lower and upper confidence bounds
#' @keywords internal
#' @importFrom lme4 lmer glmer VarCorr
bootstrap_vpc <- function(model, data, formula, n_boot, conf_level) {
  # Initialise to NA so iterations whose refit() throws — and never reach the
  # assignment inside the tryCatch body — stay NA rather than the numeric() default of 0.
  # The error handler runs in its own scope and cannot write back to this vector,
  # so the initial value is what survives a failure.
  vpc_boot <- rep(NA_real_, n_boot)
  sim_data <- stats::simulate(model, nsim = n_boot)

  for (i in 1:n_boot) {
    tryCatch({
      boot_model <- lme4::refit(model, newresp = sim_data[[i]])

      # Calculate VPC
      vc <- lme4::VarCorr(boot_model)
      var_random <- maihda_stratum_variance_lme4(boot_model)
      var_total_random <- maihda_total_random_variance_lme4(boot_model)
      var_other_random <- max(0, var_total_random - var_random)
      var_residual <- maihda_residual_variance_lme4(boot_model, vc)

      vpc_boot[i] <- var_random / (var_random + var_other_random + var_residual)
    }, error = function(e) NULL)
  }

  # Reduce to an interval, requiring a minimum number of successful refits.
  ci <- maihda_bootstrap_ci(vpc_boot, n_boot, conf_level, "VPC")

  return(ci)
}

#' Print method for maihda_summary objects
#'
#' @param x A maihda_summary object
#' @param ... Additional arguments (not used)
#' @return No return value, called for side effects.
#' @export
print.maihda_summary <- function(x, ...) {
  cat("MAIHDA Model Summary\n")
  cat("====================\n\n")

  cat("Variance Partition Coefficient (VPC/ICC):\n")
  if (maihda_vpc_has_interval(x$vpc)) {
    cat(sprintf("  Estimate: %.4f [%.4f, %.4f]\n",
                x$vpc$estimate, x$vpc$ci_lower, x$vpc$ci_upper))
    cat("  ", maihda_vpc_interval_label(x$vpc), "\n\n", sep = "")
  } else {
    cat(sprintf("  Estimate: %.4f\n\n", x$vpc$estimate))
  }

  cat("Variance Components:\n")
  print(x$variance_components, row.names = FALSE, digits = 4)
  cat("\n")

  cat("Fixed Effects:\n")
  print(x$fixed_effects, row.names = FALSE, digits = 4)
  cat("\n")

  if (!is.null(x$stratum_estimates) && nrow(x$stratum_estimates) > 0) {
    cat("Stratum Estimates (first 10):\n")
    print(utils::head(x$stratum_estimates, 10), row.names = FALSE, digits = 4)
    if (nrow(x$stratum_estimates) > 10) {
      cat(sprintf("  ... and %d more strata\n", nrow(x$stratum_estimates) - 10))
    }
  }

  invisible(x)
}
