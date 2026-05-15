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
#' @param object A maihda_model object from \code{fit_maihda()}.
#' @param bootstrap Logical indicating whether to compute bootstrap confidence
#'   intervals for VPC/ICC. Default is FALSE.
#' @param n_boot Number of bootstrap samples if bootstrap = TRUE. Default is 1000.
#' @param conf_level Confidence level for bootstrap intervals. Default is 0.95.
#' @param ... Additional arguments (not currently used).
#'
#' @return A maihda_summary object containing:
#'   \item{vpc}{Variance Partition Coefficient (ICC) with optional CI}
#'   \item{variance_components}{Data frame of variance components}
#'   \item{stratum_estimates}{Data frame of stratum-specific random effects with labels if available}
#'   \item{fixed_effects}{Fixed effects estimates}
#'   \item{model_summary}{Original model summary}
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

  engine <- object$engine
  model <- object$model

  # Extract variance components and calculate VPC
  if (engine == "lme4") {
    # Extract variance components
    vc <- lme4::VarCorr(model)
    var_random <- maihda_stratum_variance_lme4(model)
    var_residual <- maihda_residual_variance_lme4(model, vc)

    # Calculate VPC (ICC)
    vpc <- var_random / (var_random + var_residual)

    # Create variance components data frame
    variance_components <- data.frame(
      component = c("Between-stratum (random)", "Within-stratum (residual)", "Total"),
      variance = c(var_random, var_residual, var_random + var_residual),
      sd = c(sqrt(var_random), sqrt(var_residual), sqrt(var_random + var_residual)),
      proportion = c(vpc, 1 - vpc, 1.0)
    )

    # Bootstrap confidence intervals for VPC if requested
    if (bootstrap) {
      vpc_ci <- bootstrap_vpc(model, object$data, object$formula, n_boot, conf_level)
      vpc_result <- list(
        estimate = vpc,
        ci_lower = vpc_ci[1],
        ci_upper = vpc_ci[2],
        bootstrap = TRUE
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
    # Verify brms is available
    if (!requireNamespace("brms", quietly = TRUE)) {
      stop("Package 'brms' is required to summarize brms models. Please install it with: install.packages('brms')")
    }

    var_random <- maihda_stratum_variance_brms(model)
    var_residual <- maihda_residual_variance_brms(model)

    # Calculate VPC
    vpc <- var_random / (var_random + var_residual)

    variance_components <- data.frame(
      component = c("Between-stratum (random)", "Within-stratum (residual)", "Total"),
      variance = c(var_random, var_residual, var_random + var_residual),
      sd = c(sqrt(var_random), sqrt(var_residual), sqrt(var_random + var_residual)),
      proportion = c(vpc, 1 - vpc, 1.0)
    )

    # For brms, bootstrap is not implemented the same way
    vpc_result <- list(
      estimate = vpc,
      bootstrap = FALSE
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
  vpc_boot <- numeric(n_boot)
  sim_data <- stats::simulate(model, nsim = n_boot)

  for (i in 1:n_boot) {
    tryCatch({
      boot_model <- lme4::refit(model, newresp = sim_data[[i]])

      # Calculate VPC
      vc <- lme4::VarCorr(boot_model)
      var_random <- maihda_stratum_variance_lme4(boot_model)
      var_residual <- maihda_residual_variance_lme4(boot_model, vc)

      vpc_boot[i] <- var_random / (var_random + var_residual)
    }, error = function(e) {
      vpc_boot[i] <- NA
    })
  }

  # Remove NAs
  vpc_boot <- vpc_boot[is.finite(vpc_boot)]
  if (length(vpc_boot) == 0) {
    stop("All VPC bootstrap refits failed.")
  }

  # Calculate confidence interval
  alpha <- 1 - conf_level
  ci <- stats::quantile(vpc_boot, probs = c(alpha/2, 1 - alpha/2))

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
  if (x$vpc$bootstrap) {
    cat(sprintf("  Estimate: %.4f [%.4f, %.4f]\n",
                x$vpc$estimate, x$vpc$ci_lower, x$vpc$ci_upper))
    cat("  (Bootstrap 95% CI)\n\n")
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
