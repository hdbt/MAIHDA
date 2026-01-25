#' Add Stratum Labels to Estimates
#'
#' Internal helper function to merge stratum labels into stratum estimates.
#'
#' @param stratum_estimates Data frame with stratum estimates
#' @param strata_info Data frame with stratum information including labels
#' @return Data frame with labels merged in
#' @keywords internal
add_stratum_labels <- function(stratum_estimates, strata_info) {
  if (is.null(strata_info)) {
    return(stratum_estimates)
  }
  
  # Extract just stratum and label columns
  strata_info_subset <- strata_info[, c("stratum", "label")]
  
  # Merge labels into estimates
  stratum_estimates <- merge(stratum_estimates, strata_info_subset, 
                            by.x = "stratum_id", by.y = "stratum", 
                            all.x = TRUE, sort = FALSE)
  
  # Reorder columns to put label after stratum
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
#' \dontrun{
#' model <- fit_maihda(outcome ~ age + (1 | stratum), data = data)
#' summary_result <- summary_maihda(model)
#' 
#' # With bootstrap CI
#' summary_boot <- summary_maihda(model, bootstrap = TRUE, n_boot = 500)
#' }
#'
#' @export
#' @importFrom lme4 VarCorr fixef ranef
#' @importFrom stats vcov confint
summary_maihda <- function(object, bootstrap = FALSE, n_boot = 1000, 
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
    var_random <- as.numeric(vc[[1]][1])  # Between-stratum variance
    var_residual <- attr(vc, "sc")^2       # Within-stratum variance
    
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
    
    # Extract random effects (stratum estimates)
    re <- lme4::ranef(model, condVar = TRUE)
    stratum_re <- re[[1]]
    stratum_names <- names(re)
    
    if (length(stratum_names) > 0 && stratum_names[1] == "stratum") {
      # Get conditional variances (uncertainties)
      cond_var <- attr(re[[1]], "postVar")
      if (is.array(cond_var) && length(dim(cond_var)) == 3) {
        stratum_se <- sqrt(cond_var[1, 1, ])
      } else {
        stratum_se <- rep(NA, nrow(stratum_re))
      }
      
      stratum_estimates <- data.frame(
        stratum = rownames(stratum_re),
        stratum_id = as.integer(rownames(stratum_re)),
        random_effect = stratum_re[, 1],
        se = stratum_se,
        lower_95 = stratum_re[, 1] - 1.96 * stratum_se,
        upper_95 = stratum_re[, 1] + 1.96 * stratum_se,
        stringsAsFactors = FALSE
      )
      
      # Add stratum labels if strata_info is available
      stratum_estimates <- add_stratum_labels(stratum_estimates, object$strata_info)
    } else {
      stratum_estimates <- NULL
    }
    
    # Get model summary
    model_summary <- summary(model)
    
  } else if (engine == "brms") {
    # Verify brms is available
    if (!requireNamespace("brms", quietly = TRUE)) {
      stop("Package 'brms' is required to summarize brms models. Please install it with: install.packages('brms')")
    }
    
    # Extract variance components from brms model
    vc <- brms::VarCorr(model)
    var_random <- vc[[1]]$sd[1, "Estimate"]^2
    var_residual <- vc[[2]]$sd[1, "Estimate"]^2
    
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
    
    # Extract random effects
    ranef_result <- brms::ranef(model)[[1]]
    
    # Transform brms ranef output to match lme4 structure
    stratum_estimates <- data.frame(
      stratum = rownames(ranef_result),
      stratum_id = as.integer(rownames(ranef_result)),
      random_effect = ranef_result[, "Estimate"],
      se = ranef_result[, "Est.Error"],
      lower_95 = ranef_result[, "Q2.5"],
      upper_95 = ranef_result[, "Q97.5"],
      stringsAsFactors = FALSE
    )
    
    # Add stratum labels if strata_info is available
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
  n <- nrow(data)
  
  for (i in 1:n_boot) {
    # Resample with replacement
    boot_indices <- sample(1:n, n, replace = TRUE)
    boot_data <- data[boot_indices, ]
    
    # Fit model on bootstrap sample
    tryCatch({
      if (inherits(model, "lmerMod")) {
        boot_model <- lme4::lmer(formula, data = boot_data)
      } else {
        boot_model <- lme4::glmer(formula, data = boot_data, 
                                  family = family(model))
      }
      
      # Calculate VPC
      vc <- lme4::VarCorr(boot_model)
      var_random <- as.numeric(vc[[1]][1])
      var_residual <- attr(vc, "sc")^2
      vpc_boot[i] <- var_random / (var_random + var_residual)
    }, error = function(e) {
      vpc_boot[i] <- NA
    })
  }
  
  # Remove NAs
  vpc_boot <- vpc_boot[!is.na(vpc_boot)]
  
  # Calculate confidence interval
  alpha <- 1 - conf_level
  ci <- quantile(vpc_boot, probs = c(alpha/2, 1 - alpha/2))
  
  return(ci)
}

#' Print method for maihda_summary objects
#'
#' @param x A maihda_summary object
#' @param ... Additional arguments (not used)
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
    print(head(x$stratum_estimates, 10), row.names = FALSE, digits = 4)
    if (nrow(x$stratum_estimates) > 10) {
      cat(sprintf("  ... and %d more strata\n", nrow(x$stratum_estimates) - 10))
    }
  }
  
  invisible(x)
}
