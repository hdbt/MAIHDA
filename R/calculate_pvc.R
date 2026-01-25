#' Calculate Proportional Change in Between-Stratum Variance (PVC)
#'
#' Calculates the proportional change in between-stratum variance (PVC) between 
#' two MAIHDA models. The PVC measures how much the between-stratum variance 
#' changes when moving from one model to another, and is calculated as:
#' PVC = (Var_model1 - Var_model2) / Var_model1
#'
#' @param model1 A maihda_model object from \code{fit_maihda()}. This is the 
#'   reference model (typically a simpler or baseline model).
#' @param model2 A maihda_model object from \code{fit_maihda()}. This is the 
#'   comparison model (typically a more complex model with additional predictors).
#' @param bootstrap Logical indicating whether to compute bootstrap confidence
#'   intervals for PVC. Default is FALSE.
#' @param n_boot Number of bootstrap samples if bootstrap = TRUE. Default is 1000.
#' @param conf_level Confidence level for bootstrap intervals. Default is 0.95.
#'
#' @return A list containing:
#'   \item{pvc}{The estimated proportional change in variance}
#'   \item{var_model1}{Between-stratum variance from model1}
#'   \item{var_model2}{Between-stratum variance from model2}
#'   \item{ci_lower}{Lower bound of confidence interval (if bootstrap = TRUE)}
#'   \item{ci_upper}{Upper bound of confidence interval (if bootstrap = TRUE)}
#'   \item{bootstrap}{Logical indicating if bootstrap was used}
#'
#' @details
#' The PVC is interpreted as the proportional reduction (or increase if negative) 
#' in between-stratum variance when moving from model1 to model2. A positive PVC 
#' indicates that model2 explains some of the between-stratum variance present in 
#' model1, while a negative PVC suggests that model2 has more unexplained 
#' between-stratum variance.
#' 
#' When bootstrap = TRUE, the function resamples the data with replacement and 
#' refits both models for each bootstrap sample to obtain confidence intervals 
#' for the PVC estimate.
#'
#' @examples
#' \dontrun{
#' # Fit two models
#' model1 <- fit_maihda(outcome ~ age + (1 | stratum), data = data)
#' model2 <- fit_maihda(outcome ~ age + gender + (1 | stratum), data = data)
#' 
#' # Calculate PVC without bootstrap
#' pvc_result <- calculate_pvc(model1, model2)
#' print(pvc_result$pvc)
#' 
#' # Calculate PVC with bootstrap CI
#' pvc_boot <- calculate_pvc(model1, model2, bootstrap = TRUE, n_boot = 500)
#' print(pvc_boot)
#' }
#'
#' @export
#' @importFrom lme4 lmer glmer VarCorr
calculate_pvc <- function(model1, model2, bootstrap = FALSE, 
                         n_boot = 1000, conf_level = 0.95) {
  # Input validation
  if (!inherits(model1, "maihda_model")) {
    stop("'model1' must be a maihda_model object from fit_maihda()")
  }
  
  if (!inherits(model2, "maihda_model")) {
    stop("'model2' must be a maihda_model object from fit_maihda()")
  }
  
  if (model1$engine != model2$engine) {
    stop("Both models must use the same engine (lme4 or brms)")
  }
  
  # Check that models were fit on compatible data
  if (nrow(model1$data) != nrow(model2$data)) {
    warning("Models were fit on data with different numbers of observations. Results may not be meaningful.")
  }
  
  # Extract between-stratum variance from both models
  var1 <- extract_between_variance(model1)
  var2 <- extract_between_variance(model2)
  
  # Validate variances
  if (is.na(var1) || is.na(var2)) {
    stop("Unable to extract variance components from one or both models")
  }
  
  if (var1 <= 0) {
    stop("Between-stratum variance in model1 is zero or negative. PVC cannot be calculated. ",
         "This may indicate a singular fit or no between-stratum variation.")
  }
  
  # Calculate PVC
  pvc <- (var1 - var2) / var1
  
  # Create result object
  result <- list(
    pvc = pvc,
    var_model1 = var1,
    var_model2 = var2,
    bootstrap = FALSE
  )
  
  # Bootstrap confidence intervals if requested
  if (bootstrap) {
    if (!requireNamespace("boot", quietly = TRUE)) {
      warning("Package 'boot' is suggested but not installed. Computing bootstrap without boot package.")
    }
    
    pvc_ci <- bootstrap_pvc(model1, model2, n_boot, conf_level)
    result$ci_lower <- pvc_ci[1]
    result$ci_upper <- pvc_ci[2]
    result$bootstrap <- TRUE
    result$conf_level <- conf_level
  }
  
  class(result) <- "pvc_result"
  return(result)
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
  
  if (engine == "lme4") {
    # Extract variance components
    vc <- lme4::VarCorr(fitted_model)
    var_random <- as.numeric(vc[[1]][1])  # Between-stratum variance
    return(var_random)
    
  } else if (engine == "brms") {
    # Verify brms is available
    if (!requireNamespace("brms", quietly = TRUE)) {
      stop("Package 'brms' is required to work with brms models. Please install it with: install.packages('brms')")
    }
    
    # Extract variance components from brms model
    vc <- brms::VarCorr(fitted_model)
    var_random <- vc[[1]]$sd[1, "Estimate"]^2
    return(var_random)
    
  } else {
    stop("Unsupported engine: ", engine, ". Only 'lme4' and 'brms' are supported.")
  }
}

#' Bootstrap PVC
#'
#' Internal function to compute bootstrap confidence intervals for PVC.
#'
#' @param model1 First maihda_model object
#' @param model2 Second maihda_model object
#' @param n_boot Number of bootstrap samples
#' @param conf_level Confidence level
#'
#' @return A vector with lower and upper confidence bounds
#' @keywords internal
#' @importFrom lme4 lmer glmer VarCorr
bootstrap_pvc <- function(model1, model2, n_boot, conf_level) {
  pvc_boot <- numeric(n_boot)
  
  # Use the larger dataset for bootstrap (in case they differ)
  data1 <- model1$data
  data2 <- model2$data
  n1 <- nrow(data1)
  n2 <- nrow(data2)
  
  # Check if data frames are compatible for bootstrap
  # We need to ensure we can resample indices that work for both
  if (n1 != n2) {
    warning("Models fitted on different sized datasets. Using model1's data size for bootstrap.")
    n <- n1
    # Only use model1's data for bootstrap
    use_data2 <- FALSE
  } else {
    n <- n1
    use_data2 <- TRUE
  }
  
  formula1 <- model1$formula
  formula2 <- model2$formula
  engine <- model1$engine
  
  for (i in 1:n_boot) {
    # Resample with replacement
    boot_indices <- sample(1:n, n, replace = TRUE)
    boot_data1 <- data1[boot_indices, ]
    
    if (use_data2) {
      boot_data2 <- data2[boot_indices, ]
    } else {
      boot_data2 <- boot_data1  # Fallback if different sizes
    }
    
    # Fit both models on bootstrap sample
    tryCatch({
      if (engine == "lme4") {
        if (inherits(model1$model, "lmerMod")) {
          boot_model1 <- lme4::lmer(formula1, data = boot_data1)
        } else {
          boot_model1 <- lme4::glmer(formula1, data = boot_data1, 
                                     family = model1$family)
        }
        
        if (inherits(model2$model, "lmerMod")) {
          boot_model2 <- lme4::lmer(formula2, data = boot_data2)
        } else {
          boot_model2 <- lme4::glmer(formula2, data = boot_data2,
                                     family = model2$family)
        }
        
        # Extract variances
        vc1 <- lme4::VarCorr(boot_model1)
        var1 <- as.numeric(vc1[[1]][1])
        
        vc2 <- lme4::VarCorr(boot_model2)
        var2 <- as.numeric(vc2[[1]][1])
        
        # Calculate PVC
        pvc_boot[i] <- (var1 - var2) / var1
        
      } else if (engine == "brms") {
        # For brms, this would require refitting which is computationally expensive
        # Not implementing bootstrap for brms in this version
        stop("Bootstrap is not yet implemented for brms models")
      }
    }, error = function(e) {
      pvc_boot[i] <- NA
    })
  }
  
  # Remove NAs
  pvc_boot <- pvc_boot[!is.na(pvc_boot)]
  
  if (length(pvc_boot) < n_boot * 0.5) {
    warning(sprintf("More than 50%% of bootstrap samples failed. CI may be unreliable. Only %d/%d successful.", 
                    length(pvc_boot), n_boot))
  }
  
  # Calculate confidence interval
  alpha <- 1 - conf_level
  ci <- quantile(pvc_boot, probs = c(alpha/2, 1 - alpha/2))
  
  return(ci)
}

#' Print method for pvc_result objects
#'
#' @param x A pvc_result object
#' @param ... Additional arguments (not used)
#' @export
print.pvc_result <- function(x, ...) {
  cat("Proportional Change in Variance (PVC)\n")
  cat("=====================================\n\n")
  
  if (x$bootstrap) {
    conf_pct <- if (!is.null(x$conf_level)) x$conf_level * 100 else 95
    cat(sprintf("PVC: %.4f [%.4f, %.4f]\n", 
                x$pvc, x$ci_lower, x$ci_upper))
    cat(sprintf("(Bootstrap %.0f%% CI)\n\n", conf_pct))
  } else {
    cat(sprintf("PVC: %.4f\n\n", x$pvc))
  }
  
  cat("Between-stratum variance:\n")
  cat(sprintf("  Model 1: %.6f\n", x$var_model1))
  cat(sprintf("  Model 2: %.6f\n", x$var_model2))
  cat(sprintf("  Change:  %.6f (%.2f%%)\n", 
              x$var_model1 - x$var_model2,
              x$pvc * 100))
  
  cat("\nInterpretation:\n")
  if (x$pvc > 0) {
    cat(sprintf("  Model 2 explains %.1f%% of the between-stratum variance\n", x$pvc * 100))
    cat("  present in Model 1 (variance reduction).\n")
  } else if (x$pvc < 0) {
    cat(sprintf("  Model 2 has %.1f%% more between-stratum variance than\n", abs(x$pvc) * 100))
    cat("  Model 1 (variance increase).\n")
  } else {
    cat("  No change in between-stratum variance between models.\n")
  }
  
  invisible(x)
}
