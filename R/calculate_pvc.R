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
#' When bootstrap = TRUE, the function uses a parametric bootstrap: it simulates
#' new responses from model2 and refits both models with \code{lme4::refit()} for
#' each simulated response to obtain confidence intervals for the PVC estimate.
#'
#' @examples
#' \donttest{
#' # Create strata and fit two models
#' strata_result <- make_strata(maihda_sim_data, c("gender", "race"))
#' model1 <- fit_maihda(health_outcome ~ age + (1 | stratum), data = strata_result$data)
#' model2 <- fit_maihda(health_outcome ~ age + gender + (1 | stratum), data = strata_result$data)
#'
#' # Calculate PVC without bootstrap
#' pvc_result <- calculate_pvc(model1, model2)
#' print(pvc_result$pvc)
#'
#' # Calculate PVC with bootstrap CI
#' # pvc_boot <- calculate_pvc(model1, model2, bootstrap = TRUE, n_boot = 500)
#' # print(pvc_boot)
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

  validate_pvc_models(model1, model2)

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
    maihda_validate_intercept_only_random_effects_lme4(
      fitted_model,
      context = "PVC calculations"
    )
    return(maihda_stratum_variance_lme4(fitted_model))

  } else if (engine == "brms") {
    if (!requireNamespace("brms", quietly = TRUE)) {
      stop("Package 'brms' is required to work with brms models. Please install it with: install.packages('brms')")
    }
    maihda_validate_intercept_only_random_effects_brms(
      brms::VarCorr(fitted_model),
      context = "PVC calculations"
    )
    return(maihda_stratum_variance_brms(fitted_model))

  } else {
    stop("Unsupported engine: ", engine, ". Only 'lme4' and 'brms' are supported.")
  }
}

validate_pvc_models <- function(model1, model2) {
  response1 <- paste(deparse(model1$formula[[2]]), collapse = "")
  response2 <- paste(deparse(model2$formula[[2]]), collapse = "")
  if (!identical(response1, response2)) {
    stop("PVC requires both models to use the same outcome. ",
         "Model 1 uses '", response1, "' and Model 2 uses '", response2, "'.",
         call. = FALSE)
  }

  fam1 <- maihda_family(model1$model)
  fam2 <- maihda_family(model2$model)
  fam_key1 <- c(
    family = if (!is.null(fam1$family)) fam1$family else NA_character_,
    link = if (!is.null(fam1$link)) fam1$link else NA_character_
  )
  fam_key2 <- c(
    family = if (!is.null(fam2$family)) fam2$family else NA_character_,
    link = if (!is.null(fam2$link)) fam2$link else NA_character_
  )
  if (!identical(fam_key1, fam_key2)) {
    stop("PVC requires both models to use the same model family and link. ",
         "Model 1 uses ", fam_key1[["family"]], "(", fam_key1[["link"]], ") and ",
         "Model 2 uses ", fam_key2[["family"]], "(", fam_key2[["link"]], ").",
         call. = FALSE)
  }

  n1 <- maihda_nobs(model1$model)
  n2 <- maihda_nobs(model2$model)
  if (is.finite(n1) && is.finite(n2) && n1 != n2) {
    stop("PVC requires both models to use the same analytic sample. ",
         "Model 1 used ", n1, " observations and Model 2 used ", n2, ".",
         call. = FALSE)
  }

  rows1 <- maihda_row_ids(model1$model)
  rows2 <- maihda_row_ids(model2$model)
  if (!is.null(rows1) && !is.null(rows2) && !identical(rows1, rows2)) {
    stop("PVC requires both models to use the same analytic sample in the same row order.",
         call. = FALSE)
  }

  strata1 <- unique(as.character(model1$data$stratum))
  strata2 <- unique(as.character(model2$data$stratum))
  strata1 <- sort(strata1[!is.na(strata1)])
  strata2 <- sort(strata2[!is.na(strata2)])
  if (!identical(strata1, strata2)) {
    stop("PVC requires both models to use the same stratum definitions.",
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
      stop("PVC requires both models to use the same stratum labels.",
           call. = FALSE)
    }
  }

  invisible(TRUE)
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
  engine <- model1$engine
  if (engine != "lme4") {
    stop("Bootstrap is currently only supported for lme4 models.")
  }

  pvc_boot <- numeric(n_boot)

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

      # Calculate PVC
      pvc_boot[i] <- if (is.finite(var1) && var1 > 0) (var1 - var2) / var1 else NA_real_
    }, error = function(e) {
      pvc_boot[i] <- NA
    })
  }

  # Remove NAs
  pvc_boot <- pvc_boot[is.finite(pvc_boot)]

  if (length(pvc_boot) < n_boot * 0.5) {
    warning(sprintf("More than 50%% of bootstrap samples failed. CI may be unreliable. Only %d/%d successful.",
                    length(pvc_boot), n_boot))
  }
  if (length(pvc_boot) == 0) {
    stop("All PVC bootstrap refits failed or produced zero model-1 stratum variance.")
  }

  # Calculate confidence interval
  alpha <- 1 - conf_level
  ci <- stats::quantile(pvc_boot, probs = c(alpha/2, 1 - alpha/2))

  return(ci)
}

#' Print method for PVC results
#'
#' @param x A pvc_result object
#' @param ... Additional arguments
#' @return No return value, called for side effects.
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

#' Stepwise Proportional Change in Variance (PCV)
#'
#' @description
#' Estimates the proportional change in variance (PCV) sequentially by fitting
#' intermediate (partially-adjusted) models. It adds each predictor variable
#' one-by-one to gauge its unique contribution in explaining between-stratum
#' inequalities.
#'
#' @param data Data frame with observations. Ensure `make_strata()` was run first
#'   so the `stratum` variable exists.
#' @param outcome Character string; the dependent variable.
#' @param vars Character vector; predictors (strata groupings & covariates) to
#'   add sequentially to the model.
#' @param engine Modeling engine ("lme4" or "brms"). Default is "lme4".
#' @param family Error distribution and link function. Default is "gaussian".
#'
#' @return A data.frame showing the sequential models, the between-stratum
#'   variance at each step, and both the step-specific and total PCV.
#'
#' @details
#' All models are fit on the complete cases for `outcome`, `stratum`, and all
#' variables in `vars` so that each sequential variance comparison uses the same
#' analytic sample.
#'
#' @examples
#' \donttest{
#' strata_result <- make_strata(maihda_sim_data, c("gender", "race"))
#' stepwise_pcv(strata_result$data, "health_outcome", c("gender", "race", "age"))
#' }
#' @importFrom stats as.formula
#'
#' @export
stepwise_pcv <- function(data, outcome, vars, engine = "lme4", family = "gaussian") {

  if (!"stratum" %in% names(data)) {
    stop("Variable 'stratum' not found in data. Please run make_strata() first.")
  }
  required_vars <- unique(c(outcome, "stratum", vars))
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

  results <- data.frame(
    Step = integer(length(vars) + 1),
    Model = character(length(vars) + 1),
    Added_Variable = character(length(vars) + 1),
    Variance = numeric(length(vars) + 1),
    Step_PCV = numeric(length(vars) + 1),
    Total_PCV = numeric(length(vars) + 1),
    stringsAsFactors = FALSE
  )

  # Model 0: Null Model
  null_fmla <- maihda_formula_with_stratum(outcome)
  null_mod <- fit_maihda(null_fmla, data, engine = engine, family = family)
  null_var <- extract_between_variance(null_mod)

  results[1, ] <- list(
    Step = 0,
    Model = "Null Model",
    Added_Variable = "None (Intercept only)",
    Variance = null_var,
    Step_PCV = 0,
    Total_PCV = 0
  )

  prev_var <- null_var

  # Sequentially add variables
  current_vars <- c()

  for (i in seq_along(vars)) {
    var <- vars[i]
    current_vars <- c(current_vars, var)

    fmla <- maihda_formula_with_stratum(outcome, current_vars)
    mod <- fit_maihda(fmla, data, engine = engine, family = family)

    curr_var <- extract_between_variance(mod)

    step_pcv <- if (prev_var > 0) (prev_var - curr_var) / prev_var else NA
    total_pcv <- if (null_var > 0) (null_var - curr_var) / null_var else NA

    results[i + 1, ] <- list(
      Step = i,
      Model = sprintf("Model %d", i),
      Added_Variable = var,
      Variance = curr_var,
      Step_PCV = step_pcv,
      Total_PCV = total_pcv
    )

    prev_var <- curr_var
  }

  class(results) <- c("maihda_stepwise", "data.frame")
  return(results)
}
