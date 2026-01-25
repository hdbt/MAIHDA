#' Compare MAIHDA Models
#'
#' Compares variance partition coefficients (VPC/ICC) across multiple MAIHDA models,
#' with optional bootstrap confidence intervals.
#'
#' @param ... Multiple maihda_model objects to compare.
#' @param model_names Optional character vector of names for the models.
#' @param bootstrap Logical indicating whether to compute bootstrap confidence
#'   intervals. Default is FALSE.
#' @param n_boot Number of bootstrap samples if bootstrap = TRUE. Default is 1000.
#' @param conf_level Confidence level for bootstrap intervals. Default is 0.95.
#'
#' @return A data frame comparing VPC/ICC across models with optional confidence intervals.
#'
#' @examples
#' \dontrun{
#' model1 <- fit_maihda(outcome ~ age + (1 | stratum), data = data1)
#' model2 <- fit_maihda(outcome ~ age + gender + (1 | stratum), data = data2)
#' 
#' # Compare without bootstrap
#' comparison <- compare_maihda(model1, model2, 
#'                             model_names = c("Base", "With Gender"))
#' 
#' # Compare with bootstrap CI
#' comparison_boot <- compare_maihda(model1, model2,
#'                                  model_names = c("Base", "With Gender"),
#'                                  bootstrap = TRUE, n_boot = 500)
#' }
#'
#' @export
compare_maihda <- function(..., model_names = NULL, bootstrap = FALSE, 
                          n_boot = 1000, conf_level = 0.95) {
  models <- list(...)
  
  # Validate inputs
  if (length(models) == 0) {
    stop("At least one model must be provided")
  }
  
  for (i in seq_along(models)) {
    if (!inherits(models[[i]], "maihda_model")) {
      stop("All arguments must be maihda_model objects")
    }
  }
  
  # Create model names if not provided
  if (is.null(model_names)) {
    model_names <- paste0("Model", seq_along(models))
  } else {
    if (length(model_names) != length(models)) {
      stop("Length of model_names must match number of models")
    }
  }
  
  # Calculate VPC for each model
  comparison_list <- lapply(seq_along(models), function(i) {
    summary_obj <- summary_maihda(models[[i]], bootstrap = bootstrap, 
                                 n_boot = n_boot, conf_level = conf_level)
    
    if (bootstrap && summary_obj$vpc$bootstrap) {
      data.frame(
        model = model_names[i],
        vpc = summary_obj$vpc$estimate,
        ci_lower = summary_obj$vpc$ci_lower,
        ci_upper = summary_obj$vpc$ci_upper,
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        model = model_names[i],
        vpc = summary_obj$vpc$estimate,
        stringsAsFactors = FALSE
      )
    }
  })
  
  # Combine results
  comparison_df <- do.call(rbind, comparison_list)
  
  return(comparison_df)
}

#' Plot Model Comparison
#'
#' Creates a plot comparing VPC/ICC across multiple models.
#'
#' @param comparison_df A data frame from \code{compare_maihda()}.
#'
#' @return A ggplot2 object.
#'
#' @examples
#' \dontrun{
#' comparison <- compare_maihda(model1, model2, bootstrap = TRUE)
#' plot_comparison(comparison)
#' }
#'
#' @export
#' @import ggplot2
plot_comparison <- function(comparison_df) {
  if (!is.data.frame(comparison_df) || !"vpc" %in% names(comparison_df)) {
    stop("comparison_df must be a data frame with a 'vpc' column")
  }
  
  has_ci <- all(c("ci_lower", "ci_upper") %in% names(comparison_df))
  
  p <- ggplot(comparison_df, aes(x = model, y = vpc)) +
    geom_point(size = 4, color = "#0072B2") +
    labs(
      title = "Comparison of Variance Partition Coefficients",
      x = "Model",
      y = "VPC (ICC)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    ylim(0, 1)
  
  if (has_ci) {
    p <- p + geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), 
                          width = 0.2, color = "#0072B2")
  }
  
  return(p)
}
