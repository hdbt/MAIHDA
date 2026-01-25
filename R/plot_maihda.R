#' Plot MAIHDA Model Results
#'
#' Creates various plots for visualizing MAIHDA model results including caterpillar
#' plots, variance partition coefficient comparisons, observed vs. shrunken estimates,
#' and predicted subgroup values with confidence intervals.
#'
#' @param object A maihda_model object from \code{fit_maihda()}.
#' @param type Character string specifying plot type:
#'   \itemize{
#'     \item "caterpillar": Caterpillar plot of stratum random effects
#'     \item "vpc": Variance partition coefficient visualization
#'     \item "obs_vs_shrunken": Observed vs. shrunken stratum means
#'     \item "predicted": Predicted values for each stratum with confidence intervals
#'   }
#' @param summary_obj Optional maihda_summary object from \code{summary_maihda()}.
#'   If NULL, will be computed.
#' @param n_strata Maximum number of strata to display in caterpillar plot or predicted plot.
#'   Default is 50. Use NULL for all strata.
#' @param ... Additional arguments (not currently used).
#'
#' @return A ggplot2 object.
#'
#' @examples
#' \dontrun{
#' model <- fit_maihda(outcome ~ age + (1 | stratum), data = data)
#'
#' # Caterpillar plot
#' plot_maihda(model, type = "caterpillar")
#'
#' # VPC plot
#' plot_maihda(model, type = "vpc")
#'
#' # Observed vs shrunken plot
#' plot_maihda(model, type = "obs_vs_shrunken")
#'
#' # Predicted values with confidence intervals
#' plot_maihda(model, type = "predicted")
#' }
#'
#' @export
#' @import ggplot2
#' @importFrom dplyr arrange
plot_maihda <- function(object, type = c("caterpillar", "vpc", "obs_vs_shrunken", "predicted"),
                       summary_obj = NULL, n_strata = 50, ...) {
  if (!inherits(object, "maihda_model")) {
    stop("'object' must be a maihda_model object from fit_maihda()")
  }

  type <- match.arg(type)

  # Get summary if not provided
  if (is.null(summary_obj)) {
    summary_obj <- summary_maihda(object)
  }

  if (type == "caterpillar") {
    plot <- plot_caterpillar(summary_obj, n_strata)
  } else if (type == "vpc") {
    plot <- plot_vpc(summary_obj)
  } else if (type == "obs_vs_shrunken") {
    plot <- plot_obs_vs_shrunken(object, summary_obj)
  } else if (type == "predicted") {
    plot <- plot_predicted_strata(object, summary_obj, n_strata)
  }

  return(plot)
}

#' Caterpillar Plot
#'
#' @param summary_obj A maihda_summary object
#' @param n_strata Maximum number of strata to display
#' @return A ggplot2 object
#' @keywords internal
#' @import ggplot2
#' @importFrom dplyr arrange slice
plot_caterpillar <- function(summary_obj, n_strata) {
  stratum_est <- summary_obj$stratum_estimates

  if (is.null(stratum_est) || nrow(stratum_est) == 0) {
    stop("No stratum estimates available for plotting")
  }

  # Sort by random effect
  stratum_est <- dplyr::arrange(stratum_est, random_effect)

  # Limit number of strata if requested
  if (!is.null(n_strata) && nrow(stratum_est) > n_strata) {
    indices <- as.integer(seq(1, nrow(stratum_est), length.out = n_strata))
    stratum_est <- dplyr::slice(stratum_est, indices)
  }

  # Create rank variable for plotting
  stratum_est$rank <- 1:nrow(stratum_est)

  # Create plot
  p <- ggplot(stratum_est, aes(x = rank, y = random_effect)) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = lower_95, ymax = upper_95), width = 0.2, alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    labs(
      title = "Caterpillar Plot of Stratum Random Effects",
      x = "Stratum Rank",
      y = "Random Effect"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.minor = element_blank()
    )

  return(p)
}

#' VPC Visualization Plot
#'
#' @param summary_obj A maihda_summary object
#' @return A ggplot2 object
#' @keywords internal
#' @import ggplot2
plot_vpc <- function(summary_obj) {
  vpc_data <- summary_obj$variance_components[1:2, ]

  # Create plot
  p <- ggplot(vpc_data, aes(x = "", y = proportion, fill = component)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_flip() +
    scale_fill_manual(values = c("Between-stratum (random)" = "#E69F00",
                                  "Within-stratum (residual)" = "#56B4E9")) +
    labs(
      title = sprintf("Variance Partition Coefficient (VPC/ICC) = %.3f",
                     summary_obj$vpc$estimate),
      x = "",
      y = "Proportion of Variance",
      fill = "Component"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank()
    ) +
    geom_text(aes(label = sprintf("%.1f%%", proportion * 100)),
              position = position_stack(vjust = 0.5),
              color = "white", fontface = "bold", size = 5)

  return(p)
}

#' Observed vs. Shrunken Estimates Plot
#'
#' @param object A maihda_model object
#' @param summary_obj A maihda_summary object
#' @return A ggplot2 object
#' @keywords internal
#' @import ggplot2
#' @importFrom dplyr group_by summarise
#' @importFrom stats formula
plot_obs_vs_shrunken <- function(object, summary_obj) {
  data <- object$data

  # Get outcome variable name from formula
  formula_obj <- object$formula
  outcome_var <- all.vars(formula_obj)[1]

  # Check if outcome and stratum exist
  if (!outcome_var %in% names(data)) {
    stop("Outcome variable not found in data")
  }
  if (!"stratum" %in% names(data)) {
    stop("'stratum' variable not found in data. Make sure to use data from make_strata()")
  }

  # Calculate observed stratum means
  obs_means <- data |>
    dplyr::group_by(stratum) |>
    dplyr::summarise(
      observed = mean(.data[[outcome_var]], na.rm = TRUE),
      n = dplyr::n(),
      .groups = "drop"
    )
  
  # Convert stratum to character for merging (to match stratum_estimates)
  obs_means$stratum <- as.character(obs_means$stratum)

  # Get fixed effects to center the random effects
  fixed_intercept <- summary_obj$fixed_effects$estimate[1]

  # Merge with random effects (shrunken estimates)
  stratum_est <- summary_obj$stratum_estimates
  if (!is.null(stratum_est)) {
    plot_data <- merge(obs_means, stratum_est, by = "stratum")
    plot_data$shrunken <- fixed_intercept + plot_data$random_effect

    # Create plot
    p <- ggplot(plot_data, aes(x = observed, y = shrunken)) +
      geom_point(aes(size = n), alpha = 0.6, color = "#0072B2") +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
      labs(
        title = "Observed vs. Shrunken Stratum Estimates",
        x = "Observed Stratum Mean",
        y = "Shrunken Estimate (with Random Effect)",
        size = "Sample Size"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right"
      )

    return(p)
  } else {
    stop("No stratum estimates available for plotting")
  }
}

#' Plot Predicted Stratum Values with Confidence Intervals
#'
#' @param object A maihda_model object
#' @param summary_obj A maihda_summary object
#' @param n_strata Maximum number of strata to display
#' @return A ggplot2 object
#' @keywords internal
#' @import ggplot2
#' @importFrom dplyr arrange slice
#' @importFrom lme4 fixef
plot_predicted_strata <- function(object, summary_obj, n_strata) {
  # Get fixed effects intercept
  if (object$engine == "lme4") {
    fixed_intercept <- lme4::fixef(object$model)[1]
  } else if (object$engine == "brms") {
    if (!requireNamespace("brms", quietly = TRUE)) {
      stop("Package 'brms' is required. Please install it with: install.packages('brms')")
    }
    fixed_intercept <- brms::fixef(object$model)[1, "Estimate"]
  } else {
    stop("Unsupported engine: ", object$engine)
  }
  
  # Get stratum estimates
  stratum_est <- summary_obj$stratum_estimates
  
  if (is.null(stratum_est) || nrow(stratum_est) == 0) {
    stop("No stratum estimates available for plotting")
  }
  
  # Calculate predicted values (fixed effect + random effect)
  stratum_est$predicted <- fixed_intercept + stratum_est$random_effect
  stratum_est$lower <- fixed_intercept + stratum_est$lower_95
  stratum_est$upper <- fixed_intercept + stratum_est$upper_95
  
  # Keep original order (no sorting)
  # Limit number of strata if requested
  if (!is.null(n_strata) && nrow(stratum_est) > n_strata) {
    indices <- as.integer(seq(1, nrow(stratum_est), length.out = n_strata))
    stratum_est <- dplyr::slice(stratum_est, indices)
  }
  
  # Use labels if available, otherwise use numeric stratum IDs
  if ("label" %in% names(stratum_est) && !all(is.na(stratum_est$label))) {
    # Use the meaningful labels for the x-axis
    stratum_est$display_label <- stratum_est$label
  } else {
    # Fall back to stratum IDs
    stratum_est$display_label <- stratum_est$stratum
  }
  
  # Create factor to preserve order for plotting
  stratum_est$display_label <- factor(stratum_est$display_label, levels = stratum_est$display_label)
  
  # Create plot
  p <- ggplot(stratum_est, aes(x = display_label, y = predicted)) +
    geom_point(size = 2, color = "#0072B2") +
    geom_errorbar(aes(ymin = lower, ymax = upper), 
                  width = 0.2, alpha = 0.5, color = "#0072B2") +
    geom_hline(yintercept = fixed_intercept, linetype = "dashed", color = "red", alpha = 0.7) +
    labs(
      title = "Predicted Subgroup Values with 95% Confidence Intervals",
      x = "Stratum",
      y = "Predicted Value",
      caption = "Dashed line represents overall mean (fixed effect)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.caption = element_text(hjust = 0.5, face = "italic", size = 9),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    )
  
  return(p)
}
