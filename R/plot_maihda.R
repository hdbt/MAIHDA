#' Plot MAIHDA Model Results
#'
#' Creates various plots for visualizing MAIHDA model results including
#' variance partition coefficient comparisons, observed vs. shrunken estimates,
#' and predicted subgroup values with confidence intervals.
#'
#' @param x A maihda_model object from \code{fit_maihda()}.
#' @param type Character string specifying plot type:
#'   \itemize{
#'     \item "vpc": Variance partition coefficient visualization
#'     \item "obs_vs_shrunken": Observed vs. shrunken stratum means
#'     \item "predicted": Predicted values for each stratum with confidence intervals
#'     \item "risk_vs_effect": Quadrant scatterplot comparing overall risk to intersectional effect
#'     \item "effect_decomp": Visualizes additive vs intersectional deviation from global mean
#'     \item "ternary": Ternary plot analyzing the dimensional breakdown of variance
#'     \item "prediction_deviation": Detailed deviation panels for individuals or strata
#'     \item "all": Generate all available plots (default if not specified)
#'   }
#' @param summary_obj Optional maihda_summary object from \code{summary()}.
#'   If NULL, will be computed.
#' @param n_strata Maximum number of strata to display in predicted plot.
#'   Default is 50. Use NULL for all strata.
#' @param ... Additional arguments (not currently used).
#'
#' @return A ggplot2 object, or a list of ggplot2 objects if type = "all".
#'
#' @examples
#' \donttest{
#' strata_result <- make_strata(maihda_sim_data, vars = c("gender", "race"))
#' model <- fit_maihda(health_outcome ~ age + (1 | stratum), data = strata_result$data)
#'
#' # VPC plot
#' plot(model, type = "vpc")
#'
#' # Generate all plots
#' plots <- plot(model)
#' }
#'
#' @export
#' @import ggplot2
#' @importFrom dplyr arrange
plot.maihda_model <- function(x, type = c("all", "vpc", "obs_vs_shrunken", "predicted", "risk_vs_effect", "effect_decomp", "ternary", "prediction_deviation"),
                       summary_obj = NULL, n_strata = 50, ...) {
  if (!inherits(x, "maihda_model")) {
    stop("'x' must be a maihda_model object from fit_maihda()")
  }

  object <- x


  if (missing(type)) {
    type <- "all"
  } else {
    type <- match.arg(type)
  }

  # Get summary if not provided
  if (is.null(summary_obj)) {
    summary_obj <- summary(object)
  }

  if (type == "all") {
    plots <- list()

    plots$vpc <- plot_vpc(summary_obj)

    # Try obs_vs_shrunken
    if ("stratum" %in% names(object$data)) {
      plots$obs_vs_shrunken <- tryCatch(plot_obs_vs_shrunken(object, summary_obj), error = function(e) NULL)
    }

    plots$predicted <- tryCatch(plot_predicted_strata(object, summary_obj, n_strata), error = function(e) NULL)

    top_n_labels <- if (is.null(n_strata)) 10 else min(10, n_strata)
    plots$risk_vs_effect <- tryCatch(plot_risk_vs_effect(object, summary_obj, top_n_labels), error = function(e) NULL)

    plots$effect_decomp <- tryCatch(plot_effect_decomposition(object, summary_obj, top_n_labels), error = function(e) NULL)

    ternary_out <- tryCatch(maihda_ternary_plot(object)$plot, error = function(e) NULL)
    if (!is.null(ternary_out)) plots$ternary <- ternary_out

    plots$prediction_deviation <- tryCatch(plot_prediction_deviation_panels(object$model, data = object$data, type = "auto", strata_info = object$strata_info), error = function(e) NULL)

    # print them
    for (p in plots[!sapply(plots, is.null)]) { print(p) }
    return(invisible(plots))
  } else {
    if (type == "vpc") {
      plot <- plot_vpc(summary_obj)
    } else if (type == "obs_vs_shrunken") {
      plot <- plot_obs_vs_shrunken(object, summary_obj)
    } else if (type == "predicted") {
      plot <- plot_predicted_strata(object, summary_obj, n_strata)
    } else if (type == "risk_vs_effect") {
      top_n_labels <- if (is.null(n_strata)) 10 else min(10, n_strata)
      plot <- plot_risk_vs_effect(object, summary_obj, top_n_labels)
    } else if (type == "effect_decomp") {
      top_n_labels <- if (is.null(n_strata)) 10 else min(10, n_strata)
      plot <- plot_effect_decomposition(object, summary_obj, top_n_labels)
    } else if (type == "ternary") {
      plot <- maihda_ternary_plot(object)$plot
    } else if (type == "prediction_deviation") {
      plot <- plot_prediction_deviation_panels(object$model, data = object$data, type = "auto", strata_info = object$strata_info)
    }

    return(plot)
  }
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
  p <- ggplot(vpc_data, aes(x = "", y = .data$proportion, fill = .data$component)) +
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
    geom_text(aes(label = sprintf("%.1f%%", .data$proportion * 100)),
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
    dplyr::group_by(.data$stratum) |>
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
    p <- ggplot(plot_data, aes(x = .data$observed, y = .data$shrunken)) +
      geom_point(aes(size = .data$n), alpha = 0.6, color = "#0072B2") +
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
  p <- ggplot(stratum_est, aes(x = .data$display_label, y = .data$predicted)) +
    geom_point(size = 2, color = "#0072B2") +
    geom_errorbar(aes(ymin = .data$lower, ymax = .data$upper),
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

#' Risk vs. Intersectional Effect Plot
#'
#' Creates a quadrant scatterplot comparing overall marginal predicted risk against
#' pure intersectional effects (shrunken residuals). Points represent strata.
#'
#' @param object A maihda_model object
#' @param summary_obj A maihda_summary object
#' @param top_n_labels Number of most extreme strata to label (by absolute effect size)
#' @return A ggplot2 object
#' @keywords internal
#' @import ggplot2
#' @importFrom dplyr group_by summarise n arrange desc
#' @importFrom utils head
#' @importFrom stats predict
plot_risk_vs_effect <- function(object, summary_obj, top_n_labels = 10) {
  data <- object$data

  if (!"stratum" %in% names(data)) {
    stop("'stratum' variable not found in data. Make sure to use data from make_strata().")
  }

  # 1. Retrieve the predicted values strictly across all cases in the data
  # Safe approach matching what's used in plot_prediction_deviation_panels
  model_type <- object$family$family

  if (model_type %in% c("binomial", "quasibinomial")) {
    preds <- tryCatch(predict(object$model, type = "response"), error = function(e) predict(object$model, type = "response"))
  } else if (inherits(object$model, "polr") || inherits(object$model, "clm") || inherits(object$model, "ordinal")) {
    probs <- tryCatch(predict(object$model, type = "probs"), error = function(e) predict(object$model, type = "p"))
    if (is.matrix(probs) || is.data.frame(probs)) {
      k_seq <- seq_len(ncol(probs))
      preds <- rowSums(probs * matrix(k_seq, nrow = nrow(probs), ncol = ncol(probs), byrow = TRUE))
    } else {
      preds <- rep(NA, nrow(data))
    }
  } else {
    preds <- tryCatch(predict(object$model), error = function(e) list(fit = predict(object$model)))
    if (is.list(preds) && "fit" %in% names(preds)) preds <- preds$fit
  }

  if (is.null(preds) || (is.numeric(preds) && length(preds) != nrow(data))) {
      preds <- tryCatch(predict(object$model, se.fit=FALSE), error = function(e) rep(NA, nrow(data)))
  }

  # Assign to dataframe and collapse to strata level average
  data$pred_val <- preds

  stratum_means <- data |>
    dplyr::group_by(.data$stratum) |>
    dplyr::summarize(
      mean_predicted = mean(.data$pred_val, na.rm = TRUE),
      n = dplyr::n(),
      .groups = "drop"
    )

  stratum_means$stratum <- as.character(stratum_means$stratum)

  # 2. Extract intersectional shrunken residuals
  stratum_est <- summary_obj$stratum_estimates
  if (is.null(stratum_est)) stop("No stratum estimates available for plotting")
  stratum_est$stratum <- as.character(stratum_est$stratum)

  # Merge Risk (pred) + Effect (random)
  plot_data <- merge(stratum_means, stratum_est, by = "stratum")

  # Map appropriate text labels to dots
  if (!is.null(object$strata_info) && "label" %in% names(object$strata_info)) {
    id_map <- setNames(object$strata_info$label, object$strata_info$stratum)
    plot_data$label <- id_map[plot_data$stratum]
  } else {
    plot_data$label <- paste("Stratum", plot_data$stratum)
  }

  # Compute centers
  global_mean <- mean(plot_data$mean_predicted, na.rm = TRUE)
  x_title <- "Mean Predicted Value (Overall Risk)"
  if (model_type %in% c("binomial", "quasibinomial")) x_title <- "Mean Predicted Probability (Risk)"
  if (inherits(object$model, "polr") || inherits(object$model, "clm") || inherits(object$model, "ordinal")) x_title <- "Average Expected Category Score"

  # Label the ones with largest intersectional residuals (positive or negative)
  label_data <- plot_data |>
    dplyr::arrange(dplyr::desc(abs(.data$random_effect))) |>
    utils::head(top_n_labels)

  # Create quadrant plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$mean_predicted, y = .data$random_effect)) +
    ggplot2::geom_vline(xintercept = global_mean, linetype = "dashed", color = "gray50") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::geom_point(ggplot2::aes(size = .data$n), alpha = 0.6, color = "#0072B2") +
    ggrepel::geom_label_repel(data = label_data, ggplot2::aes(label = .data$label), size = 3, min.segment.length = 0) +
    ggplot2::labs(
      title = "Risk vs. Intersectional Effect",
      subtitle = "Marginal predicted scores vs pure intersectional random effects\nTop-Right: Double Penalty (High Risk + Unique Penalty factor) \nBottom-Right: High Risk but fully explained by additive characteristics",
      x = x_title,
      y = "Random Effect (Intersectional Penalty/Advantage)",
      size = "Sample Size"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, face = "italic", size = 9),
      legend.position = "right"
    )

  return(p)
}
#' Effect Decomposition Plot
#'
#' Decomposes the total deviation from the overall mean into the additive (fixed) component
#' and the intersectional (random) component for each stratum.
#'
#' @param object A maihda_model object
#' @param summary_obj A maihda_summary object
#' @param top_n_labels Number of most extreme strata to label
#' @return A ggplot2 object
#' @keywords internal
#' @import ggplot2
#' @importFrom dplyr group_by summarise n arrange desc mutate row_number
#' @importFrom utils head
#' @importFrom stats predict setNames fitted
plot_effect_decomposition <- function(object, summary_obj, top_n_labels = 10) {
  data <- object$data

  if (!"stratum" %in% names(data)) {
    stop("'stratum' variable not found in data. Make sure to use data from make_strata().")
  }

  # Calculate both Full predictions and Fixed-effect ONLY predictions

  if (object$engine == "lme4") {
    preds_total <- tryCatch(predict(object$model, type = "response"), error = function(e) rep(NA, nrow(data)))
    preds_fixed <- tryCatch(predict(object$model, type = "response", re.form = NA), error = function(e) rep(NA, nrow(data)))
  } else if (object$engine == "brms") {
    preds_total <- tryCatch(fitted(object$model)[, "Estimate"], error = function(e) rep(NA, nrow(data)))
    preds_fixed <- tryCatch(fitted(object$model, re_formula = NA)[, "Estimate"], error = function(e) rep(NA, nrow(data)))
  } else {
    stop("Engine not supported for effect decomposition.")
  }

  data$pred_total <- preds_total
  data$pred_fixed <- preds_fixed

  global_mean <- mean(data$pred_total, na.rm = TRUE)

  # Aggregate to stratum level
  stratum_means <- data |>
    dplyr::group_by(.data$stratum) |>
    dplyr::summarize(
      mean_total = mean(.data$pred_total, na.rm = TRUE),
      mean_fixed = mean(.data$pred_fixed, na.rm = TRUE),
      .groups = "drop"
    )

  stratum_means$stratum <- as.character(stratum_means$stratum)

  # Map appropriate text labels
  if (!is.null(object$strata_info) && "label" %in% names(object$strata_info)) {
    id_map <- stats::setNames(object$strata_info$label, as.character(object$strata_info$stratum))
    stratum_means$label <- id_map[stratum_means$stratum]
  } else {
    stratum_means$label <- paste("Stratum", stratum_means$stratum)
  }

  # Calculate components
  # Total Dev = Additive Dev + Intersectional Dev
  stratum_means <- stratum_means |>
    dplyr::mutate(
      total_dev = .data$mean_total - global_mean,
      additive_dev = .data$mean_fixed - global_mean,
      intersectional_dev = .data$mean_total - .data$mean_fixed,
      abs_total_dev = abs(.data$total_dev)
    ) |>
    dplyr::arrange(.data$total_dev) |>
    dplyr::mutate(rank = dplyr::row_number())

  # Create segment definitions for stacking
  # Additive goes from 0 -> additive_dev
  # Intersectional goes from additive_dev -> total_dev
  seg_data <- rbind(
    data.frame(
      rank = stratum_means$rank,
      label = stratum_means$label,
      Component = "Additive Effect (Demographics)",
      y_start = 0,
      y_end = stratum_means$additive_dev,
      abs_total_dev = stratum_means$abs_total_dev
    ),
    data.frame(
      rank = stratum_means$rank,
      label = stratum_means$label,
      Component = "Intersectional Effect (Penalty/Advantage)",
      y_start = stratum_means$additive_dev,
      y_end = stratum_means$total_dev,
      abs_total_dev = stratum_means$abs_total_dev
    )
  )

  # Set component ordering so Additive is handled first
  seg_data$Component <- factor(seg_data$Component, levels = c("Additive Effect (Demographics)", "Intersectional Effect (Penalty/Advantage)"))

  # Label the most extreme overall cases
  label_data <- stratum_means |>
    dplyr::arrange(dplyr::desc(.data$abs_total_dev)) |>
    utils::head(top_n_labels)

  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    # Draw segments stacked directly simulating waterfall
    ggplot2::geom_segment(data = seg_data, ggplot2::aes(x = .data$rank, xend = .data$rank, y = .data$y_start, yend = .data$y_end, color = .data$Component), linewidth = 3, alpha = 0.8) +
    # Draw a point at the final Total Deviation
    ggplot2::geom_point(data = stratum_means, ggplot2::aes(x = .data$rank, y = .data$total_dev), size = 1.5, color = "black") +
    # Label extremes
    ggrepel::geom_label_repel(data = label_data, ggplot2::aes(x = .data$rank, y = .data$total_dev, label = .data$label), size = 3, min.segment.length = 0) +
    ggplot2::scale_color_manual(values = c("Additive Effect (Demographics)" = "gray60", "Intersectional Effect (Penalty/Advantage)" = "#D55E00")) +
    ggplot2::labs(
      title = "Deviation Decomposition: Additive vs. Intersectional Effects",
      subtitle = "Visualizing how much of the stratum's deviation from the global mean is due to additive risk factors vs. true intersectionality.\nThe black dot represents the Total Marginal Deviation from the mean.",
      x = "Stratum Rank (Ordered by Total Predicted Deviation)",
      y = "Deviation from Global Mean",
      color = "Effect Component"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, face = "italic", size = 9),
      legend.position = "bottom"
    )

  return(p)
}
