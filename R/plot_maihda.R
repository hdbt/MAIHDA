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
#'     \item "obs_vs_shrunken": Observed vs. shrunken stratum means. The y-axis
#'       (model-based estimate) includes the fixed effects, so for a
#'       covariate-adjusted model the distance from the diagonal reflects both
#'       shrinkage \emph{and} covariate adjustment, not shrinkage alone; it is a
#'       pure shrinkage view only for an intercept-only (null) model
#'     \item "predicted": Predicted values for each stratum with confidence intervals
#'     \item "risk_vs_effect": Quadrant scatterplot of each stratum's mean predicted outcome against its random effect
#'     \item "effect_decomp": Visualizes additive vs intersectional deviation from global mean
#'     \item "ternary": Ternary diagnostic of the relative additive, intersectional, and uncertainty signals per stratum (a normalized-magnitude diagnostic, not a variance decomposition)
#'     \item "prediction_deviation": Detailed deviation panels for individuals or strata
#'     \item "all": Generate all available plots (default if not specified)
#'   }
#' @param summary_obj Optional maihda_summary object from \code{summary()}.
#'   If NULL, will be computed.
#' @param n_strata Maximum number of strata to display in the predicted plot.
#'   When there are more strata than this, the first \code{n_strata} (in stratum
#'   order) are shown and the plot caption notes how many were omitted. Default
#'   is 50. Use NULL for all strata.
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

    plots$prediction_deviation <- tryCatch(plot_prediction_deviation_panels(object, type = "auto"), error = function(e) NULL)

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
      plot <- plot_prediction_deviation_panels(object, type = "auto")
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
  vpc_data <- summary_obj$variance_components[
    summary_obj$variance_components$component != "Total", , drop = FALSE
  ]
  component_colors <- c(
    "Between-stratum (random)" = "#E69F00",
    "Other random effects" = "#009E73",
    "Within-stratum (residual)" = "#56B4E9"
  )
  missing_colors <- setdiff(vpc_data$component, names(component_colors))
  if (length(missing_colors) > 0) {
    component_colors[missing_colors] <- "#999999"
  }

  # Create plot
  p <- ggplot(vpc_data, aes(x = "", y = .data$proportion, fill = .data$component)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_flip() +
    scale_fill_manual(values = component_colors) +
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
#' @details The x-axis is each stratum's raw observed mean; the y-axis is the
#'   model-based stratum estimate, which includes the fixed-effect contribution.
#'   For an intercept-only (null) model the vertical distance from the diagonal is
#'   pure shrinkage toward the grand mean. For a covariate-adjusted model the model
#'   estimate also moves with the stratum's covariate profile, so distance from the
#'   diagonal reflects \emph{both} shrinkage and covariate adjustment and should
#'   not be read as shrinkage alone. The caption notes which case applies.
#'
#' @param object A maihda_model object
#' @param summary_obj A maihda_summary object
#' @return A ggplot2 object
#' @keywords internal
#' @import ggplot2
#' @importFrom dplyr group_by summarise
#' @importFrom stats formula terms
plot_obs_vs_shrunken <- function(object, summary_obj) {
  data <- object$data

  observed_response <- maihda_observed_response_from_model_frame(data, object$formula)
  if (!"stratum" %in% names(data)) {
    stop("'stratum' variable not found in data. Make sure to use data from make_strata()")
  }

  observed_outcome <- maihda_observed_outcome_for_plot(observed_response, object$family)

  # Calculate observed stratum means
  obs_data <- data
  obs_data$.maihda_observed_numerator <- observed_outcome$numerator
  obs_data$.maihda_observed_denominator <- observed_outcome$denominator
  obs_data$.maihda_prior_weight <- maihda_prior_weights(object)
  obs_means <- obs_data |>
    dplyr::group_by(.data$stratum) |>
    dplyr::summarise(
      observed = maihda_observed_weighted_mean(
        .data$.maihda_observed_numerator,
        .data$.maihda_observed_denominator,
        .data$.maihda_prior_weight
      ),
      n = maihda_observed_sample_size(
        .data$.maihda_observed_numerator,
        .data$.maihda_observed_denominator
      ),
      .groups = "drop"
    )

  # Convert stratum to character for merging (to match stratum_estimates)
  obs_means$stratum <- as.character(obs_means$stratum)

  # Merge with random effects (shrunken estimates)
  stratum_est <- summary_obj$stratum_estimates
  if (!is.null(stratum_est)) {
    pred_data <- if (object$engine == "lme4") {
      maihda_stratum_predictions_lme4(object, summary_obj, scale = "response")
    } else if (object$engine == "brms") {
      maihda_stratum_predictions_brms(object, summary_obj, scale = "response")
    } else {
      stop("Unsupported engine: ", object$engine)
    }

    plot_data <- merge(obs_means, stratum_est, by = "stratum")
    pred_idx <- match(as.character(plot_data$stratum), as.character(pred_data$stratum))
    plot_data$shrunken <- pred_data$predicted_row[pred_idx]

    # The y-axis (model estimate) includes the fixed effects, so for an adjusted
    # model the vertical gap from the diagonal mixes shrinkage with covariate
    # adjustment; only an intercept-only model gives a pure shrinkage view. Flag
    # which case applies in the caption rather than letting it be misread.
    fixed_terms <- tryCatch(
      attr(stats::terms(reformulas::nobars(object$formula)), "term.labels"),
      error = function(e) character(0)
    )
    interpretation_caption <- if (length(fixed_terms) > 0) {
      paste("Adjusted model: the y-axis includes fixed effects, so distance from",
            "the diagonal reflects both shrinkage and covariate adjustment.")
    } else {
      paste("Null model: vertical distance from the diagonal is shrinkage of the",
            "stratum mean toward the grand mean.")
    }

    # Create plot
    p <- ggplot(plot_data, aes(x = .data$observed, y = .data$shrunken)) +
      geom_point(aes(size = .data$n), alpha = 0.6, color = "#0072B2") +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
      labs(
        title = "Observed vs. Shrunken Stratum Estimates",
        x = "Observed Stratum Mean",
        y = "Shrunken Estimate (with Random Effect)",
        size = "Sample Size",
        caption = interpretation_caption
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

maihda_observed_response_from_model_frame <- function(data, formula_obj) {
  response <- tryCatch(stats::model.response(data), error = function(e) NULL)
  if (!is.null(response)) {
    return(response)
  }

  outcome_var <- all.vars(formula_obj)[1]
  if (!outcome_var %in% names(data)) {
    stop("Outcome variable not found in data")
  }

  data[[outcome_var]]
}

maihda_observed_plot_values <- function(numerator, denominator = NULL) {
  numerator <- as.numeric(numerator)
  if (is.null(denominator)) {
    denominator <- rep(1, length(numerator))
  }
  data.frame(
    numerator = numerator,
    denominator = as.numeric(denominator)
  )
}

maihda_observed_complete <- function(numerator, denominator) {
  is.finite(numerator) & is.finite(denominator) & denominator > 0
}

maihda_observed_weighted_mean <- function(numerator, denominator, w = NULL) {
  keep <- maihda_observed_complete(numerator, denominator)
  if (!any(keep)) {
    return(NA_real_)
  }

  # Incorporate the model's prior/precision weights so the observed stratum mean is
  # on the same weighted footing as the weighted shrunken estimate. These are lme4
  # prior/precision weights, not a complex survey design -- no design-based
  # (e.g. Taylor-linearised) variance is computed -- so results are not
  # survey-representative. With unit weights this is the previous
  # sum(numerator)/sum(denominator).
  if (is.null(w)) {
    w <- rep(1, length(numerator))
  }
  w <- as.numeric(w)
  w[!is.finite(w)] <- 0

  sum(w[keep] * numerator[keep]) / sum(w[keep] * denominator[keep])
}

maihda_observed_sample_size <- function(numerator, denominator) {
  keep <- maihda_observed_complete(numerator, denominator)
  if (!any(keep)) {
    return(0)
  }

  sum(denominator[keep])
}

maihda_observed_outcome_for_plot <- function(x, family = NULL) {
  fam_name <- if (!is.null(family) && !is.null(family$family)) family$family else NULL
  is_binomial <- !is.null(fam_name) && fam_name %in% c("binomial", "quasibinomial")

  if ((is.matrix(x) || is.data.frame(x)) && is_binomial && ncol(x) == 2) {
    x_mat <- as.matrix(x)
    if (!all(vapply(seq_len(ncol(x_mat)), function(j) is.numeric(x_mat[, j]), logical(1)))) {
      stop("Observed-vs-shrunken plots require numeric success/failure counts for matrix binomial outcomes.",
           call. = FALSE)
    }
    totals <- rowSums(x_mat, na.rm = FALSE)
    numerator <- x_mat[, 1]
    numerator[!is.finite(totals) | totals <= 0] <- NA_real_
    return(maihda_observed_plot_values(numerator, totals))
  }

  if (is.numeric(x)) {
    return(maihda_observed_plot_values(x))
  }
  if (is.logical(x)) {
    return(maihda_observed_plot_values(x))
  }
  if (is.factor(x)) {
    if (is_binomial && nlevels(x) == 2) {
      return(maihda_observed_plot_values(x == levels(x)[2]))
    }
    stop("Observed-vs-shrunken plots require a numeric outcome, or a two-level factor for binomial models.",
         call. = FALSE)
  }
  if (is.character(x) && is_binomial && length(unique(stats::na.omit(x))) == 2) {
    levels_x <- sort(unique(stats::na.omit(x)))
    return(maihda_observed_plot_values(x == levels_x[2]))
  }

  stop("Observed-vs-shrunken plots require a numeric outcome, or a binary outcome that can be converted to 0/1.",
       call. = FALSE)
}

#' Plot Predicted Stratum Values with Confidence Intervals
#'
#' @param object A maihda_model object
#' @param summary_obj A maihda_summary object
#' @param n_strata Maximum number of strata to display (the first n_strata, in stratum order)
#' @param scale Prediction scale: "response" (default) or "link"
#' @return A ggplot2 object
#' @keywords internal
#' @import ggplot2
#' @importFrom dplyr arrange slice
plot_predicted_strata <- function(object, summary_obj, n_strata, scale = c("response", "link")) {
  scale <- match.arg(scale)

  pred_data <- if (object$engine == "lme4") {
    maihda_stratum_predictions_lme4(object, summary_obj, scale = scale)
  } else if (object$engine == "brms") {
    maihda_stratum_predictions_brms(object, summary_obj, scale = scale)
  } else {
    stop("Unsupported engine: ", object$engine)
  }

  # Weight the across-strata reference by each stratum's summed prior weights
  # (w_sum), which equals the row count for an unweighted model.
  ref_weights <- if ("w_sum" %in% names(pred_data)) pred_data$w_sum else pred_data$n
  fixed_reference <- stats::weighted.mean(pred_data$fixed_row, ref_weights, na.rm = TRUE)

  # Get stratum estimates
  stratum_est <- summary_obj$stratum_estimates

  if (is.null(stratum_est) || nrow(stratum_est) == 0) {
    stop("No stratum estimates available for plotting")
  }

  pred_idx <- match(as.character(stratum_est$stratum), as.character(pred_data$stratum))
  stratum_est$predicted <- pred_data$predicted_row[pred_idx]
  stratum_est$lower <- pred_data$lower_row[pred_idx]
  stratum_est$upper <- pred_data$upper_row[pred_idx]

  # Keep original order (no sorting). n_strata is a MAXIMUM: when there are more
  # strata than that, show the first n_strata in stratum order rather than
  # thinning an evenly-spaced subset across all of them. The old stride sampling
  # silently dropped strata from the middle while implying full coverage; the
  # caption below records how many were omitted so the cap is not silent.
  n_total_strata <- nrow(stratum_est)
  truncated_strata <- !is.null(n_strata) && n_total_strata > n_strata
  if (truncated_strata) {
    stratum_est <- utils::head(stratum_est, n_strata)
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
    geom_hline(yintercept = fixed_reference, linetype = "dashed", color = "red", alpha = 0.7) +
    labs(
      title = "Predicted Subgroup Values with Conditional 95% Intervals",
      x = "Stratum",
      y = "Predicted Value",
      caption = paste0(
        "Intervals reflect random-effect (conditional) uncertainty only, ",
        "not fixed-effect uncertainty.\nDashed line is the mean fixed-only prediction.",
        if (truncated_strata) {
          sprintf("\nShowing the first %d of %d strata (n_strata = %d).",
                  n_strata, n_total_strata, n_strata)
        } else {
          ""
        }
      )
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

#' Mean Prediction vs. Stratum Random Effect Plot
#'
#' Creates a quadrant scatterplot comparing each stratum's mean predicted outcome
#' against its stratum random effect (shrunken between-stratum deviation). Points
#' represent strata. Whether a higher predicted value is "worse" or "better"
#' depends on the outcome, so the axes are not framed as risk. The random effect
#' equals the \emph{pure} intersectional (interaction) component only when the
#' additive main effects of the strata variables are included in the model;
#' otherwise it also absorbs those omitted main effects.
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

  if (object$engine == "brms" || inherits(object$model, "brmsfit")) {
    if (!requireNamespace("brms", quietly = TRUE)) {
      stop("Package 'brms' is required to plot the mean prediction vs. stratum ",
           "random effect for brms models.", call. = FALSE)
    }
    preds <- stats::fitted(object$model, newdata = data, re_formula = NA, summary = TRUE)[, "Estimate"]
  } else if (model_type %in% c("binomial", "quasibinomial")) {
    preds <- tryCatch(
      predict(object$model, newdata = data, type = "response", re.form = NA),
      error = function(e) predict(object$model, type = "response", re.form = NA)
    )
  } else if (inherits(object$model, "polr") || inherits(object$model, "clm") || inherits(object$model, "ordinal")) {
    probs <- tryCatch(
      predict(object$model, newdata = data, type = "probs"),
      error = function(e) predict(object$model, newdata = data, type = "p")
    )
    if (is.matrix(probs) || is.data.frame(probs)) {
      k_seq <- seq_len(ncol(probs))
      preds <- rowSums(probs * matrix(k_seq, nrow = nrow(probs), ncol = ncol(probs), byrow = TRUE))
    } else {
      preds <- rep(NA, nrow(data))
    }
  } else {
    preds <- tryCatch(
      predict(object$model, newdata = data, type = "response", re.form = NA),
      error = function(e) list(fit = predict(object$model, type = "response", re.form = NA))
    )
    if (is.list(preds) && "fit" %in% names(preds)) preds <- preds$fit
  }

  if (is.null(preds) || (is.numeric(preds) && length(preds) != nrow(data))) {
      preds <- tryCatch(
        predict(object$model, newdata = data, type = "response", re.form = NA, se.fit = FALSE),
        error = function(e) rep(NA, nrow(data))
      )
  }
  if (is.matrix(preds) || is.data.frame(preds)) {
    preds <- preds[, 1]
  }
  preds <- as.numeric(preds)
  if (length(preds) != nrow(data)) {
    stop("Could not compute one prediction per analytic row.", call. = FALSE)
  }

  # Assign to dataframe and collapse to strata level average. The model's
  # prior/precision weights make the per-stratum mean (and the reference centre
  # below) reflect the weighted fit; these are lme4 prior weights, not a complex
  # survey design, so the result is not survey-representative. For an unweighted
  # model the weights are all 1 and this reduces to the previous plain means.
  data$pred_val <- preds
  data$.maihda_w <- maihda_prior_weights(object)

  stratum_means <- data |>
    dplyr::group_by(.data$stratum) |>
    dplyr::summarize(
      mean_predicted = stats::weighted.mean(.data$pred_val, .data$.maihda_w, na.rm = TRUE),
      n = dplyr::n(),
      w_sum = sum(.data$.maihda_w, na.rm = TRUE),
      .groups = "drop"
    )

  stratum_means$stratum <- as.character(stratum_means$stratum)

  # 2. Extract intersectional shrunken residuals
  stratum_est <- summary_obj$stratum_estimates
  if (is.null(stratum_est)) stop("No stratum estimates available for plotting")
  stratum_est$stratum <- as.character(stratum_est$stratum)

  # Merge mean prediction + stratum random effect
  plot_data <- merge(stratum_means, stratum_est, by = "stratum")

  # Map appropriate text labels to dots
  if (!is.null(object$strata_info) && "label" %in% names(object$strata_info)) {
    id_map <- setNames(object$strata_info$label, object$strata_info$stratum)
    plot_data$label <- id_map[plot_data$stratum]
  } else {
    plot_data$label <- paste("Stratum", plot_data$stratum)
  }

  # Compute the reference centre as the population mean, weighting each stratum by
  # its summed prior weights (w_sum, = stratum size for an unweighted model), so
  # common and rare strata are represented in proportion to their weight -- matching
  # the weighted reference line in plot_predicted_strata().
  ref_w <- if ("w_sum" %in% names(plot_data)) plot_data$w_sum else plot_data$n
  global_mean <- if (any(is.finite(ref_w))) {
    stats::weighted.mean(plot_data$mean_predicted, ref_w, na.rm = TRUE)
  } else {
    mean(plot_data$mean_predicted, na.rm = TRUE)
  }
  x_title <- "Mean Predicted Value"
  if (model_type %in% c("binomial", "quasibinomial")) x_title <- "Mean Predicted Probability"
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
      title = "Mean Prediction vs. Stratum Random Effect",
      subtitle = paste0(
        "Mean predicted outcome per stratum vs the stratum random effect. ",
        "Whether a higher predicted value is 'worse' or 'better' depends on the ",
        "outcome.\nThe random effect is the pure intersectional (interaction) ",
        "effect only when the strata main effects are in the model; otherwise it ",
        "also includes those additive main effects."
      ),
      x = x_title,
      y = "Stratum random effect (between-stratum deviation)",
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

  # Compute full and fixed-only predictions on the LINK scale. The additive
  # decomposition (total = additive + intersectional) is only exact on the model
  # scale: eta = X*beta + u_stratum. On the response scale, for non-identity links
  # (logit/log) the split is not additive. For Gaussian/identity the link scale
  # equals the response scale, so this is unchanged there.
  if (object$engine == "lme4") {
    preds_total <- tryCatch(predict(object$model, type = "link"), error = function(e) rep(NA, nrow(data)))
    preds_fixed <- tryCatch(predict(object$model, type = "link", re.form = NA), error = function(e) rep(NA, nrow(data)))
  } else if (object$engine == "brms") {
    preds_total <- tryCatch(brms::posterior_linpred(object$model, summary = TRUE)[, "Estimate"], error = function(e) rep(NA, nrow(data)))
    preds_fixed <- tryCatch(brms::posterior_linpred(object$model, re_formula = NA, summary = TRUE)[, "Estimate"], error = function(e) rep(NA, nrow(data)))
  } else {
    stop("Engine not supported for effect decomposition.")
  }

  data$pred_total <- preds_total
  data$pred_fixed <- preds_fixed
  # The model's prior/precision weights so the per-stratum and global means reflect
  # the weighted fit (the stratum random-effect component below is the
  # weight-invariant BLUP). These are lme4 prior weights, not a complex survey
  # design, so the result is not survey-representative. Unit weights reproduce the
  # previous unweighted means exactly.
  data$.maihda_w <- maihda_prior_weights(object)

  global_mean <- stats::weighted.mean(data$pred_total, data$.maihda_w, na.rm = TRUE)

  # Aggregate to stratum level
  stratum_means <- data |>
    dplyr::group_by(.data$stratum) |>
    dplyr::summarize(
      mean_total = stats::weighted.mean(.data$pred_total, .data$.maihda_w, na.rm = TRUE),
      mean_fixed = stats::weighted.mean(.data$pred_fixed, .data$.maihda_w, na.rm = TRUE),
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

  # Calculate components: total_dev = additive_dev + intersectional_dev.
  # The intersectional (stratum) component is the stratum random effect (BLUP)
  # itself, taken from the summary, NOT total-minus-fixed. With additional random
  # effects (e.g. (1 | site)) total-minus-fixed would also absorb those, wrongly
  # attributing them to the stratum; using the stratum random effect isolates the
  # intersectional component. For the canonical single-stratum model the two are
  # identical. Strata absent from the random-effect table contribute 0.
  re_map <- stats::setNames(
    as.numeric(summary_obj$stratum_estimates$random_effect),
    as.character(summary_obj$stratum_estimates$stratum)
  )
  stratum_means$intersectional_dev <- unname(re_map[stratum_means$stratum])
  stratum_means$intersectional_dev[is.na(stratum_means$intersectional_dev)] <- 0

  stratum_means <- stratum_means |>
    dplyr::mutate(
      additive_dev = .data$mean_fixed - global_mean,
      total_dev = .data$additive_dev + .data$intersectional_dev,
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
      Component = "Fixed-effect component",
      y_start = 0,
      y_end = stratum_means$additive_dev,
      abs_total_dev = stratum_means$abs_total_dev
    ),
    data.frame(
      rank = stratum_means$rank,
      label = stratum_means$label,
      Component = "Stratum random-effect component",
      y_start = stratum_means$additive_dev,
      y_end = stratum_means$total_dev,
      abs_total_dev = stratum_means$abs_total_dev
    )
  )

  # Set component ordering so Additive is handled first
  seg_data$Component <- factor(seg_data$Component, levels = c("Fixed-effect component", "Stratum random-effect component"))

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
    ggplot2::scale_color_manual(values = c("Fixed-effect component" = "gray60", "Stratum random-effect component" = "#D55E00")) +
    ggplot2::labs(
      title = "Deviation Decomposition: Fixed vs. Stratum-Random Components",
      subtitle = paste0(
        "Stratum deviation split into the fixed-effect component and the stratum ",
        "random effect (BLUP), on the model (link) scale.\nThe black dot is their ",
        "sum; any other random effects (e.g. (1 | site)) are not included. The ",
        "stratum component is the pure intersectional (interaction) effect only when ",
        "the strata main effects are in the model; otherwise it also absorbs them."
      ),
      x = "Stratum Rank (Ordered by Total Predicted Deviation)",
      y = "Deviation from Global Mean (link scale)",
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
