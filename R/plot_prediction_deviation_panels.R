maihda_binomial_observed_01 <- function(x, n) {
  if (is.null(x) || length(x) != n || !maihda_is_binary_vector(x)) {
    return(rep(NA_integer_, n))
  }

  maihda_binary_to_01(x)
}

maihda_binomial_abs_deviance_residual <- function(obs_outcome_01, fitted) {
  out <- rep(0, length(fitted))
  known_obs <- !is.na(obs_outcome_01) &
    obs_outcome_01 %in% c(0L, 1L) &
    is.finite(fitted)

  if (!any(known_obs)) {
    return(out)
  }

  p <- pmin(pmax(fitted[known_obs], .Machine$double.eps), 1 - .Machine$double.eps)
  y <- obs_outcome_01[known_obs]
  dev <- ifelse(y == 1L, -2 * log(p), -2 * log1p(-p))
  out[known_obs] <- sqrt(pmax(dev, 0))
  out
}

maihda_prediction_panel_auto_type <- function(model) {
  if (inherits(model, "polr") || inherits(model, "clm") || inherits(model, "ordinal")) {
    return("ordinal")
  }

  fam <- maihda_family(model)
  fam_name <- if (!is.null(fam) && !is.null(fam$family)) fam$family else NULL
  if (!is.null(fam_name) && fam_name %in% c("binomial", "quasibinomial", "bernoulli")) {
    return("binomial")
  }
  if (!is.null(fam_name) && fam_name %in% c("cumulative", "sratio", "cratio", "acat", "ordinal")) {
    return("ordinal")
  }

  "gaussian"
}

maihda_prediction_panel_fitted <- function(model, data, type) {
  if (inherits(model, "brmsfit")) {
    if (!requireNamespace("brms", quietly = TRUE)) {
      stop("Package 'brms' is required to plot prediction deviations from brms models.",
           call. = FALSE)
    }
    fit <- stats::fitted(model, newdata = data, summary = TRUE)
    if (is.null(dim(fit)) || !"Estimate" %in% colnames(fit)) {
      stop("Could not extract fitted estimates from brms model.", call. = FALSE)
    }
    # Derive the interval half-width from the posterior 2.5/97.5% quantiles so the
    # downstream `estimate +/- 1.96 * se` reflects the actual posterior spread
    # rather than assuming Est.Error (the posterior SD) describes a normal
    # interval. (This still renders a symmetric bar; full asymmetric posterior
    # intervals would require carrying the quantiles through the aggregation.)
    se <- if (all(c("Q2.5", "Q97.5") %in% colnames(fit))) {
      (fit[, "Q97.5"] - fit[, "Q2.5"]) / (2 * stats::qnorm(0.975))
    } else if ("Est.Error" %in% colnames(fit)) {
      fit[, "Est.Error"]
    } else {
      # NA (not 0) so downstream CI bars are omitted rather than collapsed.
      rep(NA_real_, nrow(data))
    }
    return(list(fit = as.numeric(fit[, "Estimate"]), se.fit = as.numeric(se)))
  }

  # SE fallbacks below are NA_real_ — not 0 — because predict() for lme4::merMod
  # (and some other mixed-model classes) does not implement se.fit. Returning 0
  # produced ci_lower == ci_upper == fitted, i.e. fake zero-width "95% CI" bars.
  # NA propagates through fitted +/- 1.96 * se and ggplot drops the geom_errorbar
  # layer for those rows, which honestly communicates "no SE available".
  if (type == "binomial") {
    preds <- tryCatch(
      predict(model, newdata = data, type = "response", se.fit = TRUE),
      error = function(e) list(
        fit = predict(model, newdata = data, type = "response"),
        se.fit = rep(NA_real_, nrow(data))
      )
    )
  } else {
    preds <- tryCatch(
      predict(model, newdata = data, se.fit = TRUE),
      error = function(e) list(
        fit = predict(model, newdata = data),
        se.fit = rep(NA_real_, nrow(data))
      )
    )
  }

  if (is.numeric(preds)) {
    preds <- list(fit = preds, se.fit = rep(NA_real_, nrow(data)))
  }
  preds$fit <- as.numeric(preds$fit)
  if (length(preds$fit) != nrow(data)) {
    stop("Predictions must have one fitted value per row in 'data'. ",
         "Use the original model frame or provide prediction data compatible with the fitted model.",
         call. = FALSE)
  }
  if (is.null(preds$se.fit) || length(preds$se.fit) != nrow(data)) {
    preds$se.fit <- rep(NA_real_, nrow(data))
  } else {
    preds$se.fit <- as.numeric(preds$se.fit)
  }
  preds
}

maihda_prediction_panel_ordinal_probs <- function(model, data) {
  probs <- tryCatch(
    predict(model, newdata = data, type = "probs"),
    error = function(e) NULL
  )
  if (is.null(probs)) {
    probs <- tryCatch(
      predict(model, newdata = data, type = "p"),
      error = function(e) NULL
    )
  }
  if (is.null(probs)) {
    probs <- tryCatch(
      predict(model, newdata = data, type = "prob"),
      error = function(e) NULL
    )
  }

  if (is.list(probs) && !is.matrix(probs) && !is.data.frame(probs)) {
    if (!is.null(probs$fit)) {
      probs <- probs$fit
    } else if (!is.null(probs$prob)) {
      probs <- probs$prob
    }
  }

  if (is.null(probs) || (!is.matrix(probs) && !is.data.frame(probs))) {
    stop("Could not extract probability matrix from ordinal model.", call. = FALSE)
  }

  probs <- as.data.frame(probs)
  if (nrow(probs) != nrow(data)) {
    stop("Ordinal predictions must have one probability row per row in 'data'. ",
         "Use the original model frame or provide prediction data compatible with the fitted model.",
         call. = FALSE)
  }

  probs[] <- lapply(probs, as.numeric)
  probs
}

maihda_prediction_panel_binomial_residuals <- function(model, data, fitted, obs_outcome_01) {
  aligned_obs <- length(obs_outcome_01) == length(fitted)
  aligned_resids <- if (aligned_obs) {
    maihda_binomial_abs_deviance_residual(obs_outcome_01, fitted)
  } else {
    rep(0, length(fitted))
  }
  has_aligned_obs <- aligned_obs && any(
    !is.na(obs_outcome_01) &
      obs_outcome_01 %in% c(0L, 1L) &
      is.finite(fitted)
  )

  if (has_aligned_obs) {
    return(aligned_resids)
  }

  if (inherits(model, "brmsfit")) {
    return(aligned_resids)
  }

  model_resids <- tryCatch(abs(residuals(model, type = "deviance")), error = function(e) NULL)
  if (is.numeric(model_resids) && length(model_resids) == nrow(data)) {
    return(model_resids)
  }

  aligned_resids
}

#' Plot Prediction Deviation Panels
#'
#' @description Creates an advanced, publication-ready two-panel dashboard for visualizing
#' predicted values and identifying deviant cases in linear, binomial, or ordinal models.
#'
#' @param model A fitted model object (e.g., from `lm()`, `glm()`, `MASS::polr()`, or `lme4::glmer()`).
#' @param data The original data frame used to fit the model. If `NULL`, attempts to extract from the model.
#' @param type Model type: "auto" (default), "gaussian", "binomial", or "ordinal".
#' @param ordinal_mode For ordinal models: "surprise" (default, based on observation probability) or "expected_score".
#' @param top_n_labels Number of extreme/deviant cases to label on the plot. Default is 5.
#' @param strata_info Optional data frame of strata labels, generally extracted from `maihda_model` objects.
#'
#' @return A `patchwork` object containing two `ggplot2` panels.
#' @importFrom rlang check_installed .data
#' @importFrom stats predict formula residuals model.frame
#' @importFrom utils head
#' @import ggplot2
#' @import patchwork
#' @import dplyr
#' @import tidyr
#' @import ggrepel
#' @export
#'
plot_prediction_deviation_panels <- function(model, data = NULL,
                                             type = c("auto", "gaussian", "binomial", "ordinal"),
                                             ordinal_mode = c("surprise", "expected_score"),
                                             top_n_labels = 5,
                                             strata_info = NULL) {

  rlang::check_installed(c("ggplot2", "patchwork", "dplyr", "tidyr", "ggrepel"))

  type <- match.arg(type)
  ordinal_mode <- match.arg(ordinal_mode)

  # Check if model is a maihda_model
  if (inherits(model, "maihda_model")) {
    if (is.null(data)) data <- model$data
    strata_info <- model$strata_info
    model <- model$model
  }

  if (is.null(data)) {
    data <- tryCatch(
      {
        if (inherits(model, "merMod")) {
          model@frame
        } else {
          model.frame(model)
        }
      },
      error = function(e) stop("Please provide the original 'data' argument, could not extract from model.")
    )
  }

  # Auto-detect model type if requested
  if (type == "auto") {
    type <- maihda_prediction_panel_auto_type(model)
  }

  get_extreme_labels <- function(df, metric_col, n) {
    df |> dplyr::arrange(dplyr::desc(abs(.data[[metric_col]]))) |> utils::head(n)
  }

  if (type == "gaussian") {
    # GAUSSIAN / LINEAR LOGIC
    # Approximate predict, some packages handle se.fit differently, so wrap safely
    preds <- maihda_prediction_panel_fitted(model, data, "gaussian")

    df <- data |>
      dplyr::mutate(
        id = dplyr::row_number(),
        fitted = preds$fit,
        se = preds$se.fit
      )

    if ("stratum" %in% names(df)) {
      df <- df |>
        dplyr::group_by(.data$stratum) |>
        dplyr::summarize(
          fitted = mean(.data$fitted, na.rm = TRUE),
          se = mean(.data$se, na.rm = TRUE),
          .groups = "drop"
        )

      if (!is.null(strata_info) && "label" %in% names(strata_info)) {
        id_map <- setNames(strata_info$label, strata_info$stratum)
        df$id <- id_map[as.character(df$stratum)]
      } else {
        df$id <- paste("Stratum", df$stratum)
      }
      x_label <- "Stratum Rank"
    } else {
      x_label <- "Case Rank"
    }

    df <- df |>
      dplyr::mutate(
        ci_lower = .data$fitted - 1.96 * .data$se,
        ci_upper = .data$fitted + 1.96 * .data$se,
        mean_fitted = mean(.data$fitted, na.rm = TRUE),
        deviation = .data$fitted - .data$mean_fitted,
        abs_deviation = abs(.data$deviation),
        direction = ifelse(.data$deviation > 0, "Above Mean", "Below Mean")
      ) |>
      dplyr::arrange(.data$fitted) |>
      dplyr::mutate(rank = dplyr::row_number())

    label_df <- get_extreme_labels(df, "deviation", top_n_labels)

    p1 <- ggplot2::ggplot(df, ggplot2::aes(x = .data$fitted)) +
      ggplot2::geom_density(fill = "gray80", alpha = 0.5) +
      ggplot2::geom_vline(ggplot2::aes(xintercept = .data$mean_fitted[1]), linetype = "dashed", color = "black") +
      ggplot2::geom_rug(data = label_df, color = "red", linewidth = 1) +
      ggplot2::labs(title = "Distribution of Fitted Values", x = NULL, y = "Density") +
      ggplot2::theme_minimal()

    p2 <- ggplot2::ggplot(df, ggplot2::aes(x = .data$rank, y = .data$fitted)) +
      ggplot2::geom_segment(ggplot2::aes(xend = .data$rank, yend = .data$mean_fitted), color = "gray60") +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$ci_lower, ymax = .data$ci_upper), width = 0, color = "gray50", alpha = 0.5) +
      ggplot2::geom_point(ggplot2::aes(color = .data$direction, size = .data$abs_deviation)) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = .data$mean_fitted[1]), linetype = "dashed") +
      ggrepel::geom_label_repel(data = label_df, ggplot2::aes(label = .data$id), size = 3, min.segment.length = 0) +
      ggplot2::scale_color_manual(values = c("Above Mean" = "#0072B2", "Below Mean" = "#D55E00")) +
      ggplot2::labs(
        x = x_label, y = "Fitted Value", color = "Direction", size = "Deviation\nMagnitude",
        caption = if (identical(x_label, "Stratum Rank")) {
          paste("Stratum intervals are approximate: the mean of the individual",
                "prediction SEs, not the SE of the stratum-mean prediction.")
        } else {
          NULL
        }
      ) +
      ggplot2::theme_minimal()

    return(patchwork::wrap_plots(p1, p2, ncol = 1, heights = c(1, 2)))

  } else if (type == "binomial") {
    # BINOMIAL / LOGISTIC LOGIC
    preds <- maihda_prediction_panel_fitted(model, data, "binomial")

    # Try to extract response variable
    form <- tryCatch(formula(model), error = function(e) NULL)
    obs_outcome <- NULL
    obs_outcome_01 <- rep(NA_integer_, nrow(data))
    if (!is.null(form)) {
      resp_name <- as.character(form)[2]
      if (resp_name %in% names(data)) {
        raw_outcome <- data[[resp_name]]
        obs_outcome <- as.factor(raw_outcome)
        obs_outcome_01 <- maihda_binomial_observed_01(raw_outcome, nrow(data))
      }
    }
    if (is.null(obs_outcome)) {
      obs_outcome <- factor(rep(NA, nrow(data)))
    }

    resids <- maihda_prediction_panel_binomial_residuals(
      model, data, preds$fit, obs_outcome_01
    )

    df <- data |>
      dplyr::mutate(
        id = dplyr::row_number(),
        obs_outcome = obs_outcome,
        obs_outcome_01 = obs_outcome_01,
        fitted = preds$fit,
        se = preds$se.fit,
        abs_res_dev = resids
      )

is_aggregated <- "stratum" %in% names(df)

    if (is_aggregated) {
      df <- df |>
        dplyr::group_by(.data$stratum) |>
        dplyr::summarize(
          fitted = mean(.data$fitted, na.rm = TRUE),
          se = mean(.data$se, na.rm = TRUE),
          abs_res_dev = mean(.data$abs_res_dev, na.rm = TRUE),
          .groups = "drop"
        )

      if (!is.null(strata_info) && "label" %in% names(strata_info)) {
        id_map <- setNames(strata_info$label, strata_info$stratum)
        df$id <- id_map[as.character(df$stratum)]
      } else {
        df$id <- paste("Stratum", df$stratum)
      }
      x_label <- "Stratum Rank"
    } else {
      x_label <- "Case Rank"
    }

    df <- df |>
      dplyr::mutate(
        ci_lower = pmax(0, .data$fitted - 1.96 * .data$se),
        ci_upper = pmin(1, .data$fitted + 1.96 * .data$se),
        mean_fitted = mean(.data$fitted, na.rm = TRUE),
        deviation = .data$fitted - .data$mean_fitted,
        direction = ifelse(.data$deviation > 0, "Above Mean", "Below Mean")
      )

    if (!is_aggregated) {
      wrong <- rep(NA_character_, nrow(df))
      known_obs <- !is.na(df$obs_outcome_01)
      wrong[known_obs] <- ifelse(
        (df$fitted[known_obs] > 0.5 & df$obs_outcome_01[known_obs] == 0) |
          (df$fitted[known_obs] < 0.5 & df$obs_outcome_01[known_obs] == 1),
        "Wrong",
        "Correct"
      )
      df$wrong <- factor(wrong, levels = c("Correct", "Wrong"))
    }

    df <- df |>
      dplyr::arrange(.data$fitted) |>
      dplyr::mutate(rank = dplyr::row_number())

    label_df <- df |> dplyr::arrange(dplyr::desc(.data$abs_res_dev)) |> utils::head(top_n_labels)

    p1 <- ggplot2::ggplot(df, ggplot2::aes(x = .data$fitted)) +
      ggplot2::geom_density(fill = "gray80", alpha = 0.5) +
      ggplot2::geom_vline(ggplot2::aes(xintercept = .data$mean_fitted[1]), linetype = "dashed", color = "black") +
      ggplot2::geom_rug(data = label_df, color = "red", linewidth = 1) +
      ggplot2::labs(title = "Distribution of Predicted Probabilities", x = NULL, y = "Density") +
      ggplot2::theme_minimal()

    p2 <- ggplot2::ggplot(df, ggplot2::aes(x = .data$rank, y = .data$fitted)) +
      ggplot2::geom_segment(ggplot2::aes(xend = .data$rank, yend = .data$mean_fitted), color = "gray60", alpha = 0.5) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$ci_lower, ymax = .data$ci_upper), width = 0, color = "gray70", alpha = 0.3)

    if (is_aggregated) {
      p2 <- p2 +
        ggplot2::geom_point(ggplot2::aes(color = .data$direction, size = .data$abs_res_dev), alpha = 0.8) +
        ggplot2::labs(
          x = x_label, y = "Predicted Probability", color = "Direction", size = "|Deviance\nResidual|",
          caption = paste("Stratum intervals are approximate: the mean of the individual",
                          "prediction SEs, not the SE of the stratum-mean prediction.")
        )
    } else {
      p2 <- p2 +
        ggplot2::geom_point(ggplot2::aes(color = .data$direction, size = .data$abs_res_dev, shape = .data$obs_outcome), alpha = 0.8)

      if (any(df$wrong == "Wrong", na.rm = TRUE)) {
        p2 <- p2 + ggplot2::geom_point(data = dplyr::filter(df, .data$wrong == "Wrong"), shape = 1, color = "red", ggplot2::aes(size = .data$abs_res_dev + 0.5))
      }
      p2 <- p2 + ggplot2::labs(x = x_label, y = "Predicted Probability", color = "Direction", size = "|Deviance\nResidual|", shape = "Observed")
    }

    p2 <- p2 +
      ggplot2::geom_hline(ggplot2::aes(yintercept = .data$mean_fitted[1]), linetype = "dashed") +
        ggrepel::geom_label_repel(data = label_df, ggplot2::aes(label = .data$id), size = 3, min.segment.length = 0) +
      ggplot2::scale_color_manual(values = c("Above Mean" = "#0072B2", "Below Mean" = "#D55E00")) +
      ggplot2::theme_minimal()

    return(patchwork::wrap_plots(p1, p2, ncol = 1, heights = c(1, 2)))

  } else if (type == "ordinal") {
    # ORDINAL LOGIC
    probs <- maihda_prediction_panel_ordinal_probs(model, data)
    prob_mat <- as.matrix(probs)
    prob_cols <- colnames(probs)

    form <- tryCatch(formula(model), error = function(e) NULL)
    obs_cat <- rep(NA, nrow(data))
    if (!is.null(form)) {
      resp_name <- as.character(form)[2]
      if (resp_name %in% names(data)) obs_cat <- as.character(data[[resp_name]])
    }

    if (ordinal_mode == "surprise") {
      df <- as.data.frame(probs) |>
        dplyr::mutate(
          id = dplyr::row_number(),
          obs_cat = obs_cat
        )

      k_seq <- seq_len(ncol(prob_mat))
      df$expected_score <- rowSums(prob_mat * matrix(k_seq, nrow = nrow(prob_mat), ncol = ncol(prob_mat), byrow = TRUE))

      # Probability of observed category
      df$observed_prob <- NA
      for (i in seq_len(nrow(df))) {
        col_idx <- match(df$obs_cat[i], prob_cols)
        if (!is.na(col_idx)) {
          df$observed_prob[i] <- prob_mat[i, col_idx]
        }
      }

      if ("stratum" %in% names(data)) {
        df$stratum <- data$stratum
        df <- df |>
          dplyr::group_by(.data$stratum) |>
          dplyr::summarize(
            dplyr::across(tidyselect::all_of(prob_cols), \(x) mean(x, na.rm = TRUE)),
            expected_score = mean(.data$expected_score, na.rm = TRUE),
            observed_prob = mean(.data$observed_prob, na.rm = TRUE),
            .groups = "drop"
          )

        if (!is.null(strata_info) && "label" %in% names(strata_info)) {
          id_map <- setNames(strata_info$label, strata_info$stratum)
          df$id <- id_map[as.character(df$stratum)]
        } else {
          df$id <- paste("Stratum", df$stratum)
        }
        x_label <- "Stratum Rank (Ordered by Expected Category Score)"
      } else {
        x_label <- "Case Rank (Ordered by Expected Category Score)"
      }

      df$surprise <- -log(df$observed_prob)

      df <- df |>
        dplyr::arrange(.data$expected_score) |>
        dplyr::mutate(rank = dplyr::row_number())

      df_long <- df |>
        tidyr::pivot_longer(cols = tidyselect::all_of(prob_cols), names_to = "Category", values_to = "Probability") |>
        dplyr::mutate(Category = factor(.data$Category, levels = prob_cols))

      label_df <- df |> dplyr::arrange(dplyr::desc(.data$surprise)) |> utils::head(top_n_labels)

      p1 <- ggplot2::ggplot(df_long, ggplot2::aes(x = .data$rank, y = .data$Probability, fill = .data$Category)) +
        ggplot2::geom_area(alpha = 0.8) +
        ggplot2::scale_fill_viridis_d(option = "magma") +
        ggplot2::labs(title = "Predicted Category Probability Structure", x = NULL, y = "Probability") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_blank())

      p2 <- ggplot2::ggplot(df, ggplot2::aes(x = .data$rank, y = .data$surprise)) +
        ggplot2::geom_segment(ggplot2::aes(xend = .data$rank, yend = 0), color = "gray50") +
        ggplot2::geom_point(ggplot2::aes(color = .data$surprise, size = .data$surprise)) +
        ggplot2::scale_color_viridis_c(option = "inferno") +
        ggrepel::geom_label_repel(data = label_df, ggplot2::aes(label = .data$id), size = 3) +
        ggplot2::labs(x = x_label, y = "Surprise\n(-log(P(Observed)))", color = "Surprise", size = "Surprise") +
        ggplot2::theme_minimal()

      return(patchwork::wrap_plots(p1, p2, ncol = 1, heights = c(1, 2)))

    } else {
      # expected_score
      k_seq <- seq_len(ncol(prob_mat))
      exp_scores <- rowSums(prob_mat * matrix(k_seq, nrow = nrow(prob_mat), ncol = ncol(prob_mat), byrow = TRUE))

      df <- data |>
        dplyr::mutate(
          id = dplyr::row_number(),
          fitted = exp_scores
        )

      if ("stratum" %in% names(df)) {
        df <- df |>
          dplyr::group_by(.data$stratum) |>
          dplyr::summarize(
            fitted = mean(.data$fitted, na.rm = TRUE),
            .groups = "drop"
          )

        if (!is.null(strata_info) && "label" %in% names(strata_info)) {
          id_map <- setNames(strata_info$label, strata_info$stratum)
          df$id <- id_map[as.character(df$stratum)]
        } else {
          df$id <- paste("Stratum", df$stratum)
        }
        x_label <- "Stratum Rank"
      } else {
        x_label <- "Case Rank"
      }

      df <- df |>
        dplyr::mutate(
          mean_fitted = mean(.data$fitted, na.rm = TRUE),
          deviation = .data$fitted - .data$mean_fitted,
          abs_deviation = abs(.data$deviation),
          direction = ifelse(.data$deviation > 0, "Above Mean", "Below Mean")
        ) |>
        dplyr::arrange(.data$fitted) |>
        dplyr::mutate(rank = dplyr::row_number())

      label_df <- get_extreme_labels(df, "deviation", top_n_labels)

      p1 <- ggplot2::ggplot(df, ggplot2::aes(x = .data$fitted)) +
        ggplot2::geom_density(fill = "gray80", alpha = 0.5) +
        ggplot2::geom_vline(ggplot2::aes(xintercept = .data$mean_fitted[1]), linetype = "dashed", color = "black") +
        ggplot2::geom_rug(data = label_df, color = "red", linewidth = 1) +
        ggplot2::labs(title = "Distribution of Expected Category Scores", x = NULL, y = "Density") +
        ggplot2::theme_minimal()

      p2 <- ggplot2::ggplot(df, ggplot2::aes(x = .data$rank, y = .data$fitted)) +
        ggplot2::geom_segment(ggplot2::aes(xend = .data$rank, yend = .data$mean_fitted), color = "gray60") +
        ggplot2::geom_point(ggplot2::aes(color = .data$direction, size = .data$abs_deviation)) +
        ggplot2::geom_hline(ggplot2::aes(yintercept = .data$mean_fitted[1]), linetype = "dashed") +
        ggrepel::geom_label_repel(data = label_df, ggplot2::aes(label = .data$id), size = 3, min.segment.length = 0) +
        ggplot2::scale_color_manual(values = c("Above Mean" = "#0072B2", "Below Mean" = "#D55E00")) +
        ggplot2::labs(x = x_label, y = "Expected Score", color = "Direction", size = "Deviation\nMagnitude") +
        ggplot2::theme_minimal()

      return(patchwork::wrap_plots(p1, p2, ncol = 1, heights = c(1, 2)))
    }
  }
}
