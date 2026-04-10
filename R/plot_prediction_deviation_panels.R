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
  strata_info <- NULL
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
    if (inherits(model, "polr") || inherits(model, "clm") || inherits(model, "ordinal")) {
      type <- "ordinal"
    } else if (inherits(model, "glm") && model$family$family %in% c("binomial", "quasibinomial")) {
      type <- "binomial"
    } else if (inherits(model, "glmerMod") && model@resp$family$family %in% c("binomial", "quasibinomial")) {
      type <- "binomial"
    } else {
      type <- "gaussian"
    }
  }

  get_extreme_labels <- function(df, metric_col, n) {
    df |> dplyr::arrange(dplyr::desc(abs(.data[[metric_col]]))) |> utils::head(n)
  }

  if (type == "gaussian") {
    # GAUSSIAN / LINEAR LOGIC
    # Approximate predict, some packages handle se.fit differently, so wrap safely
    preds <- tryCatch(predict(model, se.fit = TRUE), error = function(e) list(fit = predict(model), se.fit = rep(0, nrow(data))))
    if (is.numeric(preds)) preds <- list(fit = preds, se.fit = rep(0, nrow(data)))

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
      ggplot2::labs(x = x_label, y = "Fitted Value", color = "Direction", size = "Deviation\nMagnitude") +
      ggplot2::theme_minimal()

    return(patchwork::wrap_plots(p1, p2, ncol = 1, heights = c(1, 2)))

  } else if (type == "binomial") {
    # BINOMIAL / LOGISTIC LOGIC
    preds <- tryCatch(predict(model, type = "response", se.fit = TRUE),
                      error = function(e) list(fit = predict(model, type = "response"), se.fit = rep(0, nrow(data))))
    if (is.numeric(preds)) preds <- list(fit = preds, se.fit = rep(0, nrow(data)))

    # Try to extract response variable
    form <- tryCatch(formula(model), error = function(e) NULL)
    obs_outcome <- NULL
    if (!is.null(form)) {
      resp_name <- as.character(form)[2]
      if (resp_name %in% names(data)) {
        obs_outcome <- as.factor(data[[resp_name]])
      }
    }
    if (is.null(obs_outcome)) {
      obs_outcome <- factor(rep(NA, nrow(data)))
    }

    resids <- tryCatch(abs(residuals(model, type = "deviance")), error = function(e) rep(0, nrow(data)))

    df <- data |>
      dplyr::mutate(
        id = dplyr::row_number(),
        obs_outcome = obs_outcome,
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
      df <- df |> dplyr::mutate(
        wrong = factor(ifelse((.data$fitted > 0.5 & as.numeric(as.character(.data$obs_outcome)) == 0) |
                                (.data$fitted < 0.5 & as.numeric(as.character(.data$obs_outcome)) == 1),
                              "Wrong", "Correct"))
      )
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
        ggplot2::labs(x = x_label, y = "Predicted Probability", color = "Direction", size = "|Deviance\nResidual|")
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
    probs <- tryCatch(predict(model, type = "probs"), error = function(e) predict(model, type = "p"))
    if (!is.matrix(probs) && !is.data.frame(probs)) stop("Could not extract probability matrix from ordinal model.")

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

      k_seq <- seq_len(ncol(probs))
      df$expected_score <- rowSums(probs * matrix(k_seq, nrow = nrow(probs), ncol = ncol(probs), byrow = TRUE))

      # Probability of observed category
      df$observed_prob <- NA
      for (i in seq_len(nrow(df))) {
        col_idx <- match(df$obs_cat[i], colnames(probs))
        if (!is.na(col_idx)) {
          df$observed_prob[i] <- probs[i, col_idx]
        }
      }

      if ("stratum" %in% names(data)) {
        df$stratum <- data$stratum
        df <- df |>
          dplyr::group_by(.data$stratum) |>
          dplyr::summarize(
            dplyr::across(tidyselect::all_of(colnames(probs)), \(x) mean(x, na.rm = TRUE)),
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
        tidyr::pivot_longer(cols = tidyselect::all_of(colnames(probs)), names_to = "Category", values_to = "Probability") |>
        dplyr::mutate(Category = factor(.data$Category, levels = colnames(probs)))

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
      k_seq <- seq_len(ncol(probs))
      exp_scores <- rowSums(probs * matrix(k_seq, nrow = nrow(probs), ncol = ncol(probs), byrow = TRUE))

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
