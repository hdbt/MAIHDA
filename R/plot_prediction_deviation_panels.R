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

# Prediction weights aligned to `data`'s rows, used to make the per-stratum
# aggregation a weighted mean for weighted fits (consistent with the weighted VPC
# and the other stratum-level plots); for an aggregated-binomial fit each row is
# weighted by its binomial TRIAL count, matching the trial-weighting the raw-model
# fallback below already applies via weights(type = "prior"). Falls back to unit
# weights -- so the weighted means reduce EXACTLY to plain means -- when the model
# is unweighted, the weights cannot be recovered, or they do not align with `data`
# (e.g. user-supplied prediction data). These are prior/precision (and trial)
# weights, not a complex survey design (no design-based variance is computed).
maihda_prediction_panel_prior_weights <- function(maihda_obj, model, data) {
  n <- nrow(data)
  w <- NULL
  if (!is.null(maihda_obj)) {
    w <- tryCatch(maihda_prediction_weights(maihda_obj), error = function(e) NULL)
  }
  if (is.null(w)) {
    w <- tryCatch(stats::weights(model, type = "prior"), error = function(e) NULL)
  }
  if (is.null(w) || !is.numeric(w) || length(w) != n) {
    return(rep(1, n))
  }
  w <- as.numeric(w)
  w[!is.finite(w)] <- NA_real_
  w
}

maihda_prediction_panel_auto_type <- function(model) {
  if (inherits(model, "polr") || inherits(model, "clm") ||
      inherits(model, "clmm") || inherits(model, "ordinal")) {
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
  # Count models must predict on the response (count) scale: routing them through
  # the Gaussian branch would plot link-scale (log) predictions under Gaussian
  # labels and calculations.
  if (!is.null(fam_name) && fam_name %in% c("poisson", "quasipoisson", "negbinomial")) {
    return("poisson")
  }

  "gaussian"
}

maihda_prediction_panel_fitted <- function(model, data, type, fitted_data = FALSE) {
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

  # lme4: when predicting the model's OWN fitted rows (no external newdata), reuse the
  # stored linear predictor by calling predict() WITHOUT newdata, so any offset is
  # retained -- predicting with newdata = the stored model frame drops an external
  # offset= and errors on a formula offset() term (its raw variable lives only as the
  # frame's derived "offset(...)"/"(offset)" column). A genuine external newdata cannot
  # reconstruct an external offset, so reject it rather than return a silently wrong
  # prediction (mirroring predict_maihda()'s individual-prediction path).
  if (inherits(model, "merMod")) {
    if (isTRUE(fitted_data)) {
      fit <- if (type == "binomial" || type == "poisson") {
        stats::predict(model, type = "response")
      } else {
        stats::predict(model)
      }
      return(list(fit = as.numeric(fit), se.fit = rep(NA_real_, length(fit))))
    }
    if (maihda_mermod_has_external_offset(model)) {
      stop("This model was fit with an external offset (offset = ... passed to ",
           "fit_maihda()), which cannot be reconstructed for the supplied prediction ",
           "data. Refit with the offset written into the formula (e.g. ",
           "... + offset(log(exposure))), or omit 'data' to use the fitted rows.",
           call. = FALSE)
    }
  }

  # SE fallbacks below are NA_real_ -- not 0 -- for model classes whose predict()
  # does not implement se.fit (returning 0 produced ci_lower == ci_upper == fitted,
  # i.e. fake zero-width "95% CI" bars; NA instead propagates through fitted +/-
  # 1.96 * se and ggplot drops the geom_errorbar layer, honestly communicating "no
  # SE available"). Where predict.merMod DOES return an (approximate) se.fit it is
  # used only for case-level bars; the per-stratum panels never average these row
  # SEs into a stratum SE (that is not a valid SE of a weighted mean) -- see
  # maihda_prediction_panel_attach_ci() / maihda_prediction_stratum_interval_brms().
  if (type == "binomial" || type == "poisson") {
    # Count and binary models are summarised on the response scale (expected count
    # / probability), not the link scale that predict() returns by default for a
    # GLM(M).
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
  if (inherits(model, "clmm")) {
    # predict.clmm does not exist: rebuild the location eta = x'beta + u from
    # the stored components (fixed-effects-only $terms, $xlevels, $beta, and the
    # stratum conditional modes) and difference the cumulative probabilities.
    # Including the random effect matches the other branches of this panel,
    # whose predict() calls include random effects by default.
    maihda_require_ordinal()
    tt <- stats::delete.response(model$terms)
    mf <- stats::model.frame(tt, data, xlev = model$xlevels,
                             na.action = stats::na.pass)
    X <- stats::model.matrix(tt, mf)
    beta <- model$beta
    eta <- if (is.null(beta) || length(beta) == 0) {
      rep(0, nrow(data))
    } else {
      missing_cols <- setdiff(names(beta), colnames(X))
      if (length(missing_cols) > 0) {
        stop("Could not rebuild the clmm design matrix; missing column(s): ",
             paste(missing_cols, collapse = ", "), call. = FALSE)
      }
      drop(X[, names(beta), drop = FALSE] %*% beta)
    }
    # A formula offset term is part of the latent location clmm fits but is not a
    # column of X, so add it explicitly or the rebuilt probabilities omit it.
    off <- stats::model.offset(mf)
    if (!is.null(off)) {
      eta <- eta + off
    }
    re_list <- tryCatch(ordinal::ranef(model), error = function(e) NULL)
    if (!is.null(re_list) && "stratum" %in% names(re_list) &&
        "stratum" %in% names(data)) {
      tab <- re_list[["stratum"]]
      re_col <- intersect(c("(Intercept)", "Intercept"), colnames(tab))
      if (length(re_col) > 0) {
        u <- stats::setNames(as.numeric(tab[[re_col[1]]]), rownames(tab))
        u <- u[as.character(data$stratum)]
        u[is.na(u)] <- 0
        eta <- eta + unname(u)
      }
    }
    probs <- maihda_ordinal_category_probs(eta, model$alpha, model$link)
    return(as.data.frame(probs))
  }

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

# Correct per-stratum response-scale interval from posterior draws: aggregate the
# predictions WITHIN each stratum per draw (a weighted mean over the stratum's rows)
# and take posterior quantiles across draws. Unlike a weighted mean of the row-level
# SEs, this respects the fixed/random-effect uncertainty the rows in a stratum share
# (the SE of a weighted mean is a'Sigma a, not the mean of the component SEs), and it
# keeps the interval asymmetric (quantiles, not a symmetric pseudo-SE). `ep` is the
# ndraws x nobs posterior expected-value matrix on the response scale, `strat` the
# length-nobs stratum labels, `w` the length-nobs prior/precision weights. Returns a
# data frame keyed by sort(unique(strat)) with `fitted` (posterior mean of the
# stratum mean) and central `level` quantiles `ci_lower`/`ci_upper`. Pure/Stan-free.
maihda_stratum_interval_from_epred <- function(ep, strat, w = NULL, level = 0.95) {
  ep <- as.matrix(ep)
  strat <- as.character(strat)
  nobs <- ncol(ep)
  if (length(strat) != nobs) {
    stop("`strat` must have one label per column of `ep`.", call. = FALSE)
  }
  if (is.null(w) || length(w) != nobs) {
    w <- rep(1, nobs)
  }
  w <- as.numeric(w)
  a <- (1 - level) / 2
  groups <- sort(unique(strat))
  rows <- lapply(groups, function(g) {
    sel <- which(strat == g)
    ww <- w[sel]
    ww[!is.finite(ww) | ww < 0] <- 0
    if (sum(ww) <= 0) {
      ww <- rep(1, length(sel))  # degenerate weights -> equal weighting
    }
    # Per-draw weighted stratum mean, then posterior quantiles across draws.
    m <- as.vector(ep[, sel, drop = FALSE] %*% (ww / sum(ww)))
    q <- stats::quantile(m, c(a, 1 - a), names = FALSE, na.rm = TRUE)
    c(mean(m, na.rm = TRUE), q[1], q[2])
  })
  M <- do.call(rbind, rows)
  data.frame(stratum = groups, fitted = M[, 1], ci_lower = M[, 2],
             ci_upper = M[, 3], stringsAsFactors = FALSE)
}

# brms wrapper for maihda_stratum_interval_from_epred(): draws the response-scale
# posterior expected values on the plotting `data` and aggregates them within each
# stratum per draw. Returns NULL (interval omitted) for a non-brms model, an
# unavailable posterior_epred, or a non-2-D epred (e.g. an ordinal/categorical
# family), so the caller falls back to no stratum error bars rather than a
# statistically invalid averaged SE.
maihda_prediction_stratum_interval_brms <- function(model, data, weight = NULL,
                                                    level = 0.95) {
  if (!inherits(model, "brmsfit") || !requireNamespace("brms", quietly = TRUE)) {
    return(NULL)
  }
  if (is.null(data) || !"stratum" %in% names(data)) {
    return(NULL)
  }
  ep <- tryCatch(brms::posterior_epred(model, newdata = data),
                 error = function(e) NULL)
  if (is.null(ep) || length(dim(ep)) != 2L || ncol(ep) != nrow(data)) {
    return(NULL)
  }
  maihda_stratum_interval_from_epred(ep, data$stratum, weight, level = level)
}

# Attach `ci_lower`/`ci_upper` to a prediction-panel data frame, clamped to
# [lower, upper]. For a stratum-aggregated panel the interval comes from the
# draw-based `strat_ci` (brms) or is omitted as NA (frequentist) -- never from a
# weighted mean of row-level SEs. For a case-level panel it is the symmetric normal
# interval `fitted +/- 1.96 * se`, NA where the model exposes no se.fit (so ggplot
# drops those bars).
maihda_prediction_panel_attach_ci <- function(df, aggregated, strat_ci,
                                              lower = -Inf, upper = Inf) {
  n <- nrow(df)
  if (aggregated) {
    if (!is.null(strat_ci)) {
      idx <- match(as.character(df$stratum), as.character(strat_ci$stratum))
      ci_lower <- strat_ci$ci_lower[idx]
      ci_upper <- strat_ci$ci_upper[idx]
    } else {
      ci_lower <- rep(NA_real_, n)
      ci_upper <- rep(NA_real_, n)
    }
  } else {
    ci_lower <- df$fitted - 1.96 * df$se
    ci_upper <- df$fitted + 1.96 * df$se
  }
  df$ci_lower <- pmax(lower, pmin(upper, ci_lower))
  df$ci_upper <- pmax(lower, pmin(upper, ci_upper))
  df
}

#' Plot Prediction Deviation Panels
#'
#' @description Creates an advanced, publication-ready two-panel dashboard for
#' visualizing predicted values and highlighting the most notable cases or strata.
#' What "notable" means depends on the model type, and the labelled points are
#' \emph{not} statistical outliers in the regression-diagnostic sense:
#' \itemize{
#'   \item Gaussian and Poisson (and the ordinal \code{"expected_score"} mode):
#'     the cases/strata whose prediction sits furthest from the mean prediction
#'     (largest deviation), ranked by absolute deviation.
#'   \item Binomial: the cases/strata with the largest absolute deviance residual,
#'     i.e. where the observed 0/1 outcome is least consistent with the fitted
#'     probability (worst-fit points), ranked by \eqn{|deviance residual|}.
#'   \item Ordinal \code{"surprise"} mode: the cases/strata with the highest
#'     surprise \eqn{-\log P(\text{observed category})}, i.e. the least probable
#'     observations under the model.
#' }
#'
#' @param model A fitted model object (e.g., from `lm()`, `glm()`, `MASS::polr()`, or `lme4::glmer()`).
#' @param data The original data frame used to fit the model. If `NULL`, attempts to extract from the model.
#' @param type Model type: "auto" (default), "gaussian", "poisson", "binomial", or "ordinal".
#' @param ordinal_mode For ordinal models: "surprise" (default, based on observation probability) or "expected_score".
#' @param top_n_labels Number of points to label on the plot. The ranking metric
#'   depends on the model type (see Description): deviation from the mean
#'   prediction for Gaussian/Poisson and the ordinal expected-score mode, absolute
#'   deviance residual for binomial, and surprise for the ordinal surprise mode.
#'   Default is 5.
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
                                             type = c("auto", "gaussian", "poisson", "binomial", "ordinal"),
                                             ordinal_mode = c("surprise", "expected_score"),
                                             top_n_labels = 5,
                                             strata_info = NULL) {

  rlang::check_installed(c("ggplot2", "patchwork", "dplyr", "tidyr", "ggrepel"))

  type <- match.arg(type)
  ordinal_mode <- match.arg(ordinal_mode)

  # Whether the caller supplied external prediction data. When they did NOT, the
  # panels predict the model's own fitted rows, and those predictions must be taken
  # WITHOUT newdata so lme4 reuses its stored linear predictor (which includes any
  # offset); passing the stored model frame back as newdata drops an external offset=
  # and errors on a formula offset() term. See maihda_prediction_panel_fitted().
  data_supplied <- !is.null(data)

  # Check if model is a maihda_model. Keep the wrapper so prior/precision weights
  # can be recovered for the weighted stratum aggregation before unwrapping.
  maihda_obj <- NULL
  if (inherits(model, "maihda_model")) {
    maihda_obj <- model
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

  # Prior/precision weights for the per-stratum aggregation (unit weights for an
  # unweighted fit, so the weighted means below reduce to plain means).
  prior_w <- maihda_prediction_panel_prior_weights(maihda_obj, model, data)

  # Auto-detect model type if requested
  if (type == "auto") {
    type <- maihda_prediction_panel_auto_type(model)
  }

  get_extreme_labels <- function(df, metric_col, n) {
    df |> dplyr::arrange(dplyr::desc(abs(.data[[metric_col]]))) |> utils::head(n)
  }

  if (type == "gaussian" || type == "poisson") {
    # GAUSSIAN / LINEAR (and POISSON / COUNT) LOGIC. Both rank strata/cases by how
    # far their prediction sits from the mean prediction; counts are summarised on
    # the response (expected-count) scale with count labels, and the symmetric
    # interval is clamped at 0.
    is_count <- type == "poisson"
    preds <- maihda_prediction_panel_fitted(model, data, type,
                                            fitted_data = !data_supplied)

    value_dist_title <- if (is_count) "Distribution of Predicted Counts" else "Distribution of Fitted Values"
    value_axis_label <- if (is_count) "Predicted Count" else "Fitted Value"

    df <- data |>
      dplyr::mutate(
        id = dplyr::row_number(),
        fitted = preds$fit,
        se = preds$se.fit,
        weight = prior_w
      )

    aggregated <- "stratum" %in% names(df)
    strat_ci <- NULL
    if (aggregated) {
      # Aggregate only the point estimate. Row-level SEs are NOT averaged into a
      # stratum SE (the SE of a weighted mean is a'Sigma a, not the weighted mean
      # of the component SEs); the stratum interval is computed correctly from the
      # posterior draws for brms and omitted otherwise.
      strat_ci <- maihda_prediction_stratum_interval_brms(model, data, prior_w)
      df <- maihda_weighted_stratum_aggregate(df, "fitted")

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

    df <- maihda_prediction_panel_attach_ci(
      df, aggregated, strat_ci, lower = if (is_count) 0 else -Inf
    )

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
      ggplot2::labs(title = value_dist_title, x = NULL, y = "Density") +
      theme_maihda()

    p2 <- ggplot2::ggplot(df, ggplot2::aes(x = .data$rank, y = .data$fitted)) +
      ggplot2::geom_segment(ggplot2::aes(xend = .data$rank, yend = .data$mean_fitted), color = "gray60") +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$ci_lower, ymax = .data$ci_upper), width = 0, color = "gray50", alpha = 0.5) +
      ggplot2::geom_point(ggplot2::aes(color = .data$direction, size = .data$abs_deviation)) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = .data$mean_fitted[1]), linetype = "dashed") +
      ggrepel::geom_label_repel(data = label_df, ggplot2::aes(label = .data$id), size = 3, min.segment.length = 0) +
      ggplot2::scale_color_manual(values = c("Above Mean" = "#0072B2", "Below Mean" = "#D55E00")) +
      ggplot2::labs(
        x = x_label, y = value_axis_label, color = "Direction", size = "Deviation\nMagnitude"
      ) +
      theme_maihda()

    return(patchwork::wrap_plots(p1, p2, ncol = 1, heights = c(1, 2)))

  } else if (type == "binomial") {
    # BINOMIAL / LOGISTIC LOGIC
    preds <- maihda_prediction_panel_fitted(model, data, "binomial",
                                            fitted_data = !data_supplied)

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
        abs_res_dev = resids,
        weight = prior_w
      )

    is_aggregated <- "stratum" %in% names(df)
    strat_ci <- NULL

    if (is_aggregated) {
      # Aggregate the point estimate and the absolute deviance residual (a genuine
      # per-stratum mean diagnostic). Row-level SEs are NOT averaged into a stratum
      # SE; the stratum probability interval is computed correctly from the
      # posterior draws for brms and omitted otherwise.
      strat_ci <- maihda_prediction_stratum_interval_brms(model, data, prior_w)
      df <- maihda_weighted_stratum_aggregate(
        df, c("fitted", "abs_res_dev")
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

    df <- maihda_prediction_panel_attach_ci(
      df, is_aggregated, strat_ci, lower = 0, upper = 1
    )

    df <- df |>
      dplyr::mutate(
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
      theme_maihda()

    p2 <- ggplot2::ggplot(df, ggplot2::aes(x = .data$rank, y = .data$fitted)) +
      ggplot2::geom_segment(ggplot2::aes(xend = .data$rank, yend = .data$mean_fitted), color = "gray60", alpha = 0.5) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$ci_lower, ymax = .data$ci_upper), width = 0, color = "gray70", alpha = 0.3)

    if (is_aggregated) {
      p2 <- p2 +
        ggplot2::geom_point(ggplot2::aes(color = .data$direction, size = .data$abs_res_dev), alpha = 0.8) +
        ggplot2::labs(
          x = x_label, y = "Predicted Probability", color = "Direction", size = "|Deviance\nResidual|"
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
      theme_maihda()

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

      # Per-observation surprise (negative log-likelihood of the observed
      # category). The stratum-level value is the MEAN of this -- average surprise
      # / log loss = mean(-log(p)). Collapsing probabilities first and taking
      # -log(mean(p)) is a different (smaller, by Jensen) quantity that can change
      # the stratum ranking, so surprise is computed per row and then averaged.
      df$surprise <- -log(df$observed_prob)

      if ("stratum" %in% names(data)) {
        df$stratum <- data$stratum
        df$weight <- prior_w
        # Prior-weight-weighted per-stratum means of the category probabilities and
        # the surprise/score summaries (the stratum surprise stays the average of
        # the per-row -log P, now weighted); reduces to plain means when unweighted.
        df <- maihda_weighted_stratum_aggregate(
          df, c(prob_cols, "expected_score", "observed_prob", "surprise")
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

      # df$surprise is already the per-observation value (case-level) or its
      # per-stratum mean (stratum-level); do not recompute it from a collapsed
      # probability here.

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
        theme_maihda() +
        ggplot2::theme(axis.text.x = ggplot2::element_blank())

      p2 <- ggplot2::ggplot(df, ggplot2::aes(x = .data$rank, y = .data$surprise)) +
        ggplot2::geom_segment(ggplot2::aes(xend = .data$rank, yend = 0), color = "gray50") +
        ggplot2::geom_point(ggplot2::aes(color = .data$surprise, size = .data$surprise)) +
        ggplot2::scale_color_viridis_c(option = "inferno") +
        ggrepel::geom_label_repel(data = label_df, ggplot2::aes(label = .data$id), size = 3) +
        ggplot2::labs(x = x_label, y = "Surprise\n(-log(P(Observed)))", color = "Surprise", size = "Surprise") +
        theme_maihda()

      return(patchwork::wrap_plots(p1, p2, ncol = 1, heights = c(1, 2)))

    } else {
      # expected_score
      k_seq <- seq_len(ncol(prob_mat))
      exp_scores <- rowSums(prob_mat * matrix(k_seq, nrow = nrow(prob_mat), ncol = ncol(prob_mat), byrow = TRUE))

      df <- data |>
        dplyr::mutate(
          id = dplyr::row_number(),
          fitted = exp_scores,
          weight = prior_w
        )

      if ("stratum" %in% names(df)) {
        # Prior-weight-weighted per-stratum mean expected score; reduces to the
        # plain mean when the fit is unweighted.
        df <- maihda_weighted_stratum_aggregate(df, c("fitted"))

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
        theme_maihda()

      p2 <- ggplot2::ggplot(df, ggplot2::aes(x = .data$rank, y = .data$fitted)) +
        ggplot2::geom_segment(ggplot2::aes(xend = .data$rank, yend = .data$mean_fitted), color = "gray60") +
        ggplot2::geom_point(ggplot2::aes(color = .data$direction, size = .data$abs_deviation)) +
        ggplot2::geom_hline(ggplot2::aes(yintercept = .data$mean_fitted[1]), linetype = "dashed") +
        ggrepel::geom_label_repel(data = label_df, ggplot2::aes(label = .data$id), size = 3, min.segment.length = 0) +
        ggplot2::scale_color_manual(values = c("Above Mean" = "#0072B2", "Below Mean" = "#D55E00")) +
        ggplot2::labs(x = x_label, y = "Expected Score", color = "Direction", size = "Deviation\nMagnitude") +
        theme_maihda()

      return(patchwork::wrap_plots(p1, p2, ncol = 1, heights = c(1, 2)))
    }
  }
}
