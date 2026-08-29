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
#'     \item "predicted": Predicted values for each stratum with confidence
#'       intervals. By default the strata are ordered highest-to-lowest predicted
#'       value (a ranked caterpillar plot); see \code{order_by} to change or
#'       disable the ordering
#'     \item "upset": UpSet-style alternative to \code{"predicted"} -- an
#'       intersection-size bar, a category matrix encoding each stratum's level on
#'       every dimension, and the predicted-value panel, all sharing one column
#'       order. Replaces the long intersectional axis labels with the matrix.
#'       Binary 0/1 dimensions show as a single present/absent row; multi-level
#'       factors get one row per level. Columns are ordered by intersection size
#'       (largest first) by default; \code{order_by = "predicted_desc"} gives the
#'       ranked caterpillar of \code{"predicted"} with the matrix in place of the
#'       text labels
#'     \item "effect_decomp": Visualizes additive vs intersectional deviation from global mean
#'     \item "prediction_deviation": Detailed deviation panels for individuals or strata
#'     \item "context_vpc": Stratum vs. context variance bars for a contextual
#'       cross-classified fit (\code{fit_maihda(context = )}); errors otherwise
#'     \item "vpc_trajectory": Time-varying VPC/ICC curve for a \strong{longitudinal}
#'       fit (\code{fit_maihda(id =, time =)}); errors otherwise. For a longitudinal
#'       model \code{"vpc"} and \code{"all"} also route here
#'     \item "trajectories": Predicted per-stratum mean trajectories over time
#'       (longitudinal fits only)
#'     \item "all": Generate all available plots (default if not specified)
#'   }
#' @param summary_obj Optional maihda_summary object from \code{summary()}.
#'   If NULL, will be computed.
#' @param n_strata Maximum number of strata to display in the predicted plot.
#'   When there are more strata than this, the first \code{n_strata} (in stratum
#'   order) are shown and the plot caption notes how many were omitted. Default
#'   is 50. Use NULL for all strata.
#' @param highlight_interactions Highlight the strata that carry a credibly
#'   non-zero intersectional interaction (from \code{\link{maihda_interactions}})
#'   on the BLUP-based views (\code{"effect_decomp"}, \code{"predicted"},
#'   \code{"obs_vs_shrunken"}); other views ignore it. \code{FALSE} (default) off;
#'   \code{TRUE} computes the flags with \code{maihda_interactions()} defaults; or
#'   pass a multiple-testing method such as \code{"BH"} or a
#'   \code{maihda_interactions} object to reuse a specific \code{conf_level}/
#'   \code{adjust}. For the pure-interaction reading the model should be the
#'   adjusted (or crossed-dimensions) model -- e.g. via
#'   \code{plot()} on a \code{\link{maihda}} analysis, which routes these views to
#'   the adjusted model automatically. Which column of the screen drives the
#'   highlight set is governed by \code{highlight_by}.
#' @param highlight_by Which interaction-screen column defines the highlighted
#'   strata: \code{"flag"} (default) uses the zero-centred \code{flagged} column
#'   (credibly non-zero interaction), preserving the historical behaviour;
#'   \code{"rope"} uses the equivalence \code{decision} column, highlighting the
#'   strata classified \code{"relevant"} (interaction interval entirely outside the
#'   region of practical equivalence). \code{"rope"} requires a screen carrying a
#'   \code{decision} column: either pass \code{rope}, or supply a
#'   \code{maihda_interactions} object built with \code{rope}; otherwise it errors.
#' @param rope Equivalence region (a "smallest interaction of interest") forwarded
#'   to \code{\link{maihda_interactions}} when the screen is computed here (i.e.
#'   when \code{highlight_interactions} is \code{TRUE} or a \code{p.adjust} method
#'   name). \code{NULL} (default) adds no equivalence classification; a single
#'   positive \code{d} means the symmetric region \code{c(-d, d)} on the latent
#'   (link) scale, or supply \code{c(lower, upper)}. Needed for
#'   \code{highlight_by = "rope"} unless the supplied \code{maihda_interactions}
#'   object already carries a \code{decision} column. Ignored when a precomputed
#'   \code{maihda_interactions} object is passed (its own \code{rope} is used).
#' @param only_flagged Show \emph{only} the flagged strata rather than dimming the
#'   rest. \code{FALSE} (default) keeps every stratum (flagged ones highlighted);
#'   \code{TRUE} restricts the \code{"predicted"} and \code{"obs_vs_shrunken"}
#'   views to the highlighted strata (those carrying a credibly non-zero
#'   interaction, or -- under \code{highlight_by = "rope"} -- those classified
#'   ROPE-\code{"relevant"}), so a highlighted stratum is never hidden by the
#'   \code{n_strata} cap. When \code{TRUE} and
#'   \code{highlight_interactions} is left \code{FALSE}, the flags are computed
#'   with \code{\link{maihda_interactions}} defaults; pass
#'   \code{highlight_interactions} to choose the \code{conf_level}/\code{adjust}.
#'   A captioned empty panel is returned when no stratum is flagged. It does not
#'   apply to \code{"effect_decomp"} (whose waterfall exists to show each flagged
#'   stratum's place in the \emph{full} distribution); that view stays highlighted.
#'   Independently of this argument, whenever interactions are highlighted the
#'   \code{n_strata} cap on \code{"predicted"} becomes flag-aware: every flagged
#'   stratum is kept and the remaining slots are filled according to \code{select}.
#' @param select When the \code{n_strata} cap must drop strata, which to keep:
#'   \code{"order"} (default; the first n_strata in stratum order, the historical
#'   behaviour) or \code{"deviation"} (the n_strata furthest from the reference
#'   line -- largest \code{|predicted - reference|}, so the most extreme strata in
#'   \emph{both} directions). Applies to \code{"predicted"} and, for a longitudinal
#'   fit, \code{"trajectories"} (where it keeps the strata whose trajectories swing
#'   furthest from the population curve). Flagged strata are always kept; this
#'   governs the fill and the unflagged case. \code{select} changes \emph{which}
#'   strata appear; their left-to-right display order is a separate choice governed
#'   by \code{order_by}.
#' @param order_by For \code{type = "predicted"} and \code{type = "upset"}, how to
#'   order the strata that are displayed (\strong{display-only} -- it does not
#'   change \emph{which} strata are shown, that is \code{n_strata}/\code{select},
#'   nor the predicted values, intervals, reference line, or highlighted set):
#'   \code{"predicted_desc"} orders from the highest predicted value to the lowest,
#'   \code{"stratum"} keeps the native stratum order, \code{"predicted_asc"} orders
#'   from lowest to highest, \code{"deviation"} orders by largest absolute deviation
#'   from the reference line (\code{|predicted - reference|}), and \code{"size"}
#'   orders by intersection (stratum) size, largest first. The default is
#'   \strong{per view}: \code{"predicted_desc"} for \code{"predicted"} (labels run
#'   from the highest predicted value down) and \code{"size"} for \code{"upset"}
#'   (the UpSet convention, which also makes its intersection-size bar monotone).
#'   On the \code{"upset"} view the three value-based orders sort on the quantity
#'   the bottom panel actually shows, so they follow \code{quantity} -- with
#'   \code{quantity = "interaction"} they order by the stratum random effect and
#'   \code{"deviation"} measures the distance from zero. Combining
#'   \code{order_by = "predicted_desc"} with \code{type = "upset"} gives the
#'   ranked caterpillar of the \code{"predicted"} view drawn against the UpSet
#'   category matrix instead of long text labels. Ignored by the other plot types.
#' @param quantity For \code{type = "upset"}, which quantity the bottom panel
#'   shows: \code{"predicted"} (default) the stratum's predicted value (fixed +
#'   random effect) against the across-strata reference line, or
#'   \code{"interaction"} the stratum random effect (the BLUP) against zero --
#'   the deviation from the model's fixed prediction, which is the \emph{pure}
#'   intersectional interaction when the dimension main effects are in the model
#'   (the adjusted model). Ignored by the other plot types.
#' @param ... Additional arguments (not currently used).
#'
#' @return For a single \code{type}, a \pkg{ggplot2} object that you can extend
#'   with the usual \code{+} grammar (themes, \code{\link[ggplot2]{labs}()},
#'   added layers, or a replacement fill/colour scale). Some types return a
#'   richer object: \code{"prediction_deviation"} returns a \pkg{patchwork} of
#'   two panels and \code{"upset"} a \pkg{patchwork} of three panels (theme
#'   every panel at once with \code{& theme_*()}).
#'   \code{type = "all"} returns a named list of ggplot2 objects.
#'
#' @examples
#' \donttest{
#' strata_result <- make_strata(maihda_sim_data, vars = c("gender", "race"))
#' model <- fit_maihda(health_outcome ~ age + (1 | stratum), data = strata_result$data)
#'
#' # VPC plot
#' plot(model, type = "vpc")
#'
#' # Single-type plots are ggplot objects -- restyle them with ggplot2:
#' plot(model, type = "vpc") +
#'   ggplot2::theme_classic() +
#'   ggplot2::labs(title = "Variance partition, restyled")
#'
#' # Generate all plots (a named list); pick one out to restyle it:
#' plots <- plot(model)
#' plots$predicted + ggplot2::theme_bw()
#' }
#'
#' @export
#' @import ggplot2
#' @importFrom dplyr arrange
plot.maihda_model <- function(x, type = c("all", "vpc", "obs_vs_shrunken", "predicted", "upset", "effect_decomp", "prediction_deviation", "context_vpc", "vpc_trajectory", "trajectories"),
                       summary_obj = NULL, n_strata = 50, highlight_interactions = FALSE, only_flagged = FALSE, highlight_by = c("flag", "rope"), rope = NULL, select = c("order", "deviation"), order_by = c("predicted_desc", "stratum", "predicted_asc", "deviation", "size"), quantity = c("predicted", "interaction"), ...) {
  if (!inherits(x, "maihda_model")) {
    stop("'x' must be a maihda_model object from fit_maihda()")
  }

  object <- x

  # Restricting the view to flagged strata needs flags to restrict by; if the
  # caller asked for only_flagged but no highlight, fall back to the default
  # interaction screen so only_flagged works on its own.
  only_flagged <- maihda_validate_only_flagged(only_flagged)
  if (only_flagged && isFALSE(highlight_interactions)) {
    highlight_interactions <- TRUE
  }
  # Which column of the interaction screen defines the highlighted set: the
  # zero-centred flag ("flag") or the ROPE equivalence decision ("rope").
  highlight_by <- match.arg(highlight_by)
  # Which strata survive the n_strata cap on the predicted / trajectory views.
  select <- match.arg(select)
  # How the stratum views order the strata they display (display-only; does not
  # change which strata are shown -- that is n_strata / select). The two views
  # have DIFFERENT sensible defaults (highest predicted first for "predicted";
  # largest intersection first for "upset", the UpSet convention its size bar
  # relies on), so an unsupplied order_by forwards as NULL -- "let the view
  # choose" -- rather than as this formal's first choice. NULL also arrives from
  # plot.maihda_analysis(), which forwards the same sentinel; a supplied value is
  # still validated here so a typo errors at the entry point.
  order_by_supplied <- !missing(order_by) && !is.null(order_by)
  order_by <- if (order_by_supplied) match.arg(order_by) else NULL
  # Which quantity the upset view's estimate panel shows: the predicted value or
  # the stratum random effect (interaction).
  quantity <- match.arg(quantity)


  if (missing(type)) {
    type <- "all"
  } else {
    type <- match.arg(type)
  }

  # Get summary if not provided
  if (is.null(summary_obj)) {
    summary_obj <- summary(object)
  }

  # Resolve the set of strata to highlight (NULL = no highlight). The BLUP-based
  # views (effect_decomp / predicted / obs_vs_shrunken) mark them; other views
  # ignore it. highlight_by chooses the screen column: the zero-centred flag, or
  # the ROPE "relevant" classification (which needs a `rope`).
  highlight_ids <- maihda_resolve_highlight(object, highlight_interactions,
                                            highlight_by = highlight_by, rope = rope)

  # Longitudinal (growth-curve) models: the time-varying VPC and the stratum mean
  # trajectories replace the cross-sectional VPC bar (whose single proportion stack
  # is undefined when the between-stratum variance varies with time). type = "vpc"
  # redirects to the trajectory, and "all" yields the two trajectory views.
  if (!is.null(object$longitudinal_info)) {
    if (type == "all") {
      plots <- list(
        vpc_trajectory = plot_vpc_trajectory(summary_obj),
        trajectories = tryCatch(
          plot_stratum_trajectories(object, summary_obj, n_strata, select = select),
          error = function(e) NULL)
      )
      for (p in plots[!vapply(plots, is.null, logical(1))]) print(p)
      return(invisible(plots))
    }
    if (type %in% c("vpc", "vpc_trajectory")) {
      return(plot_vpc_trajectory(summary_obj))
    }
    if (type == "trajectories") {
      return(plot_stratum_trajectories(object, summary_obj, n_strata, select = select))
    }
    # Every remaining view (predicted, obs_vs_shrunken, effect_decomp,
    # prediction_deviation) is a cross-sectional BLUP scalar per stratum, which
    # misrepresents a growth model's trajectory estimand. Refuse them and point
    # to the trajectory views above.
    maihda_stop_longitudinal_scalar(
      paste0("plot(type = \"", type, "\")"),
      stratum_slope = maihda_object_stratum_slope(object))
  } else if (type %in% c("vpc_trajectory", "trajectories")) {
    stop("type = \"", type, "\" is only available for a longitudinal MAIHDA ",
         "(fit_maihda(id = , time = )).", call. = FALSE)
  }

  if (type == "all") {
    plots <- list()

    plots$vpc <- plot_vpc(summary_obj)

    # Try obs_vs_shrunken. A panel that cannot be built for this model degrades
    # to NULL (the montage still prints the rest) but warns with the underlying
    # reason, so a silently missing panel is never mistaken for a complete set.
    if ("stratum" %in% names(object$data)) {
      plots$obs_vs_shrunken <- maihda_try_optional(plot_obs_vs_shrunken(object, summary_obj, highlight = highlight_ids, only_flagged = only_flagged), "Plot panel 'obs_vs_shrunken'")
    }

    plots$predicted <- maihda_try_optional(plot_predicted_strata(object, summary_obj, n_strata, highlight = highlight_ids, only_flagged = only_flagged, select = select, order_by = order_by), "Plot panel 'predicted'")

    top_n_labels <- if (is.null(n_strata)) 10 else min(10, n_strata)
    plots$effect_decomp <- maihda_try_optional(plot_effect_decomposition(object, summary_obj, top_n_labels, highlight = highlight_ids), "Plot panel 'effect_decomp'")

    plots$prediction_deviation <- maihda_try_optional(plot_prediction_deviation_panels(object, type = "auto"), "Plot panel 'prediction_deviation'")

    if (!is.null(object$context_info)) {
      plots$context_vpc <- maihda_try_optional(plot_context_vpc(summary_obj), "Plot panel 'context_vpc'")
    }

    # print them
    for (p in plots[!sapply(plots, is.null)]) { print(p) }
    return(invisible(plots))
  } else {
    if (type == "vpc") {
      plot <- plot_vpc(summary_obj)
    } else if (type == "context_vpc") {
      plot <- plot_context_vpc(summary_obj)
    } else if (type == "obs_vs_shrunken") {
      plot <- plot_obs_vs_shrunken(object, summary_obj, highlight = highlight_ids, only_flagged = only_flagged)
    } else if (type == "predicted") {
      plot <- plot_predicted_strata(object, summary_obj, n_strata, highlight = highlight_ids, only_flagged = only_flagged, select = select, order_by = order_by)
    } else if (type == "upset") {
      plot <- plot_upset_strata(object, summary_obj, n_strata, highlight = highlight_ids, only_flagged = only_flagged, select = select, order_by = order_by, quantity = quantity)
    } else if (type == "effect_decomp") {
      # The waterfall's value IS the full-distribution context, so filtering it
      # away defeats the view; keep it highlighted and say so rather than no-op.
      if (only_flagged) {
        message("plot(): 'only_flagged' does not apply to type = \"effect_decomp\" -- ",
                "its waterfall exists to show each flagged stratum's place in the full ",
                "distribution, so the flagged strata are highlighted in context instead.")
      }
      top_n_labels <- if (is.null(n_strata)) 10 else min(10, n_strata)
      plot <- plot_effect_decomposition(object, summary_obj, top_n_labels, highlight = highlight_ids)
    } else if (type == "prediction_deviation") {
      plot <- plot_prediction_deviation_panels(object, type = "auto")
    }

    return(plot)
  }
}

# Resolve the `highlight_interactions` plot argument to a character vector of
# highlighted stratum ids (or NULL = no highlight). Accepts FALSE/NULL (off), TRUE
# (compute the screen with maihda_interactions() defaults), a p.adjust method name
# such as "BH" (compute the screen with that adjustment), or a precomputed
# maihda_interactions object (so callers can set conf_level/adjust/rope once and
# reuse). `highlight_by` selects which screen column defines the set: the
# zero-centred `flagged` column ("flag") or the equivalence `decision` column
# ("rope", the strata classified "relevant"). `rope` is forwarded to
# maihda_interactions() when the screen is computed here so the decision column
# exists; it is ignored when a precomputed object is supplied (that object's own
# rope governs its decision column).
maihda_resolve_highlight <- function(model, highlight_interactions,
                                     highlight_by = c("flag", "rope"), rope = NULL) {
  highlight_by <- match.arg(highlight_by)
  if (is.null(highlight_interactions) || isFALSE(highlight_interactions)) {
    return(NULL)
  }
  flags <- if (inherits(highlight_interactions, "maihda_interactions")) {
    highlight_interactions
  } else if (isTRUE(highlight_interactions)) {
    maihda_interactions(model, rope = rope)
  } else if (is.character(highlight_interactions) && length(highlight_interactions) == 1L) {
    choices <- c("none", stats::p.adjust.methods)
    if (!highlight_interactions %in% choices) {
      stop("'highlight_interactions' must be FALSE, TRUE, a multiple-comparison ",
           "method name (e.g. \"BH\"), or a maihda_interactions object from ",
           "maihda_interactions().", call. = FALSE)
    }
    maihda_interactions(model, adjust = highlight_interactions, rope = rope)
  } else {
    stop("'highlight_interactions' must be FALSE, TRUE, a multiple-comparison ",
         "method name (e.g. \"BH\"), or a maihda_interactions object from ",
         "maihda_interactions().", call. = FALSE)
  }
  ids <- maihda_highlight_ids(flags, highlight_by)
  # Carry the screen's parameters along so downstream views can caption an
  # only_flagged subset honestly (e.g. "95% interval, BH-adjusted", or the ROPE
  # region for highlight_by = "rope"). Attributes ride harmlessly through the
  # `as.character(stratum) %in% highlight` membership checks that consume `ids`
  # elsewhere.
  attr(ids, "conf_level") <- attr(flags, "conf_level")
  attr(ids, "adjust") <- attr(flags, "adjust")
  attr(ids, "engine") <- attr(flags, "engine")
  attr(ids, "highlight_by") <- highlight_by
  attr(ids, "rope") <- attr(flags, "rope")
  ids
}

# Extract the highlighted stratum ids from a resolved interaction screen per
# `highlight_by`: the zero-centred `flagged` column ("flag", the historical
# default) or the equivalence `decision` column ("rope", the strata classified
# "relevant"). The ROPE path needs a screen that was computed with a `rope`, so it
# errors with an actionable message when the `decision` column is absent.
maihda_highlight_ids <- function(flags, highlight_by) {
  if (identical(highlight_by, "rope")) {
    if (!"decision" %in% names(flags)) {
      stop("highlight_by = \"rope\" highlights the ROPE-relevant strata, but the ",
           "interaction screen has no 'decision' column (it was computed without a ",
           "'rope'). Pass 'rope' -- a positive half-width d for the region ",
           "c(-d, d), or c(lower, upper) on the latent (link) scale -- e.g. ",
           "plot(fit, highlight_interactions = TRUE, highlight_by = \"rope\", ",
           "rope = 0.4); or supply a maihda_interactions object built with ",
           "maihda_interactions(fit, rope = ...).", call. = FALSE)
    }
    # NA decisions (interval not computable) are never "relevant"; %in% drops them.
    as.character(flags$stratum[flags$decision %in% "relevant"])
  } else {
    as.character(flags$stratum[flags$flagged %in% TRUE])
  }
}

# The noun naming the highlighted set, for honest captions/labels: "flagged"
# (zero-centred screen) or "ROPE-relevant" (equivalence screen). Reads the
# highlight_by attribute attached by maihda_resolve_highlight(); defaults to
# "flagged" when absent (a bare id vector), preserving the historical wording.
maihda_highlight_noun <- function(highlight) {
  if (identical(attr(highlight, "highlight_by"), "rope")) "ROPE-relevant" else "flagged"
}

# Validate the `only_flagged` plot argument: NULL -> FALSE, otherwise a single
# non-NA logical. Anything else is a usage error.
maihda_validate_only_flagged <- function(only_flagged) {
  if (is.null(only_flagged)) return(FALSE)
  if (!is.logical(only_flagged) || length(only_flagged) != 1L || is.na(only_flagged)) {
    stop("'only_flagged' must be a single TRUE or FALSE.", call. = FALSE)
  }
  only_flagged
}

# Human-readable description of the interaction screen behind a highlight set, for
# an honest only_flagged caption. Mirrors the basis line of
# maihda_print_interactions_line(): a credible interval for brms, otherwise the
# conf_level interval with the multiplicity stance actually used. For a ROPE
# highlight (highlight_by = "rope") it instead names the equivalence region the
# "relevant" classification is read against. Reads the attributes attached by
# maihda_resolve_highlight(); defaults to a 95% interval when they are absent
# (e.g. a bare character vector of ids).
maihda_highlight_screen_label <- function(highlight) {
  conf <- attr(highlight, "conf_level")
  conf_pct <- if (is.null(conf)) 95 else conf * 100

  # ROPE highlight: the "relevant" decision is "conf_level interval entirely
  # outside the region", so caption the region (symmetric -> |effect| > d) and the
  # interval level, on the latent (link) scale the interaction lives on.
  if (identical(attr(highlight, "highlight_by"), "rope")) {
    rope <- attr(highlight, "rope")
    region <- if (is.null(rope)) {
      "the equivalence region"
    } else if (isTRUE(all.equal(rope[1], -rope[2]))) {
      sprintf("|effect| > %g", rope[2])
    } else {
      sprintf("interval outside [%g, %g]", rope[1], rope[2])
    }
    return(sprintf("%s on the latent (link) scale, %.0f%% interval",
                   region, conf_pct))
  }

  adjust <- attr(highlight, "adjust")
  engine <- attr(highlight, "engine")
  if (identical(engine, "brms")) {
    sprintf("%.0f%% credible interval", conf_pct)
  } else if (is.null(adjust) || identical(adjust, "none")) {
    sprintf("%.0f%% interval, unadjusted", conf_pct)
  } else {
    sprintf("%.0f%% interval, %s-adjusted", conf_pct, adjust)
  }
}

# Placeholder for an only_flagged view when no stratum is highlighted: an empty
# panel carrying the same title plus an explanatory annotation, so the filtered
# view degrades gracefully instead of erroring or drawing a bare axis. Mirrors the
# print method's "No strata show interaction credibly different from zero". The
# lead phrase tracks the highlight mode (flag vs ROPE) read from `highlight`.
maihda_no_flagged_plot <- function(title, screen_label, highlight = NULL) {
  body <- if (identical(attr(highlight, "highlight_by"), "rope")) {
    paste0("No strata classified as ROPE-relevant\n(", screen_label, ").")
  } else {
    paste0("No strata flagged as carrying a credibly\n",
           "non-zero interaction (", screen_label, ").")
  }
  ggplot2::ggplot() +
    ggplot2::annotate(
      "text", x = 0, y = 0, size = 4, color = "grey30",
      label = body) +
    ggplot2::scale_x_continuous(limits = c(-1, 1)) +
    ggplot2::scale_y_continuous(limits = c(-1, 1)) +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    theme_maihda() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank())
}

# Choose which `k` of the candidate row indices `cand_idx` to keep when a capped
# view must drop strata, honouring the selection rule. "order" keeps the first k
# (stratum order, the historical behaviour); "deviation" keeps the k most extreme
# by |deviation| -- magnitude, so BOTH tails are pulled in, ranked by distance --
# where `mag` is a per-row magnitude aligned to the full table (non-finite -> -Inf
# so incomputable rows sort last). The returned indices are sorted ascending, so
# the DISPLAYED x-axis stays in stratum order regardless of how the survivors were
# chosen: selection and display order are deliberately separate.
maihda_pick_strata <- function(cand_idx, k, select, mag) {
  if (length(cand_idx) <= k) return(sort(cand_idx))
  chosen <- if (identical(select, "deviation")) {
    m <- mag[cand_idx]
    m[!is.finite(m)] <- -Inf
    cand_idx[order(m, decreasing = TRUE)[seq_len(k)]]
  } else {
    utils::head(cand_idx, k)
  }
  sort(chosen)
}

# Append a star to the strata flagged for highlighting, for use in plot labels.
maihda_highlight_label <- function(label, stratum, highlight) {
  flagged <- as.character(stratum) %in% highlight
  ifelse(flagged, paste0(label, " *"), as.character(label))
}

# Discrete aesthetics for the focus-by-contrast highlight: non-flagged strata
# dimmed (neutral grey, low opacity), flagged strata solid in the accent colour.
# Named by the logical flag so they map directly onto an aes(colour/alpha =
# .maihda_flag). Replaces the old open-circle "ring" overlay, which added geometry
# on top of already-busy plots and, on the effect-decomposition view, sat at the
# total deviation rather than at the interaction it flagged.
maihda_highlight_palette <- function() c(`FALSE` = "#9AA0A6", `TRUE` = "#D55E00")
maihda_highlight_alpha   <- function() c(`FALSE` = 0.30, `TRUE` = 1.00)

#' VPC Visualization Plot
#'
#' @param summary_obj A maihda_summary object
#' @return A ggplot2 object
#' @keywords internal
#' @import ggplot2
plot_vpc <- function(summary_obj) {
  vc <- summary_obj$variance_components
  vpc_data <- vc[vc$component != "Total", , drop = FALSE]
  is_cc <- identical(attr(vc, "kind"), "cross_classified")
  is_ctx <- identical(attr(vc, "kind"), "contextual")

  if (is_cc) {
    # Crossed-dimensions split: one slice per dimension (additive), one for the
    # interaction, any contextual random intercepts, then the residual. Colour the
    # additive dimensions from a qualitative palette, the interaction in orange (it
    # is the "new between"), contexts in green, the residual in blue. The component
    # order in the table drives the stack.
    add_comps <- vpc_data$component[grepl("^Additive: ", vpc_data$component)]
    ctx_comps <- vpc_data$component[grepl("^Context: ", vpc_data$component)]
    dim_palette <- c("#CC79A7", "#009E73", "#0072B2", "#D55E00", "#117733",
                     "#882255", "#44AA99", "#332288")
    component_colors <- stats::setNames(rep("#999999", nrow(vpc_data)),
                                        vpc_data$component)
    if (length(add_comps) > 0) {
      component_colors[add_comps] <-
        dim_palette[((seq_along(add_comps) - 1) %% length(dim_palette)) + 1]
    }
    if (length(ctx_comps) > 0) {
      ctx_palette <- c("#117733", "#44AA99", "#999933", "#DDCC77")
      component_colors[ctx_comps] <-
        ctx_palette[((seq_along(ctx_comps) - 1) %% length(ctx_palette)) + 1]
    }
    component_colors["Intersectional interaction"] <- "#E69F00"
    component_colors["Within-stratum (residual)"] <- "#56B4E9"
    plot_title <- sprintf("Variance Partition (crossed-dimensions), VPC/ICC = %.3f",
                          summary_obj$vpc$estimate)
    # Keep the table's component ordering (additive dims, interaction, residual).
    vpc_data$component <- factor(vpc_data$component, levels = vpc_data$component)
  } else if (is_ctx) {
    # Contextual cross-classified split: the between-stratum (intersectional)
    # slice, one slice per context (the general contextual effects), any other
    # random effects, and the residual. Stratum keeps the canonical orange so the
    # plot reads like the single-stratum VPC bar with the context broken out.
    component_colors <- maihda_vpc_component_colors(vpc_data$component)
    plot_title <- sprintf("Variance Partition (stratum x context), VPC/ICC = %.3f",
                          summary_obj$vpc$estimate)
    vpc_data$component <- factor(vpc_data$component, levels = vpc_data$component)
  } else {
    component_colors <- maihda_vpc_component_colors(vpc_data$component)
    plot_title <- sprintf("Variance Partition Coefficient (VPC/ICC) = %.3f",
                          summary_obj$vpc$estimate)
  }

  # Create plot
  p <- ggplot(vpc_data, aes(x = "", y = .data$proportion, fill = .data$component)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_flip() +
    scale_fill_manual(values = component_colors) +
    labs(
      title = plot_title,
      x = "",
      y = "Proportion of Variance",
      fill = "Component"
    ) +
    theme_maihda() +
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

# Map a set of variance-partition component names to fill colours, shared by the
# single stacked VPC bar (plot_vpc, its standard and contextual branches) and the
# null-vs-adjusted change bar (plot_vpc_change) so the two views never drift.
# Between-stratum stays the canonical orange, each "Context: <var>" slice draws from
# the green context palette, the other-random-effects and residual slices keep their
# hues, and anything unmapped falls back to grey.
maihda_vpc_component_colors <- function(components) {
  components <- as.character(components)
  colors <- stats::setNames(rep("#999999", length(components)), components)
  ctx <- components[grepl("^Context: ", components)]
  if (length(ctx) > 0) {
    ctx_palette <- c("#117733", "#44AA99", "#999933", "#DDCC77")
    colors[ctx] <- ctx_palette[((seq_along(ctx) - 1) %% length(ctx_palette)) + 1]
  }
  fixed <- c(
    "Between-stratum (random)"  = "#E69F00",
    "Other random effects"      = "#009E73",
    "Within-stratum (residual)" = "#56B4E9"
  )
  present <- intersect(names(fixed), components)
  colors[present] <- fixed[present]
  colors
}

#' Null-vs-Adjusted Variance-Partition Change Plot
#'
#' Draws the null and adjusted models' variance partitions as two stacked bars on a
#' shared axis and annotates the PCV, so the drop in the between-stratum share after
#' accounting for the dimensions' additive main effects is visible in one figure.
#' Backs \code{plot(<maihda_analysis>, type = "vpc", model = "both")} for a
#' cross-sectional analysis.
#'
#' @param null_summary,adjusted_summary The \code{maihda_summary} objects of the
#'   null and adjusted models of a \code{\link{maihda}} analysis.
#' @param pcv Optional \code{pcv_result} (the analysis's \code{$pcv}); its estimate
#'   is shown in the subtitle. \code{NULL} or a non-finite estimate drops the PCV
#'   clause (e.g. a boundary fit where the PCV is undefined).
#' @return A ggplot2 object.
#' @keywords internal
#' @import ggplot2
plot_vpc_change <- function(null_summary, adjusted_summary, pcv = NULL) {
  one <- function(s, label) {
    vc <- s$variance_components
    d <- vc[vc$component != "Total", c("component", "proportion"), drop = FALSE]
    d$model <- label
    d
  }
  dat <- rbind(one(null_summary, "Null model"),
               one(adjusted_summary, "Adjusted model"))

  # Component fill order: the null model's component order first (between-stratum
  # -> ... -> residual), then any adjusted-only components, so the stack and legend
  # read consistently across the two bars.
  comp_levels <- unique(as.character(dat$component))
  dat$component <- factor(dat$component, levels = comp_levels)
  # After coord_flip() the first factor level sits at the bottom, so order the
  # models adjusted-then-null to put the null bar on top (natural null -> adjusted
  # top-to-bottom reading).
  dat$model <- factor(dat$model, levels = c("Adjusted model", "Null model"))

  component_colors <- maihda_vpc_component_colors(comp_levels)

  # Subtitle: the between-stratum VPC/ICC shift, then the PCV (the additive share of
  # the between-stratum variance) when it is available.
  null_vpc <- null_summary$vpc$estimate
  adj_vpc  <- adjusted_summary$vpc$estimate
  subtitle <- sprintf("Between-stratum VPC/ICC: %.1f%% (null) \u2192 %.1f%% (adjusted)",
                      null_vpc * 100, adj_vpc * 100)
  pcv_val <- if (!is.null(pcv)) pcv$pcv else NULL
  if (!is.null(pcv_val) && is.finite(pcv_val)) {
    subtitle <- paste0(subtitle, sprintf("  |  PCV = %.1f%%", pcv_val * 100))
  }

  ggplot(dat, aes(x = .data$model, y = .data$proportion, fill = .data$component)) +
    geom_bar(stat = "identity", width = 0.7, color = "white") +
    coord_flip() +
    scale_fill_manual(values = component_colors) +
    geom_text(aes(label = sprintf("%.1f%%", .data$proportion * 100)),
              position = position_stack(vjust = 0.5),
              color = "white", fontface = "bold", size = 4) +
    labs(
      title = "Variance Partition: Null vs. Adjusted Model",
      subtitle = subtitle,
      x = NULL,
      y = "Proportion of Variance",
      fill = "Component",
      caption = "PCV: proportional change in the between-stratum variance (null -> adjusted)."
    ) +
    theme_maihda() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      panel.grid = element_blank()
    )
}

#' Null-vs-Adjusted VPC-over-time Change Plot (longitudinal MAIHDA)
#'
#' Overlays the null and adjusted models' time-varying VPC/ICC curves on one axis,
#' so the reduction in the between-stratum share after adjustment is visible over
#' the whole time range. Backs \code{plot(<maihda_analysis>, type = "vpc",
#' model = "both")} for a longitudinal analysis; complements
#' \code{\link{plot_pcv_trajectory}} (the additive-share PCV(t)).
#'
#' @param null_summary,adjusted_summary Longitudinal \code{maihda_summary} objects.
#' @return A ggplot2 object.
#' @keywords internal
#' @import ggplot2
plot_vpc_trajectory_change <- function(null_summary, adjusted_summary) {
  lng_n <- null_summary$longitudinal
  lng_a <- adjusted_summary$longitudinal
  if (is.null(lng_n) || is.null(lng_a)) {
    stop("plot_vpc_trajectory_change() needs two longitudinal summaries.",
         call. = FALSE)
  }
  vt_n <- lng_n$vpc_t; vt_n$model <- "Null model"
  vt_a <- lng_a$vpc_t; vt_a$model <- "Adjusted model"
  dat <- rbind(vt_n[, c("time", "estimate", "model")],
               vt_a[, c("time", "estimate", "model")])
  dat$model <- factor(dat$model, levels = c("Null model", "Adjusted model"))

  ggplot(dat, aes(x = .data$time, y = .data$estimate, color = .data$model)) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 2) +
    geom_vline(xintercept = lng_n$ref_time, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c("Null model" = "#E69F00",
                                  "Adjusted model" = "#0072B2")) +
    labs(
      title = "Time-varying VPC/ICC: Null vs. Adjusted Model",
      subtitle = sprintf("Dashed line: reference time %s = %g",
                         lng_n$time, lng_n$ref_time),
      x = lng_n$time,
      y = "VPC/ICC (between-stratum share of variance)",
      color = "Model"
    ) +
    theme_maihda() +
    theme(plot.title = element_text(face = "bold"))
}

#' Stratum vs. Context Variance Plot (contextual cross-classified MAIHDA)
#'
#' One bar per variance component -- the between-stratum (intersectional)
#' variance, each context's variance, any other random effects, and the residual
#' -- on the variance scale, with each component's share of the total printed
#' above its bar. Complements \code{plot_vpc()}'s stacked proportion bar by
#' showing the \emph{magnitudes} the shares are computed from.
#'
#' @param summary_obj A \code{maihda_summary} from a contextual
#'   cross-classified fit (\code{fit_maihda(context = )}).
#' @return A ggplot2 object.
#' @keywords internal
#' @import ggplot2
plot_context_vpc <- function(summary_obj) {
  vc <- summary_obj$variance_components
  if (!identical(attr(vc, "kind"), "contextual") || is.null(summary_obj$context)) {
    stop("No contextual partition is available. Fit the model with ",
         "fit_maihda(context = ) (or maihda(context = )) to plot the stratum ",
         "vs. context variances.", call. = FALSE)
  }

  bar_data <- vc[vc$component != "Total", , drop = FALSE]
  bar_data$component <- factor(bar_data$component, levels = bar_data$component)

  ctx_comps <- levels(bar_data$component)[grepl("^Context: ", levels(bar_data$component))]
  component_colors <- stats::setNames(rep("#999999", nrow(bar_data)),
                                      levels(bar_data$component))
  component_colors["Between-stratum (random)"] <- "#E69F00"
  if (length(ctx_comps) > 0) {
    ctx_palette <- c("#117733", "#44AA99", "#999933", "#DDCC77")
    component_colors[ctx_comps] <-
      ctx_palette[((seq_along(ctx_comps) - 1) %% length(ctx_palette)) + 1]
  }
  component_colors["Other random effects"] <- "#009E73"
  component_colors["Within-stratum (residual)"] <- "#56B4E9"

  caption <- paste(
    "Contextual cross-classified MAIHDA: individuals are cross-classified by their",
    "intersectional stratum and the higher-level context(s).",
    "The between-stratum variance is conditional on the context random effect(s);",
    "the context variance is the between-context component of unexplained variance.",
    sep = "\n")

  ggplot(bar_data, aes(x = .data$component, y = .data$variance,
                       fill = .data$component)) +
    geom_col(color = "white") +
    geom_text(aes(label = sprintf("%.1f%%", .data$proportion * 100)),
              vjust = -0.4, fontface = "bold", size = 4) +
    scale_fill_manual(values = component_colors, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(
      title = sprintf("Stratum vs. Context Variance (VPC/ICC = %.3f)",
                      summary_obj$vpc$estimate),
      x = NULL,
      y = "Variance",
      caption = caption
    ) +
    theme_maihda() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 20, hjust = 1),
      panel.grid.major.x = element_blank()
    )
}

#' Observed vs. Shrunken Estimates Plot
#'
#' @details The x-axis is each stratum's raw observed mean; the y-axis is the
#'   model-based stratum estimate, which includes the fixed-effect contribution.
#'   For an intercept-only (null) model the vertical distance from the diagonal is
#'   pure shrinkage toward the grand mean. For a covariate-adjusted model the model
#'   estimate also moves with the stratum's covariate profile, so distance from the
#'   diagonal reflects \emph{both} shrinkage and covariate adjustment and should
#'   not be read as shrinkage alone.
#'
#' @param object A maihda_model object
#' @param summary_obj A maihda_summary object
#' @param highlight Optional character vector of highlighted stratum ids (flagged,
#'   or ROPE-relevant under \code{highlight_by = "rope"}), with the
#'   interaction-screen parameters attached as attributes.
#' @param only_flagged When TRUE, show only the highlighted strata; a captioned
#'   empty panel is returned if none are.
#' @return A ggplot2 object
#' @keywords internal
#' @import ggplot2
#' @importFrom dplyr group_by summarise
#' @importFrom stats formula terms
plot_obs_vs_shrunken <- function(object, summary_obj, highlight = NULL, only_flagged = FALSE) {
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
    } else if (object$engine == "wemix") {
      maihda_stratum_predictions_wemix(object, summary_obj, scale = "response")
    } else if (object$engine == "ordinal") {
      # Response scale = expected category score, matching the observed mean
      # category score computed above for an ordered-factor outcome.
      maihda_stratum_predictions_ordinal(object, summary_obj, scale = "response")
    } else {
      stop("Unsupported engine: ", object$engine)
    }

    plot_data <- merge(obs_means, stratum_est, by = "stratum")
    pred_idx <- match(as.character(plot_data$stratum), as.character(pred_data$stratum))
    plot_data$shrunken <- pred_data$predicted_row[pred_idx]
    plot_data$.maihda_flag <- as.character(plot_data$stratum) %in% highlight

    # The y = x diagonal is the only reference and is filtering-invariant, so an
    # only_flagged subset is safe. With nothing flagged, degrade to a captioned
    # empty panel rather than a bare diagonal.
    caption_txt <- NULL
    if (only_flagged) {
      n_total <- nrow(plot_data)
      n_flagged <- sum(plot_data$.maihda_flag)
      screen_label <- maihda_highlight_screen_label(highlight)
      noun <- maihda_highlight_noun(highlight)
      if (n_flagged == 0) {
        return(maihda_no_flagged_plot("Observed vs. Shrunken Stratum Estimates",
                                      screen_label, highlight))
      }
      plot_data <- plot_data[plot_data$.maihda_flag, , drop = FALSE]
      caption_txt <- sprintf("Showing the %d %s strata (%s) of %d total.",
                             n_flagged, noun, screen_label, n_total)
    }

    has_hl <- any(plot_data$.maihda_flag)

    # Create plot. When interactions are highlighted, focus by contrast -- flagged
    # strata solid in the accent colour, the rest dimmed -- instead of ringing them.
    point_layer <- if (has_hl) {
      geom_point(aes(size = .data$n, color = .data$.maihda_flag,
                     alpha = .data$.maihda_flag))
    } else {
      geom_point(aes(size = .data$n), alpha = 0.6, color = "#0072B2")
    }
    p <- ggplot(plot_data, aes(x = .data$observed, y = .data$shrunken)) +
      point_layer +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
      labs(
        title = "Observed vs. Shrunken Stratum Estimates",
        x = "Observed Stratum Mean",
        y = "Shrunken Estimate (with Random Effect)",
        size = "Sample Size",
        caption = caption_txt
      ) +
      theme_maihda() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.caption = element_text(hjust = 0.5, face = "italic", size = 9),
        legend.position = "right"
      )
    if (has_hl) {
      p <- p +
        scale_color_manual(values = maihda_highlight_palette(), guide = "none") +
        scale_alpha_manual(values = maihda_highlight_alpha(), guide = "none")
    }

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
    is_cumulative <- !is.null(fam_name) &&
      maihda_normalize_family_name(fam_name) == "cumulative"
    if (is_cumulative) {
      # Cumulative (ordinal) outcome: the observed value is the category score
      # (1..K in level order), whose stratum mean is the observed counterpart
      # of the model's expected category score.
      return(maihda_observed_plot_values(as.integer(x)))
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
#' @param highlight Optional character vector of highlighted stratum ids -- the
#'   flagged strata, or the ROPE-relevant strata under
#'   \code{highlight_by = "rope"} -- with the interaction-screen parameters
#'   attached as attributes (including \code{highlight_by} and \code{rope}).
#' @param only_flagged When TRUE, show only the highlighted strata (those in
#'   \code{highlight}); a captioned empty panel is returned if none are.
#' @param select When the \code{n_strata} cap drops strata, which to keep:
#'   \code{"order"} (default) the first n_strata in stratum order, or
#'   \code{"deviation"} the n_strata furthest from the reference line (largest
#'   \code{|predicted - reference|}, so both tails). Flagged strata are kept
#'   regardless; this governs the fill / the unflagged case. It controls
#'   \emph{which} strata are shown, separately from how they are ordered for
#'   display.
#' @param order_by Display order of the shown strata (\strong{display-only}; does
#'   not change which strata are shown -- that is \code{n_strata}/\code{select} --
#'   nor the predicted values, intervals, or reference line): \code{"predicted_desc"}
#'   (default) highest predicted at the top, \code{"stratum"} native stratum order,
#'   \code{"predicted_asc"} lowest at the top, \code{"deviation"} largest
#'   \code{|predicted - reference|} at the top, or \code{"size"} largest stratum
#'   first. \code{NULL} means the caller expressed no preference and takes this
#'   view's default (\code{"predicted_desc"}).
#' @return A ggplot2 object
#' @keywords internal
#' @import ggplot2
#' @importFrom dplyr arrange slice
plot_predicted_strata <- function(object, summary_obj, n_strata, scale = c("response", "link"), highlight = NULL, only_flagged = FALSE, select = c("order", "deviation"), order_by = c("predicted_desc", "stratum", "predicted_asc", "deviation", "size")) {
  # NULL is the "no explicit order asked for" sentinel the dispatchers forward;
  # this view's default is the ranked caterpillar (highest predicted at the top).
  order_by <- if (is.null(order_by)) "predicted_desc" else match.arg(order_by)
  prep <- maihda_prepare_predicted_strata(
    object, summary_obj, n_strata, scale = scale,
    highlight = highlight, only_flagged = only_flagged, select = select)
  if (isTRUE(prep$no_flagged)) {
    return(maihda_no_flagged_plot(
      "Predicted Subgroup Values with Conditional 95% Intervals", prep$screen_label,
      highlight))
  }
  stratum_est <- prep$stratum_est
  fixed_reference <- prep$fixed_reference
  caption_txt <- prep$caption_txt

  # Display order of the kept strata. This is purely cosmetic: the strata were
  # already selected (n_strata / select) and the reference line fixed from ALL
  # strata, so reordering here changes neither the shown set nor any value. order()
  # sinks any NA predictions to the bottom. Because the factor levels below are
  # reversed and the view is coord_flip()ped, the first row renders at the top --
  # so decreasing = TRUE puts the largest value on top.
  ord <- switch(order_by,
    stratum        = seq_len(nrow(stratum_est)),
    predicted_desc = order(stratum_est$predicted, decreasing = TRUE),
    predicted_asc  = order(stratum_est$predicted, decreasing = FALSE),
    deviation      = order(abs(stratum_est$predicted - fixed_reference), decreasing = TRUE),
    size           = order(stratum_est$n, decreasing = TRUE))
  stratum_est <- stratum_est[ord, , drop = FALSE]

  # Create factor to preserve order for plotting. Levels are reversed so that
  # after coord_flip() the first stratum sits at the top of the axis (natural
  # top-to-bottom reading order) rather than the bottom.
  stratum_est$display_label <- factor(stratum_est$display_label, levels = rev(stratum_est$display_label))

  has_hl <- any(stratum_est$.maihda_flag)

  # Create plot. Highlighted: flagged strata solid in the accent colour, the rest
  # dimmed (focus by contrast rather than ringing); flagged labels are starred.
  pt_layers <- if (has_hl) {
    list(
      geom_point(aes(color = .data$.maihda_flag, alpha = .data$.maihda_flag), size = 2),
      geom_errorbar(aes(ymin = .data$lower, ymax = .data$upper,
                        color = .data$.maihda_flag, alpha = .data$.maihda_flag),
                    width = 0.2)
    )
  } else {
    list(
      geom_point(size = 2, color = "#0072B2"),
      geom_errorbar(aes(ymin = .data$lower, ymax = .data$upper),
                    width = 0.2, alpha = 0.5, color = "#0072B2")
    )
  }
  p <- ggplot(stratum_est, aes(x = .data$display_label, y = .data$predicted)) +
    pt_layers +
    geom_hline(yintercept = fixed_reference, linetype = "dashed", color = "red", alpha = 0.7) +
    labs(
      title = "Predicted Subgroup Values with Conditional 95% Intervals",
      x = "Stratum",
      y = "Predicted Value",
      caption = caption_txt
    ) +
    coord_flip() +
    theme_maihda() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.caption = element_text(hjust = 0.5, face = "italic", size = 9),
      panel.grid.minor = element_blank()
    )
  if (has_hl) {
    p <- p +
      scale_color_manual(values = maihda_highlight_palette(), guide = "none") +
      scale_alpha_manual(values = maihda_highlight_alpha(), guide = "none")
  }

  return(p)
}

# Shared data preparation for the two stratum-prediction views: the text
# `predicted` plot and the `upset` composite. Computes each stratum's predicted
# value + conditional interval, the across-strata reference (the dashed line),
# the flag/highlight membership, and applies the only_flagged filter / n_strata
# cap / `select` rule -- everything both views need before they diverge on
# layout. Returns a list with `stratum_est` (the kept strata, carrying
# predicted/lower/upper/n/.maihda_flag/display_label), `fixed_reference`,
# `caption_txt`, `screen_label`, and `no_flagged` (TRUE when only_flagged found
# nothing flagged, so the caller returns its own titled empty panel).
maihda_prepare_predicted_strata <- function(object, summary_obj, n_strata,
                                            scale = c("response", "link"),
                                            highlight = NULL, only_flagged = FALSE,
                                            select = c("order", "deviation")) {
  scale <- match.arg(scale)
  select <- match.arg(select)

  pred_data <- if (object$engine == "lme4") {
    maihda_stratum_predictions_lme4(object, summary_obj, scale = scale)
  } else if (object$engine == "brms") {
    maihda_stratum_predictions_brms(object, summary_obj, scale = scale)
  } else if (object$engine == "wemix") {
    maihda_stratum_predictions_wemix(object, summary_obj, scale = scale)
  } else if (object$engine == "ordinal") {
    maihda_stratum_predictions_ordinal(object, summary_obj, scale = scale)
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
  # Carry the stratum size through for the upset view's intersection-size bars;
  # the text view ignores it.
  stratum_est$n <- pred_data$n[pred_idx]

  # Mark strata flagged as carrying a credibly non-zero interaction BEFORE any
  # truncation, so the cap can be made flag-aware (a flagged stratum past the cap
  # must not be silently dropped). The reference line above is already computed
  # from ALL strata, so neither filtering nor capping below shifts it.
  stratum_est$.maihda_flag <- as.character(stratum_est$stratum) %in% highlight
  n_total_strata <- nrow(stratum_est)
  n_flagged_total <- sum(stratum_est$.maihda_flag)
  screen_label <- maihda_highlight_screen_label(highlight)
  # Noun for the highlighted set in captions: "flagged" or "ROPE-relevant".
  noun <- maihda_highlight_noun(highlight)

  # Per-stratum magnitude = visual distance from the reference line. select =
  # "deviation" keeps the most extreme by this; "order" ignores it. The reference
  # is computed from ALL strata above, so the ranking is stable under capping.
  dev_mag <- abs(stratum_est$predicted - fixed_reference)
  # Phrase a kept subset of size k to match the rule actually used, so a
  # deviation-selected view is never read as a stratum-order one. Word order
  # differs ("first 2" vs "2 most extreme"), so build the whole phrase here.
  sel_phrase <- function(k) {
    if (identical(select, "deviation")) sprintf("%d most extreme", k) else sprintf("first %d", k)
  }

  caption_txt <- ""
  if (only_flagged) {
    # Restrict to flagged strata. With none flagged, signal the caller to return
    # its own titled empty panel rather than erroring or drawing a bare axis.
    if (n_flagged_total == 0) {
      return(list(no_flagged = TRUE, screen_label = screen_label))
    }
    flagged_idx <- which(stratum_est$.maihda_flag)
    # A cap still applies for readability when MANY strata are flagged; which
    # flagged strata survive it follows `select`.
    k <- if (is.null(n_strata)) length(flagged_idx) else n_strata
    keep_idx <- maihda_pick_strata(flagged_idx, k, select, dev_mag)
    capped <- length(keep_idx) < n_flagged_total
    stratum_est <- stratum_est[keep_idx, , drop = FALSE]
    caption_txt <- if (capped) {
      sprintf("\nShowing the %s of %d %s strata (%s); %d strata total.",
              sel_phrase(length(keep_idx)), n_flagged_total, noun, screen_label, n_total_strata)
    } else {
      sprintf("\nShowing the %d %s strata (%s) of %d total.",
              n_flagged_total, noun, screen_label, n_total_strata)
    }
  } else if (!is.null(n_strata) && n_total_strata > n_strata) {
    # Cap exceeded. n_strata is a MAXIMUM: keep the first n_strata in stratum order
    # (select = "order") or the n_strata furthest from the reference (select =
    # "deviation") -- never an evenly-spaced stride, which silently dropped strata
    # from the middle while implying full coverage.
    if (n_flagged_total > 0) {
      # Flag-aware cap: keep EVERY flagged stratum, then fill the remaining slots
      # per `select`. A flagged stratum beyond the cap is the signal the highlight
      # exists to surface, so it is never dropped -- even if that means showing
      # more than n_strata.
      flagged_idx <- which(stratum_est$.maihda_flag)
      nonflag_idx <- which(!stratum_est$.maihda_flag)
      n_fill <- max(0, n_strata - length(flagged_idx))
      fill_idx <- maihda_pick_strata(nonflag_idx, n_fill, select, dev_mag)
      keep_idx <- sort(union(flagged_idx, fill_idx))
      stratum_est <- stratum_est[keep_idx, , drop = FALSE]
      caption_txt <- if (length(flagged_idx) > n_strata) {
        sprintf(paste0("\nShowing all %d %s strata of %d (n_strata = %d ",
                       "exceeded to keep every %s stratum)."),
                length(flagged_idx), noun, n_total_strata, n_strata, noun)
      } else {
        sprintf(paste0("\nShowing %d of %d strata: all %d %s plus the %s ",
                       "others (n_strata = %d)."),
                length(keep_idx), n_total_strata, length(flagged_idx), noun,
                sel_phrase(n_fill), n_strata)
      }
    } else {
      keep_idx <- maihda_pick_strata(seq_len(n_total_strata), n_strata, select, dev_mag)
      stratum_est <- stratum_est[keep_idx, , drop = FALSE]
      caption_txt <- if (identical(select, "deviation")) {
        sprintf("\nShowing the %d strata furthest from the reference, of %d (n_strata = %d).",
                n_strata, n_total_strata, n_strata)
      } else {
        sprintf("\nShowing the first %d of %d strata (n_strata = %d).",
                n_strata, n_total_strata, n_strata)
      }
    }
  }

  # Use labels if available, otherwise use numeric stratum IDs
  if ("label" %in% names(stratum_est) && !all(is.na(stratum_est$label))) {
    # Use the meaningful labels for the x-axis
    stratum_est$display_label <- stratum_est$label
  } else {
    # Fall back to stratum IDs
    stratum_est$display_label <- stratum_est$stratum
  }

  # Star the flagged strata's axis labels so the highlight survives the (possibly
  # truncated) view. In only_flagged mode every shown stratum is flagged, so a
  # star on each would be redundant noise -- skip it there.
  if (!only_flagged && any(stratum_est$.maihda_flag)) {
    stratum_est$display_label <- maihda_highlight_label(
      stratum_est$display_label, stratum_est$stratum, highlight)
  }

  list(
    stratum_est = stratum_est,
    fixed_reference = fixed_reference,
    caption_txt = caption_txt,
    screen_label = screen_label,
    no_flagged = FALSE
  )
}

#' UpSet-style Predicted Stratum Plot
#'
#' Composite alternative to the text-labelled \code{"predicted"} view that
#' replaces the long intersectional axis labels with an UpSet-style category
#' matrix. Three panels share one column order: a top bar of intersection
#' (stratum) sizes, a middle matrix encoding each stratum's category on every
#' dimension, and a bottom panel of predicted values with conditional intervals.
#' Columns are ordered by intersection size (largest first) by default, or by any
#' other \code{order_by} rule -- \code{"predicted_desc"} turns the view into the
#' ranked caterpillar of \code{plot_predicted_strata()} drawn against the category
#' matrix instead of long text labels. Binary 0/1 (or logical) dimensions collapse
#' to a single present/absent row; multi-level factors get one row per level, and
#' each column lights exactly one dot per dimension.
#'
#' @inheritParams plot_predicted_strata
#' @param order_by Left-to-right column order (\strong{display-only}; does not
#'   change which strata are shown, nor any value): \code{"size"} (default,
#'   largest intersection first -- the UpSet convention), \code{"stratum"} native
#'   stratum order, or \code{"predicted_desc"} / \code{"predicted_asc"} /
#'   \code{"deviation"}, which sort on the quantity the estimate panel actually
#'   shows and so follow \code{quantity}. \code{NULL} takes this view's default.
#' @return A \pkg{patchwork} object stacking the three panels.
#' @keywords internal
#' @import ggplot2
plot_upset_strata <- function(object, summary_obj, n_strata, scale = c("response", "link"),
                              highlight = NULL, only_flagged = FALSE,
                              select = c("order", "deviation"),
                              order_by = c("size", "predicted_desc", "stratum",
                                           "predicted_asc", "deviation"),
                              quantity = c("predicted", "interaction")) {
  scale <- match.arg(scale)
  select <- match.arg(select)
  # NULL is the "no explicit order asked for" sentinel the dispatchers forward;
  # this view's default is the UpSet convention, largest intersection first.
  order_by <- if (is.null(order_by)) "size" else match.arg(order_by)
  quantity <- match.arg(quantity)
  is_interaction <- quantity == "interaction"
  plot_title <- if (is_interaction) {
    "Stratum Random Effects by Intersection"
  } else {
    "Predicted Subgroup Values by Intersection"
  }

  # The matrix encodes each stratum's category on every dimension, so it needs
  # the per-dimension stratum table from make_strata(); refuse a model that
  # carries only a bare stratum id.
  strata_vars <- object$strata_vars
  strata_info <- object$strata_info
  if (is.null(strata_vars) || length(strata_vars) == 0 || is.null(strata_info) ||
      !all(strata_vars %in% names(strata_info))) {
    stop("type = \"upset\" needs the per-dimension stratum table from ",
         "make_strata(); this model does not carry one. Use type = \"predicted\".",
         call. = FALSE)
  }

  prep <- maihda_prepare_predicted_strata(
    object, summary_obj, n_strata, scale = scale,
    highlight = highlight, only_flagged = only_flagged, select = select)
  if (isTRUE(prep$no_flagged)) {
    return(maihda_no_flagged_plot(plot_title, prep$screen_label, highlight))
  }
  stratum_est <- prep$stratum_est
  fixed_reference <- prep$fixed_reference
  caption_txt <- prep$caption_txt

  # The interaction view plots the stratum random effect (the BLUP) against a
  # zero line -- the deviation from the model's fixed prediction, which is the
  # pure interaction only when the dimension main effects are in the model. It
  # needs the random-effect interval the summary attaches.
  if (is_interaction &&
      !all(c("random_effect", "lower_95", "upper_95") %in% names(stratum_est))) {
    stop("quantity = \"interaction\" needs the stratum random-effect interval, ",
         "which this fit does not expose. Use quantity = \"predicted\".",
         call. = FALSE)
  }

  # Which column the estimate panel shows: the predicted value against the
  # across-strata reference, or the random effect (interaction) against zero.
  # Resolved HERE, above the ordering, because the value-based orders sort on the
  # quantity actually plotted -- ordering by `predicted` while drawing the random
  # effect would render as an unsorted panel.
  est_y   <- if (is_interaction) "random_effect" else "predicted"
  est_lo  <- if (is_interaction) "lower_95" else "lower"
  est_hi  <- if (is_interaction) "upper_95" else "upper"
  est_ref <- if (is_interaction) 0 else fixed_reference
  est_ylab <- if (is_interaction) "Stratum random effect" else "Predicted Value"

  # Left-to-right order of the kept strata. Purely cosmetic: `select` already
  # chose WHICH strata survive the cap and the reference line was fixed from ALL
  # strata, so this changes neither the shown set nor any value. The default is
  # intersection size (largest first) -- the UpSet convention, the most useful
  # here since the largest strata carry the most reliably estimated effects, and
  # the only order under which the top size bar decreases monotonically.
  # order() sinks any NA estimates to the right-hand end.
  ord <- switch(order_by,
    size           = order(stratum_est$n, decreasing = TRUE),
    stratum        = seq_len(nrow(stratum_est)),
    predicted_desc = order(stratum_est[[est_y]], decreasing = TRUE),
    predicted_asc  = order(stratum_est[[est_y]], decreasing = FALSE),
    deviation      = order(abs(stratum_est[[est_y]] - est_ref), decreasing = TRUE))
  stratum_est <- stratum_est[ord, , drop = FALSE]
  stratum_est$rank <- seq_len(nrow(stratum_est))
  k <- nrow(stratum_est)
  info_idx <- match(as.character(stratum_est$stratum), as.character(strata_info$stratum))

  # Lay out the matrix rows: a binary 0/1 (or logical) dimension is a single
  # present/absent row, while a multi-level factor expands to one row per level.
  # Each row records the level it lights up for (`on_level`) and a display label
  # (the variable name for an indicator, "var: level" for a factor level). Rows
  # stack top-to-bottom in `strata_vars` order; within a factor, levels keep
  # their natural order.
  row_specs <- list()
  for (v in strata_vars) {
    col <- strata_info[[v]]
    if (maihda_dim_is_indicator(col)) {
      row_specs[[length(row_specs) + 1L]] <- list(
        var = v, on_level = maihda_dim_levels(col)[2], label = v)
    } else {
      for (lev in maihda_dim_levels(col)) {
        row_specs[[length(row_specs) + 1L]] <- list(
          var = v, on_level = lev, label = paste0(v, ": ", lev))
      }
    }
  }
  n_rows <- length(row_specs)
  # Descending y so the first row sits at the top of the panel.
  for (i in seq_len(n_rows)) row_specs[[i]]$y <- n_rows - i + 1L
  y_breaks <- vapply(row_specs, function(s) s$y, numeric(1))
  y_labels <- vapply(row_specs, function(s) s$label, character(1))

  # Long matrix: one entry per (stratum, matrix-row). A cell is "on" when the
  # stratum's value on that dimension equals the row's level, so each column
  # lights exactly one dot per dimension (or none, for an absent indicator).
  mat <- do.call(rbind, lapply(row_specs, function(s) {
    vals <- strata_info[[s$var]][info_idx]
    data.frame(
      rank = stratum_est$rank,
      dim_y = s$y,
      on = as.character(vals) == as.character(s$on_level),
      stringsAsFactors = FALSE
    )
  }))

  # Vertical connector spanning the "on" dots within each column (the UpSet line).
  on_mat <- mat[mat$on, , drop = FALSE]
  seg <- do.call(rbind, lapply(split(on_mat, on_mat$rank), function(d) {
    if (nrow(d) < 2) return(NULL)
    data.frame(rank = d$rank[1], y0 = min(d$dim_y), y1 = max(d$dim_y))
  }))

  on_col <- "#2B2D42"
  off_col <- "#D9D9D9"
  x_scale <- scale_x_continuous(limits = c(0.5, k + 0.5), expand = c(0, 0))
  base_theme <- theme_maihda() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(2, 6, 2, 6)
    )

  # Panel 1: intersection (stratum) sizes.
  p_bar <- ggplot(stratum_est, aes(x = .data$rank, y = .data$n)) +
    geom_col(fill = on_col, width = 0.7) +
    x_scale +
    scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
    labs(x = NULL, y = "Intersection\nsize") +
    base_theme

  # Panel 2: the dot matrix -- this replaces the long text axis labels. Segments
  # are drawn first so the dots sit on top of the connector.
  matrix_layers <- list(
    geom_point(data = mat,
               aes(x = .data$rank, y = .data$dim_y, color = .data$on), size = 2.3)
  )
  if (!is.null(seg)) {
    matrix_layers <- c(
      list(geom_segment(data = seg,
                        aes(x = .data$rank, xend = .data$rank,
                            y = .data$y0, yend = .data$y1),
                        color = on_col, linewidth = 0.7)),
      matrix_layers)
  }
  p_matrix <- ggplot() +
    matrix_layers +
    scale_color_manual(values = c(`TRUE` = on_col, `FALSE` = off_col), guide = "none") +
    x_scale +
    scale_y_continuous(breaks = y_breaks, labels = y_labels,
                       limits = c(0.5, n_rows + 0.5), expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    base_theme +
    theme(panel.grid.major.y = element_blank())

  # Panel 3: the per-stratum estimate (inherits the highlight). The column shown
  # was resolved from `quantity` above, alongside the display order.
  has_hl <- any(stratum_est$.maihda_flag)
  est_layers <- if (has_hl) {
    list(
      geom_point(aes(color = .data$.maihda_flag, alpha = .data$.maihda_flag), size = 2),
      geom_errorbar(aes(ymin = .data[[est_lo]], ymax = .data[[est_hi]],
                        color = .data$.maihda_flag, alpha = .data$.maihda_flag),
                    width = 0.25)
    )
  } else {
    list(
      geom_point(size = 2, color = "#0072B2"),
      geom_errorbar(aes(ymin = .data[[est_lo]], ymax = .data[[est_hi]]),
                    width = 0.25, alpha = 0.5, color = "#0072B2")
    )
  }
  # The caption states the order actually drawn -- it is no longer fixed, and a
  # figure that misreports its own column order is worse than one that omits it.
  est_word <- if (is_interaction) "stratum random effect" else "predicted value"
  order_note <- switch(order_by,
    size           = "intersection size (largest first)",
    stratum        = "stratum order",
    predicted_desc = paste0(est_word, " (highest first)"),
    predicted_asc  = paste0(est_word, " (lowest first)"),
    deviation      = paste0("absolute deviation from ",
                            if (is_interaction) "zero" else "the reference",
                            " (largest first)"))
  note <- paste0("Dark dot = the stratum's category on each dimension; columns ",
                 "ordered by ", order_note, ".")
  if (is_interaction) {
    note <- paste0(note, "\nLower panel: stratum random effect on the link scale ",
                   "(deviation from the model's fixed prediction; the pure ",
                   "interaction when the dimension main effects are in the model).")
  }
  p_est <- ggplot(stratum_est, aes(x = .data$rank, y = .data[[est_y]])) +
    geom_hline(yintercept = est_ref, linetype = "dashed", color = "red", alpha = 0.7) +
    est_layers +
    x_scale +
    labs(x = NULL, y = est_ylab, caption = paste0(note, caption_txt)) +
    base_theme +
    theme(plot.caption = element_text(hjust = 0.5, face = "italic", size = 9))
  if (has_hl) {
    p_est <- p_est +
      scale_color_manual(values = maihda_highlight_palette(), guide = "none") +
      scale_alpha_manual(values = maihda_highlight_alpha(), guide = "none")
  }

  patchwork::wrap_plots(p_bar, p_matrix, p_est, ncol = 1,
                        heights = c(2.2, 0.45 * n_rows + 0.4, 3.0)) +
    patchwork::plot_annotation(
      title = plot_title,
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))
}

#' Recommended Figure Size for the UpSet Stratum Plot
#'
#' Computes sensible \code{width} and \code{height} (in inches) for
#' \code{plot(object, type = "upset")}, so a knitr chunk or
#' \code{\link[ggplot2]{ggsave}()} call can size the figure to its content. The
#' UpSet composite grows \emph{taller} with the number of matrix rows (one per
#' binary 0/1 dimension, one per level of a multi-level factor) and \emph{wider}
#' with the number of strata columns shown, so a single fixed size tends to crop
#' or stretch it -- particularly for multi-level designs (many rows) or a large
#' \code{n_strata} (many columns; UpSet is an inherently wide format).
#'
#' @param object A \code{maihda_model} from \code{\link{fit_maihda}} or a
#'   \code{maihda} analysis from \code{\link{maihda}}; it must carry the
#'   per-dimension stratum table from \code{\link{make_strata}}.
#' @param n_strata Maximum number of strata the plot will show -- pass the same
#'   value you give \code{plot()}. \code{NULL} means all strata. Default 50.
#' @return A list with numeric \code{width} and \code{height} (inches) plus the
#'   \code{rows} (matrix rows) and \code{cols} (strata shown) they derive from.
#' @examples
#' \donttest{
#' strata <- make_strata(maihda_sim_data, vars = c("gender", "race"))
#' model <- fit_maihda(health_outcome ~ (1 | stratum), data = strata$data)
#' sz <- maihda_upset_size(model, n_strata = 30)
#' ggplot2::ggsave(
#'   tempfile(fileext = ".png"),
#'   plot(model, type = "upset", n_strata = 30),
#'   width = sz$width, height = sz$height)
#' }
#' @seealso \code{\link{plot.maihda_model}}
#' @export
maihda_upset_size <- function(object, n_strata = 50) {
  model <- if (inherits(object, "maihda_analysis")) object$model else object
  strata_vars <- model$strata_vars
  strata_info <- model$strata_info
  if (is.null(strata_vars) || length(strata_vars) == 0 || is.null(strata_info) ||
      !all(strata_vars %in% names(strata_info))) {
    stop("maihda_upset_size() needs the per-dimension stratum table from ",
         "make_strata(); this object does not carry one.", call. = FALSE)
  }
  if (!is.null(n_strata) &&
      (!is.numeric(n_strata) || length(n_strata) != 1 ||
       is.na(n_strata) || n_strata < 1)) {
    stop("'n_strata' must be a single positive number or NULL.", call. = FALSE)
  }

  # Matrix rows: one for a binary 0/1 (or logical) indicator, one per level for a
  # multi-level factor -- matching the layout plot_upset_strata() builds, so the
  # row count is exactly what will be drawn.
  n_rows <- sum(vapply(strata_vars, function(v) {
    col <- strata_info[[v]]
    if (maihda_dim_is_indicator(col)) 1L else length(maihda_dim_levels(col))
  }, integer(1)))

  n_total <- nrow(strata_info)
  n_cols <- if (is.null(n_strata)) n_total else min(n_total, as.integer(n_strata))

  # ~0.28 in per column and ~0.25 in per row, atop fixed gutters for the y-axis
  # labels (width) and the bar + estimate panels + caption (height); floored so
  # tiny designs still get a usable canvas.
  list(
    width  = max(6, round(2 + 0.28 * n_cols, 1)),
    height = max(4.5, round(4 + 0.25 * n_rows, 1)),
    rows   = n_rows,
    cols   = n_cols
  )
}

#' Effect Decomposition Plot
#'
#' Decomposes the total deviation from the overall mean into the additive (fixed) component
#' and the intersectional (random) component for each stratum.
#'
#' @param object A maihda_model object
#' @param summary_obj A maihda_summary object
#' @param top_n_labels Number of most extreme strata to label
#' @param highlight Optional character vector of stratum ids to highlight. When
#'   supplied, labels are restricted to these strata rather than the most extreme
#'   overall deviations.
#' @return A ggplot2 object
#' @keywords internal
#' @import ggplot2
#' @importFrom dplyr group_by summarise n arrange desc mutate row_number
#' @importFrom utils head
#' @importFrom stats predict setNames fitted
plot_effect_decomposition <- function(object, summary_obj, top_n_labels = 10, highlight = NULL) {
  data <- object$data

  if (!"stratum" %in% names(data)) {
    stop("'stratum' variable not found in data. Make sure to use data from make_strata().")
  }

  # Cross-classified model: the additive part is carried by the dimension random
  # effects (not the fixed effects), so the additive component is computed from the
  # total deviation minus the stratum (interaction) random effect, rather than from
  # the fixed-only prediction.
  cc_mode <- !is.null(object$cc_info)

  # Compute full and fixed-only predictions on the LINK scale. The additive
  # decomposition (total = additive + intersectional) is only exact on the model
  # scale: eta = X*beta + u_stratum. On the response scale, for non-identity links
  # (logit/log) the split is not additive. For Gaussian/identity the link scale
  # equals the response scale, so this is unchanged there.
  if (object$engine == "lme4") {
    preds_total <- tryCatch(predict(object$model, type = "link"), error = function(e) rep(NA, nrow(data)))
    preds_fixed <- tryCatch(predict(object$model, type = "link", re.form = NA), error = function(e) rep(NA, nrow(data)))
  } else if (object$engine == "brms") {
    preds_total <- tryCatch(maihda_brms_linpred_mean(object$model), error = function(e) rep(NA, nrow(data)))
    preds_fixed <- tryCatch(maihda_brms_linpred_mean(object$model, re_formula = NA), error = function(e) rep(NA, nrow(data)))
  } else if (object$engine == "wemix") {
    preds_total <- tryCatch(maihda_wemix_linpred(object, include_re = TRUE), error = function(e) rep(NA, nrow(data)))
    preds_fixed <- tryCatch(maihda_wemix_linpred(object, include_re = FALSE), error = function(e) rep(NA, nrow(data)))
  } else if (object$engine == "ordinal") {
    # The latent location eta = x'beta + u: the additive/intersectional split is
    # exact on this (link) scale, exactly as for the other engines.
    preds_total <- tryCatch(maihda_clmm_linpred(object, include_re = TRUE), error = function(e) rep(NA, nrow(data)))
    preds_fixed <- tryCatch(maihda_clmm_linpred(object, include_re = FALSE), error = function(e) rep(NA, nrow(data)))
  } else {
    stop("Engine not supported for effect decomposition.")
  }

  data$pred_total <- preds_total
  data$pred_fixed <- preds_fixed
  # The model's prediction weights so the per-stratum and global means reflect the
  # weighted fit, trial-weighting each row of an aggregated-binomial fit (the stratum
  # random-effect component below is the weight-invariant BLUP). These are
  # prior/precision (and trial) weights, not a complex survey design, so the result
  # is not survey-representative. Unit weights reproduce the previous unweighted
  # means exactly.
  data$.maihda_w <- maihda_prediction_weights(object)

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

  # additive_dev = total deviation minus the intersectional (stratum) component.
  # Two-model: the additive part is the fixed-effect deviation (mean_fixed - global).
  # Cross-classified: the dimension main effects are random, so the additive part is
  # the total stratum deviation (mean_total - global) net of the interaction RE; this
  # absorbs the dimension REs (plus any covariate deviation), keeping
  # total = additive + interaction in both modes.
  stratum_means <- stratum_means |>
    dplyr::mutate(
      additive_dev = if (cc_mode) {
        .data$mean_total - global_mean - .data$intersectional_dev
      } else {
        .data$mean_fixed - global_mean
      },
      total_dev = .data$additive_dev + .data$intersectional_dev,
      abs_total_dev = abs(.data$total_dev)
    ) |>
    dplyr::arrange(.data$total_dev) |>
    dplyr::mutate(rank = dplyr::row_number())

  if (is.null(highlight) || isFALSE(highlight)) {
    highlight <- NULL
  } else {
    highlight <- as.character(highlight)
  }
  highlight_requested <- !is.null(highlight)

  # Mark strata flagged as carrying a credibly non-zero interaction.
  stratum_means$.maihda_flag <- as.character(stratum_means$stratum) %in% highlight

  additive_label <- if (cc_mode) "Additive (dimension random effects)" else "Fixed-effect component"
  interaction_label <- if (cc_mode) "Intersectional interaction" else "Stratum random-effect component"

  # Create segment definitions for stacking
  # Additive goes from 0 -> additive_dev
  # Intersectional goes from additive_dev -> total_dev
  seg_data <- rbind(
    data.frame(
      rank = stratum_means$rank,
      label = stratum_means$label,
      Component = additive_label,
      y_start = 0,
      y_end = stratum_means$additive_dev,
      abs_total_dev = stratum_means$abs_total_dev,
      flag = stratum_means$.maihda_flag
    ),
    data.frame(
      rank = stratum_means$rank,
      label = stratum_means$label,
      Component = interaction_label,
      y_start = stratum_means$additive_dev,
      y_end = stratum_means$total_dev,
      abs_total_dev = stratum_means$abs_total_dev,
      flag = stratum_means$.maihda_flag
    )
  )

  # Set component ordering so Additive is handled first
  seg_data$Component <- factor(seg_data$Component, levels = c(additive_label, interaction_label))

  has_hl <- any(stratum_means$.maihda_flag)

  # Without an interaction screen, label the most extreme overall deviations.
  # With a screen, label exactly the highlighted strata (e.g. BH survivors),
  # including the zero-row case when no stratum survives, so the labels track the
  # chosen multiplicity rule rather than unadjusted extremes.
  label_data <- if (highlight_requested) {
    stratum_means[stratum_means$.maihda_flag, , drop = FALSE]
  } else {
    stratum_means |>
      dplyr::arrange(dplyr::desc(.data$abs_total_dev)) |>
      utils::head(top_n_labels)
  }
  if (has_hl) {
    label_data$label <- maihda_highlight_label(
      label_data$label, label_data$stratum, highlight)
  }

  seg_colors <- stats::setNames(c("gray60", "#D55E00"), c(additive_label, interaction_label))

  plot_title <- if (cc_mode) {
    "Deviation Decomposition: Additive vs. Interaction (crossed-dimensions)"
  } else {
    "Deviation Decomposition: Fixed vs. Stratum-Random Components"
  }
  # Here colour already encodes the component (additive vs interaction), so the
  # focus-by-contrast highlight rides the opacity channel: flagged strata at full
  # opacity, the rest dimmed -- no ring overlay.
  seg_layer <- if (has_hl) {
    ggplot2::geom_segment(data = seg_data, ggplot2::aes(x = .data$rank, xend = .data$rank, y = .data$y_start, yend = .data$y_end, color = .data$Component, alpha = .data$flag), linewidth = 3)
  } else {
    ggplot2::geom_segment(data = seg_data, ggplot2::aes(x = .data$rank, xend = .data$rank, y = .data$y_start, yend = .data$y_end, color = .data$Component), linewidth = 3, alpha = 0.8)
  }
  total_layer <- if (has_hl) {
    ggplot2::geom_point(data = stratum_means, ggplot2::aes(x = .data$rank, y = .data$total_dev, alpha = .data$.maihda_flag), size = 1.5, color = "black")
  } else {
    ggplot2::geom_point(data = stratum_means, ggplot2::aes(x = .data$rank, y = .data$total_dev), size = 1.5, color = "black")
  }

  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    # Draw segments stacked directly simulating waterfall
    seg_layer +
    # Draw a point at the final Total Deviation
    total_layer +
    # Label extremes
    ggrepel::geom_label_repel(data = label_data, ggplot2::aes(x = .data$rank, y = .data$total_dev, label = .data$label), size = 3, min.segment.length = 0) +
    ggplot2::scale_color_manual(values = seg_colors) +
    ggplot2::labs(
      title = plot_title,
      x = "Stratum Rank (Ordered by Total Predicted Deviation)",
      y = "Deviation from Global Mean (link scale)",
      color = "Effect Component"
    ) +
    theme_maihda() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom"
    )
  if (has_hl) {
    p <- p + ggplot2::scale_alpha_manual(values = maihda_highlight_alpha(), guide = "none")
  }

  return(p)
}

#' Time-varying VPC trajectory plot (longitudinal MAIHDA)
#'
#' The between-stratum share of variance (VPC/ICC) as a function of time, with a
#' confidence/credible ribbon when available. The headline reference-time VPC is
#' marked. For a longitudinal MAIHDA the VPC is not a single number -- the
#' between-stratum variance is a random intercept + slope on time -- so this curve
#' replaces the cross-sectional VPC bar.
#'
#' @param summary_obj A \code{maihda_summary} from a longitudinal model.
#' @return A ggplot2 object.
#' @keywords internal
#' @import ggplot2
plot_vpc_trajectory <- function(summary_obj) {
  lng <- summary_obj$longitudinal
  if (is.null(lng)) {
    stop("plot_vpc_trajectory() needs a longitudinal summary (fit_maihda(id = , ",
         "time = )).", call. = FALSE)
  }
  vt <- lng$vpc_t
  has_ribbon <- any(is.finite(vt$lower) & is.finite(vt$upper))

  p <- ggplot(vt, aes(x = .data$time, y = .data$estimate))
  if (has_ribbon) {
    p <- p + geom_ribbon(aes(ymin = .data$lower, ymax = .data$upper),
                         fill = "#E69F00", alpha = 0.20)
  }
  p <- p +
    geom_line(color = "#E69F00", linewidth = 1.1) +
    geom_point(color = "#E69F00", size = 2) +
    geom_vline(xintercept = lng$ref_time, linetype = "dashed",
               color = "grey50") +
    labs(
      title = "Time-varying VPC/ICC (between-stratum share)",
      subtitle = sprintf("Dashed line: reference time %s = %g (baseline VPC = %.3f)",
                         lng$time, lng$ref_time, summary_obj$vpc$estimate),
      x = lng$time,
      y = "VPC/ICC (between-stratum share of variance)"
    ) +
    theme_maihda() +
    theme(plot.title = element_text(face = "bold"))
  p
}

#' Stratum mean-trajectory plot (longitudinal MAIHDA)
#'
#' One predicted line per stratum over time -- the fixed-part trajectory plus each
#' stratum's random intercept and slope (BLUPs) -- the longitudinal analogue of the
#' predicted-strata caterpillar. Shows how the intersectional groups fan out (or
#' converge) over time.
#'
#' @param object A longitudinal \code{maihda_model}.
#' @param summary_obj Its \code{maihda_summary}.
#' @param n_strata Maximum number of strata to draw; the rest are noted in the
#'   caption.
#' @param select When the cap drops strata, which to keep: \code{"order"}
#'   (default) the first n_strata in stratum order, or \code{"deviation"} the
#'   n_strata whose trajectory swings furthest from the population curve (largest
#'   peak \code{|random deviation|} over the time grid, either direction).
#' @return A ggplot2 object.
#' @keywords internal
#' @import ggplot2
plot_stratum_trajectories <- function(object, summary_obj, n_strata = 50, select = c("order", "deviation")) {
  select <- match.arg(select)
  lng <- summary_obj$longitudinal
  if (is.null(lng)) {
    stop("plot_stratum_trajectories() needs a longitudinal model.", call. = FALSE)
  }
  grid <- lng$time_grid
  # The stratum coefficient vectors are in the model's (possibly centered)
  # coordinates; the displayed grid stays on the original time scale, so the
  # per-stratum polynomial is evaluated at grid - center.
  grid_c <- grid - maihda_lng_time_center(object$longitudinal_info)
  # Per-stratum random intercept + slope (BLUPs) on the time polynomial.
  re <- maihda_longitudinal_stratum_re(object)
  strata <- re$stratum
  omitted <- 0L
  if (!is.null(n_strata) && length(strata) > n_strata) {
    omitted <- length(strata) - n_strata
    keep_rows <- if (identical(select, "deviation")) {
      # Peak absolute random deviation over the grid -- max_t |sum_j coef_j t^j|,
      # the stratum's own departure from the fixed trajectory -- so the strata whose
      # trajectories diverge most (either direction) survive the cap, not the first
      # by stratum id. A scalar BLUP would miss a small-intercept/large-slope fan-out.
      mag <- vapply(seq_len(nrow(re)), function(i) {
        a <- vapply(grid_c, function(t) sum(re$coef[[i]] * t^(0:(length(re$coef[[i]]) - 1))),
                    numeric(1))
        max(abs(a))
      }, numeric(1))
      order(mag, decreasing = TRUE)[seq_len(n_strata)]
    } else {
      seq_len(n_strata)
    }
    re <- re[sort(keep_rows), , drop = FALSE]
    strata <- re$stratum
  }

  # Fixed-part trajectory at the population mean covariate profile: predict on a
  # one-row-per-grid-time frame holding covariates at the data means, RE excluded.
  eta_fixed <- maihda_longitudinal_fixed_trajectory(object, grid)

  rows <- do.call(rbind, lapply(seq_len(nrow(re)), function(i) {
    a <- vapply(grid_c, function(t) sum(re$coef[[i]] * t^(0:(length(re$coef[[i]]) - 1))),
                numeric(1))
    data.frame(stratum = re$stratum[i],
               label = if (!is.null(re$label)) re$label[i] else re$stratum[i],
               time = grid, value = eta_fixed + a, stringsAsFactors = FALSE)
  }))

  cap <- if (omitted > 0) {
    sprintf("%d of %d strata shown%s; %d omitted.",
            length(strata), length(strata) + omitted,
            if (identical(select, "deviation")) " (most extreme by trajectory deviation)" else "",
            omitted)
  } else NULL

  ggplot(rows, aes(x = .data$time, y = .data$value, group = .data$stratum,
                   color = .data$label)) +
    geom_line(alpha = 0.8, linewidth = 0.7) +
    labs(
      title = "Predicted stratum trajectories",
      subtitle = "Fixed-part trajectory + each stratum's random intercept & slope",
      x = lng$time, y = "Predicted outcome (link scale)", color = "Stratum",
      caption = cap
    ) +
    theme_maihda() +
    theme(plot.title = element_text(face = "bold"),
          legend.position = if (length(strata) > 12) "none" else "right")
}

#' Time-specific PCV plot (longitudinal MAIHDA)
#'
#' The additive share -- the proportional change in the between-stratum (trajectory)
#' variance from the null to the adjusted model -- as a function of time. A high,
#' flat curve means intersectional trajectory inequalities are "mostly additive".
#'
#' @param pcv A \code{maihda_long_pcv} from a longitudinal \code{maihda()} pair.
#' @return A ggplot2 object.
#' @keywords internal
#' @import ggplot2
plot_pcv_trajectory <- function(pcv) {
  if (!inherits(pcv, "maihda_long_pcv")) {
    stop("plot_pcv_trajectory() needs a maihda_long_pcv object.", call. = FALSE)
  }
  d <- pcv$pcv_t
  ggplot(d, aes(x = .data$time, y = .data$pcv)) +
    geom_line(color = "#0072B2", linewidth = 1.1) +
    geom_point(color = "#0072B2", size = 2) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
    labs(
      title = "Additive share of between-stratum (trajectory) variance over time",
      subtitle = sprintf("PCV(t) = (Var_null(t) - Var_adjusted(t)) / Var_null(t); time = %s",
                         pcv$time),
      x = pcv$time, y = "PCV(t) (additive share)"
    ) +
    theme_maihda() +
    theme(plot.title = element_text(face = "bold"))
}

# Per-stratum random-effect coefficient vector (intercept, slope, ...) for the
# stratum grouping of a longitudinal fit, as a data frame with a list-column
# `coef`. Used by plot_stratum_trajectories(). Engine-aware (lme4 ranef /
# brms ranef posterior means).
maihda_longitudinal_stratum_re <- function(object) {
  lng <- object$longitudinal_info
  if (identical(object$engine, "lme4")) {
    re <- lme4::ranef(object$model)[["stratum"]]
    coefs <- lapply(seq_len(nrow(re)), function(i) as.numeric(re[i, ]))
    ids <- rownames(re)
  } else if (identical(object$engine, "brms")) {
    arr <- brms::ranef(object$model)[["stratum"]]
    ids <- dimnames(arr)[[1]]
    coefs <- lapply(seq_along(ids), function(i) as.numeric(arr[i, "Estimate", ]))
  } else {
    stop("Longitudinal trajectories are available for lme4/brms only.", call. = FALSE)
  }
  # Labels via the shared helper (row order preserved), then attach the per-stratum
  # coefficient vectors (intercept, slope, ...) as a list-column aligned by row.
  out <- add_stratum_labels(
    data.frame(stratum = ids, stratum_id = suppressWarnings(as.integer(ids)),
               random_effect = vapply(coefs, `[`, numeric(1), 1),
               stringsAsFactors = FALSE),
    object$strata_info)
  out$coef <- coefs
  out
}

# Fixed-part trajectory (NO random effects, re.form = NA) at the mean covariate
# profile, over a time grid. Builds a prediction frame holding every non-time
# covariate at its mean (numeric) or modal (factor) value and varying only time.
maihda_longitudinal_fixed_trajectory <- function(object, grid) {
  lng <- object$longitudinal_info
  time_term <- maihda_lng_time_term(lng)
  data <- object$data
  # RHS vars only: an addition-term response (y | trials(n)) would leave its
  # trials/weights variable in an all.vars(formula)[-1] extraction.
  fixed_vars <- all.vars(maihda_nobars(object$formula)[[3]])
  nd <- data[rep(1L, length(grid)), , drop = FALSE]
  for (v in intersect(fixed_vars, names(nd))) {
    if (v %in% c(lng$time, time_term)) next
    col <- data[[v]]
    nd[[v]] <- if (is.numeric(col)) mean(col, na.rm = TRUE) else {
      tb <- sort(table(col), decreasing = TRUE)
      rep(names(tb)[1], length(grid))
    }
  }
  nd[[lng$time]] <- grid
  # A centered fit's formula references the derived centered column, not the
  # original time; keep it aligned with the (original-scale) grid.
  if (!identical(time_term, lng$time)) {
    nd[[time_term]] <- grid - maihda_lng_time_center(lng)
  }

  if (identical(object$engine, "lme4")) {
    # predict.merMod(newdata = nd) drops an external offset= and errors on a formula
    # offset() term, so build the fixed trajectory directly. A FORMULA offset such as
    # offset(0.5 * time) is RE-EVALUATED on the trajectory grid (nd carries the grid
    # times, raw and centered), so a time-dependent offset tracks the trajectory instead
    # of being frozen; an EXTERNAL offset= vector has no expression to re-evaluate and is
    # held at its mean -- the representative-exposure analogue of holding covariates at
    # their means. Either part is NULL when absent; a no-offset fit reproduces
    # predict(re.form = NA) exactly.
    off_ext <- maihda_fitted_offset_external(object$model)
    off <- if (is.null(off_ext)) NULL else mean(off_ext, na.rm = TRUE)
    if (!is.null(attr(stats::terms(maihda_nobars(stats::formula(object$model))),
                      "offset"))) {
      # fallback = "mean": an offset term whose raw variable is absent from nd (stored
      # only as the derived "offset(...)" column, necessarily time-invariant) is held at
      # its mean, like the external offset above.
      fo <- maihda_lme4_formula_offset_at(object$model, nd, fallback = "mean")
      off <- if (is.null(off)) fo else off + fo
    }
    as.numeric(maihda_lme4_fixed_link(object$model, nd, offset = off))
  } else if (identical(object$engine, "brms")) {
    maihda_brms_linpred_mean(object$model, newdata = nd, re_formula = NA)
  } else {
    stop("Longitudinal trajectories are available for lme4/brms only.", call. = FALSE)
  }
}
