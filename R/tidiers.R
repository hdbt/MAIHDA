# Tidy / glance methods (broom interface) for the MAIHDA classes.
#
# These import the generics from the lightweight `generics` package -- the same
# generics broom/broom.mixed re-export -- so `tidy()` / `glance()` dispatch to the
# methods below whether the user loads broom, generics, or just MAIHDA. We do NOT
# depend on broom itself.
#
# Design: the methods only repackage quantities the summary already computes; they
# add no new statistics and fit no models. The unique, non-redundant content is the
# one-row `glance()` of the MAIHDA headline (VPC + PCV + AUC/MOR + additive/
# interaction shares), uniform across all engines -- nothing in broom.mixed/easystats
# produces it from the underlying fit. The raw fixed-effect / per-row tidying that
# broom.mixed already does on the underlying model is intentionally NOT reimplemented.

#' @importFrom generics tidy
#' @export
generics::tidy

#' @importFrom generics glance
#' @export
generics::glance

# NA-safe scalar extractor: NULL/empty -> NA_real_, else the first value as numeric.
maihda_num <- function(v) {
  if (is.null(v) || length(v) == 0L) return(NA_real_)
  as.numeric(v)[1]
}

# NA-safe scalar character extractor.
maihda_chr <- function(v) {
  if (is.null(v) || length(v) == 0L) return(NA_character_)
  as.character(v)[1]
}

# Normalise the engine-specific fixed_effects slot to a broom-shaped tibble:
#   term / estimate / std.error / statistic / p.value / conf.low / conf.high
# The lme4, wemix and ordinal summaries store the full Wald table (term / estimate /
# se / statistic / p_value / lower / upper) built by maihda_fixed_effects_table();
# brms::fixef() gives a matrix with Estimate / Est.Error and two quantile columns
# (named for the quantiles themselves, e.g. Q2.5 / Q97.5), and a posterior summary
# carries no test statistic or p-value, so those are NA there.
#
# The estimate/se fallbacks below cover a summary built before those columns existed
# (a stored maihda_summary object from an older version).
maihda_tidy_fixed <- function(fe) {
  if (is.null(fe)) return(tibble::tibble())

  if (is.matrix(fe)) {
    cn <- colnames(fe)
    # Interval columns are the two extreme quantile columns, whatever level the
    # fit was summarised at (Q2.5/Q97.5 at the 95% default).
    qcol <- grep("^Q[0-9.]+$", cn, value = TRUE)
    qval <- suppressWarnings(as.numeric(sub("^Q", "", qcol)))
    has_q <- length(qcol) >= 2L && !anyNA(qval)
    return(tibble::tibble(
      term      = rownames(fe),
      estimate  = if ("Estimate" %in% cn) fe[, "Estimate"] else fe[, 1],
      std.error = if ("Est.Error" %in% cn) fe[, "Est.Error"] else NA_real_,
      statistic = NA_real_,
      df        = NA_real_,
      p.value   = NA_real_,
      conf.low  = if (has_q) fe[, qcol[which.min(qval)]] else NA_real_,
      conf.high = if (has_q) fe[, qcol[which.max(qval)]] else NA_real_
    ))
  }

  est <- as.numeric(fe$estimate)
  se  <- if ("se" %in% names(fe)) as.numeric(fe$se) else NA_real_
  se[!is.na(se) & se <= 0] <- NA_real_
  z <- stats::qnorm(0.975)

  tibble::tibble(
    term      = as.character(fe$term),
    estimate  = est,
    std.error = se,
    statistic = if ("statistic" %in% names(fe)) as.numeric(fe$statistic) else est / se,
    # NA on every engine but the Gaussian lme4 one, which reports containment
    # (between-within) denominator degrees of freedom; broom's lmerTest tidier
    # carries the same column.
    df        = if ("df" %in% names(fe)) as.numeric(fe$df) else NA_real_,
    p.value   = if ("p_value" %in% names(fe)) {
      as.numeric(fe$p_value)
    } else {
      2 * stats::pnorm(-abs(est / se))
    },
    conf.low  = if ("lower" %in% names(fe)) as.numeric(fe$lower) else est - z * se,
    conf.high = if ("upper" %in% names(fe)) as.numeric(fe$upper) else est + z * se
  )
}

#' Tidy a MAIHDA summary, model, or analysis
#'
#' \code{\link[generics]{tidy}} methods that return the MAIHDA estimates as a tidy
#' \code{tibble}, ready for downstream tables (\code{gt}, \code{flextable}) and
#' \code{ggplot2}. They read the slots that \code{\link{summary.maihda_model}}
#' already computes and add no new statistics.
#'
#' @param x A \code{maihda_summary} (from \code{\link[=summary.maihda_model]{summary}}),
#'   a \code{maihda_model} (from \code{\link{fit_maihda}}), or a
#'   \code{maihda_analysis} (from \code{\link{maihda}}).
#' @param component Which estimates to return:
#'   \describe{
#'     \item{\code{"strata"}}{(default) the stratum (intersection) random-effect
#'       estimates -- one row per stratum, with \code{estimate}, \code{std.error}
#'       and \code{conf.low}/\code{conf.high}, plus the human-readable intersectional
#'       \code{label} when available.}
#'     \item{\code{"variance"}}{the variance-components table (between-stratum, any
#'       other random effects, residual, and total) with each component's variance,
#'       SD and proportion.}
#'     \item{\code{"fixed"}}{the fixed-effect estimates in broom's
#'       \code{term}/\code{estimate}/\code{std.error}/\code{statistic}/
#'       \code{p.value}/\code{conf.low}/\code{conf.high} shape.}
#'   }
#' @param which For a \code{maihda_analysis}, whether to tidy the \code{"null"}
#'   (default) or \code{"adjusted"} model's summary.
#' @param ... Unused, for S3 consistency.
#'
#' @return A \code{tibble}. For \code{component = "strata"}: columns
#'   \code{stratum}, \code{label}, \code{estimate}, \code{std.error},
#'   \code{conf.low}, \code{conf.high}. For \code{"variance"}: \code{component},
#'   \code{variance}, \code{sd}, \code{proportion}. For \code{"fixed"}: \code{term},
#'   \code{estimate}, \code{std.error}, \code{statistic}, \code{df}, \code{p.value},
#'   \code{conf.low}, \code{conf.high}.
#'
#' @section Fixed-effect statistics: For the lme4, WeMix and ordinal engines
#'   \code{statistic} is the Wald statistic (\code{estimate / std.error}),
#'   \code{p.value} its two-sided p-value, \code{conf.low}/\code{conf.high} the
#'   matching Wald interval at the summary's \code{conf_level}, and \code{df}
#'   the denominator degrees of freedom they were built on.
#'
#'   A Gaussian \pkg{lme4} fit reports containment (between-within) \code{df}
#'   and a \eqn{t}; every other engine leaves \code{df} \code{NA} and uses a z.
#'   Under \code{summary(df_method = "bootstrap")} \code{p.value} and the
#'   interval come from the parametric bootstrap and \code{df} is \code{NA}.
#'   See \code{\link{summary.maihda_model}}. For brms the estimate is the
#'   posterior mean with its \code{Est.Error} and credible interval, and
#'   \code{statistic}/\code{df}/\code{p.value} are \code{NA}. Standard errors
#'   are \code{NA} where they are undefined (a boundary fit).
#'
#' @seealso \code{\link{glance.maihda_analysis}} for the one-row model summary.
#'
#' @examples
#' data("maihda_health_data")
#' m <- fit_maihda(BMI ~ Age + (1 | Gender:Race:Education), data = maihda_health_data)
#' tidy(m)                       # stratum estimates
#' tidy(m, component = "variance")
#' tidy(m, component = "fixed")
#'
#' @name maihda_tidiers
#' @method tidy maihda_summary
#' @export
tidy.maihda_summary <- function(x, component = c("strata", "variance", "fixed"), ...) {
  component <- match.arg(component)

  if (component == "variance") {
    if (is.null(x$variance_components)) return(tibble::tibble())
    return(tibble::as_tibble(x$variance_components))
  }

  if (component == "fixed") {
    return(maihda_tidy_fixed(x$fixed_effects))
  }

  # component == "strata"
  se <- x$stratum_estimates
  if (is.null(se) || nrow(se) == 0L) return(tibble::tibble())
  tibble::tibble(
    stratum   = as.character(se$stratum),
    label     = if ("label" %in% names(se)) as.character(se$label) else NA_character_,
    estimate  = as.numeric(se$random_effect),
    std.error = if ("se" %in% names(se)) as.numeric(se$se) else NA_real_,
    conf.low  = if ("lower_95" %in% names(se)) as.numeric(se$lower_95) else NA_real_,
    conf.high = if ("upper_95" %in% names(se)) as.numeric(se$upper_95) else NA_real_
  )
}

#' @rdname maihda_tidiers
#' @method tidy maihda_model
#' @export
tidy.maihda_model <- function(x, component = c("strata", "variance", "fixed"), ...) {
  tidy(summary(x), component = component, ...)
}

#' @rdname maihda_tidiers
#' @method tidy maihda_analysis
#' @export
tidy.maihda_analysis <- function(x, component = c("strata", "variance", "fixed"),
                                 which = c("null", "adjusted"), ...) {
  component <- match.arg(component)
  which <- match.arg(which)
  s <- if (which == "adjusted") x$summary_adjusted else x$summary
  if (is.null(s)) {
    stop("No '", which, "' summary is available on this maihda_analysis (mode = '",
         maihda_chr(x$mode), "').", call. = FALSE)
  }
  tidy(s, component = component, ...)
}

#' Glance at a MAIHDA model or analysis
#'
#' \code{\link[generics]{glance}} methods that return the MAIHDA headline as a
#' one-row \code{tibble}: the variance partition coefficient (VPC/ICC), and -- for a
#' \code{maihda_analysis} -- the proportional change in variance (PCV), plus the
#' additive/interaction shares and the discriminatory accuracy (AUC, MOR) for a
#' binomial outcome. The layout is uniform across the lme4, brms, WeMix and ordinal
#' engines. No other package emits this row from the underlying fit, because PCV
#' needs the null+adjusted pair that only a \code{maihda_analysis} holds.
#'
#' @param x A \code{maihda_summary}, \code{maihda_model}, or \code{maihda_analysis}.
#' @param ... Unused, for S3 consistency.
#'
#' @return A one-row \code{tibble}. \code{glance.maihda_analysis} adds \code{pcv}
#'   (with \code{pcv.conf.low}/\code{pcv.conf.high} when bootstrapped or from a brms
#'   posterior), the adjusted-model \code{auc.adjusted}, \code{nobs}, \code{family}
#'   and \code{mode} to the columns produced for a single summary.
#'
#' @seealso \code{\link{maihda_tidiers}} for the per-estimate \code{tidy()} methods.
#'
#' @examples
#' data("maihda_health_data")
#' a <- maihda(BMI ~ Age + Gender + Race + Education + (1 | Gender:Race:Education),
#'             data = maihda_health_data)
#' glance(a)
#'
#' @name maihda_glance
#' @method glance maihda_summary
#' @export
glance.maihda_summary <- function(x, ...) {
  # Locals must not share a name with any output column: tibble() exposes
  # already-built columns to later expressions, so a local `vpc` would be shadowed
  # by the numeric `vpc` column (turning `vpc$ci_lower` into a `$`-on-atomic error).
  vp  <- x$vpc
  dec <- x$decomposition
  da  <- x$discriminatory_accuracy
  tibble::tibble(
    vpc               = maihda_num(vp$estimate),
    vpc.conf.low      = maihda_num(vp$ci_lower),
    vpc.conf.high     = maihda_num(vp$ci_upper),
    additive.share    = maihda_num(dec$additive_share),
    interaction.share = maihda_num(dec$interaction_share),
    auc               = maihda_num(da$auc),
    mor               = maihda_num(da$mor),
    n_strata          = if (is.null(x$stratum_estimates)) NA_integer_ else nrow(x$stratum_estimates),
    engine            = maihda_chr(x$engine)
  )
}

#' @rdname maihda_glance
#' @method glance maihda_model
#' @export
glance.maihda_model <- function(x, ...) {
  g <- glance(summary(x))
  g$nobs <- tryCatch(nrow(x$data), error = function(e) NA_integer_)
  g$family <- tryCatch(maihda_chr(x$family$family), error = function(e) NA_character_)
  g[, c("vpc", "vpc.conf.low", "vpc.conf.high", "additive.share", "interaction.share",
        "auc", "mor", "n_strata", "nobs", "engine", "family")]
}

#' @rdname maihda_glance
#' @method glance maihda_analysis
#' @export
glance.maihda_analysis <- function(x, ...) {
  g <- glance(x$summary)

  pcv <- x$pcv
  g$pcv          <- maihda_num(pcv$pcv)
  g$pcv.conf.low <- maihda_num(pcv$ci_lower)
  g$pcv.conf.high <- maihda_num(pcv$ci_upper)

  da_adj <- x$summary_adjusted$discriminatory_accuracy
  g$auc.adjusted <- maihda_num(da_adj$auc)

  g$nobs   <- tryCatch(nrow(x$model$data), error = function(e) NA_integer_)
  g$family <- tryCatch(maihda_chr(x$model$family$family), error = function(e) NA_character_)
  g$mode   <- maihda_chr(x$mode)

  g[, c("vpc", "vpc.conf.low", "vpc.conf.high",
        "pcv", "pcv.conf.low", "pcv.conf.high",
        "additive.share", "interaction.share",
        "auc", "auc.adjusted", "mor",
        "n_strata", "nobs", "engine", "family", "mode")]
}
