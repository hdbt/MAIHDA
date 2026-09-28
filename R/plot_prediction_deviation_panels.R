maihda_binomial_observed_01 <- function(x, n) {
  if (is.null(x) || length(x) != n || !maihda_is_binary_vector(x)) {
    return(rep(NA_integer_, n))
  }

  maihda_binary_to_01(x)
}

# The 0/1 coding of the plotted rows' outcome `x` AS THE MODEL CODED IT, on the fitted
# rows as on supplied ones: 1 for the event, 0 for the reference, NA for a missing
# value or one the fit never saw. Keyed
# on the fit, never on `x` itself, whose factor levels may be declared in another
# order and whose rows may hold only one of the two values: maihda_binomial_observed_01()
# re-derived the coding from the plotted rows, so the same rows with their factor
# levels reversed had every residual inverted and every "Wrong" / "Correct" swapped,
# and rows holding a single value could not be coded at all. Sources, in order:
#   * a fit_maihda() model's $response_recoding, the original label -> 0/1 record of
#     the analytic rows. That record is by label, so a column of `data` is read in the
#     original coding -- except a numeric 0/1 column that can only be the recoded one
#     (the stored `object$data` handed back), taken as it stands: one whose values are
#     not all original labels, or any, when the original cannot have been numeric
#     (its recorded levels are not ascending numbers). A numeric column holding
#     nothing but a value that is both (only 1s, on a 1/2-coded outcome) is read in
#     the original coding, as documented for `data`.
#   * otherwise the fitted response -- the model frame's, or for a wemix model, whose
#     stored rows are plain data, its outcome expression evaluated on them: a
#     factor's first level is the failure, as in R's binomial family and in brms
#     (whose binary families code by position in the fit's response levels); a
#     logical's TRUE is the event; 0/1 numbers are themselves.
#   * with neither (a bare WeMix fit exposes no response), the plotted rows' own two
#     values, as before.
# A value the fit never saw (another label, a stray code) is NA like a missing one; the
# rows it affects are returned as the attribute "unmatched" -- when the coding came
# from the fit -- so the panel can say so, as the ordinal surprise panel does for a
# category the fit does not have.
maihda_prediction_panel_fit_01 <- function(maihda_obj, model, x) {
  n <- length(x)
  known <- !is.na(x)
  xc <- as.character(x)
  flag <- function(out) {
    out[!known] <- NA_integer_
    structure(out, unmatched = which(known & is.na(out)))
  }
  rr <- if (!is.null(maihda_obj)) maihda_obj$response_recoding
  if (is.data.frame(rr) && all(c("level", "value") %in% names(rr)) && nrow(rr) == 2L) {
    lev <- as.character(rr$level)
    # TRUE / FALSE are 1 / 0 to a 0/1-coded outcome: read a logical column that way when
    # its own labels are not the recorded ones.
    if (is.logical(x) && !all(xc[known] %in% lev)) {
      xc <- as.character(as.integer(x))
    }
    # A numeric original is recorded in ascending order (maihda_binary_levels() sorts
    # it), so levels that are not ascending numbers -- words, FALSE/TRUE, or "1"/"0"
    # declared in that order as a factor -- came from a column a numeric one cannot
    # be in the original coding of: a numeric column is then the recoded one. Reading
    # it by label instead inverted object$data for a factor(y, levels = c(1, 0)).
    lev_num <- suppressWarnings(as.numeric(lev))
    numeric_original <- !anyNA(lev_num) && lev_num[1] < lev_num[2]
    recoded <- is.numeric(x) && all(x[known] %in% c(0, 1)) &&
      (!numeric_original || !all(xc[known] %in% lev))
    out <- if (recoded) {
      as.integer(x)
    } else {
      as.integer(rr$value[match(xc, lev)])
    }
    return(flag(out))
  }

  frame <- maihda_prediction_panel_fit_rows(maihda_obj, model)
  y <- if (is.data.frame(frame)) {
    tryCatch(stats::model.response(frame), error = function(e) NULL)
  }
  # A wemix model stores its analytic rows as plain data, not a model frame; its
  # fitted response is its outcome expression evaluated on those rows (an expression
  # outcome carries no $response_recoding, only bare columns being recoded).
  if (is.null(y) && !is.null(maihda_obj) && is.data.frame(frame)) {
    y <- tryCatch(maihda_response_from_frame(frame, maihda_obj$formula),
                  error = function(e) NULL)
  }
  if (!is.null(y) && is.null(dim(y))) {
    out <- NULL
    if (is.factor(y)) {
      lev <- levels(y)
      out <- ifelse(xc == lev[1], 0L, ifelse(xc %in% lev[-1], 1L, NA_integer_))
    } else if (is.logical(y)) {
      xl <- if (is.logical(x)) {
        x
      } else if (is.numeric(x)) {
        ifelse(x %in% c(0, 1), x == 1, NA)
      } else {
        ifelse(xc %in% c("TRUE", "FALSE"), xc == "TRUE", NA)
      }
      out <- as.integer(xl)
    } else if (is.numeric(y) && all(y[!is.na(y)] %in% c(0, 1))) {
      xn <- if (is.logical(x)) as.numeric(x) else suppressWarnings(as.numeric(xc))
      out <- ifelse(xn %in% c(0, 1), as.integer(xn), NA_integer_)
    }
    if (!is.null(out)) {
      return(flag(as.integer(out)))
    }
  }

  maihda_binomial_observed_01(x, n)
}

# |Bernoulli deviance residual| of a 0/1 outcome at probability `fitted`. NA where a
# row cannot be scored -- its outcome is unknown or its probability is not finite --
# so a stratum mean leaves it out. It was 0, which counted such a row as a perfect
# fit: plotted on the original data of an unweighted fit, each stratum's mean residual
# shrank by exactly its share of missing outcomes.
maihda_binomial_abs_deviance_residual <- function(obs_outcome_01, fitted) {
  out <- rep(NA_real_, length(fitted))
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

# |Binomial deviance residual| of `successes` out of `trials` at probability `p`:
# sqrt(2 * [s log(s / (n p)) + (n - s) log((n - s) / (n (1 - p)))]), a 0 * log(0)
# term counting as 0. This is lme4's residuals(type = "deviance") for a
# cbind(successes, failures) fit with the sign dropped, and at one trial it is the
# 0/1 residual above. NA where a row cannot be scored (missing counts, no trials,
# successes outside 0..trials, a non-finite probability), so stratum means leave it
# out.
maihda_binomial_abs_deviance_residual_agg <- function(successes, trials, p) {
  out <- rep(NA_real_, length(p))
  ok <- is.finite(successes) & is.finite(trials) & is.finite(p) & trials > 0 &
    successes >= 0 & successes <= trials
  if (!any(ok)) {
    return(out)
  }
  s <- successes[ok]
  n <- trials[ok]
  f <- n - s
  pp <- pmin(pmax(p[ok], .Machine$double.eps), 1 - .Machine$double.eps)
  term_s <- ifelse(s > 0, s * log(s / (n * pp)), 0)
  term_f <- ifelse(f > 0, f * log(f / (n * (1 - pp))), 0)
  out[ok] <- sqrt(pmax(2 * (term_s + term_f), 0))
  out
}

# Prediction weights aligned to `data`'s rows, used to make the per-stratum
# aggregation a weighted mean for weighted fits (consistent with the weighted VPC
# and the other stratum-level plots); for an aggregated-binomial fit each row is
# weighted by its binomial TRIAL count, matching the trial-weighting the raw-model
# fallback below already applies via weights(type = "prior"). These are
# prior/precision (and trial, or sampling) weights, not a complex survey design (no
# design-based variance is computed).
#
# `fitted_data` says the plotted rows ARE the rows the model was fitted to -- the
# panel was called without `data` -- so the fit's stored weights align with them by
# construction and are used as they always were -- unit weights when the model is
# unweighted (the weighted means then reduce EXACTLY to plain means), a bare brms
# `y | trials(n)` fit (no weights() method) weighted by its trial counts -- except that
# a fit whose weights() has no method but whose frame carries "(weights)" (a bare
# MASS::polr, ordinal::clm or ordinal::clmm) is now weighted by them, as its supplied
# rows are.
#
# Any other `data` goes to maihda_prediction_panel_row_weights(), which finds each
# row's own weight. The fit's stored vector used to be paired with those rows by
# position whenever the LENGTHS matched, which says nothing about which row is which:
# on the fixture of test-audit-2026-09-27b.R, reordering the fitted rows moved a
# precision-weighted Gaussian stratum mean by 2.57 and a cbind() stratum probability by
# 0.082; new rows of the same count took the training rows' weights instead of their
# own, and a subset was silently unweighted.
# `values` are the plotted values of `data`'s rows and `type` the panel type, used
# to recognise fitted rows.
maihda_prediction_panel_prior_weights <- function(maihda_obj, model, data,
                                                  fitted_data = FALSE,
                                                  values = NULL, type = NULL) {
  n <- nrow(data)
  if (!isTRUE(fitted_data)) {
    return(maihda_prediction_panel_row_weights(maihda_obj, model, data,
                                               values = values, type = type))
  }
  w <- NULL
  if (!is.null(maihda_obj)) {
    w <- tryCatch(maihda_prediction_weights(maihda_obj), error = function(e) NULL)
  }
  if (is.null(w)) {
    w <- tryCatch(maihda_fit_rows_weights(model), error = function(e) NULL)
  }
  if (is.null(w)) {
    # weights() has no method for a MASS::polr or a bare ordinal::clm / clmm fit, but
    # its model frame -- which is `data` here -- carries the weights it was given.
    # Without this a weighted polr was unweighted on its own rows while its supplied
    # rows, read by maihda_prediction_panel_row_prior(), were weighted: one fit, two
    # answers.
    w <- maihda_prediction_panel_fit_prior(data)
  }
  if (is.null(w) || !is.numeric(w) || length(w) != n) {
    w <- tryCatch(maihda_prediction_panel_brms_trials(model, data),
                  error = function(e) NULL)
    if (is.null(w)) {
      return(rep(1, n))
    }
  }
  w <- as.numeric(w)
  w[!is.finite(w)] <- NA_real_
  w
}

# The weight of each row of `data` for the stratum summaries when `data` is not the
# fitted rows: the fit's own weighting evaluated on THAT row --
#   * a design-weighted maihda_model: the row's sampling weight, from its column;
#   * a glm / lm / lme4 fit with prior weights: the row's prior weight
#     (maihda_prediction_panel_row_prior());
#   * times the row's trial count for an aggregated binomial (cbind(), brms trials()).
# A fit without weights gives every row 1. When the weight of any row with a finite
# plotted value cannot be found, EVERY row is weighted equally and a warning says so:
# the stratum summaries are then the plain means over the plotted rows, which is a
# defined quantity, where a guessed weight would not be. The warning's
# "residuals_warned" attribute tells the binomial panel it has also explained the rows
# whose deviance residuals that leaves undefined.
maihda_prediction_panel_row_weights <- function(maihda_obj, model, data, values = NULL,
                                                type = NULL) {
  n <- nrow(data)
  frame <- maihda_prediction_panel_fit_rows(maihda_obj, model)
  kind <- character(0)
  base <- NULL
  sw <- if (!is.null(maihda_obj)) maihda_obj$sampling_weights
  if (!is.null(sw)) {
    # A brms fit's stored frame carries the weights normalised to mean 1 under their
    # own name; the scale does not change a weighted mean.
    col <- if (sw %in% names(data)) {
      sw
    } else if (.maihda_brms_weights_col %in% names(data)) {
      .maihda_brms_weights_col
    }
    base <- if (!is.null(col)) maihda_prediction_panel_clean_weights(data[[col]], n)
    if (is.null(base)) base <- rep(NA_real_, n)
    kind <- "sampling weights"
  } else {
    base <- maihda_prediction_panel_row_prior(model, data, frame, values, type)
    if (!is.null(base)) kind <- "prior weights"
  }
  trials <- maihda_prediction_panel_row_trials(model, data, frame)
  if (!is.null(trials)) kind <- c(kind, "binomial trial counts")

  w <- rep(1, n)
  if (!is.null(base)) w <- w * base
  if (!is.null(trials)) w <- w * trials
  plotted <- if (is.numeric(values) && length(values) == n) {
    is.finite(values)
  } else {
    rep(TRUE, n)
  }
  unknown <- plotted & !is.finite(w)
  if (!any(unknown)) {
    return(w)
  }
  also <- ""
  if (identical(type, "binomial") && "prior weights" %in% kind) {
    also <- paste0("; they are also this binomial fit's trial counts, so those rows' ",
                   "deviance residuals are undefined and they are left out of the ",
                   "residual summaries and labels")
  }
  found_in <- c(
    `sampling weights` = sprintf("a sampling weight from the '%s' column of `data`",
                                 if (is.null(sw)) "" else sw),
    `prior weights` = paste0("a prior weight from the weights column of `data` or from ",
                             "the fitted row with the same row name and prediction"),
    `binomial trial counts` = paste0("a trial count from its observed successes and ",
                                     "failures, or its trials() term, in `data`")
  )[kind]
  warning(sprintf(paste0(
    "plot_prediction_deviation_panels(): the model's %s are not known for %d of the %d ",
    "rows of `data`, so every row is weighted equally in the stratum summaries%s. A row ",
    "takes %s; omit `data` to summarise the fitted rows with their weights."),
    paste(kind, collapse = " and "), sum(unknown), n, also,
    paste(found_in, collapse = ", and ")), call. = FALSE)
  structure(rep(1, n), residuals_warned = nzchar(also))
}

# The prior (precision) weight of each row of `data`, for a glm / lm / lme4 fit given
# prior weights; NULL for a fit without them. Read from a "(weights)" column of `data`
# (the stored frame handed back), else from the fit's `weights =` expression evaluated
# on the rows, else from the fitted row each row is shown to be
# (maihda_prediction_panel_row_identity()). NA where none of these can say.
maihda_prediction_panel_row_prior <- function(model, data, frame, values = NULL,
                                              type = NULL) {
  fit_w <- maihda_prediction_panel_fit_prior(frame)
  if (is.null(fit_w)) {
    return(NULL)
  }
  n <- nrow(data)
  w <- if ("(weights)" %in% names(data)) {
    maihda_prediction_panel_clean_weights(data[["(weights)"]], n)
  } else {
    maihda_prediction_panel_call_weights(model, data)
  }
  if (!is.null(w)) {
    return(w)
  }
  idx <- maihda_prediction_panel_row_identity(
    frame, data, maihda_prediction_panel_fit_row_values(model, type), values)
  if (is.null(idx)) rep(NA_real_, n) else fit_w[idx]
}

# The rows the model was fitted to, as stored: a maihda_model's analytic frame (lme4,
# brms) or analytic data (wemix, ordinal), otherwise the model's own frame; NULL when
# there is neither.
maihda_prediction_panel_fit_rows <- function(maihda_obj, model) {
  if (!is.null(maihda_obj) && is.data.frame(maihda_obj$data)) {
    return(maihda_obj$data)
  }
  fr <- maihda_model_frame(model)
  if (is.data.frame(fr)) fr else NULL
}

# The prior weights a fit was given, one per fitted row -- the "(weights)" column of a
# glm / lm / lme4 (or polr / clm) model frame -- or NULL when they are all 1 or absent.
# For a cbind() response these are the weights BEFORE R's binomial family multiplies
# them by each row's trials.
maihda_prediction_panel_fit_prior <- function(frame) {
  w <- if (is.data.frame(frame)) frame[["(weights)"]]
  if (!is.numeric(w) || !is.null(dim(w))) {
    return(NULL)
  }
  w <- as.numeric(w)
  if (isTRUE(all(abs(w - 1) < sqrt(.Machine$double.eps)))) {
    return(NULL)
  }
  w
}

# Weights read for `n` rows; NULL unless one number per row. A missing, non-finite or
# negative weight is unknown (NA); a zero weight is kept, and counts for nothing.
maihda_prediction_panel_clean_weights <- function(w, n) {
  if (!is.numeric(w) || !is.null(dim(w)) || length(w) != n) {
    return(NULL)
  }
  w <- as.numeric(w)
  w[!is.finite(w) | w < 0] <- NA_real_
  w
}

# A fit's `weights =` argument evaluated on the rows of `data` -- a bare glm() or
# glmer() fitted with weights = n -- when every variable it names is a column of
# `data`, else NULL. The formula's environment supplies only the FUNCTIONS the
# expression applies (1 / v, sqrt(n)), never a variable: a same-named object reachable
# from there has nothing to do with these rows. That includes a fit_maihda() model's
# own record, weights = .maihda_arg_weights, which is bound in its formula's
# environment to the TRAINING rows' weights: evaluated there, it handed new rows of the
# same count the training rows' weights, the very pairing this replaces.
maihda_prediction_panel_call_weights <- function(model, data) {
  cl <- tryCatch(stats::getCall(model), error = function(e) NULL)
  expr <- if (is.call(cl)) cl$weights
  if (is.null(expr)) {
    return(NULL)
  }
  vars <- all.vars(expr)
  if (length(vars) == 0L || !all(vars %in% names(data))) {
    return(NULL)
  }
  env <- tryCatch(environment(stats::formula(model)), error = function(e) NULL)
  if (!is.environment(env)) env <- baseenv()
  w <- tryCatch(eval(expr, data, env), error = function(e) NULL)
  maihda_prediction_panel_clean_weights(w, nrow(data))
}

# The value the panel plots on each FITTED row, on the scale it plots it: what
# maihda_prediction_panel_fitted() returns when there is no `data`. Only for the fits
# whose prior weights live on their fitted rows (lme4, glm, lm); NULL otherwise.
maihda_prediction_panel_fit_row_values <- function(model, type) {
  if (!inherits(model, c("merMod", "lm"))) {
    return(NULL)
  }
  v <- tryCatch(
    if (isTRUE(type %in% c("binomial", "poisson"))) {
      maihda_fit_rows_predict(model, type = "response")
    } else {
      maihda_fit_rows_predict(model)
    },
    error = function(e) NULL)
  if (is.numeric(v)) as.numeric(v) else NULL
}

# Which fitted row each row of `data` is -- an index into `frame`, NA for a row that is
# none of them -- or NULL when that cannot be established. A row name is only a claim:
# new data numbered 1..n carries the row names of a fit to OTHER rows numbered 1..n.
# So a named row counts only when it reproduces that fitted row -- the same plotted
# value (`values` against `fit_values`, the model's own summary of its covariates,
# stratum and offset) and, where both have one, the same stratum -- and one named row
# that does not refutes the names of every row. A row whose name is no fitted row (a
# new row, or one the fit dropped) is NA. The check cannot tell apart two fitted rows
# with the same prediction and stratum; exchanging their weights would leave every
# stratum's weighted mean of the predictions unchanged.
maihda_prediction_panel_row_identity <- function(frame, data, fit_values, values) {
  if (!is.data.frame(frame) || is.null(fit_values) || is.null(values) ||
      length(fit_values) != nrow(frame) || length(values) != nrow(data)) {
    return(NULL)
  }
  idx <- match(rownames(data), rownames(frame))
  hit <- which(!is.na(idx))
  if (length(hit) == 0L) {
    return(idx)
  }
  a <- as.numeric(values[hit])
  b <- as.numeric(fit_values[idx[hit]])
  same <- is.finite(a) & is.finite(b) & abs(a - b) <= 1e-8 * pmax(1, abs(b))
  if ("stratum" %in% names(data) && "stratum" %in% names(frame)) {
    sa <- as.character(data$stratum[hit])
    sb <- as.character(frame$stratum[idx[hit]])
    same <- same & !is.na(sa) & !is.na(sb) & sa == sb
  }
  if (!all(same)) {
    return(NULL)
  }
  idx
}

# The binomial trial count of each row of `data` when the fit's response carries one:
# a brms `y | trials(n)` term evaluated on the rows, or the successes plus failures of
# a cbind() response (NA where the row's outcome is missing). NULL for any other fit.
maihda_prediction_panel_row_trials <- function(model, data, frame, resp = NULL) {
  brms_trials <- tryCatch(maihda_prediction_panel_brms_trials(model, data),
                          error = function(e) NA_real_)
  if (!is.null(brms_trials)) {
    return(rep_len(as.numeric(brms_trials), nrow(data)))
  }
  fit_resp <- if (is.data.frame(frame)) {
    tryCatch(stats::model.response(frame), error = function(e) NULL)
  }
  if (is.null(dim(fit_resp))) {
    return(NULL)
  }
  if (is.null(resp)) resp <- maihda_prediction_panel_row_response(model, data)
  if (is.matrix(resp) && ncol(resp) == 2L && nrow(resp) == nrow(data)) {
    return(as.numeric(resp[, 1]) + as.numeric(resp[, 2]))
  }
  rep(NA_real_, nrow(data))
}

# The observed response of each row of `data` as the model's formula defines it: a
# column of `data`, or the fitted expression -- cbind(s, f), I(y > 3) -- evaluated on
# its rows (maihda_response_from_frame()). NULL when `data` cannot supply it.
maihda_prediction_panel_row_response <- function(model, data) {
  f <- tryCatch(stats::formula(model), error = function(e) NULL)
  if (is.null(f)) {
    return(NULL)
  }
  r <- tryCatch(maihda_response_from_frame(data, f), error = function(e) NULL)
  if (is.null(r)) {
    return(NULL)
  }
  if (is.null(dim(r))) unname(r) else r
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

# The name of the observed-response column the panels look up in `data`. formula()
# of a brmsfit is a brmsformula -- a list, whose as.character() deparses its
# elements, so element 2 is "list()" and never named a column: the brms ordinal
# surprise was undefined for every row and the brms binomial deviance residuals all
# 0. Unwrap it and follow brms addition terms (y | weights(w), y | thres(K)) to the
# response on their left. A trials() addition makes the response an aggregated
# count, not an observation-level 0/1 outcome, so NULL is returned for it and the
# panels treat the response as unobserved, as before. Every other model keeps the
# as.character(formula)[2] lookup.
maihda_prediction_panel_response_name <- function(model) {
  form <- tryCatch(formula(model), error = function(e) NULL)
  if (is.null(form)) {
    return(NULL)
  }
  if (inherits(form, "brmsformula")) {
    f <- form$formula
    if (!inherits(f, "formula") || length(f) != 3L) {
      return(NULL)
    }
    lhs <- f[[2]]
    if (is.call(lhs) && identical(lhs[[1]], as.name("|")) &&
        !is.null(maihda_find_trials_expr(lhs[[3]]))) {
      return(NULL)
    }
    return(paste(deparse(maihda_describe_response_expr(f)), collapse = " "))
  }
  as.character(form)[2]
}

# TRUE when the binomial panel can draw the observed outcome `x` as point shapes: a
# vector -- not the two-column cbind(successes, failures) matrix -- of categories, as a
# 0/1 or other categorical outcome is (whole numbers when numeric, the rule
# maihda_is_binary_vector() applies), and no more of them than ggplot's shape palette
# holds (six). A numeric outcome with a fractional value is a proportion with trial
# weights -- an aggregated outcome, which no shape stands for.
maihda_prediction_panel_shape_outcome <- function(x) {
  if (is.null(x) || !is.null(dim(x))) {
    return(FALSE)
  }
  seen <- unique(x[!is.na(x)])
  length(seen) <= 6L &&
    (!is.numeric(seen) || all(is.finite(seen) & abs(seen - round(seen)) < 1e-8))
}

# The formula of a brms fit, unwrapped from its brmsformula; NULL for any other model.
maihda_prediction_panel_brms_formula <- function(model) {
  if (!inherits(model, "brmsfit")) {
    return(NULL)
  }
  f <- tryCatch(formula(model), error = function(e) NULL)
  if (inherits(f, "brmsformula")) {
    f <- f$formula
  }
  if (!inherits(f, "formula") || length(f) != 3L) {
    return(NULL)
  }
  f
}

# Trial counts of a brms `y | trials(n)` fit, one per row of `data`; NULL when the
# model is not a brms fit with a trials() term. brms's fitted() and posterior_epred()
# return the expected success COUNT (trials * p) for such a fit, so the binomial
# panel divides by these to plot probabilities, as predict_maihda() does. The term is
# evaluated on `data`, the rows the panel predicts; one that cannot be evaluated there
# is an error rather than a silent return of counts.
maihda_prediction_panel_brms_trials <- function(model, data) {
  f <- maihda_prediction_panel_brms_formula(model)
  if (is.null(f)) {
    return(NULL)
  }
  lhs <- f[[2]]
  if (!is.call(lhs) || !identical(lhs[[1]], as.name("|")) ||
      is.null(maihda_find_trials_expr(lhs[[3]]))) {
    return(NULL)
  }
  trials <- maihda_trials_from_formula(f, data)
  if (is.null(trials)) {
    stop("Could not evaluate the trials() term of this brms binomial fit on `data`, ",
         "so its expected success counts cannot be put on the probability scale.",
         call. = FALSE)
  }
  trials
}

# Success counts of a brms `y | trials(n)` fit -- its response, left of the addition
# terms -- evaluated on `data`; NULL when `data` cannot supply them.
maihda_prediction_panel_brms_successes <- function(model, data) {
  f <- maihda_prediction_panel_brms_formula(model)
  if (is.null(f)) {
    return(NULL)
  }
  s <- tryCatch(eval(maihda_describe_response_expr(f), envir = data, enclos = baseenv()),
                error = function(e) NULL)
  if (!is.numeric(s) || !is.null(dim(s)) || length(s) != nrow(data)) {
    return(NULL)
  }
  as.numeric(s)
}

# Response-scale values divided by their rows' trial counts. A row without a positive
# trial count has no per-trial probability and becomes NA.
maihda_prediction_panel_per_trial <- function(x, trials) {
  ok <- is.finite(trials) & trials > 0
  out <- rep(NA_real_, length(x))
  out[ok] <- x[ok] / trials[ok]
  out
}

maihda_prediction_panel_fitted <- function(model, data, type, fitted_data = FALSE,
                                           trials = NULL, maihda_obj = NULL) {
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
    est <- as.numeric(fit[, "Estimate"])
    se <- as.numeric(se)
    if (!is.null(trials)) {
      # A `y | trials(n)` fit: the estimate and its interval are success counts.
      est <- maihda_prediction_panel_per_trial(est, trials)
      se <- maihda_prediction_panel_per_trial(se, trials)
    }
    return(list(fit = est, se.fit = se))
  }

  # wemix: predict.WeMixResults() takes the response off the CALL's formula with
  # as.name(form[[2]]), which is only defined when the response is a bare variable --
  # handed a call it raises "'language' object cannot be coerced to type 'symbol'". So
  # every wemix fit with a transformed response (log(y), sqrt(y), I(y > 0)) errored
  # here, and plot(type = "all") reported the panel as uncomputable and dropped it,
  # while the same fit spelled with a pre-computed column worked: one fit, two
  # spellings, two answers. Rebuild the linear predictor from the wrapper instead --
  # the route predict_maihda() already takes for this engine, which also re-uses the
  # FITTED transformation basis and factor coding rather than re-deriving them from the
  # prediction rows -- and leave predict() only to a bare fit, which carries no wrapper.
  # It agrees with predict.WeMixResults() to the last bit on a bare-response fit's own
  # rows, so nothing that worked before changes value. WeMix exposes no prediction SE
  # on either route, so se.fit stays NA and the case-level bars are omitted rather than
  # faked (see the NA_real_ note below).
  if (inherits(model, "WeMixResults")) {
    on_response <- type == "binomial" || type == "poisson"
    if (!is.null(maihda_obj)) {
      eta <- as.numeric(maihda_wemix_linpred(
        maihda_obj,
        newdata = if (isTRUE(fitted_data)) NULL else data,
        include_re = TRUE))
      fit <- if (on_response) {
        as.numeric(maihda_linkinv(maihda_model_family(maihda_obj))(eta))
      } else {
        eta
      }
      return(list(fit = fit, se.fit = rep(NA_real_, length(fit))))
    }
    if (isTRUE(fitted_data)) {
      # The fit's own rows: predict() WITHOUT newdata returns the stored linear
      # predictor / mean and never reaches the as.name() above, so a transformed
      # response is fine here.
      fit <- as.numeric(stats::predict(model,
                                       type = if (on_response) "response" else "link"))
      return(list(fit = fit, se.fit = rep(NA_real_, length(fit))))
    }
    # Prediction data with no wrapper to rebuild the design from leaves
    # predict.WeMixResults() as the only route -- the one that cannot take a
    # transformed response. A bare WeMix fit also cannot reach the fitted-rows branch
    # above, because model.frame() of one does not evaluate, so 'data' is not
    # optional here: name the two remedies that exist rather than pass WeMix's
    # as.name() error on.
    lhs <- tryCatch(stats::formula(stats::getCall(model))[[2]], error = function(e) NULL)
    if (!is.null(lhs) && !is.symbol(lhs)) {
      stop("WeMix's predict() method cannot predict from a fit whose response is an ",
           "expression (", paste(deparse(lhs), collapse = " "), ") rather than a bare ",
           "variable. Pass the fit_maihda() model object instead of the bare WeMix fit, ",
           "or refit with the transformed response stored as a column of 'data'.",
           call. = FALSE)
    }
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
        maihda_fit_rows_predict(model, type = "response")
      } else {
        maihda_fit_rows_predict(model)
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

maihda_prediction_panel_ordinal_probs <- function(model, data, coding = NULL) {
  if (inherits(model, "clmm")) {
    # predict.clmm does not exist: rebuild the location eta = x'beta + u from
    # the stored components (fixed-effects-only $terms, $beta, the fit's factor
    # coding, and the stratum conditional modes) and difference the cumulative
    # probabilities. Including the random effect matches the other branches of this
    # panel, whose predict() calls include random effects by default. `coding` is
    # the record a maihda_model carries; for a bare clmm it is read off the fit's
    # own $xlevels / $contrasts (see maihda_engine_fixed_coding()).
    maihda_require_ordinal()
    tt <- stats::delete.response(model$terms)
    if (is.null(coding)) {
      coding <- maihda_engine_fixed_coding(model, model$terms, data)
    }
    design <- maihda_fixed_design(tt, data, coding)
    mf <- design$frame
    X <- design$X
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
    probs <- maihda_ordinal_category_probs(eta, maihda_clmm_cutpoints(model),
                                           model$link)
    return(as.data.frame(probs))
  }

  probs <- NULL
  if (inherits(model, "brmsfit")) {
    # The posterior MEAN of each category probability, from fitted(). predict()
    # only estimates it by tabulating simulated responses, so the panel changed
    # between calls (by 0.14 in probability on a 300-draw fit) and a category no
    # draw produced got probability 0, which would be an infinite surprise. These are
    # also the probabilities predict_maihda() turns into expected category scores.
    if (!requireNamespace("brms", quietly = TRUE)) {
      stop("Package 'brms' is required to plot prediction deviations from brms models.",
           call. = FALSE)
    }
    fit <- stats::fitted(model, newdata = data, summary = TRUE)
    dims <- dim(fit)
    if (length(dims) != 3L || !"Estimate" %in% dimnames(fit)[[2]]) {
      stop("Could not extract category probabilities from the brms ordinal model.",
           call. = FALSE)
    }
    # An explicit nobs x ncat matrix: fit[, "Estimate", ] drops the unit row margin
    # of a one-row `data` (see maihda_brms_fitted_array_scores()).
    probs <- matrix(fit[, "Estimate", ], nrow = dims[1], ncol = dims[3],
                    dimnames = list(NULL, dimnames(fit)[[3]]))
  }

  if (is.null(probs)) {
    probs <- tryCatch(
      predict(model, newdata = data, type = "probs"),
      error = function(e) NULL
    )
  }
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

# The fitted response categories of an ordinal model, one per column of its
# probability matrix `probs` and in column order: a clmm fit's y.levels; the levels
# of a brms fit's stored response (brms numbers an ordered factor by its position in
# those levels, and an integer response is its own category number); otherwise the
# column names the model's predict() method gave, which polr sets to its levels.
# Read from the FIT, never from the prediction `data`, whose response may lack a
# category or declare its levels in another order. A clmm or brms fit whose
# categories cannot be read is an error, not a fall-back to its column names: those
# are positions (clmm's 1..K), and a label match against them is the defect this
# replaced.
maihda_prediction_panel_ordinal_categories <- function(model, probs) {
  k <- ncol(probs)
  cats <- if (inherits(model, "clmm")) {
    model$y.levels
  } else if (inherits(model, "brmsfit")) {
    resp <- maihda_prediction_panel_response_name(model)
    y <- if (!is.null(resp) && is.data.frame(model$data) &&
             resp %in% names(model$data)) {
      model$data[[resp]]
    }
    if (is.factor(y)) {
      levels(y)
    } else if (is.numeric(y)) {
      seq_len(k)
    }
  } else {
    colnames(probs)
  }
  cats <- as.character(cats)
  if (length(cats) != k || anyNA(cats) || anyDuplicated(cats) > 0L) {
    stop(sprintf(paste0("Could not match the ordinal model's %d category ",
                        "probabilities to its fitted response categories."), k),
         call. = FALSE)
  }
  cats
}

# Probability of each row's observed category. The category is found by its
# POSITION among `categories` -- the fitted categories, in the column order of
# `prob_mat` -- never by matching its label to a column name: the engines name those
# columns 1..K (the rebuilt clmm matrix), "P(Y = <category>)" (brms) or by level
# (polr), so on a clmm fit a label match scored no row of a low/mid/high outcome, lost
# category 0 of a 0/1/2 one and read its other rows one category off, and swapped the
# end categories of a 3/2/1 one. `obs_cat` is the observed response (NA
# where missing, which stays NA); `resp_name` / `resp_found` name its column and say
# whether `data` has it. A row whose category is not a fitted one cannot be scored
# and is warned about, as is a `data` in which no row can be scored, so the
# surprise panel is never silently empty.
maihda_prediction_panel_observed_prob <- function(prob_mat, obs_cat, categories,
                                                  resp_name = NULL,
                                                  resp_found = TRUE) {
  obs_cat <- as.character(obs_cat)
  pos <- match(obs_cat, categories)
  hit <- which(!is.na(pos))
  out <- rep(NA_real_, length(obs_cat))
  out[hit] <- prob_mat[cbind(hit, pos[hit])]

  unmatched <- !is.na(obs_cat) & is.na(pos)
  fitted_txt <- paste(categories, collapse = ", ")
  if (length(hit) == 0L && length(obs_cat) > 0L) {
    reason <- if (!resp_found) {
      if (is.null(resp_name)) {
        "the model's response could not be identified"
      } else {
        sprintf("`data` has no column '%s' with the observed response", resp_name)
      }
    } else if (!any(unmatched)) {
      "every observed response in `data` is missing"
    } else {
      sprintf("no observed category in `data` is one of the model's fitted categories (%s)",
              fitted_txt)
    }
    warning("plot_prediction_deviation_panels(): no row can be scored, so the ordinal ",
            "surprise panel is empty: ", reason, ". Use ordinal_mode = ",
            "\"expected_score\" to rank by the predictions alone.", call. = FALSE)
  } else if (any(unmatched)) {
    warning(sprintf(paste0("plot_prediction_deviation_panels(): %d row(s) of `data` have ",
                           "an observed category that is not one of the model's fitted ",
                           "categories (%s), e.g. '%s'; their surprise is undefined and ",
                           "they are left out of the surprise panel."),
                    sum(unmatched), fitted_txt, obs_cat[unmatched][1]),
            call. = FALSE)
  }
  out
}

# TRUE when a fitted binomial response codes ONE trial per row: a factor, a logical,
# or numbers that are all 0 or 1. A proportion strictly between them is a share of
# several trials.
maihda_prediction_panel_is_01 <- function(y) {
  if (is.factor(y) || is.logical(y)) {
    return(TRUE)
  }
  is.numeric(y) && all(y[!is.na(y)] %in% c(0, 1))
}

# The successes and trials of each plotted row of a glm / lme4 binomial fit whose rows
# are not single, unit-weight Bernoulli trials -- a cbind(successes, failures)
# response, a proportion, or a 0/1 outcome with prior weights -- and NULL for a fit
# whose rows are, and for every other engine; those keep the Bernoulli residual.
#
# For the binomial family a prior weight IS a trial count (see
# maihda_da_aggregated_counts()): R forms a row's deviance from its proportion y and
# weight w as 2 w [y log(y / p) + (1 - y) log((1 - y) / (1 - p))], the deviance of
# y w successes out of w trials, whole numbers or not, and a cbind() response enters
# as y = s / (s + f) with w = s + f times any prior weight. On the fitted rows that is
# the engine's own residuals(type = "deviance"), returned as `residuals`. On any other
# rows it is rebuilt from THOSE rows' response and prior weights, and is unknown where
# either is. The panel used to borrow the fitted rows' residuals whenever `data` had
# as many rows -- up to 2.83 off on a row once the fitted rows of
# test-audit-2026-09-27b.R's fixture were reordered -- and to score a 0/1 outcome as a
# single trial whatever its weight, exactly five times too small at 25 trials a row.
#
# `resp` is the rows' observed response (maihda_prediction_panel_row_response()),
# `obs_outcome_01` its 0/1 coding where it has one, and `prior` the rows' prior
# weights (NULL when the fit has none).
maihda_prediction_panel_binomial_counts <- function(model, data, frame, fitted_data,
                                                    resp, obs_outcome_01, prior) {
  if (!inherits(model, c("merMod", "glm")) || !is.data.frame(frame)) {
    return(NULL)
  }
  fit_resp <- tryCatch(stats::model.response(frame), error = function(e) NULL)
  if (is.null(fit_resp)) {
    return(NULL)
  }
  aggregated <- !is.null(dim(fit_resp))
  if (!aggregated && is.null(maihda_prediction_panel_fit_prior(frame)) &&
      maihda_prediction_panel_is_01(fit_resp)) {
    return(NULL)
  }
  n <- nrow(data)
  if (isTRUE(fitted_data)) {
    r <- tryCatch(abs(maihda_fit_rows_residuals(model, type = "deviance")),
                  error = function(e) NULL)
    if (is.numeric(r) && length(r) == n) {
      return(list(residuals = r))
    }
  }
  if (is.null(prior)) prior <- rep(1, n)
  if (aggregated) {
    ok <- is.matrix(resp) && ncol(resp) == 2L && nrow(resp) == n
    s <- if (ok) as.numeric(resp[, 1]) else rep(NA_real_, n)
    trials <- if (ok) s + as.numeric(resp[, 2]) else rep(NA_real_, n)
  } else {
    # The 0/1 coding where the outcome has one; otherwise the numeric proportion, which
    # also scores a 0/1 outcome that takes a single value on these rows.
    s <- as.numeric(obs_outcome_01)
    num <- if (is.numeric(resp) && is.null(dim(resp)) && length(resp) == n) {
      as.numeric(resp)
    } else {
      rep(NA_real_, n)
    }
    fill <- is.na(s) & is.finite(num) & num >= 0 & num <= 1
    s[fill] <- num[fill]
    trials <- rep(1, n)
  }
  list(successes = prior * s, trials = prior * trials)
}

# |Binomial deviance residual| of each plotted row, NA where the row cannot be scored.
maihda_prediction_panel_binomial_residuals <- function(model, data, fitted, obs_outcome_01,
                                                       trials = NULL, counts = NULL) {
  # A brms `y | trials(n)` fit: each row's successes out of its trials, the residual
  # lme4's residuals() gives a cbind() fit. `fitted` is already per trial.
  if (!is.null(trials)) {
    successes <- maihda_prediction_panel_brms_successes(model, data)
    if (!is.null(successes)) {
      return(maihda_binomial_abs_deviance_residual_agg(successes, trials, fitted))
    }
  }
  # A glm / lme4 fit whose rows carry several trials, or a proportion
  # (maihda_prediction_panel_binomial_counts()).
  if (!is.null(counts$residuals)) {
    return(counts$residuals)
  }
  if (!is.null(counts)) {
    return(maihda_binomial_abs_deviance_residual_agg(counts$successes, counts$trials,
                                                     fitted))
  }
  if (length(obs_outcome_01) != length(fitted)) {
    return(rep(NA_real_, length(fitted)))
  }
  maihda_binomial_abs_deviance_residual(obs_outcome_01, fitted)
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
                                                    level = 0.95, trials = NULL) {
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
  if (!is.null(trials)) {
    # A `y | trials(n)` fit draws success counts: divide each row's draws by its
    # trials. A row without a positive trial count is zeroed and given no weight, so
    # it cannot turn its stratum's per-draw mean into NA.
    ok <- is.finite(trials) & trials > 0
    ep <- sweep(ep, 2L, ifelse(ok, trials, 1), "/")
    ep[, !ok] <- 0
    if (is.null(weight) || length(weight) != length(ok)) {
      weight <- rep(1, length(ok))
    }
    weight <- ifelse(ok, weight, 0)
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
#'     i.e. where the observed 0/1 outcome -- for an aggregated binomial
#'     (\code{cbind(successes, failures)}, a proportion with its trial counts as
#'     prior weights, or a brms \code{y | trials(n)} fit), the observed successes
#'     out of trials -- is least consistent with the fitted probability (worst-fit
#'     points), ranked by \eqn{|deviance residual|}. As in R's binomial family, a
#'     prior weight counts as that many trials, so a 0/1 outcome with weights is
#'     scored as that many successes or failures. Every row's outcome is coded as
#'     the model coded it -- by label, whatever the order of a factor's levels in
#'     \code{data} -- and a factor with more than two levels, as in R's binomial
#'     family, as its first level against all the others. Predictions are per-trial
#'     probabilities for every binomial fit. A row whose outcome or trial count is
#'     unknown has no residual: it is left out of the stratum residual means, is
#'     never labelled, and is drawn as a hollow point.
#'   \item Ordinal \code{"surprise"} mode: the cases/strata with the highest
#'     surprise \eqn{-\log P(\text{observed category})}, i.e. the least probable
#'     observations under the model. It needs the observed response in
#'     \code{data}: a row whose response is missing is left out, and one whose
#'     category is not among the model's fitted categories is left out with a
#'     warning.
#' }
#'
#' @param model A fitted model object (e.g., from `lm()`, `glm()`, `MASS::polr()`, or `lme4::glmer()`).
#' @param data The original data frame used to fit the model. If `NULL`, attempts
#'   to extract from the model. When supplied, everything is computed for its own
#'   rows: their predictions, their binomial residuals (from their own outcomes and
#'   trial counts) and their weights in the stratum summaries. A row's weight is
#'   read from its column in `data` -- the sampling-weight column, or the column a
#'   bare model's `weights =` argument names (`glm()`, `lmer()`, `MASS::polr()`,
#'   ...) -- or else from the fitted row with the same row name and prediction (a
#'   `fit_maihda()` model stores its weights as values). If any plotted row's
#'   weight cannot be found, every row is weighted equally, with a warning.
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
  # unweighted fit, so the weighted means below reduce to plain means) are taken in
  # each branch once its predictions exist: on supplied rows a fitted row is
  # recognised by its prediction (maihda_prediction_panel_prior_weights()).

  # The rows handed to brms's own prediction calls. brms validates them against its
  # fitted data, so supplied rows lose an outcome still in its original coding
  # (maihda_brms_newdata(), the copy predict_maihda() hands brms too); `data` itself,
  # original outcome included, stays for the panel's own summaries. Every other model
  # predicts from `data` as it stands.
  pred_data <- if (data_supplied && inherits(model, "brmsfit")) {
    maihda_brms_newdata(maihda_obj, data)
  } else {
    data
  }

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
    preds <- maihda_prediction_panel_fitted(model, pred_data, type,
                                            fitted_data = !data_supplied,
                                            maihda_obj = maihda_obj)
    prior_w <- as.vector(maihda_prediction_panel_prior_weights(
      maihda_obj, model, data, fitted_data = !data_supplied,
      values = preds$fit, type = type))

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
      strat_ci <- maihda_prediction_stratum_interval_brms(model, pred_data, prior_w)
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
    # A brms `y | trials(n)` fit predicts expected success COUNTS; its trial counts
    # put them, and their posterior draws, on the probability scale.
    trials <- maihda_prediction_panel_brms_trials(model, data)
    preds <- maihda_prediction_panel_fitted(model, pred_data, "binomial",
                                            fitted_data = !data_supplied,
                                            trials = trials, maihda_obj = maihda_obj)
    prior_w <- maihda_prediction_panel_prior_weights(
      maihda_obj, model, data, fitted_data = !data_supplied,
      values = preds$fit, type = "binomial")
    residuals_warned <- isTRUE(attr(prior_w, "residuals_warned"))
    prior_w <- as.vector(prior_w)

    # Try to extract response variable. A cbind(successes, failures) response sits in
    # the model's data as a two-column matrix under its deparsed name: an aggregated
    # count, not one 0/1 outcome per row (as.factor() of it has two entries per row),
    # so it has no 0/1 coding and is scored from its successes and trials below.
    resp_name <- maihda_prediction_panel_response_name(model)
    obs_outcome <- NULL
    obs_outcome_01 <- rep(NA_integer_, nrow(data))
    unmatched <- integer(0)
    raw_outcome <- if (!is.null(resp_name) && resp_name %in% names(data)) {
      data[[resp_name]]
    }
    # An outcome written as an expression -- I(y > 3), cbind(s, f) -- is a column of
    # the model frame but not of the original data: evaluate it on the plotted rows.
    if (is.null(raw_outcome) && is.null(trials)) {
      raw_outcome <- maihda_prediction_panel_row_response(model, data)
    }
    if (!is.null(raw_outcome) && is.null(dim(raw_outcome))) {
      obs_outcome <- as.factor(raw_outcome)
      # Every row is coded as the model coded its outcome, not by the plotted rows' own
      # levels or values -- the fitted rows too. R's binomial family takes a factor's
      # first level as the failure and EVERY other level as the event, which a factor
      # with three or more levels cannot be read back from by its values
      # (maihda_binomial_observed_01() codes only two): an unweighted glm, glmer or
      # fit_maihda() fit of one had no row scored on its own rows, where the same rows
      # passed as `data` were scored in full.
      obs_outcome_01 <- if (length(raw_outcome) == nrow(data)) {
        maihda_prediction_panel_fit_01(maihda_obj, model, raw_outcome)
      } else {
        maihda_binomial_observed_01(raw_outcome, nrow(data))
      }
      unmatched <- as.integer(attr(obs_outcome_01, "unmatched"))
      obs_outcome_01 <- as.vector(obs_outcome_01)
    }
    if (is.null(obs_outcome)) {
      obs_outcome <- factor(rep(NA, nrow(data)))
    }

    # Rows carrying several trials, or a proportion, are scored from their own
    # successes and trials (maihda_prediction_panel_binomial_counts()).
    fit_rows <- maihda_prediction_panel_fit_rows(maihda_obj, model)
    counts <- maihda_prediction_panel_binomial_counts(
      model, data, fit_rows, fitted_data = !data_supplied, resp = raw_outcome,
      obs_outcome_01 = obs_outcome_01,
      prior = if (data_supplied) {
        maihda_prediction_panel_row_prior(model, data, fit_rows, preds$fit, "binomial")
      } else {
        maihda_prediction_panel_fit_prior(fit_rows)
      })
    resids <- maihda_prediction_panel_binomial_residuals(
      model, data, preds$fit, obs_outcome_01, trials = trials, counts = counts
    )
    # A supplied outcome the fit never saw (another label, a stray code) cannot be
    # coded; say so rather than drop its rows unannounced. Only the rows that are in
    # fact left unscored: a weighted fit still scores a proportion as successes out of
    # its trials (maihda_prediction_panel_binomial_counts()).
    unmatched <- unmatched[!is.finite(resids[unmatched])]
    unmatched_txt <- if (length(unmatched)) {
      sprintf("the outcome of %d row(s) of `data` is not one the model was fitted on (e.g. '%s')",
              length(unmatched), as.character(raw_outcome[unmatched[1]]))
    }
    if (nrow(data) > 0L && !any(is.finite(resids)) && !residuals_warned) {
      outcome_found <- !is.null(raw_outcome) ||
        (!is.null(trials) && !is.null(maihda_prediction_panel_brms_successes(model, data)))
      warning("plot_prediction_deviation_panels(): no row of `data` can be scored, so no ",
              "deviance residual is shown and no point is labelled: ",
              if (!outcome_found) {
                "`data` does not contain the model's observed outcome"
              } else if (!is.null(unmatched_txt)) {
                unmatched_txt
              } else {
                "no row has an observed outcome the model can score beside a finite prediction"
              }, ".", call. = FALSE)
    } else if (!is.null(unmatched_txt)) {
      warning("plot_prediction_deviation_panels(): ", unmatched_txt, "; those rows cannot ",
              "be scored and are left out of the residual summaries and labels.",
              call. = FALSE)
    }

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
      strat_ci <- maihda_prediction_stratum_interval_brms(model, pred_data, prior_w,
                                                          trials = trials)
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

    # A row or stratum without a deviance residual (its outcome or trials unknown) is
    # never labelled as a worst fit, and its point -- when its prediction is known -- is
    # drawn hollow at a fixed size rather than dropped for want of a size. The sized
    # layer is narrowed to the scored rows only when there is such a point to draw
    # hollow, so a panel without one is drawn exactly as before.
    scored <- is.finite(df$abs_res_dev)
    unscored_df <- df[!scored & is.finite(df$fitted), , drop = FALSE]
    scored_df <- if (nrow(unscored_df) > 0L) df[scored, , drop = FALSE]
    label_df <- (if (all(scored)) df else df[scored, , drop = FALSE]) |>
      dplyr::arrange(dplyr::desc(.data$abs_res_dev)) |> utils::head(top_n_labels)

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
        ggplot2::geom_point(data = scored_df, ggplot2::aes(color = .data$direction, size = .data$abs_res_dev), alpha = 0.8) +
        ggplot2::labs(
          x = x_label, y = "Predicted Probability", color = "Direction", size = "|Deviance\nResidual|"
        )
    } else {
      # The point shape shows each case's observed outcome. An aggregated binomial has
      # none to show -- successes out of trials in a cbind() matrix, a proportion with
      # trial weights, a brms trials() count -- and neither has a case whose outcome
      # `data` does not carry. Mapped to shape, those came out NA (a proportion's values
      # past the palette's six likewise) and ggplot dropped every such point, so they
      # take a fixed shape: the default when no case shows an outcome, and a hollow
      # ring beside cases that do, so as not to pass for an observed category.
      shape_df <- if (is.null(scored_df)) df else scored_df
      shown <- maihda_prediction_panel_shape_outcome(raw_outcome) & !is.na(shape_df$obs_outcome)
      if (all(shown)) {
        p2 <- p2 +
          ggplot2::geom_point(data = scored_df, ggplot2::aes(color = .data$direction, size = .data$abs_res_dev, shape = .data$obs_outcome), alpha = 0.8)
      } else if (!any(shown)) {
        p2 <- p2 +
          ggplot2::geom_point(data = scored_df, ggplot2::aes(color = .data$direction, size = .data$abs_res_dev), alpha = 0.8)
      } else {
        p2 <- p2 +
          ggplot2::geom_point(data = shape_df[shown, , drop = FALSE], ggplot2::aes(color = .data$direction, size = .data$abs_res_dev, shape = .data$obs_outcome), alpha = 0.8) +
          ggplot2::geom_point(data = shape_df[!shown, , drop = FALSE], ggplot2::aes(color = .data$direction, size = .data$abs_res_dev),
                              shape = 1, alpha = 0.8, show.legend = FALSE)
      }

      if (any(df$wrong == "Wrong", na.rm = TRUE)) {
        p2 <- p2 + ggplot2::geom_point(data = dplyr::filter(df, .data$wrong == "Wrong", is.finite(.data$abs_res_dev)), shape = 1, color = "red", ggplot2::aes(size = .data$abs_res_dev + 0.5))
      }
      # No shape legend to title when no case shows an outcome (labs() drops a waiver();
      # ggplot announces a label for an unmapped aesthetic as unknown).
      p2 <- p2 + ggplot2::labs(x = x_label, y = "Predicted Probability", color = "Direction", size = "|Deviance\nResidual|",
                               shape = if (any(shown)) "Observed" else ggplot2::waiver())
    }
    if (nrow(unscored_df) > 0L) {
      p2 <- p2 + ggplot2::geom_point(data = unscored_df, ggplot2::aes(color = .data$direction),
                                     shape = 1, size = 1.5, show.legend = FALSE)
    }

    p2 <- p2 +
      ggplot2::geom_hline(ggplot2::aes(yintercept = .data$mean_fitted[1]), linetype = "dashed") +
        ggrepel::geom_label_repel(data = label_df, ggplot2::aes(label = .data$id), size = 3, min.segment.length = 0) +
      ggplot2::scale_color_manual(values = c("Above Mean" = "#0072B2", "Below Mean" = "#D55E00")) +
      theme_maihda()

    return(patchwork::wrap_plots(p1, p2, ncol = 1, heights = c(1, 2)))

  } else if (type == "ordinal") {
    # ORDINAL LOGIC
    probs <- maihda_prediction_panel_ordinal_probs(
      model, pred_data,
      coding = if (!is.null(maihda_obj) && inherits(model, "clmm")) {
        maihda_object_fixed_coding(maihda_obj)
      })
    prob_mat <- as.matrix(probs)
    # No fitted-row recognition here (it needs the lme4 / glm predictions): of the
    # ordinal fits only a bare polr / clm / clmm carries prior weights, read from its
    # "(weights)" frame column on its own rows and from its `weights =` column on
    # supplied rows (the engine = "ordinal" maihda_model refuses weights).
    prior_w <- as.vector(maihda_prediction_panel_prior_weights(
      maihda_obj, model, data, fitted_data = !data_supplied, type = "ordinal"))

    if (ordinal_mode == "surprise") {
      # The fitted categories label the probability columns for display and locate
      # each row's observed category by position. The columns themselves are keyed
      # prob_1..prob_K: category labels are arbitrary strings, and as column names a
      # category called "n", "weight" or "id" was overwritten by this panel's own
      # columns (a polr category "n" was drawn as each stratum's row count).
      categories <- maihda_prediction_panel_ordinal_categories(model, probs)
      prob_cols <- paste0("prob_", seq_len(ncol(prob_mat)))
      colnames(prob_mat) <- prob_cols

      resp_name <- maihda_prediction_panel_response_name(model)
      resp_found <- !is.null(resp_name) && resp_name %in% names(data)
      obs_cat <- rep(NA, nrow(data))
      if (resp_found) obs_cat <- as.character(data[[resp_name]])

      df <- as.data.frame(prob_mat) |>
        dplyr::mutate(
          id = dplyr::row_number(),
          obs_cat = obs_cat
        )

      k_seq <- seq_len(ncol(prob_mat))
      df$expected_score <- rowSums(prob_mat * matrix(k_seq, nrow = nrow(prob_mat), ncol = ncol(prob_mat), byrow = TRUE))

      # Probability of the observed category, located by its position among the
      # fitted categories (see maihda_prediction_panel_observed_prob()).
      df$observed_prob <- maihda_prediction_panel_observed_prob(
        prob_mat, obs_cat, categories, resp_name = resp_name, resp_found = resp_found
      )

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
        dplyr::mutate(Category = factor(.data$Category, levels = prob_cols,
                                        labels = categories))

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
