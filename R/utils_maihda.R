# Internal helpers shared across model summaries, PCV, predictions, and plots.

# Map engine-specific family labels onto the package's canonical names so every
# downstream `fam$family ==` / `%in%` comparison sees one spelling per family.
# The one label that genuinely needs this is lme4's negative binomial: a
# glmer.nb() fit reports family "Negative Binomial(<theta>)" with the estimated
# theta embedded in the STRING, so the label differs between any two fits (e.g.
# a null and an adjusted model) and never matches a fixed name. brms already
# reports "negbinomial", which is adopted as the canonical name. "ordinal" is
# accepted as an alias of brms's "cumulative". Unknown names pass through.
maihda_normalize_family_name <- function(name) {
  if (is.null(name) || !is.character(name) || length(name) != 1 || is.na(name)) {
    return(name)
  }
  if (grepl("^Negative Binomial\\(", name)) {
    return("negbinomial")
  }
  if (name == "ordinal") {
    return("cumulative")
  }
  name
}

maihda_family <- function(model) {
  fam <- tryCatch(stats::family(model), error = function(e) NULL)
  if (is.null(fam) && inherits(model, "brmsfit")) {
    fam <- tryCatch(model$family, error = function(e) NULL)
  }
  # Canonicalise the family name (see maihda_normalize_family_name) so callers
  # can compare against fixed names; the link and linkinv are left untouched.
  if (!is.null(fam) && !is.null(fam$family)) {
    fam$family <- maihda_normalize_family_name(fam$family)
  }
  fam
}

# TRUE when the model's negative-binomial dispersion theta was FIXED by the user
# (a MASS::negative.binomial(theta) family object passed to glmer), rather than
# ESTIMATED (lme4::glmer.nb or the brms 'shape' parameter). The wrapper records
# the family exactly as it resolved it at fit time: an estimated fit carries the
# canonical "negbinomial" marker, while a fixed-theta fit carries MASS's
# "Negative Binomial(<theta>)" label with theta embedded in the string. The
# fitted object's family reads "Negative Binomial(<theta>)" in BOTH cases, so the
# wrapper-recorded family is the only reliable signal of which one it was.
maihda_negbin_theta_is_fixed <- function(model) {
  fam <- model$family
  if (!is.list(fam) || is.null(fam$family) || !is.character(fam$family) ||
      length(fam$family) != 1 || is.na(fam$family)) {
    return(FALSE)
  }
  grepl("^Negative Binomial\\(", fam$family)
}

# "family(link)" key for a maihda_model, used to decide whether two models are
# comparable (same family and link). Prefers the fitted object's family and
# falls back to the family the wrapper recorded at fit time -- stats::family()
# is undefined for engines like wemix (WeMixResults). Names are canonical via
# maihda_family()/maihda_normalize_family_name(), so e.g. two glmer.nb() fits
# with different ESTIMATED thetas still compare equal.
#
# The negative binomial needs special care: a FIXED, user-specified theta
# (MASS::negative.binomial(theta)) is part of the model SPECIFICATION -- two such
# fits with different thetas assume different dispersion and are NOT comparable --
# so theta is kept in the key. An estimated theta differs between fits only by
# estimation noise and is normalized away (to plain "negbinomial").
maihda_model_family_key <- function(model) {
  fam <- maihda_family(model$model)
  if (is.null(fam) && is.list(model$family)) {
    fam <- model$family
    if (!is.null(fam$family)) {
      fam$family <- maihda_normalize_family_name(fam$family)
    }
  }
  fam_name <- if (!is.null(fam$family)) fam$family else NA_character_
  link <- if (!is.null(fam$link)) fam$link else NA_character_
  if (identical(fam_name, "negbinomial") && maihda_negbin_theta_is_fixed(model)) {
    # Restore the user's fixed-theta label ("Negative Binomial(<theta>)") so two
    # different fixed thetas yield different keys and fail the comparability check.
    fam_name <- model$family$family
  }
  paste0(fam_name, "(", link, ")")
}

maihda_linkinv <- function(fam) {
  if (!is.null(fam) && !is.null(fam$linkinv)) {
    return(fam$linkinv)
  }

  link <- if (!is.null(fam) && !is.null(fam$link)) fam$link else "identity"
  switch(link,
         identity = function(eta) eta,
         log = exp,
         logit = stats::plogis,
         probit = stats::pnorm,
         cloglog = function(eta) 1 - exp(-exp(eta)),
         inverse = function(eta) 1 / eta,
         stop("Unsupported link function for response-scale transformation: ", link, call. = FALSE))
}

# Capture fit-quality diagnostics (singular fit, non-convergence) so they can be
# surfaced on demand. lme4 returns singular fits silently and only warns about
# convergence once, at fit time, so we re-read both from the fitted object for
# reporting in print()/summary(). brms (Stan) convergence is diagnosed elsewhere
# via Rhat, so for brmsfit objects only the engine is recorded here.
maihda_fit_diagnostics <- function(model) {
  diagnostics <- list(
    engine = NA_character_,
    singular = NA,
    converged = NA,
    messages = character(0)
  )

  if (inherits(model, "merMod")) {
    diagnostics$engine <- "lme4"
    diagnostics$singular <- tryCatch(isTRUE(lme4::isSingular(model)),
                                     error = function(e) NA)
    msgs <- tryCatch(model@optinfo$conv$lme4$messages,
                     error = function(e) NULL)
    msgs <- as.character(msgs)
    msgs <- msgs[nzchar(msgs)]
    diagnostics$messages <- msgs
    diagnostics$converged <- length(msgs) == 0
  } else if (inherits(model, "WeMixResults")) {
    diagnostics$engine <- "wemix"
    # WeMix raises a hard error when optimisation fails, so a returned fit has
    # converged; flag a boundary (zero between-stratum variance) fit the same way
    # lme4's isSingular() does, since it makes the VPC unreliable in exactly the
    # same way.
    diagnostics$converged <- TRUE
    diagnostics$singular <- tryCatch({
      vd <- model$varDF
      s <- as.numeric(vd$vcov[vd$grp == "stratum"][1])
      is.finite(s) && s < 1e-8
    }, error = function(e) NA)
  } else if (inherits(model, "clmm")) {
    diagnostics$engine <- "ordinal"
    # clmm stores the optimiser result: convergence code 0 means converged, and
    # the message carries the optimiser's own wording otherwise. A boundary
    # (zero between-stratum variance) fit is flagged like lme4's isSingular().
    conv <- tryCatch(model$optRes$convergence, error = function(e) NULL)
    diagnostics$converged <- if (is.numeric(conv) && length(conv) == 1) {
      conv == 0
    } else {
      NA
    }
    if (isFALSE(diagnostics$converged)) {
      msg <- tryCatch(as.character(model$optRes$message), error = function(e) character(0))
      diagnostics$messages <- msg[nzchar(msg)]
    }
    diagnostics$singular <- tryCatch({
      vc <- ordinal::VarCorr(model)
      s <- as.numeric(vc[["stratum"]][1, 1])
      is.finite(s) && s < 1e-8
    }, error = function(e) NA)
  } else if (inherits(model, "brmsfit")) {
    diagnostics$engine <- "brms"
    # Stan/brms convergence is not flagged at fit time the way lme4 surfaces
    # singular/convergence warnings, so check the standard HMC diagnostics here:
    # the maximum Rhat (chain mixing; > 1.01 is the common threshold) and the
    # number of divergent transitions (a divergence signals the sampler could not
    # explore the posterior, so estimates may be biased). Both are wrapped so a
    # missing brms or an unusual fit degrades to "no diagnostics" rather than error.
    if (requireNamespace("brms", quietly = TRUE)) {
      msgs <- character(0)
      max_rhat <- tryCatch(suppressWarnings(max(brms::rhat(model), na.rm = TRUE)),
                           error = function(e) NA_real_)
      if (is.finite(max_rhat) && max_rhat > 1.01) {
        msgs <- c(msgs, sprintf(
          "Max Rhat = %.3f (> 1.01): the chains may not have converged.", max_rhat))
      }
      n_div <- tryCatch({
        np <- brms::nuts_params(model)
        sum(np$Value[np$Parameter == "divergent__"], na.rm = TRUE)
      }, error = function(e) NA_real_)
      if (is.finite(n_div) && n_div > 0) {
        msgs <- c(msgs, sprintf(
          "%d divergent transition(s): the posterior may be biased; consider increasing adapt_delta.",
          as.integer(n_div)))
      }
      diagnostics$messages <- msgs
      # "No warning" only means "converged" when a diagnostic was actually
      # available to raise one. A one-chain fit, a variational/approximate
      # algorithm, or an unusual object can leave BOTH the maximum Rhat and the
      # divergence count non-finite; there is then no evidence either way, so
      # leave converged = NA rather than silently reporting a clean fit.
      any_diag <- is.finite(max_rhat) || is.finite(n_div)
      diagnostics$converged <- if (any_diag) length(msgs) == 0 else NA
    }
  }

  structure(diagnostics, class = "maihda_fit_diagnostics")
}

# Format fit diagnostics as a character vector of report lines (empty when the fit
# looks clean), shared by the maihda_model and maihda_summary print methods.
maihda_format_fit_diagnostics <- function(diagnostics) {
  if (is.null(diagnostics)) {
    return(character(0))
  }

  out <- character(0)
  if (isTRUE(diagnostics$singular)) {
    out <- c(
      out,
      "Singular fit: at least one variance component is estimated at (or near) zero.",
      "  The between-stratum variance and any VPC/PCV derived from it may be unreliable."
    )
  }
  if (isFALSE(diagnostics$converged) && length(diagnostics$messages) > 0) {
    header <- switch(
      as.character(diagnostics$engine),
      brms = "MCMC convergence diagnostics (brms):",
      ordinal = "Convergence warnings reported by ordinal::clmm:",
      lme4 = "Convergence warnings reported by lme4:",
      "Convergence warnings:"
    )
    out <- c(out, header, paste0("  - ", diagnostics$messages))
  }
  out
}

# Print the fit-diagnostics block (nothing is printed for a clean fit).
maihda_print_fit_diagnostics <- function(diagnostics) {
  diag_lines <- maihda_format_fit_diagnostics(diagnostics)
  if (length(diag_lines) == 0) {
    return(invisible(NULL))
  }
  pal <- maihda_palette()
  cat(pal$warn("Fit diagnostics:"), "\n", sep = "")        # genuine caveat -> warn
  cat(pal$warn(paste0("  ", diag_lines)), sep = "\n")
  cat("\n\n")
  invisible(NULL)
}

# Run `expr`, returning its value; on error return NULL but re-emit the failure
# as a warning carrying the ORIGINAL condition's message, so an output a caller
# asked for is never dropped without a trace. `what` names the output for the
# message. Used for the optional summary / plot companions that must not break
# the core result but also must not fail silently (a bare tryCatch(., NULL)
# hides the reason a requested calculation is missing).
maihda_try_optional <- function(expr, what) {
  tryCatch(
    expr,
    error = function(e) {
      warning(sprintf("%s could not be computed and was omitted: %s",
                      what, conditionMessage(e)), call. = FALSE)
      NULL
    }
  )
}

maihda_validate_conf_level <- function(conf_level) {
  if (!is.numeric(conf_level) || length(conf_level) != 1 ||
      is.na(conf_level) || !is.finite(conf_level) ||
      conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be a single number between 0 and 1.", call. = FALSE)
  }

  as.numeric(conf_level)
}

# Reduce successful bootstrap draws to a central interval, requiring a minimum
# number of successful refits (so an interval is never returned from one or a
# handful of draws) and warning when the failure rate is high.
maihda_bootstrap_ci <- function(values, n_boot, conf_level, what = "VPC") {
  values <- values[is.finite(values)]
  n_ok <- length(values)
  failed <- n_boot - n_ok
  min_ok <- 10L

  if (n_ok == 0) {
    stop("All ", what, " bootstrap refits failed; no interval can be computed.",
         call. = FALSE)
  }
  if (n_ok < min_ok) {
    stop(sprintf(paste0("Only %d of %d %s bootstrap refits succeeded; at least %d ",
                        "are required to form an interval. Increase n_boot or check ",
                        "for singular/failing fits."),
                 n_ok, n_boot, what, min_ok), call. = FALSE)
  }
  if (failed > n_boot * 0.5) {
    warning(sprintf("%d of %d %s bootstrap refits failed (%.0f%%); the interval may be unreliable.",
                    failed, n_boot, what, 100 * failed / n_boot), call. = FALSE)
  }

  alpha <- 1 - conf_level
  ci <- stats::quantile(values, probs = c(alpha / 2, 1 - alpha / 2), names = FALSE)
  # n_ok is the effective number of successful refits; mc_se is the Monte Carlo
  # standard error OF THE BOOTSTRAP MEAN (sd / sqrt(n_ok)) -- a coarse gauge of
  # bootstrap noise, NOT the uncertainty of the percentile interval's endpoints,
  # which is a different (larger) order-statistic quantity. The print methods
  # label it accordingly.
  attr(ci, "n_ok") <- n_ok
  attr(ci, "mc_se") <- if (n_ok > 1) stats::sd(values) / sqrt(n_ok) else NA_real_
  ci
}

maihda_validate_bootstrap_args <- function(n_boot, conf_level) {
  # Forming an interval needs at least maihda_bootstrap_ci()'s minimum number of
  # successful refits, so reject n_boot below that floor here -- failing fast with
  # a clear message rather than only erroring later, after the bootstrap has run.
  min_boot <- 10L
  if (!is.numeric(n_boot) || length(n_boot) != 1 ||
      is.na(n_boot) || !is.finite(n_boot) ||
      n_boot < min_boot || n_boot != floor(n_boot)) {
    stop("'n_boot' must be a single whole number >= ", min_boot,
         " (at least ", min_boot, " successful refits are required to form an interval).",
         call. = FALSE)
  }

  conf_level <- maihda_validate_conf_level(conf_level)

  # The hard floor above (10) exists so the test suite can bootstrap quickly; the
  # default is 1000. Between the two, a percentile interval's tail endpoints are
  # order statistics estimated from very few draws, so they are unstable -- warn
  # (without blocking) when n_boot is well below the count needed for a dependable
  # 2.5/97.5% interval. Published guidance uses hundreds to ~1000+ replications.
  rec_boot <- 200L
  if (n_boot < rec_boot) {
    warning(sprintf(paste0("n_boot = %d is low for a %g%% percentile interval: the ",
                           "tail endpoints are order statistics from few draws and are ",
                           "unstable below ~%d replications. Increase n_boot (the ",
                           "default, 1000, is dependable) or read the interval as ",
                           "indicative only."),
                    as.integer(n_boot), 100 * conf_level, rec_boot),
            call. = FALSE)
  }

  list(n_boot = as.integer(n_boot), conf_level = conf_level)
}

maihda_quote_name <- function(name) {
  if (!is.character(name) || length(name) != 1 || is.na(name) || name == "") {
    stop("Variable names must be non-empty character strings.", call. = FALSE)
  }

  paste(deparse(as.name(name), backtick = TRUE), collapse = "")
}

# Drop random-effect bar terms from a formula, keeping the response intact.
# reformulas::nobars() descends into the LEFT-hand side too, so a brms
# addition-term response such as `y | trials(n)` is itself treated as a bar
# term and stripped -- `y | trials(n) ~ x + (1 | g)` becomes `x ~ 1`, silently
# promoting a predictor to the response. Stripping only the right-hand side
# preserves the response expression exactly (also covering call-valued
# responses like cbind(successes, failures), for which nobars() on a bars-only
# RHS returns a bare call rather than a formula).
maihda_nobars <- function(formula) {
  if (!inherits(formula, "formula")) {
    return(reformulas::nobars(formula))
  }
  out <- formula
  rhs <- reformulas::nobars(formula[[length(formula)]])
  if (is.null(rhs)) rhs <- 1
  out[[length(out)]] <- rhs
  out
}

# Detect fixed-effect interaction term(s) among the stratum dimensions in a formula's
# fixed part -- the bug behind a corrupt MAIHDA decomposition. The shorthand
# `var1 * var2` expands to `var1 + var2 + var1:var2`, so the dimensions' main-effect
# detection (which only looks for the additive terms) classifies the formula as
# fully-specified and strips ONLY the main effects, leaving the fixed `var1:var2`
# interaction in both the derived null and adjusted formulas. That fixed cell-means
# term duplicates the intersectional stratum random intercept (1 | var1:var2): it
# absorbs the between-stratum variance into fixed effects, pinning the stratum RE at a
# singular boundary and collapsing the PCV (and the pure-interaction BLUP diagnostic).
#
# Returns the offending fixed-effect term labels (character(0) when there are none).
# Detection uses the terms() factors matrix (rows = variables, columns = terms) rather
# than string-splitting the ":"-joined labels, so it is robust to non-syntactic names
# (e.g. `gender var`) and to variable order within the interaction. A term is flagged
# when its interaction order is >= 2 and EVERY variable it involves is a stratum
# dimension; a dimension-by-covariate interaction (e.g. age:gender, a legitimate
# covariate adjustment) is left alone. `strata_vars` are the raw dimension names and
# `dim_terms` (optional) the adjusted-model main-effect terms -- the reconstructed
# `.maihda_dim_*` factor for an auto-binned numeric dimension -- both matched in raw
# and backtick-quoted form (terms() quotes non-syntactic variable labels).
maihda_dimension_interaction_terms <- function(formula, strata_vars,
                                               dim_terms = character(0)) {
  if (is.null(strata_vars) || length(strata_vars) < 2) {
    return(character(0))
  }
  tt <- tryCatch(stats::terms(maihda_nobars(formula)),
                 error = function(e) NULL)
  if (is.null(tt)) return(character(0))
  factors <- attr(tt, "factors")
  if (is.null(factors) || length(factors) == 0) return(character(0))
  ord <- attr(tt, "order")
  term_labels <- colnames(factors)
  var_labels <- rownames(factors)

  dim_names <- unique(c(strata_vars, dim_terms))
  dim_quoted <- vapply(dim_names, maihda_quote_name, character(1))
  is_dim_var <- var_labels %in% dim_names | var_labels %in% dim_quoted

  flagged <- character(0)
  for (j in seq_along(term_labels)) {
    if (ord[j] < 2L) next
    involved <- factors[, j] != 0
    if (any(involved) && all(is_dim_var[involved])) {
      flagged <- c(flagged, term_labels[j])
    }
  }
  flagged
}

# Stop with a single, actionable message when a formula carries a fixed interaction
# among the stratum dimensions (see maihda_dimension_interaction_terms()). Called by
# the workflow entry points (maihda(), compare_maihda_groups()) before they derive
# the null/adjusted formulas, so the user is rejected up front rather than handed a
# silently corrupt (NA/NULL) PCV. Rejection -- not silent stripping -- is the safe
# choice: the MAIHDA adjusted model is defined to carry only the dimensions' ADDITIVE
# main effects, with the intersection estimated by the stratum random effect.
maihda_check_no_dimension_interaction <- function(formula, strata_vars,
                                                  dim_terms = character(0),
                                                  fn = "maihda") {
  flagged <- maihda_dimension_interaction_terms(formula, strata_vars, dim_terms)
  if (length(flagged) == 0) {
    return(invisible(NULL))
  }
  stop(fn, "(): the formula's fixed part contains the interaction term(s) ",
       paste(flagged, collapse = ", "), " between the stratum-defining dimensions (",
       paste(strata_vars, collapse = ", "), "). A fixed interaction among the ",
       "stratum dimensions duplicates the intersectional stratum random intercept ",
       "(1 | ", paste(strata_vars, collapse = ":"), "): it absorbs the ",
       "between-stratum variance into fixed cell means, which pins the stratum ",
       "variance at a singular boundary and makes the PCV invalid. The MAIHDA ",
       "adjusted model takes only the dimensions' ADDITIVE main effects -- write ",
       "`", paste(strata_vars, collapse = " + "), "` (a `*` expands to the additive ",
       "terms PLUS this interaction). The intersection itself is estimated by the ",
       "stratum random effect; use decomposition = \"crossed-dimensions\" or ",
       "maihda_interactions() to quantify it.", call. = FALSE)
}

maihda_formula_with_stratum <- function(outcome, vars = character()) {
  if (!is.character(vars)) {
    stop("'vars' must be a character vector.", call. = FALSE)
  }

  fixed_terms <- vapply(vars, maihda_quote_name, character(1))
  random_term <- paste0("(1 | ", maihda_quote_name("stratum"), ")")
  rhs <- c(if (length(fixed_terms) > 0) fixed_terms else "1", random_term)

  stats::reformulate(rhs, response = maihda_quote_name(outcome))
}

maihda_is_binary_vector <- function(x) {
  if (!is.null(dim(x))) {
    return(FALSE)
  }

  x <- x[!is.na(x)]
  length(unique(x)) == 2
}

maihda_binary_levels <- function(x) {
  x <- x[!is.na(x)]

  if (is.factor(x)) {
    return(levels(droplevels(x)))
  }
  if (is.logical(x)) {
    return(c(FALSE, TRUE))
  }
  if (is.numeric(x)) {
    return(sort(unique(x)))
  }

  levels(factor(x))
}

maihda_binary_to_01 <- function(x) {
  levels_x <- maihda_binary_levels(x)
  key <- as.character(levels_x)
  out <- match(as.character(x), key) - 1L
  out[is.na(x)] <- NA_integer_
  as.integer(out)
}

# All levels present in a stratum dimension, in display order (factor order,
# FALSE/TRUE for logical, sorted unique for numeric). Generalises
# maihda_binary_levels() to multi-level factors so the upset matrix can lay out
# one row per level.
maihda_dim_levels <- function(x) {
  x <- x[!is.na(x)]
  if (is.factor(x)) {
    return(levels(droplevels(x)))
  }
  if (is.logical(x)) {
    return(c(FALSE, TRUE))
  }
  if (is.numeric(x)) {
    return(sort(unique(x)))
  }
  levels(factor(x))
}

# TRUE when a dimension is a clean present/absent indicator -- logical, or numeric
# coded exactly {0, 1} -- so the upset matrix can collapse it to a single
# present/absent row (the classic UpSet look). Factors (including two-level ones,
# where neither level is an "absence") and other codings get one row per level.
maihda_dim_is_indicator <- function(x) {
  if (is.logical(x)) {
    return(length(maihda_dim_levels(x)) == 2)
  }
  if (is.numeric(x)) {
    lv <- maihda_dim_levels(x)
    return(length(lv) == 2 && all(as.character(lv) %in% c("0", "1")))
  }
  FALSE
}

# Sampling/design weights are dropped by the weighted engines (WeMix, the brms
# pseudo-likelihood) whenever they are non-finite or non-positive -- the
# is.finite(w) & w > 0 rule in maihda_fit_wemix() / maihda_prepare_brms_sampling_weights().
# lme4 PRIOR weights, by contrast, only drop NA. The shared analytic-frame row mask
# (maihda_row_mask) implements the NA-only rule, so to make binary/ordinal detection
# and 0/1 recoding see exactly the rows a weighted engine fits, map the engine-excluded
# weights (non-finite or <= 0) to NA here before handing the vector to the detection
# chain. Without this a zero-weight row carrying an out-of-sample outcome value could
# flip detection (e.g. fit gaussian instead of binomial) relative to the fitted rows.
maihda_sampling_weight_mask <- function(w) {
  w <- as.numeric(w)
  w[!is.finite(w) | w <= 0] <- NA_real_
  w
}

# lme4 PRECISION weights (weights=) that are zero -- or negative / non-finite --
# carry no information, but lmer does NOT drop such a row: it returns a degenerate
# fit (logLik -Inf, NA gradient) while the residual-variance helper silently
# discards the zero weight and still reports a finite VPC, so the variance
# estimates differ materially from fitting after removing the row. Map those
# weights to NA -- the case lme4 DOES drop, reproducing the row-removed fit
# exactly -- so binary/ordinal (family) detection, strata auto-binning, the
# analytic-size checks (min_group_n), and the fit all key off the SAME analytic
# sample. Only the lme4 path takes precision weights; wemix/brms/ordinal reject
# them, so this is a no-op there (and for a NULL / non-numeric / wrong-length
# value). Returns the (possibly remapped) weights and warns once -- tagged with
# `fn` -- when any row is remapped. Shared by fit_maihda(), compare_maihda_groups()
# and maihda() so the three cannot drift apart: the high-level APIs had each
# re-derived the analytic sample WITHOUT this normalization, letting an invalid
# precision weight flip the detected family/engine and inflate a group's analytic n.
maihda_normalize_precision_weights <- function(weights_value, n, engine, fn) {
  if (identical(engine, "wemix") || identical(engine, "brms") ||
      !is.numeric(weights_value) || length(weights_value) != n) {
    return(weights_value)
  }
  bad_w <- !is.na(weights_value) & !(is.finite(weights_value) & weights_value > 0)
  if (any(bad_w)) {
    warning(sprintf(paste0("%s: dropped %d row(s) with a non-positive or ",
                           "non-finite precision weight ('weights') before ",
                           "fitting (a zero weight otherwise degenerates the ",
                           "lme4 fit)."), fn, sum(bad_w)), call. = FALSE)
    weights_value[bad_w] <- NA_real_
  }
  weights_value
}

# Logical keep-mask over `n` rows reproducing the row selection lme4/brms apply
# before fitting: an (already-evaluated) `subset` value -- logical (recycled,
# NA -> drop), positive/negative numeric indices, or character row names -- the
# removal of rows whose (already-evaluated) prior `weights` is NA, and the removal
# of rows whose (already-evaluated) external `offset` is NA. lme4 puts the offset in
# its model frame, so na.omit drops offset-NA rows exactly as it drops weight-NA and
# response/covariate-NA rows; folding the offset in here keeps binary/ordinal family
# detection, 0/1 recoding, strata auto-binning and longitudinal centering keyed off
# the SAME analytic sample the engine fits. Without it an offset-NA row carrying an
# out-of-sample outcome value flips the family (e.g. gaussian instead of binomial),
# or an out-of-range occasion mis-anchors the growth time-centering. Subset, weights
# and offset arrive as VALUES (resolved with the data mask in fit_maihda) so this
# helper never needs to evaluate user expressions.
maihda_row_mask <- function(data, subset = NULL, weights = NULL, offset = NULL) {
  n <- nrow(data)
  keep <- rep(TRUE, n)
  if (!is.null(subset)) {
    if (is.logical(subset)) {
      s <- if (length(subset) == n) subset else rep_len(subset, n)
      keep <- keep & !is.na(s) & s
    } else if (is.numeric(subset)) {
      sel <- logical(n)
      idx <- tryCatch(seq_len(n)[subset], error = function(e) integer(0))
      idx <- idx[!is.na(idx) & idx >= 1L & idx <= n]
      sel[idx] <- TRUE
      keep <- keep & sel
    } else if (is.character(subset)) {
      keep <- keep & (rownames(data) %in% subset)
    }
  }
  if (!is.null(weights) && length(weights) == n) {
    keep <- keep & !is.na(weights)
  }
  if (!is.null(offset) && length(offset) == n) {
    keep <- keep & !is.na(offset)
  }
  keep
}

# Expand an (already-evaluated) `subset` value into a full-length logical mask over
# `n` rows. Positional subsets -- numeric row indices (e.g. c(1:10, 31:40)) or a
# recycled logical mask -- are GLOBAL: index k means the k-th row of the data the
# subset was written against. They must be turned into a full-length mask BEFORE
# being sliced per group, otherwise the same vector is silently reinterpreted
# relative to each subgroup and selects the wrong rows. Mirrors the subset
# handling in maihda_row_mask(). Character (row-name) subsets are name-based, not
# positional -- base `[` preserves row names in a group's slice, so they are
# returned unchanged and matched per group by name; NULL stays NULL.
maihda_normalize_subset <- function(subset, n) {
  if (is.null(subset)) {
    return(NULL)
  }
  if (is.logical(subset)) {
    mask <- if (length(subset) == n) subset else rep_len(subset, n)
    mask[is.na(mask)] <- FALSE
    return(mask)
  }
  if (is.numeric(subset)) {
    mask <- logical(n)
    # Let base R validate the numeric subscripts rather than swallowing its error.
    # `seq_len(n)[subset]` only errors on a genuinely invalid subset -- mixing
    # positive and negative indices (e.g. c(-1, 2)), or mixing NA with negative
    # indices; out-of-range indices merely yield NA, which is filtered below. The old
    # tryCatch(..., error = integer(0)) collapsed those invalid subsets to an empty
    # mask, so compare_maihda_groups()/maihda() silently fit zero rows (every group
    # reported as skipped, n = 0) instead of erroring like base R and a bare fit do.
    idx <- tryCatch(
      seq_len(n)[subset],
      error = function(e) stop("invalid 'subset': ", conditionMessage(e),
                               call. = FALSE))
    idx <- idx[!is.na(idx) & idx >= 1L & idx <= n]
    mask[idx] <- TRUE
    return(mask)
  }
  subset
}

# The analytic model frame lme4/brms will actually fit: response and fixed-effect
# transformations are applied, the grouping variables are retained, any `subset` is
# honoured, rows with a missing prior weight are dropped, and rows missing after
# those transformations are dropped (na.omit). This is what the binary
# auto-detection and the small-sample guards key off, so they see exactly the rows
# the model uses rather than the raw columns -- which ignore transformations (e.g.
# log(x) of a non-positive value), subsetting, and weight-based row removal.
# Returns NULL when the frame cannot be built (callers fall back to a raw check).
maihda_analytic_model_frame <- function(formula, data, subset = NULL,
                                        weights = NULL, offset = NULL) {
  keep <- maihda_row_mask(data, subset = subset, weights = weights, offset = offset)
  data <- data[keep, , drop = FALSE]

  fr_form <- tryCatch(reformulas::subbars(formula), error = function(e) NULL)
  if (is.null(fr_form)) {
    return(NULL)
  }
  environment(fr_form) <- environment(formula)
  tryCatch(
    stats::model.frame(fr_form, data = data, na.action = stats::na.omit),
    error = function(e) NULL
  )
}

# Logical keep-mask (over the rows of `data`) selecting exactly the rows that survive
# into the analytic model frame above: rows kept by the subset / prior-weight row mask
# AND not dropped by na.omit once the formula's response and fixed-effect
# transformations are applied. This is maihda_analytic_model_frame()'s row selection
# exposed as a position-based mask (the na.action indices are mapped back to original
# rows, so it does not rely on rownames being preserved), letting the wemix/ordinal
# engines prefilter `data` to the rows their fit actually uses. A raw
# complete.cases() over the formula's columns misses NAs that only appear after a
# transformed term (e.g. log(x) of x <= 0), so the stored frame and the engine's
# fitted rows could disagree. Returns NULL when the frame cannot be built (callers
# fall back to a raw complete.cases()).
maihda_analytic_keep_mask <- function(formula, data, subset = NULL, weights = NULL,
                                      offset = NULL) {
  premask <- maihda_row_mask(data, subset = subset, weights = weights, offset = offset)
  fr_form <- tryCatch(reformulas::subbars(formula), error = function(e) NULL)
  if (is.null(fr_form)) {
    return(NULL)
  }
  environment(fr_form) <- environment(formula)
  mf <- tryCatch(
    stats::model.frame(fr_form, data = data[premask, , drop = FALSE],
                       na.action = stats::na.omit),
    error = function(e) NULL
  )
  if (is.null(mf)) {
    return(NULL)
  }
  keep <- premask
  omit <- attr(mf, "na.action")
  if (!is.null(omit)) {
    # `omit` indexes rows WITHIN data[premask, ]; map those positions back to the
    # corresponding original rows before clearing them.
    kept_idx <- which(premask)
    keep[kept_idx[as.integer(omit)]] <- FALSE
  }
  keep
}

# The model response over the analytic sample (post-transformation, post-NA,
# post-subset, post-weight-NA and post-offset-NA). Only plain-symbol responses
# qualify as a Bernoulli candidate, so a transformed or aggregated response
# (log(y), cbind(s, f), `y | trials(n)`) yields NULL -- "not a single two-level
# response".
maihda_analytic_response <- function(formula, data, subset = NULL,
                                     weights = NULL, offset = NULL) {
  if (length(formula) != 3L || !is.symbol(formula[[2]])) {
    return(NULL)
  }
  fr <- maihda_analytic_model_frame(formula, data, subset = subset,
                                    weights = weights, offset = offset)
  if (is.null(fr)) {
    return(NULL)
  }
  tryCatch(stats::model.response(fr), error = function(e) NULL)
}

# Recode a vector to 0/1 using exactly two reference levels: the first level
# becomes 0, the second becomes 1, and anything else (including NA) becomes NA.
maihda_recode_to_01 <- function(x, levels_2) {
  key <- as.character(levels_2)
  as.integer(match(as.character(x), key) - 1L)
}

maihda_prepare_binomial_response <- function(data, formula, subset = NULL,
                                             weights = NULL, offset = NULL) {
  response <- formula[[2]]
  if (!is.symbol(response)) {
    return(data)
  }

  outcome <- as.character(response)
  # Recode against the analytic sample lme4/glmer actually fits (transformations
  # applied; rows excluded by `subset`, a missing weight, a missing offset, or
  # missingness dropped), matching the binary detection in fit_maihda(). A
  # character/factor outcome whose
  # third value appears only on excluded rows is still recoded to 0/1; the
  # out-of-sample value becomes NA (and is dropped) rather than left as a stray
  # level that glmer() would reject.
  resp <- maihda_analytic_response(formula, data, subset = subset,
                                   weights = weights, offset = offset)
  if (!outcome %in% names(data) || is.null(resp) ||
      !maihda_is_binary_vector(resp)) {
    return(data)
  }

  analytic_levels <- maihda_binary_levels(resp)
  data[[outcome]] <- maihda_recode_to_01(data[[outcome]], analytic_levels)

  # Record (and surface) which level became the modeled event (= 1). The mapping
  # follows base glm/glmer convention -- for a character outcome the levels are
  # alphabetical, for a factor the declared level order -- so "case"/"control"
  # could map either way depending on the column type. That is easy to invert
  # silently, especially when family = "binomial" is passed explicitly (no
  # auto-detect warning fires), so we attach the mapping to the data (fit_maihda()
  # stores it on the model as $response_recoding) and emit one informational
  # message. A response already coded 0/1 is a no-op and stays silent.
  recoding <- data.frame(
    level = as.character(analytic_levels),
    value = c(0L, 1L),
    role = c("reference", "event"),
    stringsAsFactors = FALSE
  )
  attr(data, "response_recoding") <- recoding
  if (!identical(as.character(analytic_levels), c("0", "1"))) {
    message(sprintf(
      "Binary outcome '%s' recoded to 0/1: '%s' = 0 (reference), '%s' = 1 (modeled event). %s",
      outcome, analytic_levels[1], analytic_levels[2],
      "Set the factor levels (or supply a 0/1 outcome) to control which level is the event."
    ))
  }
  data
}

# TRUE only when the model response is a single two-level (Bernoulli) vector,
# i.e. a plain symbol naming a binary column. Aggregated binomial responses such
# as cbind(success, failure) or `y | trials(n)` are calls, not symbols, so they
# return FALSE and must remain a binomial() model.
#
# The check is evaluated on the analytic sample lme4/brms actually fits -- the
# model frame after transformations, missing-row dropping, any `subset`, and the
# removal of rows with a missing prior weight -- so a response that is 0/1 only
# once excluded rows are removed is still recognised as Bernoulli.
maihda_response_is_binary <- function(formula, data, subset = NULL,
                                      weights = NULL, offset = NULL) {
  resp <- maihda_analytic_response(formula, data, subset = subset,
                                   weights = weights, offset = offset)
  if (is.null(resp)) {
    return(FALSE)
  }
  maihda_is_binary_vector(resp)
}

# TRUE when the analytic response is an ORDERED factor with 3+ observed levels
# -- the signature of a cumulative (ordinal) outcome. Used by fit_maihda()'s
# default-family auto-detection, after the binary check (a 2-level ordered
# factor is a binomial model). An unordered factor stays FALSE: its level order
# is not declared meaningful, so silently treating it as ordinal would be wrong.
maihda_response_is_ordinal <- function(formula, data, subset = NULL,
                                       weights = NULL, offset = NULL) {
  resp <- maihda_analytic_response(formula, data, subset = subset,
                                   weights = weights, offset = offset)
  if (is.null(resp)) {
    return(FALSE)
  }
  is.ordered(resp) && nlevels(droplevels(resp)) >= 3
}

# TRUE if a random-effect grouping expression is a plain variable or a colon
# interaction of plain variables (e.g. a, a:b, a:b:c) -- the only forms automatic
# strata creation understands. Function-call terms such as interaction(a, b),
# paste(a, b) or cut(age, 3) are FALSE: their semantics would be silently ignored
# (all.vars() would just extract a and b), so they must go through make_strata().
maihda_is_colon_interaction <- function(expr) {
  if (is.symbol(expr)) {
    return(TRUE)
  }
  if (is.call(expr) && identical(expr[[1]], as.name(":"))) {
    return(all(vapply(as.list(expr)[-1], maihda_is_colon_interaction, logical(1))))
  }
  FALSE
}

maihda_model_frame <- function(model, fallback = NULL) {
  out <- tryCatch(stats::model.frame(model), error = function(e) NULL)
  if (is.null(out) && inherits(model, "merMod")) {
    out <- tryCatch(model@frame, error = function(e) NULL)
  }
  if (is.null(out)) {
    out <- fallback
  }
  out
}

maihda_nobs <- function(model) {
  tryCatch(stats::nobs(model), error = function(e) {
    frame <- maihda_model_frame(model)
    if (is.null(frame)) NA_integer_ else nrow(frame)
  })
}

# Content fingerprint of a model's analytic response vector. Two models fitted to
# the same data share a fingerprint even if the rows were reordered or carry
# default 1:n names; models fitted to unrelated data do not. Used to catch
# comparisons/PCV across different datasets that happen to share n, row names and
# stratum ids.
maihda_response_fingerprint <- function(model) {
  frame <- maihda_model_frame(model)
  if (is.null(frame)) {
    return(NA_character_)
  }
  resp <- tryCatch(stats::model.response(frame), error = function(e) NULL)
  if (is.null(resp)) {
    return(NA_character_)
  }
  resp <- unname(resp)
  if (is.numeric(resp)) {
    paste(formatC(resp, format = "g", digits = 12), collapse = "\r")
  } else {
    paste(as.character(resp), collapse = "\r")
  }
}

# Fingerprint of a model's prior weights, so PCV / model comparisons do not
# silently combine fits that share an outcome, sample and strata but used DIFFERENT
# prior weights (which change the variance estimates). An unweighted fit and an
# explicit weights = rep(1, n) fit both map to "unit", so they compare as equal;
# engines where prior weights are not recoverable (e.g. brms) also degrade to
# "unit" rather than erroring, leaving their current behaviour unchanged.
maihda_weight_fingerprint <- function(model) {
  w <- tryCatch(stats::weights(model, type = "prior"), error = function(e) NULL)
  if (is.null(w) || length(w) == 0) {
    return("unit")
  }
  w <- as.numeric(w)
  if (all(is.finite(w)) && isTRUE(all(abs(w - 1) < sqrt(.Machine$double.eps)))) {
    return("unit")
  }
  paste(formatC(w, format = "g", digits = 12), collapse = "\r")
}

maihda_row_ids <- function(model) {
  frame <- maihda_model_frame(model)
  if (is.null(frame)) {
    return(NULL)
  }
  row.names(frame)
}

# ---- wrapper-level analytic-sample helpers ---------------------------------
# The nobs / row-id / response-fingerprint helpers above read the RAW fitted
# object's model frame, which is undefined for some engines -- notably
# WeMixResults, whose stats::nobs()/model.frame() are not implemented. These
# wrapper-level companions take the maihda_model and fall back to its stored
# analytic $data (the exact rows the engine fit) so the PCV and model-comparison
# sample-identity checks are NOT silently skipped for those engines (two WeMix
# fits with the same formula/n/strata but different outcome values must still be
# distinguished).

maihda_wrapper_nobs <- function(model) {
  n <- maihda_nobs(model$model)
  if (!is.finite(n) && is.data.frame(model$data)) {
    n <- nrow(model$data)
  }
  n
}

maihda_wrapper_row_ids <- function(model) {
  ids <- maihda_row_ids(model$model)
  if (!is.null(ids)) {
    return(ids)
  }
  if (is.data.frame(model$data)) {
    return(row.names(model$data))
  }
  NULL
}

maihda_wrapper_response_fingerprint <- function(model) {
  fp <- maihda_response_fingerprint(model$model)
  if (!is.na(fp)) {
    return(fp)
  }
  # Engine exposes no usable model frame (WeMixResults): fingerprint the response
  # from the wrapper's analytic $data via the (two-sided) formula instead. The
  # encoding matches maihda_response_fingerprint() exactly so a value computed
  # this way compares equal to one read from a model frame on the same data.
  if (is.data.frame(model$data) && !is.null(model$formula) &&
      length(model$formula) == 3L) {
    resp <- tryCatch(eval(model$formula[[2]], envir = model$data),
                     error = function(e) NULL)
    if (!is.null(resp)) {
      resp <- unname(resp)
      return(if (is.numeric(resp)) {
        paste(formatC(resp, format = "g", digits = 12), collapse = "\r")
      } else {
        paste(as.character(resp), collapse = "\r")
      })
    }
  }
  NA_character_
}

maihda_infer_strata_vars <- function(strata_info) {
  if (is.null(strata_info)) {
    return(NULL)
  }

  vars <- setdiff(names(strata_info), c("stratum", "label", "n"))
  if (length(vars) == 0) {
    return(NULL)
  }

  vars
}

maihda_refresh_strata_counts <- function(strata_info, data) {
  if (is.null(strata_info) ||
      !"stratum" %in% names(strata_info) ||
      is.null(data) ||
      !"stratum" %in% names(data)) {
    return(strata_info)
  }

  counts <- table(as.character(data$stratum), useNA = "no")
  refreshed_n <- as.integer(counts[match(as.character(strata_info$stratum), names(counts))])
  refreshed_n[is.na(refreshed_n)] <- 0L
  strata_info$n <- refreshed_n

  strata_info
}

maihda_match_strata_rows <- function(data, lookup, vars) {
  if (length(vars) == 0) {
    return(rep(NA_integer_, nrow(data)))
  }
  if (nrow(data) == 0) {
    return(integer())
  }
  if (is.null(lookup) || nrow(lookup) == 0) {
    return(rep(NA_integer_, nrow(data)))
  }

  # Encode each row as a composite key built from per-column integer codes joined
  # by "\r". The codes come from a shared per-variable level table (the distinct
  # values seen in `lookup`), so a missing value, or a value not present in
  # `lookup`, gets no code and the row matches no stratum. Integer codes cannot
  # contain "\r", so distinct value combinations always produce distinct keys even
  # when the values themselves contain the display separator. This preserves the
  # exact semantics of the previous row-by-row matcher (first matching lookup row,
  # NA when nothing matches) while replacing its O(rows * strata * vars) scan with a
  # single vectorised match(), which scales to large survey data.
  levels_by_var <- lapply(vars, function(v) unique(as.character(lookup[[v]])))
  names(levels_by_var) <- vars

  encode <- function(df) {
    code_cols <- lapply(vars, function(v) {
      match(as.character(df[[v]]), levels_by_var[[v]])
    })
    any_missing <- Reduce(`|`, lapply(code_cols, is.na))
    keys <- do.call(paste, c(code_cols, list(sep = "\r")))
    keys[any_missing] <- NA_character_
    keys
  }

  match(encode(data), encode(lookup))
}

maihda_non_intercept_effects <- function(effect_names) {
  if (is.null(effect_names)) {
    return(character())
  }

  setdiff(effect_names, c("(Intercept)", "Intercept"))
}

maihda_stop_for_random_slopes <- function(offending, context) {
  if (length(offending) == 0) {
    return(invisible(TRUE))
  }

  details <- paste(
    sprintf(
      "%s (%s)",
      names(offending),
      vapply(offending, paste, collapse = ", ", FUN.VALUE = character(1))
    ),
    collapse = "; "
  )

  stop(
    context,
    " currently supports intercept-only random effects. Random slopes were found in: ",
    details,
    ". Fit random-intercept MAIHDA models for VPC/ICC summaries.",
    call. = FALSE
  )
}

maihda_validate_intercept_only_random_effects_lme4 <- function(model, context = "MAIHDA variance calculations") {
  vc <- lme4::VarCorr(model)
  offending <- list()

  for (group in names(vc)) {
    group_mat <- as.matrix(vc[[group]])
    non_intercepts <- maihda_non_intercept_effects(rownames(group_mat))
    if (length(non_intercepts) > 0) {
      offending[[group]] <- non_intercepts
    }
  }

  maihda_stop_for_random_slopes(offending, context)
  invisible(TRUE)
}

maihda_stratum_variance_lme4 <- function(model, group = "stratum") {
  vc <- lme4::VarCorr(model)
  if (!group %in% names(vc)) {
    stop("No '", group, "' random-effect variance found in the model.")
  }

  group_vc <- as.matrix(vc[[group]])
  effect_names <- rownames(group_vc)
  intercept_name <- intersect(c("(Intercept)", "Intercept"), effect_names)
  if (length(intercept_name) == 0) {
    stop("The '", group, "' random effect must include an intercept for MAIHDA variance calculations.")
  }

  as.numeric(group_vc[intercept_name[1], intercept_name[1]])
}

maihda_total_random_variance_lme4 <- function(model) {
  maihda_validate_intercept_only_random_effects_lme4(model)

  vc <- lme4::VarCorr(model)
  variances <- unlist(lapply(vc, function(group_vc) {
    group_mat <- as.matrix(group_vc)
    if (is.null(dim(group_mat))) {
      return(numeric())
    }
    vals <- as.numeric(diag(group_mat))
    vals[is.finite(vals)]
  }), use.names = FALSE)

  sum(variances, na.rm = TRUE)
}

maihda_stratum_variance_brms <- function(model, group = "stratum") {
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("Package 'brms' is required to work with brms models. Please install it with: install.packages('brms')")
  }

  # Posterior-mean between-stratum variance E[sd^2] from the SD draws. Using the
  # draws (rather than the posterior summary SD squared, E[sd]^2) keeps
  # calculate_pcv()/stepwise_pcv() consistent with the draws-based VPC reported by
  # summary.maihda_model(). Falls back to the summary SD if draws are unavailable.
  draws <- tryCatch(maihda_posterior_draws_brms(model), error = function(e) NULL)
  if (!is.null(draws)) {
    rv <- tryCatch(maihda_random_variance_draws_brms(draws, group = group),
                   error = function(e) NULL)
    if (!is.null(rv)) {
      return(mean(rv$stratum))
    }
  }

  vc <- brms::VarCorr(model)
  if (!group %in% names(vc)) {
    stop("No '", group, "' random-effect variance found in the brms model.")
  }

  sd_tab <- vc[[group]]$sd
  if (is.null(dim(sd_tab))) {
    stop("Could not extract '", group, "' standard deviations from the brms model.")
  }

  effect_names <- rownames(sd_tab)
  idx <- match(TRUE, effect_names %in% c("(Intercept)", "Intercept"))
  if (is.na(idx)) {
    if (nrow(sd_tab) == 1) {
      idx <- 1
    } else {
      stop("The '", group, "' random effect must include an intercept for MAIHDA variance calculations.")
    }
  }

  as.numeric(sd_tab[idx, "Estimate"]^2)
}

maihda_validate_intercept_only_random_effects_brms <- function(vc, context = "MAIHDA variance calculations") {
  random_groups <- setdiff(names(vc), c("residual__", "sigma"))
  offending <- list()

  for (group in random_groups) {
    sd_tab <- vc[[group]]$sd
    if (is.null(dim(sd_tab))) {
      next
    }

    effect_names <- rownames(sd_tab)
    if (is.null(effect_names) && nrow(sd_tab) == 1) {
      next
    }

    non_intercepts <- maihda_non_intercept_effects(effect_names)
    if (length(non_intercepts) > 0) {
      offending[[group]] <- non_intercepts
    }
  }

  maihda_stop_for_random_slopes(offending, context)
  invisible(TRUE)
}

maihda_total_random_variance_brms <- function(model) {
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("Package 'brms' is required to work with brms models. Please install it with: install.packages('brms')")
  }

  vc <- brms::VarCorr(model)
  maihda_validate_intercept_only_random_effects_brms(vc)

  # Posterior-mean total random-effect variance E[sum sd^2] from the SD draws,
  # consistent with maihda_stratum_variance_brms(); falls back to summary SDs.
  draws <- tryCatch(maihda_posterior_draws_brms(model), error = function(e) NULL)
  if (!is.null(draws)) {
    sd_cols <- grep("^sd_", names(draws), value = TRUE)
    if (length(sd_cols) > 0) {
      total_per_draw <- Reduce(`+`, lapply(sd_cols, function(cn) as.numeric(draws[[cn]])^2))
      return(mean(total_per_draw, na.rm = TRUE))
    }
  }

  random_groups <- setdiff(names(vc), c("residual__", "sigma"))
  variances <- unlist(lapply(random_groups, function(group) {
    sd_tab <- vc[[group]]$sd
    if (is.null(dim(sd_tab)) || !"Estimate" %in% colnames(sd_tab)) {
      return(numeric())
    }
    vals <- as.numeric(sd_tab[, "Estimate"])^2
    vals[is.finite(vals)]
  }), use.names = FALSE)

  sum(variances, na.rm = TRUE)
}

maihda_variance_components_table <- function(var_stratum, var_other_random, var_residual) {
  var_other_random <- max(0, var_other_random, na.rm = TRUE)
  total_variance <- var_stratum + var_other_random + var_residual

  components <- "Between-stratum (random)"
  variances <- var_stratum

  if (is.finite(var_other_random) && var_other_random > sqrt(.Machine$double.eps)) {
    components <- c(components, "Other random effects")
    variances <- c(variances, var_other_random)
  }

  components <- c(components, "Within-stratum (residual)")
  variances <- c(variances, var_residual)

  proportions <- if (is.finite(total_variance) && total_variance > 0) {
    variances / total_variance
  } else {
    rep(NA_real_, length(variances))
  }

  data.frame(
    component = c(components, "Total"),
    variance = c(variances, total_variance),
    sd = sqrt(c(variances, total_variance)),
    proportion = c(proportions, 1.0)
  )
}

# ---- Cross-classified (group-aware) variance partition ---------------------
# These generalise the single-"stratum" extractors above to a model with several
# crossed intercept-only random effects, as fitted by the cross-classified MAIHDA
# decomposition: outcome ~ covars + (1 | dim1) + ... + (1 | stratum). Each dimension
# RE variance is that dimension's ADDITIVE main-effect variance; the intersection
# ("stratum") RE variance is the INTERACTION beyond additive. The default
# single-stratum path is untouched -- these are only reached when maihda() tags a
# model with $cc_info.

# Named numeric vector of every grouping factor's intercept variance (lme4). Unlike
# maihda_stratum_variance_lme4(), which returns only the "stratum" component, this
# returns one entry per random effect, keyed by the grouping-factor name.
maihda_random_variances_lme4 <- function(model) {
  vc <- lme4::VarCorr(model)
  vapply(names(vc), function(g) {
    group_mat <- as.matrix(vc[[g]])
    nm <- rownames(group_mat)
    icn <- intersect(c("(Intercept)", "Intercept"), nm)
    if (length(icn) == 0) NA_real_ else as.numeric(group_mat[icn[1], icn[1]])
  }, numeric(1))
}

# Per-draw intercept variances for every grouping factor (brms). Returns a named
# list of equal-length numeric vectors (one per random effect). brms names a
# group-level SD column sd_<group>__<coef>; the group name is the text before the
# first "__" (so dimension names containing "_" such as .maihda_dim_age are handled),
# and an intercept-only MAIHDA model contributes exactly one column per group.
maihda_group_variance_draws_brms <- function(draws) {
  if (!is.data.frame(draws)) {
    stop("'draws' must be a data frame of posterior draws.", call. = FALSE)
  }
  sd_cols <- grep("^sd_", names(draws), value = TRUE)
  if (length(sd_cols) == 0) {
    stop("No random-effect standard-deviation draws (sd_*) found in the brms posterior.")
  }
  bodies <- sub("^sd_", "", sd_cols)
  # brms names each group-level SD column sd_<group>__<coef>. The grouping-factor
  # name may itself contain "__" (e.g. a user context or dimension column literally
  # named site__id), so recover the group by splitting on the LAST "__": the group is
  # everything before it and the coefficient is the final segment. Splitting on the
  # FIRST "__" -- sub("__.*$", "", bodies) -- truncated site__id to site, so that
  # grouping was reported missing (the by-name lookups in maihda_cc_summary_brms())
  # or silently merged under the wrong key. brms coefficient names never contain
  # "__" (interactions use ":", powers use I()), so the last "__" is unambiguously
  # the group/coefficient separator, and two coefficients of one group still share a
  # group key -- preserving the random-slope rejection below.
  last_sep <- vapply(gregexpr("__", bodies, fixed = TRUE),
                     function(pos) if (pos[1L] == -1L) NA_integer_ else max(pos),
                     integer(1))
  groups <- ifelse(is.na(last_sep), bodies, substr(bodies, 1L, last_sep - 1L))
  out <- list()
  for (g in unique(groups)) {
    cols_g <- sd_cols[groups == g]
    if (length(cols_g) > 1L) {
      stop("The '", g, "' random effect must include only an intercept for MAIHDA ",
           "variance calculations (found multiple sd_", g,
           "__* draws, suggesting random slopes).", call. = FALSE)
    }
    out[[g]] <- as.numeric(draws[[cols_g]])^2
  }
  out
}

# Named numeric vector of posterior-mean intercept variances per group (brms),
# mirroring maihda_random_variances_lme4() for the point summary.
maihda_random_variances_brms <- function(model) {
  draws <- maihda_posterior_draws_brms(model)
  gv <- maihda_group_variance_draws_brms(draws)
  vapply(gv, mean, numeric(1))
}

# Split a named variance vector into additive (sum of the dimension REs), the
# interaction (the intersection RE), and the per-dimension vector relabelled by the
# human dimension name. dim_groups maps each strata_var to its RE grouping-factor
# name; interaction_group is the intersection RE name (usually "stratum"). Works on a
# named numeric vector (lme4 point) -- the per-draw brms path computes the same split
# directly from the draws list.
maihda_cc_variance_split <- function(var_named, dim_groups, interaction_group) {
  dim_re <- unname(dim_groups)
  needed <- c(dim_re, interaction_group)
  missing_re <- setdiff(needed, names(var_named))
  if (length(missing_re) > 0) {
    stop("Cross-classified variance partition is missing the random effect(s): ",
         paste(missing_re, collapse = ", "),
         ". Expected one intercept per dimension plus the intersection effect.",
         call. = FALSE)
  }
  per_dim <- var_named[dim_re]
  names(per_dim) <- names(dim_groups)
  list(
    per_dim = per_dim,
    additive = sum(per_dim),
    interaction = unname(var_named[interaction_group])
  )
}

# The crossed-dimensions partition arithmetic. Elementwise, so it serves the lme4
# point estimate (scalars), the lme4 bootstrap, and the brms posterior (per-draw
# vectors) identically. between = additive + interaction (total between-strata
# variance); vpc = between / (between + within + other); additive/interaction shares
# are the split of the between-strata variance (the crossed-dimensions analogue of
# the PCV). 'other' carries any random-effect variance outside the stratum structure
# -- e.g. a contextual (1 | school) intercept from fit_maihda(context = ) -- so the
# VPC stays the between-strata share of ALL unexplained variance, consistent with
# the single-stratum summary; it defaults to 0 for the plain decomposition.
maihda_cc_partition <- function(additive, interaction, within, other = 0) {
  between <- additive + interaction
  total <- between + within + other
  # The additive / interaction shares split the between-strata variance, so they are
  # undefined when there is no between-strata variance (between == 0): additive / 0
  # and interaction / 0 are 0 / 0 = NaN, a degenerate fit (every stratum's additive
  # and interaction variance estimated at exactly zero). Report NA rather than leaking
  # NaN to the public summary. Likewise the VPC is undefined when total variance is 0.
  safe_ratio <- function(num, den) ifelse(den > 0, num / den, NA_real_)
  list(
    additive = additive,
    interaction = interaction,
    between = between,
    within = within,
    other = other,
    total = total,
    vpc = safe_ratio(between, total),
    additive_share = safe_ratio(additive, between),
    interaction_share = safe_ratio(interaction, between)
  )
}

# Variance-components table for a crossed-dimensions summary: one (non-overlapping)
# row per dimension, then the intersection, any contextual random intercepts, the
# within-stratum (residual) term, and the total, so the proportions sum to 1 and
# plot_vpc() can read it directly. The additive subtotal and the shares live on the
# summary's $decomposition, not here, to avoid double counting. The table is tagged
# so plot_vpc() colours it as a crossed-dimensions split. (The attr value keeps the
# historical "cross_classified" spelling for compatibility with stored objects.)
maihda_cc_components_table <- function(per_dim, interaction_var, within_var,
                                       per_context = NULL) {
  components <- c(sprintf("Additive: %s", names(per_dim)),
                  "Intersectional interaction")
  variances <- c(unname(per_dim), interaction_var)
  if (!is.null(per_context) && length(per_context) > 0) {
    components <- c(components, sprintf("Context: %s", names(per_context)))
    variances <- c(variances, unname(per_context))
  }
  components <- c(components, "Within-stratum (residual)")
  variances <- c(variances, within_var)
  total_variance <- sum(variances)
  proportions <- if (is.finite(total_variance) && total_variance > 0) {
    variances / total_variance
  } else {
    rep(NA_real_, length(variances))
  }
  out <- data.frame(
    component = c(components, "Total"),
    variance = c(variances, total_variance),
    sd = sqrt(c(variances, total_variance)),
    proportion = c(proportions, 1.0),
    stringsAsFactors = FALSE
  )
  attr(out, "kind") <- "cross_classified"
  attr(out, "n_dimensions") <- length(per_dim)
  attr(out, "n_contexts") <- if (is.null(per_context)) 0L else length(per_context)
  out
}

# ---- Contextual cross-classified MAIHDA (stratum x place/institution) -------
# These support fit_maihda(context = ): individuals cross-classified by their
# intersectional stratum AND one or more higher-level contexts (school, hospital,
# region, ...), outcome ~ covars + (1 | stratum) + (1 | context). This is the
# literature's "cross-classified MAIHDA" (e.g. hospitals in Wemrell & Merlo's AMI
# study; schools in Prior et al.'s London students study) -- distinct from the
# crossed-DIMENSIONS decomposition above, which crosses the stratum dimensions'
# own main effects. The partition splits the unexplained variance into
# between-stratum vs. between-context vs. residual.

# Validate the `context` argument against the data: a character vector of existing
# column names, none of which may be (or collide with) the stratum machinery.
# Returns the validated (de-duplicated) vector, or NULL.
maihda_validate_context <- function(context, data) {
  if (is.null(context)) {
    return(NULL)
  }
  if (!is.character(context) || length(context) == 0 || anyNA(context) ||
      any(!nzchar(context))) {
    stop("'context' must be a character vector of column names in 'data' (e.g. ",
         "context = \"school\").", call. = FALSE)
  }
  context <- unique(context)
  missing_cols <- setdiff(context, names(data))
  if (length(missing_cols) > 0) {
    stop("Context variable(s) not found in data: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  if ("stratum" %in% context) {
    stop("'context' cannot include \"stratum\": the stratum random intercept is the ",
         "intersectional effect itself. Name the higher-level context column(s) ",
         "(e.g. school, hospital, region) instead.", call. = FALSE)
  }
  context
}

# The contextual partition arithmetic. Elementwise like maihda_cc_partition(), so
# it serves the lme4 point estimate (scalars), the lme4 bootstrap, and the brms
# posterior (per-draw vectors) identically. var_context is a NAMED list (or vector)
# with one element per context variable; var_other carries any further random
# effects outside stratum/context so the total still spans all unexplained variance.
# vpc_stratum is the headline MAIHDA VPC/ICC (the between-stratum share);
# vpc_context_total is the contexts' share (the general contextual effect).
maihda_context_partition <- function(var_stratum, var_context, var_residual,
                                     var_other = 0) {
  context_list <- if (is.list(var_context)) var_context else as.list(var_context)
  context_total <- Reduce(`+`, context_list)
  total <- var_stratum + context_total + var_residual + var_other
  list(
    stratum = var_stratum,
    context = context_list,
    context_total = context_total,
    residual = var_residual,
    other = var_other,
    total = total,
    vpc_stratum = var_stratum / total,
    vpc_context = lapply(context_list, function(v) v / total),
    vpc_context_total = context_total / total
  )
}

# Variance-components table for a contextual summary: the between-stratum row, one
# "Context: <name>" row per context, any other random effects, the residual, and
# the total. Keeps the canonical "Between-stratum (random)" / "Within-stratum
# (residual)" labels so existing consumers (compare_maihda_groups, plot_vpc) that
# match on them keep working; tagged kind = "contextual" for the plot layer.
maihda_context_components_table <- function(var_stratum, per_context,
                                            var_other_random, var_residual) {
  var_other_random <- max(0, var_other_random, na.rm = TRUE)
  components <- c("Between-stratum (random)",
                  sprintf("Context: %s", names(per_context)))
  variances <- c(var_stratum, unname(per_context))
  if (is.finite(var_other_random) && var_other_random > sqrt(.Machine$double.eps)) {
    components <- c(components, "Other random effects")
    variances <- c(variances, var_other_random)
  }
  components <- c(components, "Within-stratum (residual)")
  variances <- c(variances, var_residual)
  total_variance <- sum(variances)
  proportions <- if (is.finite(total_variance) && total_variance > 0) {
    variances / total_variance
  } else {
    rep(NA_real_, length(variances))
  }
  out <- data.frame(
    component = c(components, "Total"),
    variance = c(variances, total_variance),
    sd = sqrt(c(variances, total_variance)),
    proportion = c(proportions, 1.0),
    stringsAsFactors = FALSE
  )
  attr(out, "kind") <- "contextual"
  attr(out, "n_contexts") <- length(per_context)
  out
}

# Normalize the user-facing `decomposition` argument shared by maihda() and
# compare_maihda_groups(). The value "cross-classified" historically named the
# crossed-DIMENSIONS decomposition; it was renamed "crossed-dimensions" when the
# contextual cross-classified model (fit_maihda(context = ), the literature's
# cross-classified MAIHDA) was added, to free the term. The old value still works
# as a deprecated alias, with a one-time warning per call.
maihda_resolve_decomposition <- function(decomposition) {
  # A caller may pass its own multi-element formal default (e.g.
  # c("two-model", "crossed-dimensions")); take the first as the default, mirroring
  # match.arg()'s behaviour, so adding "longitudinal" to the choices below does not
  # break callers whose default vector no longer equals the full choice set.
  if (length(decomposition) > 1) {
    decomposition <- decomposition[1]
  }
  if (identical(decomposition, "cross-classified")) {
    warning("decomposition = \"cross-classified\" has been renamed ",
            "\"crossed-dimensions\" (it crosses the stratum dimensions' main ",
            "effects as random intercepts). The old value still works but is ",
            "deprecated -- \"cross-classified MAIHDA\" now refers to the ",
            "contextual stratum-by-place model fitted via the 'context' argument.",
            call. = FALSE)
    decomposition <- "crossed-dimensions"
  }
  match.arg(decomposition, c("two-model", "crossed-dimensions", "longitudinal"))
}

maihda_stratum_ranef_lme4 <- function(model, group = "stratum") {
  re <- lme4::ranef(model, condVar = TRUE)
  if (!group %in% names(re)) {
    stop("No '", group, "' random effects found in the model.")
  }

  group_re <- re[[group]]
  effect_names <- colnames(group_re)
  intercept_name <- intersect(c("(Intercept)", "Intercept"), effect_names)
  if (length(intercept_name) == 0) {
    stop("The '", group, "' random effect must include an intercept for MAIHDA stratum estimates.")
  }

  effect_idx <- match(intercept_name[1], effect_names)
  cond_var <- attr(group_re, "postVar")
  if (is.array(cond_var) && length(dim(cond_var)) == 3) {
    se <- sqrt(cond_var[effect_idx, effect_idx, ])
  } else {
    se <- rep(NA_real_, nrow(group_re))
  }

  random_effect <- group_re[[effect_idx]]
  data.frame(
    stratum = rownames(group_re),
    stratum_id = suppressWarnings(as.integer(rownames(group_re))),
    random_effect = random_effect,
    se = se,
    lower_95 = random_effect - 1.96 * se,
    upper_95 = random_effect + 1.96 * se,
    stringsAsFactors = FALSE
  )
}

maihda_stratum_ranef_brms <- function(model, group = "stratum") {
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("Package 'brms' is required to work with brms models. Please install it with: install.packages('brms')")
  }

  re <- brms::ranef(model, summary = TRUE)
  if (!group %in% names(re)) {
    stop("No '", group, "' random effects found in the brms model.")
  }

  group_re <- re[[group]]
  if (length(dim(group_re)) != 3) {
    stop("Unexpected brms random-effects shape; expected levels x summaries x effects.")
  }

  effect_names <- dimnames(group_re)[[3]]
  idx <- match(TRUE, effect_names %in% c("(Intercept)", "Intercept"))
  if (is.na(idx)) {
    if (length(effect_names) == 1) {
      idx <- 1
    } else {
      stop("The '", group, "' random effect must include an intercept for MAIHDA stratum estimates.")
    }
  }

  data.frame(
    stratum = dimnames(group_re)[[1]],
    stratum_id = suppressWarnings(as.integer(dimnames(group_re)[[1]])),
    random_effect = group_re[, "Estimate", idx],
    se = group_re[, "Est.Error", idx],
    lower_95 = group_re[, "Q2.5", idx],
    upper_95 = group_re[, "Q97.5", idx],
    stringsAsFactors = FALSE
  )
}

# Gaussian level-1 (residual) variance for VPC. For an unweighted lmer this is the
# usual sigma^2. With prior weights w_i the observation-level residual variance is
# sigma^2 / w_i, so there is no single residual variance; we report the mean
# conditional residual variance, mean_i(sigma^2 / w_i) = sigma^2 * mean(1 / w_i),
# as the level-1 variance (see the VPC note in summary.maihda_model()). Reduces to
# sigma^2 when all weights are 1.
maihda_gaussian_residual_variance_lme4 <- function(model, vc = lme4::VarCorr(model)) {
  sigma2 <- attr(vc, "sc")^2
  w <- tryCatch(stats::weights(model, type = "prior"), error = function(e) NULL)
  if (is.null(w) || length(w) == 0) {
    return(sigma2)
  }
  w <- as.numeric(w)
  w <- w[is.finite(w) & w > 0]
  if (length(w) == 0 || isTRUE(all(abs(w - 1) < sqrt(.Machine$double.eps)))) {
    return(sigma2)
  }
  sigma2 * mean(1 / w)
}

# Prior (precision/frequency) weights of a RAW fitted model object (glmerMod /
# brmsfit) as a plain numeric vector, or NULL when the fit is unweighted or the
# weights are unavailable. lme4 exposes them through stats::weights(type = "prior")
# (all 1 for an unweighted fit, so the count-VPC average reduces exactly to the
# unweighted mean); anything else falls through the tryCatch to NULL. Distinct from
# maihda_prior_weights(), which takes the maihda_model WRAPPER and aligns row weights
# (sampling/aggregated-binomial aware) for the plot/summary aggregations.
#
# brms exposes NO weights.brmsfit accessor, so stats::weights() would return NULL and
# a sampling-weighted brms fit would be read as unweighted -- biasing the population
# count-family (Poisson/NB) VPC, whose level-1 variance averages a per-row term that
# MUST be weighted (see maihda_residual_variance_draws_brms()). A sampling-weighted
# MAIHDA brms fit instead carries the (mean-1 normalized) likelihood weights in the
# reserved `.maihda_sw` data column, aligned to the same rows stats::fitted() returns,
# so read them straight off brmsfit$data; an unweighted brms fit has no such column
# and correctly degrades to NULL (the unweighted mean).
maihda_fit_prior_weights <- function(model) {
  if (inherits(model, "brmsfit")) {
    d <- tryCatch(model$data, error = function(e) NULL)
    if (is.data.frame(d) && ".maihda_sw" %in% names(d)) {
      w <- suppressWarnings(as.numeric(d[[".maihda_sw"]]))
      if (length(w) > 0) {
        return(w)
      }
    }
    return(NULL)
  }
  w <- tryCatch(stats::weights(model, type = "prior"), error = function(e) NULL)
  if (is.null(w) || length(w) == 0) {
    return(NULL)
  }
  as.numeric(w)
}

# Mean of a per-observation latent-scale quantity, weighted by prior weights when
# present. The count-family (Poisson/negative-binomial) level-1 VPC averages a
# per-row latent variance log1p(1/lambda + ...); with frequency weights w_i that average
# must be sum(w_i v_i) / sum(w_i) so a weighted fit reproduces the plain mean over the
# equivalent duplicated-row data. Reduces to mean(v, na.rm = TRUE) when weights are
# absent, length-mismatched, or all 1.
maihda_weighted_obs_mean <- function(v, w = NULL) {
  v <- as.numeric(v)
  if (is.null(w) || length(w) != length(v)) {
    return(mean(v, na.rm = TRUE))
  }
  w <- as.numeric(w)
  ok <- is.finite(v) & is.finite(w) & w > 0
  if (!any(ok) || isTRUE(all(abs(w[ok] - 1) < sqrt(.Machine$double.eps)))) {
    return(mean(v, na.rm = TRUE))
  }
  sum(w[ok] * v[ok]) / sum(w[ok])
}

# Theta (the negative-binomial size/dispersion parameter) of an lme4 NB fit.
# glmer.nb() stores its ML estimate retrievably via getME(); a fixed-theta
# glmer(family = MASS::negative.binomial(theta)) fit has no such slot, so fall
# back to parsing theta out of the family label "Negative Binomial(<theta>)"
# (read RAW via stats::family() -- maihda_family() would normalize the label
# away). Var(Y|u) = mu + mu^2/theta, so larger theta means less overdispersion.
maihda_negbin_theta_lme4 <- function(model) {
  th <- tryCatch(lme4::getME(model, "glmer.nb.theta"), error = function(e) NULL)
  if (is.numeric(th) && length(th) == 1 && is.finite(th) && th > 0) {
    return(th)
  }
  raw <- tryCatch(stats::family(model)$family, error = function(e) NULL)
  if (is.character(raw) && length(raw) == 1) {
    th <- suppressWarnings(
      as.numeric(sub("^Negative Binomial\\(([^)]+)\\)$", "\\1", raw)))
    if (is.finite(th) && th > 0) {
      return(th)
    }
  }
  stop("Could not recover the negative-binomial theta from the lme4 fit.",
       call. = FALSE)
}

# ---- marginal (population-averaged) count mean for the latent-scale VPC -----
# The count-family level-1 variances below -- Stryhn et al. (2006)
# log1p(1/lambda) for Poisson, Nakagawa, Johnson & Schielzeth (2017)
# log1p(1/lambda + 1/theta) for the negative binomial -- are defined at the
# MARGINAL expectation of the response, lambda_i = E[y_i] =
# exp(x_i'beta + v_i/2), with v_i the row's total latent random-effect variance
# (the log-normal mean correction), NOT at the conditional fitted mean
# exp(x_i'beta + u_hat) that stats::fitted() returns: the BLUPs would make the
# level-1 variance depend on the realized (and shrunken) random effects, which
# none of the cited references do. For the canonical null model every row's
# lambda is exactly Nakagawa's plug-in exp(beta0 + sigma^2/2); with covariates
# the per-row lambda_i extends it Austin-et-al.-style across the sample, and the
# callers average the per-row terms (prior-weighted where applicable).

# Per-observation latent-scale random-effect variance: sum over every
# random-effect term of z_i' Sigma_g z_i, where z_i is the row's RE design
# vector for that term and Sigma_g the term's covariance block. For the
# canonical intercept-only MAIHDA structures this is the constant sum of the
# intercept variances; for a longitudinal growth model it varies with the row's
# time value (z_i = (1, t_i, ...)). `bars` are reformulas::findbars() terms and
# `blocks` a parallel list of covariance blocks whose dimnames carry the design
# column names ("(Intercept)", the time term, ...). Pure (no fit object), so it
# serves lme4 and brms alike and is unit-testable Stan-free.
maihda_latent_re_variance_rows <- function(bars, blocks, data) {
  v <- numeric(nrow(data))
  for (k in seq_along(bars)) {
    Sg <- as.matrix(blocks[[k]])
    lhs_txt <- paste(deparse(bars[[k]][[2]]), collapse = " ")
    X <- stats::model.matrix(stats::as.formula(paste("~", lhs_txt)), data = data)
    if (nrow(X) != nrow(data)) {
      stop("Could not build the random-effect design rows for the marginal ",
           "count mean (rows were dropped).", call. = FALSE)
    }
    want <- rownames(Sg)
    idx <- match(want, colnames(X))
    if (anyNA(idx)) {
      stop("Could not align the random-effect design column(s) ",
           paste(want[is.na(idx)], collapse = ", "),
           " for the marginal count mean.", call. = FALSE)
    }
    Xk <- X[, idx, drop = FALSE]
    v <- v + unname(rowSums((Xk %*% Sg) * Xk))
  }
  v
}

# The grouping-factor name of a random-effect bar term, as VarCorr()/brms name
# it (no backticks around a non-syntactic column).
maihda_bar_group_name <- function(bar) {
  g <- bar[[3]]
  if (is.name(g)) as.character(g) else paste(deparse(g), collapse = " ")
}

# Marginal count mean lambda_i = exp(x_i'beta + v_i/2) of an lme4 count fit
# (see the section header). eta is the fixed-part linear predictor (offset
# included); v_i comes from the VarCorr blocks over the fitted model frame, so
# growth-model random slopes contribute their time-varying share.
maihda_count_marginal_mu_lme4 <- function(model) {
  eta <- as.numeric(stats::predict(model, re.form = NA, type = "link"))
  bars <- reformulas::findbars(stats::formula(model))
  if (is.null(bars) || length(bars) == 0) {
    stop("Could not read the random-effect terms off the lme4 formula for the ",
         "marginal count mean.", call. = FALSE)
  }
  vc <- lme4::VarCorr(model)
  blocks <- lapply(bars, function(b) {
    blk <- vc[[maihda_bar_group_name(b)]]
    if (is.null(blk)) {
      stop("No '", maihda_bar_group_name(b), "' variance block found for the ",
           "marginal count mean.", call. = FALSE)
    }
    blk
  })
  v <- maihda_latent_re_variance_rows(bars, blocks, maihda_model_frame(model))
  pmax(exp(eta + v / 2), .Machine$double.eps)
}

# Posterior-mean covariance block of a brms random-effect term, keyed by the
# DESIGN column names ("(Intercept)", the slope column, ...) so it aligns with
# maihda_latent_re_variance_rows(). Diagonals are E[sd^2]; off-diagonals
# E[cor * sd_i * sd_j] (brms orders cor_<g>__<a>__<b> by parameter order, so
# both orders are tried). Pure function of the draws, unit-testable Stan-free.
maihda_brms_re_block_mean <- function(draws, group, design_names) {
  to_par <- function(nm) if (identical(nm, "(Intercept)")) "Intercept" else nm
  k <- length(design_names)
  sds <- lapply(design_names, function(nm) {
    col <- paste0("sd_", group, "__", to_par(nm))
    if (!col %in% names(draws)) {
      stop("No 'sd_", group, "__", to_par(nm), "' draws found in the brms ",
           "posterior for the marginal count mean.", call. = FALSE)
    }
    as.numeric(draws[[col]])
  })
  S <- matrix(0, k, k, dimnames = list(design_names, design_names))
  for (i in seq_len(k)) {
    S[i, i] <- mean(sds[[i]]^2)
    for (j in seq_len(i - 1L)) {
      c1 <- paste0("cor_", group, "__", to_par(design_names[j]), "__",
                   to_par(design_names[i]))
      c2 <- paste0("cor_", group, "__", to_par(design_names[i]), "__",
                   to_par(design_names[j]))
      rho <- if (c1 %in% names(draws)) {
        as.numeric(draws[[c1]])
      } else if (c2 %in% names(draws)) {
        as.numeric(draws[[c2]])
      } else {
        0
      }
      S[i, j] <- S[j, i] <- mean(rho * sds[[i]] * sds[[j]])
    }
  }
  S
}

# Posterior-mean linear predictor of a brms fit, one value per prediction row.
# brms's posterior_linpred() returns the ndraws x nobs DRAWS matrix and has no
# summary argument: a `summary = TRUE` passed by the caller is silently
# swallowed by prepare_predictions() (brms hard-codes summary = FALSE in its
# internal posterior_epred() call), so the historical
# posterior_linpred(..., summary = TRUE)[, "Estimate"] idiom fails on the
# unnamed draws matrix ("no 'dimnames' attribute"). Collapse the draws to their
# posterior mean here -- the same estimand that idiom intended -- in one shared
# place. `...` is forwarded to posterior_linpred (newdata, re_formula,
# allow_new_levels, ...). If a brms version ever returns a summary matrix
# again, its "Estimate" column is that same posterior mean, so it is accepted.
maihda_brms_linpred_mean <- function(model, ...) {
  draws <- brms::posterior_linpred(model, ...)
  # A categorical() likelihood returns a 3-D ndraws x nobs x ncat array (one
  # linear predictor per response category), for which a single posterior-mean
  # linear predictor per row is undefined -- and colMeans() below would silently
  # flatten it to a wrong-length nobs * ncat vector. No family the MAIHDA
  # summaries support produces more than 2 dims; reject rather than return
  # garbage.
  if (length(dim(draws)) > 2L) {
    stop("brms::posterior_linpred() returned a ", length(dim(draws)),
         "-dimensional draws array (one linear predictor per response ",
         "category, as for a categorical() likelihood), so a single ",
         "posterior-mean linear predictor per row is not defined for this ",
         "model's family.", call. = FALSE)
  }
  if (!is.null(colnames(draws)) && "Estimate" %in% colnames(draws)) {
    return(as.numeric(draws[, "Estimate"]))
  }
  if (is.null(dim(draws))) {
    # Defensive: a dimensionless return can only be the draws of a single
    # prediction row (brms never returns one draw), so its mean is that row's
    # posterior-mean linear predictor.
    return(mean(as.numeric(draws)))
  }
  as.numeric(colMeans(draws))
}

# brms counterpart of maihda_count_marginal_mu_lme4(): lambda_i =
# exp(eta_i + v_i/2) with eta the posterior-mean fixed-part linear predictor
# and v_i built from the posterior-mean covariance blocks -- point-estimate
# plug-ins, consistent with the "mu held at its posterior mean" treatment the
# per-draw residual-variance path documents.
maihda_count_marginal_mu_brms <- function(model,
                                          draws = maihda_posterior_draws_brms(model)) {
  eta <- maihda_brms_linpred_mean(model, re_formula = NA)
  f <- model$formula
  if (inherits(f, "brmsformula") && inherits(f$formula, "formula")) {
    f <- f$formula
  }
  bars <- tryCatch(reformulas::findbars(f), error = function(e) NULL)
  if (is.null(bars) || length(bars) == 0) {
    stop("Could not read the random-effect terms off the brms formula for the ",
         "marginal count mean.", call. = FALSE)
  }
  data <- model$data
  blocks <- lapply(bars, function(b) {
    lhs_txt <- paste(deparse(b[[2]]), collapse = " ")
    design_names <- colnames(
      stats::model.matrix(stats::as.formula(paste("~", lhs_txt)), data = data))
    maihda_brms_re_block_mean(draws, maihda_bar_group_name(b), design_names)
  })
  v <- maihda_latent_re_variance_rows(bars, blocks, data)
  pmax(exp(eta + v / 2), .Machine$double.eps)
}

maihda_residual_variance_lme4 <- function(model, vc = lme4::VarCorr(model)) {
  fam <- maihda_family(model)
  if (is.null(fam)) {
    stop("Unable to determine model family for residual variance calculation.")
  }

  latent_families <- c("binomial", "bernoulli", "quasibinomial", "cumulative", "sratio", "cratio", "acat", "ordinal")
  if (fam$family == "gaussian") {
    # Only the identity link yields a coherent VPC: attr(vc, "sc")^2 is the
    # residual variance on the response scale, while the between-stratum variance
    # is on the linear-predictor (link) scale. With a non-identity link (e.g.
    # gaussian(link = "log")) the two live on different scales, so their ratio is
    # not a valid variance partition.
    maihda_stop_gaussian_non_identity_vpc(fam$link)
    return(maihda_gaussian_residual_variance_lme4(model, vc))
  }
  if (fam$family %in% latent_families && fam$link == "logit") {
    return((pi^2) / 3)
  }
  if (fam$family %in% latent_families && fam$link == "probit") {
    return(1)
  }
  if (fam$family == "poisson" && fam$link == "log") {
    # Stryhn et al. (2006) latent-scale level-1 variance approximation
    # log(1 + 1/lambda), evaluated at the MARGINAL expected count lambda_i =
    # exp(x_i'beta + v_i/2) (see maihda_count_marginal_mu_lme4) -- exactly
    # Nakagawa et al.'s exp(beta0 + sigma^2/2) plug-in for the null model --
    # not at the conditional fitted mean, whose BLUPs would tie the level-1
    # variance to the realized random effects. The simpler 1/lambda form is the
    # first-order Taylor expansion and matches log1p(1/lambda) only when
    # 1/lambda is small; for low-count outcomes it overestimates residual
    # variance and biases VPC downward. With prior weights the per-row
    # variances are averaged by those weights (see maihda_weighted_obs_mean()).
    mu <- maihda_count_marginal_mu_lme4(model)
    return(maihda_weighted_obs_mean(log1p(1 / mu), maihda_fit_prior_weights(model)))
  }
  if (fam$family == "negbinomial" && fam$link == "log") {
    # Negative-binomial analogue of the Poisson approximation above: the
    # lognormal latent-scale level-1 variance ln(1 + 1/lambda + 1/theta) of
    # Nakagawa, Johnson & Schielzeth (2017, J R Soc Interface 14:20170213),
    # likewise at the marginal expected count. The extra 1/theta term carries
    # the overdispersion, so it reduces to the Stryhn Poisson form as
    # theta -> Inf. Prior weights average the per-row terms.
    mu <- maihda_count_marginal_mu_lme4(model)
    theta <- maihda_negbin_theta_lme4(model)
    return(maihda_weighted_obs_mean(log1p(1 / mu + 1 / theta),
                                    maihda_fit_prior_weights(model)))
  }

  stop("VPC residual variance is not implemented for family '", fam$family,
       "' with link '", fam$link, "'.")
}

# A Gaussian model with a non-identity link mixes a response-scale residual
# variance with a link-scale between-stratum variance, so no valid VPC/ICC exists.
# Raise one clear, shared error rather than silently returning an invalid ratio.
maihda_stop_gaussian_non_identity_vpc <- function(link) {
  if (identical(link, "identity")) {
    return(invisible(NULL))
  }
  stop("VPC/ICC is not available for a Gaussian model with a non-identity link ",
       "(here '", link, "'): the residual variance is on the response scale while ",
       "the between-stratum variance is on the link scale, so their ratio is not a ",
       "valid variance partition. Refit with the identity link, or transform the ",
       "outcome and use gaussian(link = \"identity\").", call. = FALSE)
}

maihda_residual_variance_brms <- function(model) {
  fam <- maihda_family(model)
  if (is.null(fam)) {
    stop("Unable to determine brms model family for residual variance calculation.")
  }

  latent_families <- c("binomial", "bernoulli", "quasibinomial", "cumulative", "sratio", "cratio", "acat", "ordinal")
  if (fam$family == "gaussian") {
    # See maihda_residual_variance_lme4(): a non-identity Gaussian link mixes a
    # response-scale residual with a link-scale random effect, so the VPC is invalid.
    maihda_stop_gaussian_non_identity_vpc(fam$link)
    sigma_est <- tryCatch(stats::sigma(model), error = function(e) NA_real_)
    if (length(sigma_est) > 0 && is.finite(sigma_est[1])) {
      return(as.numeric(sigma_est[1]^2))
    }
    vc <- brms::VarCorr(model)
    residual_name <- intersect(c("residual__", "sigma"), names(vc))
    if (length(residual_name) > 0) {
      return(as.numeric(vc[[residual_name[1]]]$sd[1, "Estimate"]^2))
    }
  }
  if (fam$family %in% latent_families && fam$link == "logit") {
    return((pi^2) / 3)
  }
  if (fam$family %in% latent_families && fam$link == "probit") {
    return(1)
  }
  if (fam$family == "poisson" && fam$link == "log") {
    # Stryhn et al. (2006) latent-scale level-1 variance approximation
    # log(1 + 1/lambda), at the MARGINAL expected count (posterior-mean fixed
    # part + lognormal RE-variance correction; see maihda_count_marginal_mu_brms
    # and the lme4 sibling for the rationale). Likelihood weights average the
    # per-row terms (see maihda_weighted_obs_mean()).
    mu <- maihda_count_marginal_mu_brms(model)
    return(maihda_weighted_obs_mean(log1p(1 / mu), maihda_fit_prior_weights(model)))
  }
  if (fam$family == "negbinomial" && fam$link == "log") {
    # Nakagawa, Johnson & Schielzeth (2017) lognormal latent-scale level-1
    # variance ln(1 + 1/lambda + 1/theta); brms's 'shape' parameter IS theta.
    # Evaluated at the marginal expected counts (posterior-mean plug-ins) and
    # posterior-mean shape (see maihda_residual_variance_draws_brms for the
    # per-draw treatment).
    draws <- maihda_posterior_draws_brms(model)
    mu <- maihda_count_marginal_mu_brms(model, draws)
    if (!"shape" %in% names(draws)) {
      stop("Could not extract the negative-binomial 'shape' (theta) draws from ",
           "the brms posterior.")
    }
    shape <- mean(as.numeric(draws[["shape"]]), na.rm = TRUE)
    return(maihda_weighted_obs_mean(log1p(1 / mu + 1 / shape),
                                    maihda_fit_prior_weights(model)))
  }

  stop("VPC residual variance is not implemented for brms family '", fam$family,
       "' with link '", fam$link, "'.")
}

# ---- brms VPC/ICC from posterior draws -------------------------------------
# The helpers above collapse the posterior to summary SDs and so report E[sd]^2;
# the helpers below instead work draw-by-draw (E[sd^2]) and yield a credible
# interval for the VPC/ICC. summary.maihda_model() uses the draws-based path.

# Posterior draws of a brms model as a data frame. Column names follow the
# standard Stan/brms convention: group-level SDs are `sd_<group>__<coef>`
# (e.g. `sd_stratum__Intercept`), the Gaussian residual SD is `sigma`, and
# population-level effects are `b_<coef>`. This is equivalent to
# posterior::as_draws_df(model) for our purposes but relies only on the
# brms-exported as.data.frame() method, so it adds no hard dependency.
maihda_posterior_draws_brms <- function(model) {
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("Package 'brms' is required to work with brms models. Please install it with: install.packages('brms')")
  }

  draws <- tryCatch(as.data.frame(model), error = function(e) NULL)
  if (is.null(draws) || !is.data.frame(draws) || nrow(draws) == 0) {
    stop("Could not extract posterior draws from the brms model.")
  }

  draws
}

# Per-draw random-effect variances from a brms posterior draws data frame.
# Returns a list of equal-length numeric vectors:
#   $stratum - between-stratum variance per draw (sd_<group>__Intercept squared)
#   $total   - summed variance across every group-level SD per draw
#   $other   - total minus stratum (variance of any non-stratum random effects)
# This is a pure function of the draws data frame, so it can be unit tested with
# a hand-built data frame without fitting a Stan model.
maihda_random_variance_draws_brms <- function(draws, group = "stratum") {
  if (!is.data.frame(draws)) {
    stop("'draws' must be a data frame of posterior draws.", call. = FALSE)
  }

  sd_cols <- grep("^sd_", names(draws), value = TRUE)
  if (length(sd_cols) == 0) {
    stop("No random-effect standard-deviation draws (sd_*) found in the brms posterior.")
  }

  # Stratum-group SD columns, matched by an exact "<group>__" prefix rather than
  # a regex so that group names containing metacharacters are treated literally.
  # Intercept-only models contribute exactly one column (sd_<group>__Intercept);
  # more than one means random slopes, which MAIHDA VPC/ICC does not support.
  bodies <- sub("^sd_", "", sd_cols)
  group_cols <- sd_cols[startsWith(bodies, paste0(group, "__"))]
  if (length(group_cols) == 0L) {
    stop("No '", group, "' random-effect SD draws (sd_", group,
         "__*) found in the brms posterior.")
  }
  if (length(group_cols) > 1L) {
    stop("The '", group, "' random effect must include only an intercept for MAIHDA ",
         "variance calculations (found multiple sd_", group,
         "__* draws, suggesting random slopes).")
  }

  var_stratum <- as.numeric(draws[[group_cols]])^2
  var_total <- Reduce(`+`, lapply(sd_cols, function(cn) as.numeric(draws[[cn]])^2))
  var_other <- pmax(0, var_total - var_stratum)

  list(stratum = var_stratum, total = var_total, other = var_other)
}

# Per-draw level-1 (residual) variance for a brms model. Mirrors the latent /
# distributional choices of maihda_residual_variance_brms() but returns a
# per-draw vector:
#   gaussian       -> sigma draws squared (exact)
#   logit latent   -> pi^2 / 3 (constant across draws)
#   probit latent  -> 1        (constant across draws)
#   poisson(log)   -> mean(log1p(1 / lambda)) at the MARGINAL expected counts
#                     lambda_{d,i} = exp(eta_{d,i} + v_d / 2). For the
#                     intercept-only VPC structures (strata / crossed / contextual)
#                     this is computed PER DRAW from the fixed-part linear-predictor
#                     draws and the draw's total random-intercept variance
#                     (maihda_count_residual_variance_draws_brms_perdraw), so the
#                     interval reflects fixed-effect and RE-variance uncertainty and
#                     their correlation with the level-1 term. A random-slope
#                     (longitudinal) structure, whose row design the fast path does
#                     not carry, falls back to the posterior-mean plug-in marginal
#                     mean held constant across draws (maihda_count_marginal_mu_brms).
#   negbinomial(log) -> mean(log1p(1 / lambda + 1 / shape_d)) per draw d: the
#                     'shape' (= theta) draws are always propagated; lambda is
#                     computed per draw for the intercept-only path and held at the
#                     posterior-mean plug-in for the random-slope fallback, exactly
#                     as in the poisson case.
# Takes the model (for the family and, for poisson/negbinomial, the marginal
# means) and the draws data frame (for the sigma/shape draws and the draw count).
# Per-draw count-family latent level-1 (residual) variance from the fixed-part
# linear-predictor DRAWS. For each posterior draw d it forms the marginal expected
# count lambda_{d,i} = exp(eta_{d,i} + v_d / 2) -- eta the fixed-part (re_formula =
# NA) link-scale linear predictor and v_d the draw's total random-INTERCEPT variance
# -- then returns the (optionally prior-weighted) mean over observations of
# log1p(1 / lambda_{d,i} [+ extra_d]). `extra_d` carries the negative-binomial
# 1 / shape term (NULL for Poisson). lambda is clamped away from 0 exactly as the
# posterior-mean helpers do (pmax(., .Machine$double.eps)). Propagating lambda per
# draw -- rather than holding it at a posterior-mean plug-in -- lets the count VPC
# credible interval reflect fixed-effect and RE-variance uncertainty and their
# correlation with the level-1 term. Pure (eta_link is a plain ndraws x nobs
# matrix), so it is unit-testable without a Stan fit.
maihda_count_resid_var_from_linpred <- function(eta_link, var_total, w = NULL,
                                                extra = NULL) {
  eta_link <- as.matrix(eta_link)
  nd <- nrow(eta_link)
  if (length(var_total) != nd) {
    stop("var_total must have one value per posterior draw (row of eta_link).",
         call. = FALSE)
  }
  # Marginal mean per draw/row, clamped away from 0; then 1 / lambda.
  mu <- exp(sweep(eta_link, 1L, var_total / 2, "+"))
  mu[mu < .Machine$double.eps] <- .Machine$double.eps
  inv_mu <- 1 / mu
  if (!is.null(extra)) {
    if (length(extra) != nd) {
      stop("extra must have one value per posterior draw.", call. = FALSE)
    }
    inv_mu <- sweep(inv_mu, 1L, extra, "+")
  }
  terms <- log1p(inv_mu)                                  # ndraws x nobs
  if (is.null(w)) {
    return(rowMeans(terms, na.rm = TRUE))
  }
  vapply(seq_len(nd),
         function(d) maihda_weighted_obs_mean(terms[d, ], w),
         numeric(1))
}

# brms gatherer for maihda_count_resid_var_from_linpred(): returns the per-draw
# count residual variance when the random effects are INTERCEPT-ONLY (the strata /
# crossed-dimensions / contextual VPC structures), else NULL to signal the caller to
# fall back to the documented posterior-mean plug-in. A random-slope growth model
# needs the row-varying random-effect design (z_i' Sigma z_i) that this fast path
# does not carry, so it is deliberately excluded here. `extra` is a length-ndraws
# vector (negative-binomial 1 / shape) or NULL (Poisson). Returns NULL rather than
# erroring on any structural mismatch so the caller's plug-in path always remains
# available.
maihda_count_residual_variance_draws_brms_perdraw <- function(model, draws,
                                                              extra = NULL) {
  f <- model$formula
  if (inherits(f, "brmsformula") && inherits(f$formula, "formula")) {
    f <- f$formula
  }
  bars <- tryCatch(reformulas::findbars(f), error = function(e) NULL)
  if (is.null(bars) || length(bars) == 0) {
    return(NULL)
  }
  is_intercept_only <- all(vapply(bars, function(b) {
    identical(paste(deparse(b[[2]]), collapse = " "), "1")
  }, logical(1)))
  if (!is_intercept_only) {
    return(NULL)
  }

  sd_cols <- grep("^sd_", names(draws), value = TRUE)
  if (length(sd_cols) == 0) {
    return(NULL)
  }
  var_total <- Reduce(`+`, lapply(sd_cols, function(cn) as.numeric(draws[[cn]])^2))

  eta_link <- tryCatch(brms::posterior_linpred(model, re_formula = NA),
                       error = function(e) NULL)
  if (is.null(eta_link)) {
    return(NULL)
  }
  eta_link <- as.matrix(eta_link)
  # posterior_linpred returns ndraws x nobs in canonical draw order, matching the
  # as.data.frame() draw order that produced `draws`; require the row counts agree
  # before trusting that alignment, else fall back.
  if (length(dim(eta_link)) != 2L || nrow(eta_link) != length(var_total)) {
    return(NULL)
  }

  w <- maihda_fit_prior_weights(model)
  maihda_count_resid_var_from_linpred(eta_link, var_total, w = w, extra = extra)
}

maihda_residual_variance_draws_brms <- function(model, draws) {
  fam <- maihda_family(model)
  if (is.null(fam)) {
    stop("Unable to determine brms model family for residual variance calculation.")
  }

  n <- nrow(draws)
  latent_families <- c("binomial", "bernoulli", "quasibinomial", "cumulative",
                       "sratio", "cratio", "acat", "ordinal")

  if (fam$family == "gaussian") {
    # A non-identity Gaussian link mixes scales (see maihda_residual_variance_lme4()).
    maihda_stop_gaussian_non_identity_vpc(fam$link)
    if ("sigma" %in% names(draws)) {
      return(as.numeric(draws[["sigma"]])^2)
    }
    sigma_est <- tryCatch(stats::sigma(model), error = function(e) NA_real_)
    if (length(sigma_est) > 0 && is.finite(sigma_est[1])) {
      return(rep(as.numeric(sigma_est[1])^2, n))
    }
    stop("Could not extract residual SD draws ('sigma') from the gaussian brms model.")
  }
  if (fam$family %in% latent_families && fam$link == "logit") {
    return(rep((pi^2) / 3, n))
  }
  if (fam$family %in% latent_families && fam$link == "probit") {
    return(rep(1, n))
  }
  if (fam$family == "poisson" && fam$link == "log") {
    # Per-draw marginal mean for the intercept-only VPC path (propagates
    # fixed-effect + RE-variance uncertainty); NULL means fall back to the
    # documented posterior-mean plug-in (random-slope / longitudinal structure).
    rv <- maihda_count_residual_variance_draws_brms_perdraw(model, draws)
    if (!is.null(rv)) {
      return(rv)
    }
    mu <- maihda_count_marginal_mu_brms(model, draws)
    w <- maihda_fit_prior_weights(model)
    return(rep(maihda_weighted_obs_mean(log1p(1 / mu), w), n))
  }
  if (fam$family == "negbinomial" && fam$link == "log") {
    # Nakagawa et al. (2017) ln(1 + 1/lambda + 1/theta), theta = brms 'shape'.
    # The shape draws are always propagated; lambda is computed per draw for the
    # intercept-only path and held at the posterior-mean plug-in otherwise (see
    # the header comment). Likelihood weights average the per-row terms.
    if (!"shape" %in% names(draws)) {
      stop("Could not extract the negative-binomial 'shape' (theta) draws from ",
           "the brms posterior.")
    }
    shape_d <- as.numeric(draws[["shape"]])
    rv <- maihda_count_residual_variance_draws_brms_perdraw(model, draws,
                                                            extra = 1 / shape_d)
    if (!is.null(rv)) {
      return(rv)
    }
    mu <- maihda_count_marginal_mu_brms(model, draws)
    w <- maihda_fit_prior_weights(model)
    return(vapply(shape_d,
                  function(s) maihda_weighted_obs_mean(log1p(1 / mu + 1 / s), w),
                  numeric(1)))
  }

  stop("VPC residual variance is not implemented for brms family '", fam$family,
       "' with link '", fam$link, "'.")
}

# Summarise a posterior VPC/ICC from per-draw variance components. var_stratum,
# var_other_random and var_residual are per-draw numeric vectors (length-1
# scalars are recycled). Returns the posterior point estimate (median by
# default) and a central credible interval at conf_level. Pure and Stan-free.
maihda_vpc_posterior_summary <- function(var_stratum, var_other_random,
                                         var_residual, conf_level = 0.95,
                                         point = c("median", "mean")) {
  point <- match.arg(point)

  n <- max(length(var_stratum), length(var_other_random), length(var_residual))
  recycle <- function(v) if (length(v) == 1L) rep(v, n) else v
  var_stratum <- recycle(var_stratum)
  var_other_random <- recycle(var_other_random)
  var_residual <- recycle(var_residual)
  if (length(var_stratum) != n || length(var_other_random) != n ||
      length(var_residual) != n) {
    stop("Per-draw variance vectors must share a common length (or be scalars).",
         call. = FALSE)
  }

  total <- var_stratum + var_other_random + var_residual
  vpc <- var_stratum / total
  vpc <- vpc[is.finite(vpc)]
  if (length(vpc) == 0) {
    stop("No finite VPC draws were available to summarise the posterior VPC/ICC.",
         call. = FALSE)
  }

  alpha <- 1 - conf_level
  estimate <- if (point == "median") stats::median(vpc) else mean(vpc)
  ci <- stats::quantile(vpc, probs = c(alpha / 2, 1 - alpha / 2), names = FALSE)

  list(
    estimate = estimate,
    ci_lower = ci[1],
    ci_upper = ci[2],
    conf_level = conf_level,
    n_draws = length(vpc)
  )
}

# Orchestrates the draws-based brms VPC/ICC used by summary.maihda_model():
# validates the intercept-only random-effect structure, extracts per-draw
# variances, and returns both the posterior VPC summary (estimate + credible
# interval) and the posterior-mean variance components for the components table.
maihda_vpc_draws_brms <- function(model, conf_level = 0.95, group = "stratum",
                                  point = c("median", "mean")) {
  point <- match.arg(point)
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("Package 'brms' is required to work with brms models. Please install it with: install.packages('brms')")
  }

  maihda_validate_intercept_only_random_effects_brms(brms::VarCorr(model))

  draws <- maihda_posterior_draws_brms(model)
  rv <- maihda_random_variance_draws_brms(draws, group = group)
  var_residual <- maihda_residual_variance_draws_brms(model, draws)

  vpc <- maihda_vpc_posterior_summary(
    rv$stratum, rv$other, var_residual,
    conf_level = conf_level, point = point
  )

  list(
    vpc = vpc,
    var_stratum = mean(rv$stratum),
    var_other_random = mean(rv$other),
    var_residual = mean(var_residual)
  )
}

# TRUE when a $vpc list carries a finite lower/upper interval (an lme4 bootstrap
# CI or a brms posterior credible interval). The print methods use this to
# decide whether to show an interval alongside the point estimate.
maihda_vpc_has_interval <- function(vpc) {
  !is.null(vpc) &&
    !is.null(vpc$ci_lower) && !is.null(vpc$ci_upper) &&
    length(vpc$ci_lower) == 1L && length(vpc$ci_upper) == 1L &&
    is.finite(vpc$ci_lower) && is.finite(vpc$ci_upper)
}

# Human-readable label for the VPC interval, distinguishing the lme4 parametric
# bootstrap CI from the brms posterior credible interval.
maihda_vpc_interval_label <- function(vpc) {
  conf_pct <- if (!is.null(vpc$conf_level)) vpc$conf_level * 100 else 95
  if (identical(vpc$method, "posterior")) {
    sprintf("(Posterior %.0f%% credible interval)", conf_pct)
  } else {
    sprintf("(Bootstrap %.0f%% CI)", conf_pct)
  }
}

# Prior/precision weights aligned to a model's analytic rows (object$data is the
# model frame, so its row order matches weights(model, type = "prior")). Returns
# rep(1, n) when the model is unweighted or the weights cannot be recovered (e.g.
# brms), so weighted aggregations below reduce EXACTLY to the unweighted ones for
# the common case. Used so stratum-level plot summaries honour the lme4 prior
# weights, consistent with the weighted Gaussian VPC. These are prior/precision
# weights, NOT a complex survey design -- no design-based variance is computed and
# results are not survey-representative.
maihda_prior_weights <- function(object) {
  n <- nrow(object$data)
  # Design-weighted fits (sampling_weights supplied): aggregate with the SAMPLING
  # weights, so stratum-level summaries are population-representative under the
  # weights. The wemix analytic frame keeps the original weight column; a brms
  # fit's model frame instead carries the normalized .maihda_sw column (mean-1
  # scaling does not change weighted means).
  if (!is.null(object$sampling_weights)) {
    sw_col <- if (object$sampling_weights %in% names(object$data)) {
      object$sampling_weights
    } else if (".maihda_sw" %in% names(object$data)) {
      ".maihda_sw"
    } else {
      NULL
    }
    if (!is.null(sw_col)) {
      w <- as.numeric(object$data[[sw_col]])
      if (length(w) == n) {
        w[!is.finite(w)] <- NA_real_
        return(w)
      }
    }
  }
  # Aggregated binomial responses (an lme4 cbind(success, failure) matrix, or a
  # brms `y | trials(n)`) carry the binomial TRIALS per row. For an OBSERVED
  # stratum summary that aggregation is numerator/denominator based (sum of
  # successes / sum of trials), so the trials are already in the denominator and
  # multiplying the rows by them would double-count -- such models are unit-weighted
  # HERE. (An lme4 cbind() fit's two-column matrix response is recognised below; a
  # brms trials() fit has a single-column response and no weights.brmsfit, so it
  # falls through to unit weights too -- the same correct answer for this path.)
  # Averaging fitted probabilities / linear predictors is different: a row with more
  # trials carries more information and must be trial-weighted, so the prediction
  # path uses maihda_prediction_weights() (below) instead of this.
  resp <- tryCatch(stats::model.response(object$data), error = function(e) NULL)
  if (!is.null(resp) && !is.null(dim(resp))) {
    return(rep(1, n))
  }
  w <- tryCatch(stats::weights(object$model, type = "prior"), error = function(e) NULL)
  if (is.null(w) || length(w) != n) {
    return(rep(1, n))
  }
  w <- as.numeric(w)
  w[!is.finite(w)] <- NA_real_
  w
}

# The argument expression of the first `trials(...)` call inside a brms response
# addition term (the RHS of `y | ...`), or NULL when there is none. Recurses so a
# combined term such as `trials(n) + weights(w)` is handled regardless of order.
maihda_find_trials_expr <- function(expr) {
  if (!is.call(expr)) {
    return(NULL)
  }
  if (identical(expr[[1]], as.name("trials")) && length(expr) >= 2L) {
    return(expr[[2]])
  }
  for (i in seq_along(expr)[-1]) {
    found <- maihda_find_trials_expr(expr[[i]])
    if (!is.null(found)) {
      return(found)
    }
  }
  NULL
}

# Binomial TRIAL counts of a brms `y | trials(n)` fit, or NULL when the fit is not
# a brms aggregated-binomial model. By default the counts are aligned to
# object$data; pass `data` (e.g. a prediction newdata) to evaluate the trials()
# term on a different frame instead, so a response-scale prediction can be
# normalised by the trial counts of the rows it was made on. brms exposes
# model.frame.brmsfit but NO weights.brmsfit, so -- unlike an lme4 cbind() fit,
# whose trials come through stats::weights(type = "prior") -- the counts are not
# recoverable that way and a brms trials() fit would silently fall back to unit
# row weights. Parse the trials() addition term off the stored formula instead and
# evaluate it on the supplied frame, so brms trials() fits are trial-weighted like
# their lme4 cbind() counterparts.
maihda_brms_trial_counts <- function(object, data = object$data) {
  if (!inherits(object$model, "brmsfit")) {
    return(NULL)
  }
  f <- object$formula
  if (inherits(f, "brmsformula") && inherits(f$formula, "formula")) {
    f <- f$formula
  }
  if (!inherits(f, "formula") || length(f) < 3L) {
    return(NULL)
  }
  lhs <- f[[2]]
  if (!is.call(lhs) || !identical(lhs[[1]], as.name("|"))) {
    return(NULL)
  }
  trials_expr <- maihda_find_trials_expr(lhs[[3]])
  if (is.null(trials_expr)) {
    return(NULL)
  }
  vals <- tryCatch(eval(trials_expr, envir = data, enclos = baseenv()),
                   error = function(e) NULL)
  if (is.null(vals)) {
    return(NULL)
  }
  as.numeric(vals)
}

# Convert a brms response-scale prediction to the probability scale for an
# aggregated-binomial `y | trials(n)` fit. brms's fitted()/posterior_epred()
# return the expected COUNT of successes (trials * p) for such a fit, whereas
# lme4's type = "response" for a cbind(success, failure) fit returns the per-trial
# probability; divide by the per-row trial counts (evaluated on the prediction
# `newdata`) so both engines return a probability on the response scale. A
# Bernoulli / non-binomial fit has no trials() term, so the counts are NULL and the
# estimate is returned unchanged. The counts must be finite and positive to divide;
# otherwise the estimate is left as-is rather than producing Inf/NaN.
maihda_brms_response_to_prob <- function(object, est, newdata) {
  trials <- maihda_brms_trial_counts(object, newdata)
  if (is.null(trials) || length(trials) != length(est) ||
      !all(is.finite(trials)) || any(trials <= 0)) {
    return(est)
  }
  est / trials
}

# Weights for averaging PREDICTIONS (fitted probabilities / linear predictors) over
# the rows of a stratum, as distinct from the observed numerator/denominator
# summaries that maihda_prior_weights() serves. Identical to maihda_prior_weights()
# EXCEPT for aggregated-binomial responses (cbind(success, failure) / `y |
# trials(n)`), where each row aggregates a different number of Bernoulli trials: the
# row's fitted probability/linear predictor must be weighted by the binomial TRIAL
# counts, not by unit weights, or a 2-trial row and a 200-trial row would count
# equally in the stratum mean. The counts come from stats::weights(type = "prior")
# for an lme4 cbind() fit; for a brms `y | trials(n)` fit (no weights.brmsfit) they
# are parsed from the formula's trials() term by maihda_brms_trial_counts(). For
# non-binomial fits the two functions coincide, so unweighted and prior/sampling-
# weighted models are unaffected. Trial counts are folded on top of any sampling
# weights (a sampling-weighted aggregated-binomial row represents sw population
# units, each contributing `trials` draws). These remain prior/precision (and,
# where present, sampling) weights -- NOT a complex survey design; no design-based
# variance is computed.
maihda_prediction_weights <- function(object) {
  w <- maihda_prior_weights(object)
  n <- length(w)

  # brms `y | trials(n)`: the response is a single column (not a cbind() matrix)
  # and there is no weights.brmsfit, so the lme4 trials path below misses it. Read
  # the trial counts off the formula and fold them onto any sampling weights,
  # exactly as the lme4 cbind() branch folds prior-weight trials below.
  brms_trials <- maihda_brms_trial_counts(object)
  if (!is.null(brms_trials) && length(brms_trials) == n &&
      all(is.finite(brms_trials))) {
    return(w * brms_trials)
  }

  resp <- tryCatch(stats::model.response(object$data), error = function(e) NULL)
  if (is.null(resp) || is.null(dim(resp))) {
    return(w)
  }
  trials <- tryCatch(stats::weights(object$model, type = "prior"),
                     error = function(e) NULL)
  if (is.null(trials) || length(trials) != n || any(!is.finite(trials))) {
    return(w)
  }
  w * as.numeric(trials)
}

# Per-stratum prior-weight-weighted means of the named columns of `pred_df` (which
# must carry `stratum` and `weight` columns). Returns a data frame with `stratum`,
# the weighted columns, integer `n` (row count) and `w_sum` (sum of weights),
# ordered by stratum. Reduces exactly to unweighted means when all weights are
# equal, so unweighted models are unaffected.
maihda_weighted_stratum_aggregate <- function(pred_df, cols) {
  groups <- sort(unique(pred_df$stratum))
  wmean <- function(col) {
    vapply(groups, function(g) {
      sel <- pred_df$stratum == g
      stats::weighted.mean(pred_df[[col]][sel], pred_df$weight[sel], na.rm = TRUE)
    }, numeric(1))
  }
  out <- data.frame(stratum = groups, stringsAsFactors = FALSE)
  for (col in cols) {
    out[[col]] <- wmean(col)
  }
  out$n <- as.integer(vapply(groups, function(g) sum(pred_df$stratum == g), integer(1)))
  out$w_sum <- vapply(groups,
                      function(g) sum(pred_df$weight[pred_df$stratum == g], na.rm = TRUE),
                      numeric(1))
  out
}

# Build a random-effects prediction formula (~ (lhs | g) + ...) that keeps ONLY
# the grouping factors named in `groups`, for a re.form-/re_formula-scoped
# predict(). Returns NA (predict with the fixed effects only) when no term
# matches or no group is requested. Group names are compared in deparsed form and
# also against their backtick-quoted form, so a non-syntactic grouping name
# matches either way (mirrors maihda_da_re_scopes()). Used to build the
# crossed-dimensions stratum-prediction baseline from the additive dimension REs
# alone -- excluding the interaction stratum RE (added back per row) and any
# non-intersectional random effect (e.g. a contextual/site RE), which must not
# leak into the ranked per-stratum prediction.
maihda_re_form_for_groups <- function(formula, groups) {
  bars <- reformulas::findbars(formula)
  groups <- groups[!is.na(groups) & nzchar(groups)]
  if (is.null(bars) || length(bars) == 0 || length(groups) == 0) {
    return(NA)
  }
  quoted <- vapply(groups, maihda_quote_name, character(1))
  want <- unique(c(groups, quoted))
  keep <- vapply(bars, function(b) {
    paste(deparse(b[[3]]), collapse = "") %in% want
  }, logical(1))
  if (!any(keep)) {
    return(NA)
  }
  re_txt <- vapply(bars[keep], function(b) {
    paste0("(", paste(deparse(b), collapse = " "), ")")
  }, character(1))
  stats::as.formula(paste("~", paste(re_txt, collapse = " + ")))
}

# Per-fitted-row offset of an lme4 fit, read off its model frame: lme4 stores an
# external offset= as the "(offset)" column and a formula offset() term as its
# "offset(...)" column, so both are recovered (and summed, matching
# stats::model.offset()). Returns NULL when the fit carries no offset. The TOTAL stored
# offset is only valid for predictions on the UNMODIFIED fitted rows; the longitudinal
# grid paths, which move the time column, instead split it into
# maihda_fitted_offset_external() (frozen -- no expression to re-evaluate) plus
# maihda_lme4_formula_offset_at() (re-evaluated, so a time-dependent formula offset
# tracks the modified time).
maihda_fitted_offset <- function(model) {
  mf <- maihda_model_frame(model)
  if (is.null(mf)) {
    return(NULL)
  }
  cols <- grep("^\\(offset\\)$|^offset\\(", names(mf), value = TRUE)
  if (length(cols) == 0) {
    return(NULL)
  }
  Reduce(`+`, lapply(cols, function(cc) as.numeric(mf[[cc]])))
}

# TRUE when a raw lme4 fit carries an EXTERNAL offset (the offset= fitting argument)
# rather than a formula offset() term: lme4 names the external-offset model-frame
# column "(offset)", whereas a formula offset() is stored under "offset(...)". The
# distinction matters for newdata predictions -- predict.merMod re-evaluates a formula
# offset() from newdata but cannot recover an external one. (predict_maihda.R has the
# maihda_model-wrapper analogue maihda_lme4_has_external_offset(), which reads the same
# "(offset)" column off object$data.)
maihda_mermod_has_external_offset <- function(model) {
  mf <- tryCatch(maihda_model_frame(model), error = function(e) NULL)
  !is.null(mf) && "(offset)" %in% names(mf)
}

# The EXTERNAL offset= contribution only -- the "(offset)" model-frame column, per fitted
# row; NULL when the fit used no external offset. Unlike a formula offset() term it has no
# stored expression to re-evaluate on a modified-time grid, so the longitudinal
# trajectories freeze it (a constant-exposure assumption). The formula-offset counterpart
# is maihda_lme4_formula_offset_at(), which DOES re-evaluate; splitting the two lets a
# time-dependent formula offset track the trajectory time while an external offset stays
# fixed.
maihda_fitted_offset_external <- function(model) {
  mf <- maihda_model_frame(model)
  if (is.null(mf) || !"(offset)" %in% names(mf)) {
    return(NULL)
  }
  as.numeric(mf[["(offset)"]])
}

# Re-evaluate a fit's FORMULA offset() term(s) on `newdata`, returning the summed per-row
# offset (NULL when the fit carries no formula offset). Unlike the stored offset
# (maihda_fitted_offset()), this tracks a time-dependent offset such as offset(0.5 * time)
# to the times in `newdata`: the longitudinal prediction grids vary time (and hold other
# covariates at observed or representative values), so a frozen offset would mis-place
# every count / link-scale trajectory. offset() itself is a no-op wrapper (identity), so
# it is bound to identity in a private evaluation frame; the raw variables resolve from
# `newdata`, everything else from the formula's environment. An offset term whose raw
# variables are NOT all present in `newdata` cannot be re-evaluated: the model frame
# stores such a term only under its derived "offset(...)" column -- e.g. offset(logE)
# where logE is not otherwise a predictor -- and it is necessarily time-invariant (a
# time-referencing term would find its time column on these grids). It falls back to the
# fit's stored per-row values instead: as-is when `newdata` aligns with the fitted rows
# (fallback = "rows", the marginalizing VPC(t) grid), or their mean for a
# representative-profile grid (fallback = "mean", the fixed trajectory).
maihda_lme4_formula_offset_at <- function(model, newdata,
                                          fallback = c("rows", "mean")) {
  fallback <- match.arg(fallback)
  tt <- stats::delete.response(stats::terms(maihda_nobars(stats::formula(model))))
  off_idx <- attr(tt, "offset")
  if (is.null(off_idx) || length(off_idx) == 0L) {
    return(NULL)
  }
  varlist <- attr(tt, "variables")
  env <- environment(stats::formula(model))
  off_env <- new.env(parent = if (is.null(env)) baseenv() else env)
  off_env$offset <- function(x) x
  mf <- maihda_model_frame(model)
  n <- nrow(newdata)
  vals <- lapply(off_idx, function(i) {
    expr <- varlist[[i + 1L]]
    if (all(all.vars(expr) %in% names(newdata))) {
      v <- as.numeric(eval(expr, newdata, off_env))
      return(if (length(v) == 1L) rep(v, n) else v)
    }
    stored <- if (is.null(mf)) NULL else mf[[paste(deparse(expr), collapse = " ")]]
    if (is.null(stored)) {
      stop("Cannot evaluate the formula offset term '",
           paste(deparse(expr), collapse = " "), "' on the prediction grid: its ",
           "variable(s) are absent from the data and no stored offset column was ",
           "found.", call. = FALSE)
    }
    stored <- as.numeric(stored)
    if (identical(fallback, "mean")) rep(mean(stored, na.rm = TRUE), n) else stored
  })
  Reduce(`+`, vals)
}

# Fixed-effects (re.form = NA) link-scale linear predictor of an lme4 fit on arbitrary
# newdata, WITH the supplied offset added. predict.merMod(newdata = ) drops an external
# offset (the offset= fitting argument) and ERRORS on a formula offset() term whose raw
# variable is absent from newdata -- which is exactly the case for the package's own
# prediction grids, built from the stored model frame (where the offset lives only as
# the derived "offset(...)"/"(offset)" column). So build X %*% beta from the fixed-
# effect design directly, excluding the offset from the model matrix (model.matrix()
# never adds an offset column), and add the caller-supplied offset. The model's terms
# object is used (not a rebuilt formula) so predvars-dependent terms such as poly()
# evaluate on the fitted basis; a formula offset()'s raw variable is filled with a
# harmless dummy where absent (it never enters the model matrix). Matches
# predict(., re.form = NA, type = "link") exactly for a no-offset fit. `offset` is a
# per-row vector (recycled if scalar) or NULL.
maihda_lme4_fixed_link <- function(model, newdata, offset = NULL) {
  beta <- lme4::fixef(model)
  tt <- stats::delete.response(stats::terms(maihda_nobars(stats::formula(model))))
  off_idx <- attr(tt, "offset")
  if (!is.null(off_idx) && length(off_idx) > 0) {
    varlist <- attr(tt, "variables")
    off_vars <- unlist(lapply(off_idx, function(i) all.vars(varlist[[i + 1L]])))
    for (v in setdiff(off_vars, names(newdata))) newdata[[v]] <- 1
  }
  mfr <- maihda_model_frame(model)
  term_vars <- all.vars(tt)
  fac <- names(mfr)[vapply(mfr, is.factor, logical(1))]
  xlev <- lapply(mfr[intersect(fac, term_vars)], levels)
  mf <- stats::model.frame(tt, newdata, na.action = stats::na.pass, xlev = xlev)
  contr <- attr(lme4::getME(model, "X"), "contrasts")
  X <- stats::model.matrix(tt, mf, contrasts.arg = contr)
  missing_cols <- setdiff(names(beta), colnames(X))
  if (length(missing_cols) > 0) {
    stop("Could not rebuild the fixed-effects design for prediction; missing ",
         "column(s): ", paste(missing_cols, collapse = ", "), ".", call. = FALSE)
  }
  eta <- as.numeric(X[, names(beta), drop = FALSE] %*% beta)
  if (!is.null(offset)) {
    eta <- eta + if (length(offset) == 1L) rep(offset, length(eta)) else offset
  }
  eta
}

maihda_stratum_predictions_lme4 <- function(object, summary_obj, scale = c("response", "link")) {
  scale <- match.arg(scale)
  if (!is.null(object$longitudinal_info)) {
    # Backstop: a single per-stratum prediction adds only the random INTERCEPT and
    # drops the random slope, so it is a cross-sectional value, not the trajectory.
    maihda_stop_longitudinal_scalar("A single per-stratum prediction")
  }
  data <- object$data
  if (!"stratum" %in% names(data)) {
    stop("'stratum' variable not found in fitted model data.")
  }

  model <- object$model
  fam <- maihda_family(model)
  linkinv <- maihda_linkinv(fam)
  prior_w <- maihda_prediction_weights(object)

  # Predict on the FITTED rows WITHOUT newdata so lme4 reuses its stored linear
  # predictor, which INCLUDES any offset. Passing newdata = data (the stored model
  # frame) would drop an external offset= and error on a formula offset() term (its
  # raw variable lives only as the frame's derived "offset(...)"/"(offset)" column).
  # object$data IS the model frame, so predict()'s output is in the same row order as
  # `data` / `key` below.
  eta_fixed <- stats::predict(model, re.form = NA, type = "link")
  stratum_est <- summary_obj$stratum_estimates
  if (is.null(stratum_est) || nrow(stratum_est) == 0) {
    stop("No stratum estimates available.")
  }

  key <- as.character(data$stratum)
  re_key <- as.character(stratum_est$stratum)
  idx <- match(key, re_key)

  transform_eta <- function(eta) {
    if (scale == "response") linkinv(eta) else eta
  }

  # Baseline the intersection (stratum) random effect is added to. For the canonical
  # single-stratum model this is the fixed-only linear predictor. For a
  # cross-classified model it must also carry the ADDITIVE dimension random effects
  # so the stratum prediction is fixed effects + dimension REs + interaction RE. It
  # must NOT carry any non-intersectional random effect (a contextual/site RE, which
  # a crossed-dimensions fit can also include): predicting all REs and subtracting
  # only u_stratum would leave that contextual effect in the baseline, so the ranked
  # stratum prediction would depend on each stratum's observed context composition
  # rather than the intended intersectional scope. Scope the baseline to the
  # dimension REs only (re.form); the interaction RE is added per row below.
  eta_base <- eta_fixed
  if (!is.null(object$cc_info)) {
    re_form <- maihda_re_form_for_groups(stats::formula(model),
                                         unname(object$cc_info$dim_groups))
    # As above: no newdata, so the stored offset is retained and the rows stay aligned.
    eta_base <- as.numeric(stats::predict(model, re.form = re_form, type = "link"))
  }

  pred_df <- data.frame(
    stratum = key,
    predicted_row = transform_eta(eta_base + stratum_est$random_effect[idx]),
    lower_row = transform_eta(eta_base + stratum_est$lower_95[idx]),
    upper_row = transform_eta(eta_base + stratum_est$upper_95[idx]),
    fixed_row = transform_eta(eta_base),
    weight = prior_w,
    stringsAsFactors = FALSE
  )

  # Prediction-weighted per-stratum means: a weighted fit's stratum predictions
  # reflect the prior/precision (or sampling) weights, and an aggregated-binomial
  # fit weights each row by its binomial TRIAL count so rows with more trials count
  # proportionally more (consistent with the weighted VPC). Identical to the plain
  # unweighted means for an ordinary unweighted model. These are prior/precision (and
  # trial) weights, not a complex survey design (not survey-representative).
  maihda_weighted_stratum_aggregate(
    pred_df, c("predicted_row", "lower_row", "upper_row", "fixed_row")
  )
}

maihda_apply_autobin_info <- function(strata_data, autobin_info) {
  if (is.null(autobin_info) || length(autobin_info) == 0) {
    return(strata_data)
  }

  for (v in intersect(names(autobin_info), names(strata_data))) {
    info <- autobin_info[[v]]
    if (!is.null(info$breaks) && !is.null(info$labels)) {
      strata_data[[v]] <- cut(strata_data[[v]], breaks = info$breaks,
                              include.lowest = TRUE, labels = info$labels)
    }
  }

  strata_data
}

# Recreate the .maihda_dim_<var> tertile factors that maihda_adjusted_terms() adds
# when a numeric stratum dimension was auto-binned. The adjusted / cross-classified
# model formulae reference these columns as fixed effects, so prediction newdata must
# carry them too; rebuild each from the stored breaks/labels (leaving the original
# numeric column intact, exactly as maihda_adjusted_terms() does). Out-of-range
# numeric values cut() to NA -- maihda_autobin_out_of_range() reports them separately.
maihda_add_binned_dim_columns <- function(data, autobin_info) {
  if (is.null(autobin_info) || length(autobin_info) == 0) {
    return(data)
  }

  for (v in intersect(names(autobin_info), names(data))) {
    info <- autobin_info[[v]]
    if (!is.null(info$breaks) && !is.null(info$labels)) {
      # Overwriting any '.maihda_dim_<v>' already in newdata is intentional: it is
      # the model's internal term, reserved at fit time (make_strata() rejects a
      # user column of this name), so predicting on the fitted data -- which already
      # carries the package's own copy -- rebuilds the identical factor.
      data[[maihda_dim_col(v)]] <- cut(
        data[[v]], breaks = info$breaks,
        include.lowest = TRUE, labels = info$labels
      )
    }
  }

  data
}

maihda_autobin_out_of_range <- function(data, autobin_info) {
  if (is.null(autobin_info) || length(autobin_info) == 0) {
    return(character())
  }

  out <- character()
  for (v in intersect(names(autobin_info), names(data))) {
    info <- autobin_info[[v]]
    if (is.null(info$breaks) || !is.numeric(data[[v]])) {
      next
    }

    x <- data[[v]]
    breaks <- range(info$breaks, na.rm = TRUE)
    bad <- !is.na(x) & (x < breaks[1] | x > breaks[2])
    if (any(bad)) {
      out <- c(out, paste0(v, " outside [", signif(breaks[1], 6), ", ",
                           signif(breaks[2], 6), "]"))
    }
  }

  out
}

# Row-wise display labels from stratum-defining columns, joining each row's
# values with `sep`. Each column is converted with as.character() individually --
# NOT via apply(), whose as.matrix() coercion runs numeric columns through
# format() in a mixed-type frame and pads them to a common width (values 1 and
# 10 would label as " 1" / "10", leaving a stray space in "m \u00d7  1"). Using
# as.character() also matches the value encoding of the stratum matcher
# (maihda_match_strata_rows), so labels and matching agree. unname() keeps a
# stratum variable literally named "sep"/"collapse" from being captured by
# paste()'s own arguments.
maihda_paste_label_rows <- function(df, sep) {
  cols <- unname(lapply(df, as.character))
  do.call(paste, c(cols, list(sep = sep)))
}

maihda_stratum_labels <- function(data, vars, sep = " \u00d7 ", autobin_info = NULL) {
  strata_data <- data[, vars, drop = FALSE]
  strata_data <- maihda_apply_autobin_info(strata_data, autobin_info)

  has_missing <- apply(strata_data, 1, function(x) any(is.na(x)))
  labels <- rep(NA_character_, nrow(strata_data))
  labels[!has_missing] <- maihda_paste_label_rows(
    strata_data[!has_missing, , drop = FALSE], sep
  )
  labels
}

maihda_stratum_lookup <- function(data, strata_info, vars, sep = " \u00d7 ",
                                  autobin_info = NULL) {
  strata_data <- data[, vars, drop = FALSE]
  strata_data <- maihda_apply_autobin_info(strata_data, autobin_info)

  has_missing <- apply(strata_data, 1, function(x) any(is.na(x)))
  out <- rep(NA_character_, nrow(strata_data))

  if (!all(vars %in% names(strata_info))) {
    labels <- maihda_stratum_labels(data, vars, sep, autobin_info)
    stratum_map <- stats::setNames(as.character(strata_info$stratum), strata_info$label)
    out <- unname(stratum_map[labels])
    return(out)
  }

  complete_idx <- which(!has_missing)
  if (length(complete_idx) == 0) {
    return(out)
  }

  matches <- maihda_match_strata_rows(
    strata_data[complete_idx, , drop = FALSE],
    strata_info[, vars, drop = FALSE],
    vars
  )
  matched <- !is.na(matches)
  out[complete_idx[matched]] <- as.character(strata_info$stratum[matches[matched]])
  out
}

# Strata the fitted model actually saw: the unique, non-missing stratum values in
# the analytic data (the levels for which a random effect was estimated). A
# stratum outside this set is "new" and has no estimated effect. Returns NULL when
# the analytic data carries no stratum column, in which case callers skip the
# unseen-stratum check (rather than reject everything).
maihda_known_strata <- function(object) {
  if (is.null(object$data) || !"stratum" %in% names(object$data)) {
    return(NULL)
  }
  s <- unique(as.character(object$data$stratum))
  s[!is.na(s)]
}

# Error if `stratum` names any level the model never saw. `type` only tailors the
# hint: an individual-level prediction can fall back to the population average via
# allow_new_levels = TRUE, but a stratum-level prediction has no random effect to
# report for an unseen stratum, so no such fallback is offered there.
maihda_check_known_strata <- function(stratum, known, type = "individual") {
  if (is.null(known)) {
    return(invisible(NULL))
  }
  hint <- if (identical(type, "individual")) {
    " Pass allow_new_levels = TRUE for a population-average (fixed-effects-only) prediction."
  } else {
    ""
  }
  # A missing (NA) stratum has no estimated random effect. For an INDIVIDUAL
  # prediction the WeMix/ordinal linpred helpers would silently map it to zero (a
  # population-average prediction), whereas lme4 rejects an NA grouping level
  # outright -- so reject it here to keep every engine consistent, unless the
  # caller opted into the population-average fallback (that path skips this check
  # entirely; see maihda_prepare_prediction_data()). A STRATUM-level prediction
  # instead simply yields no row for an NA stratum (an empty result, matching
  # lme4), so it is left to fall through.
  if (identical(type, "individual") && anyNA(stratum)) {
    stop("newdata has row(s) with a missing ('NA') stratum, so no stratum random ",
         "effect can be applied.", hint, call. = FALSE)
  }
  wanted <- unique(as.character(stratum))
  wanted <- wanted[!is.na(wanted)]
  unknown <- setdiff(wanted, known)
  if (length(unknown) == 0) {
    return(invisible(NULL))
  }
  stop("newdata contains strata not present in the fitted model: ",
       paste(utils::head(unknown, 5), collapse = ", "),
       if (length(unknown) > 5) ", ..." else "", ".", hint,
       call. = FALSE)
}

maihda_prepare_prediction_data <- function(object, newdata, type = "individual",
                                           allow_new_levels = FALSE) {
  if (!is.data.frame(newdata)) {
    stop("'newdata' must be a data frame", call. = FALSE)
  }
  # The adjusted / cross-classified model formula references reconstructed
  # .maihda_dim_<var> tertile factors for any auto-binned numeric dimension; recreate
  # them from the stored breaks/labels so newdata carries the same fixed-effect columns
  # as the fitting data. Done before the early return so it applies even when the
  # caller already supplies a 'stratum' column (as long as the numeric source columns
  # are present).
  newdata <- maihda_add_binned_dim_columns(newdata, object$strata_autobin_info)

  # A longitudinal fit on internally centered time references the derived
  # .maihda_ctime column; rebuild it from the original time column so newdata in
  # original time units predicts correctly. Always recomputed (idempotent when
  # predicting on the fitted data, which already carries the package's own copy),
  # so a caller can never double-center by passing rows copied from object$data.
  lng <- object$longitudinal_info
  if (!is.null(lng)) {
    time_term <- maihda_lng_time_term(lng)
    if (!identical(time_term, lng$time)) {
      if (lng$time %in% names(newdata)) {
        newdata[[time_term]] <- newdata[[lng$time]] - maihda_lng_time_center(lng)
      } else if (!time_term %in% names(newdata)) {
        stop("Cannot predict from this longitudinal model: newdata must contain ",
             "the time column '", lng$time, "'.", call. = FALSE)
      }
    }
  }

  # allow_new_levels only relaxes individual-level predictions (which fall back to
  # the population average); stratum-level predictions always require known strata.
  permit_new <- isTRUE(allow_new_levels) && identical(type, "individual")

  if ("stratum" %in% names(newdata)) {
    # A caller-supplied 'stratum' column previously skipped validation entirely,
    # so a misspelled or genuinely new stratum silently flowed through to a
    # fixed-only prediction (the WeMix/ordinal helpers map an unseen stratum's
    # random effect to 0). Validate it here unless the caller opted into
    # population-average predictions.
    if (!permit_new) {
      maihda_check_known_strata(newdata$stratum, maihda_known_strata(object), type)
    }
    return(newdata)
  }

  strata_info <- object$strata_info
  strata_vars <- object$strata_vars
  if (is.null(strata_vars) || length(strata_vars) == 0) {
    strata_vars <- maihda_infer_strata_vars(strata_info)
  }
  if (is.null(strata_vars) || length(strata_vars) == 0) {
    return(newdata)
  }

  missing_vars <- setdiff(strata_vars, names(newdata))
  if (length(missing_vars) > 0) {
    stop("Cannot rebuild 'stratum' for prediction. Missing grouping variables in newdata: ",
         paste(missing_vars, collapse = ", "), call. = FALSE)
  }

  if (is.null(strata_info) || !all(c("stratum", "label") %in% names(strata_info))) {
    stop("Cannot rebuild 'stratum' for prediction because training strata labels were not stored.",
         call. = FALSE)
  }

  sep <- object$strata_sep
  if (is.null(sep) || length(sep) != 1) {
    sep <- " \u00d7 "
  }

  out_of_range <- maihda_autobin_out_of_range(newdata, object$strata_autobin_info)
  if (length(out_of_range) > 0) {
    stop("Cannot rebuild 'stratum' for prediction because numeric grouping values ",
         "fall outside the training auto-bin ranges: ",
         paste(out_of_range, collapse = "; "), ".",
         call. = FALSE)
  }

  labels <- maihda_stratum_labels(newdata, strata_vars, sep, object$strata_autobin_info)
  newdata$stratum <- maihda_stratum_lookup(
    newdata,
    strata_info,
    strata_vars,
    sep,
    object$strata_autobin_info
  )

  # A row missing a stratum-defining dimension yields no stratum label at all
  # (labels[i] is NA), so it is neither a known nor a "new" stratum -- its stratum
  # stays NA. For an INDIVIDUAL prediction the WeMix/ordinal helpers would map that
  # (absent) random effect to zero, silently returning a population-average
  # prediction where lme4 rejects an NA grouping level; reject it here unless the
  # caller opted into the fallback (permit_new keeps the NA stratum, which the
  # engines then treat as fixed-effects-only). A stratum-level prediction yields no
  # row for such a stratum downstream, so it is left to fall through.
  if (identical(type, "individual") && !permit_new && anyNA(labels)) {
    stop(sprintf(paste0("newdata has %d row(s) whose stratum-defining variable(s) ",
                        "are missing (NA), so their stratum cannot be determined."),
                 sum(is.na(labels))),
         " Pass allow_new_levels = TRUE for a population-average ",
         "(fixed-effects-only) prediction.", call. = FALSE)
  }
  unknown <- !is.na(labels) & is.na(newdata$stratum)
  if (any(unknown)) {
    if (permit_new) {
      # Keep each new combination as its own stratum label so the engine maps it
      # to a zero random effect (a population-average prediction) rather than
      # erroring; mirrors the supplied-'stratum' branch above.
      newdata$stratum[unknown] <- labels[unknown]
    } else {
      unknown_labels <- unique(labels[unknown])
      hint <- if (identical(type, "individual")) {
        " Pass allow_new_levels = TRUE for a population-average (fixed-effects-only) prediction."
      } else {
        ""
      }
      stop("newdata contains strata combinations that were not present when the model was fit: ",
           paste(utils::head(unknown_labels, 5), collapse = ", "),
           if (length(unknown_labels) > 5) ", ..." else "", ".", hint,
           call. = FALSE)
    }
  }

  newdata
}

maihda_stratum_predictions_brms <- function(object, summary_obj, scale = c("response", "link")) {
  scale <- match.arg(scale)
  if (!is.null(object$longitudinal_info)) {
    # Backstop: as in the lme4 helper, the scalar per-stratum prediction drops the
    # random slope and so misrepresents a growth model's trajectory estimand.
    maihda_stop_longitudinal_scalar("A single per-stratum prediction")
  }
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("Package 'brms' is required. Please install it with: install.packages('brms')")
  }

  data <- object$data
  if (!"stratum" %in% names(data)) {
    stop("'stratum' variable not found in fitted model data.")
  }

  model <- object$model
  fam <- maihda_family(model)
  linkinv <- maihda_linkinv(fam)
  prior_w <- maihda_prediction_weights(object)
  eta_fixed <- maihda_brms_linpred_mean(model, newdata = data, re_formula = NA)

  # A cumulative (ordinal) brms fit's response scale is the EXPECTED CATEGORY
  # SCORE sum_k k * P(Y = k) in [1, K], NOT the scalar inverse link (which would
  # return a single cumulative probability in [0, 1]). Mirror the clmm path and
  # the individual brms path in predict_maihda(): map the latent location to the
  # score with the shared cumulative helpers. posterior_linpred() returns the
  # location mu (thresholds excluded), so P(Y <= k) = g^-1(alpha_k - eta) holds
  # with the posterior-mean thresholds. As on the clmm path this plugs point
  # estimates into the nonlinear score (the stratum prediction is built from the
  # stratum effect's point summaries anyway), rather than averaging the score
  # over the posterior the way the per-individual fitted() array does.
  is_ordinal <- maihda_family_is_ordinal(fam)
  if (is_ordinal) {
    ord_thresholds <- maihda_brms_ordinal_thresholds(model)
    ord_link <- fam$link
  }

  stratum_est <- summary_obj$stratum_estimates
  if (is.null(stratum_est) || nrow(stratum_est) == 0) {
    stop("No stratum estimates available.")
  }

  key <- as.character(data$stratum)
  re_key <- as.character(stratum_est$stratum)
  idx <- match(key, re_key)

  transform_eta <- function(eta) {
    if (scale != "response") {
      return(eta)
    }
    if (is_ordinal) {
      maihda_ordinal_eta_to_score(eta, ord_thresholds, ord_link)
    } else {
      linkinv(eta)
    }
  }

  # See the lme4 sibling: in a cross-classified model the stratum prediction is
  # fixed effects + ADDITIVE dimension REs + interaction RE, and must exclude any
  # non-intersectional random effect (a contextual/site RE the fit may also carry),
  # which predicting all REs and subtracting only u_stratum would leave in the
  # baseline and leak into the ranked stratum prediction. Scope the baseline to the
  # dimension REs only (re_formula); the interaction RE is added per row below.
  eta_base <- eta_fixed
  if (!is.null(object$cc_info)) {
    f <- model$formula
    if (inherits(f, "brmsformula") && inherits(f$formula, "formula")) {
      f <- f$formula
    }
    re_form <- maihda_re_form_for_groups(f, unname(object$cc_info$dim_groups))
    eta_base <- maihda_brms_linpred_mean(model, newdata = data,
                                         re_formula = re_form)
  }

  pred_df <- data.frame(
    stratum = key,
    predicted_row = transform_eta(eta_base + stratum_est$random_effect[idx]),
    lower_row = transform_eta(eta_base + stratum_est$lower_95[idx]),
    upper_row = transform_eta(eta_base + stratum_est$upper_95[idx]),
    fixed_row = transform_eta(eta_base),
    weight = prior_w,
    stringsAsFactors = FALSE
  )

  # Prediction-weighted per-stratum means: a weighted fit's stratum predictions
  # reflect the prior/precision (or sampling) weights, and an aggregated-binomial
  # fit weights each row by its binomial TRIAL count so rows with more trials count
  # proportionally more (consistent with the weighted VPC). Identical to the plain
  # unweighted means for an ordinary unweighted model. These are prior/precision (and
  # trial) weights, not a complex survey design (not survey-representative).
  maihda_weighted_stratum_aggregate(
    pred_df, c("predicted_row", "lower_row", "upper_row", "fixed_row")
  )
}

maihda_add_strata_columns <- function(data, strata_info) {
  if (is.null(strata_info) || !"stratum" %in% names(strata_info)) {
    return(data)
  }

  idx <- match(as.character(data$stratum), as.character(strata_info$stratum))
  extra_cols <- setdiff(names(strata_info), "stratum")
  for (col in extra_cols) {
    if (!col %in% names(data)) {
      data[[col]] <- strata_info[[col]][idx]
    }
  }
  data
}
