# Internal helpers shared across model summaries, PVC, predictions, and plots.

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

# "family(link)" key for a maihda_model, used to decide whether two models are
# comparable (same family and link). Prefers the fitted object's family and
# falls back to the family the wrapper recorded at fit time -- stats::family()
# is undefined for engines like wemix (WeMixResults). Names are canonical via
# maihda_family()/maihda_normalize_family_name(), so e.g. two glmer.nb() fits
# with different estimated thetas still compare equal.
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
      diagnostics$converged <- length(msgs) == 0
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
    header <- if (identical(diagnostics$engine, "brms")) {
      "MCMC convergence diagnostics (brms):"
    } else {
      "Convergence warnings reported by lme4:"
    }
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
  cat("Fit diagnostics:\n")
  cat(paste0("  ", diag_lines), sep = "\n")
  cat("\n\n")
  invisible(NULL)
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
  # Report the Monte Carlo error of the interval so users can judge whether enough
  # draws were used. n_ok is the effective number of successful refits; mc_se is the
  # Monte Carlo standard error of the bootstrap distribution's centre
  # (sd / sqrt(n_ok)), a simple proxy for how much the interval would move with a
  # different bootstrap seed. The print methods surface these.
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

  list(n_boot = as.integer(n_boot), conf_level = maihda_validate_conf_level(conf_level))
}

maihda_quote_name <- function(name) {
  if (!is.character(name) || length(name) != 1 || is.na(name) || name == "") {
    stop("Variable names must be non-empty character strings.", call. = FALSE)
  }

  paste(deparse(as.name(name), backtick = TRUE), collapse = "")
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

# Logical keep-mask over `n` rows reproducing the row selection lme4/brms apply
# before fitting: an (already-evaluated) `subset` value -- logical (recycled,
# NA -> drop), positive/negative numeric indices, or character row names -- and
# the removal of rows whose (already-evaluated) prior `weights` is NA. Subset and
# weights arrive as VALUES (resolved with the data mask in fit_maihda) so this
# helper never needs to evaluate user expressions.
maihda_row_mask <- function(data, subset = NULL, weights = NULL) {
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
    idx <- tryCatch(seq_len(n)[subset], error = function(e) integer(0))
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
                                        weights = NULL) {
  keep <- maihda_row_mask(data, subset = subset, weights = weights)
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

# The model response over the analytic sample (post-transformation, post-NA,
# post-subset and post-weight-NA). Only plain-symbol responses qualify as a
# Bernoulli candidate, so a transformed or aggregated response (log(y),
# cbind(s, f), `y | trials(n)`) yields NULL -- "not a single two-level response".
maihda_analytic_response <- function(formula, data, subset = NULL,
                                     weights = NULL) {
  if (length(formula) != 3L || !is.symbol(formula[[2]])) {
    return(NULL)
  }
  fr <- maihda_analytic_model_frame(formula, data, subset = subset,
                                    weights = weights)
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
                                             weights = NULL) {
  response <- formula[[2]]
  if (!is.symbol(response)) {
    return(data)
  }

  outcome <- as.character(response)
  # Recode against the analytic sample lme4/glmer actually fits (transformations
  # applied; rows excluded by `subset`, a missing weight, or missingness dropped),
  # matching the binary detection in fit_maihda(). A character/factor outcome whose
  # third value appears only on excluded rows is still recoded to 0/1; the
  # out-of-sample value becomes NA (and is dropped) rather than left as a stray
  # level that glmer() would reject.
  resp <- maihda_analytic_response(formula, data, subset = subset,
                                   weights = weights)
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
                                      weights = NULL) {
  resp <- maihda_analytic_response(formula, data, subset = subset,
                                   weights = weights)
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
                                       weights = NULL) {
  resp <- maihda_analytic_response(formula, data, subset = subset,
                                   weights = weights)
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
# comparisons/PVC across different datasets that happen to share n, row names and
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

# Fingerprint of a model's prior weights, so PVC / model comparisons do not
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
  # calculate_pvc()/stepwise_pcv() consistent with the draws-based VPC reported by
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
  groups <- sub("__.*$", "", bodies)
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
  list(
    additive = additive,
    interaction = interaction,
    between = between,
    within = within,
    other = other,
    total = total,
    vpc = between / total,
    additive_share = additive / between,
    interaction_share = interaction / between
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
  if (identical(decomposition, "cross-classified")) {
    warning("decomposition = \"cross-classified\" has been renamed ",
            "\"crossed-dimensions\" (it crosses the stratum dimensions' main ",
            "effects as random intercepts). The old value still works but is ",
            "deprecated -- \"cross-classified MAIHDA\" now refers to the ",
            "contextual stratum-by-place model fitted via the 'context' argument.",
            call. = FALSE)
    decomposition <- "crossed-dimensions"
  }
  match.arg(decomposition, c("two-model", "crossed-dimensions"))
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
    # Stryhn et al. (2006) latent-scale level-1 variance approximation: log(1 + 1/mu).
    # The simpler 1/mu form is the first-order Taylor expansion and matches log1p(1/mu)
    # only when 1/mu is small (i.e. mu large); for low-count Poisson outcomes (mu < ~2)
    # it overestimates residual variance and biases VPC downward.
    mu <- stats::fitted(model)
    mu <- pmax(as.numeric(mu), .Machine$double.eps)
    return(mean(log1p(1 / mu), na.rm = TRUE))
  }
  if (fam$family == "negbinomial" && fam$link == "log") {
    # Negative-binomial analogue of the Poisson approximation above: the
    # lognormal latent-scale level-1 variance ln(1 + 1/mu + 1/theta) of
    # Nakagawa, Johnson & Schielzeth (2017, J R Soc Interface 14:20170213).
    # The extra 1/theta term carries the overdispersion, so it reduces to the
    # Stryhn Poisson form as theta -> Inf.
    mu <- stats::fitted(model)
    mu <- pmax(as.numeric(mu), .Machine$double.eps)
    theta <- maihda_negbin_theta_lme4(model)
    return(mean(log1p(1 / mu + 1 / theta), na.rm = TRUE))
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
    # Stryhn et al. (2006) latent-scale level-1 variance approximation: log(1 + 1/mu).
    # See the lme4 sibling for the rationale; 1/mu alone biases VPC downward at low mu.
    mu <- stats::fitted(model, summary = TRUE)[, "Estimate"]
    mu <- pmax(as.numeric(mu), .Machine$double.eps)
    return(mean(log1p(1 / mu), na.rm = TRUE))
  }
  if (fam$family == "negbinomial" && fam$link == "log") {
    # Nakagawa, Johnson & Schielzeth (2017) lognormal latent-scale level-1
    # variance ln(1 + 1/mu + 1/theta); brms's 'shape' parameter IS theta.
    # Evaluated at the posterior-mean fitted means and posterior-mean shape
    # (see maihda_residual_variance_draws_brms for the per-draw treatment).
    mu <- stats::fitted(model, summary = TRUE)[, "Estimate"]
    mu <- pmax(as.numeric(mu), .Machine$double.eps)
    draws <- maihda_posterior_draws_brms(model)
    if (!"shape" %in% names(draws)) {
      stop("Could not extract the negative-binomial 'shape' (theta) draws from ",
           "the brms posterior.")
    }
    shape <- mean(as.numeric(draws[["shape"]]), na.rm = TRUE)
    return(mean(log1p(1 / mu + 1 / shape), na.rm = TRUE))
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
#   poisson(log)   -> mean(log1p(1 / mu)) at the POSTERIOR-MEAN fitted means
#                     (constant across draws). A per-draw version would need an
#                     ndraws x nobs posterior_epred() matrix, which is
#                     prohibitively expensive; the random-effect SD draws (the
#                     dominant source of VPC uncertainty) are still propagated --
#                     only the small level-1 term is held at its posterior mean.
#   negbinomial(log) -> mean(log1p(1 / mu + 1 / shape_d)) per draw d: the
#                     'shape' (= theta) draws are propagated, while mu is held
#                     at the posterior-mean fitted means for the same cost
#                     reason as the poisson case.
# Takes the model (for the family and, for poisson/negbinomial, the fitted
# means) and the draws data frame (for the sigma/shape draws and the draw count).
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
    mu <- stats::fitted(model, summary = TRUE)[, "Estimate"]
    mu <- pmax(as.numeric(mu), .Machine$double.eps)
    return(rep(mean(log1p(1 / mu), na.rm = TRUE), n))
  }
  if (fam$family == "negbinomial" && fam$link == "log") {
    # Nakagawa et al. (2017) ln(1 + 1/mu + 1/theta), theta = brms 'shape'.
    # The shape draws are propagated; mu is held at the posterior-mean fitted
    # means (see the header comment).
    mu <- stats::fitted(model, summary = TRUE)[, "Estimate"]
    mu <- pmax(as.numeric(mu), .Machine$double.eps)
    if (!"shape" %in% names(draws)) {
      stop("Could not extract the negative-binomial 'shape' (theta) draws from ",
           "the brms posterior.")
    }
    shape_d <- as.numeric(draws[["shape"]])
    return(vapply(shape_d,
                  function(s) mean(log1p(1 / mu + 1 / s), na.rm = TRUE),
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
  # Aggregated binomial responses (cbind(success, failure) / `y | trials(n)`)
  # report the binomial TRIALS through weights(type = "prior"); those are already
  # encoded in the response denominator, so multiplying by them would double-count.
  # Treat such models as unweighted for stratum-level aggregation.
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

maihda_stratum_predictions_lme4 <- function(object, summary_obj, scale = c("response", "link")) {
  scale <- match.arg(scale)
  data <- object$data
  if (!"stratum" %in% names(data)) {
    stop("'stratum' variable not found in fitted model data.")
  }

  model <- object$model
  fam <- maihda_family(model)
  linkinv <- maihda_linkinv(fam)
  prior_w <- maihda_prior_weights(object)

  eta_fixed <- stats::predict(model, newdata = data, re.form = NA, type = "link")
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
  # cross-classified model it must also carry the dimension random effects (the
  # additive part) so the stratum prediction includes ALL random effects, not just
  # the interaction: eta_base = eta(all REs) - u_stratum.
  eta_base <- eta_fixed
  if (!is.null(object$cc_info)) {
    eta_allre <- stats::predict(model, newdata = data, type = "link")
    u_row <- stratum_est$random_effect[idx]
    u_row[is.na(u_row)] <- 0
    eta_base <- as.numeric(eta_allre) - u_row
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

  # Prior-weight-weighted per-stratum means, so a weighted fit's stratum
  # predictions reflect the prior/precision weights (consistent with the weighted
  # VPC); identical to the previous unweighted means for an unweighted model. These
  # are lme4 prior weights, not a complex survey design (not survey-representative).
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

maihda_stratum_labels <- function(data, vars, sep = " \u00d7 ", autobin_info = NULL) {
  strata_data <- data[, vars, drop = FALSE]
  strata_data <- maihda_apply_autobin_info(strata_data, autobin_info)

  has_missing <- apply(strata_data, 1, function(x) any(is.na(x)))
  labels <- rep(NA_character_, nrow(strata_data))
  labels[!has_missing] <- apply(
    strata_data[!has_missing, , drop = FALSE],
    1,
    function(x) paste(x, collapse = sep)
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

maihda_prepare_prediction_data <- function(object, newdata) {
  if (!is.data.frame(newdata)) {
    stop("'newdata' must be a data frame", call. = FALSE)
  }
  if ("stratum" %in% names(newdata)) {
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

  unknown <- !is.na(labels) & is.na(newdata$stratum)
  if (any(unknown)) {
    unknown_labels <- unique(labels[unknown])
    stop("newdata contains strata combinations that were not present when the model was fit: ",
         paste(utils::head(unknown_labels, 5), collapse = ", "),
         if (length(unknown_labels) > 5) ", ..." else "",
         call. = FALSE)
  }

  newdata
}

maihda_stratum_predictions_brms <- function(object, summary_obj, scale = c("response", "link")) {
  scale <- match.arg(scale)
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
  prior_w <- maihda_prior_weights(object)
  eta_fixed <- brms::posterior_linpred(model, newdata = data, re_formula = NA, summary = TRUE)[, "Estimate"]

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

  # See the lme4 sibling: in a cross-classified model the stratum prediction must
  # also include the dimension random effects (the additive part). eta_base =
  # eta(all REs) - u_stratum, so adding the stratum RE recovers the full prediction.
  eta_base <- eta_fixed
  if (!is.null(object$cc_info)) {
    eta_allre <- brms::posterior_linpred(model, newdata = data,
                                         summary = TRUE)[, "Estimate"]
    u_row <- stratum_est$random_effect[idx]
    u_row[is.na(u_row)] <- 0
    eta_base <- as.numeric(eta_allre) - u_row
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

  # Prior-weight-weighted per-stratum means, so a weighted fit's stratum
  # predictions reflect the prior/precision weights (consistent with the weighted
  # VPC); identical to the previous unweighted means for an unweighted model. These
  # are lme4 prior weights, not a complex survey design (not survey-representative).
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
