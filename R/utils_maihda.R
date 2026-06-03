# Internal helpers shared across model summaries, PVC, predictions, and plots.

maihda_family <- function(model) {
  fam <- tryCatch(stats::family(model), error = function(e) NULL)
  if (is.null(fam) && inherits(model, "brmsfit")) {
    fam <- tryCatch(model$family, error = function(e) NULL)
  }
  fam
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

maihda_validate_conf_level <- function(conf_level) {
  if (!is.numeric(conf_level) || length(conf_level) != 1 ||
      is.na(conf_level) || !is.finite(conf_level) ||
      conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be a single number between 0 and 1.", call. = FALSE)
  }

  as.numeric(conf_level)
}

maihda_validate_bootstrap_args <- function(n_boot, conf_level) {
  if (!is.numeric(n_boot) || length(n_boot) != 1 ||
      is.na(n_boot) || !is.finite(n_boot) ||
      n_boot < 1 || n_boot != floor(n_boot)) {
    stop("'n_boot' must be a single positive whole number.", call. = FALSE)
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

maihda_prepare_binomial_response <- function(data, formula) {
  response <- formula[[2]]
  if (!is.symbol(response)) {
    return(data)
  }

  outcome <- as.character(response)
  if (!outcome %in% names(data) || !maihda_is_binary_vector(data[[outcome]])) {
    return(data)
  }

  data[[outcome]] <- maihda_binary_to_01(data[[outcome]])
  data
}

# TRUE only when the model response is a single two-level (Bernoulli) vector,
# i.e. a plain symbol naming a binary column. Aggregated binomial responses such
# as cbind(success, failure) or `y | trials(n)` are calls, not symbols, so they
# return FALSE and must remain a binomial() model.
#
# The check is evaluated on the analytic sample (complete cases over the model
# variables), matching what lme4/brms actually fit: a response that is 0/1 once
# rows with missing covariates are dropped is still recognised as Bernoulli.
maihda_response_is_binary <- function(formula, data) {
  if (length(formula) != 3L) {
    return(FALSE)
  }
  response <- formula[[2]]
  if (!is.symbol(response)) {
    return(FALSE)
  }
  outcome <- as.character(response)
  if (!outcome %in% names(data)) {
    return(FALSE)
  }

  model_vars <- intersect(all.vars(formula), names(data))
  keep <- if (length(model_vars) > 0) {
    stats::complete.cases(data[, model_vars, drop = FALSE])
  } else {
    rep(TRUE, nrow(data))
  }
  maihda_is_binary_vector(data[[outcome]][keep])
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

  out <- rep(NA_integer_, nrow(data))
  for (i in seq_len(nrow(data))) {
    matches <- rep(TRUE, nrow(lookup))
    for (var in vars) {
      value <- data[[var]][i]
      matches <- matches & !is.na(value) & !is.na(lookup[[var]])
      matches <- matches & as.character(lookup[[var]]) == as.character(value)
    }
    idx <- which(matches)
    if (length(idx) > 0) {
      out[i] <- idx[1]
    }
  }

  out
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

maihda_residual_variance_lme4 <- function(model, vc = lme4::VarCorr(model)) {
  fam <- maihda_family(model)
  if (is.null(fam)) {
    stop("Unable to determine model family for residual variance calculation.")
  }

  latent_families <- c("binomial", "bernoulli", "quasibinomial", "cumulative", "sratio", "cratio", "acat", "ordinal")
  if (fam$family == "gaussian") {
    return(attr(vc, "sc")^2)
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

  stop("VPC residual variance is not implemented for family '", fam$family,
       "' with link '", fam$link, "'.")
}

maihda_residual_variance_brms <- function(model) {
  fam <- maihda_family(model)
  if (is.null(fam)) {
    stop("Unable to determine brms model family for residual variance calculation.")
  }

  latent_families <- c("binomial", "bernoulli", "quasibinomial", "cumulative", "sratio", "cratio", "acat", "ordinal")
  if (fam$family == "gaussian") {
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
# Takes the model (for the family and, for poisson, the fitted means) and the
# draws data frame (for the sigma draws and the draw count).
maihda_residual_variance_draws_brms <- function(model, draws) {
  fam <- maihda_family(model)
  if (is.null(fam)) {
    stop("Unable to determine brms model family for residual variance calculation.")
  }

  n <- nrow(draws)
  latent_families <- c("binomial", "bernoulli", "quasibinomial", "cumulative",
                       "sratio", "cratio", "acat", "ordinal")

  if (fam$family == "gaussian") {
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

maihda_stratum_predictions_lme4 <- function(object, summary_obj, scale = c("response", "link")) {
  scale <- match.arg(scale)
  data <- object$data
  if (!"stratum" %in% names(data)) {
    stop("'stratum' variable not found in fitted model data.")
  }

  model <- object$model
  fam <- maihda_family(model)
  linkinv <- maihda_linkinv(fam)

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

  pred_df <- data.frame(
    stratum = key,
    predicted_row = transform_eta(eta_fixed + stratum_est$random_effect[idx]),
    lower_row = transform_eta(eta_fixed + stratum_est$lower_95[idx]),
    upper_row = transform_eta(eta_fixed + stratum_est$upper_95[idx]),
    fixed_row = transform_eta(eta_fixed),
    stringsAsFactors = FALSE
  )

  out <- stats::aggregate(
    pred_df[, c("predicted_row", "lower_row", "upper_row", "fixed_row")],
    by = list(stratum = pred_df$stratum),
    FUN = mean,
    na.rm = TRUE
  )
  out$n <- as.integer(stats::aggregate(
    pred_df$predicted_row,
    by = list(stratum = pred_df$stratum),
    FUN = length
  )$x)
  out
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

  pred_df <- data.frame(
    stratum = key,
    predicted_row = transform_eta(eta_fixed + stratum_est$random_effect[idx]),
    lower_row = transform_eta(eta_fixed + stratum_est$lower_95[idx]),
    upper_row = transform_eta(eta_fixed + stratum_est$upper_95[idx]),
    fixed_row = transform_eta(eta_fixed),
    stringsAsFactors = FALSE
  )

  out <- stats::aggregate(
    pred_df[, c("predicted_row", "lower_row", "upper_row", "fixed_row")],
    by = list(stratum = pred_df$stratum),
    FUN = mean,
    na.rm = TRUE
  )
  out$n <- as.integer(stats::aggregate(
    pred_df$predicted_row,
    by = list(stratum = pred_df$stratum),
    FUN = length
  )$x)
  out
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
