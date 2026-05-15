# Internal helpers shared across model summaries, PVC, predictions, and plots.

maihda_family <- function(model) {
  fam <- tryCatch(stats::family(model), error = function(e) NULL)
  if (is.null(fam) && inherits(model, "brmsfit") && requireNamespace("brms", quietly = TRUE)) {
    fam <- brms::family(model)
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

maihda_stratum_variance_brms <- function(model, group = "stratum") {
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("Package 'brms' is required to work with brms models. Please install it with: install.packages('brms')")
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

  latent_families <- c("binomial", "quasibinomial", "cumulative", "sratio", "cratio", "acat", "ordinal")
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
    mu <- stats::fitted(model)
    mu <- pmax(as.numeric(mu), .Machine$double.eps)
    return(mean(1 / mu, na.rm = TRUE))
  }

  stop("VPC residual variance is not implemented for family '", fam$family,
       "' with link '", fam$link, "'.")
}

maihda_residual_variance_brms <- function(model) {
  fam <- maihda_family(model)
  if (is.null(fam)) {
    stop("Unable to determine brms model family for residual variance calculation.")
  }

  latent_families <- c("binomial", "quasibinomial", "cumulative", "sratio", "cratio", "acat", "ordinal")
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
    mu <- brms::fitted(model, summary = TRUE)[, "Estimate"]
    mu <- pmax(as.numeric(mu), .Machine$double.eps)
    return(mean(1 / mu, na.rm = TRUE))
  }

  stop("VPC residual variance is not implemented for brms family '", fam$family,
       "' with link '", fam$link, "'.")
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
