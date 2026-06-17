#' Add Stratum Labels to Estimates
#'
#' Internal helper function to merge stratum labels into stratum estimates.
#'
#' @param stratum_estimates Data frame with stratum estimates
#' @param strata_info Data frame with stratum information including labels
#' @return Data frame with labels merged in
#' @keywords internal
add_stratum_labels <- function(stratum_estimates, strata_info) {
  if (is.null(strata_info) || !"stratum" %in% names(strata_info) || !"label" %in% names(strata_info)) {
    return(stratum_estimates)
  }

  idx <- match(as.character(stratum_estimates$stratum), as.character(strata_info$stratum))
  stratum_estimates$label <- strata_info$label[idx]

  col_order <- c("stratum", "stratum_id", "label", "random_effect", "se", "lower_95", "upper_95")
  stratum_estimates <- stratum_estimates[, col_order[col_order %in% names(stratum_estimates)]]

  return(stratum_estimates)
}

#' Summarize MAIHDA Model
#'
#' Provides a summary of a MAIHDA model including variance partition coefficients
#' (VPC/ICC) and stratum-specific estimates.
#'
#' @section Interpreting the VPC/ICC: The VPC is the between-stratum variance
#'   divided by the total \emph{unexplained} variance. For the canonical
#'   single-stratum model that denominator is between-stratum + residual, but if the
#'   model includes additional random effects (e.g. \code{(1 | site)}) their
#'   variance is included in the denominator too (between-stratum + other random
#'   effects + residual), so the VPC is the between-stratum \emph{share} of all
#'   unexplained variance. It is a conditional/residual ICC that excludes variance
#'   captured by the fixed effects, so for models with covariates it is conditional
#'   on them. It is most commonly read from the null model
#'   \code{outcome ~ 1 + (1 | stratum)}, where it is the total between-stratum
#'   share. For non-Gaussian families the level-1 (residual) variance uses a
#'   latent/distributional approximation (\eqn{\pi^2/3} for logistic,
#'   \eqn{\log(1 + 1/\mu)} for Poisson per Stryhn et al. 2006, and
#'   \eqn{\log(1 + 1/\mu + 1/\theta)} for the negative binomial per Nakagawa,
#'   Johnson & Schielzeth 2017), so the
#'   VPC is on that latent scale; for a \emph{weighted} Gaussian model the level-1
#'   variance is the mean conditional residual variance,
#'   \eqn{\bar{\sigma^2 / w_i}}, since the per-observation residual variance is
#'   \eqn{\sigma^2 / w_i}. The stratum random effects represent the total
#'   between-stratum deviation; they equal the \emph{pure} intersectional
#'   (interaction) component only when the additive main effects of the strata
#'   variables are included in the model.
#'
#' @param object A maihda_model object from \code{fit_maihda()}.
#' @param bootstrap Logical indicating whether to compute parametric bootstrap
#'   confidence intervals for VPC/ICC. Default is FALSE. Supported for lme4
#'   models only; \code{brms} models always return a posterior credible interval
#'   (see Details), so \code{bootstrap = TRUE} is rejected for them.
#'   For a negative-binomial model (\code{glmer.nb}) the bootstrap refits via
#'   \code{lme4::refit()}, which holds the dispersion parameter theta fixed at
#'   its original estimate, so the interval is conditional on the estimated
#'   theta (theta's own sampling uncertainty is not propagated). The
#'   \code{ordinal} (clmm) engine has no simulate/refit machinery, so
#'   \code{bootstrap = TRUE} is rejected there (use \code{engine = "brms"} for
#'   interval estimates).
#' @param n_boot Number of bootstrap samples if bootstrap = TRUE. Default is 1000.
#' @param conf_level Confidence level for the VPC/ICC interval -- the lme4
#'   bootstrap CI or the brms posterior credible interval. Default is 0.95.
#' @param response_vpc Logical; for a binomial (lme4) model, also compute the
#'   response-scale VPC (\code{\link{maihda_vpc_response}}) and attach it as the
#'   \code{vpc_response} slot. It is estimated by simulation, so it is opt-in (default
#'   \code{FALSE}) and uses \code{seed} for reproducibility. Ignored for other
#'   families/engines.
#' @param seed Optional integer seed for the response-scale VPC simulation when
#'   \code{response_vpc = TRUE}.
#' @param ... Additional arguments (not currently used).
#'
#' @return A maihda_summary object containing:
#'   \item{vpc}{Variance Partition Coefficient (ICC); for lme4 with
#'     \code{bootstrap = TRUE} and for all brms models this includes
#'     \code{ci_lower}/\code{ci_upper}/\code{conf_level}. For a contextual
#'     cross-classified fit this is the \emph{between-stratum} share of all
#'     unexplained variance (net of the context)}
#'   \item{variance_components}{Data frame of variance components. For a
#'     contextual cross-classified fit (\code{fit_maihda(context = )}) each
#'     context appears as its own \code{Context: <name>} row}
#'   \item{context}{For a contextual cross-classified fit, the stratum vs.
#'     context partition: per-context variances and shares, the contexts' total
#'     share (\code{vpc_context_total}, with an interval when bootstrapped or for
#'     brms), and the between-stratum share (\code{vpc_stratum}); \code{NULL}
#'     otherwise}
#'   \item{discriminatory_accuracy}{For a binomial/Bernoulli outcome, the
#'     \code{maihda_da} object (AUC + MOR) from
#'     \code{\link{maihda_discriminatory_accuracy}}; \code{NULL} otherwise (and for a
#'     cross-classified fit, whose single-stratum between-variance the MOR needs is
#'     not defined across crossed random effects)}
#'   \item{vpc_response}{The response-scale VPC (\code{maihda_vpc_response}) when
#'     \code{response_vpc = TRUE} for a binomial lme4 model; \code{NULL} otherwise}
#'   \item{stratum_estimates}{Data frame of stratum-specific random effects with labels if available}
#'   \item{fixed_effects}{Fixed effects estimates}
#'   \item{thresholds}{For a cumulative (ordinal) clmm fit, the threshold (cut
#'     point) estimates with standard errors -- the cumulative model's
#'     "intercepts"; NULL otherwise}
#'   \item{model_summary}{Original model summary}
#'   \item{diagnostics}{Fit-quality diagnostics (singular fit / convergence)
#'     carried over from the fitted model and reported by the print method}
#'
#' @note
#' For \code{lme4} models a VPC/ICC interval is obtained from a parametric
#' bootstrap (\code{bootstrap = TRUE}). For \code{brms} models the VPC/ICC is
#' summarised directly from the posterior draws: the reported estimate is the
#' posterior median of the per-draw VPC (\eqn{E[\sigma^2]}-based, not the biased
#' \eqn{E[\sigma]^2}) and the interval is a central credible interval at
#' \code{conf_level} (default 95\%), so no \code{bootstrap} argument is needed.
#' The variance-components table reports the posterior-mean variance components,
#' so the stratum proportion shown there may differ slightly from the headline
#' VPC because the median of a ratio is not the ratio of means. For non-Gaussian
#' \code{brms} families the level-1 (residual) variance uses the usual
#' latent-scale approximation; for \code{poisson(log)} it is evaluated at the
#' posterior-mean fitted values rather than per draw to avoid an expensive
#' \eqn{ndraws \times nobs} computation.
#'
#' @examples
#' \donttest{
#' strata_result <- make_strata(maihda_sim_data, vars = c("gender", "race"))
#' model <- fit_maihda(health_outcome ~ age + (1 | stratum), data = strata_result$data)
#' summary_result <- summary(model)
#'
#' # With bootstrap CI
#' # summary_boot <- summary(model, bootstrap = TRUE, n_boot = 50)
#' }
#'
#' @export
#' @importFrom lme4 VarCorr fixef ranef
summary.maihda_model <- function(object, bootstrap = FALSE, n_boot = 1000,
                          conf_level = 0.95, response_vpc = FALSE, seed = NULL, ...) {
  if (!inherits(object, "maihda_model")) {
    stop("'object' must be a maihda_model object from fit_maihda()")
  }

  if (!is.logical(bootstrap) || length(bootstrap) != 1 || is.na(bootstrap)) {
    stop("'bootstrap' must be TRUE or FALSE.", call. = FALSE)
  }
  if (bootstrap) {
    bootstrap_args <- maihda_validate_bootstrap_args(n_boot, conf_level)
    n_boot <- bootstrap_args$n_boot
    conf_level <- bootstrap_args$conf_level
  }

  engine <- object$engine
  model <- object$model
  # A crossed-dimensions model (tagged by maihda(decomposition =
  # "crossed-dimensions")) has several crossed REs: each dimension carries its
  # additive main-effect variance and the intersection ("stratum") RE the
  # interaction. A contextual cross-classified model (fit_maihda(context = )) has
  # the stratum RE crossed with one or more higher-level context REs, and the two
  # tags can co-occur. When neither is present, the variance path below is
  # identical to the historical single-stratum summary.
  cc <- object$cc_info
  ctx <- object$context_info
  # A longitudinal (3-level growth) fit (fit_maihda(id =, time =)) has random
  # slopes on time at the stratum and individual levels, so the between-stratum
  # variance -- and the VPC -- is time-varying. It routes to the longitudinal
  # path below, which reads the full random-effect covariance blocks instead of a
  # single intercept variance (and skips the intercept-only guard, kept for every
  # other model).
  lng <- object$longitudinal_info
  decomposition <- NULL
  context_summary <- NULL
  longitudinal <- NULL
  thresholds <- NULL

  # Extract variance components and calculate VPC
  if (engine == "lme4") {
    # Extract variance components
    vc <- lme4::VarCorr(model)
    if (!is.null(lng)) {
      lng_res <- maihda_longitudinal_summary_lme4(object, bootstrap, n_boot,
                                                  conf_level)
      variance_components <- lng_res$variance_components
      vpc_result <- lng_res$vpc_result
      longitudinal <- lng_res$longitudinal
    } else if (!is.null(cc)) {
      cc_res <- maihda_cc_summary_lme4(object, cc, vc, bootstrap, n_boot, conf_level)
      variance_components <- cc_res$variance_components
      vpc_result <- cc_res$vpc_result
      decomposition <- cc_res$decomposition
      context_summary <- cc_res$context
    } else if (!is.null(ctx)) {
      ctx_res <- maihda_context_summary_lme4(object, ctx, vc, bootstrap, n_boot,
                                             conf_level)
      variance_components <- ctx_res$variance_components
      vpc_result <- ctx_res$vpc_result
      context_summary <- ctx_res$context
    } else {
      var_random <- maihda_stratum_variance_lme4(model)
      var_total_random <- maihda_total_random_variance_lme4(model)
      var_other_random <- max(0, var_total_random - var_random)
      var_residual <- maihda_residual_variance_lme4(model, vc)

      # Calculate VPC (ICC)
      vpc <- var_random / (var_random + var_other_random + var_residual)

      # Create variance components data frame
      variance_components <- maihda_variance_components_table(
        var_random, var_other_random, var_residual
      )

      # Bootstrap confidence intervals for VPC if requested
      if (bootstrap) {
        vpc_ci <- bootstrap_vpc(model, object$data, object$formula, n_boot, conf_level)
        vpc_result <- list(
          estimate = vpc,
          ci_lower = vpc_ci[1],
          ci_upper = vpc_ci[2],
          conf_level = conf_level,
          bootstrap = TRUE,
          method = "bootstrap",
          n_boot_ok = attr(vpc_ci, "n_ok"),
          mc_se = attr(vpc_ci, "mc_se")
        )
      } else {
        vpc_result <- list(
          estimate = vpc,
          bootstrap = FALSE
        )
      }
    }

    # Extract fixed effects
    fixed_effects <- data.frame(
      term = names(lme4::fixef(model)),
      estimate = lme4::fixef(model),
      row.names = NULL
    )

    # Stratum (intersection) random-effect estimates -- the interaction residuals in
    # the cross-classified model (the named interaction group), or the single stratum
    # RE otherwise.
    interaction_group <- if (!is.null(cc)) cc$interaction_group else "stratum"
    stratum_estimates <- maihda_stratum_ranef_lme4(model, group = interaction_group)
    stratum_estimates <- add_stratum_labels(stratum_estimates, object$strata_info)

    # Get model summary
    model_summary <- summary(model)

  } else if (engine == "wemix") {
    if (bootstrap) {
      stop("Bootstrap VPC intervals are not available for the wemix engine: the ",
           "parametric bootstrap relies on lme4's simulate()/refit(), and a ",
           "design-based interval would require replicate weights (not yet ",
           "implemented). The design-weighted VPC is reported as a point ",
           "estimate; for interval estimates refit with engine = \"brms\" ",
           "(pseudo-posterior, with caveats).", call. = FALSE)
    }
    if (!is.null(cc) || !is.null(ctx)) {
      stop("Crossed-dimensions and contextual partitions are not available for ",
           "the wemix engine (WeMix fits no crossed random effects).",
           call. = FALSE)
    }

    # Canonical single-stratum partition from the pseudo-ML variance components.
    # For a binomial-logit fit the level-1 variance is the latent-scale pi^2/3,
    # exactly as in the lme4/brms summaries, so VPCs are comparable across engines.
    vars <- maihda_wemix_variances(object)
    vpc <- vars$stratum / (vars$stratum + vars$residual)
    variance_components <- maihda_variance_components_table(
      vars$stratum, 0, vars$residual
    )
    vpc_result <- list(estimate = vpc, bootstrap = FALSE)

    # WeMix fixed-effect standard errors are design-consistent (sandwich), the
    # main payoff of the design-weighted fit -- include them alongside the
    # estimates (the lme4 table reports estimates only).
    fixed_effects <- data.frame(
      term = names(object$model$coef),
      estimate = as.numeric(object$model$coef),
      se = as.numeric(object$model$SE[names(object$model$coef)]),
      row.names = NULL
    )

    stratum_estimates <- maihda_wemix_stratum_ranef(object)
    stratum_estimates <- add_stratum_labels(stratum_estimates, object$strata_info)

    model_summary <- tryCatch(summary(model), error = function(e) NULL)

  } else if (engine == "ordinal") {
    if (bootstrap) {
      stop("Bootstrap VPC intervals are not available for the ordinal engine: ",
           "the parametric bootstrap relies on lme4's simulate()/refit(), which ",
           "do not exist for ordinal::clmm fits. The VPC is reported as a point ",
           "estimate; for interval estimates refit with engine = \"brms\" ",
           "(posterior credible intervals).", call. = FALSE)
    }
    if (!is.null(cc) || !is.null(ctx)) {
      stop("Crossed-dimensions and contextual partitions are not available for ",
           "the ordinal engine (the clmm path fits the canonical single ",
           "(1 | stratum) structure only); use engine = \"brms\".", call. = FALSE)
    }

    # Canonical single-stratum partition on the latent scale: the level-1
    # variance is pi^2/3 (logit) or 1 (probit), the same latent treatment the
    # binomial summaries use, so cumulative VPCs are comparable to them.
    vars <- maihda_clmm_variances(object)
    vpc <- vars$stratum / (vars$stratum + vars$residual)
    variance_components <- maihda_variance_components_table(
      vars$stratum, 0, vars$residual
    )
    vpc_result <- list(estimate = vpc, bootstrap = FALSE)

    # Location coefficients with Hessian-based SEs; the thresholds (the
    # cumulative model's "intercepts") are reported separately below.
    beta <- object$model$beta
    if (is.null(beta) || length(beta) == 0) {
      fixed_effects <- data.frame(term = character(0), estimate = numeric(0),
                                  se = numeric(0))
    } else {
      V <- tryCatch(stats::vcov(object$model), error = function(e) NULL)
      beta_se <- rep(NA_real_, length(beta))
      if (!is.null(V) && all(names(beta) %in% rownames(V))) {
        beta_se <- sqrt(pmax(diag(V)[names(beta)], 0))
      }
      fixed_effects <- data.frame(
        term = names(beta),
        estimate = as.numeric(beta),
        se = as.numeric(beta_se),
        row.names = NULL
      )
    }

    stratum_estimates <- maihda_clmm_stratum_ranef(object)
    stratum_estimates <- add_stratum_labels(stratum_estimates, object$strata_info)

    thresholds <- tryCatch(maihda_clmm_thresholds(object), error = function(e) NULL)

    model_summary <- tryCatch(summary(model), error = function(e) NULL)

  } else if (engine == "brms") {
    if (bootstrap) {
      stop("Bootstrap VPC confidence intervals are only supported for lme4 models. ",
           "brms summaries already return a posterior credible interval for the ",
           "VPC/ICC, so 'bootstrap = TRUE' is not needed.",
           call. = FALSE)
    }

    # Verify brms is available
    if (!requireNamespace("brms", quietly = TRUE)) {
      stop("Package 'brms' is required to summarize brms models. Please install it with: install.packages('brms')")
    }

    conf_level <- maihda_validate_conf_level(conf_level)

    if (!is.null(lng)) {
      lng_res <- maihda_longitudinal_summary_brms(object, conf_level)
      variance_components <- lng_res$variance_components
      vpc_result <- lng_res$vpc_result
      longitudinal <- lng_res$longitudinal
    } else if (!is.null(cc)) {
      cc_res <- maihda_cc_summary_brms(object, cc, conf_level)
      variance_components <- cc_res$variance_components
      vpc_result <- cc_res$vpc_result
      decomposition <- cc_res$decomposition
      context_summary <- cc_res$context
    } else if (!is.null(ctx)) {
      ctx_res <- maihda_context_summary_brms(object, ctx, conf_level)
      variance_components <- ctx_res$variance_components
      vpc_result <- ctx_res$vpc_result
      context_summary <- ctx_res$context
    } else {
      # Summarise the VPC/ICC from posterior draws (E[sd^2], with a credible
      # interval) rather than from the posterior summary SDs (E[sd]^2, no interval).
      vpc_draws <- maihda_vpc_draws_brms(model, conf_level = conf_level)

      # Components table reports the posterior-mean variance of each component.
      variance_components <- maihda_variance_components_table(
        vpc_draws$var_stratum, vpc_draws$var_other_random, vpc_draws$var_residual
      )

      vpc_result <- list(
        estimate = vpc_draws$vpc$estimate,
        ci_lower = vpc_draws$vpc$ci_lower,
        ci_upper = vpc_draws$vpc$ci_upper,
        conf_level = conf_level,
        bootstrap = FALSE,
        method = "posterior"
      )
    }

    # Extract fixed effects
    fixed_effects <- brms::fixef(model)

    interaction_group <- if (!is.null(cc)) cc$interaction_group else "stratum"
    stratum_estimates <- maihda_stratum_ranef_brms(model, group = interaction_group)
    stratum_estimates <- add_stratum_labels(stratum_estimates, object$strata_info)

    model_summary <- summary(model)
  }

  # Family-defined companions to the VPC for a binomial/Bernoulli outcome: the
  # discriminatory accuracy (AUC + MOR) -- the "DA" in MAIHDA -- always, and the
  # response-scale VPC on request (it is a simulation, hence opt-in and seeded).
  # Both summarise the single fitted model with no refit, mirroring how the
  # cross-classified additive/interaction `decomposition` slot above is computed in
  # this same summary layer so that fit_maihda() and maihda() share the logic.
  # Skipped for a cross-classified fit (its MOR / response-VPC need a single-stratum
  # between-variance that is not defined across crossed random effects) and wrapped so
  # a bonus summary never breaks the core VPC.
  discriminatory_accuracy <- NULL
  vpc_response <- NULL
  fam_name <- tryCatch(maihda_model_family_name(object),
                       error = function(e) NA_character_)
  if (is.null(cc) && is.null(lng) && isTRUE(fam_name %in% c("binomial", "bernoulli"))) {
    discriminatory_accuracy <- tryCatch(
      maihda_discriminatory_accuracy(object), error = function(e) NULL)
    if (isTRUE(response_vpc) && identical(engine, "lme4") &&
        identical(fam_name, "binomial")) {
      vpc_response <- tryCatch(
        maihda_vpc_response(object, seed = seed), error = function(e) NULL)
    }
  }

  # Create summary object
  result <- structure(
    list(
      vpc = vpc_result,
      variance_components = variance_components,
      decomposition = decomposition,
      context = context_summary,
      longitudinal = longitudinal,
      discriminatory_accuracy = discriminatory_accuracy,
      vpc_response = vpc_response,
      stratum_estimates = stratum_estimates,
      fixed_effects = fixed_effects,
      thresholds = thresholds,
      model_summary = model_summary,
      engine = engine,
      cc_info = cc,
      context_info = ctx,
      longitudinal_info = lng,
      diagnostics = object$diagnostics
    ),
    class = "maihda_summary"
  )

  return(result)
}

#' Crossed-dimensions variance summary (lme4)
#'
#' Internal helper for \code{\link{summary.maihda_model}} when the model is a
#' crossed-dimensions MAIHDA fit (\code{object$cc_info} set). Partitions the crossed
#' random-effect variances into the additive (sum of the dimension REs) and
#' interaction (intersection RE) components, builds the variance-components table and
#' the VPC, and -- when \code{bootstrap = TRUE} -- adds parametric-bootstrap intervals
#' for the VPC and the additive/interaction shares. When the fit also carries a
#' contextual random intercept (\code{object$context_info} set), the context
#' variance enters the VPC denominator and is reported as its own component row
#' and \code{context} element.
#'
#' @param object A \code{maihda_model} (crossed-dimensions).
#' @param cc The \code{cc_info} list (\code{dim_groups}, \code{interaction_group}).
#' @param vc The model's \code{VarCorr}.
#' @param bootstrap,n_boot,conf_level Bootstrap controls.
#' @return A list with \code{variance_components}, \code{vpc_result},
#'   \code{decomposition} and \code{context} (NULL without a context).
#' @keywords internal
maihda_cc_summary_lme4 <- function(object, cc, vc, bootstrap, n_boot, conf_level) {
  model <- object$model
  var_named <- maihda_random_variances_lme4(model)
  var_within <- maihda_residual_variance_lme4(model, vc)
  split <- maihda_cc_variance_split(var_named, cc$dim_groups, cc$interaction_group)
  ctx_vars <- if (!is.null(object$context_info)) {
    object$context_info$context_vars
  } else {
    character(0)
  }
  per_context <- if (length(ctx_vars) > 0) var_named[ctx_vars] else NULL
  var_context_total <- if (is.null(per_context)) 0 else sum(per_context)
  part <- maihda_cc_partition(split$additive, split$interaction, var_within,
                              var_context_total)

  variance_components <- maihda_cc_components_table(split$per_dim, split$interaction,
                                                   var_within, per_context)

  decomposition <- list(
    additive_var = split$additive,
    interaction_var = split$interaction,
    between_var = part$between,
    within_var = var_within,
    additive_share = part$additive_share,
    interaction_share = part$interaction_share,
    per_dim = split$per_dim,
    bootstrap = FALSE
  )

  context_summary <- NULL
  if (length(ctx_vars) > 0) {
    context_summary <- list(
      context_vars = ctx_vars,
      var_stratum = part$between,
      per_context = per_context,
      context_var_total = var_context_total,
      within_var = var_within,
      other_var = 0,
      vpc_stratum = part$vpc,
      vpc_context = per_context / part$total,
      vpc_context_total = var_context_total / part$total,
      bootstrap = FALSE
    )
  }

  if (bootstrap) {
    boot <- bootstrap_cc(model, cc, n_boot, conf_level, ctx_vars = ctx_vars)
    vpc_result <- list(
      estimate = part$vpc,
      ci_lower = boot$vpc[1],
      ci_upper = boot$vpc[2],
      conf_level = conf_level,
      bootstrap = TRUE,
      method = "bootstrap",
      n_boot_ok = attr(boot$vpc, "n_ok"),
      mc_se = attr(boot$vpc, "mc_se")
    )
    decomposition$bootstrap <- TRUE
    decomposition$conf_level <- conf_level
    decomposition$additive_share_ci <- c(boot$additive_share[1], boot$additive_share[2])
    decomposition$interaction_share_ci <- c(boot$interaction_share[1],
                                            boot$interaction_share[2])
    if (!is.null(context_summary)) {
      context_summary$bootstrap <- TRUE
      context_summary$conf_level <- conf_level
      context_summary$vpc_context_total_ci <- c(boot$context_vpc[1],
                                                boot$context_vpc[2])
    }
  } else {
    vpc_result <- list(estimate = part$vpc, bootstrap = FALSE)
  }

  list(variance_components = variance_components, vpc_result = vpc_result,
       decomposition = decomposition, context = context_summary)
}

#' Crossed-dimensions variance summary (brms)
#'
#' brms counterpart of \code{\link{maihda_cc_summary_lme4}}: computes the additive /
#' interaction partition per posterior draw and returns posterior point estimates with
#' credible intervals for the VPC and the shares (no bootstrap -- the posterior already
#' supplies the interval). A contextual random intercept
#' (\code{object$context_info} set) enters the per-draw VPC denominator and is
#' reported as its own component row and \code{context} element.
#'
#' @param object A \code{maihda_model} (crossed-dimensions, brms engine).
#' @param cc The \code{cc_info} list.
#' @param conf_level Credible-interval level.
#' @param point Posterior point estimate, "median" (default) or "mean".
#' @return A list with \code{variance_components}, \code{vpc_result},
#'   \code{decomposition} and \code{context} (NULL without a context).
#' @keywords internal
maihda_cc_summary_brms <- function(object, cc, conf_level, point = c("median", "mean")) {
  point <- match.arg(point)
  model <- object$model
  draws <- maihda_posterior_draws_brms(model)
  gv <- maihda_group_variance_draws_brms(draws)
  within_draws <- maihda_residual_variance_draws_brms(model, draws)

  dim_re <- unname(cc$dim_groups)
  ctx_vars <- if (!is.null(object$context_info)) {
    object$context_info$context_vars
  } else {
    character(0)
  }
  missing_re <- setdiff(c(dim_re, cc$interaction_group, ctx_vars), names(gv))
  if (length(missing_re) > 0) {
    stop("Crossed-dimensions brms summary is missing the random effect(s): ",
         paste(missing_re, collapse = ", "), ".", call. = FALSE)
  }

  additive_draws <- Reduce(`+`, gv[dim_re])
  interaction_draws <- gv[[cc$interaction_group]]
  context_total_draws <- if (length(ctx_vars) > 0) Reduce(`+`, gv[ctx_vars]) else 0
  part <- maihda_cc_partition(additive_draws, interaction_draws, within_draws,
                              context_total_draws)

  summ <- function(v) {
    v <- v[is.finite(v)]
    if (length(v) == 0) {
      return(list(estimate = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_))
    }
    a <- 1 - conf_level
    pt <- if (point == "median") stats::median(v) else mean(v)
    ci <- stats::quantile(v, probs = c(a / 2, 1 - a / 2), names = FALSE)
    list(estimate = pt, ci_lower = ci[1], ci_upper = ci[2])
  }
  vpc_s <- summ(part$vpc)
  add_s <- summ(part$additive_share)
  int_s <- summ(part$interaction_share)

  per_dim_mean <- vapply(dim_re, function(g) mean(gv[[g]]), numeric(1))
  names(per_dim_mean) <- names(cc$dim_groups)
  interaction_mean <- mean(interaction_draws)
  within_mean <- mean(within_draws)
  per_context_mean <- if (length(ctx_vars) > 0) {
    vapply(ctx_vars, function(g) mean(gv[[g]]), numeric(1))
  } else {
    NULL
  }

  variance_components <- maihda_cc_components_table(per_dim_mean, interaction_mean,
                                                   within_mean, per_context_mean)

  vpc_result <- list(
    estimate = vpc_s$estimate,
    ci_lower = vpc_s$ci_lower,
    ci_upper = vpc_s$ci_upper,
    conf_level = conf_level,
    bootstrap = FALSE,
    method = "posterior"
  )

  decomposition <- list(
    additive_var = mean(additive_draws),
    interaction_var = interaction_mean,
    between_var = mean(part$between),
    within_var = within_mean,
    additive_share = add_s$estimate,
    interaction_share = int_s$estimate,
    additive_share_ci = c(add_s$ci_lower, add_s$ci_upper),
    interaction_share_ci = c(int_s$ci_lower, int_s$ci_upper),
    per_dim = per_dim_mean,
    bootstrap = FALSE,
    method = "posterior",
    conf_level = conf_level
  )

  context_summary <- NULL
  if (length(ctx_vars) > 0) {
    ctx_total_s <- summ(part$other / part$total)
    context_summary <- list(
      context_vars = ctx_vars,
      var_stratum = mean(part$between),
      per_context = per_context_mean,
      context_var_total = mean(context_total_draws),
      within_var = within_mean,
      other_var = 0,
      vpc_stratum = vpc_s$estimate,
      vpc_context = vapply(ctx_vars, function(g) summ(gv[[g]] / part$total)$estimate,
                           numeric(1)),
      vpc_context_total = ctx_total_s$estimate,
      vpc_context_total_ci = c(ctx_total_s$ci_lower, ctx_total_s$ci_upper),
      bootstrap = FALSE,
      method = "posterior",
      conf_level = conf_level
    )
  }

  list(variance_components = variance_components, vpc_result = vpc_result,
       decomposition = decomposition, context = context_summary)
}

#' Bootstrap a crossed-dimensions MAIHDA partition (lme4)
#'
#' Parametric bootstrap (simulate from the fitted model, refit) of the
#' crossed-dimensions VPC and the additive / interaction shares, returning a
#' percentile interval for each via \code{maihda_bootstrap_ci}. lme4 only -- brms
#' returns posterior credible intervals directly. When \code{ctx_vars} names
#' contextual random intercepts, their variance enters each refit's VPC denominator
#' and a \code{context_vpc} interval (the contexts' total share) is returned too.
#'
#' @param model The underlying lme4 model object.
#' @param cc The \code{cc_info} list.
#' @param n_boot Number of bootstrap samples.
#' @param conf_level Confidence level.
#' @param ctx_vars Character vector of contextual grouping factors (may be empty).
#' @return A list with \code{vpc}, \code{additive_share}, \code{interaction_share}
#'   (and \code{context_vpc} when \code{ctx_vars} is non-empty), each a length-2
#'   interval carrying \code{n_ok}/\code{mc_se} attributes.
#' @keywords internal
#' @importFrom lme4 refit
bootstrap_cc <- function(model, cc, n_boot, conf_level, ctx_vars = character(0)) {
  vpc_boot <- rep(NA_real_, n_boot)
  additive_boot <- rep(NA_real_, n_boot)
  interaction_boot <- rep(NA_real_, n_boot)
  context_boot <- rep(NA_real_, n_boot)
  sim_data <- stats::simulate(model, nsim = n_boot)

  for (i in seq_len(n_boot)) {
    tryCatch({
      boot_model <- lme4::refit(model, newresp = sim_data[[i]])
      var_named <- maihda_random_variances_lme4(boot_model)
      var_within <- maihda_residual_variance_lme4(boot_model)
      split <- maihda_cc_variance_split(var_named, cc$dim_groups, cc$interaction_group)
      var_context <- if (length(ctx_vars) > 0) sum(var_named[ctx_vars]) else 0
      part <- maihda_cc_partition(split$additive, split$interaction, var_within,
                                  var_context)
      vpc_boot[i] <- part$vpc
      additive_boot[i] <- part$additive_share
      interaction_boot[i] <- part$interaction_share
      context_boot[i] <- var_context / part$total
    }, error = function(e) NULL)
  }

  out <- list(
    vpc = maihda_bootstrap_ci(vpc_boot, n_boot, conf_level, "VPC"),
    additive_share = maihda_bootstrap_ci(additive_boot, n_boot, conf_level,
                                         "additive share"),
    interaction_share = maihda_bootstrap_ci(interaction_boot, n_boot, conf_level,
                                            "interaction share")
  )
  if (length(ctx_vars) > 0) {
    out$context_vpc <- maihda_bootstrap_ci(context_boot, n_boot, conf_level,
                                           "context VPC")
  }
  out
}

#' Contextual cross-classified variance summary (lme4)
#'
#' Internal helper for \code{\link{summary.maihda_model}} when the model carries a
#' contextual random intercept (\code{object$context_info} set,
#' \code{fit_maihda(context = )}) without the crossed-dimensions decomposition.
#' Partitions the unexplained variance into between-stratum vs. between-context
#' (one share per context variable) vs. residual. The headline VPC stays the
#' between-stratum share of all unexplained variance -- numerically identical to
#' the generic single-stratum summary, which folds the context into "Other random
#' effects" -- but the context is now named, given its own component row(s), and
#' returned as a \code{context} element.
#'
#' @param object A \code{maihda_model} with \code{context_info}.
#' @param ctx The \code{context_info} list (\code{context_vars}).
#' @param vc The model's \code{VarCorr}.
#' @param bootstrap,n_boot,conf_level Bootstrap controls.
#' @return A list with \code{variance_components}, \code{vpc_result}, \code{context}.
#' @keywords internal
maihda_context_summary_lme4 <- function(object, ctx, vc, bootstrap, n_boot,
                                        conf_level) {
  model <- object$model
  var_named <- maihda_random_variances_lme4(model)
  var_within <- maihda_residual_variance_lme4(model, vc)
  ctx_vars <- ctx$context_vars
  missing_re <- setdiff(c("stratum", ctx_vars), names(var_named))
  if (length(missing_re) > 0) {
    stop("Contextual variance partition is missing the random effect(s): ",
         paste(missing_re, collapse = ", "),
         ". Expected the stratum intercept plus one intercept per context.",
         call. = FALSE)
  }
  var_stratum <- unname(var_named[["stratum"]])
  per_context <- var_named[ctx_vars]
  # Any further random effects beyond stratum + context (rare; e.g. a manual
  # extra grouping) stay in the denominator as "Other random effects".
  var_other <- max(0, sum(var_named, na.rm = TRUE) - var_stratum - sum(per_context))
  part <- maihda_context_partition(var_stratum, as.list(per_context), var_within,
                                   var_other)

  variance_components <- maihda_context_components_table(var_stratum, per_context,
                                                         var_other, var_within)

  context_summary <- list(
    context_vars = ctx_vars,
    var_stratum = var_stratum,
    per_context = per_context,
    context_var_total = part$context_total,
    within_var = var_within,
    other_var = var_other,
    vpc_stratum = part$vpc_stratum,
    vpc_context = unlist(part$vpc_context),
    vpc_context_total = part$vpc_context_total,
    bootstrap = FALSE
  )

  if (bootstrap) {
    boot <- bootstrap_context(model, ctx_vars, n_boot, conf_level)
    vpc_result <- list(
      estimate = part$vpc_stratum,
      ci_lower = boot$vpc[1],
      ci_upper = boot$vpc[2],
      conf_level = conf_level,
      bootstrap = TRUE,
      method = "bootstrap",
      n_boot_ok = attr(boot$vpc, "n_ok"),
      mc_se = attr(boot$vpc, "mc_se")
    )
    context_summary$bootstrap <- TRUE
    context_summary$conf_level <- conf_level
    context_summary$vpc_context_total_ci <- c(boot$context_vpc[1],
                                              boot$context_vpc[2])
  } else {
    vpc_result <- list(estimate = part$vpc_stratum, bootstrap = FALSE)
  }

  list(variance_components = variance_components, vpc_result = vpc_result,
       context = context_summary)
}

#' Contextual cross-classified variance summary (brms)
#'
#' brms counterpart of \code{\link{maihda_context_summary_lme4}}: computes the
#' stratum / context / residual partition per posterior draw and returns posterior
#' point estimates with credible intervals for the between-stratum VPC and the
#' contexts' total share (no bootstrap -- the posterior supplies the interval).
#'
#' @param object A \code{maihda_model} with \code{context_info} (brms engine).
#' @param ctx The \code{context_info} list.
#' @param conf_level Credible-interval level.
#' @param point Posterior point estimate, "median" (default) or "mean".
#' @return A list with \code{variance_components}, \code{vpc_result}, \code{context}.
#' @keywords internal
maihda_context_summary_brms <- function(object, ctx, conf_level,
                                        point = c("median", "mean")) {
  point <- match.arg(point)
  model <- object$model
  draws <- maihda_posterior_draws_brms(model)
  gv <- maihda_group_variance_draws_brms(draws)
  within_draws <- maihda_residual_variance_draws_brms(model, draws)

  ctx_vars <- ctx$context_vars
  missing_re <- setdiff(c("stratum", ctx_vars), names(gv))
  if (length(missing_re) > 0) {
    stop("Contextual brms summary is missing the random effect(s): ",
         paste(missing_re, collapse = ", "), ".", call. = FALSE)
  }

  stratum_draws <- gv[["stratum"]]
  context_draws <- gv[ctx_vars]
  other_groups <- setdiff(names(gv), c("stratum", ctx_vars))
  other_draws <- if (length(other_groups) > 0) Reduce(`+`, gv[other_groups]) else 0
  part <- maihda_context_partition(stratum_draws, context_draws, within_draws,
                                   other_draws)

  summ <- function(v) {
    v <- v[is.finite(v)]
    if (length(v) == 0) {
      return(list(estimate = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_))
    }
    a <- 1 - conf_level
    pt <- if (point == "median") stats::median(v) else mean(v)
    ci <- stats::quantile(v, probs = c(a / 2, 1 - a / 2), names = FALSE)
    list(estimate = pt, ci_lower = ci[1], ci_upper = ci[2])
  }
  vpc_s <- summ(part$vpc_stratum)
  ctx_total_s <- summ(part$vpc_context_total)

  per_context_mean <- vapply(ctx_vars, function(g) mean(gv[[g]]), numeric(1))
  other_mean <- if (is.numeric(other_draws) && length(other_draws) > 1) {
    mean(other_draws)
  } else {
    other_draws
  }
  variance_components <- maihda_context_components_table(
    mean(stratum_draws), per_context_mean, other_mean, mean(within_draws)
  )

  vpc_result <- list(
    estimate = vpc_s$estimate,
    ci_lower = vpc_s$ci_lower,
    ci_upper = vpc_s$ci_upper,
    conf_level = conf_level,
    bootstrap = FALSE,
    method = "posterior"
  )

  context_summary <- list(
    context_vars = ctx_vars,
    var_stratum = mean(stratum_draws),
    per_context = per_context_mean,
    context_var_total = mean(part$context_total),
    within_var = mean(within_draws),
    other_var = other_mean,
    vpc_stratum = vpc_s$estimate,
    vpc_context = vapply(ctx_vars, function(g) summ(gv[[g]] / part$total)$estimate,
                         numeric(1)),
    vpc_context_total = ctx_total_s$estimate,
    vpc_context_total_ci = c(ctx_total_s$ci_lower, ctx_total_s$ci_upper),
    bootstrap = FALSE,
    method = "posterior",
    conf_level = conf_level
  )

  list(variance_components = variance_components, vpc_result = vpc_result,
       context = context_summary)
}

#' Bootstrap a contextual cross-classified MAIHDA partition (lme4)
#'
#' Parametric bootstrap (simulate from the fitted model, refit) of the
#' between-stratum VPC and the contexts' total share for a contextual
#' cross-classified fit, returning a percentile interval for each via
#' \code{maihda_bootstrap_ci}. lme4 only -- brms returns posterior credible
#' intervals directly.
#'
#' @param model The underlying lme4 model object.
#' @param ctx_vars Character vector of context grouping factors.
#' @param n_boot Number of bootstrap samples.
#' @param conf_level Confidence level.
#' @return A list with \code{vpc} (between-stratum share) and \code{context_vpc}
#'   (contexts' total share), each a length-2 interval carrying
#'   \code{n_ok}/\code{mc_se} attributes.
#' @keywords internal
#' @importFrom lme4 refit
bootstrap_context <- function(model, ctx_vars, n_boot, conf_level) {
  vpc_boot <- rep(NA_real_, n_boot)
  context_boot <- rep(NA_real_, n_boot)
  sim_data <- stats::simulate(model, nsim = n_boot)

  for (i in seq_len(n_boot)) {
    tryCatch({
      boot_model <- lme4::refit(model, newresp = sim_data[[i]])
      var_named <- maihda_random_variances_lme4(boot_model)
      var_within <- maihda_residual_variance_lme4(boot_model)
      var_stratum <- unname(var_named[["stratum"]])
      per_context <- var_named[ctx_vars]
      var_other <- max(0, sum(var_named, na.rm = TRUE) - var_stratum -
                         sum(per_context))
      part <- maihda_context_partition(var_stratum, as.list(per_context),
                                       var_within, var_other)
      vpc_boot[i] <- part$vpc_stratum
      context_boot[i] <- part$vpc_context_total
    }, error = function(e) NULL)
  }

  list(
    vpc = maihda_bootstrap_ci(vpc_boot, n_boot, conf_level, "VPC"),
    context_vpc = maihda_bootstrap_ci(context_boot, n_boot, conf_level,
                                      "context VPC")
  )
}

#' Bootstrap VPC/ICC
#'
#' Internal function to compute bootstrap confidence intervals for VPC.
#'
#' @param model An lme4 model object
#' @param data The data used to fit the model
#' @param formula The model formula
#' @param n_boot Number of bootstrap samples
#' @param conf_level Confidence level
#'
#' @return A vector with lower and upper confidence bounds
#' @keywords internal
#' @importFrom lme4 lmer glmer VarCorr
bootstrap_vpc <- function(model, data, formula, n_boot, conf_level) {
  # Initialise to NA so iterations whose refit() throws — and never reach the
  # assignment inside the tryCatch body — stay NA rather than the numeric() default of 0.
  # The error handler runs in its own scope and cannot write back to this vector,
  # so the initial value is what survives a failure.
  vpc_boot <- rep(NA_real_, n_boot)
  sim_data <- stats::simulate(model, nsim = n_boot)

  for (i in 1:n_boot) {
    tryCatch({
      boot_model <- lme4::refit(model, newresp = sim_data[[i]])

      # Calculate VPC
      vc <- lme4::VarCorr(boot_model)
      var_random <- maihda_stratum_variance_lme4(boot_model)
      var_total_random <- maihda_total_random_variance_lme4(boot_model)
      var_other_random <- max(0, var_total_random - var_random)
      var_residual <- maihda_residual_variance_lme4(boot_model, vc)

      vpc_boot[i] <- var_random / (var_random + var_other_random + var_residual)
    }, error = function(e) NULL)
  }

  # Reduce to an interval, requiring a minimum number of successful refits.
  ci <- maihda_bootstrap_ci(vpc_boot, n_boot, conf_level, "VPC")

  return(ci)
}

#' Print the additive vs. intersectional decomposition of a crossed-dimensions summary
#'
#' @param d The \code{decomposition} list from a crossed-dimensions
#'   \code{\link{summary.maihda_model}}.
#' @return No return value, called for side effects.
#' @keywords internal
maihda_print_cc_decomposition <- function(d) {
  fmt_share <- function(est, ci) {
    if (!is.null(ci) && length(ci) == 2L && all(is.finite(ci))) {
      sprintf("%.1f%% [%.1f%%, %.1f%%]", est * 100, ci[1] * 100, ci[2] * 100)
    } else {
      sprintf("%.1f%%", est * 100)
    }
  }
  cat("Additive vs. Intersectional Decomposition (crossed-dimensions):\n")
  cat(sprintf("  Additive (sum of dimension main effects) variance: %.4f\n",
              d$additive_var))
  cat(sprintf("  Intersectional interaction variance:               %.4f\n",
              d$interaction_var))
  cat(sprintf("  Total between-strata variance:                     %.4f\n",
              d$between_var))
  cat(sprintf("  Additive share of between-strata variance:    %s\n",
              fmt_share(d$additive_share, d$additive_share_ci)))
  cat(sprintf("  Interaction share of between-strata variance: %s\n",
              fmt_share(d$interaction_share, d$interaction_share_ci)))
  per_dim <- d$per_dim
  if (!is.null(per_dim) && length(per_dim) > 0) {
    cat("  Per-dimension additive variance:\n")
    for (nm in names(per_dim)) {
      cat(sprintf("    %s: %.4f\n", nm, per_dim[[nm]]))
    }
  }
  cat("  Note: the additive share is the crossed-dimensions analogue of the PCV but\n",
      "  a different estimator; interpret the interaction share cautiously.\n\n", sep = "")
  invisible(NULL)
}

#' Print the stratum vs. context partition of a contextual cross-classified summary
#'
#' @param ctx The \code{context} list from a contextual
#'   \code{\link{summary.maihda_model}} (a model fitted with
#'   \code{fit_maihda(context = )}).
#' @return No return value, called for side effects.
#' @keywords internal
maihda_print_context_partition <- function(ctx) {
  fmt_share <- function(est, ci = NULL) {
    if (!is.null(ci) && length(ci) == 2L && all(is.finite(ci))) {
      sprintf("%.1f%% [%.1f%%, %.1f%%]", est * 100, ci[1] * 100, ci[2] * 100)
    } else {
      sprintf("%.1f%%", est * 100)
    }
  }
  cat("Contextual Cross-Classified Partition (stratum x context):\n")
  cat(sprintf("  Between-stratum (intersectional) variance: %.4f (share %s)\n",
              ctx$var_stratum, fmt_share(ctx$vpc_stratum)))
  per_context <- ctx$per_context
  vpc_context <- ctx$vpc_context
  for (nm in names(per_context)) {
    cat(sprintf("  Context '%s' variance: %.4f (share %s)\n",
                nm, per_context[[nm]], fmt_share(vpc_context[[nm]])))
  }
  if (length(per_context) > 1) {
    cat(sprintf("  All contexts combined: %.4f (share %s)\n",
                ctx$context_var_total,
                fmt_share(ctx$vpc_context_total, ctx$vpc_context_total_ci)))
  } else if (!is.null(ctx$vpc_context_total_ci)) {
    cat(sprintf("  Context share interval: %s\n",
                fmt_share(ctx$vpc_context_total, ctx$vpc_context_total_ci)))
  }
  cat("  Note: the headline VPC/ICC is the between-stratum share net of the\n",
      "  context(s) -- intersectional clustering not attributable to shared place\n",
      "  or institution. The context share is the general contextual effect.\n\n",
      sep = "")
  invisible(NULL)
}

#' Print method for maihda_summary objects
#'
#' @param x A maihda_summary object
#' @param ... Additional arguments (not used)
#' @return No return value, called for side effects.
#' @export
print.maihda_summary <- function(x, ...) {
  cat("MAIHDA Model Summary\n")
  cat("====================\n\n")

  maihda_print_fit_diagnostics(x$diagnostics)

  is_lng <- !is.null(x$longitudinal)
  if (is_lng) {
    cat(sprintf("Variance Partition Coefficient (VPC/ICC) at baseline (%s = %g):\n",
                x$longitudinal$time, x$longitudinal$ref_time))
  } else {
    cat("Variance Partition Coefficient (VPC/ICC):\n")
  }
  if (maihda_vpc_has_interval(x$vpc)) {
    cat(sprintf("  Estimate: %.4f [%.4f, %.4f]\n",
                x$vpc$estimate, x$vpc$ci_lower, x$vpc$ci_upper))
    cat("  ", maihda_vpc_interval_label(x$vpc), "\n", sep = "")
    if (!is.null(x$vpc$mc_se) && is.finite(x$vpc$mc_se)) {
      cat(sprintf("  (%d successful bootstrap draws; Monte Carlo SE %.4f)\n",
                  as.integer(x$vpc$n_boot_ok), x$vpc$mc_se))
    }
    cat("\n")
  } else {
    cat(sprintf("  Estimate: %.4f\n\n", x$vpc$estimate))
  }

  cat("Variance Components:\n")
  print(x$variance_components, row.names = FALSE, digits = 4)
  cat("\n")

  if (is_lng) {
    vt <- x$longitudinal$vpc_t
    cat(sprintf("Time-varying VPC/ICC (between-stratum share over %s):\n",
                x$longitudinal$time))
    cat(sprintf("  range %.4f to %.4f across %s in [%g, %g].\n",
                min(vt$estimate, na.rm = TRUE), max(vt$estimate, na.rm = TRUE),
                x$longitudinal$time, min(vt$time), max(vt$time)))
    cat("  The between-stratum variance is a function of time (random intercept +\n")
    cat("  slope), so the VPC varies; it depends on where time is zeroed. See\n")
    cat("  plot(type = \"vpc_trajectory\") for the full curve.\n\n")
  }

  if (!is.null(x$decomposition)) {
    maihda_print_cc_decomposition(x$decomposition)
  }

  if (!is.null(x$context)) {
    maihda_print_context_partition(x$context)
  }

  # Discriminatory accuracy (AUC + MOR) and, when requested, the response-scale VPC --
  # the binomial companions to the latent-scale VPC. Reuse their own print methods.
  if (!is.null(x$discriminatory_accuracy)) {
    print(x$discriminatory_accuracy)
    cat("\n")
  }
  if (!is.null(x$vpc_response)) {
    print(x$vpc_response)
    cat("\n")
  }

  cat("Fixed Effects:\n")
  print(x$fixed_effects, row.names = FALSE, digits = 4)
  cat("\n")

  if (!is.null(x$thresholds) && nrow(x$thresholds) > 0) {
    cat("Thresholds (cumulative cut points; they take the intercept's place):\n")
    print(x$thresholds, row.names = FALSE, digits = 4)
    cat("\n")
  }

  if (!is.null(x$stratum_estimates) && nrow(x$stratum_estimates) > 0) {
    # For a longitudinal fit the stratum random_effect is the baseline (intercept)
    # deviation only -- the random slope is not shown here -- so label it as such
    # and point to the trajectory tools rather than letting it read as a single
    # cross-sectional stratum effect.
    if (is_lng) {
      cat("Stratum baseline (intercept) deviations (first 10):\n")
    } else {
      cat("Stratum Estimates (first 10):\n")
    }
    print(utils::head(x$stratum_estimates, 10), row.names = FALSE, digits = 4)
    if (nrow(x$stratum_estimates) > 10) {
      cat(sprintf("  ... and %d more strata\n", nrow(x$stratum_estimates) - 10))
    }
    if (is_lng) {
      cat("  (random slope not shown; use predict(type = \"strata\") for the ",
          "per-stratum intercept and slope, or plot(type = \"trajectories\")).\n",
          sep = "")
    }
  }

  invisible(x)
}
