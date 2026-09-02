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

#' Fixed-effect table with Wald statistics
#'
#' Internal helper that assembles the \code{fixed_effects} slot of a
#' \code{maihda_summary} for the frequentist engines: the Wald statistic
#' (\code{estimate / se}), its two-sided p-value and the matching Wald interval
#' at \code{conf_level}. With \code{df} the reference is a \eqn{t} on those
#' degrees of freedom, otherwise the normal.
#'
#' @param term Character vector of coefficient names.
#' @param estimate Numeric vector of point estimates.
#' @param se Numeric vector of standard errors (\code{NULL} for none).
#' @param conf_level Interval level.
#' @param df Denominator degrees of freedom, named by \code{term} or in
#'   \code{term} order; \code{NULL} (default) for the normal approximation.
#' @return A data frame with \code{term}, \code{estimate}, \code{se},
#'   \code{statistic}, \code{df}, \code{p_value}, \code{lower} and \code{upper}.
#'   \code{df} is \code{NA} on the normal path.
#' @keywords internal
maihda_fixed_effects_table <- function(term, estimate, se, conf_level = 0.95,
                                       df = NULL) {
  term <- as.character(term)
  estimate <- as.numeric(estimate)
  se <- if (is.null(se)) rep(NA_real_, length(estimate)) else as.numeric(se)
  # A zero/NA SE (boundary or non-positive-definite Hessian) leaves the Wald
  # quantities undefined rather than infinite.
  se[!is.na(se) & se <= 0] <- NA_real_
  statistic <- estimate / se

  # Degrees of freedom are matched by name where they carry one, so a df vector
  # ordered differently from `term` (or missing a term) cannot silently shift.
  if (is.null(df)) {
    df <- rep(NA_real_, length(estimate))
  } else {
    df <- if (!is.null(names(df))) as.numeric(df)[match(term, names(df))]
          else as.numeric(df)[seq_along(term)]
    df[!is.na(df) & df <= 0] <- NA_real_
  }

  # A missing df falls back to the normal, which is the t's own limit, so the
  # two paths stay consistent within one table. as.numeric() keeps the empty
  # table numeric: ifelse() on a zero-length test returns logical(0).
  safe_df <- ifelse(is.na(df), 1, df)
  crit <- as.numeric(ifelse(is.na(df), stats::qnorm((1 + conf_level) / 2),
                            stats::qt((1 + conf_level) / 2, df = safe_df)))
  p_value <- as.numeric(ifelse(is.na(df), 2 * stats::pnorm(-abs(statistic)),
                               2 * stats::pt(-abs(statistic), df = safe_df)))

  data.frame(
    term      = term,
    estimate  = estimate,
    se        = se,
    statistic = statistic,
    df        = as.numeric(df),
    p_value   = p_value,
    lower     = estimate - crit * se,
    upper     = estimate + crit * se,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

#' Tag a summary with its role in a two-model analysis
#'
#' Internal helper. The null and adjusted summaries of a \code{\link{maihda}}
#' analysis are indistinguishable once pulled out of the analysis object, which
#' makes a printed summary easy to misread -- the null model's fixed effects are
#' the intercept and covariates only, because the strata dimensions are its
#' random-effect grouping rather than fixed effects. Stamping the role lets
#' \code{\link{print.maihda_summary}} say which model it is showing.
#'
#' @param s A \code{maihda_summary}, or \code{NULL}.
#' @param role \code{"null"} or \code{"adjusted"}.
#' @return \code{s} with a \code{"maihda_role"} attribute (\code{NULL} in, \code{NULL} out).
#' @keywords internal
maihda_tag_role <- function(s, role) {
  if (is.null(s)) return(NULL)
  attr(s, "maihda_role") <- role
  s
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
#'   latent/distributional approximation (\eqn{\pi^2/3} for logistic;
#'   \eqn{\log(1 + 1/\lambda)} for Poisson per Stryhn et al. 2006 and
#'   \eqn{\log(1 + 1/\lambda + 1/\theta)} for the negative binomial per Nakagawa,
#'   Johnson & Schielzeth 2017 -- their \code{"delta"} and \code{"trigamma"}
#'   alternatives are available via \code{fit_maihda(count_approximation = )}, and
#'   the choice is reported in the printed summary and in \code{count_vpc} because
#'   the three diverge materially below a marginal count of 2 -- each evaluated at a single \emph{marginal}
#'   expected count \eqn{\lambda}: the mean over the analytic sample of the
#'   row-level \eqn{\lambda_i = \exp(x_i'\beta + v_i/2)} -- the fixed-part
#'   prediction with the log-normal correction for the row's total random-effect
#'   variance \eqn{v_i}. The counts are averaged \emph{before} the transform,
#'   which is where the cited \eqn{\lambda} is defined and which reduces to
#'   Nakagawa et al.'s \eqn{\lambda = \exp(\beta_0 + \sigma^2/2)} in the null
#'   model; \emph{not} at the conditional fitted means, whose BLUPs would tie the
#'   level-1 variance to the realized random effects), so the
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
#'   interval estimates). For a Gaussian model carrying lme4 precision
#'   \code{weights}, the simulated responses draw each residual at
#'   \eqn{\sigma / \sqrt{w_i}}, so the interval rests on the same
#'   \eqn{\sigma^2 / w_i} semantics as the point estimate above.
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
#' @param df_method Reference distribution for the fixed-effect p-values and
#'   Wald intervals of a Gaussian \code{lme4} fit: \code{"between-within"}
#'   (default) a \eqn{t} on containment degrees of freedom, \code{"normal"} a z.
#'   Every other engine uses a z regardless.
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
#'   \item{longitudinal}{For a longitudinal (growth-curve) fit, the time-varying
#'     summary: \code{vpc_t} (the VPC over a reporting grid, with bootstrap or
#'     posterior bands), the per-level variances over that grid, the stratum and
#'     individual covariance blocks, and the two \emph{trajectory VPCs}
#'     \code{vpc_intercept} and \code{vpc_slope} described below. \code{NULL} for a
#'     cross-sectional fit}
#'   \item{context}{For a contextual cross-classified fit, the stratum vs.
#'     context partition: per-context variances and shares, the contexts' total
#'     share (\code{vpc_context_total}, with an interval when bootstrapped or for
#'     brms), and the between-stratum share (\code{vpc_stratum}); \code{NULL}
#'     otherwise}
#'   \item{discriminatory_accuracy}{For a binomial/Bernoulli outcome, the
#'     \code{maihda_da} object (AUC + MOR) from
#'     \code{\link{maihda_discriminatory_accuracy}}; \code{NULL} otherwise. A
#'     contextual fit (\code{fit_maihda(context = )}) is included -- its headline AUC
#'     is the intersectional-scope concordance that excludes the context random
#'     effect. \code{NULL} for a crossed-dimensions fit (whose headline here is the
#'     additive/interaction decomposition) and a longitudinal fit}
#'   \item{count_vpc}{For a log-link count model, the \code{approximation} the
#'     level-1 variance came from, the marginal count \code{lambda} (and
#'     \code{theta} / \code{lambda_effective} for the negative binomial) it was
#'     evaluated at, the \code{alternatives} all three approximations give at that
#'     \code{lambda} (\code{level1_variance} is the one used), and
#'     \code{low_count} -- \code{TRUE} below the \eqn{\lambda = 2} threshold above
#'     which Nakagawa et al. (2017) report the three agree. These are plug-in
#'     values at a single \code{lambda}: on the likelihood engines that is exactly
#'     the number in \code{variance_components}, but a \code{brms} summary works
#'     draw by draw and reports \eqn{E[\sigma^2_e]} there, which differs slightly.
#'     \code{NULL} for every other family.}
#'   \item{vpc_response}{The response-scale VPC (\code{maihda_vpc_response}) when
#'     \code{response_vpc = TRUE} for a binomial lme4 model, including a contextual
#'     fit (the context variance enters the VPC denominator); \code{NULL} otherwise
#'     (including for crossed-dimensions and longitudinal fits)}
#'   \item{stratum_estimates}{Data frame of stratum-specific random effects with labels if available}
#'   \item{fixed_effects}{Fixed-effect estimates. For the lme4, WeMix and ordinal
#'     engines a data frame with \code{term}, \code{estimate}, \code{se},
#'     \code{statistic}, \code{df}, \code{p_value} and the Wald interval
#'     \code{lower}/\code{upper} at \code{conf_level}; \code{df} is \code{NA}
#'     wherever the reference is a z. The WeMix standard errors are its sandwich
#'     (robust) ones. For brms, the \code{brms::fixef()} matrix (posterior mean,
#'     \code{Est.Error} and the credible-interval quantiles at
#'     \code{conf_level}). Available in a tidy, engine-independent shape from
#'     \code{tidy(x, component = "fixed")}}
#'   \item{conf_level}{The interval level used for the fixed effects (and, when
#'     bootstrapped or Bayesian, the VPC)}
#'   \item{df_method}{The reference distribution the \code{fixed_effects} table
#'     used, \code{"between-within"} or \code{"normal"}}
#'   \item{thresholds}{For a cumulative (ordinal) clmm fit, the threshold (cut
#'     point) estimates with standard errors -- the cumulative model's
#'     "intercepts"; NULL otherwise}
#'   \item{model_summary}{Original model summary}
#'   \item{diagnostics}{Fit-quality diagnostics (singular fit / convergence)
#'     carried over from the fitted model and reported by the print method}
#'   \item{strata_autobin_info}{The auto-binning recipe carried over from the
#'     fitted model: for each numeric stratum dimension \code{make_strata()}
#'     discretised, its \code{breaks} and \code{labels}. The cut-points are
#'     quantiles of the analytic sample, so they define the strata (and hence the
#'     estimand); the print method reports them. An empty list when nothing was
#'     binned}
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
#' latent-scale approximation. For the count families (Poisson, negative
#' binomial) the marginal expected counts are propagated \emph{per draw} -- from
#' the fixed-part linear-predictor draws and each draw's total random-intercept
#' variance -- for the intercept-only VPC structures (strata, crossed-dimensions,
#' contextual), so the credible interval reflects fixed-effect and random-effect
#' variance uncertainty; the negative-binomial \code{shape} draws are always
#' propagated. A random-slope (longitudinal) structure, whose per-row
#' random-effect design that fast path does not carry, instead holds the marginal
#' expected counts at a posterior-mean plug-in (constant across draws) to avoid an
#' expensive \eqn{ndraws \times nobs} computation.
#'
#' @section Fixed-effect degrees of freedom:
#' A Gaussian \code{lme4} fit refers each Wald statistic to a \eqn{t} on
#' \emph{containment} (between-within) degrees of freedom, reported in the
#' \code{df} column: a term absorbed by a random-effect grouping is tested
#' against that grouping's units minus the terms it absorbs, and a term absorbed
#' by none against \eqn{n} minus the random-effect levels. A random slope counts
#' as absorbing. Set \code{df_method = "normal"} for a z instead.
#'
#' A GLMM, a WeMix pseudo-ML fit and an \code{ordinal::clmm} fit have no
#' finite-sample \eqn{t} and use the Wald z; a brms summary reports the
#' posterior. For Kenward-Roger or Satterthwaite, apply \pkg{pbkrtest} or
#' \pkg{lmerTest} to \code{x$model}.
#'
#' @section Two VPCs for a longitudinal fit:
#' A longitudinal summary reports \strong{two different variance partitions}, and
#' they are not interchangeable. They differ in one term of the denominator: the
#' level-1 (occasion) residual variance \eqn{\sigma^2_e}, which in a growth model is
#' \emph{within-individual volatility} -- how far a single measurement falls from
#' that person's own smooth trajectory. It mixes measurement error, genuine
#' short-term fluctuation, and any misfit of the assumed functional form. A
#' cross-sectional MAIHDA cannot separate it from between-individual variance at
#' all; repeated measures are what split the two.
#'
#' \strong{The headline VPC} (\code{vpc}, and \code{longitudinal$vpc_t} over time)
#' keeps it:
#' \deqn{VPC_S(t) = \frac{Var_S(t)}{Var_S(t) + Var_I(t) + \sigma^2_e}.}
#' This is the discriminatory-accuracy question -- how much of an \emph{observed
#' measurement} at time \eqn{t} a stratum accounts for -- and it is the quantity
#' comparable to published cross-sectional MAIHDA VPCs. Report this one unless you
#' specifically mean the trajectory question.
#'
#' \strong{The trajectory VPCs} (\code{longitudinal$vpc_intercept} and
#' \code{longitudinal$vpc_slope}) drop it, following Bell et al. (2024), equation
#' (5):
#' \deqn{VPC_{intercept} = \frac{Var_S(t_0)}{Var_S(t_0) + Var_I(t_0)}, \qquad
#'       VPC_{slope} = \frac{SlopeVar_S(t_0)}{SlopeVar_S(t_0) + SlopeVar_I(t_0)}.}
#' These ask how intersectionally patterned people's \emph{trajectories} are i.e., what
#' share of the between-individual variation in where a trajectory starts, and in how
#' fast it changes, lies between strata. Because \eqn{\sigma^2_e} is absent neither is
#' affected by how noisy the outcome measure is, which makes them comparable across
#' studies using different instruments.
#'
#' Both are evaluated at the baseline \eqn{t_0} (\code{ref_time}, the earliest
#' observed time), pairing with \code{PCV_intercept} and \code{PCV_slope} from
#' \code{maihda(decomposition = "longitudinal")}. The intercept VPC depends on where
#' time is zeroed and the slope VPC does not, as Bell et al. note; their own examples
#' centre on mean age rather than the baseline, so an intercept VPC replicated from
#' the paper will differ from the one reported here unless the reference points are
#' aligned. \code{vpc_slope} is \code{NA} when the model was fit with
#' \code{stratum_slope = FALSE} (no between-stratum slope variance exists to take a
#' share of).
#'
#' @references
#' Bell, A., Evans, C., Holman, D., & Leckie, G. (2024). Extending intersectional
#' multilevel analysis of individual heterogeneity and discriminatory accuracy
#' (MAIHDA) to study individual longitudinal trajectories, with application to
#' mental health in the UK. \emph{Social Science & Medicine}, 351, 116955.
#' \doi{10.1016/j.socscimed.2024.116955}
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
                          conf_level = 0.95, response_vpc = FALSE, seed = NULL,
                          df_method = c("between-within", "normal"), ...) {
  if (!inherits(object, "maihda_model")) {
    stop("'object' must be a maihda_model object from fit_maihda()")
  }
  df_method <- match.arg(df_method)

  if (!is.logical(bootstrap) || length(bootstrap) != 1 || is.na(bootstrap)) {
    stop("'bootstrap' must be TRUE or FALSE.", call. = FALSE)
  }
  # Validated for every engine and both bootstrap settings: conf_level also sets
  # the fixed-effect interval below, not just the bootstrap VPC interval.
  conf_level <- maihda_validate_conf_level(conf_level)
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
      var_residual <- maihda_residual_variance_lme4(
        model, vc, approximation = maihda_count_approximation(object))

      # Calculate VPC (ICC)
      vpc <- var_random / (var_random + var_other_random + var_residual)

      # Create variance components data frame
      variance_components <- maihda_variance_components_table(
        var_random, var_other_random, var_residual
      )

      # Bootstrap confidence intervals for VPC if requested
      if (bootstrap) {
        vpc_ci <- bootstrap_vpc(model, object$data, object$formula, n_boot,
                                conf_level,
                                approximation = maihda_count_approximation(object))
        vpc_result <- list(
          estimate = vpc,
          ci_lower = vpc_ci[1],
          ci_upper = vpc_ci[2],
          conf_level = conf_level,
          bootstrap = TRUE,
          method = "bootstrap",
          n_boot_ok = attr(vpc_ci, "n_ok"),
          n_boot_nonconverged = attr(vpc_ci, "n_nonconverged"),
          mc_se = attr(vpc_ci, "mc_se")
        )
      } else {
        vpc_result <- list(
          estimate = vpc,
          bootstrap = FALSE
        )
      }
    }

    # Get model summary
    model_summary <- summary(model)

    # Fixed effects with their standard errors, Wald statistics and intervals,
    # read off the model summary's coefficient matrix (Estimate / Std. Error --
    # identical to lme4::fixef() plus sqrt(diag(vcov()))). lme4 reports no
    # p-value for a Gaussian fit by design; the one added here is the Wald p that
    # matches the interval (see maihda_fixed_effects_table). A Gaussian fit gets
    # containment (between-within) denominator degrees of freedom, so a dimension
    # main effect is tested against the number of strata rather than n; a GLMM has
    # no finite-sample t reference and stays on the normal.
    fe_coef <- stats::coef(model_summary)
    fe_df <- if (df_method == "between-within") maihda_containment_df(model) else NULL
    fixed_effects <- maihda_fixed_effects_table(
      term = rownames(fe_coef),
      estimate = fe_coef[, "Estimate"],
      se = fe_coef[, "Std. Error"],
      conf_level = conf_level,
      df = fe_df
    )

    # Stratum (intersection) random-effect estimates -- the interaction residuals in
    # the cross-classified model (the named interaction group), or the single stratum
    # RE otherwise.
    interaction_group <- if (!is.null(cc)) cc$interaction_group else "stratum"
    # A longitudinal growth block carries time slopes by design; its intercept
    # column (the deviation at the model's coefficient origin) is reported here
    # while the full trajectory lives in the `longitudinal` slot, so it opts out
    # of the extractor's intercept-only guard. Every other path reaching this
    # line is already slope-free (the ordinary branch validated the whole model
    # above; the crossed builder rejects slope terms).
    stratum_estimates <- maihda_stratum_ranef_lme4(model, group = interaction_group,
                                                   allow_slope_columns = !is.null(lng))
    stratum_estimates <- add_stratum_labels(stratum_estimates, object$strata_info)

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

    # WeMix reports sandwich (robust) fixed-effect standard errors -- include them
    # alongside the estimates (the lme4 table reports estimates only). They cover
    # the weighting and within-stratum dependence, NOT the clustering or
    # stratification of a complex sample design, which this wrapper has no way to
    # represent from a single person-weight column (see fit_maihda's
    # `sampling_weights`); do not describe them as design-consistent.
    fixed_effects <- maihda_fixed_effects_table(
      term = names(object$model$coef),
      estimate = as.numeric(object$model$coef),
      se = as.numeric(object$model$SE[names(object$model$coef)]),
      conf_level = conf_level
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
      fixed_effects <- maihda_fixed_effects_table(character(0), numeric(0),
                                                  numeric(0), conf_level)
    } else {
      V <- tryCatch(stats::vcov(object$model), error = function(e) NULL)
      beta_se <- rep(NA_real_, length(beta))
      if (!is.null(V) && all(names(beta) %in% rownames(V))) {
        beta_se <- sqrt(pmax(diag(V)[names(beta)], 0))
      }
      fixed_effects <- maihda_fixed_effects_table(
        term = names(beta),
        estimate = as.numeric(beta),
        se = as.numeric(beta_se),
        conf_level = conf_level
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
      vpc_draws <- maihda_vpc_draws_brms(
        model, conf_level = conf_level,
        approximation = maihda_count_approximation(object))

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

    # Extract fixed effects. Kept in brms::fixef()'s matrix form (Estimate /
    # Est.Error / two quantile columns); the quantiles are taken at conf_level,
    # so the interval matches the rest of this summary. A posterior summary has
    # no test statistic or p-value -- tidy() reports those as NA for brms.
    fixed_effects <- brms::fixef(
      model,
      probs = c((1 - conf_level) / 2, 1 - (1 - conf_level) / 2)
    )

    interaction_group <- if (!is.null(cc)) cc$interaction_group else "stratum"
    # Longitudinal opt-out as in the lme4 branch above: the growth block's time
    # slopes are reported via the `longitudinal` slot, not this table.
    stratum_estimates <- maihda_stratum_ranef_brms(model, group = interaction_group,
                                                   allow_slope_columns = !is.null(lng))
    stratum_estimates <- add_stratum_labels(stratum_estimates, object$strata_info)

    model_summary <- summary(model)
  }

  # Family-defined companions to the VPC for a binomial/Bernoulli outcome: the
  # discriminatory accuracy (AUC + MOR) -- the "DA" in MAIHDA -- always, and the
  # response-scale VPC on request (it is a simulation, hence opt-in and seeded).
  # Both summarise the single fitted model with no refit, mirroring how the
  # crossed-dimensions additive/interaction `decomposition` slot above is computed
  # in this same summary layer so that fit_maihda() and maihda() share the logic.
  #
  # A CONTEXTUAL fit (fit_maihda(context = )) is INCLUDED: the DA's headline AUC is
  # the intersectional-scope concordance that EXCLUDES the context random effect
  # (auc_full carries the all-effects value separately), and maihda_vpc_response()
  # integrates the context variance into the VPC denominator -- so both sit on the
  # same stratum-vs-context estimand this summary reports, not the mismatched one
  # the earlier skip guarded against. Skipped only for a crossed-dimensions fit --
  # whose headline here is the additive/interaction decomposition above, the DA /
  # response-VPC companions remaining available from the standalone helpers -- and a
  # longitudinal fit, whose VPC is time-varying. Wrapped so a bonus summary never
  # breaks the core VPC.
  discriminatory_accuracy <- NULL
  vpc_response <- NULL
  fam_name <- tryCatch(maihda_model_family_name(object),
                       error = function(e) NA_character_)
  if (is.null(cc) && is.null(lng) &&
      isTRUE(fam_name %in% c("binomial", "bernoulli"))) {
    discriminatory_accuracy <- maihda_try_optional(
      maihda_discriminatory_accuracy(object),
      "Discriminatory accuracy (AUC/MOR)")
    if (isTRUE(response_vpc) && identical(engine, "lme4") &&
        identical(fam_name, "binomial")) {
      vpc_response <- maihda_try_optional(
        maihda_vpc_response(object, seed = seed), "Response-scale VPC")
    }
  }

  # Which count level-1 approximation the VPC used, and the marginal lambda it was
  # evaluated at -- Nakagawa et al. (2017) ask that the choice be reported because
  # the three approximations diverge materially at low counts, and this is the only
  # place the user sees it. NULL for every non-count family. Wrapped like the other
  # optional summaries so a failure here never breaks the VPC.
  count_vpc <- maihda_try_optional(maihda_count_vpc_note(object),
                                   "Count VPC approximation note")

  # Create summary object
  result <- structure(
    list(
      vpc = vpc_result,
      count_vpc = count_vpc,
      variance_components = variance_components,
      decomposition = decomposition,
      context = context_summary,
      longitudinal = longitudinal,
      discriminatory_accuracy = discriminatory_accuracy,
      vpc_response = vpc_response,
      stratum_estimates = stratum_estimates,
      fixed_effects = fixed_effects,
      conf_level = conf_level,
      # Which reference distribution the fixed-effect table used. Only the
      # Gaussian lme4 path can honour "between-within"; every other engine
      # reports "normal" whatever was asked for, so the print method and any
      # downstream reader can say what was actually done.
      df_method = if (is.data.frame(fixed_effects) &&
                      "df" %in% names(fixed_effects) &&
                      any(!is.na(fixed_effects$df))) "between-within" else "normal",
      thresholds = thresholds,
      model_summary = model_summary,
      engine = engine,
      cc_info = cc,
      context_info = ctx,
      longitudinal_info = lng,
      diagnostics = object$diagnostics,
      strata_autobin_info = object$strata_autobin_info
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
  var_within <- maihda_residual_variance_lme4(
    model, vc, approximation = maihda_count_approximation(object))
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
    boot <- bootstrap_cc(model, cc, n_boot, conf_level, ctx_vars = ctx_vars,
                         approximation = maihda_count_approximation(object))
    vpc_result <- list(
      estimate = part$vpc,
      ci_lower = boot$vpc[1],
      ci_upper = boot$vpc[2],
      conf_level = conf_level,
      bootstrap = TRUE,
      method = "bootstrap",
      n_boot_ok = attr(boot$vpc, "n_ok"),
      n_boot_nonconverged = attr(boot$vpc, "n_nonconverged"),
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
  within_draws <- maihda_residual_variance_draws_brms(
    model, draws, approximation = maihda_count_approximation(object))

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
#' @param approximation The count level-1 variance approximation of the
#'   fitted model (\code{maihda_count_approximation()}); inert for every
#'   non-count family.
#' @return A list with \code{vpc}, \code{additive_share}, \code{interaction_share}
#'   (and \code{context_vpc} when \code{ctx_vars} is non-empty), each a length-2
#'   interval carrying \code{n_ok}/\code{mc_se} attributes.
#' @keywords internal
#' @importFrom lme4 refit
bootstrap_cc <- function(model, cc, n_boot, conf_level, ctx_vars = character(0),
                         approximation = "lognormal") {
  vpc_boot <- rep(NA_real_, n_boot)
  additive_boot <- rep(NA_real_, n_boot)
  interaction_boot <- rep(NA_real_, n_boot)
  context_boot <- rep(NA_real_, n_boot)
  # Count contributing draws whose refit optimiser did not converge, so the reported
  # n_boot_ok does not silently imply convergence (see bootstrap_vpc()).
  n_nonconv <- 0L
  sim_data <- maihda_simulate_lme4(model, nsim = n_boot)

  for (i in seq_len(n_boot)) {
    tryCatch({
      boot_model <- lme4::refit(model, newresp = sim_data[[i]])
      var_named <- maihda_random_variances_lme4(boot_model)
      var_within <- maihda_residual_variance_lme4(boot_model,
                                                  approximation = approximation)
      split <- maihda_cc_variance_split(var_named, cc$dim_groups, cc$interaction_group)
      var_context <- if (length(ctx_vars) > 0) sum(var_named[ctx_vars]) else 0
      part <- maihda_cc_partition(split$additive, split$interaction, var_within,
                                  var_context)
      vpc_boot[i] <- part$vpc
      additive_boot[i] <- part$additive_share
      interaction_boot[i] <- part$interaction_share
      context_boot[i] <- var_context / part$total
      if (maihda_lme4_optimizer_failed(boot_model)) n_nonconv <- n_nonconv + 1L
    }, error = function(e) NULL)
  }

  # Non-converged draws are retained but reported (see bootstrap_pcv()).
  if (n_nonconv > 0) {
    warning(sprintf(paste0(
      "%d of %d contributing crossed-dimensions decomposition bootstrap draw(s) had ",
      "an lme4 optimiser that did not converge; they are retained in the intervals. ",
      "n_boot_ok counts converged and non-converged refits alike -- interpret the ",
      "intervals accordingly."), n_nonconv, sum(is.finite(vpc_boot))), call. = FALSE)
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
  attr(out$vpc, "n_nonconverged") <- n_nonconv
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
  var_within <- maihda_residual_variance_lme4(
    model, vc, approximation = maihda_count_approximation(object))
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
    boot <- bootstrap_context(model, ctx_vars, n_boot, conf_level,
                              approximation = maihda_count_approximation(object))
    vpc_result <- list(
      estimate = part$vpc_stratum,
      ci_lower = boot$vpc[1],
      ci_upper = boot$vpc[2],
      conf_level = conf_level,
      bootstrap = TRUE,
      method = "bootstrap",
      n_boot_ok = attr(boot$vpc, "n_ok"),
      n_boot_nonconverged = attr(boot$vpc, "n_nonconverged"),
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
  within_draws <- maihda_residual_variance_draws_brms(
    model, draws, approximation = maihda_count_approximation(object))

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
#' @param approximation The count level-1 variance approximation of the
#'   fitted model (\code{maihda_count_approximation()}); inert for every
#'   non-count family.
#' @return A list with \code{vpc} (between-stratum share) and \code{context_vpc}
#'   (contexts' total share), each a length-2 interval carrying
#'   \code{n_ok}/\code{mc_se} attributes.
#' @keywords internal
#' @importFrom lme4 refit
bootstrap_context <- function(model, ctx_vars, n_boot, conf_level,
                              approximation = "lognormal") {
  vpc_boot <- rep(NA_real_, n_boot)
  context_boot <- rep(NA_real_, n_boot)
  # Count contributing draws whose refit optimiser did not converge, so the reported
  # n_boot_ok does not silently imply convergence (see bootstrap_vpc()).
  n_nonconv <- 0L
  sim_data <- maihda_simulate_lme4(model, nsim = n_boot)

  for (i in seq_len(n_boot)) {
    tryCatch({
      boot_model <- lme4::refit(model, newresp = sim_data[[i]])
      var_named <- maihda_random_variances_lme4(boot_model)
      var_within <- maihda_residual_variance_lme4(boot_model,
                                                  approximation = approximation)
      var_stratum <- unname(var_named[["stratum"]])
      per_context <- var_named[ctx_vars]
      var_other <- max(0, sum(var_named, na.rm = TRUE) - var_stratum -
                         sum(per_context))
      part <- maihda_context_partition(var_stratum, as.list(per_context),
                                       var_within, var_other)
      vpc_boot[i] <- part$vpc_stratum
      context_boot[i] <- part$vpc_context_total
      if (maihda_lme4_optimizer_failed(boot_model)) n_nonconv <- n_nonconv + 1L
    }, error = function(e) NULL)
  }

  # Non-converged draws are retained but reported (see bootstrap_pcv()).
  if (n_nonconv > 0) {
    warning(sprintf(paste0(
      "%d of %d contributing contextual decomposition bootstrap draw(s) had an lme4 ",
      "optimiser that did not converge; they are retained in the intervals. ",
      "n_boot_ok counts converged and non-converged refits alike -- interpret the ",
      "intervals accordingly."), n_nonconv, sum(is.finite(vpc_boot))), call. = FALSE)
  }

  out <- list(
    vpc = maihda_bootstrap_ci(vpc_boot, n_boot, conf_level, "VPC"),
    context_vpc = maihda_bootstrap_ci(context_boot, n_boot, conf_level,
                                      "context VPC")
  )
  attr(out$vpc, "n_nonconverged") <- n_nonconv
  out
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
#' @param approximation The count level-1 variance approximation of the fitted
#'   model (\code{maihda_count_approximation()}); inert for every non-count family.
#'
#' @return A vector with lower and upper confidence bounds
#' @keywords internal
#' @importFrom lme4 lmer glmer VarCorr
bootstrap_vpc <- function(model, data, formula, n_boot, conf_level,
                          approximation = "lognormal") {
  # Initialise to NA so iterations whose refit() throws -- and never reach the
  # assignment inside the tryCatch body -- stay NA rather than the numeric() default of 0.
  # The error handler runs in its own scope and cannot write back to this vector,
  # so the initial value is what survives a failure.
  vpc_boot <- rep(NA_real_, n_boot)
  # Count draws that contribute to the interval but whose refit optimiser did not
  # converge, so the reported n_boot_ok does not silently imply convergence.
  n_nonconv <- 0L
  sim_data <- maihda_simulate_lme4(model, nsim = n_boot)

  for (i in 1:n_boot) {
    tryCatch({
      boot_model <- lme4::refit(model, newresp = sim_data[[i]])

      # Calculate VPC
      vc <- lme4::VarCorr(boot_model)
      var_random <- maihda_stratum_variance_lme4(boot_model)
      var_total_random <- maihda_total_random_variance_lme4(boot_model)
      var_other_random <- max(0, var_total_random - var_random)
      var_residual <- maihda_residual_variance_lme4(
        boot_model, vc, approximation = approximation)

      vpc_boot[i] <- var_random / (var_random + var_other_random + var_residual)
      if (maihda_lme4_optimizer_failed(boot_model)) {
        n_nonconv <- n_nonconv + 1L
      }
    }, error = function(e) NULL)
  }

  # Non-converged draws are retained but reported (see bootstrap_pcv()).
  if (n_nonconv > 0) {
    warning(sprintf(paste0(
      "%d of %d contributing VPC bootstrap draw(s) had an lme4 optimiser that did ",
      "not converge; they are retained in the interval. n_boot_ok counts converged ",
      "and non-converged refits alike -- interpret the interval accordingly."),
      n_nonconv, sum(is.finite(vpc_boot))), call. = FALSE)
  }

  # Reduce to an interval, requiring a minimum number of successful refits.
  ci <- maihda_bootstrap_ci(vpc_boot, n_boot, conf_level, "VPC")
  attr(ci, "n_nonconverged") <- n_nonconv

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
    if (!is.finite(est)) {
      # NA when there is no between-strata variance to split (see maihda_cc_partition).
      return("NA (no between-strata variance to split)")
    }
    if (!is.null(ci) && length(ci) == 2L && all(is.finite(ci))) {
      sprintf("%.1f%% [%.1f%%, %.1f%%]", est * 100, ci[1] * 100, ci[2] * 100)
    } else {
      sprintf("%.1f%%", est * 100)
    }
  }
  cat(maihda_palette()$bold("Additive vs. Intersectional Decomposition (crossed-dimensions):"), "\n", sep = "")
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
  cat(maihda_palette()$muted(paste0(
      "  Note: the additive share is the crossed-dimensions analogue of the PCV but\n",
      "  a different estimator; interpret the interaction share cautiously.\n\n")))
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
  cat(maihda_palette()$bold("Contextual Cross-Classified Partition (stratum x context):"), "\n", sep = "")
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
  cat(maihda_palette()$muted(paste0(
      "  Note: the headline VPC/ICC is the between-stratum share conditional on\n",
      "  the context random effect(s). The context share is the between-context\n",
      "  component of the unexplained variance.\n\n")))
  invisible(NULL)
}

#' Print method for maihda_summary objects
#'
#' @param x A maihda_summary object
#' @param ... Additional arguments (not used)
#' @return No return value, called for side effects.
#' @export
print.maihda_summary <- function(x, ...) {
  pal <- maihda_palette()
  cat(pal$bold("MAIHDA Model Summary"), "\n", sep = "")
  cat("====================\n\n")

  # Say which of a two-model analysis's models this is: the null model's fixed
  # effects are the intercept and covariates only (its strata dimensions are the
  # random-effect grouping), so an unlabelled print reads as if the dimensions
  # were missing. Tagged by maihda(); absent for a bare fit_maihda() summary.
  role <- attr(x, "maihda_role")
  if (identical(role, "null")) {
    cat(pal$muted("Null model. Use which = \"adjusted\" for the adjusted model.\n\n"))
  } else if (identical(role, "adjusted")) {
    cat(pal$muted("Adjusted model.\n\n"))
  }
  maihda_print_fit_diagnostics(x$diagnostics)
  # The auto-bin recipe belongs next to the numbers it produced: a tertile-binned
  # numeric dimension makes the strata (and so the VPC/PCV below) data-dependent,
  # and make_strata()'s fit-time message is long gone by the time a saved model is
  # printed in another session.
  maihda_print_autobin_info(x$strata_autobin_info)

  is_lng <- !is.null(x$longitudinal)
  if (is_lng) {
    cat(sprintf("Variance Partition Coefficient (VPC/ICC) at baseline (%s = %g):\n",
                x$longitudinal$time, x$longitudinal$ref_time))
  } else {
    cat("Variance Partition Coefficient (VPC/ICC):\n")
  }
  if (maihda_vpc_has_interval(x$vpc)) {
    cat(sprintf("  Estimate: %s [%.4f, %.4f]\n",
                pal$accent(sprintf("%.4f", x$vpc$estimate)), x$vpc$ci_lower, x$vpc$ci_upper))
    cat("  ", pal$muted(maihda_vpc_interval_label(x$vpc)), "\n", sep = "")
    if (!is.null(x$vpc$mc_se) && is.finite(x$vpc$mc_se)) {
      cat(pal$muted(sprintf(
        "  (%d successful bootstrap draws; Monte Carlo SE of the bootstrap mean %.4f)\n",
        as.integer(x$vpc$n_boot_ok), x$vpc$mc_se)))
    }
    cat("\n")
  } else {
    cat(sprintf("  Estimate: %s\n\n", pal$accent(sprintf("%.4f", x$vpc$estimate))))
  }
  maihda_print_count_vpc_note(x$count_vpc, pal)

  cat(pal$bold("Variance Components:"), "\n", sep = "")
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
    # Bell et al. (2024) eq. (5). A DIFFERENT denominator from the VPC above -- the
    # occasion-level residual is excluded -- so they answer a different question.
    # Only the intercept VPC is ordered against the VPC above (larger, at the
    # reference time); the slope VPC is on the slope scale and is not comparable to
    # it at all. Kept to a few lines here; ?summary.maihda_model and the
    # longitudinal vignette carry the full contrast.
    tv_i <- x$longitudinal$vpc_intercept
    tv_s <- x$longitudinal$vpc_slope
    if (!is.null(tv_i) || !is.null(tv_s)) {
      fmt_tv <- function(v) if (isTRUE(is.finite(v))) sprintf("%.4f", v) else "NA"
      cat("Trajectory VPCs (Bell et al. 2024, eq. 5; occasion-level variance excluded):\n")
      cat(sprintf("  Intercept (%s = %g): %s    Slope: %s\n",
                  x$longitudinal$time, x$longitudinal$ref_time,
                  pal$accent(fmt_tv(tv_i)), pal$accent(fmt_tv(tv_s))))
      # One sentence on WHICH QUESTION each answers -- the part a reader cannot infer
      # from the header, and the one that stops the larger number being reported as
      # "the VPC". The ordering claim is deliberately split: "at the reference time"
      # is load-bearing for the intercept (it exceeds the VPC at the same time point,
      # but not necessarily its maximum over the whole trajectory), and the slope VPC
      # gets no ordering at all because its numerator is a slope variance.
      cat("  These ask how intersectionally patterned trajectories are; the VPC\n")
      cat("  above asks how much of an observed measurement a stratum explains.\n")
      cat("  Intercept exceeds it at the reference time; Slope is on a different\n")
      cat("  scale and is not comparable to it. See ?summary.maihda_model.\n\n")
    }
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

  # The frequentist engines carry Wald statistics and an interval (see
  # maihda_fixed_effects_table); brms carries the posterior summary matrix,
  # whose column names already label the quantiles.
  fe_level <- if (is.null(x$conf_level)) 95 else 100 * x$conf_level
  fe_print <- x$fixed_effects
  # A t reference is only ever used where finite-sample degrees of freedom were
  # available (a Gaussian lme4 fit); everywhere else the df column is all-NA and
  # is dropped rather than printed as a column of blanks.
  fe_has_df <- is.data.frame(fe_print) && "df" %in% names(fe_print) &&
    any(!is.na(fe_print$df))
  fe_header <- if (is.data.frame(fe_print) && "lower" %in% names(fe_print)) {
    sprintf("Fixed Effects (Wald %s, %s%% intervals):",
            if (fe_has_df) "t" else "z", format(round(fe_level, 1), trim = TRUE))
  } else {
    "Fixed Effects:"
  }
  cat(pal$bold(fe_header), "\n", sep = "")
  if (is.data.frame(fe_print) && "p_value" %in% names(fe_print)) {
    # Without this a p of 1e-200 prints as "0.00000" next to a p of 0.04.
    fe_print$p_value <- format.pval(fe_print$p_value, digits = 3, eps = 1e-16)
  }
  if (is.data.frame(fe_print) && "df" %in% names(fe_print) && !fe_has_df) {
    fe_print$df <- NULL
  }
  print(fe_print, row.names = FALSE, digits = 4)
  if (fe_has_df) {
    # One line, not a lecture: without it the small df on a dimension main effect
    # reads as a typo, but this prints on every Gaussian summary.
    cat("  df: containment (between-within).\n")
  }
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
