# Response-scale (observed-scale) VPC for binomial and count MAIHDA models.
#
# The package's default VPC for a logistic model is on the latent scale, with a
# fixed level-1 variance (pi^2/3 for the logit link). That is scale-free but not on
# an interpretable probability metric. This file adds the response-scale VPC via the
# simulation method of Goldstein, Browne & Rasbash (2002): simulate stratum random
# effects, map them through the inverse link to probabilities, and express the
# between-stratum variance of those probabilities as a share of the total. It is an
# interpretable complement to the latent-scale VPC, not a replacement (it depends on
# the overall prevalence).
#
# For a log-link count model (Poisson, negative binomial) the response scale is the
# count scale, and the partition has a closed form: the exact expressions of Leckie,
# Browne, Goldstein, Merlo & Austin (2020), evaluated at each observation's marginal
# expected count. It involves no simulation and no approximation of the level-1
# variance, which is what separates it from the latent-scale count VPC of summary().

# Number of Gauss-Hermite nodes used for the inner (non-stratum) integral. The
# inverse-link integrand is smooth and bounded in [0, 1], so quadrature converges
# fast; 80 nodes is machine-accurate for the logit/probit/cloglog links and cheap
# (an n_sim x 80 matrix). Kept as a named internal constant so the returned object
# and print method can report the effective inner design.
maihda_vpc_response_inner_nodes <- 80L

# Gauss-Hermite quadrature nodes/weights for the standard-normal measure
# (probabilists' Hermite): returns `nodes` x_k and `weights` w_k with sum(w_k) = 1
# such that sum(w_k * f(x_k)) approximates E[f(Z)] for Z ~ N(0, 1). Built from the
# Golub-Welsch algorithm -- the eigenvalues of the symmetric tridiagonal Jacobi
# matrix (zero diagonal, sqrt(1..n-1) off-diagonals) are the nodes, and the squared
# first components of its (orthonormal) eigenvectors are the weights. Base R only,
# so it adds no statmod dependency. Being deterministic, it removes the inner
# Monte-Carlo error the response-scale VPC previously carried in its fixed 500-draw
# inner sample (so the returned VPC no longer stays conditional on an undocumented
# inner draw count as n_sim grows). Pure and unit-testable.
maihda_gauss_hermite_normal <- function(n = maihda_vpc_response_inner_nodes) {
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 2L) {
    stop("Gauss-Hermite quadrature needs n >= 2 nodes.", call. = FALSE)
  }
  off <- sqrt(seq_len(n - 1L))
  jacobi <- matrix(0, n, n)
  jacobi[cbind(2:n, 1:(n - 1L))] <- off
  jacobi[cbind(1:(n - 1L), 2:n)] <- off
  ev <- eigen(jacobi, symmetric = TRUE)
  nodes <- ev$values
  weights <- ev$vectors[1, ]^2          # mass-1 normal measure => sum(weights) = 1
  ord <- order(nodes)
  list(nodes = nodes[ord], weights = weights[ord])
}

# Families whose response scale is the count scale.
maihda_count_vpc_families <- c("poisson", "negbinomial")

# Exact count-scale variance partition of a two-level random-intercept model with
# the log link (Leckie, Browne, Goldstein, Merlo & Austin 2020: eq. 8 for the
# Poisson, eq. 13-14 for the mean-dispersion negative binomial, NB2). With
# u ~ N(0, var_between), a conditional mean mu = exp(eta + u) and
# Var(y | u) = mu + alpha * mu^2, the law of total variance splits Var(y) into
#   level 2:  Var(E[y | u]) = lambda^2 * (exp(var_between) - 1)
#   level 1:  E[Var(y | u)] = lambda + alpha * lambda^2 * exp(var_between)
# where lambda = exp(eta + var_between / 2) is the MARGINAL expected count, and the
# VPC is level 2 over their sum. `alpha` is 0 for a Poisson model and 1 / theta for
# the negative binomial of lme4 and brms, whose variance is mu + mu^2 / theta.
# Vectorised over `lambda`. Pure, so it is tested against the published values.
#
# The VPC is formed after dividing numerator and denominator by
# lambda * exp(var_between): the same ratio, but finite where lambda^2 or
# exp(var_between) would overflow.
maihda_count_vpc_exact <- function(lambda, var_between, alpha = 0) {
  if (!is.numeric(lambda) || any(!is.finite(lambda)) || any(lambda < 0)) {
    stop("'lambda' must hold finite, non-negative marginal expected counts.",
         call. = FALSE)
  }
  if (!is.numeric(var_between) || length(var_between) != 1L ||
      !is.finite(var_between) || var_between < 0) {
    stop("'var_between' must be a single finite, non-negative variance.",
         call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha < 0) {
    stop("'alpha' must be a single finite, non-negative overdispersion.",
         call. = FALSE)
  }
  between <- -expm1(-var_between)            # 1 - exp(-var_between)
  list(
    lambda = lambda,
    level2 = lambda^2 * expm1(var_between),
    level1 = lambda + alpha * lambda^2 * exp(var_between),
    vpc = lambda * between / (lambda * between + exp(-var_between) + alpha * lambda)
  )
}

# Quantiles of a per-observation quantity under prior weights read as frequencies,
# the companion of maihda_weighted_obs_mean(): with whole-number weights it returns
# stats::quantile(rep(v, w), probs) -- the default (type 7) quantile of the
# duplicated-row data -- without building that vector. The k-th order statistic of
# the expanded sample is the first sorted value whose cumulative weight reaches k,
# and type 7 interpolates between the order statistics either side of
# 1 + (N - 1) * p with N = sum(w). Reduces to stats::quantile(v, probs) when the
# weights are absent, length-mismatched or all 1, on the same rule as the mean.
maihda_weighted_obs_quantile <- function(v, w = NULL, probs = c(0.25, 0.5, 0.75)) {
  v <- as.numeric(v)
  if (is.null(w) || length(w) != length(v)) {
    return(stats::quantile(v, probs = probs, names = FALSE, na.rm = TRUE))
  }
  w <- as.numeric(w)
  ok <- is.finite(v) & is.finite(w) & w > 0
  if (!any(ok) || isTRUE(all(abs(w[ok] - 1) < sqrt(.Machine$double.eps)))) {
    return(stats::quantile(v, probs = probs, names = FALSE, na.rm = TRUE))
  }
  ord <- order(v[ok])
  sorted <- v[ok][ord]
  cum_w <- cumsum(w[ok][ord])
  n_w <- cum_w[length(cum_w)]
  order_stat <- function(k) {
    sorted[pmin(findInterval(k, cum_w, left.open = TRUE) + 1L, length(sorted))]
  }
  h <- 1 + (n_w - 1) * probs
  lo <- floor(h)
  x_lo <- order_stat(pmax(lo, 1))
  x_hi <- order_stat(pmin(lo + 1, n_w))
  x_lo + (h - lo) * (x_hi - x_lo)
}

# Why the count-scale VPC is not available for `model`, as a complete error
# message, or NULL when it is. maihda_vpc_response() raises it and
# summary(response_vpc = TRUE) skips the fit on it, as it skips any other family or
# engine the response-scale VPC does not cover -- one function, so the two cannot
# disagree on what is covered.
#
# The closed form is the two-level random-intercept result: a single normal random
# intercept on the log scale. A second random effect of any kind -- a context, an
# explicit (1 | site), the dimension intercepts of a crossed-dimensions fit, a
# random slope -- changes the marginal moments, and with crossed effects the count
# scale has no additive partition to read one stratum share from (Leckie et al.
# derive nested levels and random coefficients, not cross-classified models).
maihda_count_vpc_refusal <- function(model) {
  if (!identical(model$engine, "lme4")) {
    return("maihda_vpc_response() is currently implemented for the lme4 engine.")
  }
  fam <- maihda_model_family_name(model)
  if (!isTRUE(fam %in% maihda_count_vpc_families)) {
    return(paste0("The count-scale VPC is defined for Poisson and negative-binomial ",
                  "models; this model uses family = '", fam, "'."))
  }
  link <- maihda_model_link_name(model)
  if (!identical(link, "log")) {
    return(paste0("The count-scale VPC is defined for the log link; this model uses ",
                  "link = '", link, "'."))
  }

  only <- paste0("The count-scale VPC is implemented for a model whose only random ",
                 "effect is the stratum random intercept; ")
  if (!is.null(model$longitudinal_info)) {
    return(paste0(only, "a longitudinal MAIHDA has random slopes on time."))
  }
  if (!is.null(model$cc_info)) {
    return(paste0(only, "a crossed-dimensions MAIHDA has a random intercept for ",
                  "each dimension as well."))
  }
  if (!is.null(model$context_info)) {
    return(paste0(only, "this model also has the context random intercept(s) ",
                  paste(model$context_info$context_vars, collapse = ", "), "."))
  }
  slopes <- tryCatch({
    maihda_validate_intercept_only_random_effects_lme4(
      model$model, context = "The count-scale VPC")
    NULL
  }, error = function(e) conditionMessage(e))
  if (!is.null(slopes)) {
    return(slopes)
  }
  groups <- names(lme4::VarCorr(model$model))
  other <- setdiff(groups, "stratum")
  if (length(other) > 0) {
    return(paste0(only, "this model also has the random effect(s) ",
                  paste(other, collapse = ", "), "."))
  }
  if (!"stratum" %in% groups) {
    return("No 'stratum' random-effect variance found in the model.")
  }
  NULL
}

# The count-scale VPC of a Poisson or negative-binomial fit. lambda_i is the
# marginal expected count the latent-scale count VPC is built on
# (maihda_count_marginal_mu_lme4(): the fixed-part linear predictor, offset
# included, plus half the random-intercept variance), so both VPCs of one fit refer
# to the same counts. With covariates or an offset lambda_i varies, and so does the
# VPC: it is evaluated for every observation and summarised by its mean, as Leckie
# et al. (2020) report it, with the median, quartiles and range alongside. Prior
# weights enter those summaries as frequencies, which is how the latent-scale path
# reads them (maihda_weighted_obs_mean()), so a weighted fit returns what the
# duplicated-row data would.
maihda_vpc_response_count <- function(model) {
  why <- maihda_count_vpc_refusal(model)
  if (!is.null(why)) {
    stop(why, call. = FALSE)
  }
  fitted_model <- model$model
  fam <- maihda_model_family_name(model)
  var_between <- maihda_stratum_variance_lme4(fitted_model)
  # Var(y | u) = mu + mu^2 / theta, so alpha = 1 / theta; a Poisson model is the
  # theta -> Inf limit.
  theta <- if (identical(fam, "negbinomial")) {
    maihda_negbin_theta_lme4(fitted_model)
  } else {
    Inf
  }
  lambda <- maihda_count_marginal_mu_lme4(fitted_model)
  vpc <- maihda_count_vpc_exact(lambda, var_between, alpha = 1 / theta)$vpc
  w <- maihda_fit_prior_weights(fitted_model)
  quants <- maihda_weighted_obs_quantile(vpc, w, probs = c(0.25, 0.5, 0.75))

  structure(
    list(estimate = maihda_weighted_obs_mean(vpc, w),
         scale = "response", method = "closed-form", family = fam,
         median = quants[2], quartiles = quants[c(1, 3)],
         range = range(vpc),
         lambda = maihda_weighted_obs_mean(lambda, w),
         lambda_range = range(lambda),
         theta = theta,
         var_between = var_between, var_other = 0),
    class = "maihda_vpc_response"
  )
}

#' Response-scale VPC for a binomial or count MAIHDA model
#'
#' @description
#' Computes the variance partition coefficient on the scale of the observed
#' response: the probability scale for a binomial MAIHDA model, the count scale for
#' a Poisson or negative-binomial one.
#'
#' For a \strong{binomial} model it uses the simulation method of Goldstein, Browne &
#' Rasbash (2002). Stratum random effects \eqn{u \sim N(0, \sigma^2_u)} are
#' simulated and converted to predicted probabilities \eqn{p = g^{-1}(\eta + u)}
#' (with \eqn{\eta} the fixed-part linear predictor); the VPC is then the
#' between-stratum variance of \eqn{p} as a share of the total
#' (between + the binomial within-stratum variance \eqn{\overline{p(1-p)}}).
#' Unlike the latent-scale VPC (fixed level-1 variance \eqn{\pi^2/3} for the logit),
#' the response-scale VPC depends on the overall outcome prevalence, so report it as
#' a complement to -- not a replacement for -- the latent-scale value.
#'
#' For a \strong{Poisson or negative-binomial} model with the log link it evaluates
#' the exact expressions of Leckie, Browne, Goldstein, Merlo & Austin (2020). With
#' \eqn{\lambda = \exp(\eta + \sigma^2_u / 2)} the marginal expected count, the
#' between-stratum variance of the counts is
#' \eqn{\lambda^2 (\exp(\sigma^2_u) - 1)} and their within-stratum variance is
#' \eqn{\lambda + \alpha \lambda^2 \exp(\sigma^2_u)}, where \eqn{\alpha = 1 / \theta}
#' is the negative-binomial overdispersion (0 for a Poisson model); the VPC is the
#' between-stratum variance as a share of their sum. It is the share of the variance
#' of the counts themselves that lies between strata, and it involves no
#' approximation of the level-1 variance. The latent-scale VPC of
#' \code{\link{summary.maihda_model}} refers to the log scale instead and depends on
#' the \code{count_approximation} of \code{\link{fit_maihda}}, so the two answer
#' different questions. Report the count-scale VPC when the partition should
#' describe the observed counts, or when the expected count is low and the
#' latent-scale approximations disagree.
#'
#' @details
#' \strong{Binomial models.}
#' The fixed part \eqn{\eta} is collapsed to a single value -- the mean linear
#' predictor \eqn{\bar\eta} over the analytic sample -- before the random effect is
#' simulated around it. The result is therefore a VPC \emph{evaluated at the mean
#' covariate profile} (a conditional-at-mean estimate), not one marginalised over the
#' empirical covariate distribution. For the canonical strata-only (null) model
#' \eqn{\eta} is constant (the intercept), so the two coincide and the value is
#' exact. For an \emph{adjusted} model (one with covariates) they can differ, because
#' the inverse link is nonlinear and \eqn{g^{-1}(\bar\eta) \neq \overline{g^{-1}(\eta)}}:
#' read the response-scale VPC from the null model, or interpret an adjusted value as
#' conditional on the average covariate profile rather than as a covariate-averaged
#' (marginal) VPC.
#'
#' The method is binomial-link agnostic: it maps the simulated stratum effects through
#' whichever inverse link the model uses (logit, probit, cloglog, ...), so a non-logit
#' binomial fit is computed on its own scale rather than rejected. Only the family is
#' required to be binomial.
#'
#' When the model carries random intercepts \emph{beyond} the stratum (a
#' contextual \code{(1 | school)} or an explicit \code{(1 | site)}), the
#' simulation integrates over them: the reported estimate is the
#' \emph{stratum share} \eqn{Var(E[p \mid u_{stratum}])} of the total
#' response-scale variance, where the total includes the variation the other
#' random effects induce in \eqn{p} plus the binomial within-variance.
#' Simulating the stratum effect alone would overstate the stratum share. The
#' inner integral over the combined non-stratum effect is evaluated by
#' \emph{deterministic Gauss-Hermite quadrature}, not an inner Monte-Carlo
#' sample, so it contributes no simulation error of its own and increasing
#' \code{n_sim} converges the whole estimate (the outer stratum draw is then the
#' only Monte-Carlo dimension). The number of quadrature nodes used is reported in
#' \code{inner_nodes}.
#'
#' \strong{Count models.} \eqn{\lambda} is computed for each observation of the
#' analytic sample from its fixed-part linear predictor, any offset included, so
#' with covariates or an offset the VPC differs between observations. The
#' \code{estimate} is the mean of the observation-level VPCs; their median,
#' quartiles and range are returned with it, as are the mean and range of
#' \eqn{\lambda}. Prior \code{weights} of the fit enter these summaries as
#' frequencies. A null model without an offset has a single \eqn{\lambda} and a
#' single VPC. \code{n_sim} and \code{seed} are not used.
#'
#' The expressions are those of a model whose only random effect is the stratum
#' random intercept. A count model with a \code{context} or any other additional
#' random effect, a crossed-dimensions fit, a longitudinal fit, and a link other
#' than the log are refused.
#'
#' @param model A binomial, Poisson or negative-binomial \code{maihda_model} (lme4
#'   engine) from \code{\link{fit_maihda}}.
#' @param n_sim Number of Monte-Carlo draws of the stratum random effect (>= 100)
#'   for a binomial model. Default 10000.
#' @param seed Optional integer seed for the binomial simulation.
#'
#' @return An object of class \code{maihda_vpc_response}: a list with
#'   \code{estimate}, \code{scale = "response"}, \code{method}
#'   (\code{"simulation"} for a binomial model, \code{"closed-form"} for a count
#'   model), \code{var_between} (the latent-scale between-stratum variance) and
#'   \code{var_other} (the summed latent-scale variance of any non-stratum random
#'   intercepts, 0 when there are none).
#'
#'   For a binomial model it also holds \code{n_sim}, \code{lp_fixed} (the mean
#'   fixed-part linear predictor), \code{inner_method} (\code{"gauss-hermite"} when
#'   non-stratum effects are integrated out, else \code{"none"}) and
#'   \code{inner_nodes} (the number of Gauss-Hermite quadrature nodes used for that
#'   inner integral, \code{NA} when there are no non-stratum effects).
#'
#'   For a count model it also holds \code{family}, the \code{median},
#'   \code{quartiles} and \code{range} of the observation-level VPCs (all equal to
#'   \code{estimate} when \eqn{\lambda} is the same for every observation),
#'   \code{lambda} (the mean marginal expected count), \code{lambda_range} and
#'   \code{theta} (the negative-binomial dispersion, \code{Inf} for a Poisson
#'   model).
#'
#' @references
#' Goldstein, H., Browne, W., & Rasbash, J. (2002). Partitioning variation in
#' multilevel models. \emph{Understanding Statistics}, 1(4), 223-231.
#'
#' Leckie, G., Browne, W. J., Goldstein, H., Merlo, J., & Austin, P. C. (2020).
#' Partitioning variation in multilevel models for count data.
#' \emph{Psychological Methods}, 25(6), 787-801. \doi{10.1037/met0000265}
#'
#' @seealso \code{\link{maihda_discriminatory_accuracy}}, \code{\link{summary.maihda_model}}
#'
#' @examples
#' \donttest{
#' strata <- make_strata(maihda_health_data, vars = c("Gender", "Race"))
#' d <- maihda_health_data
#' d$stratum <- strata$data$stratum
#' m <- fit_maihda(Obese ~ (1 | stratum), data = d, family = "binomial")
#' maihda_vpc_response(m, seed = 1)
#'
#' # A count outcome: the VPC on the count scale, in closed form
#' set.seed(1)
#' u <- rnorm(max(d$stratum), sd = 0.4)
#' d$visits <- rpois(nrow(d), exp(0.5 + u[d$stratum]))
#' m_count <- fit_maihda(visits ~ (1 | stratum), data = d, family = "poisson")
#' maihda_vpc_response(m_count)
#' }
#'
#' @export
maihda_vpc_response <- function(model, n_sim = 10000, seed = NULL) {
  if (!inherits(model, "maihda_model")) {
    stop("'model' must be a maihda_model object from fit_maihda().", call. = FALSE)
  }
  if (!identical(model$engine, "lme4")) {
    stop("maihda_vpc_response() is currently implemented for the lme4 engine.",
         call. = FALSE)
  }
  fam <- maihda_model_family_name(model)
  # A count model has a closed form on its own response scale. n_sim and seed
  # belong to the binomial simulation and are not read on this path, which draws
  # no random number.
  if (isTRUE(fam %in% maihda_count_vpc_families)) {
    return(maihda_vpc_response_count(model))
  }
  if (!identical(fam, "binomial")) {
    stop("The response-scale VPC is defined for binomial MAIHDA models (any ",
         "binomial link) and, on the count scale, for Poisson and negative-binomial ",
         "models with the log link; this model uses family = '", fam, "'.",
         call. = FALSE)
  }
  if (!is.numeric(n_sim) || length(n_sim) != 1 || !is.finite(n_sim) ||
      n_sim < 100 || n_sim != floor(n_sim)) {
    stop("'n_sim' must be a single whole number >= 100.", call. = FALSE)
  }
  # rnorm() silently truncates a fractional count (rnorm(100.5) draws 100) while the
  # recorded n_sim would keep the fraction; cast so the draw count and the report agree.
  n_sim <- as.integer(n_sim)

  # The simulation collapses every random effect to a single intercept variance:
  # the stratum draw uses var_between and the non-stratum effects are summed into
  # one normal (var_other, below). That is correct only for intercept-only random
  # effects. A random slope -- (x | stratum) or (x | site) -- carries a
  # slope variance, an intercept-slope covariance, and a row-specific design value
  # that this method does not integrate over, so its between-stratum variance is
  # not a scalar. Reject such a model explicitly, with the same message the scalar
  # summary/VPC path uses, rather than reading only the intercept variance and
  # silently returning a wrong (or NA) response-scale VPC. Longitudinal growth fits
  # (which carry random slopes by construction) route through the time-varying path
  # and never reach maihda_vpc_response(). lme4 is the only engine here (checked above).
  maihda_validate_intercept_only_random_effects_lme4(
    model$model,
    context = "The response-scale VPC"
  )

  # The random-intercept variances by grouping factor, read once: the
  # stratum-structure groups define the between-stratum (intersectional) variance
  # and every other group is a non-stratum effect integrated into the total
  # (var_other, below).
  all_vars <- tryCatch(maihda_random_variances_lme4(model$model),
                       error = function(e) NULL)

  # Grouping factors that jointly make up the intersectional between-stratum
  # variance. For the canonical fit this is the single "stratum" intercept. For a
  # crossed-dimensions fit ($cc_info) the additive dimension REs and the
  # interaction RE are independent crossed effects that SUM at the intersection
  # level, so all of them are between-stratum -- exactly as the MOR
  # (maihda_mor_between_variance()) and the crossed-dimensions VPC
  # (maihda_cc_partition(): between = additive + interaction) read it. Treating
  # only the interaction ("stratum") component as between-stratum -- and the
  # additive dimension REs as non-stratum (var_other) -- understates the
  # response-scale stratum share by orders of magnitude.
  stratum_groups <- if (is.null(model$cc_info)) {
    "stratum"
  } else {
    unique(c(model$cc_info$interaction_group, unname(model$cc_info$dim_groups)))
  }

  var_between <- tryCatch(
    if (is.null(model$cc_info)) {
      extract_between_variance(model)
    } else {
      missing_re <- setdiff(stratum_groups, names(all_vars))
      if (length(missing_re) > 0) {
        stop("The response-scale VPC is missing the crossed random effect(s): ",
             paste(missing_re, collapse = ", "), ".", call. = FALSE)
      }
      sum(all_vars[stratum_groups])
    },
    error = function(e) NA_real_)
  if (!is.numeric(var_between) || length(var_between) != 1 ||
      !is.finite(var_between) || var_between < 0) {
    return(structure(
      list(estimate = NA_real_, scale = "response", method = "simulation",
           n_sim = n_sim, var_between = var_between, var_other = NA_real_,
           lp_fixed = NA_real_, inner_method = NA_character_,
           inner_nodes = NA_integer_),
      class = "maihda_vpc_response"))
  }

  # Non-stratum random-intercept variances (a contextual (1 | school) or an
  # explicit (1 | site)). The simulation must integrate over them: they widen
  # the distribution of p -- entering the total variance -- without separating
  # strata, so simulating the stratum effect alone overstates the stratum
  # share. Independent random intercepts sum on the link scale, so their
  # combined effect is one normal with the summed variance. Every stratum-
  # structure group (for a crossed fit, the dimension REs and the interaction)
  # is excluded here -- it already lives in var_between above.
  var_other <- 0
  if (!is.null(all_vars)) {
    other <- all_vars[setdiff(names(all_vars), stratum_groups)]
    var_other <- sum(other[is.finite(other)])
  }

  fitted_model <- model$model
  linkinv <- stats::family(fitted_model)$linkinv
  # Fixed-part linear predictor (random effects excluded), collapsed to its sample
  # mean. For the canonical null / strata-only model this is exactly the intercept,
  # so the VPC below is exact; with covariates this is the mean fixed linear
  # predictor, so the VPC is a conditional-at-mean estimate (evaluated at the average
  # covariate profile) rather than one integrated over the covariate distribution --
  # see the @details section of the function documentation.
  lp_fixed <- mean(maihda_fit_rows_predict(fitted_model, re.form = NA, type = "link"),
                   na.rm = TRUE)

  if (!is.null(seed)) {
    # Keep reproducibility local: snapshot the caller's RNG state and restore it
    # on exit, so passing seed= does not silently reseed the session and perturb
    # the user's subsequent random draws. (Base-R equivalent of withr::with_seed;
    # withr is not a package dependency.)
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
      on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
    } else {
      # RNG was uninitialised before this call; remove the state we introduce so
      # the session is left exactly as we found it.
      on.exit(
        if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
          rm(".Random.seed", envir = globalenv())
        },
        add = TRUE
      )
    }
    set.seed(seed)
  }
  u <- stats::rnorm(n_sim, mean = 0, sd = sqrt(var_between))
  inner_nodes <- NA_integer_
  if (var_other > 0) {
    # Nested integration: E[p | u_stratum] over the combined non-stratum effect.
    # The inner 1-D Gaussian integral is done by DETERMINISTIC Gauss-Hermite
    # quadrature rather than a fixed inner Monte-Carlo sample, so it carries no
    # simulation error of its own: increasing n_sim now converges the whole
    # estimate (the outer stratum draw is the only remaining Monte-Carlo
    # dimension, and n_sim controls it). The stratum share is Var(E[p | u_stratum])
    # over the TOTAL response-scale variance: Var(p) across strata and the other
    # effects plus the binomial within-variance E[p(1-p)]. All three moments are
    # taken against the same quadrature weights so they stay mutually consistent.
    gh <- maihda_gauss_hermite_normal(maihda_vpc_response_inner_nodes)
    inner_nodes <- length(gh$nodes)
    z <- gh$nodes * sqrt(var_other)          # quadrature nodes for N(0, var_other)
    wq <- gh$weights                         # quadrature weights (sum to 1)
    p_mat <- linkinv(outer(u, z, `+`) + lp_fixed)          # n_sim x inner_nodes
    m_stratum <- as.vector(p_mat %*% wq)                   # E[p | u_i], exact inner
    v_between <- stats::var(m_stratum)
    # Var(p) over the joint distribution via the LAW OF TOTAL VARIANCE, so the
    # between piece of the denominator is the SAME estimator as the numerator.
    # Reading it instead off pooled population moments (E[p^2] - E[p]^2) divides an
    # n-1 sample variance by a total whose between part uses n, inflating the VPC by
    # exactly n/(n-1) -- 1.01% at the allowed minimum n_sim = 100. The inner term is
    # a mean of per-draw quantities, already on that same footing. This form is also
    # exactly continuous with the var_other == 0 branch below, where it vanishes.
    v_cond <- as.vector((p_mat * p_mat) %*% wq) - m_stratum^2   # Var(p | u_i)
    v_total_p <- v_between + mean(pmax(v_cond, 0))         # Var(p) over the joint
    v_within <- mean(as.vector((p_mat * (1 - p_mat)) %*% wq))  # E[p(1-p)]
    vpc <- v_between / (v_total_p + v_within)
  } else {
    p <- linkinv(lp_fixed + u)
    v_between <- stats::var(p)
    v_within <- mean(p * (1 - p))
    vpc <- v_between / (v_between + v_within)
  }

  structure(
    list(estimate = vpc, scale = "response", method = "simulation",
         n_sim = n_sim, var_between = var_between, var_other = var_other,
         lp_fixed = lp_fixed,
         inner_method = if (var_other > 0) "gauss-hermite" else "none",
         inner_nodes = inner_nodes),
    class = "maihda_vpc_response"
  )
}

#' @export
print.maihda_vpc_response <- function(x, ...) {
  pal <- maihda_palette()
  if (identical(x$method, "closed-form")) {
    cat(pal$bold("Count-scale VPC (closed form, Leckie et al. 2020)"), "\n", sep = "")
    cat(sprintf("  VPC: %s\n",
                if (is.finite(x$estimate)) pal$accent(sprintf("%.4f", x$estimate)) else "NA"))
    # One lambda (a null model without an offset) gives one VPC; otherwise the
    # estimate is a mean over observations and the spread says how far it varies.
    if (isTRUE(all.equal(x$lambda_range[1], x$lambda_range[2]))) {
      cat(pal$muted(sprintf(
        "  Marginal expected count %s; between-stratum variance %.4f (log scale).\n",
        format(x$lambda, digits = 3), x$var_between)))
    } else {
      cat(pal$muted(sprintf(paste0(
        "  Mean of the observation-level VPCs (median %.4f, quartiles %.4f and\n",
        "  %.4f, range %.4f to %.4f).\n",
        "  Marginal expected count: mean %s, range %s to %s. Between-stratum\n",
        "  variance %.4f (log scale).\n"),
        x$median, x$quartiles[1], x$quartiles[2], x$range[1], x$range[2],
        format(x$lambda, digits = 3), format(x$lambda_range[1], digits = 3),
        format(x$lambda_range[2], digits = 3), x$var_between)))
    }
    if (is.finite(x$theta)) {
      cat(pal$muted(sprintf("  Negative-binomial theta %s.\n",
                            format(x$theta, digits = 3))))
    }
    return(invisible(x))
  }
  cat(pal$bold("Response-scale VPC (simulation method)"), "\n", sep = "")
  cat(sprintf("  VPC: %s\n",
              if (is.finite(x$estimate)) pal$accent(sprintf("%.4f", x$estimate)) else "NA"))
  cat(pal$muted(sprintf("  %d simulated stratum effects; between-stratum variance %.4f (latent scale).\n",
              x$n_sim, x$var_between)))
  if (is.finite(x$var_other) && x$var_other > 0) {
    nodes_txt <- if (!is.null(x$inner_nodes) && is.finite(x$inner_nodes)) {
      sprintf(" by %d-node Gauss-Hermite quadrature", x$inner_nodes)
    } else {
      ""
    }
    cat(pal$muted(sprintf(paste0(
      "  Integrated over non-stratum random effects (latent variance %.4f)%s:\n",
      "  the estimate is the stratum share of the total response-scale variance.\n"),
      x$var_other, nodes_txt)))
  }
  invisible(x)
}
