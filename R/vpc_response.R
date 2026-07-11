# Response-scale (observed / probability-scale) VPC for binomial MAIHDA models.
#
# The package's default VPC for a logistic model is on the latent scale, with a
# fixed level-1 variance (pi^2/3 for the logit link). That is scale-free but not on
# an interpretable probability metric. This file adds the response-scale VPC via the
# simulation method of Goldstein, Browne & Rasbash (2002): simulate stratum random
# effects, map them through the inverse link to probabilities, and express the
# between-stratum variance of those probabilities as a share of the total. It is an
# interpretable complement to the latent-scale VPC, not a replacement (it depends on
# the overall prevalence).

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

#' Response-scale VPC for a binomial MAIHDA model
#'
#' @description
#' Computes the variance partition coefficient on the response (probability) scale
#' for a binomial MAIHDA model, using the simulation method of Goldstein, Browne &
#' Rasbash (2002). Stratum random effects \eqn{u \sim N(0, \sigma^2_u)} are
#' simulated and converted to predicted probabilities \eqn{p = g^{-1}(\eta + u)}
#' (with \eqn{\eta} the fixed-part linear predictor); the VPC is then the
#' between-stratum variance of \eqn{p} as a share of the total
#' (between + the binomial within-stratum variance \eqn{\overline{p(1-p)}}).
#'
#' Unlike the latent-scale VPC (fixed level-1 variance \eqn{\pi^2/3} for the logit),
#' the response-scale VPC depends on the overall outcome prevalence, so report it as
#' a complement to -- not a replacement for -- the latent-scale value.
#'
#' @details
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
#' @param model A binomial \code{maihda_model} (lme4 engine) from
#'   \code{\link{fit_maihda}}.
#' @param n_sim Number of Monte-Carlo draws of the stratum random effect (>= 100).
#'   Default 10000.
#' @param seed Optional integer seed for reproducibility.
#'
#' @return An object of class \code{maihda_vpc_response}: a list with
#'   \code{estimate}, \code{scale = "response"}, \code{method = "simulation"},
#'   \code{n_sim}, \code{var_between} (the latent-scale between-stratum variance),
#'   \code{var_other} (the summed latent-scale variance of any non-stratum random
#'   intercepts, 0 when there are none), \code{lp_fixed} (the mean fixed-part
#'   linear predictor), \code{inner_method} (\code{"gauss-hermite"} when non-stratum
#'   effects are integrated out, else \code{"none"}) and \code{inner_nodes} (the
#'   number of Gauss-Hermite quadrature nodes used for that inner integral, \code{NA}
#'   when there are no non-stratum effects).
#'
#' @references
#' Goldstein, H., Browne, W., & Rasbash, J. (2002). Partitioning variation in
#' multilevel models. \emph{Understanding Statistics}, 1(4), 223-231.
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
  if (!identical(fam, "binomial")) {
    stop("The response-scale VPC is only defined for binomial MAIHDA models ",
         "(any binomial link); this model uses family = '", fam, "'.", call. = FALSE)
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

  var_between <- tryCatch(extract_between_variance(model), error = function(e) NA_real_)
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
  # combined effect is one normal with the summed variance.
  all_vars <- tryCatch(maihda_random_variances_lme4(model$model),
                       error = function(e) NULL)
  var_other <- 0
  if (!is.null(all_vars)) {
    other <- all_vars[setdiff(names(all_vars), "stratum")]
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
  lp_fixed <- mean(stats::predict(fitted_model, re.form = NA, type = "link"), na.rm = TRUE)

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
    e_p <- mean(m_stratum)                                 # E[p]
    e_p2 <- mean(as.vector((p_mat * p_mat) %*% wq))        # E[p^2]
    v_total_p <- max(e_p2 - e_p^2, 0)                      # Var(p) over the joint
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
