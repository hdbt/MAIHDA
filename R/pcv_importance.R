# Order-invariant PCV attribution (Shapley / dominance) -- issue #63.
#
# The engine room: every attribution here is computed from ONE value function
#
#   v(S) = (V0 - V(S)) / V0
#
# where V0 is the null model's between-stratum variance and V(S) the
# between-stratum variance after adding the variable subset S as fixed effects
# -- i.e. v(S) is exactly the Total_PCV that stepwise_pcv() reports when its
# path ends at S. Subsets are encoded as bitmasks over seq_along(vars)
# (bit i-1 set <=> vars[i] in S), each subset model is fit AT MOST ONCE and
# cached, and the attribution methods are pure functions of the resulting
# v-lookup -- which is what lets the parametric bootstrap re-run the identical
# combinatorial computation on refit variances.

# Bit i-1 set <=> vars[i] in the subset. Vectorised over the bit positions.
maihda_mask_bits <- function(mask, k) {
  which(bitwAnd(mask, bitwShiftL(1L, seq_len(k) - 1L)) != 0L)
}

# Subset sizes (popcounts) for all masks 0..n_masks-1, as a vector indexed by
# mask + 1. One pass per bit, so O(k * 2^k) integer ops -- trivial at the k
# where the exact attribution is feasible at all.
maihda_mask_sizes <- function(n_masks, k) {
  sizes <- integer(n_masks)
  masks <- seq_len(n_masks) - 1L
  for (b in seq_len(k)) {
    has <- bitwAnd(masks, bitwShiftL(1L, b - 1L)) != 0L
    sizes[has] <- sizes[has] + 1L
  }
  sizes
}

# Exact Shapley values from a complete v-table (numeric vector indexed mask + 1,
# length 2^k, v[1] = v(empty) = 0):
#   phi_i = sum over S not containing i of |S|! (k-|S|-1)! / k! * (v(S u {i}) - v(S)).
# By the telescoping/efficiency property sum(phi) = v(full) exactly.
maihda_shapley_exact_phi <- function(v, k) {
  n_masks <- length(v)
  masks <- seq_len(n_masks) - 1L
  sizes <- maihda_mask_sizes(n_masks, k)
  fact <- factorial(0:k)
  # w[s + 1] = s! (k - s - 1)! / k! for s = 0..k-1
  w <- fact[seq_len(k)] * fact[seq(k, 1)] / fact[k + 1]
  phi <- numeric(k)
  for (i in seq_len(k)) {
    bit <- bitwShiftL(1L, i - 1L)
    without_i <- masks[bitwAnd(masks, bit) == 0L]
    s <- sizes[without_i + 1L]
    phi[i] <- sum(w[s + 1L] * (v[bitwOr(without_i, bit) + 1L] - v[without_i + 1L]))
  }
  phi
}

# Budescu dominance tables from a complete v-table.
#   conditional[i, s+1]: variable i's average marginal contribution
#     v(S u {i}) - v(S) over all subsets S (i not in S) of size s -- one column
#     per adjustment-set size 0..k-1.
#   general: rowMeans(conditional). Averaging within each size and then across
#     sizes reproduces the Shapley weights exactly, so general == the Shapley
#     values (LMG); the identity is pinned by tests.
#   complete[i, j]: TRUE when i's marginal contribution STRICTLY exceeds j's for
#     every subset containing neither (Budescu 1993's complete dominance; both
#     v(S u {i}) and v(S u {j}) subtract the same v(S), so comparing the raw
#     v values is equivalent). A single tie or reversal makes it FALSE;
#     diagonal NA.
maihda_dominance_tables <- function(v, k, var_names) {
  n_masks <- length(v)
  masks <- seq_len(n_masks) - 1L
  sizes <- maihda_mask_sizes(n_masks, k)
  conditional <- matrix(NA_real_, nrow = k, ncol = k,
                        dimnames = list(var_names, paste0("size_", 0:(k - 1))))
  for (i in seq_len(k)) {
    bit <- bitwShiftL(1L, i - 1L)
    without_i <- masks[bitwAnd(masks, bit) == 0L]
    marg <- v[bitwOr(without_i, bit) + 1L] - v[without_i + 1L]
    s <- sizes[without_i + 1L]
    conditional[i, ] <- vapply(0:(k - 1), function(sz) mean(marg[s == sz]),
                               numeric(1))
  }
  complete <- matrix(NA, nrow = k, ncol = k,
                     dimnames = list(var_names, var_names))
  for (i in seq_len(k)) {
    for (j in seq_len(k)) {
      if (i == j) next
      bi <- bitwShiftL(1L, i - 1L)
      bj <- bitwShiftL(1L, j - 1L)
      neither <- masks[bitwAnd(masks, bitwOr(bi, bj)) == 0L]
      complete[i, j] <- all(v[bitwOr(neither, bi) + 1L] >
                              v[bitwOr(neither, bj) + 1L])
    }
  }
  list(conditional = conditional, general = rowMeans(conditional),
       complete = complete)
}

#' Order-Invariant PCV Attribution Across Predictors (Shapley / Dominance)
#'
#' @description
#' Apportions the between-stratum variance reduction (the PCV) among a set of
#' predictors. Where \code{\link{stepwise_pcv}} adds the variables one at a
#' time -- so each variable's contribution depends on its entry order --
#' \code{pcv_importance()} treats the PCV as a value function over variable
#' subsets and attributes it \emph{fairly}: the flagship \code{"shapley"}
#' method averages each variable's marginal PCV over every possible entry
#' order, and \code{"dominance"} adds Budescu's pairwise dominance detail (its
#' general dominance weights coincide with the Shapley values). All methods
#' satisfy the \strong{efficiency} identity: the contributions sum exactly to
#' the full model's total PCV. (The order-dependent \code{"sequential"} method
#' is \strong{deprecated}; use \code{\link{stepwise_pcv}} for the sequential
#' path -- see the \code{method} argument.)
#'
#' Two attribution targets are useful in a MAIHDA. Passing the \emph{stratum
#' dimensions} (e.g. \code{c("gender", "race", "education")}) splits the
#' additive share -- the canonical null-to-adjusted PCV -- fairly across the
#' dimensions: which social dimension drives the additive between-stratum
#' inequality. Passing \emph{individual-level covariates} gives an
#' order-invariant counterpart to \code{stepwise_pcv()}: how much each
#' covariate explains of the between-stratum variance (subject to the
#' latent-scale caveat below for non-Gaussian families).
#'
#' @param data Data frame with observations. Ensure \code{\link{make_strata}}
#'   was run first so the \code{stratum} variable exists.
#' @param outcome Character string; the dependent variable.
#' @param vars Character vector of predictors (stratum dimensions and/or
#'   covariates) among which the PCV is apportioned. Order does not affect the
#'   \code{"shapley"} and \code{"dominance"} results; it defines the path for
#'   \code{"sequential"}.
#' @param method Attribution method: \code{"shapley"} (default; order-invariant
#'   Shapley values) or \code{"dominance"} (Budescu dominance analysis: general
#'   dominance -- equal to the Shapley values -- plus the conditional and
#'   complete dominance detail). \code{"sequential"} (the order-dependent
#'   one-at-a-time path) is \strong{deprecated} and will be removed in a future
#'   release: it still runs but warns, and \code{\link{stepwise_pcv}} is the
#'   supported sequential decomposition -- it additionally reports the
#'   step-specific \code{Step_PCV} and, for a binary outcome, the
#'   discriminatory-accuracy trajectory (AUC, MOR).
#' @param approx For \code{method = "shapley"} only: \code{"exact"} fits all
#'   \eqn{2^k - 1} non-empty variable subsets (plus the null); \code{"montecarlo"}
#'   samples \code{n_perm} random entry orders instead. When \code{approx} is
#'   not supplied it defaults to \code{"exact"} for up to 10 variables and
#'   switches to \code{"montecarlo"} (with a message) beyond.
#'   \code{"sequential"} needs only the \eqn{k} path models and ignores
#'   \code{approx}; \code{"dominance"} requires every subset, so it errors if
#'   \code{approx = "montecarlo"} is requested.
#' @param n_perm Number of random permutations for
#'   \code{approx = "montecarlo"}. Default 2000. Distinct subset fits are
#'   cached, so the number of \emph{model fits} is bounded by
#'   \code{min(2^k, n_perm * (k - 1) + k + 1)}; set a seed
#'   (\code{set.seed()}) for reproducible sampling.
#' @param engine Modeling engine ("lme4", "brms", "wemix", or "ordinal");
#'   default "lme4". Resolved exactly as in \code{\link{stepwise_pcv}}:
#'   switches to "wemix" when \code{sampling_weights} is supplied and to
#'   "ordinal" for an ordinal family or ordered-factor outcome. Exact
#'   attribution with \code{engine = "brms"} refits a Stan model for
#'   \emph{every} subset and is strongly discouraged beyond a handful of
#'   variables (see Details).
#' @param family Error distribution and link function. Default "gaussian";
#'   a binary or ordered-factor outcome is auto-detected when \code{family} is
#'   left unspecified, mirroring \code{\link{stepwise_pcv}}.
#' @param context Optional higher-level context column(s) (e.g.
#'   \code{"school"}), forwarded to every subset fit so each model carries the
#'   crossed \code{(1 | context)} intercept alongside the stratum effect; the
#'   attribution is then of the between-stratum PCV \strong{net of} the
#'   context, exactly as in \code{\link{stepwise_pcv}}. \code{lme4}/\code{brms}
#'   engines only.
#' @param sampling_weights Optional name of a sampling-weight column for
#'   design-weighted fits; see \code{\link{fit_maihda}}. The weight column
#'   joins the complete-case filter so every subset fit uses the same analytic
#'   sample.
#' @param bootstrap Logical; compute parametric-bootstrap confidence intervals
#'   for each contribution by refitting \emph{every} subset model on responses
#'   simulated from the full model. \strong{lme4 engine and exact attribution
#'   only}; the cost is \code{n_boot} times the number of subset models (up to
#'   \code{n_boot * 2^k} refits). Default FALSE.
#' @param n_boot Number of bootstrap draws if \code{bootstrap = TRUE}. Default 1000.
#' @param conf_level Confidence level for bootstrap intervals. Default 0.95.
#' @param estimation Variance-estimation basis for the between-stratum variances the
#'   attribution differences across subset models, \code{"fitted"} (default) or
#'   \code{"ML"}; see \code{\link{calculate_pcv}}. Affects Gaussian \code{lmer} fits
#'   only.
#'
#' @return An object of class \code{maihda_pcv_importance}: a list with
#'   \item{importance}{Data frame with one row per variable (in the order of
#'     \code{vars}): \code{Contribution} (the variable's share of the total
#'     PCV, on the PCV scale), \code{Share} (\code{Contribution / total_pcv};
#'     \code{NA} when the total PCV is zero), plus \code{MC_SE} (Monte-Carlo
#'     standard errors, \code{approx = "montecarlo"} only) and
#'     \code{CI_lower}/\code{CI_upper} (\code{bootstrap = TRUE} only).}
#'   \item{total_pcv}{The full-model PCV, \eqn{(V_0 - V_{full}) / V_0}; the
#'     contributions sum to this value (efficiency).}
#'   \item{null_variance, full_variance}{Between-stratum variance of the null
#'     and the full (all-\code{vars}) model.}
#'   \item{subsets}{Data frame of every subset model fit: the variables in the
#'     subset, its size, between-stratum variance, and PCV \eqn{v(S)}.}
#'   \item{conditional, complete_dominance}{\code{method = "dominance"} only:
#'     the conditional dominance matrix (variables x adjustment-set size) and
#'     the pairwise complete-dominance matrix.}
#'   \item{method, approx, n_perm, n_fits, n_obs, engine, family, context,
#'     bootstrap, conf_level, n_boot_ok, estimation, estimation_used}{Metadata;
#'     \code{n_fits} counts the distinct models fit (including the null),
#'     \code{estimation} is the variance-estimation basis
#'     (\code{"fitted"}/\code{"ML"}) requested for every subset model's
#'     between-stratum variance, and \code{estimation_used} is the basis actually used
#'     (\code{"mixed"} when an ML refit was skipped at the boundary, leaving a subset
#'     model on REML).}
#'
#' @details
#' \strong{Value function and efficiency.} Write \eqn{V_0} for the null model's
#' between-stratum variance and \eqn{V(S)} for the between-stratum variance
#' after adding the variable subset \eqn{S} as fixed effects (each put on the
#' \code{estimation} basis -- the default \code{"fitted"} keeps each \code{lmer}
#' fit's own REML variance, \code{"ML"} refits it with maximum likelihood --
#' exactly as in \code{\link{calculate_pcv}} and \code{\link{stepwise_pcv}}, so
#' all attributions live on the same scale as the rest of the package). The value
#' function is \eqn{v(S) = (V_0 - V(S)) / V_0} -- the total PCV of the model
#' that adds \eqn{S} -- and the Shapley contribution of variable \eqn{i}
#' averages its marginal PCV \eqn{v(S \cup \{i\}) - v(S)} over all subsets with
#' the usual Shapley weights. Because every permutation's marginals telescope
#' to \eqn{v(N)}, the contributions of \emph{every} method here (including the
#' Monte-Carlo approximation, for any \code{n_perm}) sum exactly to the
#' full-model total PCV. This is the multilevel-PCV analogue of the LMG /
#' Shorrocks-Shapley decomposition of \eqn{R^2} (Groemping 2006; Shorrocks 2013).
#'
#' \strong{Sequential method (deprecated) vs. \code{stepwise_pcv()}.} The
#' \code{"sequential"} method is \strong{deprecated} in favour of
#' \code{\link{stepwise_pcv}}, which owns the sequential trajectory and also
#' reports the \code{Step_PCV} column and the discriminatory-accuracy path.
#' While it remains, its contributions are the \emph{increments in total PCV}
#' along the entry order, \eqn{v(\{x_1..x_i\}) - v(\{x_1..x_{i-1}\})} -- i.e.
#' \code{diff(c(0, Total_PCV))} of the \code{\link{stepwise_pcv}} table, so
#' they sum to the same total. They are \emph{not} the \code{Step_PCV} column,
#' which normalises each step by the \emph{previous} step's variance rather
#' than by \eqn{V_0}.
#'
#' \strong{Exact vs. Monte-Carlo cost.} The exact attribution fits
#' \eqn{2^k - 1} non-empty subset models plus the null -- \eqn{2^k} fits in
#' total, with every fit cached and reused across all marginal differences:
#' 256 fits at \eqn{k = 8}, 1024 at
#' \eqn{k = 10}. That is feasible for \code{lme4} but quickly infeasible
#' beyond, and \emph{each} of those fits is a separate Stan run under
#' \code{engine = "brms"} -- exact attribution on brms is therefore strongly
#' discouraged except for very small \eqn{k}. The Monte-Carlo route samples
#' \code{n_perm} entry orders (unbiased for the Shapley values, since a
#' uniform random permutation reproduces the Shapley weights), reports a
#' per-variable Monte-Carlo standard error, and warns when the largest
#' \code{MC_SE} exceeds 0.01 on the PCV scale (increase \code{n_perm}).
#'
#' \strong{Latent-scale families and rescaling.} For binomial/ordinal (and, in
#' attenuated form, count) families, adding a predictor that varies
#' \emph{within} strata rescales the latent metric, so its marginal PCV mixes
#' explained variance with rescaling (see the latent-scale note in
#' \code{\link{calculate_pcv}}); Shapley values of the PCV inherit this.
#' Attribution among the stratum \emph{dimensions} -- constant within each
#' stratum -- is largely unaffected.
#'
#' \strong{Suppression and negative contributions.} The PCV can be negative,
#' so a contribution (and its \code{Share}) can be negative or exceed 100\%.
#' Efficiency still holds; contributions are reported \emph{signed} and no
#' non-negative normalisation is applied. A negative Shapley contribution
#' flags a suppressor-style variable whose inclusion tends to \emph{raise} the
#' between-stratum variance.
#'
#' \strong{Bootstrap.} With \code{bootstrap = TRUE} the whole attribution is
#' bootstrapped: responses are simulated from the full model and every subset
#' model is refit per draw (\code{lme4::refit}), giving percentile intervals
#' per contribution -- \code{n_boot} times \code{n_fits} refits, so gate it by
#' cost. Available for \code{engine = "lme4"} with exact attribution (any
#' \code{method}); the Monte-Carlo approximation already carries its own
#' sampling error and is not combined with the bootstrap. As in
#' \code{\link{calculate_pcv}}, \code{refit()} holds a \code{glmer.nb}
#' dispersion parameter fixed at its original estimate.
#'
#' @section Reproducibility:
#' \code{approx = "montecarlo"} draws random permutations with the session RNG;
#' call \code{set.seed()} first for reproducible results. The exact methods are
#' deterministic.
#'
#' @examples
#' \donttest{
#' strata <- make_strata(maihda_sim_data, c("gender", "race"))
#'
#' # Order-invariant split of the PCV across the two dimensions and age
#' imp <- pcv_importance(strata$data, "health_outcome", c("gender", "race", "age"))
#' print(imp)
#' plot(imp)
#'
#' # Dominance analysis: adds conditional / complete dominance detail
#' pcv_importance(strata$data, "health_outcome", c("gender", "race", "age"),
#'                method = "dominance")
#'
#' # Monte-Carlo approximation (set a seed for reproducibility)
#' set.seed(42)
#' pcv_importance(strata$data, "health_outcome", c("gender", "race", "age"),
#'                approx = "montecarlo", n_perm = 500)
#' }
#'
#' @seealso \code{\link{stepwise_pcv}} for the sequential path with the
#'   discriminatory-accuracy trajectory, \code{\link{calculate_pcv}} for the
#'   two-model PCV and the latent-scale caveat, and \code{\link{maihda}} for
#'   the canonical null-vs-adjusted decomposition.
#'
#' @references
#' Budescu, D. V. (1993). Dominance analysis: a new approach to the problem of
#' relative importance of predictors in multiple regression.
#' \emph{Psychological Bulletin}, 114(3), 542-551.
#'
#' Groemping, U. (2006). Relative importance for linear regression in R: the
#' package relaimpo. \emph{Journal of Statistical Software}, 17(1), 1-27.
#'
#' Shorrocks, A. F. (2013). Decomposition procedures for distributional
#' analysis: a unified framework based on the Shapley value. \emph{Journal of
#' Economic Inequality}, 11, 99-126.
#'
#' @export
pcv_importance <- function(data, outcome, vars,
                           method = c("shapley", "sequential", "dominance"),
                           approx = c("exact", "montecarlo"),
                           n_perm = 2000,
                           engine = "lme4", family = "gaussian",
                           context = NULL, sampling_weights = NULL,
                           bootstrap = FALSE, n_boot = 1000, conf_level = 0.95,
                           estimation = c("fitted", "ML")) {
  method <- match.arg(method)
  if (identical(method, "sequential")) {
    warning("pcv_importance(method = \"sequential\") is deprecated and will be ",
            "removed in a future release. Use stepwise_pcv() for the sequential ",
            "(order-dependent) path -- it also reports the step-specific PCV ",
            "and, for a binary outcome, the discriminatory-accuracy trajectory ",
            "(AUC, MOR); or use method = \"shapley\" for an order-invariant ",
            "split of the same total PCV.", call. = FALSE)
  }
  approx_missing <- missing(approx)
  approx <- match.arg(approx)
  estimation <- match.arg(estimation)

  if (!is.character(vars) || length(vars) < 1 || anyNA(vars)) {
    stop("'vars' must be a character vector naming at least one predictor.",
         call. = FALSE)
  }
  if (anyDuplicated(vars)) {
    stop("'vars' contains duplicated variable names: ",
         paste(unique(vars[duplicated(vars)]), collapse = ", "),
         ". Each predictor can enter the attribution once.", call. = FALSE)
  }
  k <- length(vars)
  # Subset bitmasks are 32-bit integers, and even the Monte-Carlo route fits
  # k models per NEW permutation -- far past any realistic MAIHDA predictor set.
  if (k > 25) {
    stop("pcv_importance() supports at most 25 variables (", k, " supplied).",
         call. = FALSE)
  }

  if (!is.logical(bootstrap) || length(bootstrap) != 1 || is.na(bootstrap)) {
    stop("'bootstrap' must be TRUE or FALSE.", call. = FALSE)
  }
  if (bootstrap) {
    bootstrap_args <- maihda_validate_bootstrap_args(n_boot, conf_level)
    n_boot <- bootstrap_args$n_boot
    conf_level <- bootstrap_args$conf_level
  }

  # Shared preamble with stepwise_pcv(): validation, engine/family handshakes,
  # ONE complete-case analytic sample, auto-binned dimension reconstruction.
  setup <- maihda_pcv_attribution_setup(
    data, outcome, vars, engine = engine, family = family, context = context,
    sampling_weights = sampling_weights,
    engine_missing = missing(engine), family_missing = missing(family),
    fn = "pcv_importance", what = "the PCV attribution"
  )
  data <- setup$data
  model_terms <- setup$model_terms
  engine <- setup$engine
  family <- setup$family
  context <- setup$context
  sampling_weights <- setup$sampling_weights

  # Resolve the attribution route. 'approx' only ever means anything for the
  # Shapley method: sequential fits just the k path models, and dominance needs
  # every subset (its conditional/complete tables average and compare over ALL
  # adjustment sets), so a Monte-Carlo dominance request is rejected rather
  # than silently degraded.
  if (method == "sequential") {
    approx <- "exact"
  } else if (method == "dominance") {
    if (!approx_missing && approx == "montecarlo") {
      stop("pcv_importance(): dominance analysis needs every subset model ",
           "(approx = \"exact\"); for a large 'vars' use method = \"shapley\" ",
           "with approx = \"montecarlo\" instead.", call. = FALSE)
    }
    approx <- "exact"
  } else if (method == "shapley" && approx_missing && k > 10) {
    approx <- "montecarlo"
    message("pcv_importance(): ", k, " variables would need ",
            bitwShiftL(1L, k) - 1L, " subset fits for the exact Shapley ",
            "attribution; using approx = \"montecarlo\" with n_perm = ",
            n_perm, ". Pass approx = \"exact\" to force the exact computation.")
  }
  if (approx == "exact" && method != "sequential" && k > 15) {
    stop("pcv_importance(): the exact attribution over ", k, " variables needs ",
         "2^", k, " - 1 subset model fits, which is not feasible. Use ",
         "method = \"shapley\" with approx = \"montecarlo\".", call. = FALSE)
  }
  if (approx == "montecarlo") {
    if (!is.numeric(n_perm) || length(n_perm) != 1 || is.na(n_perm) ||
        !is.finite(n_perm) || n_perm < 10 || n_perm > 1e6 ||
        n_perm != floor(n_perm)) {
      stop("'n_perm' must be a single whole number between 10 and 1e6.",
           call. = FALSE)
    }
    n_perm <- as.integer(n_perm)
  }

  if (bootstrap && !identical(engine, "lme4")) {
    stop("Bootstrap intervals for the PCV attribution are only available for ",
         "lme4 models (the parametric bootstrap relies on lme4's ",
         "simulate()/refit()). Call pcv_importance() with bootstrap = FALSE ",
         "for a point attribution.", call. = FALSE)
  }
  if (bootstrap && approx == "montecarlo") {
    stop("pcv_importance(): the bootstrap is available for the exact ",
         "attribution only; approx = \"montecarlo\" already carries its own ",
         "Monte-Carlo error (MC_SE). Use approx = \"exact\" (feasible up to ",
         "~10 variables) or bootstrap = FALSE.", call. = FALSE)
  }

  # Message-only fit-count bound; computed in doubles so a large n_perm * k
  # cannot overflow integer arithmetic.
  n_subsets_planned <- if (approx == "exact" && method != "sequential") {
    2^k - 1
  } else if (method == "sequential") {
    k
  } else {
    min(2^k - 1, as.numeric(n_perm) * max(k - 1, 1) + k)
  }
  if (identical(engine, "brms")) {
    message("pcv_importance(): engine = \"brms\" runs a separate Stan fit for ",
            "the null model and each distinct predictor subset (up to ",
            n_subsets_planned, " subsets here); this can be very slow. ",
            "Exact Shapley on brms is discouraged -- see ?pcv_importance.")
  } else if (approx == "exact" && method != "sequential" && k >= 8) {
    message("pcv_importance(): fitting all ", n_subsets_planned,
            " non-empty subsets of the ", k, " variables (plus the null model).")
  }

  # ---- null model / value-function machinery ---------------------------------
  # Each subset model is fit once, put on the requested variance-estimation basis
  # (estimation = "ML" refits a REML lmer fit with ML; "fitted", the default, keeps
  # it -- the variances are compared across models with different fixed effects, see
  # calculate_pcv()), and cached by bitmask. Models themselves are retained
  # only when the bootstrap needs to refit them.
  keep_models <- bootstrap
  cache <- new.env(parent = emptyenv())
  # Set when any subset model keeps its REML fit because its ML refit was skipped at the
  # boundary OR failed (maihda_pcv_refit_ml()); under estimation = "ML" either makes the
  # attribution's basis "mixed", not a pure ML comparison.
  ml_refit_skipped_any <- FALSE
  ml_refit_failed_any <- FALSE

  fit_mask <- function(mask) {
    key <- as.character(mask)
    entry <- cache[[key]]
    if (!is.null(entry)) return(entry)
    terms_s <- model_terms[maihda_mask_bits(mask, k)]
    fmla <- maihda_formula_with_stratum(outcome, terms_s)
    mod <- tryCatch(
      maihda_pcv_apply_estimation(
        fit_maihda(fmla, data, engine = engine,
                   family = family, context = context,
                   sampling_weights = sampling_weights), estimation),
      error = function(e) {
        stop("pcv_importance(): the subset model {",
             paste(vars[maihda_mask_bits(mask, k)], collapse = ", "),
             "} failed to fit: ", conditionMessage(e), call. = FALSE)
      })
    if (isTRUE(mod$ml_refit_skipped_boundary)) ml_refit_skipped_any <<- TRUE
    if (isTRUE(mod$ml_refit_failed)) ml_refit_failed_any <<- TRUE
    variance <- extract_between_variance(mod)
    if (!is.finite(variance)) {
      stop("pcv_importance(): no finite between-stratum variance for the ",
           "subset model {",
           paste(vars[maihda_mask_bits(mask, k)], collapse = ", "), "}.",
           call. = FALSE)
    }
    entry <- list(variance = variance,
                  # Flag a boundary-level (effectively singular) fit. The null (mask
                  # 0) is the PCV denominator, so its flag guards the degenerate-null
                  # stop below; the FULL model's flag is read after the attribution to
                  # warn when the shares saturate near 100% as a singular-fit artefact
                  # rather than genuine attenuation (parallel to calculate_pcv()'s
                  # adjusted-model boundary flag).
                  at_boundary = maihda_pcv_null_at_boundary(mod),
                  model = if (keep_models) mod$model else NULL,
                  family_key = maihda_model_family_key(mod),
                  family_name = maihda_model_family_name(mod))
    assign(key, entry, envir = cache)
    entry
  }

  null_entry <- fit_mask(0L)
  null_variance <- null_entry$variance
  # lme4: reject a degenerate null denominator first. v_of() divides every
  # contribution by null_variance, so a non-positive OR a strictly-positive but
  # boundary-level null variance (an effectively singular fit, e.g. 1.56e-17 --
  # which passes a plain > 0 test) would explode the whole attribution. The lme4
  # singularity guard calculate_pcv()/stepwise_pcv() apply covers both.
  if (isTRUE(null_entry$at_boundary)) {
    maihda_pcv_degenerate_null_stop("the null model")
  }
  # The other engines (brms/wemix/ordinal) are not flagged above; a non-positive
  # null variance there still has no PCV to attribute.
  if (null_variance <= 0) {
    stop("Between-stratum variance in the null model is zero or negative. ",
         "The PCV attribution cannot be calculated. This may indicate a ",
         "singular fit or no between-stratum variation.", call. = FALSE)
  }

  v_of <- function(mask) {
    if (mask == 0L) return(0)
    (null_variance - fit_mask(mask)$variance) / null_variance
  }

  full_mask <- bitwShiftL(1L, k) - 1L

  # ---- point estimates --------------------------------------------------------
  # Pure computation given a v-lookup, so the bootstrap below can re-run the
  # identical combinatorics on refit variances. Returns the contribution vector
  # (plus the dominance tables at the point estimate only).
  mc_perms <- NULL
  if (approx == "montecarlo") {
    mc_perms <- t(vapply(seq_len(n_perm), function(p) sample.int(k),
                         integer(k)))
  }
  # Prefix masks for the sequential path: the first i variables occupy bits
  # 0..i-1, so the i-th prefix mask is 2^i - 1.
  seq_masks <- if (method == "sequential") bitwShiftL(1L, seq_len(k)) - 1L else NULL

  compute_phi <- function(vfun) {
    if (method == "sequential") {
      v_path <- vapply(seq_masks, vfun, numeric(1))
      return(diff(c(0, v_path)))
    }
    if (approx == "montecarlo") {
      marg <- matrix(NA_real_, nrow = n_perm, ncol = k)
      for (p in seq_len(n_perm)) {
        mask <- 0L
        v_prev <- 0
        for (pos in seq_len(k)) {
          i <- mc_perms[p, pos]
          mask <- bitwOr(mask, bitwShiftL(1L, i - 1L))
          v_cur <- vfun(mask)
          marg[p, i] <- v_cur - v_prev
          v_prev <- v_cur
        }
      }
      return(structure(colMeans(marg),
                       mc_se = apply(marg, 2, stats::sd) / sqrt(n_perm)))
    }
    # exact shapley / dominance: needs the complete v-table
    v_table <- c(0, vapply(seq_len(bitwShiftL(1L, k) - 1L), vfun, numeric(1)))
    maihda_shapley_exact_phi(v_table, k)
  }

  phi <- compute_phi(v_of)
  mc_se <- attr(phi, "mc_se")
  phi <- as.numeric(phi)

  # Flag when the FULL model's between-stratum variance is at the singularity boundary
  # (any estimation basis): the attributed shares then sum to ~100% as a singular-fit
  # artefact, not genuine attenuation -- the pcv_importance() analogue of
  # calculate_pcv()'s adjusted-model boundary flag. This is the COMMON case for additive
  # strata (a full additive model leaves ~0 interaction variance), so it is carried as a
  # silent status attribute (full_at_boundary, attached to the result below) for
  # programmatic inspection and the print method, NOT raised as a per-call warning.
  # compute_phi() has fit the full model on every path (its mask is the last prefix /
  # permutation step / subset).
  full_entry <- cache[[as.character(full_mask)]]
  full_at_boundary <- !is.null(full_entry) && isTRUE(full_entry$at_boundary)

  dominance <- NULL
  if (method == "dominance") {
    v_table <- c(0, vapply(seq_len(full_mask), v_of, numeric(1)))
    dominance <- maihda_dominance_tables(v_table, k, vars)
    # General dominance IS the Shapley value; report it as the contribution
    # (identical to phi up to floating error, pinned by tests).
    phi <- as.numeric(dominance$general)
  }

  total_pcv <- v_of(full_mask)
  full_variance <- fit_mask(full_mask)$variance

  if (approx == "montecarlo" && any(is.finite(mc_se)) &&
      max(mc_se, na.rm = TRUE) > 0.01) {
    warning("pcv_importance(): the largest Monte-Carlo standard error is ",
            sprintf("%.4f", max(mc_se, na.rm = TRUE)),
            " on the PCV scale; the Shapley estimates may not have converged. ",
            "Increase n_perm.", call. = FALSE)
  }

  importance <- data.frame(
    Variable = vars,
    Contribution = phi,
    Share = if (abs(total_pcv) > 1e-12) phi / total_pcv else NA_real_,
    stringsAsFactors = FALSE
  )
  if (!is.null(mc_se)) {
    importance$MC_SE <- mc_se
  }

  fitted_masks <- sort(as.integer(ls(cache)))
  subsets <- data.frame(
    Variables = vapply(fitted_masks, function(m) {
      if (m == 0L) "(none: null model)" else
        paste(vars[maihda_mask_bits(m, k)], collapse = " + ")
    }, character(1)),
    Size = vapply(fitted_masks, function(m) length(maihda_mask_bits(m, k)),
                  integer(1)),
    Variance = vapply(fitted_masks, function(m) cache[[as.character(m)]]$variance,
                      numeric(1)),
    PCV = vapply(fitted_masks, v_of, numeric(1)),
    stringsAsFactors = FALSE
  )
  subsets <- subsets[order(subsets$Size, fitted_masks), , drop = FALSE]
  rownames(subsets) <- NULL

  result <- list(
    importance = importance,
    method = method,
    approx = if (method == "shapley") approx else "exact",
    n_perm = if (approx == "montecarlo") n_perm else NULL,
    total_pcv = total_pcv,
    null_variance = null_variance,
    full_variance = full_variance,
    full_at_boundary = full_at_boundary,
    subsets = subsets,
    n_fits = length(fitted_masks),
    n_obs = nrow(data),
    outcome = outcome,
    vars = vars,
    engine = engine,
    family = null_entry$family_key,
    family_name = null_entry$family_name,
    estimation = estimation,
    estimation_used = maihda_pcv_estimation_used(
      estimation, ml_refit_skipped_any || ml_refit_failed_any, engine = engine),
    context = context,
    sampling_weights = sampling_weights,
    conditional = dominance$conditional,
    complete_dominance = dominance$complete,
    bootstrap = FALSE
  )

  # ---- parametric bootstrap ---------------------------------------------------
  # Simulate responses from the FULL model (all vars -- the analogue of
  # bootstrap_pcv() simulating from model2), refit every cached subset model
  # per draw, and re-run the identical attribution on the refit variances. A
  # draw fails as a whole (NA row) if any refit errors or its null variance is
  # not positive.
  if (bootstrap) {
    boot_masks <- fitted_masks
    boot_keys <- as.character(boot_masks)
    boot_models <- lapply(boot_keys, function(key) cache[[key]]$model)
    names(boot_models) <- boot_keys
    message("pcv_importance(): parametric bootstrap -- ", n_boot, " draws x ",
            length(boot_masks), " subset refits (",
            n_boot * length(boot_masks), " refits total).")
    sim_data <- stats::simulate(boot_models[[as.character(full_mask)]],
                                nsim = n_boot)
    boot_phi <- matrix(NA_real_, nrow = n_boot, ncol = k)
    for (b in seq_len(n_boot)) {
      # Refit every subset on the simulated response, keeping the refit models so
      # the null draw's boundary can be tested before it enters a denominator.
      refits_b <- tryCatch(
        lapply(boot_keys, function(key)
          lme4::refit(boot_models[[key]], newresp = sim_data[[b]])),
        error = function(e) NULL)
      if (is.null(refits_b)) next
      names(refits_b) <- boot_keys
      # Drop the whole draw when the null refit lands on the zero boundary: a
      # non-positive OR a strictly-positive-but-degenerate (~1e-12) null variance
      # would blow up every ratio below, so exclude it exactly as the
      # point-estimate guard and calculate_pcv()'s bootstrap do.
      if (isTRUE(tryCatch(maihda_stratum_at_boundary_lme4(refits_b[["0"]]),
                          error = function(e) TRUE))) next
      variances_b <- vapply(refits_b, maihda_stratum_variance_lme4, numeric(1))
      v0_b <- variances_b[["0"]]
      if (!is.finite(v0_b) || v0_b <= 0) next
      vfun_b <- function(mask) {
        if (mask == 0L) return(0)
        (v0_b - variances_b[[as.character(mask)]]) / v0_b
      }
      # For method = "dominance" this is still the right draw statistic:
      # general dominance is mathematically the Shapley value.
      boot_phi[b, ] <- as.numeric(compute_phi(vfun_b))
    }

    # Draws fail as a whole, so the successful-draw count is shared across
    # variables; let the first column carry the low-success warning and
    # silence the identical repeats.
    ci_list <- vector("list", k)
    for (i in seq_len(k)) {
      ci_list[[i]] <- if (i == 1) {
        maihda_bootstrap_ci(boot_phi[, i], n_boot, conf_level,
                            "PCV attribution")
      } else {
        suppressWarnings(maihda_bootstrap_ci(boot_phi[, i], n_boot, conf_level,
                                             "PCV attribution"))
      }
    }
    result$importance$CI_lower <- vapply(ci_list, function(ci) ci[1], numeric(1))
    result$importance$CI_upper <- vapply(ci_list, function(ci) ci[2], numeric(1))
    result$bootstrap <- TRUE
    result$conf_level <- conf_level
    result$n_boot <- n_boot
    result$n_boot_ok <- attr(ci_list[[1]], "n_ok")
  }

  class(result) <- "maihda_pcv_importance"
  result
}

#' Print a PCV importance attribution
#'
#' @param x A \code{maihda_pcv_importance} object from \code{\link{pcv_importance}}.
#' @param digits Number of digits for the contribution columns. Default 4.
#' @param ... Additional arguments (not used).
#' @return Invisibly, \code{x}.
#' @export
print.maihda_pcv_importance <- function(x, digits = 4, ...) {
  pal <- maihda_palette()
  method_label <- switch(x$method,
    shapley = if (identical(x$approx, "montecarlo")) {
      sprintf("Shapley values (Monte-Carlo, %d permutations)", x$n_perm)
    } else {
      "Shapley values (exact)"
    },
    sequential = "Sequential increments (order-dependent)",
    dominance = "Dominance analysis (general dominance = Shapley)"
  )

  cat(pal$bold("PCV Attribution Across Predictors"), "\n", sep = "")
  cat("=================================\n\n")
  cat(sprintf("Method:  %s\n", pal$accent(method_label)))
  cat(sprintf("Outcome: %s   Engine: %s (%s)\n", x$outcome, x$engine, x$family))
  cat(sprintf("Analytic sample: %d observations; %d models fit (incl. null).\n",
              x$n_obs, x$n_fits))
  basis <- x$estimation_used
  if (is.null(basis)) basis <- x$estimation
  if (!is.null(basis)) {
    cat(sprintf("Variance basis: %s\n", maihda_pcv_basis_label(basis)))
  }
  if (!is.null(x$context)) {
    cat(sprintf("Context: %s (contributions are net of the context random intercept).\n",
                paste(x$context, collapse = ", ")))
  }
  cat(sprintf("\nBetween-stratum variance: null %.6f -> full model %.6f\n",
              x$null_variance, x$full_variance))
  cat(sprintf("Total PCV (null -> all variables): %s\n\n",
              pal$accent(sprintf(paste0("%.", digits, "f"), x$total_pcv))))
  if (isTRUE(x$full_at_boundary)) {
    cat(pal$muted(paste0(
      "Note: the full model's between-stratum variance is at the singularity boundary (a\n",
      "singular fit), so Total PCV is pinned near 100%. This is consistent with genuinely\n",
      "additive strata (no interaction beyond the main effects) as well as a degenerate\n",
      "fit; the two cannot be told apart.\n\n")))
  }

  tab <- x$importance
  fmt <- paste0("%.", digits, "f")
  disp <- data.frame(
    Variable = tab$Variable,
    Contribution = sprintf(fmt, tab$Contribution),
    Share = ifelse(is.na(tab$Share), "--",
                   sprintf("%+.1f%%", 100 * tab$Share)),
    stringsAsFactors = FALSE
  )
  if ("MC_SE" %in% names(tab)) {
    disp$MC_SE <- sprintf(fmt, tab$MC_SE)
  }
  if (all(c("CI_lower", "CI_upper") %in% names(tab))) {
    conf_pct <- if (!is.null(x$conf_level)) x$conf_level * 100 else 95
    disp[[sprintf("CI_%.0f%%", conf_pct)]] <-
      sprintf(paste0("[", fmt, ", ", fmt, "]"), tab$CI_lower, tab$CI_upper)
  }
  total_row <- disp[1, , drop = FALSE]
  total_row[1, ] <- ""
  total_row$Variable <- "Total"
  total_row$Contribution <- sprintf(fmt, sum(tab$Contribution))
  total_row$Share <- if (anyNA(tab$Share)) "--" else "100.0%"
  print(rbind(disp, total_row), row.names = FALSE, right = TRUE)

  cat(pal$muted(paste0(
    "\nContributions are shares of the null model's between-stratum variance\n",
    "explained (the PCV scale) and sum to the full-model Total PCV (efficiency).\n")))
  if (identical(x$method, "sequential")) {
    cat(pal$muted(paste0(
      "Order-dependent: each value is the increment in Total PCV when the\n",
      "variable enters IN THE ORDER GIVEN (= diff of stepwise_pcv()'s Total_PCV,\n",
      "not its Step_PCV). Use method = \"shapley\" for an order-invariant split.\n")))
  }
  if (x$bootstrap && !is.null(x$n_boot_ok)) {
    cat(pal$muted(sprintf(
      "Parametric bootstrap intervals from %d successful draws of %d.\n",
      as.integer(x$n_boot_ok), as.integer(x$n_boot))))
  }
  if (any(x$importance$Contribution < 0)) {
    cat(pal$muted(paste0(
      "Negative contributions flag suppression: the variable tends to RAISE the\n",
      "between-stratum variance; signed values are reported (no normalisation).\n")))
  }
  if (!is.na(x$family_name) &&
      !identical(x$family_name, "gaussian")) {
    cat(pal$muted(paste0(
      "Latent-scale caveat: for this family, contributions of predictors that\n",
      "vary WITHIN strata mix explained variance with latent-scale rescaling\n",
      "(see ?calculate_pcv); stratum-constant dimensions are largely unaffected.\n")))
  }

  if (identical(x$method, "dominance") && !is.null(x$conditional)) {
    cat("\n", pal$bold("Conditional dominance"), sep = "")
    cat(pal$muted(" (average marginal PCV by adjustment-set size):\n"))
    print(round(x$conditional, digits))
    comp <- x$complete_dominance
    pairs <- character(0)
    if (!is.null(comp)) {
      for (i in seq_len(nrow(comp))) {
        for (j in seq_len(ncol(comp))) {
          if (i != j && isTRUE(comp[i, j])) {
            pairs <- c(pairs, paste(rownames(comp)[i], ">", colnames(comp)[j]))
          }
        }
      }
    }
    cat(pal$bold("Complete dominance:"), " ",
        if (length(pairs) > 0) paste(pairs, collapse = ", ") else
          "none established", "\n", sep = "")
  }

  invisible(x)
}

#' Plot a PCV importance attribution
#'
#' Bar chart of the per-variable contributions to the total PCV, ordered by
#' size, with a dashed zero line (negative bars flag suppression) and, when
#' available, parametric-bootstrap confidence intervals. The bars sum to the
#' full-model Total PCV (the efficiency property).
#'
#' @param x A \code{maihda_pcv_importance} object from \code{\link{pcv_importance}}.
#' @param ... Additional arguments (not used).
#' @return A ggplot object.
#' @export
#' @import ggplot2
#' @importFrom rlang .data
plot.maihda_pcv_importance <- function(x, ...) {
  if (!inherits(x, "maihda_pcv_importance")) {
    stop("'x' must be a maihda_pcv_importance object from pcv_importance().",
         call. = FALSE)
  }
  df <- x$importance
  # Ascending order + coord_flip puts the largest contribution on top. The
  # sequential method keeps its (meaningful) entry order, top-down.
  if (identical(x$method, "sequential")) {
    df$Variable <- factor(df$Variable, levels = rev(df$Variable))
  } else {
    df <- df[order(df$Contribution), , drop = FALSE]
    df$Variable <- factor(df$Variable, levels = df$Variable)
  }
  has_ci <- all(c("CI_lower", "CI_upper") %in% names(df))

  subtitle <- switch(x$method,
    shapley = if (identical(x$approx, "montecarlo")) {
      sprintf("Shapley values (Monte-Carlo, %d permutations)", x$n_perm)
    } else {
      "Shapley values (exact, order-invariant)"
    },
    sequential = "Sequential increments in the order given (order-dependent)",
    dominance = "General dominance (= Shapley values)"
  )
  caption <- sprintf("Contributions sum to the full-model Total PCV = %.3f.",
                     x$total_pcv)

  p <- ggplot(df, aes(x = .data$Variable, y = .data$Contribution)) +
    geom_col(fill = "#0072B2") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    coord_flip() +
    labs(
      title = "PCV attribution across predictors",
      subtitle = subtitle,
      x = NULL,
      y = "Contribution to total PCV (share of null between-stratum variance)",
      caption = caption
    ) +
    theme_maihda() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      plot.caption = element_text(hjust = 0.5, face = "italic", size = 9)
    )
  if (has_ci) {
    p <- p + geom_errorbar(aes(ymin = .data$CI_lower, ymax = .data$CI_upper),
                           width = 0.2, color = "#D55E00")
  }
  p
}
