#' Fit MAIHDA Model
#'
#' Fits a multilevel model for MAIHDA (Multilevel Analysis of Individual
#' Heterogeneity and Discriminatory Accuracy) using lme4, brms, or -- for
#' design-weighted (survey) data -- WeMix.
#'
#' @param formula A formula specifying the model. Can include a random effect
#'   for stratum (e.g., \code{outcome ~ fixed_vars + (1 | stratum)}) or can
#'   directly specify the intersection variables to be used for forming strata
#'   (e.g., \code{outcome ~ fixed_vars + (1 | var1:var2:var3)}). If variables
#'   other than "stratum" are provided in the random effect, \code{\link{make_strata}}
#'   will be called internally to compute the strata and the formula will be
#'   updated.
#' @param data A data frame containing the variables in the formula.
#' @param engine Character string specifying which engine to use: "lme4"
#'   (default), "brms", "wemix" (design-weighted pseudo-maximum-likelihood via
#'   \code{WeMix::mix()}; requires \code{sampling_weights}), or "ordinal"
#'   (cumulative link mixed model via \code{ordinal::clmm()}; requires an
#'   ordinal family). When \code{sampling_weights} is supplied and \code{engine}
#'   is left at its default, the engine switches to "wemix" automatically (with
#'   a message); likewise an ordinal family (or an auto-detected ordered-factor
#'   outcome) switches the default engine to "ordinal".
#' @param family Character string, family object, or family function specifying
#'   the model family. Common options: "gaussian", "binomial", "poisson",
#'   "negbinomial". Default is "gaussian".
#'   \code{family = "negbinomial"} fits an overdispersed count model with the
#'   dispersion parameter theta \emph{estimated} from the data: lme4 via
#'   \code{lme4::glmer.nb()} and brms via its \code{shape} parameter (log link
#'   only; not supported by the wemix engine). A fixed-theta
#'   \code{MASS::negative.binomial(theta)} family object is also accepted with
#'   \code{engine = "lme4"} and is fitted with \code{glmer()}, honouring the
#'   supplied theta.
#'   \code{family = "ordinal"} (alias \code{"cumulative"}; or
#'   \code{\link{maihda_cumulative}("probit")} / \code{brms::cumulative()} for a
#'   non-logit link) fits a cumulative (proportional-odds) model for an
#'   \emph{ordered-factor} outcome: \code{ordinal::clmm()} under the automatic
#'   "ordinal" engine, \code{brms::cumulative()} under \code{engine = "brms"}.
#'   The VPC/ICC lives on the latent scale (level-1 variance \eqn{\pi^2/3}
#'   logit / 1 probit, as for binomial models) and response-scale predictions
#'   are \emph{expected category scores} (categories scored 1..K in order). An
#'   ordered-factor outcome with 3+ levels under the default family selects
#'   this model automatically, with a warning. The logit and probit links are
#'   supported; \code{sampling_weights}, \code{context}, and lme4-style
#'   \code{weights}/\code{subset}/\code{offset} arguments are not available on
#'   the clmm path.
#'   If the outcome variable appears to be binary and the default family is used,
#'   the function will automatically switch to "binomial", recode two-level
#'   responses to 0/1 for \code{glmer()}, and issue a warning.
#'   When a two-level non-0/1 response is recoded (on either the auto-detected or
#'   an explicit \code{family = "binomial"} path), the mapping follows the usual
#'   convention -- the first level becomes 0 (reference) and the second becomes 1
#'   (the modeled event), where "first/second" means alphabetical order for a
#'   character outcome and the declared order for a factor. The chosen mapping is
#'   reported via a \code{message()} and stored on the result as
#'   \code{$response_recoding}; set the factor levels (or supply a 0/1 outcome) to
#'   control which level is the event.
#'   Although any valid family object is accepted for fitting, the MAIHDA variance
#'   summaries (\code{\link{summary.maihda_model}}, VPC/ICC, PCV) are only defined
#'   for \code{gaussian("identity")}, the binomial/Bernoulli families with a logit
#'   or probit link, \code{poisson("log")}, and the negative binomial with a log
#'   link (level-1 variance \code{log(1 + 1/mu + 1/theta)}; Nakagawa, Johnson &
#'   Schielzeth 2017). Other families (for example \code{Gamma(link = "log")})
#'   will fit, but \code{summary()} and the VPC/PCV helpers will stop with an
#'   "not implemented" error because no level-1 variance is defined for them.
#' @param autobin Logical indicating whether numeric variables used only for
#'   automatic strata creation should be binned by \code{\link{make_strata}}.
#'   Default is TRUE.
#' @param context Optional character vector naming one or more higher-level
#'   \emph{context} columns in \code{data} (e.g. \code{"school"},
#'   \code{"hospital"}, \code{"region"}). Each enters the model as a crossed
#'   intercept-only random effect alongside the intersectional stratum effect --
#'   \code{outcome ~ covars + (1 | stratum) + (1 | context)} -- giving the
#'   \emph{contextual cross-classified MAIHDA} of the literature (individuals
#'   cross-classified by stratum and place/institution). \code{\link{summary.maihda_model}}
#'   then partitions the unexplained variance into between-stratum vs.
#'   between-context vs. residual, and the headline VPC/ICC remains the
#'   between-stratum share (now net of the context). A context variable may not be
#'   a stratum dimension or \code{"stratum"} itself, and may not already appear as
#'   a fixed-effect term (its variance would then be absorbed by the fixed part).
#'   A context with few levels (say < 10) weakly identifies its variance and often
#'   yields a singular lme4 fit; the \code{brms} engine handles this better.
#'   Writing the random effect directly in the formula (\code{... + (1 | school)})
#'   fits the same model but is summarised generically as "Other random effects";
#'   only \code{context =} activates the labelled contextual partition.
#'   Not supported by the \code{wemix} engine.
#' @param sampling_weights Optional single character string naming a numeric
#'   column of \code{data} holding individual \emph{sampling} (survey/design)
#'   weights, for a \strong{design-weighted MAIHDA} on complex-survey data
#'   (e.g. NHANES, PISA). Sampling weights are not the same thing as lme4's
#'   \code{weights=} (precision weights, which rescale the residual variance), so
#'   combining \code{sampling_weights} with \code{engine = "lme4"} is an error.
#'   Two engines support them:
#'   \itemize{
#'     \item \code{engine = "wemix"} (chosen automatically when \code{engine} is
#'       left at its default): weighted pseudo-maximum-likelihood via
#'       \code{WeMix::mix()} (Rabe-Hesketh & Skrondal 2006), the estimator used
#'       for NAEP/PISA analysis. The individual weights enter at level 1
#'       unchanged and the level-2 (stratum) weights are 1, because
#'       intersectional strata are exhaustive population cells included with
#'       certainty. Supports \code{gaussian(identity)} and
#'       \code{binomial(logit)} models with the canonical single
#'       \code{(1 | stratum)} random intercept. Fixed-effect standard errors are
#'       design-consistent (sandwich); the VPC/PCV are reported as point
#'       estimates (no bootstrap -- see \code{\link{summary.maihda_model}}).
#'     \item \code{engine = "brms"}: the weights enter the model as likelihood
#'       weights (\code{y | weights(w)}), normalized to mean 1, giving a
#'       \emph{pseudo-posterior}: point estimates are design-consistent but
#'       credible intervals are not design-based -- interpret them cautiously.
#'   }
#'   Rows with a missing or non-positive sampling weight are dropped with a
#'   warning. Default \code{NULL} (unweighted).
#' @param ... Additional arguments passed to \code{lmer}/\code{glmer} (lme4),
#'   \code{brm} (brms), or \code{WeMix::mix()} (wemix; e.g. \code{nQuad},
#'   \code{fast}). The lme4-style \code{weights}/\code{subset}/\code{offset}
#'   arguments are not supported by the wemix engine.
#'
#' @return A maihda_model object containing:
#'   \item{model}{The fitted model object (lme4, brms, or WeMix)}
#'   \item{engine}{The engine used ("lme4", "brms", or "wemix")}
#'   \item{sampling_weights}{The sampling-weight column name when supplied,
#'     NULL otherwise}
#'   \item{formula}{The model formula}
#'   \item{data}{The data used for fitting}
#'   \item{family}{The family used}
#'   \item{strata_info}{The strata information from make_strata() if available, NULL otherwise}
#'   \item{context_vars}{The context variable name(s) when \code{context} was
#'     supplied, NULL otherwise}
#'   \item{response_recoding}{For a recoded two-level outcome, a data frame mapping
#'     each original level to its 0/1 value and role (reference/event); NULL when no
#'     recoding occurred}
#'   \item{diagnostics}{Fit-quality diagnostics (singular fit / convergence) for
#'     lme4 models, surfaced by the print and summary methods}
#'
#' @examples
#' \donttest{
#' # Standard approach: manually create strata first
#' strata_result <- make_strata(maihda_sim_data, vars = c("gender", "race", "education"))
#' model <- fit_maihda(health_outcome ~ age + (1 | stratum),
#'                     data = strata_result$data,
#'                     engine = "lme4")
#'
#' # Simplified approach: specify stratifying variables directly in the grouping structure
#' # The function internally calls make_strata() to create intersectionals
#' model2 <- fit_maihda(health_outcome ~ age + (1 | gender:race:education),
#'                      data = maihda_sim_data,
#'                      engine = "lme4")
#'
#' # Contextual cross-classified MAIHDA: strata crossed with a higher-level
#' # context (here country) -- the literature's cross-classified MAIHDA.
#' data(maihda_country_data)
#' model3 <- fit_maihda(math ~ 1 + (1 | gender:ses),
#'                      data = maihda_country_data,
#'                      context = "country")
#' summary(model3)  # between-stratum vs. between-country vs. residual
#' }
#'
#' @export
#' @importFrom lme4 lmer glmer
#' @importFrom reformulas findbars nobars
#' @importFrom rlang enquos eval_tidy
#' @importFrom stats gaussian binomial poisson
fit_maihda <- function(formula, data, engine = "lme4", family = "gaussian",
                       autobin = TRUE, context = NULL, sampling_weights = NULL,
                       ...) {
  # Input validation
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula object")
  }

  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }

  # Sampling (design) weights select the design-weighted engine. lme4 weights are
  # PRECISION weights -- feeding survey weights to lmer/glmer maximises the wrong
  # objective -- so an explicit engine = "lme4" with sampling_weights is an error
  # rather than a silent misfit; the default engine switches to "wemix".
  if (!is.null(sampling_weights)) {
    sampling_weights <- maihda_validate_sampling_weights(sampling_weights, data)
    if (missing(engine)) {
      engine <- "wemix"
      message("fit_maihda(): 'sampling_weights' supplied; using engine = \"wemix\" ",
              "(design-weighted pseudo-maximum-likelihood via WeMix). Set 'engine' ",
              "explicitly to silence this message or to choose engine = \"brms\".")
    } else if (identical(engine, "lme4")) {
      stop("Sampling weights are not supported by engine = \"lme4\": lme4's ",
           "weights are precision weights (they rescale the residual variance), ",
           "not sampling weights, so survey weights would give invalid estimates. ",
           "Use engine = \"wemix\" (pseudo-maximum-likelihood, recommended) or ",
           "engine = \"brms\" (pseudo-posterior).", call. = FALSE)
    }
  }

  if (!is.character(engine) || length(engine) != 1 ||
      !engine %in% c("lme4", "brms", "wemix", "ordinal")) {
    stop("'engine' should be one of: lme4, brms, wemix, ordinal", call. = FALSE)
  }

  context <- maihda_validate_context(context, data)

  if (identical(engine, "wemix")) {
    if (is.null(sampling_weights)) {
      stop("engine = \"wemix\" is the design-weighted MAIHDA fit and requires ",
           "'sampling_weights' (the sampling-weight column). For an unweighted ",
           "fit use engine = \"lme4\" or \"brms\".", call. = FALSE)
    }
    if (!is.null(context)) {
      stop("engine = \"wemix\" does not support 'context' (WeMix fits no crossed ",
           "random effects). Use engine = \"lme4\" or \"brms\" for a contextual ",
           "cross-classified model.", call. = FALSE)
    }
  }

  # Capture the forwarded engine arguments as quosures (each keeps its expression
  # AND its environment) and evaluate them once, here, against the data. Plain
  # `...` forwarding turns data-masked arguments into ..1/..2 promises that bypass
  # the data mask once fit_maihda is called through maihda()/compare_maihda_groups();
  # rlang::eval_tidy() instead resolves each argument against the data columns first
  # and then the caller's scope, so weights = a_column, weights = a_caller_variable,
  # and subset = y %in% c("no", "yes") all work at any nesting depth. Evaluating the
  # subset here, against the ORIGINAL response, also makes it immune to the 0/1
  # recoding below. The resulting values feed binary detection and the engine call.
  dot_vals <- lapply(rlang::enquos(...), function(q) rlang::eval_tidy(q, data = data))
  subset_value <- dot_vals[["subset"]]
  weights_value <- dot_vals[["weights"]]

  if (!is.null(sampling_weights) && "weights" %in% names(dot_vals)) {
    stop("Supply either 'sampling_weights' (design weights) or 'weights' ",
         "(precision weights), not both.", call. = FALSE)
  }
  if (identical(engine, "wemix")) {
    unsupported_dots <- intersect(c("weights", "subset", "offset"), names(dot_vals))
    if (length(unsupported_dots) > 0) {
      stop("Argument(s) not supported by engine = \"wemix\": ",
           paste(unsupported_dots, collapse = ", "),
           ". Subset or transform the data before fitting.", call. = FALSE)
    }
  }

  # Rows with a missing sampling weight are dropped by the weighted engines, so
  # base the binary-outcome detection (and any 0/1 recoding) on the same rows by
  # passing the sampling weights through the detector's missing-weight handling.
  detect_weights <- if (!is.null(sampling_weights)) {
    data[[sampling_weights]]
  } else {
    weights_value
  }

  # Automatically switch to binomial for binary outcomes if family is default.
  # Detect on the analytic sample lme4/brms will actually fit -- the model frame
  # after transformations, NA-dropping, any `subset`, and dropping rows with a
  # missing prior weight -- so an outcome that is only 0/1 once excluded rows are
  # removed is still recognised as binary. An ORDERED factor outcome with 3+
  # levels under the default family likewise auto-switches to the cumulative
  # (ordinal) model -- it would otherwise just error inside lmer() -- with the
  # binary check taking precedence (a 2-level ordered factor is a binomial
  # model). Both switches warn so the family choice is never silent.
  if (missing(family)) {
    is_binary <- tryCatch(
      maihda_response_is_binary(formula, data, subset = subset_value,
                                weights = detect_weights),
      error = function(e) FALSE)
    if (isTRUE(is_binary)) {
      warning("The outcome variable appears to be binary. Automatically switching to family = 'binomial'. To fit a Linear Probability Model, explicitly specify family = 'gaussian'.", call. = FALSE)
      family <- "binomial"
    } else if (isTRUE(tryCatch(
      maihda_response_is_ordinal(formula, data, subset = subset_value,
                                 weights = detect_weights),
      error = function(e) FALSE))) {
      warning("The outcome variable is an ordered factor. Automatically ",
              "switching to the cumulative (ordinal) model, family = 'ordinal'. ",
              "Specify a family explicitly to override.", call. = FALSE)
      family <- "ordinal"
    }
  }

  # Family <-> engine handshake for the cumulative (ordinal) model. lme4 cannot
  # fit it, so an ordinal family with the default engine auto-switches to the
  # clmm-based "ordinal" engine (mirroring the sampling_weights -> wemix switch
  # above), an explicit engine = "lme4" is an error, and engine = "ordinal"
  # without an ordinal family is an error. An explicit engine = "wemix" falls
  # through to maihda_wemix_check_family()'s targeted rejection.
  is_ordinal <- maihda_family_is_ordinal(
    if (is.function(family)) tryCatch(family(), error = function(e) NULL) else family
  )
  if (is_ordinal) {
    if (missing(engine) && is.null(sampling_weights)) {
      engine <- "ordinal"
      message("fit_maihda(): ordinal (cumulative) family; using engine = ",
              "\"ordinal\" (ordinal::clmm). Set 'engine' explicitly to silence ",
              "this message or to choose engine = \"brms\".")
    } else if (identical(engine, "lme4")) {
      stop("lme4 cannot fit a cumulative (ordinal) model. Use engine = ",
           "\"ordinal\" (ordinal::clmm, the default for this family) or ",
           "engine = \"brms\" (brms::cumulative).", call. = FALSE)
    }
  } else if (identical(engine, "ordinal")) {
    stop("engine = \"ordinal\" fits cumulative (ordinal) models; supply ",
         "family = \"ordinal\" / maihda_cumulative() (or let an ordered-factor ",
         "outcome select it automatically).", call. = FALSE)
  }
  if (identical(engine, "ordinal")) {
    if (!is.null(sampling_weights)) {
      stop("engine = \"ordinal\" does not support 'sampling_weights'. Use ",
           "engine = \"brms\" for a sampling-weighted cumulative model ",
           "(pseudo-posterior).", call. = FALSE)
    }
    if (!is.null(context)) {
      stop("engine = \"ordinal\" does not support 'context' (the clmm path fits ",
           "the canonical single (1 | stratum) structure only). Use engine = ",
           "\"brms\" for a contextual cross-classified cumulative model.",
           call. = FALSE)
    }
    unsupported_dots <- intersect(c("weights", "subset", "offset"), names(dot_vals))
    if (length(unsupported_dots) > 0) {
      stop("Argument(s) not supported by engine = \"ordinal\": ",
           paste(unsupported_dots, collapse = ", "),
           ". Subset or transform the data before fitting.", call. = FALSE)
    }
  }

  # Parse formula to find grouping variables. Automatic strata creation is only
  # safe for the documented shorthand: one intercept-only non-stratum grouping
  # term such as (1 | gender:race). More complex random-effect structures should
  # be specified explicitly after calling make_strata().
  re_terms <- reformulas::findbars(formula)
  strata_info <- attr(data, "strata_info")
  strata_vars <- attr(data, "strata_vars")
  if (is.null(strata_vars) || length(strata_vars) == 0) {
    strata_vars <- maihda_infer_strata_vars(strata_info)
  }
  strata_sep <- attr(data, "strata_sep")
  strata_autobin_info <- attr(data, "strata_autobin_info")

  if (length(re_terms) > 0) {
    grouping_vars_by_term <- lapply(re_terms, function(x) all.vars(x[[3]]))
    grouping_vars <- unique(unlist(grouping_vars_by_term, use.names = FALSE))
    has_stratum_group <- any(vapply(grouping_vars_by_term, function(vars) {
      identical(vars, "stratum")
    }, logical(1)))

    if (!has_stratum_group) {
      if (length(re_terms) != 1) {
        stop("Automatic strata creation only supports a single intercept-only random effect, ",
             "for example (1 | gender:race). For more complex random-effects structures, ",
             "call make_strata() first and include (1 | stratum) explicitly.",
             call. = FALSE)
      }

      random_lhs <- paste(deparse(re_terms[[1]][[2]]), collapse = " ")
      if (random_lhs != "1") {
        stop("Automatic strata creation only supports intercept-only random effects, ",
             "for example (1 | gender:race).",
             call. = FALSE)
      }

      if (!maihda_is_colon_interaction(re_terms[[1]][[3]])) {
        stop("Automatic strata creation supports a single variable or a colon ",
             "interaction such as (1 | gender:race). For other grouping expressions ",
             "(e.g. interaction(), paste(), cut()), call make_strata() first and use ",
             "(1 | stratum).", call. = FALSE)
      }

      strata_vars <- grouping_vars_by_term[[1]]
      missing_grouping_vars <- setdiff(strata_vars, names(data))
      if (length(missing_grouping_vars) > 0) {
        stop("Grouping variables not found in data: ",
             paste(missing_grouping_vars, collapse = ", "), call. = FALSE)
      }
      if ("stratum" %in% names(data)) {
        stop("Automatic strata creation would overwrite an existing 'stratum' column. ",
             "Use the existing (1 | stratum) model or rename/remove that column first.",
             call. = FALSE)
      }

      strata_result <- make_strata(data, vars = strata_vars, autobin = autobin)
      data$stratum <- strata_result$data$stratum
      strata_info <- strata_result$strata_info
      strata_sep <- strata_result$sep
      strata_autobin_info <- strata_result$autobin_info
      attr(data, "strata_info") <- strata_info
      attr(data, "strata_vars") <- strata_vars
      attr(data, "strata_sep") <- strata_sep
      attr(data, "strata_autobin_info") <- strata_autobin_info

      fixed_formula <- reformulas::nobars(formula)
      formula <- stats::update(fixed_formula, . ~ . + (1 | stratum))
    }
  }

  # Contextual cross-classified MAIHDA: append the higher-level context random
  # intercept(s) AFTER the stratum random effect is resolved, so the shorthand
  # (1 | var1:var2) path and the pre-built (1 | stratum) path both end up with
  # outcome ~ covars + (1 | stratum) + (1 | context). Idempotent: a context that
  # is already a random-effect grouping (e.g. when maihda() refits a derived
  # formula that carries the context term) is validated and tagged but not
  # appended again.
  context_info <- NULL
  if (!is.null(context)) {
    re_terms_now <- reformulas::findbars(formula)
    grouping_vars_now <- unique(unlist(
      lapply(re_terms_now, function(x) all.vars(x[[3]])), use.names = FALSE))
    if (!"stratum" %in% grouping_vars_now) {
      stop("'context' adds a crossed contextual random effect alongside the ",
           "intersectional stratum effect, but the formula has no stratum random ",
           "effect. Use the shorthand (1 | var1:var2) or include (1 | stratum).",
           call. = FALSE)
    }
    clash_dims <- intersect(context, strata_vars)
    if (length(clash_dims) > 0) {
      stop("Context variable(s) ", paste(clash_dims, collapse = ", "),
           " also define the intersectional strata. A variable cannot be both a ",
           "stratum dimension and a higher-level context; remove it from one of ",
           "the two roles.", call. = FALSE)
    }
    clash_fixed <- intersect(context, all.vars(reformulas::nobars(formula)[[3]]))
    if (length(clash_fixed) > 0) {
      stop("Context variable(s) ", paste(clash_fixed, collapse = ", "),
           " already appear in the fixed part of the formula, which would absorb ",
           "the context variance the contextual partition is meant to estimate. ",
           "Supply the context only via 'context', or only as a fixed effect, ",
           "not both.", call. = FALSE)
    }
    context_to_add <- setdiff(context, grouping_vars_now)
    if (length(context_to_add) > 0) {
      re_add <- paste(
        sprintf("(1 | %s)",
                vapply(context_to_add, maihda_quote_name, character(1))),
        collapse = " + ")
      formula <- stats::update(formula,
                               stats::as.formula(paste(". ~ . +", re_add)))
    }
    context_info <- list(context_vars = context)
  }

  # Convert family to family object if it's a string or constructor function.
  # "negbinomial" (the brms spelling) resolves to a plain marker list rather
  # than a stats family object: there is no theta-free negative-binomial family
  # constructor in stats -- lme4 estimates theta itself via glmer.nb() and brms
  # via its 'shape' parameter, so no theta is needed (or wanted) here.
  if (is.character(family)) {
    family <- switch(family,
                     gaussian = gaussian(),
                     binomial = binomial(),
                     poisson = poisson(),
                     negbinomial = list(family = "negbinomial", link = "log"),
                     ordinal = ,
                     cumulative = maihda_cumulative("logit"),
                     stop("Unsupported family: ", family))
  } else if (is.function(family)) {
    family <- family()
  }

  if (!is.list(family) || is.null(family$family) || is.null(family$link)) {
    stop("'family' must be a family name, family object, or family function.",
         call. = FALSE)
  }

  # Recompute on the RESOLVED family (the handshake above peeked at the raw
  # input) and validate the cumulative link: the latent-scale VPC is defined for
  # logit and probit only. Reached for maihda_cumulative()/brms::cumulative()
  # objects as well as the "ordinal"/"cumulative" strings.
  is_ordinal <- maihda_family_is_ordinal(family)
  if (is_ordinal) {
    maihda_ordinal_check_family(family)
  }

  # Negative binomial in any accepted form: the "negbinomial" marker from the
  # string path, brms::negbinomial (function or object), or a fixed-theta
  # MASS::negative.binomial(theta) object, whose label "Negative Binomial(<theta>)"
  # the normalizer maps to the canonical name. Only the log link is supported:
  # the latent-scale level-1 variance underlying the VPC (log1p(1/mu + 1/theta))
  # is derived for the log link, mirroring the poisson("log") restriction.
  is_negbin <- identical(maihda_normalize_family_name(family$family), "negbinomial")
  if (is_negbin && !identical(family$link, "log")) {
    stop("The negative-binomial family is only supported with the log link ",
         "(the latent-scale level-1 variance behind the VPC/ICC is defined for ",
         "it); this model uses link = '", family$link, "'.", call. = FALSE)
  }

  # Recode a two-level (Bernoulli) response to 0/1 (glmer and brms bernoulli both
  # accept 0/1). Aggregated binomial responses -- cbind(success, failure) or
  # `y | trials(n)` -- are left untouched and remain binomial() models.
  is_binomial_family <- family$family %in% c("binomial", "quasibinomial")
  response_is_binary <- is_binomial_family &&
    maihda_response_is_binary(formula, data, subset = subset_value,
                              weights = detect_weights)
  response_recoding <- NULL
  if (response_is_binary) {
    data <- maihda_prepare_binomial_response(data, formula, subset = subset_value,
                                             weights = detect_weights)
    # Mapping of original outcome levels to 0/1 (which level is the modeled event),
    # captured so it is inspectable on the returned model object.
    response_recoding <- attr(data, "response_recoding")
  }

  if (identical(engine, "wemix")) {
    # WeMix supports linear and binomial-logit models with the canonical single
    # (1 | stratum) intercept; reject anything else with a targeted message before
    # touching the engine. An aggregated binomial response (cbind / trials) has no
    # WeMix representation either.
    maihda_wemix_check_family(family)
    maihda_wemix_check_formula(formula)
    if (is_binomial_family && !response_is_binary) {
      stop("engine = \"wemix\" supports a binary (Bernoulli) 0/1 outcome only; ",
           "aggregated binomial responses (cbind(success, failure) or trials) ",
           "are not supported. Use engine = \"lme4\" or \"brms\".", call. = FALSE)
    }
  }

  # Build the engine call from the already-evaluated `...` values. Each value is
  # bound in a private environment and referenced by name so the model's stored
  # call stays small and readable (e.g. weights = .maihda_arg_weights) rather than
  # embedding whole vectors, and the pre-evaluated subset is a plain logical that no
  # longer depends on the (now recoded) response. The formula's environment is
  # pointed at this env so that lme4/brms, which evaluate `weights`/`subset` against
  # the formula's environment, find the bound values (and the `data` symbol).
  if (engine == "wemix") {
    # Design-weighted pseudo-ML fit. The guard above bans data-masked engine
    # arguments (weights/subset/offset) for wemix, so the remaining dots are plain
    # values WeMix::mix() takes directly -- the fit_env machinery the other
    # engines need is unnecessary here. maihda_fit_wemix() also pre-builds the
    # analytic sample (complete cases, positive weights), so the `data` it
    # returns matches the rows actually fitted.
    maihda_require_wemix()
    wemix_fit <- maihda_fit_wemix(formula, data, family, sampling_weights, dot_vals)
    model <- wemix_fit$model
    data <- wemix_fit$data
  } else if (engine == "ordinal") {
    # Cumulative link mixed model via ordinal::clmm(). Like the wemix branch,
    # the guard above bans data-masked engine arguments (weights/subset/offset),
    # so the remaining dots are plain values clmm() takes directly, and
    # maihda_fit_clmm() pre-builds the analytic sample (complete cases) so the
    # `data` it returns matches the rows actually fitted.
    maihda_require_ordinal()
    maihda_ordinal_check_formula(formula)
    data <- maihda_ordinal_prepare_response(data, formula)
    ord_fit <- maihda_fit_clmm(formula, data, family, dot_vals)
    model <- ord_fit$model
    data <- ord_fit$data
  } else {

  fit_env <- new.env(parent = environment(formula))
  fit_env$data <- data
  dot_args <- list()
  for (nm in names(dot_vals)) {
    bind_nm <- paste0(".maihda_arg_", nm)
    assign(bind_nm, dot_vals[[nm]], envir = fit_env)
    dot_args[[nm]] <- as.name(bind_nm)
  }
  environment(formula) <- fit_env

  # A negative-binomial request without a theta (the "negbinomial" marker from
  # the family string, or brms::negbinomial) means theta is to be ESTIMATED:
  # lme4 does that in glmer.nb(), brms via its 'shape' parameter. A fixed-theta
  # MASS::negative.binomial(theta) object instead is a complete GLM family that
  # plain glmer() accepts, so it takes the ordinary glmer path below, honouring
  # the user's theta.
  negbin_estimate_theta <- is_negbin &&
    identical(family$family, "negbinomial")

  if (engine == "lme4") {
    if (negbin_estimate_theta) {
      # glmer.nb() takes NO family argument -- it fits a Poisson glmer first,
      # ML-estimates theta, and refits with negative.binomial(theta). Its `...`
      # is forwarded to glmer(), so the bound dot_args (weights/subset/offset/
      # control, plus glmer.nb's own interval/tol/nb.control) pass through
      # unchanged via the fit_call below.
      fit_fun <- quote(lme4::glmer.nb)
      fit_args <- list(formula = formula, data = quote(data))
    } else {
      # Use lmer only for Gaussian with the identity link -- lmer takes no family
      # argument and silently ignores a non-identity link. Route a non-identity
      # Gaussian (e.g. gaussian(link = "log")) through glmer() so the link is
      # actually honoured, consistent with the family reported on the result.
      use_lmer <- family$family == "gaussian" && family$link == "identity"
      fit_fun <- if (use_lmer) quote(lme4::lmer) else quote(lme4::glmer)
      fit_args <- list(formula = formula, data = quote(data))
      if (!use_lmer) {
        fit_args$family <- family
      }
    }
  } else if (engine == "brms") {
    # Check if brms is installed
    if (!requireNamespace("brms", quietly = TRUE)) {
      stop("Package 'brms' is required but not installed. Please install it with: install.packages('brms')")
    }

    # Sampling weights enter the brms model as likelihood weights, which gives a
    # PSEUDO-posterior: point estimates are design-consistent, intervals are not
    # design-based. The weights are normalized to mean 1 so expansion weights do
    # not inflate the effective sample size, and rows with missing/non-positive
    # weights are dropped here so the stored data matches the fitted rows.
    if (!is.null(sampling_weights)) {
      prep <- maihda_prepare_brms_sampling_weights(data, formula, sampling_weights)
      data <- prep$data
      formula <- prep$formula
      fit_env$data <- data
      message("fit_maihda(): sampling weights enter the brms model as likelihood ",
              "weights (normalized to mean 1), giving a pseudo-posterior: point ",
              "estimates are design-consistent, but credible intervals are not ",
              "design-based -- interpret them cautiously.")
    }

    # brms models a 0/1 response with bernoulli(); passing binomial() would
    # require a trials specification and errors on Bernoulli data. Only rewrite
    # when the response really is a two-level vector -- aggregated binomial
    # (cbind / trials) must stay binomial().
    if (response_is_binary) {
      family <- brms::bernoulli(link = family$link)
    }

    if (is_negbin) {
      if (!negbin_estimate_theta) {
        # A fixed-theta MASS::negative.binomial(theta) object has no brms
        # counterpart (brms always estimates its 'shape' = theta).
        stop("A fixed-theta negative.binomial(theta) family object is only ",
             "supported by engine = \"lme4\". For brms, use family = ",
             "\"negbinomial\" (theta is estimated as the 'shape' parameter).",
             call. = FALSE)
      }
      # Convert the plain marker list from the family-string path into the
      # proper brmsfamily; a user-supplied brms::negbinomial object already is one.
      family <- brms::negbinomial(link = family$link)
    }

    if (is_ordinal) {
      # The same response validation/coercion the clmm path applies: brms's
      # cumulative() needs an ordered factor, and the category order is
      # load-bearing either way.
      data <- maihda_ordinal_prepare_response(data, formula)
      fit_env$data <- data
      family <- brms::cumulative(link = family$link)
    }

    fit_fun <- quote(brms::brm)
    fit_args <- list(formula = formula, data = quote(data),
                     family = family)
  }

  fit_call <- as.call(c(list(fit_fun), fit_args, dot_args))
  model <- eval(fit_call, fit_env)

  }

  # Capture fit-quality diagnostics (singular fit / non-convergence) so they can
  # be reported by print()/summary(); lme4 surfaces these only once at fit time.
  diagnostics <- maihda_fit_diagnostics(model)

  # Store the actual analytic model frame so downstream calculations use the
  # same rows as lme4/brms after their NA handling. The wemix and ordinal paths
  # pre-built their analytic samples above (complete cases), so `data` already
  # IS the fitted frame -- and model.frame() is undefined for WeMixResults,
  # while clmm's frame would drop the non-model columns (e.g. the stratum
  # dimension variables) the plots and group comparisons need.
  model_data <- if (engine %in% c("wemix", "ordinal")) {
    data
  } else {
    maihda_model_frame(model, fallback = data)
  }
  strata_info <- maihda_refresh_strata_counts(strata_info, model_data)
  attr(model_data, "strata_info") <- strata_info
  attr(model_data, "strata_vars") <- strata_vars
  attr(model_data, "strata_sep") <- strata_sep
  attr(model_data, "strata_autobin_info") <- strata_autobin_info

  result <- structure(
    list(
      model = model,
      engine = engine,
      formula = formula,
      data = model_data,
      original_data = data,
      family = family,
      strata_info = strata_info,
      strata_vars = strata_vars,
      strata_sep = strata_sep,
      strata_autobin_info = strata_autobin_info,
      context_vars = context,
      context_info = context_info,
      sampling_weights = sampling_weights,
      response_recoding = response_recoding,
      diagnostics = diagnostics
    ),
    class = "maihda_model"
  )

  return(result)
}

#' Print method for maihda_model
#'
#' @param x A maihda_model object
#' @param ... Additional arguments
#' @return No return value, called for side effects.
#' @export
print.maihda_model <- function(x, ...) {
  cat("MAIHDA Model\n")
  cat("============\n\n")
  cat("Engine:", x$engine, "\n")
  cat("Family:", x$family$family, "\n")
  cat("Formula:", deparse(x$formula), "\n")
  if (!is.null(x$context_vars)) {
    cat("Context:", paste(x$context_vars, collapse = ", "),
        "(crossed contextual random intercept)\n")
  }
  if (!is.null(x$sampling_weights)) {
    cat("Sampling weights:", x$sampling_weights,
        if (identical(x$engine, "wemix")) {
          "(design-weighted pseudo-maximum-likelihood)"
        } else {
          "(likelihood weights; pseudo-posterior)"
        }, "\n")
  }
  cat("\n")
  maihda_print_fit_diagnostics(x$diagnostics)
  cat("Underlying model:\n")
  print(x$model, ...)
  invisible(x)
}
