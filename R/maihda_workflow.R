#' Run a Complete MAIHDA Analysis
#'
#' A single high-level entry point that runs the standard two-model MAIHDA workflow
#' and returns one bundled object. It fits the \strong{null} model (covariates plus
#' the intersectional random intercept, \emph{excluding} the stratum dimensions' main
#' effects) and the \strong{adjusted} model (the null plus the additive main effects of
#' the stratum-defining dimensions), summarises the variance partition (VPC/ICC) of the
#' null model, and reports the \strong{PCV} -- the proportional change in
#' between-stratum variance from the null to the adjusted model, i.e. the additive
#' share of the intersectional inequality. When a higher-level grouping variable is
#' supplied it also compares this decomposition across that variable's levels.
#'
#' \strong{Binomial companions.} For a binary outcome the model summaries also carry
#' the discriminatory accuracy (AUC / C-statistic and Median Odds Ratio) -- the "DA"
#' in MAIHDA -- automatically, so the null model's strata-only AUC sits alongside its
#' VPC; set \code{response_vpc = TRUE} to add the (simulation-based) response-scale
#' VPC as well. These are read from \code{summary(x)} and the attached
#' \code{summary_adjusted}, and the headline \code{print()} shows the null-model AUC.
#'
#' This is a convenience wrapper around \code{\link{fit_maihda}},
#' \code{\link{calculate_pvc}}, \code{\link{summary.maihda_model}} and
#' \code{\link{compare_maihda_groups}}. It is \emph{intrinsically} a two-model
#' decomposition and has no single-model mode -- for a single fit (e.g. just the
#' null-model VPC / discriminatory accuracy), call \code{\link{fit_maihda}} directly.
#'
#' \strong{The dimensions' additive main effects.} You may write them in the formula --
#' the fully-specified, lme4-native adjusted model
#' \code{outcome ~ covars + var1 + var2 + (1 | var1:var2)} -- or omit them. Either way
#' the null \emph{excludes} the dimension main effects and the adjusted \emph{includes}
#' them: when the formula already lists them it is taken as the adjusted model and the
#' null is derived by dropping them; when they are missing \code{maihda()} adds them to
#' the adjusted model and emits a \code{message()} so the decomposition stays explicit.
#' The dimensions themselves are read from the random term: the shorthand
#' \code{(1 | var1:var2)} and \code{\link{make_strata}} both record them, and a numeric
#' dimension that \code{make_strata()} auto-binned enters the adjusted model as its
#' reconstructed tertile factor (matching the strata), not as a linear term. Because
#' \code{maihda()} is intrinsically a decomposition, it \strong{errors} (rather than
#' returning a null-only result) when it cannot build the adjusted model -- when the
#' dimensions cannot be recovered (a hand-built \code{stratum} column records none) or
#' there is only one dimension (no intersection to decompose). Use \code{\link{fit_maihda}}
#' for those single-model fits.
#'
#' @param formula A model formula, using either the intersectional shorthand
#'   \code{outcome ~ covars + (1 | var1:var2)} or \code{... + (1 | stratum)} when
#'   \code{data} already has a \code{stratum} column from \code{\link{make_strata}}.
#'   The dimensions' additive main effects may be listed (the fully-specified adjusted
#'   model) or omitted (added automatically, with a message); see Details.
#' @param data A data frame with the model variables (and the \code{group}
#'   variable, if used).
#' @param group Optional character string naming a higher-level grouping variable
#'   (e.g. \code{"country"}). When supplied, \code{\link{compare_maihda_groups}}
#'   is run and attached to the result. \code{group} runs a \emph{stratified}
#'   analysis (one independent model per level); to instead model a higher level
#'   \emph{jointly}, crossed with the strata, use \code{context}. The two are
#'   different designs, so supplying both errors.
#' @param context Optional character vector naming higher-level \emph{context}
#'   column(s) in \code{data} (e.g. \code{"school"}, \code{"hospital"},
#'   \code{"region"}), forwarded to \code{\link{fit_maihda}}. Each context enters
#'   every fitted model as a crossed random intercept --
#'   \code{outcome ~ covars + (1 | stratum) + (1 | context)} -- the
#'   \emph{contextual cross-classified MAIHDA} of the literature. The summaries
#'   then partition the unexplained variance into between-stratum vs.
#'   between-context vs. residual, the headline VPC/ICC becomes the between-stratum
#'   share \emph{net of} the context, and the PCV decomposition is computed with the
#'   context partialled out (the context random intercept is carried by both the
#'   null and the adjusted model). Cannot be combined with \code{group}; a context
#'   with few levels weakly identifies its variance (consider \code{engine = "brms"}).
#' @param engine Modeling engine, "lme4" (default), "brms", "wemix" (the
#'   design-weighted pseudo-maximum-likelihood fit; requires
#'   \code{sampling_weights} and is selected automatically when they are supplied
#'   with the default engine), or "ordinal" (cumulative link mixed model via
#'   \code{ordinal::clmm()}; selected automatically for an ordinal family or an
#'   ordered-factor outcome).
#' @param family Model family. Default "gaussian". As in \code{\link{fit_maihda}},
#'   a binary outcome is auto-detected when \code{family} is left at the default,
#'   and the same resolved family is then used for the group comparison so all
#'   models agree.
#' @param decomposition How to decompose the intersectional inequality into additive
#'   and interaction parts. \code{"two-model"} (default) is the classic MAIHDA
#'   approach: a null model and an adjusted model (the dimensions' additive main
#'   effects as \emph{fixed} effects), with the additive share read from the PCV
#'   (proportional change in between-stratum variance). \code{"crossed-dimensions"}
#'   fits a \strong{single} model that enters each dimension's additive main effect as
#'   a \emph{random} intercept -- \code{outcome ~ covars + (1 | dim1) + ... +
#'   (1 | stratum)} -- so each dimension's RE variance is its additive contribution and
#'   the intersection (\code{stratum}) RE variance is the interaction beyond additive;
#'   the additive and interaction \emph{shares} of the total between-strata variance are
#'   read directly from that one fit. The two modes target the same scientific question
#'   with different estimators, so their additive shares are conceptually parallel but
#'   not numerically identical. The crossed-dimensions additive share is a partial-pooling
#'   estimate: dimensions with few levels (e.g. a binary sex variable, whose variance is
#'   estimated from two groups) are poorly identified and often give a singular lme4 fit
#'   -- the \code{brms} engine handles this better. See Details.
#'   (\code{"cross-classified"} is accepted as a deprecated alias for
#'   \code{"crossed-dimensions"}, with a warning: in the MAIHDA literature
#'   \dQuote{cross-classified} refers to the contextual stratum-by-place model,
#'   which this package fits via \code{context}.)
#'   \code{"longitudinal"} fits a 3-level \strong{growth-curve} MAIHDA (requires
#'   \code{id} and \code{time}; selected automatically when they are supplied): a
#'   null and an adjusted growth model, where the adjusted model adds the
#'   dimensions' main effects \emph{and their interactions with time}
#'   (\code{dim:time}). The between-stratum variance is then time-varying and the
#'   PCV is reported separately for the baseline (intercept) and the slope variance
#'   -- the additive-vs-multiplicative split of the intersectional trajectory
#'   inequality (Bell, Evans, Holman & Leckie 2024). See \code{\link{fit_maihda}}.
#' @param id,time,time_degree For a \strong{longitudinal} MAIHDA: the person/unit
#'   identifier column, the numeric measurement-time column, and the growth-curve
#'   polynomial degree (1 = linear). Supplying \code{id}/\code{time} selects
#'   \code{decomposition = "longitudinal"}. See \code{\link{fit_maihda}} for the
#'   model structure; \code{group}, \code{context}, and \code{sampling_weights} are
#'   not supported alongside them. Default \code{NULL} (cross-sectional).
#' @param autobin Logical passed to \code{\link{make_strata}}; tertile-bins numeric
#'   grouping variables. Default TRUE.
#' @param shared_strata Logical, forwarded to \code{\link{compare_maihda_groups}}
#'   when \code{group} is supplied: build strata once on the full data so VPCs are
#'   comparable across groups (TRUE, default) or rebuild them within each group.
#' @param min_group_n Minimum group size for the per-group comparison, forwarded
#'   to \code{\link{compare_maihda_groups}}. Default 30.
#' @param bootstrap Logical; compute parametric-bootstrap VPC confidence intervals
#'   (lme4 only) for both the overall summary and the per-group comparison.
#'   Default FALSE.
#' @param n_boot Number of bootstrap samples when \code{bootstrap = TRUE}.
#' @param conf_level Confidence level for bootstrap intervals. Default 0.95.
#' @param response_vpc Logical; for a binomial (lme4) outcome, also attach the
#'   response-scale VPC (\code{\link{maihda_vpc_response}}) to the model summaries.
#'   It is estimated by simulation, so it is opt-in (default \code{FALSE}) and uses
#'   \code{seed}. The discriminatory accuracy (AUC + MOR) is attached automatically
#'   for a binomial/Bernoulli outcome regardless of this flag (see Details).
#' @param seed Optional integer seed for the response-scale VPC simulation.
#' @param sampling_weights Optional name of a sampling-weight column for a
#'   \strong{design-weighted MAIHDA} on complex-survey data; see
#'   \code{\link{fit_maihda}} for the full semantics (engine selection, the
#'   pseudo-likelihood weighting, and what is/is not design-based). Both the null
#'   and the adjusted model (and any per-group fits) use the same weights, so the
#'   PCV is a design-weighted decomposition. Not compatible with
#'   \code{engine = "lme4"}, \code{decomposition = "crossed-dimensions"} under
#'   the wemix engine, or \code{bootstrap = TRUE}.
#' @param interactions Whether to compute the per-stratum interaction diagnostic
#'   (\code{\link{maihda_interactions}}) on the adjusted / crossed-dimensions model
#'   and attach it as the \code{interactions} slot, surfaced in \code{print()}.
#'   \code{TRUE} (default) uses the diagnostic's own default correction
#'   (\code{adjust = "BH"}, FDR -- since scanning all strata is a screening
#'   question); \code{FALSE} skips it; or pass a \code{\link[stats]{p.adjust}} method
#'   name -- including \code{"none"} for the uncorrected, per-stratum
#'   individual-testing view. Uses \code{conf_level}. Not computed for a longitudinal
#'   decomposition. The computation is cheap (it reads the stratum estimates the
#'   summary already holds; no refit).
#' @param ... Additional arguments passed to \code{\link{fit_maihda}} (and on to
#'   \code{lmer}/\code{glmer}).
#'
#' @return An object of class \code{maihda_analysis}: a list with
#'   \item{model}{the fitted \code{maihda_model} (see \code{\link{fit_maihda}}); the
#'     \strong{null} model in \code{"two-model"} mode, or the single
#'     \strong{crossed-dimensions} model in \code{"crossed-dimensions"} mode}
#'   \item{summary}{the model's \code{maihda_summary} (VPC/ICC, variance components,
#'     stratum estimates; plus the additive/interaction \code{decomposition} in
#'     crossed-dimensions mode, and the stratum-vs-context \code{context} partition
#'     when \code{context} is supplied)}
#'   \item{model_adjusted}{the fitted \strong{adjusted} \code{maihda_model}
#'     (\code{"two-model"} and \code{"longitudinal"} modes; \code{NULL} otherwise)}
#'   \item{summary_adjusted}{the adjusted model's \code{maihda_summary}, or \code{NULL}}
#'   \item{pcv}{the proportional change in variance: the \code{pvc_result} from
#'     \code{\link{calculate_pvc}} in \code{"two-model"} mode, or a
#'     \code{maihda_long_pcv} (the intercept/slope and time-specific PCV) in
#'     \code{"longitudinal"} mode; \code{NULL} otherwise}
#'   \item{decomposition}{the additive/interaction partition (additive and interaction
#'     variances and shares, per-dimension variances; \code{"crossed-dimensions"} mode
#'     only, \code{NULL} otherwise)}
#'   \item{groups}{a \code{maihda_group_comparison} when \code{group} is supplied,
#'     otherwise \code{NULL}}
#'   \item{interactions}{the \code{maihda_interactions} diagnostic (per-stratum
#'     interaction BLUPs and flags) when \code{interactions} is not \code{FALSE},
#'     otherwise \code{NULL}}
#'   \item{mode}{\code{"two-model"}, \code{"crossed-dimensions"}, or \code{"longitudinal"}}
#'   \item{context_vars}{the context variable name(s) when \code{context} was
#'     supplied, otherwise \code{NULL}}
#'   \item{formula, adjusted_formula, group_var, call}{bookkeeping for printing}
#'
#' @seealso \code{\link{fit_maihda}} for the single-model fitter,
#'   \code{\link{compare_maihda_groups}} for the group comparison, and
#'   \code{\link{summary.maihda_model}} for the variance summary.
#'
#' @examples
#' \donttest{
#' data(maihda_health_data)
#'
#' # One call: null + adjusted fit, VPC summary, and PCV decomposition. Writing the
#' # dimensions' additive main effects (Gender + Race) gives the fully-specified
#' # adjusted model; maihda() derives the null by dropping them.
#' a <- maihda(BMI ~ Age + Gender + Race + (1 | Gender:Race), data = maihda_health_data)
#' a                              # VPC (null) and PCV (null -> adjusted)
#' a$pcv                          # proportional change in between-stratum variance
#' a$formula                      # null:     BMI ~ Age + (1 | stratum)
#' a$adjusted_formula             # adjusted: null + Gender + Race main effects
#'
#' # Omitting them is equivalent -- maihda() adds them to the adjusted model and
#' # emits a message; the null and PCV are identical to the explicit form above.
#' a0 <- maihda(BMI ~ Age + (1 | Gender:Race), data = maihda_health_data)
#'
#' plot(a, type = "vpc")          # null model
#' plot(a, type = "effect_decomp")# adjusted model (additive vs intersectional)
#'
#' # Crossed-dimensions decomposition: one model, the dimensions' main effects entered
#' # as RANDOM intercepts. The additive and interaction shares of the between-strata
#' # variance are read directly from the single fit (no null/adjusted pair).
#' cc <- maihda(BMI ~ Age + (1 | Gender:Race), data = maihda_health_data,
#'              decomposition = "crossed-dimensions")
#' cc                                    # VPC and additive/interaction shares
#' cc$decomposition$additive_share       # crossed-dimensions analogue of the PCV
#' cc$formula                            # BMI ~ Age + (1|Gender) + (1|Race) + (1|stratum)
#'
#' # Add a higher-level grouping variable to also compare across its levels.
#' # maihda_country_data has a real country grouping (PISA achievement data):
#' data(maihda_country_data)
#' a2 <- maihda(math ~ 1 + (1 | gender:ses), data = maihda_country_data,
#'              group = "country")
#' a2
#' plot(a2, type = "group_vpc")
#' plot(a2, type = "group_pcv")
#'
#' # Contextual cross-classified MAIHDA: instead of one model per country (group=),
#' # model the strata CROSSED with country in a single fit. The summary partitions
#' # the unexplained variance into between-stratum vs. between-country vs. residual,
#' # and the PCV is computed with country partialled out.
#' a3 <- maihda(math ~ 1 + (1 | gender:ses), data = maihda_country_data,
#'              context = "country")
#' a3
#' a3$summary$context$vpc_context_total  # the country (general contextual) share
#' plot(a3, type = "context_vpc")
#' }
#'
#' @export
maihda <- function(formula, data, group = NULL, context = NULL, engine = "lme4",
                   family = "gaussian",
                   decomposition = c("two-model", "crossed-dimensions",
                                     "longitudinal"),
                   autobin = TRUE, shared_strata = TRUE,
                   min_group_n = 30, bootstrap = FALSE,
                   n_boot = 1000, conf_level = 0.95,
                   response_vpc = FALSE, seed = NULL,
                   sampling_weights = NULL,
                   id = NULL, time = NULL, time_degree = 1,
                   interactions = TRUE, ...) {
  call <- match.call()

  # Longitudinal (growth-curve) MAIHDA selects itself when 'time' is supplied and
  # the decomposition is left at its default; an explicit non-longitudinal
  # decomposition with id/time is a contradiction.
  if (!is.null(time) || !is.null(id)) {
    if (missing(decomposition)) {
      decomposition <- "longitudinal"
    } else if (!identical(maihda_resolve_decomposition(decomposition),
                          "longitudinal")) {
      stop("'id'/'time' request a longitudinal MAIHDA, but decomposition = \"",
           maihda_resolve_decomposition(decomposition), "\" was given. Use ",
           "decomposition = \"longitudinal\" (or omit it) for a growth-curve fit.",
           call. = FALSE)
    }
  }
  decomposition <- maihda_resolve_decomposition(decomposition)
  if (identical(decomposition, "longitudinal") && (is.null(time) || is.null(id))) {
    stop("decomposition = \"longitudinal\" requires both 'id' (the person/unit ",
         "identifier) and 'time' (the measurement-time column). See ?maihda.",
         call. = FALSE)
  }
  if (identical(decomposition, "longitudinal") &&
      (!is.null(group) || !is.null(context) || !is.null(sampling_weights))) {
    stop("A longitudinal MAIHDA does not support 'group', 'context', or ",
         "'sampling_weights' (out of scope). Drop them for the growth-curve fit.",
         call. = FALSE)
  }

  # group= (stratified: one independent model per level) and context= (joint:
  # strata crossed with the higher level in ONE model) are different designs for
  # bringing in a higher level; combining them is ambiguous, so error early.
  if (!is.null(group) && !is.null(context)) {
    stop("Supply either 'group' (stratified comparison: one independent model per ",
         "level) or 'context' (contextual cross-classified model: strata crossed ",
         "with the higher level in one fit), not both.", call. = FALSE)
  }

  # Sampling weights select the design-weighted engine, mirroring fit_maihda():
  # the default engine switches to "wemix" (with a message), an explicit lme4 is
  # an error, and the wemix-incompatible workflow options fail early with
  # targeted messages instead of mid-decomposition.
  if (!is.null(sampling_weights)) {
    sampling_weights <- maihda_validate_sampling_weights(sampling_weights, data)
    if (missing(engine)) {
      engine <- "wemix"
      message("maihda(): 'sampling_weights' supplied; using engine = \"wemix\" ",
              "(design-weighted pseudo-maximum-likelihood via WeMix). Set 'engine' ",
              "explicitly to silence this message or to choose engine = \"brms\".")
    } else if (identical(engine, "lme4")) {
      stop("Sampling weights are not supported by engine = \"lme4\" (lme4's ",
           "weights are precision weights, not sampling weights). Use ",
           "engine = \"wemix\" or \"brms\".", call. = FALSE)
    }
  }
  if (identical(engine, "wemix")) {
    if (identical(decomposition, "crossed-dimensions")) {
      stop("decomposition = \"crossed-dimensions\" needs crossed random effects, ",
           "which WeMix does not fit. Use the default two-model decomposition ",
           "with engine = \"wemix\", or engine = \"brms\" for the ",
           "crossed-dimensions form.", call. = FALSE)
    }
    if (isTRUE(bootstrap)) {
      stop("Bootstrap intervals are not available for engine = \"wemix\" (a ",
           "design-based interval would require replicate weights). Set ",
           "bootstrap = FALSE.", call. = FALSE)
    }
  }

  # Ordinal (cumulative) family <-> engine handshake, mirroring fit_maihda():
  # the wrappers pass 'engine' explicitly to every fit, so fit_maihda()'s own
  # missing(engine) auto-switch could never fire through them. An ordered-factor
  # outcome under all-default family/engine likewise selects the ordinal engine
  # here (fit_maihda() then auto-detects the family on the analytic sample).
  if (missing(family) && missing(engine) && is.null(sampling_weights) &&
      isTRUE(tryCatch(maihda_response_is_ordinal(formula, data),
                      error = function(e) FALSE))) {
    engine <- "ordinal"
    message("maihda(): the outcome is an ordered factor; using the cumulative ",
            "(ordinal) model with engine = \"ordinal\" (ordinal::clmm). Specify ",
            "'family'/'engine' explicitly to override.")
  }
  is_ordinal <- maihda_family_is_ordinal(
    if (is.function(family)) tryCatch(family(), error = function(e) NULL) else family
  )
  if (is_ordinal) {
    if (missing(engine) && is.null(sampling_weights)) {
      engine <- "ordinal"
      message("maihda(): ordinal (cumulative) family; using engine = \"ordinal\" ",
              "(ordinal::clmm). Set 'engine' explicitly to silence this message ",
              "or to choose engine = \"brms\".")
    } else if (identical(engine, "lme4")) {
      stop("lme4 cannot fit a cumulative (ordinal) model. Use engine = ",
           "\"ordinal\" (ordinal::clmm, the default for this family) or ",
           "engine = \"brms\" (brms::cumulative).", call. = FALSE)
    }
  }
  if (identical(engine, "ordinal")) {
    if (identical(decomposition, "crossed-dimensions")) {
      stop("decomposition = \"crossed-dimensions\" needs crossed random effects, ",
           "which the ordinal (clmm) engine does not fit. Use the default ",
           "two-model decomposition, or engine = \"brms\" for the ",
           "crossed-dimensions form.", call. = FALSE)
    }
    if (isTRUE(bootstrap)) {
      stop("Bootstrap intervals are not available for engine = \"ordinal\" ",
           "(ordinal::clmm has no simulate()/refit() machinery). Set ",
           "bootstrap = FALSE, or use engine = \"brms\" for posterior credible ",
           "intervals.", call. = FALSE)
    }
  }

  # Fit the supplied formula first -- this resolves the stratum dimensions (and the
  # stratum column) and the family. Depending on whether the formula already lists the
  # dimensions' additive main effects, it serves as either the null or the adjusted
  # model below. When the user leaves 'family' at the default we omit it so fit_maihda()
  # can auto-detect a binary outcome, then reuse the resolved family for every model so
  # they agree.
  if (missing(family)) {
    model <- fit_maihda(formula, data, engine = engine, autobin = autobin,
                        context = context, sampling_weights = sampling_weights,
                        id = id, time = time, time_degree = time_degree, ...)
  } else {
    model <- fit_maihda(formula, data, engine = engine, family = family,
                        autobin = autobin, context = context,
                        sampling_weights = sampling_weights,
                        id = id, time = time, time_degree = time_degree, ...)
  }
  family_used <- model$family

  # maihda() IS the two-model decomposition -- it has no single-model mode. Building
  # the adjusted model requires at least two identifiable stratum dimensions, so error
  # (rather than return a null-only result) when that is impossible, pointing the user
  # at the shorthand / make_strata() or at fit_maihda() for a single fit. The usual
  # shorthand and make_strata() paths both record the dimensions; only a hand-built
  # 'stratum' column from a custom grouping records none.
  strata_vars <- model$strata_vars
  if (is.null(strata_vars) || length(strata_vars) == 0) {
    stop("maihda() builds an adjusted model from the stratum-defining variables to ",
         "compute the PCV, but they could not be identified (a pre-built 'stratum' ",
         "column from a custom grouping records none). Use the shorthand ",
         "(1 | var1:var2) or run make_strata() so the dimensions are recorded -- or ",
         "call fit_maihda() directly for a single-model fit.", call. = FALSE)
  }
  if (length(strata_vars) < 2) {
    stop("maihda() decomposes the intersectional inequality into additive and ",
         "interaction parts, which needs at least two stratum dimensions, but the ",
         "strata are defined by a single dimension (", paste(strata_vars, collapse = ", "),
         "). With one dimension there is no intersection to decompose; use ",
         "fit_maihda() for a single-dimension random-intercept fit.", call. = FALSE)
  }

  # --- Null and adjusted models ------------------------------------------------
  # The stratum dimensions' additive main effects belong in the ADJUSTED model. They
  # may be written explicitly in the formula (the fully-specified, lme4-native form)
  # or omitted; either way the null model excludes them and the adjusted includes
  # them. maihda_adjusted_terms() gives the correct term per dimension -- an
  # auto-binned numeric dimension enters as its reconstructed tertile factor
  # (.maihda_dim_*), matching the strata, not the raw column.
  dim_terms <- maihda_adjusted_terms(strata_vars, model$strata_autobin_info,
                                     model$original_data)$terms
  supplied_fixed <- attr(stats::terms(reformulas::nobars(model$formula)), "term.labels")
  # terms() backtick-quotes non-syntactic term labels (e.g. `gender var`), so the
  # dimension main effects must be compared in their quoted form -- a raw-name
  # intersect would miss them and misclassify a fully-specified adjusted formula as
  # the null (producing identical null/adjusted formulas and a corrupt PCV).
  # present_terms keeps the RAW names (remove_terms re-quotes them) and missing_vars
  # stays parallel to strata_vars.
  dim_terms_quoted <- vapply(dim_terms, maihda_quote_name, character(1))
  dim_present <- dim_terms_quoted %in% supplied_fixed
  present_terms <- dim_terms[dim_present]
  missing_vars <- strata_vars[!dim_present]

  remove_terms <- function(f, terms) {
    # maihda_quote_name() safely backtick-quotes each term, including the rare
    # legal name containing a backtick that a manual sprintf("`%s`") would break.
    quoted <- vapply(terms, maihda_quote_name, character(1))
    stats::update(f, stats::as.formula(
      paste(". ~ . -", paste(quoted, collapse = " - "))))
  }

  # --- Crossed-dimensions decomposition -----------------------------------------
  # A single model that estimates each dimension's additive main effect as a RANDOM
  # intercept (its variance is that dimension's additive contribution) and the
  # intersection as the interaction RE. This is the crossed-dimensions alternative
  # to the two-model fixed-effects PCV below. Any dimension main effects the user
  # wrote as fixed terms are stripped first -- here they enter as random effects.
  if (decomposition == "crossed-dimensions") {
    base_formula <- if (length(present_terms) > 0) {
      remove_terms(model$formula, present_terms)
    } else {
      model$formula
    }
    cc <- maihda_cross_classified_formula(base_formula, strata_vars,
                                          model$strata_autobin_info,
                                          model$original_data)
    # The formula builder strips ALL random effects before adding the dimension +
    # intersection intercepts; passing `context` again re-appends (and re-tags) any
    # contextual random intercept, so the two structures compose in one fit.
    cc_model <- fit_maihda(cc$formula, cc$data, engine = engine,
                           family = family_used, context = context,
                           sampling_weights = sampling_weights, ...)
    # Tag the fit so summary()/plot() read the additive-vs-interaction partition.
    cc_model$cc_info <- list(dim_groups = cc$dim_groups,
                             interaction_group = cc$interaction_group,
                             dim_labels = strata_vars)
    summary_obj <- summary(cc_model, bootstrap = bootstrap, n_boot = n_boot,
                           conf_level = conf_level, response_vpc = response_vpc,
                           seed = seed)

    groups <- NULL
    if (!is.null(group)) {
      group_formula <- if (length(present_terms) > 0) {
        remove_terms(formula, present_terms)
      } else {
        formula
      }
      groups <- compare_maihda_groups(
        group_formula, data, group = group, engine = engine, family = family_used,
        shared_strata = shared_strata, min_group_n = min_group_n,
        autobin = autobin, bootstrap = bootstrap, n_boot = n_boot,
        conf_level = conf_level, decomposition = "crossed-dimensions",
        sampling_weights = sampling_weights, ...
      )
    }

    cc_analysis <- structure(
      list(
        model = cc_model,
        summary = summary_obj,
        model_adjusted = NULL,
        summary_adjusted = NULL,
        pcv = NULL,
        decomposition = summary_obj$decomposition,
        groups = groups,
        formula = cc_model$formula,
        adjusted_formula = NULL,
        group_var = group,
        mode = "crossed-dimensions",
        context_vars = context,
        interactions = NULL,
        call = call
      ),
      class = "maihda_analysis"
    )
    return(maihda_attach_interactions(cc_analysis, interactions, conf_level))
  }

  # --- Longitudinal (growth-curve) decomposition --------------------------------
  # A 3-level growth pair: the NULL growth model (no dimension main effects) and an
  # ADJUSTED growth model that adds the dimensions' main effects AND their time
  # interactions (dim:time). The stratum-level intercept/slope variance remaining in
  # the adjusted model is the intersectional interaction beyond additive, so the PCV
  # in the baseline (intercept) and slope variances is the additive-vs-multiplicative
  # split of the intersectional trajectory inequality (Bell, Evans, Holman & Leckie
  # 2024). The VPC itself is time-varying (see summary()).
  if (decomposition == "longitudinal") {
    null_model <- if (length(present_terms) > 0) {
      fit_maihda(remove_terms(model$formula, present_terms), model$original_data,
                 engine = engine, family = family_used,
                 id = id, time = time, time_degree = time_degree, ...)
    } else {
      model
    }
    laf <- maihda_longitudinal_adjusted_formula(
      null_model$formula, strata_vars, null_model$strata_autobin_info,
      null_model$original_data, time = time, time_degree = time_degree)
    adjusted_model <- fit_maihda(laf$formula, laf$data, engine = engine,
                                 family = family_used, id = id, time = time,
                                 time_degree = time_degree, ...)

    summary_obj <- summary(null_model, bootstrap = bootstrap, n_boot = n_boot,
                           conf_level = conf_level)
    summary_adj <- summary(adjusted_model, bootstrap = bootstrap, n_boot = n_boot,
                           conf_level = conf_level)
    pcv <- tryCatch(
      maihda_longitudinal_pcv(null_model, adjusted_model),
      error = function(e) {
        warning("maihda(): the longitudinal PCV could not be computed (",
                conditionMessage(e), "). Returning the fitted null and adjusted ",
                "models without a PCV.", call. = FALSE)
        NULL
      })

    return(structure(
      list(
        model = null_model,
        summary = summary_obj,
        model_adjusted = adjusted_model,
        summary_adjusted = summary_adj,
        pcv = pcv,
        decomposition = NULL,
        groups = NULL,
        formula = null_model$formula,
        adjusted_formula = adjusted_model$formula,
        group_var = NULL,
        mode = "longitudinal",
        context_vars = NULL,
        interactions = NULL,
        call = call
      ),
      class = "maihda_analysis"
    ))
  }

  # In the branches below every derived formula already carries any contextual
  # random intercept (it is part of model$formula and survives remove_terms /
  # maihda_adjusted_formula); passing `context` again is idempotent in
  # fit_maihda() -- it only re-tags the fit so the summaries report the
  # stratum-vs-context partition. The null and adjusted models therefore share the
  # same context structure and the PCV is computed with the context partialled out.
  if (length(present_terms) == 0) {
    # No dimension main effects in the formula: the supplied fit IS the null; build
    # the adjusted by adding the dimensions' additive main effects.
    null_model <- model
    af <- maihda_adjusted_formula(null_model$formula, strata_vars,
                                  null_model$strata_autobin_info,
                                  null_model$original_data)
    adjusted_model <- fit_maihda(af$formula, af$data, engine = engine,
                                 family = family_used, context = context,
                                 sampling_weights = sampling_weights, ...)
    adjusted_formula <- af$formula
  } else if (length(missing_vars) == 0) {
    # Every dimension main effect is already present (the fully-specified adjusted
    # model): the supplied fit IS the adjusted; derive the null by removing them.
    adjusted_model <- model
    adjusted_formula <- model$formula
    null_model <- fit_maihda(remove_terms(model$formula, present_terms),
                             model$original_data, engine = engine,
                             family = family_used, context = context,
                             sampling_weights = sampling_weights, ...)
  } else {
    # Some dimension main effects present, some missing: fit a clean null (covariates
    # only) and a clean adjusted (all dimension main effects).
    null_model <- fit_maihda(remove_terms(model$formula, present_terms),
                             model$original_data, engine = engine,
                             family = family_used, context = context,
                             sampling_weights = sampling_weights, ...)
    af <- maihda_adjusted_formula(null_model$formula, strata_vars,
                                  null_model$strata_autobin_info,
                                  null_model$original_data)
    adjusted_model <- fit_maihda(af$formula, af$data, engine = engine,
                                 family = family_used, context = context,
                                 sampling_weights = sampling_weights, ...)
    adjusted_formula <- af$formula
  }

  # When maihda() supplied any dimension main effects itself, say so -- the
  # decomposition stays explicit rather than silent. Listing them in the formula (the
  # lme4-native form) makes the adjusted model self-documenting and suppresses this.
  if (length(missing_vars) > 0) {
    message("maihda(): added the additive main effect(s) of the stratum dimension(s) ",
            paste(missing_vars, collapse = ", "),
            " to the adjusted model; the null model excludes them. List them in the ",
            "formula to specify the adjusted model explicitly.")
  }

  summary_adj <- summary(adjusted_model, bootstrap = bootstrap, n_boot = n_boot,
                         conf_level = conf_level, response_vpc = response_vpc,
                         seed = seed)
  # A successfully fitted pair can still leave the PCV undefined when the null model
  # has zero between-stratum variance (a boundary/singular fit); keep both models and
  # warn rather than aborting in that numerical edge case.
  pcv <- tryCatch(
    calculate_pvc(null_model, adjusted_model, bootstrap = bootstrap, n_boot = n_boot,
                  conf_level = conf_level),
    error = function(e) {
      warning("maihda(): the PCV could not be computed (", conditionMessage(e),
              "). Returning the fitted null and adjusted models without a PCV.",
              call. = FALSE)
      NULL
    })

  summary_obj <- summary(null_model, bootstrap = bootstrap, n_boot = n_boot,
                         conf_level = conf_level, response_vpc = response_vpc,
                         seed = seed)

  groups <- NULL
  if (!is.null(group)) {
    # Strip any dimension main effects the user wrote so the per-group decomposition
    # matches the overall one: each group's null excludes them and its adjusted adds
    # them. (compare_maihda_groups() builds the per-group strata and adjusted itself.)
    group_formula <- if (length(present_terms) > 0) remove_terms(formula, present_terms) else formula
    groups <- compare_maihda_groups(
      group_formula, data, group = group, engine = engine, family = family_used,
      shared_strata = shared_strata, min_group_n = min_group_n,
      autobin = autobin, bootstrap = bootstrap, n_boot = n_boot,
      conf_level = conf_level, sampling_weights = sampling_weights, ...
    )
  }

  analysis <- structure(
    list(
      model = null_model,
      summary = summary_obj,
      model_adjusted = adjusted_model,
      summary_adjusted = summary_adj,
      pcv = pcv,
      decomposition = NULL,
      groups = groups,
      formula = null_model$formula,
      adjusted_formula = adjusted_formula,
      group_var = group,
      mode = "two-model",
      context_vars = context,
      interactions = NULL,
      call = call
    ),
    class = "maihda_analysis"
  )
  maihda_attach_interactions(analysis, interactions, conf_level)
}

# Compact discriminatory-accuracy line for the analysis headline: the null /
# strata-only model's AUC + MOR (Merlo's central DA quantity), the adjusted model's
# AUC for comparison, and the response-scale VPC when present. A no-op for a
# non-binomial outcome, where the summaries carry no `discriminatory_accuracy`.
maihda_print_analysis_da <- function(summary_obj, summary_adjusted = NULL) {
  da <- summary_obj$discriminatory_accuracy
  if (is.null(da)) return(invisible(NULL))
  fmt <- function(v, d = 3) {
    if (isTRUE(is.finite(v))) formatC(v, format = "f", digits = d) else "NA"
  }
  pal <- maihda_palette()
  cat("\n", pal$bold("Discriminatory accuracy (null model):"), "\n", sep = "")
  cat(sprintf("  AUC: %s | MOR: %s | cases/controls: %d/%d\n",
              pal$accent(fmt(da$auc)), pal$accent(fmt(da$mor)), da$n_case, da$n_control))
  da_adj <- if (!is.null(summary_adjusted)) summary_adjusted$discriminatory_accuracy else NULL
  if (!is.null(da_adj) && isTRUE(is.finite(da_adj$auc))) {
    cat(sprintf("  Adjusted-model AUC: %s\n", pal$accent(fmt(da_adj$auc))))
  }
  vr <- summary_obj$vpc_response
  if (!is.null(vr) && isTRUE(is.finite(vr$estimate))) {
    cat(sprintf("  Response-scale VPC (null): %s\n", pal$accent(fmt(vr$estimate, 4))))
  }
  invisible(NULL)
}

#' Print a MAIHDA Analysis
#'
#' @param x A \code{maihda_analysis} object from \code{\link{maihda}}.
#' @param ... Additional arguments (not used).
#' @return No return value, called for side effects.
#' @export
print.maihda_analysis <- function(x, ...) {
  pal <- maihda_palette()
  cat(pal$bold("MAIHDA Analysis"), "\n", sep = "")
  cat("===============\n\n")

  # Crossed-dimensions mode: one model, additive (dimension REs) vs interaction
  # (intersection RE) read directly. Distinct layout from the two-model PCV below.
  if (identical(x$mode, "crossed-dimensions")) {
    cat("Decomposition:   crossed-dimensions (single model)\n")
    cat("Formula:         ", paste(deparse(x$formula), collapse = " "), "\n", sep = "")
    cat("Engine: ", x$model$engine, " | Family: ", x$model$family$family, "\n", sep = "")
    if (!is.null(x$context_vars)) {
      cat("Context:         ", paste(x$context_vars, collapse = ", "),
          " (crossed contextual random intercept)\n", sep = "")
    }
    maihda_print_fit_diagnostics(x$model$diagnostics)

    vpc <- x$summary$vpc
    if (maihda_vpc_has_interval(vpc)) {
      cat(sprintf("VPC/ICC: %s [%.4f, %.4f]\n",
                  pal$accent(sprintf("%.4f", vpc$estimate)), vpc$ci_lower, vpc$ci_upper))
    } else {
      cat(sprintf("VPC/ICC: %s\n", pal$accent(sprintf("%.4f", vpc$estimate))))
    }
    if (!is.null(x$decomposition)) {
      cat("\n")
      maihda_print_cc_decomposition(x$decomposition)
    }
    if (!is.null(x$summary$context)) {
      maihda_print_context_partition(x$summary$context)
    }
    maihda_print_analysis_da(x$summary, x$summary_adjusted)
    if (!is.null(x$summary$stratum_estimates)) {
      cat("Strata: ", nrow(x$summary$stratum_estimates), "\n", sep = "")
    }
    maihda_print_interactions_line(x$interactions)
    if (!is.null(x$groups)) {
      cat(sprintf("\nGroup comparison by '%s':\n", x$group_var))
      print(x$groups)
    }
    cat(pal$muted("\nUse summary() for variance components and plot(type = ...) for figures.\n"))
    return(invisible(x))
  }

  # Longitudinal mode: a null/adjusted growth pair with a time-varying VPC and the
  # additive-vs-multiplicative PCV (a maihda_long_pcv, printed by its own method).
  if (identical(x$mode, "longitudinal")) {
    lng <- x$summary$longitudinal
    cat("Decomposition:   longitudinal (3-level growth curve)\n")
    cat("Null formula:    ", paste(deparse(x$formula), collapse = " "), "\n", sep = "")
    cat("Adjusted formula:", paste(deparse(x$adjusted_formula), collapse = " "),
        "\n", sep = "")
    cat("Engine: ", x$model$engine, " | Family: ", x$model$family$family, "\n",
        sep = "")
    maihda_print_fit_diagnostics(x$model$diagnostics)

    vpc <- x$summary$vpc
    if (maihda_vpc_has_interval(vpc)) {
      cat(sprintf("VPC/ICC at baseline (%s = %g): %s [%.4f, %.4f]\n",
                  lng$time, lng$ref_time, pal$accent(sprintf("%.4f", vpc$estimate)),
                  vpc$ci_lower, vpc$ci_upper))
    } else {
      cat(sprintf("VPC/ICC at baseline (%s = %g): %s\n",
                  lng$time, lng$ref_time, pal$accent(sprintf("%.4f", vpc$estimate))))
    }
    if (!is.null(lng)) {
      vt <- lng$vpc_t
      cat(sprintf("VPC/ICC over %s: %.4f to %.4f (time-varying).\n", lng$time,
                  min(vt$estimate, na.rm = TRUE), max(vt$estimate, na.rm = TRUE)))
    }
    if (!is.null(x$pcv)) {
      cat("\n")
      print(x$pcv)
    }
    if (!is.null(x$summary$stratum_estimates)) {
      cat("\nStrata: ", nrow(x$summary$stratum_estimates), "\n", sep = "")
    }
    cat(pal$muted(paste0("\nUse summary() for variance components and ",
        "plot(type = \"vpc_trajectory\" / \"trajectories\") for figures.\n")))
    return(invisible(x))
  }

  cat("Null formula:    ", paste(deparse(x$formula), collapse = " "), "\n", sep = "")
  if (!is.null(x$adjusted_formula)) {
    cat("Adjusted formula:", paste(deparse(x$adjusted_formula), collapse = " "), "\n", sep = "")
  }
  cat("Engine: ", x$model$engine, " | Family: ", x$model$family$family, "\n", sep = "")
  if (!is.null(x$context_vars)) {
    cat("Context: ", paste(x$context_vars, collapse = ", "),
        " (crossed contextual random intercept in the null and adjusted models)\n",
        sep = "")
  }

  vpc <- x$summary$vpc
  if (maihda_vpc_has_interval(vpc)) {
    cat(sprintf("VPC/ICC (null): %s [%.4f, %.4f]\n",
                pal$accent(sprintf("%.4f", vpc$estimate)), vpc$ci_lower, vpc$ci_upper))
  } else {
    cat(sprintf("VPC/ICC (null): %s\n", pal$accent(sprintf("%.4f", vpc$estimate))))
  }
  if (!is.null(x$summary$context)) {
    ctx <- x$summary$context
    cat(sprintf("Context share (null): %.4f (between-%s share of unexplained variance)\n",
                ctx$vpc_context_total, paste(ctx$context_vars, collapse = "/")))
  }

  if (!is.null(x$pcv)) {
    pcv <- x$pcv
    if (isTRUE(pcv$bootstrap) && !is.null(pcv$ci_lower)) {
      cat(sprintf("PCV (null -> adjusted): %s [%.4f, %.4f]\n",
                  pal$accent(sprintf("%.4f", pcv$pvc)), pcv$ci_lower, pcv$ci_upper))
    } else {
      cat(sprintf("PCV (null -> adjusted): %s\n", pal$accent(sprintf("%.4f", pcv$pvc))))
    }
    cat(sprintf("Between-stratum variance: %.4f (null) -> %.4f (adjusted)\n",
                pcv$var_model1, pcv$var_model2))
    if (pcv$pvc >= 0) {
      cat(sprintf(paste0("  ~%.1f%% of the between-stratum variance is additive (the ",
                         "dimensions' main\n  effects); the remainder is the between-stratum ",
                         "variance remaining after the\n  additive main effects -- a ",
                         "model-dependent quantity\n"),
                  pcv$pvc * 100))
    } else {
      cat("  PCV < 0: the additive main effects do not account for the between-stratum\n",
          "  variance (possible suppression/rescaling).\n", sep = "")
    }
  }

  maihda_print_analysis_da(x$summary, x$summary_adjusted)

  if (!is.null(x$summary$stratum_estimates)) {
    cat("Strata: ", nrow(x$summary$stratum_estimates), "\n", sep = "")
  }

  maihda_print_interactions_line(x$interactions)

  if (!is.null(x$groups)) {
    cat(sprintf("\nGroup comparison by '%s':\n", x$group_var))
    print(x$groups)
  }

  cat(pal$muted("\nUse summary() for variance components and plot(type = ...) for figures.\n"))
  invisible(x)
}

#' Summarize a MAIHDA Analysis
#'
#' Returns the variance summary (VPC/ICC, variance components, stratum estimates)
#' of the fitted model. The per-group comparison, when present, is attached as the
#' \code{"groups"} attribute.
#'
#' @param object A \code{maihda_analysis} object from \code{\link{maihda}}.
#' @param ... Additional arguments (not used).
#' @return The \code{maihda_summary} for the fitted model.
#' @export
summary.maihda_analysis <- function(object, ...) {
  out <- object$summary
  attr(out, "groups") <- object$groups
  attr(out, "pcv") <- object$pcv
  attr(out, "adjusted") <- object$summary_adjusted
  attr(out, "interactions") <- object$interactions
  out
}

#' Plot a MAIHDA Analysis
#'
#' Dispatches each \code{type} to the model it is valid on. The VPC and shrinkage
#' views (\code{"vpc"}, \code{"obs_vs_shrunken"}, \code{"predicted"}) use the
#' \strong{null} model. The additive-vs-intersectional views (\code{"risk_vs_effect"},
#' \code{"effect_decomp"}, \code{"ternary"}, \code{"prediction_deviation"}) use the
#' \strong{adjusted} model, whose fixed effects carry the dimensions' additive part so
#' the stratum random effect is the pure interaction; with fewer than two dimensions
#' (no adjusted model) they fall back to the null model. Group types
#' (\code{"group_vpc"}, \code{"group_components"}, \code{"group_between_variance"},
#' \code{"group_pcv"}, \code{"group_additive_share"}) use the group comparison when
#' \code{\link{maihda}} was called with a \code{group}.
#'
#' @param x A \code{maihda_analysis} object from \code{\link{maihda}}.
#' @param type One of the model types ("all", "vpc", "obs_vs_shrunken", "predicted",
#'   "risk_vs_effect", "effect_decomp", "ternary", "prediction_deviation"), the
#'   contextual type ("context_vpc", a stratum-vs-context variance bar; requires
#'   \code{maihda(context = )}), a longitudinal type ("vpc_trajectory",
#'   "trajectories", "pcv_trajectory"; requires \code{decomposition =
#'   "longitudinal"}), or a group type ("group_vpc", "group_components",
#'   "group_between_variance", "group_pcv", "group_additive_share"). Default "all".
#'   For a longitudinal analysis "all" shows the VPC-over-time, the stratum
#'   trajectories, and the time-specific PCV.
#' @param highlight_interactions Highlight strata with a credibly non-zero
#'   intersectional interaction on the BLUP-based views (see
#'   \code{\link{maihda_interactions}} and \code{\link[=plot.maihda_model]{plot}}).
#'   \code{FALSE} (default), \code{TRUE} (computed from this analysis's adjusted /
#'   crossed-dimensions model), a multiple-testing method such as \code{"BH"}, or
#'   a \code{maihda_interactions} object. The flags are computed once from the
#'   correct (adjusted) model and reused across views.
#' @param ... Additional arguments passed to the underlying plot method.
#' @return A ggplot2 object, or (for \code{type = "all"}) an invisible list of them.
#' @export
plot.maihda_analysis <- function(x, type = "all", highlight_interactions = FALSE, ...) {
  type <- match.arg(type, c(
    "all", "vpc", "obs_vs_shrunken", "predicted", "risk_vs_effect",
    "effect_decomp", "ternary", "prediction_deviation", "context_vpc",
    "vpc_trajectory", "trajectories", "pcv_trajectory",
    "group_vpc", "group_components", "group_between_variance", "group_pcv",
    "group_additive_share"
  ))

  # Longitudinal analysis: the trajectory views replace the cross-sectional ones.
  # "all" shows the VPC-over-time and the stratum mean trajectories; "pcv_trajectory"
  # plots the time-specific additive share from the null/adjusted pair.
  if (identical(x$mode, "longitudinal")) {
    if (type == "pcv_trajectory") {
      if (is.null(x$pcv)) {
        stop("No longitudinal PCV is available for this analysis.", call. = FALSE)
      }
      return(plot_pcv_trajectory(x$pcv))
    }
    if (type == "all") {
      plots <- list(
        vpc_trajectory = plot(x$model, type = "vpc_trajectory",
                              summary_obj = x$summary, ...),
        trajectories = tryCatch(plot(x$model, type = "trajectories",
                                     summary_obj = x$summary, ...),
                                error = function(e) NULL),
        pcv_trajectory = tryCatch(plot_pcv_trajectory(x$pcv),
                                  error = function(e) NULL)
      )
      for (p in plots[!vapply(plots, is.null, logical(1))]) print(p)
      return(invisible(plots))
    }
    if (type %in% c("vpc_trajectory", "trajectories", "vpc")) {
      t <- if (type == "vpc") "vpc_trajectory" else type
      return(plot(x$model, type = t, summary_obj = x$summary, ...))
    }
    # Other cross-sectional types fall through to the adjusted growth model below.
  }

  if (type == "context_vpc" && is.null(x$context_vars)) {
    stop("No contextual partition is available. Call maihda() with a 'context' ",
         "argument (the contextual cross-classified model).", call. = FALSE)
  }

  group_types <- c("group_vpc", "group_components", "group_between_variance",
                   "group_pcv", "group_additive_share")
  if (type %in% group_types) {
    if (is.null(x$groups)) {
      stop("No group comparison is available. Call maihda() with a 'group' argument.",
           call. = FALSE)
    }
    gtype <- sub("^group_", "", type)
    return(plot(x$groups, type = gtype))
  }

  # The additive-vs-intersectional views are only interpretable on the adjusted
  # model -- its fixed effects (two-model) or dimension random effects
  # (crossed-dimensions) carry the dimensions' additive part, so the stratum random
  # effect is the pure interaction. The VPC, shrinkage, and context views belong to
  # the null / crossed-dimensions model. In crossed-dimensions mode there is no
  # separate adjusted model, so all views use x$model (which already carries the
  # additive structure).
  adjusted_types <- c("risk_vs_effect", "effect_decomp", "ternary", "prediction_deviation")
  adj_model <- if (!is.null(x$model_adjusted)) x$model_adjusted else x$model
  adj_summary <- if (!is.null(x$model_adjusted)) x$summary_adjusted else x$summary

  # Resolve the interaction highlight ONCE, from the analysis (so it reads the
  # adjusted / crossed-dimensions model and never trips the null-model guardrail),
  # then forward the resolved maihda_interactions object to every model view -- the
  # null-model shrinkage views included -- so the same flagged strata are marked
  # everywhere without recomputation.
  hl <- maihda_resolve_analysis_highlight(x, highlight_interactions)

  if (type == "all") {
    null_plots <- list(vpc = plot(x$model, type = "vpc", summary_obj = x$summary, ...))
    null_plots$obs_vs_shrunken <- tryCatch(
      plot(x$model, type = "obs_vs_shrunken", summary_obj = x$summary,
           highlight_interactions = hl, ...),
      error = function(e) NULL)
    null_plots$predicted <- tryCatch(
      plot(x$model, type = "predicted", summary_obj = x$summary,
           highlight_interactions = hl, ...),
      error = function(e) NULL)
    if (!is.null(x$context_vars)) {
      null_plots$context_vpc <- tryCatch(
        plot(x$model, type = "context_vpc", summary_obj = x$summary, ...),
        error = function(e) NULL)
    }

    adj_plots <- list()
    for (t in adjusted_types) {
      adj_plots[[t]] <- tryCatch(
        plot(adj_model, type = t, summary_obj = adj_summary,
             highlight_interactions = hl, ...),
        error = function(e) NULL)
    }

    model_plots <- c(null_plots, adj_plots)
    if (!is.null(x$groups)) {
      group_plots <- list(
        group_vpc = plot(x$groups, type = "vpc"),
        group_components = plot(x$groups, type = "components"),
        group_between_variance = plot(x$groups, type = "between_variance")
      )
      group_plots$group_pcv <- tryCatch(plot(x$groups, type = "pcv"),
                                         error = function(e) NULL)
      model_plots <- c(model_plots, group_plots)
    }
    for (p in model_plots[!vapply(model_plots, is.null, logical(1))]) print(p)
    return(invisible(model_plots))
  }

  if (type %in% adjusted_types) {
    return(plot(adj_model, type = type, summary_obj = adj_summary,
                highlight_interactions = hl, ...))
  }

  # vpc, obs_vs_shrunken, predicted -> null model
  plot(x$model, type = type, summary_obj = x$summary,
       highlight_interactions = hl, ...)
}

# Resolve the highlight argument for a maihda_analysis: FALSE/NULL stays FALSE;
# a maihda_interactions object passes through; TRUE reuses the interaction
# diagnostic already stored on the analysis (so the plot matches the printed
# headline exactly), falling back to computing it when the analysis was built with
# interactions = FALSE. A p.adjust method name such as "BH" computes adjusted
# flags from the analysis's adjusted/crossed-dimensions model. Either way the
# downstream model plots receive a ready object and neither recompute nor warn.
maihda_resolve_analysis_highlight <- function(x, highlight_interactions) {
  if (is.null(highlight_interactions) || isFALSE(highlight_interactions)) {
    return(FALSE)
  }
  if (inherits(highlight_interactions, "maihda_interactions")) {
    return(highlight_interactions)
  }
  if (isTRUE(highlight_interactions)) {
    if (inherits(x$interactions, "maihda_interactions")) {
      return(x$interactions)
    }
    return(maihda_interactions(x))
  }
  if (is.character(highlight_interactions) && length(highlight_interactions) == 1L) {
    choices <- c("none", stats::p.adjust.methods)
    if (!highlight_interactions %in% choices) {
      stop("'highlight_interactions' must be FALSE, TRUE, a multiple-comparison ",
           "method name (e.g. \"BH\"), or a maihda_interactions object from ",
           "maihda_interactions().", call. = FALSE)
    }
    return(maihda_interactions(x, adjust = highlight_interactions))
  }
  stop("'highlight_interactions' must be FALSE, TRUE, a multiple-comparison ",
       "method name (e.g. \"BH\"), or a maihda_interactions object from ",
       "maihda_interactions().", call. = FALSE)
}
