#' Fit MAIHDA Model
#'
#' Fits a multilevel model for MAIHDA (Multilevel Analysis of Individual
#' Heterogeneity and Discriminatory Accuracy) using lme4, brms, WeMix (for
#' design-weighted (survey) data), or -- for an ordered-factor outcome -- a
#' cumulative link mixed model via \code{ordinal::clmm()}.
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
#'   for \code{gaussian("identity")}, the binomial/Bernoulli families with a logit,
#'   probit, or complementary log-log (\code{cloglog}) link (latent level-1 variance
#'   \eqn{\pi^2/3}, 1, and \eqn{\pi^2/6} respectively), \code{poisson("log")}, and the negative binomial with a log
#'   link (level-1 variance \code{log(1 + 1/lambda + 1/theta)} at the
#'   \emph{marginal} expected count \code{lambda}; Nakagawa, Johnson & Schielzeth
#'   2017 -- see \code{\link{summary.maihda_model}} for how \code{lambda} is
#'   computed). Other families (for example \code{Gamma(link = "log")})
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
#'   weights, for a \strong{weighted MAIHDA} on survey data. Sampling weights are
#'   not the same thing as lme4's \code{weights=} (precision weights, which
#'   rescale the residual variance), so combining \code{sampling_weights} with
#'   \code{engine = "lme4"} is an error.
#'
#'   \strong{What this does and does not buy you.} A single person-level weight
#'   column is all this argument accepts: there is no representation of the
#'   sampling hierarchy -- no primary sampling units, no sampling strata, no
#'   higher-stage or conditional weights, no finite-population corrections, and no
#'   replicate weights. Weighting the likelihood makes the point estimates target
#'   the population-weighted estimand, which is the main reason to use it. The
#'   standard errors, however, are \emph{not} general design-based standard errors
#'   for a multistage or stratified sample, because none of the design features
#'   that drive them are supplied. Treat the intervals as approximate unless your
#'   design genuinely matches the assumptions below; for full design-based
#'   inference on such data, use replicate weights with a survey-specific package.
#'
#'   Two engines support them:
#'   \itemize{
#'     \item \code{engine = "wemix"} (chosen automatically when \code{engine} is
#'       left at its default): weighted pseudo-maximum-likelihood via
#'       \code{WeMix::mix()} (Rabe-Hesketh & Skrondal 2006), the estimator used
#'       for NAEP/PISA analysis. The individual weights enter at level 1
#'       unchanged and the level-2 (stratum) weights are 1, because
#'       intersectional strata are treated as exhaustive population cells
#'       included with certainty -- so this is a \emph{single-stage} weighted
#'       model, and it is that assumption the inference rests on. Supports
#'       \code{gaussian(identity)} and \code{binomial(logit)} models with the
#'       canonical single \code{(1 | stratum)} random intercept. Fixed-effect
#'       standard errors are the sandwich (robust) errors WeMix reports, which
#'       account for the weighting and for dependence within the model's own
#'       grouping (the intersectional strata) -- but not for clustering or
#'       stratification induced by the sample design, which the strata do not
#'       represent. The VPC/PCV are reported as point estimates (no bootstrap --
#'       see \code{\link{summary.maihda_model}}).
#'     \item \code{engine = "brms"}: the weights enter the model as likelihood
#'       weights (\code{y | weights(w)}), normalized to mean 1, giving a
#'       \emph{pseudo-posterior}: point estimates target the population-weighted
#'       estimand but credible intervals are not design-based -- interpret them
#'       cautiously. Full design-based (replicate-weight / linearised) variances
#'       for the variance components are not computed.
#'   }
#'   Rows with a missing or non-positive sampling weight are dropped with a
#'   warning. The column names \code{.maihda_sw} (brms likelihood weight) and
#'   \code{.maihda_l2wt} (WeMix level-2 weight) are \strong{reserved} for the
#'   design-weighted engines: they are written into the analytic data internally,
#'   so a model that supplies \code{sampling_weights} may not also reference a
#'   variable of either name in its formula (fitting would error). Default
#'   \code{NULL} (unweighted).
#' @param id Optional single character string naming a person/unit identifier
#'   column for a \strong{longitudinal (growth-curve) MAIHDA} on long-format data
#'   (one row per measurement occasion). Id values must be \emph{globally}
#'   unique to a person -- ids numbered within a site or group (person "1" in
#'   every site) would merge different people's trajectories, and an id
#'   appearing in more than one stratum is rejected with an error. Supplied
#'   together with \code{time}, it
#'   makes the model a 3-level growth curve -- occasions within individuals
#'   (\code{id}) within intersectional strata -- with a random intercept and slope
#'   on \code{time} at \emph{both} the individual and stratum levels. The growth
#'   random effects are added automatically: write the strata shorthand
#'   \code{(1 | var1:var2)} (or \code{(1 | stratum)}) only, not the slopes. The
#'   between-stratum variance (and hence the VPC) then becomes a function of time;
#'   \code{\link{summary.maihda_model}} reports the time-varying VPC. Longitudinal
#'   fits are supported by \code{engine = "lme4"}/\code{"brms"} only (not
#'   \code{wemix}/\code{ordinal}), and are incompatible with \code{context} and
#'   \code{sampling_weights}. Default \code{NULL} (cross-sectional). See
#'   Bell, Evans, Holman & Leckie (2024).
#' @param time Optional single character string naming a numeric measurement-time
#'   column (e.g. wave 0, 1, 2, ... or age), required for a longitudinal MAIHDA;
#'   see \code{id}. When the time axis does not start at 0 (age, calendar year,
#'   waves coded 10, 11, ...), the growth terms are fit on internally
#'   \emph{centered} time (\code{time - min(time)}, with a message): the raw
#'   polynomial basis over an offset range is ill-conditioned and can silently
#'   converge to a wrong solution. All results (the time-varying VPC, the
#'   PCV, plots, predictions) are reported on the original \code{time} scale;
#'   the column name \code{.maihda_ctime} is reserved for the internal centered
#'   variable. Default \code{NULL}.
#' @param time_degree Polynomial degree of the growth curve when \code{time} is
#'   supplied: 1 (default) linear, 2 quadratic, etc. The brms engine supports
#'   degree 1 only.
#' @param stratum_slope Longitudinal only: keep the stratum-level random slope(s)
#'   on \code{time}? \code{TRUE} (default) fits the canonical Bell et al. (2024)
#'   structure, \code{(time | id) + (time | stratum)}, in which the between-stratum
#'   variance is a function of time. \code{FALSE} fits \code{(time | id) +
#'   (1 | stratum)}: the individual level keeps its growth block, but strata differ
#'   in \emph{level} only, so the between-stratum variance -- and the numerator of
#'   the VPC -- is constant over time and no \code{PCV_slope} is defined. Use it
#'   when the stratum slope variance is at the singularity boundary: with few
#'   strata, few occasions per stratum, or irregular measurement times, a
#'   \code{(time | stratum)} block routinely collapses to a perfect intercept-slope
#'   correlation, and the trajectory decomposition it supports is then a boundary
#'   artefact. The VPC still varies with time through the person-level slope
#'   variance and the residual, so this is a time-constant between-stratum
#'   \emph{variance}, not a time-constant VPC.
#' @param interactions Opt-in per-stratum interaction diagnostic
#'   (\code{\link{maihda_interactions}}), attached as the \code{interactions} slot
#'   and shown by \code{print()}. \code{FALSE} (default) skips it; \code{TRUE}
#'   computes it with the diagnostic's default correction (\code{adjust = "BH"}); or
#'   pass a \code{\link[stats]{p.adjust}} method name, including \code{"none"} for the
#'   uncorrected view. It is meaningful only on an \emph{adjusted} model (the
#'   dimensions' main effects in the fixed part); on a null model
#'   \code{maihda_interactions} warns. This is the single-fit parallel to the
#'   default-on \code{interactions} of \code{\link{maihda}}.
#' @param count_approximation Which latent-scale level-1 (observation) variance
#'   approximation the VPC/ICC of a \emph{log-link count} model uses, from table 1
#'   of Nakagawa, Johnson & Schielzeth (2017): \code{"lognormal"} (default,
#'   \eqn{\ln(1 + 1/\lambda\ [+\ 1/\theta])}), \code{"delta"}
#'   (\eqn{1/\lambda\ [+\ 1/\theta]}), or \code{"trigamma"}
#'   (\eqn{\psi_1(\lambda)} for Poisson,
#'   \eqn{\psi_1((1/\lambda + 1/\theta)^{-1})} for the negative binomial).
#'   Recorded on the fit, so \code{summary()}, the bootstrap intervals and the
#'   longitudinal VPC(t) all use the same one. Inert for every other family.
#'
#'   The three agree above a marginal count of about \eqn{\lambda = 2} and diverge
#'   sharply below it -- at \eqn{\lambda = 0.34} they give level-1 variances of
#'   1.37, 2.94 and 9.76, so the VPC moves by a factor of six -- which is why
#'   \code{summary()} reports the method and the \eqn{\lambda} it was evaluated at,
#'   and warns below the threshold. The default is the same \emph{approximation}
#'   \code{insight::get_variance()} and \code{performance::icc()} use by default,
#'   but the two are not numerically identical on an adjusted model: they evaluate
#'   it at a different \eqn{\lambda}. MAIHDA averages the row-wise marginal
#'   expected counts \eqn{\exp(x_i'\beta + v_i/2)} over the analytic sample --
#'   the global-\eqn{\lambda} form of Nakagawa et al. -- while \code{insight}
#'   plugs in a single \eqn{\exp(\beta_0 + \sigma^2/2)} taken from an
#'   intercept-only null model. On a \emph{null} model, the MAIHDA headline, the
#'   two coincide exactly. On an adjusted model Jensen's inequality makes the
#'   MAIHDA \eqn{\lambda} the larger, so its level-1 variance is the smaller and
#'   its VPC slightly the larger; the gap grows with the spread of the fitted
#'   means, from well under 1\% of the level-1 variance for a weak covariate to
#'   tens of per cent for a very strong one. Nakagawa et al. themselves
#'   recommend the trigamma form; note that it is the least conservative at low
#'   counts (it drives the VPC toward zero), because it is the variance of
#'   \eqn{\log X} for \eqn{X \sim \mathrm{Gamma}(\lambda, 1)} and that
#'   approximation is weakest exactly where most counts are zero.
#' @param ... Additional arguments passed to \code{lmer}/\code{glmer} (lme4),
#'   \code{brm} (brms), or \code{WeMix::mix()} (wemix; e.g. \code{nQuad},
#'   \code{fast}). The lme4-style \code{weights} (precision weights),
#'   \code{subset}, and \code{offset} arguments are honoured only by the
#'   \code{lme4} engine, which applies them directly. The \code{wemix},
#'   \code{ordinal}, and \code{brms} engines reject them: none takes them as a
#'   top-level fitting argument (brms in particular expects weighting/offset as
#'   formula addition terms, \code{weights(.)} / \code{offset(.)}, and design
#'   weights via \code{sampling_weights}). Prefilter \code{data} instead of using
#'   \code{subset} on those engines.
#'
#' @return A maihda_model object containing:
#'   \item{model}{The fitted model object (lme4, brms, WeMix, or ordinal::clmm)}
#'   \item{engine}{The engine used ("lme4", "brms", "wemix", or "ordinal")}
#'   \item{sampling_weights}{The sampling-weight column name when supplied,
#'     NULL otherwise}
#'   \item{formula}{The model formula}
#'   \item{data}{The data used for fitting}
#'   \item{family}{The family used}
#'   \item{strata_info}{The strata information from make_strata() if available, NULL otherwise}
#'   \item{context_vars}{The context variable name(s) when \code{context} was
#'     supplied, NULL otherwise}
#'   \item{interactions}{The \code{maihda_interactions} diagnostic when
#'     \code{interactions} is not \code{FALSE}, NULL otherwise}
#'   \item{response_recoding}{For a recoded two-level outcome, a data frame mapping
#'     each original level to its 0/1 value and role (reference/event); NULL when no
#'     recoding occurred}
#'   \item{diagnostics}{Fit-quality diagnostics, surfaced by the print and
#'     summary methods: singular fit / convergence for lme4 and WeMix, MCMC
#'     convergence (maximum Rhat, divergent transitions) for brms, and the
#'     optimizer convergence code for an ordinal (clmm) fit. An lme4 fit also
#'     carries likelihood-adequacy caveats -- count overdispersion and zero
#'     inflation, stratum random-effect non-normality and longitudinal residual
#'     autocorrelation -- reported only when a conservative threshold is crossed,
#'     since the VPC/PCV and interaction estimates are conditional on the
#'     likelihood holding. Every one of those checks is lme4-only, so a cumulative
#'     (clmm) fit raises no adequacy caveat at all: it stores a single descriptive
#'     fixed-only proportional-odds proxy statistic, never flagged because it
#'     cannot be separated from stratum heterogeneity; use
#'     \code{\link{maihda_proportional_odds_test}} to test that assumption}
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
                       id = NULL, time = NULL, time_degree = 1,
                       stratum_slope = TRUE, interactions = FALSE,
                       count_approximation = c("lognormal", "delta", "trigamma"),
                       ...) {
  # Input validation
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula object")
  }

  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }

  # Validated up front (not lazily at summary() time) so a typo fails at the fit
  # rather than several minutes of MCMC later. Inert for every family whose
  # level-1 variance is a distribution constant; only the log-link count branches
  # read it.
  count_approximation <- maihda_check_count_approximation(count_approximation)

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
  dot_quos <- rlang::enquos(...)
  # Internal fast-path flag set only by maihda()'s crossed-dimensions preliminary
  # pass, which needs the resolved strata/family metadata but discards the fit
  # itself. Pull it out of the engine dots BEFORE they are evaluated, so it is
  # neither eval'd against the data nor forwarded to lme4/brms; a full fit never
  # sees it (default FALSE). It is intentionally not a formal argument -- it is an
  # internal switch, not part of the public fit_maihda() interface.
  metadata_only <- isTRUE(tryCatch(
    rlang::eval_tidy(dot_quos[[".metadata_only"]]), error = function(e) FALSE))
  dot_quos[[".metadata_only"]] <- NULL
  dot_vals <- lapply(dot_quos, function(q) rlang::eval_tidy(q, data = data))
  subset_value <- dot_vals[["subset"]]
  # Resolve a character (row-name) subset to a positional logical mask against this
  # data's own row names, up front, so the engine receives a logical vector rather
  # than a raw character subset. base `[`'s character row-indexing partial-matches
  # and inserts phantom NA rows, so a raw character subset can make the fitted sample
  # diverge from the analytic sample family detection, 0/1 recoding and strata
  # binning key off (all of which resolve names by exact %in% via maihda_row_mask).
  # A no-op for logical/numeric subsets and when there is no `subset`.
  if (is.character(subset_value)) {
    subset_value <- rownames(data) %in% subset_value
    dot_vals[["subset"]] <- subset_value
  }
  weights_value <- dot_vals[["weights"]]
  # The external offset= (lme4's only NA-dropping fitting argument besides weights)
  # is forwarded to the engine below, where lme4 puts it in its model frame and
  # na.omit drops offset-NA rows. Extract it here so every analytic mask -- family
  # detection, 0/1 recoding, strata auto-binning, longitudinal validation/centering
  # -- is keyed off the SAME rows the engine fits: an offset-NA row carrying an
  # out-of-sample outcome value would otherwise flip the auto-detected family, and an
  # out-of-range occasion would mis-anchor the growth time-centering. Only the lme4
  # engine accepts an offset (wemix/brms/ordinal reject it above/below), so this is
  # NULL on every other path.
  offset_value <- dot_vals[["offset"]]

  if (!is.null(sampling_weights) && "weights" %in% names(dot_vals)) {
    stop("Supply either 'sampling_weights' (design weights) or 'weights' ",
         "(precision weights), not both.", call. = FALSE)
  }

  # Normalize invalid lme4 PRECISION weights (zero / negative / non-finite -> NA,
  # the case lme4 actually drops) up front, so binary/ordinal detection, strata
  # auto-binning, and the fit all use the same analytic sample. Shared with
  # compare_maihda_groups() and maihda() via maihda_normalize_precision_weights()
  # so the three cannot drift apart. A no-op on the wemix/brms paths (they reject
  # precision weights below) and when no numeric `weights` was supplied.
  if (is.numeric(weights_value)) {
    weights_value <- maihda_normalize_precision_weights(
      weights_value, nrow(data), engine, "fit_maihda()")
    dot_vals[["weights"]] <- weights_value
  }

  if (identical(engine, "wemix")) {
    unsupported_dots <- intersect(c("weights", "subset", "offset"), names(dot_vals))
    if (length(unsupported_dots) > 0) {
      stop("Argument(s) not supported by engine = \"wemix\": ",
           paste(unsupported_dots, collapse = ", "),
           ". Subset or transform the data before fitting.", call. = FALSE)
    }
  }
  if (identical(engine, "brms")) {
    # brms does not take lme4's data-masked fitting arguments as top-level
    # arguments: it has no `subset`, and expects weighting / offset information as
    # formula ADDITION terms (weights(), offset()), not as `weights=` / `offset=`.
    # brms::brm() would silently absorb them into `...` and ignore them -- while
    # family detection and strata auto-binning above DID honour them -- so
    # preprocessing would describe one analytic sample and brms fit another
    # (silently changing coefficients, variance components, VPC, and PCV). Reject
    # them with guidance rather than fit the wrong model. (Design weights come
    # through 'sampling_weights', which brms supports as likelihood weights.)
    unsupported_dots <- intersect(c("weights", "subset", "offset"), names(dot_vals))
    if (length(unsupported_dots) > 0) {
      stop("Argument(s) not supported by engine = \"brms\": ",
           paste(unsupported_dots, collapse = ", "),
           ". brms takes weighting and offset information as formula addition ",
           "terms -- put offset(.) or weights(.) in the model formula, pass ",
           "design weights via 'sampling_weights', and prefilter 'data' instead ",
           "of using 'subset'.", call. = FALSE)
    }
  }

  # The weighted engines drop rows whose sampling weight is non-finite or <= 0
  # (not just NA), so base the binary/ordinal detection AND any 0/1 recoding on the
  # same analytic sample: maihda_sampling_weight_mask() maps those engine-excluded
  # weights to NA, which the shared row mask (maihda_row_mask) already drops.
  # Precision (lme4 prior) weights keep their NA-only handling.
  detect_weights <- if (!is.null(sampling_weights)) {
    maihda_sampling_weight_mask(data[[sampling_weights]])
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
                                weights = detect_weights, offset = offset_value),
      error = function(e) FALSE)
    if (isTRUE(is_binary)) {
      warning("The outcome variable appears to be binary. Automatically switching to family = 'binomial'. To fit a Linear Probability Model, explicitly specify family = 'gaussian'.", call. = FALSE)
      family <- "binomial"
    } else if (isTRUE(tryCatch(
      maihda_response_is_ordinal(formula, data, subset = subset_value,
                                 weights = detect_weights, offset = offset_value),
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

  # Longitudinal (3-level growth) MAIHDA: when 'time' is supplied, validate the
  # id/time specification now (engine and the wemix/ordinal/context/sampling
  # restrictions are already resolved). The growth formula is built AFTER strata
  # resolution below, so this only records the validated spec.
  lng_spec <- NULL
  if (!is.null(time) || !is.null(id)) {
    lng_spec <- maihda_validate_longitudinal(id, time, time_degree, data,
                                             engine = engine,
                                             sampling_weights = sampling_weights,
                                             context = context,
                                             formula = formula,
                                             stratum_slope = stratum_slope)
  }

  # Parse formula to find grouping variables and resolve the strata shorthand.
  # Shared with maihda_describe(), so the pre-model description is built from
  # exactly the same parsing/strata machinery as the fit.
  # Auto-bin the strata cut-points on EXACTLY the rows the fit uses, so
  # fit_maihda(data, subset = keep, ...) gives the same tertile boundaries -- and
  # therefore the same stratum membership and VPC/PCV -- as fitting the already
  # filtered analytic data. The analytic sample drops rows removed by `subset`, by a
  # missing outcome or covariate (na.omit over the model frame), by a missing
  # precision weight, OR by an invalid sampling weight (non-finite / <= 0, mapped to
  # NA in `detect_weights` above). An earlier version masked on `subset` alone, on
  # the assumption that the other exclusions are na.rm'd identically whether or not
  # the data is pre-filtered -- they are NOT. A row dropped from the fit for one of
  # those reasons can still carry an extreme but non-missing value of the binning
  # variable, which then shifts the quantiles and silently redefines OTHER rows'
  # strata. maihda_analytic_keep_mask() reproduces model.frame()'s row selection (so
  # transformed terms such as log(x) are handled too); it returns NULL if the frame
  # cannot be built, in which case fall back to the subset-only mask.
  strata_bin_rows <- tryCatch(
    maihda_analytic_keep_mask(formula, data, subset = subset_value,
                              weights = detect_weights, offset = offset_value),
    error = function(e) NULL)
  if (is.null(strata_bin_rows) && (!is.null(subset_value) || !is.null(offset_value))) {
    strata_bin_rows <- maihda_row_mask(data, subset = subset_value,
                                       offset = offset_value)
  }
  strata_res <- maihda_resolve_strata_formula(formula, data, autobin,
                                              bin_rows = strata_bin_rows)
  formula <- strata_res$formula
  data <- strata_res$data
  strata_info <- strata_res$strata_info
  strata_vars <- strata_res$strata_vars
  strata_sep <- strata_res$strata_sep
  strata_autobin_info <- strata_res$autobin_info

  # Longitudinal growth structure: now that the stratum grouping is resolved,
  # replace the random part with the canonical 3-level growth blocks
  # (time... | id) + (time... | stratum) and ensure the time polynomial is in the
  # fixed part. The fit then flows through the unchanged lme4/brms branches (they
  # already pass random slopes to the engine); $longitudinal_info tags the model so
  # summary()/predict()/plot() route to the time-varying path.
  longitudinal_info <- NULL
  if (!is.null(lng_spec)) {
    has_stratum_re <- any(vapply(reformulas::findbars(formula),
      function(b) "stratum" %in% all.vars(b[[3]]), logical(1)))
    if (!has_stratum_re) {
      stop("A longitudinal MAIHDA needs a stratum random effect. Use the shorthand ",
           "(1 | var1:var2) or include (1 | stratum); the id/time growth slopes are ",
           "added automatically (do not write them in the formula).", call. = FALSE)
    }
    # The analytic sample -- exactly the rows the engine fits, after `subset`, a
    # missing precision weight, and missing-value dropping over the resolved
    # formula's variables plus id/time. Every ROW-SENSITIVE longitudinal check and
    # the time-centering below key off THIS mask, not the full input: validating or
    # centering on rows the fit drops is what let an excluded row (an out-of-subset
    # occasion, a missing outcome/covariate, a zero/NA-weight row) spuriously fail
    # the cross-stratum id check, slip past the repeated-measures gate, or mis-anchor
    # the centering. Falls back to keeping every row if the mask cannot be built.
    analytic_keep <- tryCatch({
      # Transformation-aware: model.frame + na.omit over the resolved formula, so a
      # row whose fixed-effect transformation is non-finite (e.g. log(x) of x <= 0)
      # is dropped here exactly as lme4 will drop it. A raw complete.cases() over the
      # source columns instead keeps such a row (the raw x is non-NA), leaving the
      # longitudinal checks and the time-centering below keyed off rows the fit never
      # sees -- mis-anchoring the centre and spuriously failing the cross-stratum id
      # check on an excluded (id, stratum) pairing.
      k <- maihda_analytic_keep_mask(formula, data, subset = subset_value,
                                     weights = weights_value, offset = offset_value)
      if (is.null(k)) {
        # Model frame could not be built: fall back to the raw check over the
        # formula's source variables (id/time are appended below either way).
        k <- maihda_row_mask(data, subset = subset_value, weights = weights_value,
                             offset = offset_value)
        fvars <- intersect(all.vars(formula), names(data))
        if (length(fvars) > 0) {
          k <- k & stats::complete.cases(data[, fvars, drop = FALSE])
        }
      }
      # id/time are not part of the (pre-growth) formula yet, so guard them
      # explicitly. They are plain grouping/time columns with no transformation, for
      # which a raw complete.cases is exact.
      idtime <- intersect(c(lng_spec$id, lng_spec$time), names(data))
      if (length(idtime) > 0) {
        k <- k & stats::complete.cases(data[, idtime, drop = FALSE])
      }
      k
    }, error = function(e) rep(TRUE, nrow(data)))
    if (!any(analytic_keep)) analytic_keep <- rep(TRUE, nrow(data))
    analytic_data <- data[analytic_keep, , drop = FALSE]

    # Genuinely repeated measures, re-checked on the analytic sample: at least one id
    # must recur among the FITTED rows, else the (time | id) person effects are
    # unidentified and this is not a longitudinal design. maihda_validate_longitudinal()
    # already applied this to the full input; repeating it here catches data whose only
    # repeats sit on rows the fit drops (an out-of-subset/missing/zero-weight occasion).
    lng_ids <- analytic_data[[lng_spec$id]]
    if (!any(duplicated(lng_ids[!is.na(lng_ids)]))) {
      stop("The data do not look longitudinal: every '", lng_spec$id, "' value is ",
           "unique, so there are no repeated measurements to model. Supply ",
           "long-format data (one row per measurement occasion).", call. = FALSE)
    }
    # Ids must identify people globally, not within a site/group: an id spanning more
    # than one stratum would merge different people's trajectories in (time | id).
    # Checked on the analytic sample so an excluded row cannot inject a spurious
    # (id, stratum) pairing the fit never sees.
    maihda_check_longitudinal_ids(analytic_data, lng_spec$id)
    # Growth-term identifiability, re-checked on the analytic sample: the input
    # passed maihda_validate_longitudinal(), but the rows dropped above can push
    # the FITTED rows below the thresholds -- fewer than time_degree + 1 distinct
    # times, or no person left with two distinct times -- leaving lme4 to return
    # arbitrary slope variances with only a singularity note.
    maihda_check_longitudinal_times(analytic_data, lng_spec$id, lng_spec$time,
                                    lng_spec$time_degree)

    # Fit the growth terms on internally CENTERED time whenever the time axis does
    # not start at 0 (age, calendar year, waves coded 10, 11, ...): the raw
    # polynomial basis over an offset range is near-collinear, and lme4 can converge
    # to a false optimum WITHOUT flagging a convergence failure, silently corrupting
    # the time-varying VPC and the PCV. Centering makes the fit the same optimization
    # problem as the equivalent 0-anchored coding. Every user-facing time (ref_time,
    # reporting grids, plots) stays on the original scale, and prediction newdata
    # rebuilds the derived column from the original time column (see
    # maihda_prepare_prediction_data). Center on the analytic rows (analytic_keep):
    # deriving the centre from excluded early/late waves would leave the fitted times
    # off-centre and defeat the stability protection.
    tv <- data[[lng_spec$time]]
    center <- maihda_longitudinal_center(tv[analytic_keep])
    time_term <- lng_spec$time
    if (center != 0) {
      # A formula already referencing the derived column is a package-derived
      # refit (maihda()'s null/adjusted pair re-entering with the first fit's
      # original_data): recomputing the column from the original time column is
      # idempotent and needs no repeat message. A fresh user column of the
      # reserved name was rejected by maihda_validate_longitudinal() above.
      derived_refit <- .maihda_ctime_col %in% all.vars(formula)
      data[[.maihda_ctime_col]] <- tv - center
      time_term <- .maihda_ctime_col
      if (!derived_refit) {
        message("fit_maihda(): the time axis '", lng_spec$time, "' starts at ",
                center, ", not 0, so the growth terms are fit on internally ",
                "centered time (", lng_spec$time, " - ", center, ") for ",
                "numerical stability. All results are reported on the original '",
                lng_spec$time, "' scale.")
      }
    }
    formula <- maihda_longitudinal_formula(formula, lng_spec$id, time_term,
                                           lng_spec$time_degree,
                                           orig_time = lng_spec$time,
                                           stratum_slope = lng_spec$stratum_slope)
    # time_range/ref_time here are provisional (pre-fit): they are recomputed
    # from the fitted analytic frame after the engine drops rows (see below),
    # so a baseline wave lost to missing outcomes does not anchor the VPC/PCV.
    longitudinal_info <- list(id = lng_spec$id, time = lng_spec$time,
                              time_degree = lng_spec$time_degree,
                              stratum_slope = lng_spec$stratum_slope,
                              time_term = time_term, time_center = center,
                              time_range = range(tv, na.rm = TRUE),
                              ref_time = min(tv, na.rm = TRUE))
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
    formula <- maihda_apply_context_formula(formula, context, strata_vars)
    context_info <- list(context_vars = context)
  }

  # Convert family to family object if it's a string or constructor function
  # (shared with maihda_describe(), so both resolve a family spec identically).
  family <- maihda_resolve_family_spec(family)

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
                              weights = detect_weights, offset = offset_value)
  response_recoding <- NULL
  if (response_is_binary) {
    data <- maihda_prepare_binomial_response(data, formula, subset = subset_value,
                                             weights = detect_weights,
                                             offset = offset_value)
    # Mapping of original outcome levels to 0/1 (which level is the modeled event),
    # captured so it is inspectable on the returned model object.
    response_recoding <- attr(data, "response_recoding")
  }

  # The pre-fit descriptive frame stored as $original_data: the FULL data (after
  # strata resolution and any 0/1 recoding) BEFORE the engine's analytic-sample drop,
  # so maihda_describe() can report the original total, missingness, and
  # excluded-weight counts against the analytic $data. Left NULL here and defaulted to
  # the path's `data` just before the result is built; only the design-weighted brms
  # path (which prunes `data` in place below) sets it explicitly to the pre-prune frame.
  # The wemix/ordinal paths intentionally leave it NULL -> it defaults to their
  # pre-built analytic sample, matching maihda_describe()'s documented contract
  # (totals == analytic there).
  original_data <- NULL

  # Metadata-only fast path (maihda()'s crossed-dimensions preliminary pass). That
  # decomposition builds and fits a SINGLE cross-classified model from the resolved
  # strata/family metadata; the preliminary fit of the supplied formula was only
  # ever read for that metadata and then discarded, so fitting it doubled the cost
  # (worst for brms: an extra compile + MCMC). Return the metadata -- resolved above
  # by the SAME strata / family / response-recoding machinery a full fit uses, so it
  # cannot drift from the fitted path -- and skip the engine dispatch entirely.
  # $original_data mirrors the default applied at the result-building tail (the full
  # pre-fit frame, after strata resolution and any 0/1 recoding). No class of
  # "maihda_model" so a stray summary()/predict() on it errors instead of silently
  # treating an unfitted stub as a model.
  if (isTRUE(metadata_only)) {
    return(structure(
      list(
        formula = formula,
        family = family,
        count_approximation = count_approximation,
        strata_info = strata_info,
        strata_vars = strata_vars,
        strata_sep = strata_sep,
        strata_autobin_info = strata_autobin_info,
        original_data = data,
        context_vars = context,
        context_info = context_info,
        longitudinal_info = longitudinal_info,
        metadata_only = TRUE
      ),
      class = "maihda_model_metadata"
    ))
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
    # PSEUDO-posterior: point estimates target the population-weighted estimand,
    # intervals are not design-based. The weights are normalized to mean 1 so
    # expansion weights do not inflate the effective sample size, and rows with
    # missing/non-positive weights are dropped here so the stored data matches the
    # fitted rows.
    if (!is.null(sampling_weights)) {
      # Keep the pre-prune frame as $original_data (the descriptive total). The drop
      # below is the brms engine's analytic-sample selection -- the pseudo-likelihood
      # counterpart of lme4's na.omit -- so the fitted $data reflects it, but
      # $original_data must stay the full pre-fit frame (as on the lme4 path) or
      # maihda_describe() cannot report the original total / missingness / excluded
      # invalid-weight rows and would make total == analytic. A derived null/adjusted
      # refit re-enters here and re-prunes this same frame, so the fitted models are
      # unchanged; only the recorded provenance differs.
      original_data <- data
      prep <- maihda_prepare_brms_sampling_weights(data, formula, sampling_weights)
      data <- prep$data
      formula <- prep$formula
      fit_env$data <- data
      # No row-aligned `...` need re-slicing after this drop: the brms engine
      # rejects the lme4-style weights/subset/offset arguments (the only row-length
      # forwarded values) up front, so every remaining dot is a scalar / non-row
      # object brms takes as-is.
      message("fit_maihda(): sampling weights enter the brms model as likelihood ",
              "weights (normalized to mean 1), giving a pseudo-posterior: point ",
              "estimates target the population-weighted estimand, but credible ",
              "intervals are not design-based -- interpret them cautiously.")
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
      # Re-check the category count on the analytic sample the same way the clmm
      # path does, but BEFORE the (expensive) Stan fit: a category present only on
      # rows brms drops for missingness would otherwise silently reduce the model
      # order. sampling_weights, if any, already pruned `data` above, so the model
      # frame here reflects the rows brms will actually use.
      resp_name <- all.vars(formula[[2]])[1]
      ord_keep <- maihda_analytic_keep_mask(formula, data)
      ord_y <- if (is.null(ord_keep)) data[[resp_name]] else data[[resp_name]][ord_keep]
      maihda_ordinal_assert_min_levels(droplevels(ord_y), resp_name)
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

  # Capture fit-quality diagnostics (singular fit / non-convergence) and likelihood-
  # adequacy caveats so they can be reported by print()/summary(); lme4 surfaces the
  # convergence warnings only once at fit time. longitudinal_info (id / time) lets the
  # adequacy pass compute within-unit residual autocorrelation.
  diagnostics <- maihda_fit_diagnostics(model, longitudinal_info)

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

  # Recompute the longitudinal reference time and time range from the FITTED
  # analytic frame, not the pre-fit data: lme4/brms drop rows with missing
  # outcomes (or predictors), which can remove an entire baseline wave. The
  # earlier capture (from `data`) recorded ref_time = min(time) over rows that
  # may not all survive the fit; the VPC/PCV summaries then report a baseline at
  # a time not represented in the fitted sample (an extrapolation). model_data
  # holds the rows the engine actually used, so its min/range are the analytic
  # baseline and span.
  if (!is.null(longitudinal_info)) {
    # A centered fit's model frame carries only the derived centered column;
    # re-attach the original time (centered + offset) so every consumer reading
    # object$data[[time]] -- the reporting grids, the ref_time recomputation
    # below, the trajectory plots -- stays on the user's original time scale.
    tt <- longitudinal_info$time_term
    if (!identical(tt, longitudinal_info$time) &&
        !longitudinal_info$time %in% names(model_data) &&
        tt %in% names(model_data)) {
      model_data[[longitudinal_info$time]] <-
        model_data[[tt]] + longitudinal_info$time_center
    }
    tv_fit <- model_data[[longitudinal_info$time]]
    longitudinal_info$time_range <- range(tv_fit, na.rm = TRUE)
    longitudinal_info$ref_time <- min(tv_fit, na.rm = TRUE)
  }

  # Default the descriptive frame to this path's `data` unless the design-weighted
  # brms path already captured the pre-prune frame above. For lme4 / unweighted brms
  # this is the full pre-fit frame; for wemix / ordinal it is the pre-built analytic
  # sample (totals == analytic, as documented).
  if (is.null(original_data)) original_data <- data

  result <- structure(
    list(
      model = model,
      engine = engine,
      formula = formula,
      data = model_data,
      original_data = original_data,
      family = family,
      count_approximation = count_approximation,
      strata_info = strata_info,
      strata_vars = strata_vars,
      strata_sep = strata_sep,
      strata_autobin_info = strata_autobin_info,
      context_vars = context,
      context_info = context_info,
      sampling_weights = sampling_weights,
      longitudinal_info = longitudinal_info,
      response_recoding = response_recoding,
      diagnostics = diagnostics
    ),
    class = "maihda_model"
  )

  # Opt-in per-stratum interaction diagnostic (parallel to maihda(), which computes
  # it by default). For a single fit this is meaningful only on the *adjusted*
  # model -- maihda_interactions() warns if the formula looks like a null model.
  if (!isFALSE(interactions)) {
    result <- maihda_attach_interactions(result, interactions)
  }

  return(result)
}

# Resolve the intersectional-strata shorthand in a MAIHDA formula. Automatic
# strata creation is only safe for the documented shorthand: one intercept-only
# non-stratum grouping term such as (1 | gender:race); more complex
# random-effect structures must be specified explicitly after calling
# make_strata(). When the shorthand is present, make_strata() builds the strata,
# `data` gains the stratum column (plus the strata_* attributes), and the
# formula's grouping term is rewritten to (1 | stratum); a formula that already
# groups on `stratum` passes through, picking up any strata_* attributes that
# make_strata() attached to `data`. Shared by fit_maihda() and
# maihda_describe() so the pre-model description and the fit build their strata
# from the same machinery -- identical IDs, labels, counts, and validation.
# TRUE when a random-effect term's left-hand side carries a random INTERCEPT --
# whether written explicitly (`1`, `1 + x`, `0 + 1`) or implicitly (`x`, which lme4
# expands to `1 + x`). This is TRUE for any term that CONTAINS an intercept (not only a
# bare intercept), so it catches a compound intercept-plus-slope term. maihda_resolve_strata_formula()
# uses it to reject more than one intercept-bearing term on 'stratum': (1 | stratum) +
# (1 + x | stratum), and (1 | stratum) + (x | stratum), each put two intercepts on the
# same grouping factor, which lme4 splits arbitrarily across 'stratum' and 'stratum.1'
# (non-identifiable). A slope-only second term (0 + x | stratum) carries no intercept and
# stays identifiable. Reads the formula's intercept attribute -- exactly how lme4 itself
# decides whether a term has a random intercept -- so it needs no data and cannot disagree
# with the fitted structure. Returns FALSE on any parse error (conservative: no over-reject).
maihda_re_lhs_has_intercept <- function(lhs) {
  fm <- tryCatch(stats::as.formula(paste("~", paste(deparse(lhs), collapse = " "))),
                 error = function(e) NULL)
  if (is.null(fm)) return(FALSE)
  int <- tryCatch(attr(stats::terms(fm), "intercept"), error = function(e) 0L)
  isTRUE(int == 1L)
}

# TRUE when a random-effect term's design matrix -- built on the analytic `data` the
# way lme4 builds it, model.matrix(~ lhs, data) -- carries the all-ones (intercept)
# vector in its COLUMN SPACE. This catches an intercept that maihda_re_lhs_has_intercept()
# (which reads only the formula intercept attribute) misses because it is realized
# NUMERICALLY rather than written as `1`: a constant slope column such as (0 + one |
# stratum) with one == 1, or a full dummy encoding (0 + f | stratum) whose columns sum
# to the intercept. Two stratum terms that each span the intercept put two intercepts on
# the same grouping factor, which lme4 fits as separate 'stratum'/'stratum.1' variance
# components and splits the between-stratum variance arbitrarily across -- the same
# non-identifiable defect the attribute test guards against for explicit duplicates.
# Membership is tested by QR rank: the all-ones vector is in colspace(Z) exactly when
# rank([Z, 1]) == rank(Z). Falls back to the attribute test when the design cannot be
# built or ranked (an absent/all-NA column, a parse error): conservative -- it never
# silently drops the guard, and never over-rejects on an un-buildable term.
maihda_re_lhs_spans_intercept <- function(lhs, data) {
  fm <- tryCatch(stats::as.formula(paste("~", paste(deparse(lhs), collapse = " "))),
                 error = function(e) NULL)
  if (!is.null(fm) && is.data.frame(data) && nrow(data) > 0L) {
    z <- tryCatch(
      stats::model.matrix(fm, stats::model.frame(fm, data = data,
                                                 na.action = stats::na.omit)),
      error = function(e) NULL)
    if (!is.null(z) && nrow(z) > 0L && ncol(z) > 0L) {
      r_z <- tryCatch(qr(z)$rank, error = function(e) NA_integer_)
      r_aug <- tryCatch(qr(cbind(z, 1))$rank, error = function(e) NA_integer_)
      if (!is.na(r_z) && !is.na(r_aug)) {
        return(isTRUE(r_aug == r_z))
      }
    }
  }
  maihda_re_lhs_has_intercept(lhs)
}

maihda_resolve_strata_formula <- function(formula, data, autobin = TRUE,
                                          bin_rows = NULL) {
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
    is_stratum_term <- vapply(grouping_vars_by_term, function(vars) {
      identical(vars, "stratum")
    }, logical(1))
    has_stratum_group <- any(is_stratum_term)

    # Reject more than one random-effect term that puts an INTERCEPT on 'stratum'.
    # lme4 fits `... + (1 | stratum) + (1 | stratum)` -- and equally the compound
    # `(1 | stratum) + (1 + x | stratum)`, or `(1 | stratum) + (x | stratum)` whose
    # second intercept is implicit -- as two SEPARATE stratum variance components
    # ("stratum" and "stratum.1") and splits the between-stratum variance arbitrarily
    # between them. That partition is non-identifiable and makes the VPC and PCV
    # ill-defined: the summary counts only the first component as between-stratum
    # variance and misclassifies the rest as "other random effects". Test each stratum
    # term for an intercept by whether its random-effect design matrix -- built on the
    # ANALYTIC rows lme4 fits (bin_rows) -- spans the intercept vector
    # (maihda_re_lhs_spans_intercept), not merely by its formula intercept attribute.
    # That catches every spelling -- 1, 0 + 1, 1 + x, and a bare slope x with its implicit
    # intercept -- AND an intercept realized NUMERICALLY: a constant slope column such as
    # (0 + one | stratum) with one == 1, or a full dummy encoding whose columns sum to 1,
    # both of which the attribute test reports as "no intercept" yet lme4 fits as a second
    # stratum intercept variance. A single compound (1 + x | stratum) term (one correlated
    # intercept+slope) and an uncorrelated pair (1 | stratum) + (0 + x | stratum) (second
    # term a genuinely varying slope, no intercept in its column space) are both
    # identifiable and are left to the downstream intercept-only validation / the fit.
    if (sum(is_stratum_term) > 1) {
      # The analytic sample (subset / missing / weight-NA rows dropped) -- a column can be
      # constant there yet vary in the raw data, so the intercept-span test must see the
      # rows the fit uses. bin_rows is that keep mask (NULL when no filtering applies).
      analytic_data <- if (!is.null(bin_rows) && is.logical(bin_rows) &&
                           length(bin_rows) == nrow(data)) {
        data[bin_rows, , drop = FALSE]
      } else {
        data
      }
      has_intercept <- vapply(re_terms[is_stratum_term], function(term) {
        maihda_re_lhs_spans_intercept(term[[2]], analytic_data)
      }, logical(1))
      if (sum(has_intercept) > 1) {
        stop("The model formula includes ", sum(has_intercept), " random-effect terms ",
             "that each put an intercept on 'stratum' (e.g. (1 | stratum) and ",
             "(1 + x | stratum), or (x | stratum) whose intercept is implicit; ",
             "(0 + 1 | stratum) also builds an intercept column). Duplicate stratum ",
             "intercepts are non-identifiable: lme4 splits the between-stratum variance ",
             "arbitrarily across them (fitted as 'stratum' and 'stratum.1'), so the VPC ",
             "and PCV are not well defined. Use a single stratum intercept -- combine ",
             "them into one (1 + x | stratum) term, or make the extra term slope-only, ",
             "(0 + x | stratum).", call. = FALSE)
      }
    }

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

      strata_result <- make_strata(data, vars = strata_vars, autobin = autobin,
                                   bin_rows = bin_rows)
      data$stratum <- strata_result$data$stratum
      strata_info <- strata_result$strata_info
      strata_sep <- strata_result$sep
      strata_autobin_info <- strata_result$autobin_info
      attr(data, "strata_info") <- strata_info
      attr(data, "strata_vars") <- strata_vars
      attr(data, "strata_sep") <- strata_sep
      attr(data, "strata_autobin_info") <- strata_autobin_info

      fixed_formula <- maihda_nobars(formula)
      formula <- stats::update(fixed_formula, . ~ . + (1 | stratum))
    }
  }

  list(formula = formula, data = data, strata_info = strata_info,
       strata_vars = strata_vars, strata_sep = strata_sep,
       autobin_info = strata_autobin_info)
}

# Append the contextual cross-classified random intercept(s) to a formula whose
# stratum grouping is already resolved, validating the roles first: the formula
# must carry a stratum random effect, and a context variable may not double as a
# stratum dimension or as a fixed-effect term (either would absorb the variance
# the contextual partition is meant to estimate). Idempotent: a context that is
# already a random-effect grouping (e.g. a maihda() refit of a derived formula
# that carries the context term) is validated but not appended again. Shared by
# fit_maihda() and maihda_describe().
maihda_apply_context_formula <- function(formula, context, strata_vars) {
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
  clash_fixed <- intersect(context, all.vars(maihda_nobars(formula)[[3]]))
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
  formula
}

# Resolve a user family specification (name string, family/constructor function,
# or family-like list) to the family object the engines and summaries consume.
# "negbinomial" (the brms spelling) resolves to a plain marker list rather than
# a stats family object: there is no theta-free negative-binomial family
# constructor in stats -- lme4 estimates theta itself via glmer.nb() and brms
# via its 'shape' parameter, so no theta is needed (or wanted) here. Shared by
# fit_maihda() and maihda_describe().
maihda_resolve_family_spec <- function(family) {
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
  family
}

#' Print method for maihda_model
#'
#' @param x A maihda_model object
#' @param ... Additional arguments
#' @return No return value, called for side effects.
#' @export
print.maihda_model <- function(x, ...) {
  cat(maihda_palette()$bold("MAIHDA Model"), "\n", sep = "")
  cat("============\n\n")
  cat("Engine:", x$engine, "\n")
  cat("Family:", x$family$family, "\n")
  cat("Formula:", deparse(x$formula), "\n")
  if (!is.null(x$context_vars)) {
    cat("Context:", paste(x$context_vars, collapse = ", "),
        "(crossed contextual random intercept)\n")
  }
  if (!is.null(x$longitudinal_info)) {
    lng <- x$longitudinal_info
    ct <- maihda_lng_time_center(lng)
    cat(sprintf("Longitudinal: id = %s, time = %s%s, degree = %d (3-level growth)\n",
                lng$id, lng$time,
                if (ct != 0) sprintf(" (growth terms internally centered at %g)", ct) else "",
                lng$time_degree))
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
  maihda_print_interactions_line(x$interactions)
  cat("Underlying model:\n")
  print(x$model, ...)
  invisible(x)
}
