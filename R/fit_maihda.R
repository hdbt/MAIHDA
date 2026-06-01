#' Fit MAIHDA Model
#'
#' Fits a multilevel model for MAIHDA (Multilevel Analysis of Individual
#' Heterogeneity and Discriminatory Accuracy) using either lme4 or brms.
#'
#' @param formula A formula specifying the model. Can include a random effect
#'   for stratum (e.g., \code{outcome ~ fixed_vars + (1 | stratum)}) or can
#'   directly specify the intersection variables to be used for forming strata
#'   (e.g., \code{outcome ~ fixed_vars + (1 | var1:var2:var3)}). If variables
#'   other than "stratum" are provided in the random effect, \code{\link{make_strata}}
#'   will be called internally to compute the strata and the formula will be
#'   updated.
#' @param data A data frame containing the variables in the formula.
#' @param engine Character string specifying which engine to use: "lme4" (default)
#'   or "brms".
#' @param family Character string, family object, or family function specifying
#'   the model family. Common options: "gaussian", "binomial", "poisson".
#'   Default is "gaussian".
#'   If the outcome variable appears to be binary and the default family is used,
#'   the function will automatically switch to "binomial", recode two-level
#'   responses to 0/1 for \code{glmer()}, and issue a warning.
#' @param autobin Logical indicating whether numeric variables used only for
#'   automatic strata creation should be binned by \code{\link{make_strata}}.
#'   Default is TRUE.
#' @param ... Additional arguments passed to \code{lmer}/\code{glmer} (lme4) or
#'   \code{brm} (brms).
#'
#' @return A maihda_model object containing:
#'   \item{model}{The fitted model object (lme4 or brms)}
#'   \item{engine}{The engine used ("lme4" or "brms")}
#'   \item{formula}{The model formula}
#'   \item{data}{The data used for fitting}
#'   \item{family}{The family used}
#'   \item{strata_info}{The strata information from make_strata() if available, NULL otherwise}
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
#' }
#'
#' @export
#' @importFrom lme4 lmer glmer
#' @importFrom reformulas findbars nobars
#' @importFrom stats gaussian binomial poisson
fit_maihda <- function(formula, data, engine = "lme4", family = "gaussian",
                       autobin = TRUE, ...) {
  # Input validation
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula object")
  }

  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }

  if (!is.character(engine) || length(engine) != 1 || !engine %in% c("lme4", "brms")) {
    stop("'engine' should be one of: lme4, brms", call. = FALSE)
  }

  # Automatically switch to binomial for binary outcomes if family is default
  if (missing(family)) {
    tryCatch({
      if (length(formula) == 3) {
        # Extract the response variable evaluated in the data context
        outcome_vals <- eval(formula[[2]], envir = data)
        outcome_vals <- stats::na.omit(outcome_vals)

        if (is.null(dim(outcome_vals)) && length(unique(outcome_vals)) == 2) {
          warning("The outcome variable appears to be binary. Automatically switching to family = 'binomial'. To fit a Linear Probability Model, explicitly specify family = 'gaussian'.", call. = FALSE)
          family <- "binomial"
        }
      }
    }, error = function(e) {
      # Silently proceed if formula extraction fails
    })
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

  # Convert family to family object if it's a string or constructor function
  if (is.character(family)) {
    family <- switch(family,
                     gaussian = gaussian(),
                     binomial = binomial(),
                     poisson = poisson(),
                     stop("Unsupported family: ", family))
  } else if (is.function(family)) {
    family <- family()
  }

  if (!is.list(family) || is.null(family$family) || is.null(family$link)) {
    stop("'family' must be a family name, family object, or family function.",
         call. = FALSE)
  }

  # Recode a two-level (Bernoulli) response to 0/1 (glmer and brms bernoulli both
  # accept 0/1). Aggregated binomial responses -- cbind(success, failure) or
  # `y | trials(n)` -- are left untouched and remain binomial() models.
  is_binomial_family <- family$family %in% c("binomial", "quasibinomial")
  response_is_binary <- is_binomial_family && maihda_response_is_binary(formula, data)
  if (response_is_binary) {
    data <- maihda_prepare_binomial_response(data, formula)
  }

  # Fit model based on engine
  if (engine == "lme4") {
    # Check if it's a Gaussian family (use lmer) or other (use glmer)
    if (family$family == "gaussian") {
      model <- lme4::lmer(formula, data = data, ...)
    } else {
      model <- lme4::glmer(formula, data = data, family = family, ...)
    }
  } else if (engine == "brms") {
    # Check if brms is installed
    if (!requireNamespace("brms", quietly = TRUE)) {
      stop("Package 'brms' is required but not installed. Please install it with: install.packages('brms')")
    }

    # brms models a 0/1 response with bernoulli(); passing binomial() would
    # require a trials specification and errors on Bernoulli data. Only rewrite
    # when the response really is a two-level vector -- aggregated binomial
    # (cbind / trials) must stay binomial().
    if (response_is_binary) {
      family <- brms::bernoulli(link = family$link)
    }

    model <- brms::brm(formula, data = data, family = family, ...)
  }

  # Store the actual analytic model frame so downstream calculations use the
  # same rows as lme4/brms after their NA handling.
  model_data <- maihda_model_frame(model, fallback = data)
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
      strata_autobin_info = strata_autobin_info
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
  cat("Formula:", deparse(x$formula), "\n\n")
  cat("Underlying model:\n")
  print(x$model, ...)
  invisible(x)
}
