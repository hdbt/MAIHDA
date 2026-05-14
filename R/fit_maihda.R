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
#' @param family Character string or family object specifying the model family.
#'   Common options: "gaussian", "binomial", "poisson". Default is "gaussian".
#'   If the outcome variable appears to be binary (0/1) and the default family is used,
#'   the function will automatically switch to "binomial" and issue a warning.
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
fit_maihda <- function(formula, data, engine = "lme4", family = "gaussian", ...) {
  # Input validation
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula object")
  }

  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }

  engine <- match.arg(engine, c("lme4", "brms"))

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

  # Parse formula to find grouping variables
  # Check if "stratum" is not already the grouping variable
  # reformulas::findbars/nobars extracts the random effect terms (preferred)
  re_terms <- reformulas::findbars(formula)
  if (!is.null(re_terms)) {
    # Extract the names of all grouping variables
    grouping_vars <- unique(unlist(lapply(re_terms, function(x) all.vars(x[[3]]))))

    # If the grouping variables are not just "stratum", create strata
    if (length(grouping_vars) > 0 && !all(grouping_vars == "stratum")) {
      # Keep variables that exist in the data
      valid_vars <- intersect(grouping_vars, names(data))

      if (length(valid_vars) > 0) {
        strata_result <- make_strata(data, vars = valid_vars)
        data <- strata_result$data
        attr(data, "strata_info") <- strata_result$strata_info

        # Rewrite the formula to use (1 | stratum) instead of the original random effects
        # We need to drop all the original random effects and add (1 | stratum)
        fixed_formula <- reformulas::nobars(formula)
        formula <- stats::update(fixed_formula, . ~ . + (1 | stratum))
      }
    }
  }

  # Convert family to family object if it's a string
  if (is.character(family)) {
    family <- switch(family,
                     gaussian = gaussian(),
                     binomial = binomial(),
                     poisson = poisson(),
                     stop("Unsupported family: ", family))
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

    model <- brms::brm(formula, data = data, family = family, ...)
  }

  # Create maihda_model object
  # Capture strata_info if it exists as an attribute on the data
  strata_info <- attr(data, "strata_info")

  result <- structure(
    list(
      model = model,
      engine = engine,
      formula = formula,
      data = data,
      family = family,
      strata_info = strata_info
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
