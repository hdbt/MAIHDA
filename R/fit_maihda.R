#' Fit MAIHDA Model
#'
#' Fits a multilevel model for MAIHDA (Multilevel Analysis of Individual
#' Heterogeneity and Discriminatory Accuracy) using either lme4 or brms.
#'
#' @param formula A formula specifying the model. Should include random effect
#'   for stratum (e.g., \code{outcome ~ fixed_vars + (1 | stratum)}).
#' @param data A data frame containing the variables in the formula.
#' @param engine Character string specifying which engine to use: "lme4" (default)
#'   or "brms".
#' @param family Character string or family object specifying the model family.
#'   Common options: "gaussian", "binomial", "poisson". Default is "gaussian".
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
#' \dontrun{
#' # Create strata
#' strata_result <- make_strata(data, vars = c("gender", "race"))
#' 
#' # Fit model with lme4
#' model <- fit_maihda(outcome ~ age + (1 | stratum),
#'                     data = strata_result$data,
#'                     engine = "lme4")
#' 
#' # Fit model with brms (if installed)
#' model_brms <- fit_maihda(outcome ~ age + (1 | stratum),
#'                          data = strata_result$data,
#'                          engine = "brms")
#' }
#'
#' @export
#' @importFrom lme4 lmer glmer
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

#' Print method for maihda_model objects
#'
#' @param x A maihda_model object
#' @param ... Additional arguments passed to print method of underlying model
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
