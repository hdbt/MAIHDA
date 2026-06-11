#' Predict from MAIHDA Model
#'
#' Makes predictions from a fitted MAIHDA model, either at the stratum level
#' or individual level.
#'
#' @param object A maihda_model object from \code{fit_maihda()}.
#' @param newdata Optional data frame for making predictions. If NULL, uses the
#'   original data from model fitting.
#' @param type Character string specifying prediction type:
#'   \itemize{
#'     \item "individual": Individual-level predictions including random effects
#'     \item "strata": Stratum-level predictions (random effects only)
#'   }
#'   For backward compatibility, "link" or "response" may also be passed here
#'   and will be interpreted as individual-level predictions on that scale.
#' @param scale Character string specifying the prediction scale for
#'   individual-level predictions: "response" (default) or "link".
#' @param ... Additional arguments passed to predict method of underlying model.
#'
#' @return Depending on type:
#'   \itemize{
#'     \item For "individual": A numeric vector of predicted values on the
#'       requested scale
#'     \item For "strata": A data frame with stratum ID and predicted random
#'       effect. When \code{newdata} is supplied, the result is restricted to the
#'       strata present in \code{newdata} (and a stratum the model never saw is an
#'       error, as for "individual"); when \code{newdata} is \code{NULL}, every
#'       training stratum is returned.
#'   }
#'
#' @examples
#' \donttest{
#' strata_result <- make_strata(maihda_sim_data, vars = c("gender", "race"))
#' model <- fit_maihda(health_outcome ~ age + (1 | stratum), data = strata_result$data)
#' 
#' # Individual predictions
#' pred_ind <- predict_maihda(model, type = "individual")
#' 
#' # Stratum predictions
#' pred_strata <- predict_maihda(model, type = "strata")
#' }
#'
#' @export
#' @importFrom stats predict
#' @importFrom lme4 ranef
predict_maihda <- function(object, newdata = NULL,
                           type = c("individual", "strata", "response", "link"),
                           scale = c("response", "link"), ...) {
  if (!inherits(object, "maihda_model")) {
    stop("'object' must be a maihda_model object from fit_maihda()")
  }
  
  type <- match.arg(type)
  if (type %in% c("response", "link")) {
    scale <- type
    type <- "individual"
  } else {
    scale <- match.arg(scale)
  }
  engine <- object$engine
  model <- object$model
  
  if (is.null(newdata)) {
    newdata <- object$data
  } else {
    newdata <- maihda_prepare_prediction_data(object, newdata)
  }
  
  if (engine == "lme4") {
    if (type == "individual") {
      # Individual-level predictions including random effects
      predictions <- predict(model, newdata = newdata, type = scale, ...)
      return(predictions)
      
    } else if (type == "strata") {
      # Stratum-level predictions (random effects)
      result <- maihda_stratum_ranef_lme4(model)
      result$predicted <- result$random_effect
      result <- result[, c("stratum", "predicted", "se", "lower_95", "upper_95")]
      return(maihda_filter_strata_predictions(result, newdata))
    }
    
  } else if (engine == "brms") {
    # Verify brms is available
    if (!requireNamespace("brms", quietly = TRUE)) {
      stop("Package 'brms' is required to predict from brms models. Please install it with: install.packages('brms')")
    }
    
    if (type == "individual") {
      # Individual-level predictions
      predictions <- if (scale == "response") {
        stats::fitted(model, newdata = newdata, summary = TRUE, ...)[, "Estimate"]
      } else {
        brms::posterior_linpred(model, newdata = newdata, summary = TRUE, ...)[, "Estimate"]
      }
      return(predictions)
      
    } else if (type == "strata") {
      # Stratum-level predictions
      result <- maihda_stratum_ranef_brms(model)
      result$predicted <- result$random_effect
      result <- result[, c("stratum", "predicted", "se", "lower_95", "upper_95")]
      return(maihda_filter_strata_predictions(result, newdata))
    }
  }
}

# Restrict a per-stratum prediction table to the strata present in `newdata` so
# type = "strata" respects newdata the way type = "individual" does (instead of
# always returning every training stratum). A stratum in newdata that the model
# never saw is an error, matching the individual path. When newdata carries no
# stratum column the table is returned unchanged. When it carries a stratum column
# whose values are all missing it names no training stratum to keep, so the result
# is empty -- not silently every training stratum.
maihda_filter_strata_predictions <- function(result, newdata) {
  if (is.null(newdata) || !"stratum" %in% names(newdata)) {
    return(result)
  }
  wanted <- unique(as.character(newdata$stratum))
  wanted <- wanted[!is.na(wanted)]
  if (length(wanted) == 0) {
    return(result[0, , drop = FALSE])
  }
  known <- as.character(result$stratum)
  unknown <- setdiff(wanted, known)
  if (length(unknown) > 0) {
    stop("newdata contains strata not present in the fitted model: ",
         paste(utils::head(unknown, 5), collapse = ", "),
         if (length(unknown) > 5) ", ..." else "",
         call. = FALSE)
  }
  result[known %in% wanted, , drop = FALSE]
}
