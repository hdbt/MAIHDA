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
#' @param ... Additional arguments passed to predict method of underlying model.
#'
#' @return Depending on type:
#'   \itemize{
#'     \item For "individual": A numeric vector of predicted values
#'     \item For "strata": A data frame with stratum ID and predicted random effect
#'   }
#'
#' @examples
#' \dontrun{
#' model <- fit_maihda(outcome ~ age + (1 | stratum), data = data)
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
predict_maihda <- function(object, newdata = NULL, type = c("individual", "strata"), ...) {
  if (!inherits(object, "maihda_model")) {
    stop("'object' must be a maihda_model object from fit_maihda()")
  }
  
  type <- match.arg(type)
  engine <- object$engine
  model <- object$model
  
  if (is.null(newdata)) {
    newdata <- object$data
  }
  
  if (engine == "lme4") {
    if (type == "individual") {
      # Individual-level predictions including random effects
      predictions <- predict(model, newdata = newdata, ...)
      return(predictions)
      
    } else if (type == "strata") {
      # Stratum-level predictions (random effects)
      re <- lme4::ranef(model, condVar = TRUE)
      
      if (length(re) > 0 && "stratum" %in% names(re)) {
        stratum_re <- re[["stratum"]]
        
        # Get conditional variances
        cond_var <- attr(stratum_re, "postVar")
        if (is.array(cond_var) && length(dim(cond_var)) == 3) {
          stratum_se <- sqrt(cond_var[1, 1, ])
        } else {
          stratum_se <- rep(NA, nrow(stratum_re))
        }
        
        result <- data.frame(
          stratum = as.integer(rownames(stratum_re)),
          predicted = stratum_re[, 1],
          se = stratum_se,
          lower_95 = stratum_re[, 1] - 1.96 * stratum_se,
          upper_95 = stratum_re[, 1] + 1.96 * stratum_se
        )
        
        return(result)
      } else {
        stop("No stratum random effects found in model")
      }
    }
    
  } else if (engine == "brms") {
    # Verify brms is available
    if (!requireNamespace("brms", quietly = TRUE)) {
      stop("Package 'brms' is required to predict from brms models. Please install it with: install.packages('brms')")
    }
    
    if (type == "individual") {
      # Individual-level predictions
      predictions <- predict(model, newdata = newdata, ...)
      return(predictions)
      
    } else if (type == "strata") {
      # Stratum-level predictions
      re <- brms::ranef(model)
      
      if (length(re) > 0 && "stratum" %in% names(re)) {
        result <- re[["stratum"]]
        return(result)
      } else {
        stop("No stratum random effects found in model")
      }
    }
  }
}
