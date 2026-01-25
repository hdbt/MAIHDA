\name{predict_maihda}
\alias{predict_maihda}
\title{Predict from MAIHDA Model}
\description{
Makes predictions from a fitted MAIHDA model, either at the stratum level
or individual level.
}
\usage{
predict_maihda(object, newdata = NULL, type = c("individual", "strata"), ...)
}
\arguments{
  \item{object}{A maihda_model object from \code{fit_maihda()}.}
  \item{newdata}{Optional data frame for making predictions. If NULL, uses the original data from model fitting.}
  \item{type}{Character string specifying prediction type: "individual" for individual-level predictions including random effects, or "strata" for stratum-level predictions (random effects only).}
  \item{...}{Additional arguments passed to predict method of underlying model.}
}
\value{
Depending on type:
\itemize{
  \item For "individual": A numeric vector of predicted values
  \item For "strata": A data frame with stratum ID and predicted random effect
}
}
\examples{
\dontrun{
model <- fit_maihda(outcome ~ age + (1 | stratum), data = data)

# Individual predictions
pred_ind <- predict_maihda(model, type = "individual")

# Stratum predictions
pred_strata <- predict_maihda(model, type = "strata")
}
}
