\name{fit_maihda}
\alias{fit_maihda}
\title{Fit MAIHDA Model}
\description{
Fits a multilevel model for MAIHDA (Multilevel Analysis of Individual
Heterogeneity and Discriminatory Accuracy) using either lme4 or brms.
}
\usage{
fit_maihda(formula, data, engine = "lme4", family = "gaussian", ...)
}
\arguments{
  \item{formula}{A formula specifying the model. Should include random effect for stratum (e.g., \code{outcome ~ fixed_vars + (1 | stratum)}).}
  \item{data}{A data frame containing the variables in the formula.}
  \item{engine}{Character string specifying which engine to use: "lme4" (default) or "brms".}
  \item{family}{Character string or family object specifying the model family. Common options: "gaussian", "binomial", "poisson". Default is "gaussian".}
  \item{...}{Additional arguments passed to \code{lmer}/\code{glmer} (lme4) or \code{brm} (brms).}
}
\value{
A maihda_model object containing:
\describe{
  \item{model}{The fitted model object (lme4 or brms)}
  \item{engine}{The engine used ("lme4" or "brms")}
  \item{formula}{The model formula}
  \item{data}{The data used for fitting}
  \item{family}{The family used}
}
}
\examples{
\dontrun{
# Create strata
strata_result <- make_strata(data, vars = c("gender", "race"))

# Fit model with lme4
model <- fit_maihda(outcome ~ age + (1 | stratum),
                    data = strata_result$data,
                    engine = "lme4")

# Fit model with brms (if installed)
model_brms <- fit_maihda(outcome ~ age + (1 | stratum),
                         data = strata_result$data,
                         engine = "brms")
}
}
