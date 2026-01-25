\name{summary_maihda}
\alias{summary_maihda}
\title{Summarize MAIHDA Model}
\description{
Provides a summary of a MAIHDA model including variance partition coefficients
(VPC/ICC) and stratum-specific estimates.
}
\usage{
summary_maihda(object, bootstrap = FALSE, n_boot = 1000, 
               conf_level = 0.95, ...)
}
\arguments{
  \item{object}{A maihda_model object from \code{fit_maihda()}.}
  \item{bootstrap}{Logical indicating whether to compute bootstrap confidence intervals for VPC/ICC. Default is FALSE.}
  \item{n_boot}{Number of bootstrap samples if bootstrap = TRUE. Default is 1000.}
  \item{conf_level}{Confidence level for bootstrap intervals. Default is 0.95.}
  \item{...}{Additional arguments (not currently used).}
}
\value{
A maihda_summary object containing:
\describe{
  \item{vpc}{Variance Partition Coefficient (ICC) with optional CI}
  \item{variance_components}{Data frame of variance components}
  \item{stratum_estimates}{Data frame of stratum-specific random effects}
  \item{fixed_effects}{Fixed effects estimates}
  \item{model_summary}{Original model summary}
}
}
\examples{
\dontrun{
model <- fit_maihda(outcome ~ age + (1 | stratum), data = data)
summary_result <- summary_maihda(model)

# With bootstrap CI
summary_boot <- summary_maihda(model, bootstrap = TRUE, n_boot = 500)
}
}
