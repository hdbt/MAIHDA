\name{plot_maihda}
\alias{plot_maihda}
\title{Plot MAIHDA Model Results}
\description{
Creates various plots for visualizing MAIHDA model results including caterpillar
plots, variance partition coefficient comparisons, observed vs. shrunken estimates,
and predicted subgroup values with confidence intervals.
}
\usage{
plot_maihda(object, type = c("caterpillar", "vpc", "obs_vs_shrunken", "predicted"),
            summary_obj = NULL, n_strata = 50, ...)
}
\arguments{
  \item{object}{A maihda_model object from \code{fit_maihda()}.}
  \item{type}{Character string specifying plot type: "caterpillar" for caterpillar plot of stratum random effects, "vpc" for variance partition coefficient visualization, "obs_vs_shrunken" for observed vs. shrunken stratum means, or "predicted" for predicted values for each stratum with confidence intervals.}
  \item{summary_obj}{Optional maihda_summary object from \code{summary_maihda()}. If NULL, will be computed.}
  \item{n_strata}{Maximum number of strata to display in caterpillar plot or predicted plot. Default is 50. Use NULL for all strata.}
  \item{...}{Additional arguments (not currently used).}
}
\value{
A ggplot2 object.
}
\examples{
\dontrun{
model <- fit_maihda(outcome ~ age + (1 | stratum), data = data)

# Caterpillar plot
plot_maihda(model, type = "caterpillar")

# VPC plot
plot_maihda(model, type = "vpc")

# Observed vs shrunken plot
plot_maihda(model, type = "obs_vs_shrunken")

# Predicted values with confidence intervals
plot_maihda(model, type = "predicted")
}
}
