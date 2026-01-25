\name{maihda_sim_data}
\alias{maihda_sim_data}
\docType{data}
\title{
Simulated Health Data for MAIHDA Examples
}
\description{
A simulated dataset containing health outcomes and demographic variables for 500 individuals.
This dataset is designed to demonstrate intersectional health inequalities suitable for
MAIHDA analysis. The data includes main effects and intersectional effects between gender,
race, and education.
}
\usage{data("maihda_sim_data")}
\format{
  A data frame with 500 observations on the following 6 variables:
  \describe{
    \item{\code{id}}{Integer identifier for each individual}
    \item{\code{gender}}{Character variable indicating gender ("Male" or "Female")}
    \item{\code{race}}{Character variable indicating race/ethnicity ("White", "Black", "Hispanic", or "Asian")}
    \item{\code{education}}{Character variable indicating education level ("High School", "Some College", "College", or "Graduate")}
    \item{\code{age}}{Numeric variable for age in years (range: 18-80)}
    \item{\code{health_outcome}}{Numeric variable representing a health score (higher is better)}
  }
}
\details{
The health outcome was simulated with:
\itemize{
  \item Main effects of gender, race, education, and age
  \item Intersectional effects (e.g., Black women, men with high school education)
  \item Random noise with standard deviation of 10
}

The data demonstrates typical patterns in health inequalities research where outcomes vary
both by individual characteristics and their intersections.
}
\examples{
data(maihda_sim_data)

# View structure
str(maihda_sim_data)

# Summary statistics
summary(maihda_sim_data)

\dontrun{
# Example MAIHDA analysis
library(MAIHDA)

# Create strata
strata_result <- make_strata(maihda_sim_data, vars = c("gender", "race"))

# Fit model
model <- fit_maihda(health_outcome ~ age + (1 | stratum), 
                   data = strata_result$data)

# Summarize
summary_maihda(model)

# Visualize
plot_maihda(model, type = "caterpillar")
}
}
\keyword{datasets}
