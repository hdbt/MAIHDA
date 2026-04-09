library(MAIHDA)
library(dplyr)
library(ggtern)

# Load example data
data("maihda_health_data")

# Create strata
strata_result <- make_strata(maihda_health_data, vars = c("Gender", "Race", "Education"))

# Fit a MAIHDA model
model <- fit_maihda(BMI ~ Age + Gender + Race + Education + Poverty + (1 | stratum),
                    data = strata_result$data,
                    engine = "lme4")

# Generate the ternary plot using ggtern
out <- maihda_ternary_plot(model)

# Print the plot
print(out$plot)
glimpse(out$data)
