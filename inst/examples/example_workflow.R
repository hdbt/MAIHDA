# MAIHDA Package Example Workflow
# This script demonstrates the complete MAIHDA workflow

# Load the package
library(MAIHDA)

# Load example data
data("maihda_sim_data")

# View the data
head(maihda_sim_data)
summary(maihda_sim_data)

# ============================================================================
# Step 1: Create intersectional strata
# ============================================================================
strata_result <- make_strata(
  data = maihda_sim_data, 
  vars = c("gender", "race"),
  min_n = 5  # Exclude strata with fewer than 5 observations
)

# View stratum information
print(strata_result)

# ============================================================================
# Step 2: Fit MAIHDA model
# ============================================================================
model <- fit_maihda(
  formula = health_outcome ~ age + (1 | stratum),
  data = strata_result$data,
  engine = "lme4",
  family = "gaussian"
)

# View model
print(model)

# ============================================================================
# Step 3: Summarize the model
# ============================================================================

# Basic summary
summary_result <- summary_maihda(model)
print(summary_result)

# Summary with bootstrap confidence intervals (takes longer)
# summary_boot <- summary_maihda(model, bootstrap = TRUE, n_boot = 500)
# print(summary_boot)

# ============================================================================
# Step 4: Make predictions
# ============================================================================

# Individual-level predictions
pred_individual <- predict_maihda(model, type = "individual")
head(pred_individual)

# Stratum-level predictions
pred_strata <- predict_maihda(model, type = "strata")
print(pred_strata)

# ============================================================================
# Step 5: Visualize results
# ============================================================================

# Caterpillar plot of stratum random effects
plot_caterpillar <- plot_maihda(model, type = "caterpillar")
print(plot_caterpillar)

# Variance partition coefficient visualization
plot_vpc <- plot_maihda(model, type = "vpc")
print(plot_vpc)

# Observed vs. shrunken estimates
plot_obs_shrunk <- plot_maihda(model, type = "obs_vs_shrunken")
print(plot_obs_shrunk)

# ============================================================================
# Step 6: Compare models (optional)
# ============================================================================

# Fit a second model with education added
strata_result2 <- make_strata(
  data = maihda_sim_data,
  vars = c("gender", "race", "education")
)

model2 <- fit_maihda(
  formula = health_outcome ~ age + (1 | stratum),
  data = strata_result2$data,
  engine = "lme4"
)

# Compare models
comparison <- compare_maihda(
  model, model2,
  model_names = c("Gender x Race", "Gender x Race x Education"),
  bootstrap = FALSE  # Set to TRUE for bootstrap CI
)

print(comparison)

# Plot comparison
plot_comp <- plot_comparison(comparison)
print(plot_comp)
