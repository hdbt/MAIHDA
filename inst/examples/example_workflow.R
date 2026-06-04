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
summary_result <- summary(model)
print(summary_result)

# Summary with bootstrap confidence intervals (takes longer)
# summary_boot <- summary(model, bootstrap = TRUE, n_boot = 500)
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

# Predicted stratum values with confidence intervals
plot_predicted <- plot(model, type = "predicted")
print(plot_predicted)

# Variance partition coefficient visualization
plot_vpc <- plot(model, type = "vpc")
print(plot_vpc)

# Observed vs. shrunken estimates
plot_obs_shrunk <- plot(model, type = "obs_vs_shrunken")
print(plot_obs_shrunk)

# ============================================================================
# Step 6: Compare models (optional)
# ============================================================================

# Fit an adjusted model on the SAME data and strata that adds the additive main
# effects of the strata variables. Comparing the unadjusted vs additive-adjusted
# model is the canonical MAIHDA contrast: it shows how much of the between-stratum
# variance is explained by the additive (main-effect) components, leaving the
# remainder as the intersectional interaction. (Comparing models built on
# DIFFERENT strata definitions would not be a valid PCV.)
model2 <- fit_maihda(
  formula = health_outcome ~ age + gender + race + (1 | stratum),
  data = strata_result$data,
  engine = "lme4"
)

# Compare nested models on the same data and strata.
comparison <- compare_maihda(
  model, model2,
  model_names = c("Unadjusted", "Additive-adjusted"),
  bootstrap = FALSE  # Set to TRUE for bootstrap CI
)

print(comparison)

# Plot comparison
plot_comp <- plot(comparison)
print(plot_comp)
