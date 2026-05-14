test_that("summary handles binomial and gaussian residual variance correctly", {
  # Setup data
  strata_result <- make_strata(maihda_sim_data, vars = c("gender", "race"))

  # Create a gaussian model
  mod_gauss <- fit_maihda(health_outcome ~ age + (1 | stratum),
                         data = strata_result$data,
                         engine = "lme4", family = "gaussian")

  # Create a binomial model
  strata_result$data$binary_outcome <- ifelse(strata_result$data$health_outcome > mean(strata_result$data$health_outcome), 1, 0)
  mod_binom <- fit_maihda(binary_outcome ~ age + (1 | stratum),
                         data = strata_result$data,
                         engine = "lme4", family = "binomial")

  # Check summaries
  summ_gauss <- summary(mod_gauss)
  summ_binom <- summary(mod_binom)

  # Gaussian residual variance should not be pi^2/3
  gauss_resid_var <- summ_gauss$variance_components$variance[summ_gauss$variance_components$component == "Within-stratum (residual)"]
  expect_true(abs(gauss_resid_var - (pi^2 / 3)) > 0.1) # Shouldn't match pi^2/3

  # Binomial residual variance should be exactly pi^2/3 for logit link
  binom_resid_var <- summ_binom$variance_components$variance[summ_binom$variance_components$component == "Within-stratum (residual)"]
  expect_equal(binom_resid_var, (pi^2) / 3, tolerance = 1e-6)
})
