test_that("Ternary plot functions work", {
  skip_if_not_installed("lme4")

  # Minimal dummy data
  set.seed(123)
  df <- data.frame(
    stratum = rep(letters[1:5], each = 10),
    y = rnorm(50),
    x = rnorm(50)
  )

  fit <- lme4::lmer(y ~ x + (1 | stratum), data = df)

  model <- list(
    model = fit,
    engine = "lme4",
    data = df
  )
  class(model) <- "maihda_model"

  # Test compute function
  td <- compute_maihda_ternary_data(model)
  expect_s3_class(td, "tbl_df")
  expect_true(all(c("stratum", "additive_prop", "interaction_prop", "uncertainty_prop") %in% names(td)))
  expect_equal(rowSums(td[, c("additive_prop", "interaction_prop", "uncertainty_prop")]), rep(1, 5), tolerance = 1e-4)

  td_ci <- compute_maihda_ternary_data(model, uncertainty_method = "ci_width")
  expect_equal(td_ci$uncertainty, td$uncertainty * 3.92, tolerance = 1e-8)
  expect_error(
    compute_maihda_ternary_data(model, uncertainty_method = "posterior_sd"),
    "only available for brms"
  )

  # Test wrapper function
  skip_if_not_installed("ggtern")
  out <- maihda_ternary_plot(model)
  expect_type(out, "list")
  expect_true(!is.null(out$plot))
  expect_s3_class(out$data, "tbl_df")
})

test_that("response-scale ternary interaction signal uses response-scale differences", {
  skip_if_not_installed("lme4")

  set.seed(124)
  df <- data.frame(
    stratum = factor(rep(seq_len(8), each = 30)),
    x = rnorm(240)
  )
  stratum_u <- rnorm(8, sd = 0.9)[df$stratum]
  df$y <- rbinom(nrow(df), 1, stats::plogis(-0.2 + 0.6 * df$x + stratum_u))

  model <- fit_maihda(y ~ x + (1 | stratum), data = df, family = "binomial")
  td <- suppressWarnings(compute_maihda_ternary_data(model, scale = "response", verbose = FALSE))

  fe_link <- stats::predict(model$model, newdata = model$data, re.form = NA, type = "link")
  u_by_row <- td$u_j[match(as.character(model$data$stratum), as.character(td$stratum))]
  expected <- stats::aggregate(
    list(
      additive_only = stats::plogis(fe_link),
      full_prediction = stats::plogis(fe_link + u_by_row)
    ),
    by = list(stratum = as.character(model$data$stratum)),
    FUN = mean
  )
  expected$interaction_signal <- abs(expected$full_prediction - expected$additive_only)

  idx <- match(as.character(td$stratum), expected$stratum)
  expect_equal(td$additive_only, expected$additive_only[idx], tolerance = 1e-8)
  expect_equal(td$interaction_signal, expected$interaction_signal[idx], tolerance = 1e-8)
  expect_false(isTRUE(all.equal(td$interaction_signal, abs(td$u_j), tolerance = 1e-4)))
})

test_that("ternary diagnostic blocks brms cumulative on both scales (Stan-free)", {
  skip_if_not_installed("brms")

  # Regression for the audit finding: compute_maihda_ternary_data() blocked only
  # engine == "ordinal" (clmm), so a brms::cumulative() fit slipped through to the
  # scalar linkinv() path and produced cumulative probabilities in [0, 1] on the
  # response scale instead of the expected category score in [1, K] -- the same
  # wrong-scale bug fixed in maihda_stratum_predictions_brms(). The ternary
  # decomposition is a latent-scale concept the package does not define a coherent
  # response-scale version of for a cumulative model, so BOTH cumulative engines
  # are now blocked consistently. The guard reads only model$engine and
  # model$family, so a bare fake fit exercises it with no Stan.
  m <- structure(
    list(
      model  = structure(list(), class = "brmsfit"),
      engine = "brms",
      family = brms::cumulative("logit")
    ),
    class = "maihda_model"
  )
  expect_error(compute_maihda_ternary_data(m, scale = "link"),
               "not supported for cumulative \\(ordinal\\)")
  # On the response scale the block fires BEFORE the "most coherent on the link
  # scale" warning, so it errors cleanly rather than warning first.
  expect_error(compute_maihda_ternary_data(m, scale = "response"),
               "not supported for cumulative \\(ordinal\\)")

  # A non-cumulative brms family is NOT caught by the guard (it proceeds past it):
  # the family predicate the guard uses distinguishes the two.
  expect_true(maihda_family_is_ordinal(brms::cumulative("logit")))
  expect_false(maihda_family_is_ordinal(brms::brmsfamily("gaussian")))
})
