test_that("plot_maihda creates predicted plot with lme4", {
  # Create test data
  set.seed(123)
  data <- data.frame(
    stratum = rep(1:10, each = 10),
    age = rnorm(100),
    outcome = rnorm(100)
  )

  # Fit model
  model <- fit_maihda(outcome ~ age + (1 | stratum),
                     data = data,
                     engine = "lme4")

  # Create predicted plot
  plot <- plot(model, type = "predicted")

  # Check structure
  expect_true(inherits(plot, "ggplot"))
  expect_true("predicted" %in% names(plot$data))
  expect_true("lower" %in% names(plot$data))
  expect_true("upper" %in% names(plot$data))
})

test_that("plot_maihda predicted handles n_strata parameter", {
  # Create test data with more strata
  set.seed(456)
  data <- data.frame(
    stratum = rep(1:20, each = 10),
    age = rnorm(200),
    outcome = rnorm(200)
  )

  # Fit model
  model <- fit_maihda(outcome ~ age + (1 | stratum),
                     data = data,
                     engine = "lme4")

  # Create predicted plot with limited strata
  plot <- plot(model, type = "predicted", n_strata = 10)

  # Check that it limits to specified number
  expect_true(inherits(plot, "ggplot"))
  expect_lte(nrow(plot$data), 10)
})

test_that("plot_maihda predicted validates inputs", {
  # Create test data
  set.seed(789)
  data <- data.frame(
    stratum = rep(1:5, each = 10),
    age = rnorm(50),
    outcome = rnorm(50)
  )

  # Fit model
  model <- fit_maihda(outcome ~ age + (1 | stratum),
                     data = data,
                     engine = "lme4")

  # Check plot is created
  plot <- plot(model, type = "predicted")
  expect_true(inherits(plot, "ggplot"))
})

test_that("plot_maihda predicted order_by = 'stratum' preserves native order", {
  # Create test data with stratum labels
  set.seed(999)
  data <- data.frame(
    stratum = factor(rep(c("A", "B", "C", "D"), each = 10),
                     levels = c("A", "B", "C", "D")),
    age = rnorm(40),
    outcome = rnorm(40)
  )

  # Fit model
  model <- fit_maihda(outcome ~ age + (1 | stratum),
                     data = data,
                     engine = "lme4")

  # order_by = "stratum" opts out of the ranked default and keeps native order
  plot <- plot(model, type = "predicted", order_by = "stratum")

  # Check structure
  expect_true(inherits(plot, "ggplot"))
  expect_true("display_label" %in% names(plot$data))

  # Check that display_label is a factor and order is preserved
  expect_true(is.factor(plot$data$display_label))

  # The displayed rows run in the native stratum order (that of the summary's
  # stratum estimates), NOT sorted by predicted value.
  se <- summary(model)$stratum_estimates
  expect_equal(as.character(plot$data$stratum), as.character(se$stratum))
  expect_true(all(c("A", "B", "C", "D") %in% levels(plot$data$display_label)))
})

test_that("plot_maihda uses meaningful stratum labels from make_strata", {
  # Create test data with meaningful categorical variables
  set.seed(555)
  data <- data.frame(
    gender = rep(c("Male", "Female"), each = 20),
    race = rep(c("White", "Black"), times = 20),
    age = rnorm(40),
    outcome = rnorm(40)
  )

  # Use make_strata to create labeled strata
  strata_result <- make_strata(data, vars = c("gender", "race"))

  # Verify strata_info has labels
  expect_true("label" %in% names(strata_result$strata_info))
  expect_true(any(grepl("\u00d7", strata_result$strata_info$label)))

  # Fit model using data from make_strata
  model <- fit_maihda(outcome ~ age + (1 | stratum),
                     data = strata_result$data,
                     engine = "lme4")

  # Verify model has strata_info
  expect_false(is.null(model$strata_info))
  expect_true("label" %in% names(model$strata_info))

  # Get summary
  summary_obj <- summary(model)

  # Verify summary has labels
  expect_true("label" %in% names(summary_obj$stratum_estimates))

  # Create predicted plot
  plot <- plot(model, type = "predicted")

  # Check structure
  expect_true(inherits(plot, "ggplot"))
  expect_true("display_label" %in% names(plot$data))

  # Check that meaningful labels are used (should contain multiplication signs from gender \u00d7 race)
  display_labels <- as.character(plot$data$display_label)
  expect_true(any(grepl("\u00d7", display_labels)),
             info = "Plot should use meaningful labels like 'Male \u00d7 White', not numeric IDs")
})

test_that("observed vs shrunken handles binary factor outcomes", {
  set.seed(556)
  strata_result <- make_strata(maihda_sim_data, vars = c("gender", "race"))

  expect_warning(
    model <- fit_maihda(binary_outcome ~ age + (1 | stratum),
                        data = strata_result$data,
                        engine = "lme4"),
    "Automatically switching to family = 'binomial'",
    fixed = TRUE
  )

  plot <- plot(model, type = "obs_vs_shrunken")

  expect_true(inherits(plot, "ggplot"))
  expect_true(all(is.finite(plot$data$observed)))
  expect_true(all(plot$data$observed >= 0 & plot$data$observed <= 1))
})

test_that("observed vs shrunken handles two-column binomial outcomes", {
  set.seed(558)
  data <- data.frame(
    stratum = factor(rep(seq_len(8), each = 25)),
    age = rnorm(200)
  )
  stratum_effect <- rnorm(8, sd = 0.7)[data$stratum]
  probability <- plogis(-0.3 + 0.5 * data$age + stratum_effect)
  trials <- sample(1:25, nrow(data), replace = TRUE)
  data$success <- rbinom(nrow(data), size = trials, prob = probability)
  data$failure <- trials - data$success

  model <- fit_maihda(cbind(success, failure) ~ age + (1 | stratum),
                      data = data,
                      family = "binomial")
  plot <- plot(model, type = "obs_vs_shrunken")

  expected <- stats::aggregate(
    cbind(success = data$success, total = data$success + data$failure),
    by = list(stratum = as.character(data$stratum)),
    FUN = sum
  )
  expected$observed <- expected$success / expected$total
  idx <- match(as.character(plot$data$stratum), expected$stratum)

  expect_true(inherits(plot, "ggplot"))
  expect_true(all(is.finite(plot$data$observed)))
  expect_equal(plot$data$observed, expected$observed[idx], tolerance = 1e-8)
})

test_that("select = 'deviation' keeps the most extreme strata (both tails), not the first", {
  set.seed(303)
  K <- 20L; per <- 30L
  stratum <- rep(seq_len(K), each = per)
  # Stratum mean rises with id; the overall mean sits at the centre (~10.5), so the
  # most extreme strata by |predicted - reference| are the lowest AND highest ids.
  outcome <- stratum + rnorm(K * per, sd = 0.3)
  dat <- data.frame(stratum = stratum, outcome = outcome)
  m <- fit_maihda(outcome ~ 1 + (1 | stratum), data = dat, engine = "lme4")

  # Pin order_by = "stratum" so display stays in native order -- this test is
  # about which strata `select` keeps, not how `order_by` arranges them.
  full_order <- as.integer(as.character(
    plot(m, type = "predicted", n_strata = NULL, order_by = "stratum")$data$stratum))
  p_ord <- plot(m, type = "predicted", n_strata = 6, select = "order", order_by = "stratum")
  p_dev <- plot(m, type = "predicted", n_strata = 6, select = "deviation", order_by = "stratum")
  ord_shown <- as.integer(as.character(p_ord$data$stratum))
  dev_shown <- as.integer(as.character(p_dev$data$stratum))

  # order keeps the first 6 in the native (stratum) order
  expect_equal(ord_shown, utils::head(full_order, 6))
  # deviation keeps the 6 furthest from the centre -- three from each tail
  expect_setequal(dev_shown, c(1, 2, 3, 18, 19, 20))
  # selection != display order: the survivors still appear in native stratum order
  expect_equal(dev_shown, full_order[full_order %in% dev_shown])
  # captions name the rule actually used
  expect_match(p_dev$labels$caption, "furthest from the reference")
  expect_match(p_ord$labels$caption, "first 6 of 20")
})

test_that("select defaults to 'order' and is validated", {
  set.seed(404)
  dat <- data.frame(stratum = rep(1:8, each = 10), outcome = rnorm(80))
  m <- fit_maihda(outcome ~ 1 + (1 | stratum), data = dat, engine = "lme4")
  # default == explicit "order"
  d_default <- as.character(plot(m, type = "predicted", n_strata = 4)$data$stratum)
  d_order   <- as.character(plot(m, type = "predicted", n_strata = 4, select = "order")$data$stratum)
  expect_equal(d_default, d_order)
  expect_error(plot(m, type = "predicted", select = "nope"))
})

# ---- order_by: prediction-based display ordering (issue #57) -----------------

# A model whose stratum predictions rise monotonically with the stratum id, so the
# ranked orderings are unambiguous.
maihda_monotone_pred_model <- function(seed, K = 12L, per = 30L) {
  set.seed(seed)
  stratum <- rep(seq_len(K), each = per)
  outcome <- stratum + rnorm(K * per, sd = 0.3)
  fit_maihda(outcome ~ 1 + (1 | stratum),
             data = data.frame(stratum = stratum, outcome = outcome),
             engine = "lme4")
}

test_that("order_by defaults to 'predicted_desc' (ranked caterpillar)", {
  m <- maihda_monotone_pred_model(1201)
  p_default <- plot(m, type = "predicted")
  p_desc    <- plot(m, type = "predicted", order_by = "predicted_desc")

  # the default is predicted_desc
  expect_equal(as.character(p_default$data$stratum), as.character(p_desc$data$stratum))
  # rows run from highest predicted (top of the flipped axis) to lowest
  expect_equal(p_default$data$predicted,
               sort(p_default$data$predicted, decreasing = TRUE))
})

test_that("order_by = 'predicted_asc' reverses the default order", {
  m <- maihda_monotone_pred_model(1202)
  p_asc  <- plot(m, type = "predicted", order_by = "predicted_asc")
  p_desc <- plot(m, type = "predicted", order_by = "predicted_desc")

  expect_equal(p_asc$data$predicted,
               sort(p_asc$data$predicted, decreasing = FALSE))
  # asc is the exact reverse of desc
  expect_equal(as.character(p_asc$data$stratum),
               rev(as.character(p_desc$data$stratum)))
})

test_that("order_by = 'stratum' keeps the native summary order", {
  m <- maihda_monotone_pred_model(1203)
  p <- plot(m, type = "predicted", order_by = "stratum")
  se <- summary(m)$stratum_estimates
  expect_equal(as.character(p$data$stratum), as.character(se$stratum))
})

test_that("order_by = 'deviation' ranks by distance from the reference line", {
  m <- maihda_monotone_pred_model(1204)
  p <- plot(m, type = "predicted", order_by = "deviation")

  # recover the dashed reference line (the geom_hline yintercept) from the build
  bl <- ggplot2::ggplot_build(p)
  yint <- NULL
  for (ld in bl$data) if ("yintercept" %in% names(ld)) yint <- ld$yintercept[1]
  expect_false(is.null(yint))

  d <- abs(p$data$predicted - yint)
  expect_equal(d, sort(d, decreasing = TRUE))
})

test_that("order_by is display-only: same strata and values, different order", {
  m <- maihda_monotone_pred_model(1205, K = 20L, per = 25L)

  # cap so selection is exercised; select governs WHICH strata, order_by the order
  p_str <- plot(m, type = "predicted", n_strata = 8, select = "order", order_by = "stratum")
  p_dsc <- plot(m, type = "predicted", n_strata = 8, select = "order", order_by = "predicted_desc")

  # identical set of strata -- order_by does not change selection
  expect_setequal(as.character(p_str$data$stratum), as.character(p_dsc$data$stratum))
  # identical predicted values / intervals per stratum -- only the row order moved
  idx <- match(as.character(p_str$data$stratum), as.character(p_dsc$data$stratum))
  expect_equal(p_str$data$predicted, p_dsc$data$predicted[idx])
  expect_equal(p_str$data$lower, p_dsc$data$lower[idx])
  expect_equal(p_str$data$upper, p_dsc$data$upper[idx])
  # the display order genuinely differs
  expect_false(identical(as.character(p_str$data$stratum),
                         as.character(p_dsc$data$stratum)))
})

test_that("order_by keeps the meaningful make_strata labels", {
  set.seed(1206)
  data <- data.frame(
    gender = rep(c("Male", "Female"), each = 20),
    race = rep(c("White", "Black"), times = 20),
    age = rnorm(40),
    outcome = rnorm(40)
  )
  sr <- make_strata(data, vars = c("gender", "race"))
  m <- fit_maihda(outcome ~ age + (1 | stratum), data = sr$data, engine = "lme4")
  p <- plot(m, type = "predicted", order_by = "predicted_desc")
  labs <- as.character(p$data$display_label)
  times <- intToUtf8(0x00d7)  # the multiplication-sign separator in make_strata() labels
  expect_true(any(grepl(times, labs, fixed = TRUE)))
})

test_that("order_by is validated", {
  set.seed(1207)
  dat <- data.frame(stratum = rep(1:8, each = 10), outcome = rnorm(80))
  m <- fit_maihda(outcome ~ 1 + (1 | stratum), data = dat, engine = "lme4")
  expect_error(plot(m, type = "predicted", order_by = "nope"))
})
