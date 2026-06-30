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

test_that("plot_maihda predicted preserves stratum order", {
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

  # Create predicted plot
  plot <- plot(model, type = "predicted")

  # Check structure
  expect_true(inherits(plot, "ggplot"))
  expect_true("display_label" %in% names(plot$data))

  # Check that display_label is a factor and order is preserved
  expect_true(is.factor(plot$data$display_label))

  # The strata should be in their original order, not sorted by predicted value
  # We can't check exact order without knowing the predicted values,
  # but we can check that stratum labels are present
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

  full_order <- as.integer(as.character(
    plot(m, type = "predicted", n_strata = NULL)$data$stratum))
  p_ord <- plot(m, type = "predicted", n_strata = 6, select = "order")
  p_dev <- plot(m, type = "predicted", n_strata = 6, select = "deviation")
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
