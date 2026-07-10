make_describe_data <- function(seed = 8001, n = 300) {
  set.seed(seed)
  d <- data.frame(
    gender = sample(c("F", "M"), n, replace = TRUE),
    race = sample(c("X", "Y", "Z"), n, replace = TRUE),
    age = rnorm(n)
  )
  sk <- interaction(d$gender, d$race, drop = TRUE)
  d$y <- 1 + 0.3 * d$age + rnorm(nlevels(sk), sd = 0.8)[sk] + rnorm(n, sd = 0.5)
  d$bin <- rbinom(n, 1, plogis(-0.2 + 0.5 * d$age))
  d
}

test_that("maihda_describe() builds strata identical to make_strata()", {
  d <- make_describe_data()
  desc <- maihda_describe(y ~ age + (1 | gender:race), data = d)

  expect_s3_class(desc, "maihda_describe")
  expect_s3_class(desc$overview, "data.frame")
  expect_s3_class(desc$strata, "data.frame")

  s <- make_strata(d, vars = c("gender", "race"))
  obs <- desc$strata[!desc$strata$empty, , drop = FALSE]
  expect_identical(obs$stratum, s$strata_info$stratum)
  expect_identical(obs$label, s$strata_info$label)
  expect_identical(obs$n, s$strata_info$n)
  expect_identical(as.character(obs$gender), as.character(s$strata_info$gender))

  # Complete data: everything is analytic, nothing empty.
  expect_equal(desc$overview$n_total, nrow(d))
  expect_equal(desc$overview$n_missing_outcome, 0)
  expect_equal(desc$overview$n_rows_missing_dimensions, 0)
  expect_equal(desc$overview$n_analytic, nrow(d))
  expect_equal(desc$overview$n_strata_expected, 6)   # 2 x 3
  expect_equal(desc$overview$n_empty_strata, 0)

  # Dimension counts sum to the sample and match table().
  expect_equal(sum(desc$dimensions$n[desc$dimensions$dimension == "gender"]),
               nrow(d))
  g <- desc$dimensions[desc$dimensions$dimension == "gender", ]
  expect_equal(setNames(g$n, g$level), setNames(as.integer(table(d$gender)),
                                                names(table(d$gender))))

  # Gaussian summary applied, and the overall summary matches base R.
  expect_true(all(c("outcome_mean", "outcome_sd", "outcome_median") %in%
                    names(desc$strata)))
  expect_equal(desc$outcome_overall$outcome_mean, mean(d$y), tolerance = 1e-12)
  expect_match(desc$outcome_summary, "gaussian")
})

test_that("missing-data accounting separates dims, outcome, covariates and matches the fit", {
  d <- make_describe_data(8002)
  d$race[1:5] <- NA        # missing dimension -> outside strata
  d$y[6:12] <- NA          # missing outcome
  d$age[13:15] <- NA       # missing covariate (analytic drop only)

  desc <- maihda_describe(y ~ age + (1 | gender:race), data = d)
  expect_equal(desc$overview$n_rows_missing_dimensions, 5)
  expect_equal(desc$overview$n_missing_outcome, 7)
  expect_equal(desc$overview$n_analytic, nrow(d) - 15)

  # The analytic sample equals the rows a subsequent fit actually uses,
  # per stratum and in total.
  m <- fit_maihda(y ~ age + (1 | gender:race), data = d)
  expect_equal(desc$overview$n_analytic, nrow(m$data))
  obs <- desc$strata[!desc$strata$empty, , drop = FALSE]
  expect_identical(obs$n_analytic, m$strata_info$n)

  # Per-variable missingness table.
  expect_equal(desc$missingness$n_missing[desc$missingness$variable == "y"], 7)
  expect_equal(desc$missingness$n_missing[desc$missingness$variable == "race"], 5)
  expect_equal(desc$missingness$n_missing[desc$missingness$variable == "gender"], 0)
  expect_true("missing_dimensions" %in% desc$warnings$check)
  expect_match(desc$warnings$message[desc$warnings$check == "missing_dimensions"],
               "race 5")

  # Per-stratum outcome missingness adds up to the overall count.
  expect_equal(sum(obs$n_missing_outcome), 7)
})

test_that("outcome summaries are family-aware (binomial vs gaussian override)", {
  d <- make_describe_data(8003)

  # 0/1 outcome auto-detects binomial: proportions, never a Gaussian SD.
  desc_bin <- maihda_describe(bin ~ (1 | gender:race), data = d)
  expect_identical(desc_bin$family, "binomial")
  expect_true(desc_bin$family_detected)
  expect_true(all(c("outcome_events", "outcome_trials", "outcome_proportion") %in%
                    names(desc_bin$strata)))
  expect_false(any(c("outcome_sd", "outcome_mean") %in% names(desc_bin$strata)))
  expect_equal(desc_bin$outcome_overall$outcome_proportion, mean(d$bin),
               tolerance = 1e-12)
  expect_equal(nrow(desc_bin$outcome_levels), 2)
  expect_identical(desc_bin$event_level, "1")

  # Explicit family = "gaussian" (a linear probability model) is honoured.
  desc_lpm <- maihda_describe(bin ~ (1 | gender:race), data = d,
                              family = "gaussian")
  expect_identical(desc_lpm$family, "gaussian")
  expect_false(desc_lpm$family_detected)
  expect_true("outcome_mean" %in% names(desc_lpm$strata))
  expect_false("outcome_proportion" %in% names(desc_lpm$strata))

  # A factor outcome under an explicit gaussian family cannot be summarised.
  d$binf <- factor(ifelse(d$bin == 1, "yes", "no"))
  expect_error(maihda_describe(binf ~ (1 | gender:race), data = d,
                               family = "gaussian"),
               "cannot summarise")
})

test_that("aggregated cbind() binomial outcomes are summarised by events/trials", {
  d <- make_describe_data(8004, n = 200)
  d$s <- rbinom(nrow(d), 5, 0.3)
  d$f <- 5L - d$s
  desc <- maihda_describe(cbind(s, f) ~ (1 | gender:race), data = d,
                          family = "binomial")
  expect_identical(desc$family, "binomial")
  expect_equal(desc$outcome_overall$outcome_events, sum(d$s))
  expect_equal(desc$outcome_overall$outcome_trials, 5 * nrow(d))
  expect_null(desc$outcome_levels)
  expect_equal(desc$overview$n_missing_outcome, 0)

  # The same shorthand used to error inside fit_maihda() as well (nobars()
  # returns a bare call for a cbind response with a bars-only RHS); the shared
  # resolver fix covers the fit too, and the description matches it.
  m <- fit_maihda(cbind(s, f) ~ (1 | gender:race), data = d,
                  family = "binomial")
  obs <- desc$strata[!desc$strata$empty, , drop = FALSE]
  expect_identical(obs$label, m$strata_info$label)
  expect_equal(desc$overview$n_analytic, nrow(m$data))
})

test_that("an ordered-factor outcome is described on the category-score scale", {
  d <- make_describe_data(8005, n = 240)
  d$rating <- factor(sample(c("Low", "Mid", "High"), nrow(d), replace = TRUE,
                            prob = c(0.3, 0.5, 0.2)),
                     levels = c("Low", "Mid", "High"), ordered = TRUE)
  desc <- maihda_describe(rating ~ (1 | gender:race), data = d)
  expect_identical(desc$family, "cumulative")
  expect_true(desc$family_detected)
  expect_true(all(c("outcome_mean_score", "outcome_median_category") %in%
                    names(desc$strata)))
  obs <- desc$strata[!desc$strata$empty, , drop = FALSE]
  expect_true(all(obs$outcome_median_category %in% c("Low", "Mid", "High")))
  expect_equal(nrow(desc$outcome_levels), 3)
  expect_identical(desc$outcome_levels$level, c("Low", "Mid", "High"))
  expect_equal(desc$outcome_overall$outcome_mean_score,
               mean(as.integer(d$rating)), tolerance = 1e-12)
})

test_that("context units are tabulated and weak identification is flagged", {
  data(maihda_country_data)
  desc <- maihda_describe(math ~ escs + (1 | gender:ses),
                          data = maihda_country_data, context = "country")
  expect_s3_class(desc$context, "data.frame")
  expect_equal(nrow(desc$context), 6)
  expect_equal(sum(desc$context$n), nrow(maihda_country_data))
  expect_true(all(desc$context$n_analytic <= desc$context$n))
  expect_true("context_few_levels" %in% desc$warnings$check)

  # The same role-clash validation as fit_maihda().
  expect_error(maihda_describe(math ~ (1 | gender:ses),
                               data = maihda_country_data, context = "gender"),
               "also define the intersectional strata")
  expect_error(maihda_describe(math ~ escs + (1 | gender:ses),
                               data = maihda_country_data, context = "escs"),
               "fixed part")
})

test_that("empty and small strata are enumerated and flagged, never dropped", {
  data(maihda_sparse_data)
  # Remove the smallest observed cell entirely so the strata space has a real
  # hole (every dimension level survives in other cells).
  key <- with(maihda_sparse_data,
              paste(gender, ethnicity, education, age_group, sep = "\r"))
  drop_key <- names(sort(table(key)))[1]
  sp <- maihda_sparse_data[key != drop_key, ]
  desc <- maihda_describe(y ~ (1 | gender:ethnicity:education:age_group),
                          data = sp, flag_stratum_n = 5)

  expect_equal(desc$overview$n_strata_expected, 36)
  expect_true(desc$overview$n_empty_strata >= 1)
  expect_equal(sum(desc$strata$empty), desc$overview$n_empty_strata)
  expect_equal(desc$overview$n_strata_observed +
                 desc$overview$n_empty_strata, 36)

  # Empty rows carry no stratum ID (IDs stay identical to a subsequent fit)
  # and zero counts.
  empt <- desc$strata[desc$strata$empty, , drop = FALSE]
  expect_true(all(is.na(empt$stratum)))
  expect_true(all(empt$n == 0))

  # Small strata are flagged, not dropped.
  obs <- desc$strata[!desc$strata$empty, , drop = FALSE]
  expect_identical(obs$small, obs$n <= 5)
  expect_true("small_strata" %in% desc$warnings$check)
  expect_true("empty_strata" %in% desc$warnings$check)

  # include_empty_strata = FALSE keeps the counts but drops the rows.
  desc2 <- maihda_describe(y ~ (1 | gender:ethnicity:education:age_group),
                           data = sp, include_empty_strata = FALSE)
  expect_false(any(desc2$strata$empty))
  expect_equal(desc2$overview$n_empty_strata, desc$overview$n_empty_strata)
})

test_that("a fitted maihda_model and a maihda_analysis are described post hoc", {
  d <- make_describe_data(8006)
  d$y[1:6] <- NA
  desc_f <- maihda_describe(y ~ age + (1 | gender:race), data = d)
  m <- fit_maihda(y ~ age + (1 | gender:race), data = d)

  desc_m <- maihda_describe(m)
  expect_identical(desc_m$source, "maihda_model")
  expect_identical(desc_m$engine, "lme4")
  obs_f <- desc_f$strata[!desc_f$strata$empty, , drop = FALSE]
  obs_m <- desc_m$strata[!desc_m$strata$empty, , drop = FALSE]
  expect_identical(obs_m$label, obs_f$label)
  expect_identical(obs_m$n, obs_f$n)
  expect_identical(obs_m$n_analytic, obs_f$n_analytic)
  expect_equal(desc_m$overview$n_analytic, desc_f$overview$n_analytic)
  expect_equal(desc_m$overview$n_missing_outcome, 6)

  # The bundled analysis from maihda() is accepted too.
  a <- suppressMessages(maihda(y ~ age + gender + race + (1 | gender:race),
                               data = d))
  desc_a <- maihda_describe(a)
  expect_identical(desc_a$source, "maihda_analysis")
  expect_identical(desc_a$strata$label, desc_m$strata$label)

  # A recoded binary fit surfaces the recoding.
  d$binf <- factor(ifelse(d$bin == 1, "yes", "no"))
  m2 <- suppressWarnings(suppressMessages(
    fit_maihda(binf ~ (1 | gender:race), data = d)))
  desc_m2 <- maihda_describe(m2)
  expect_identical(desc_m2$family, "binomial")
  expect_s3_class(desc_m2$response_recoding, "data.frame")

  # A fitted input already carries its data/family/context/weights.
  expect_error(maihda_describe(m, data = d), "do not also supply")
  expect_error(maihda_describe(m, family = "binomial"), "do not also supply")
})

test_that("sampling weights add weighted summaries and drop invalid weights", {
  d <- make_describe_data(8007)
  set.seed(99)
  d$w <- runif(nrow(d), 0.5, 2)
  d$w[1:3] <- 0
  d$w[4] <- NA

  desc <- maihda_describe(y ~ age + (1 | gender:race), data = d,
                          sampling_weights = "w")
  expect_equal(desc$overview$n_invalid_weights, 4)
  # The weighted engines drop missing/non-positive weights from the fit.
  expect_equal(desc$overview$n_analytic, nrow(d) - 4)
  expect_true(all(c("n_weighted", "outcome_mean_weighted") %in%
                    names(desc$strata)))
  valid <- is.finite(d$w) & d$w > 0
  expect_equal(desc$overview$sum_weights, sum(d$w[valid]), tolerance = 1e-12)
  expect_equal(desc$outcome_overall$outcome_mean_weighted,
               sum(d$w[valid] * d$y[valid]) / sum(d$w[valid]),
               tolerance = 1e-12)
  expect_true("n_weighted" %in% names(desc$dimensions))

  expect_error(maihda_describe(y ~ (1 | gender:race), data = d,
                               sampling_weights = "nope"),
               "not found")
})

test_that("(1 | stratum) input recovers make_strata() dims; a bare column degrades gracefully", {
  d <- make_describe_data(8008)

  s <- make_strata(d, vars = c("gender", "race"))
  desc <- maihda_describe(y ~ age + (1 | stratum), data = s$data)
  expect_identical(desc$strata_vars, c("gender", "race"))
  obs <- desc$strata[!desc$strata$empty, , drop = FALSE]
  expect_identical(obs$label, s$strata_info$label)
  expect_identical(obs$n, s$strata_info$n)

  # Hand-built stratum column without make_strata(): no dimensions, no
  # expected-strata Cartesian, but the per-stratum table still works.
  d2 <- d
  d2$stratum <- as.integer(interaction(d2$gender, d2$race, drop = TRUE))
  desc2 <- maihda_describe(y ~ age + (1 | stratum), data = d2)
  expect_null(desc2$dimensions)
  expect_null(desc2$strata_vars)
  expect_true(is.na(desc2$overview$n_strata_expected))
  expect_equal(desc2$overview$n_strata_observed, 6)
  expect_equal(sum(desc2$strata$n), nrow(d2))
  expect_output(print(desc2), "not recorded")

  # No stratum term at all is a clear error.
  expect_error(maihda_describe(y ~ age, data = d),
               "no intersectional stratum term")
})

test_that("numeric dimensions are flagged: autobin, linear, and ID-like", {
  d <- make_describe_data(8009, n = 200)

  # >10 unique values: make_strata() auto-bins (message) and the description
  # matches a subsequent fit on the same formula.
  expect_message(
    desc <- maihda_describe(y ~ (1 | gender:age), data = d),
    "auto-binned")
  expect_true("autobinned_dimension" %in% desc$warnings$check)
  m <- suppressMessages(fit_maihda(y ~ (1 | gender:age), data = d))
  obs <- desc$strata[!desc$strata$empty, , drop = FALSE]
  expect_identical(obs$label, m$strata_info$label)

  # A few distinct numeric values: kept raw, flagged as entering adjusted
  # models linearly.
  d$code <- sample(1:4, nrow(d), replace = TRUE)
  desc2 <- maihda_describe(y ~ (1 | gender:code), data = d)
  expect_true("linear_numeric_dimension" %in% desc2$warnings$check)

  # A high-cardinality dimension looks like an identifier.
  d$id <- seq_len(nrow(d))
  desc3 <- suppressMessages(
    maihda_describe(y ~ (1 | gender:id), data = d, autobin = FALSE))
  expect_true("id_like_dimension" %in% desc3$warnings$check)
})

test_that("print() renders the description and plot() returns ggplot objects", {
  d <- make_describe_data(8010)
  d$y[1:10] <- NA
  desc <- maihda_describe(y ~ age + (1 | gender:race), data = d,
                          flag_stratum_n = 30)

  expect_output(ret <- print(desc), "MAIHDA Sample Description")
  expect_output(print(desc), "Analytic sample")
  expect_output(print(desc), "Observed outcome")
  expect_output(print(desc), "full table in \\$strata")
  expect_identical(ret, desc)   # print() returns its input invisibly

  for (tp in c("stratum_size", "outcome", "missingness")) {
    p <- plot(desc, type = tp)
    expect_s3_class(p, "ggplot")
  }
  # Binomial outcome plot uses the category distribution.
  desc_bin <- maihda_describe(bin ~ (1 | gender:race), data = d)
  expect_s3_class(plot(desc_bin, type = "outcome"), "ggplot")
  expect_error(plot(desc, type = "nope"))
})

test_that("input validation", {
  d <- make_describe_data(8011, n = 60)
  expect_error(maihda_describe(42), "must be a model formula")
  expect_error(maihda_describe(y ~ (1 | gender:race)), "'data' is required")
  expect_error(maihda_describe(y ~ (1 | gender:race), data = d,
                               flag_stratum_n = -1),
               "non-negative")
  expect_error(maihda_describe(y ~ (1 | gender:race), data = d,
                               include_empty_strata = NA),
               "TRUE or FALSE")
  expect_error(maihda_describe(y ~ (1 | gender:race), data = d,
                               family = "nope"),
               "Unsupported family")
})
