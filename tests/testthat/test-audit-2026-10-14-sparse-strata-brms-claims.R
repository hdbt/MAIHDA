# Audit 2026-10-14 (sparse-strata-brms-claims): several places said that brms handles,
# regularises or calibrates sparse strata and sparse growth data. Measured against a
# known truth (30 to 1,000 data sets per design):
# - in the design of maihda_sparse_data (36 strata, 240 rows, cells of median 6) the
#   two-model VPC is close to the truth with either engine (Gaussian: lme4 0.26 and
#   brms 0.27 for 0.26; binary: lme4 0.28 and brms 0.34 for 0.31; lme4 never singular;
#   brms intervals cover the truth in 96 to 99%);
# - in a crossed-dimensions fit of that design lme4 is singular in about 70% of data
#   sets, mostly through a dimension's variance, and understates the interaction share
#   (0.29 to 0.34 for 0.42); brms's interval for that share covers the truth in 23%
#   (Gaussian) and 70% (binary) of data sets under the default prior and in 81% / 90%
#   under the vignette's priors, so it is not calibrated;
# - with tiny cells (60 rows over 50 strata) lme4 puts the stratum variance at zero
#   in 41 to 56% of data sets and brms's posterior median VPC is about 0.05 to 0.07
#   whether the true VPC is 0, 0.02 or 0.05;
# - in a growth design with few waves (8 strata of 5 persons, 3 waves) lme4 is
#   singular in about 58% of fits, through a correlation at +-1, while its slope
#   variance is about unbiased; brms's posterior median is about as close, and the
#   upper limits of its intervals are a median of ten times the truth.
# The help pages, maihda_describe()'s note and the vignettes now say what each engine
# does. The blocks below pin the runtime note and, from the package tree, the reworded
# text.

test_that("the empty-strata note of maihda_describe() names no engine", {
  data(maihda_sparse_data, envir = environment())
  # Remove the smallest observed cell so the strata space has a hole.
  key <- with(maihda_sparse_data, paste(gender, ethnicity, education, age_group, sep = "\r"))
  sp <- maihda_sparse_data[key != names(sort(table(key)))[1], ]
  desc <- maihda_describe(y ~ (1 | gender:ethnicity:education:age_group), data = sp)
  msg <- desc$warnings$message[desc$warnings$check == "empty_strata"]
  expect_identical(msg, paste0(
    "1 of 36 expected strata are empty (no observations): the strata space is sparse; ",
    "small strata are shrunk toward the mean and their effects are imprecise with ",
    "either engine."))
  expect_false(grepl("brms", msg, fixed = TRUE))
})

test_that("the help pages and vignettes no longer say brms repairs sparse cells", {
  skip_on_cran()
  root <- testthat::test_path("..", "..")
  # Present when the suite runs from the package tree; R CMD check copies only tests/.
  skip_if_not(dir.exists(file.path(root, "man")), "the package tree is not available")
  flat <- function(...) {
    p <- file.path(root, ...)
    if (!file.exists(p)) return(NA_character_)
    gsub("[[:space:]]+", " ", paste(readLines(p, warn = FALSE), collapse = " "))
  }
  sparse_rd <- flat("man", "maihda_sparse_data.Rd")
  for (old in c("calibrated", "collapses to a singular fit", "purely a small-cell artifact",
                "honest uncertainty")) {
    expect_false(grepl(old, sparse_rd, fixed = TRUE), info = old)
  }
  expect_true(grepl("can be singular", sparse_rd, fixed = TRUE))
  describe_rd <- flat("man", "maihda_describe.Rd")
  expect_false(grepl("exactly what motivates", describe_rd, fixed = TRUE))
  expect_true(grepl("imprecise with either engine", describe_rd, fixed = TRUE))

  planning <- flat("vignettes", "planning_a_study.Rmd")
  if (!is.na(planning)) {
    expect_false(grepl("regularise the variance off the boundary", planning, fixed = TRUE))
    expect_false(grepl("which sparse cells call for", planning, fixed = TRUE))
    expect_true(grepl("`engine = \"brms\"` does not repair this", planning, fixed = TRUE))
  }
  longitudinal <- flat("vignettes", "longitudinal.Rmd")
  if (!is.na(longitudinal)) {
    expect_false(grepl("the brms engine handles this better", longitudinal, fixed = TRUE))
    expect_true(grepl("reports an interval instead", longitudinal, fixed = TRUE))
  }
  sparse_vignette <- flat("vignettes", "bayesian_sparse_maihda.Rmd")
  if (!is.na(sparse_vignette)) {
    expect_true(grepl("read with the priors it rests on", sparse_vignette, fixed = TRUE))
  }
})
