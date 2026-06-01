test_that("maihda_country_data is available and has correct structure", {
  data(maihda_country_data, envir = environment())

  expect_true(is.data.frame(maihda_country_data))
  expect_equal(nrow(maihda_country_data), 3600)
  expect_equal(ncol(maihda_country_data), 7)

  expected_cols <- c("country", "gender", "ses", "escs", "math", "reading", "low_math")
  expect_equal(names(maihda_country_data), expected_cols)

  expect_true(is.factor(maihda_country_data$country))
  expect_true(is.factor(maihda_country_data$gender))
  expect_true(is.factor(maihda_country_data$ses))
  expect_true(is.numeric(maihda_country_data$escs))
  expect_true(is.numeric(maihda_country_data$math))
  expect_true(is.numeric(maihda_country_data$reading))
  expect_true(is.factor(maihda_country_data$low_math))

  # Balanced: 6 countries x 600, strata = gender x ses (6 strata)
  expect_equal(length(levels(maihda_country_data$country)), 6)
  expect_true(all(table(maihda_country_data$country) == 600))
  expect_equal(levels(maihda_country_data$ses), c("Low", "Medium", "High"))
  expect_equal(levels(maihda_country_data$gender), c("female", "male"))

  # No missing values in the variables used for the showcase model
  expect_false(anyNA(maihda_country_data[, c("country", "gender", "ses", "math")]))
})

test_that("maihda_country_data showcases differing intersectional VPC across countries", {
  data(maihda_country_data, envir = environment())

  cmp <- compare_maihda_groups(
    math ~ 1 + (1 | gender:ses),
    data = maihda_country_data,
    group = "country"
  )

  # One estimable row per country, each with the full 6 (gender x ses) strata
  expect_equal(nrow(cmp), 6)
  expect_true(all(cmp$status == "ok"))
  expect_true(all(cmp$n_strata == 6))
  expect_true(all(is.finite(cmp$vpc) & cmp$vpc >= 0 & cmp$vpc <= 1))

  # The point of the dataset: VPC genuinely varies across countries.
  expect_gt(diff(range(cmp$vpc)), 0.02)

  # Deterministic orderings for this fixed sample (robust to minor numeric drift):
  # Germany shows the largest gender-by-SES inequality, the UK the smallest.
  expect_equal(as.character(cmp$group[which.max(cmp$vpc)]), "Germany")
  expect_equal(as.character(cmp$group[which.min(cmp$vpc)]), "United Kingdom")
})
