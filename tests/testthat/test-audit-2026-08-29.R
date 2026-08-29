# Audit 2026-08-29 -- a stratum dimension written in TRANSFORMED form in the fixed part
# (factor(a), scale(a), I(a), poly(a, 2), splines) was not recognised as that dimension's
# additive main effect: only a bare column name was. The transform therefore survived into
# the derived NULL model (deflating its between-stratum variance and invalidating the PCV),
# the adjusted model added the bare name alongside it, the crossed-dimensions model entered
# the dimension as a fixed effect AND a random intercept, and `factor(a) * b` walked past
# the fixed-cell-interaction guard entirely. The decomposition entry points now reject a
# transformed dimension up front.

maihda_tf_data <- function(seed = 1, n = 600) {
  set.seed(seed)
  a <- sample(1:3, n, TRUE)
  b <- sample(c("x", "y", "z"), n, TRUE)
  st <- paste(a, b, sep = "_")
  st_eff <- stats::setNames(stats::rnorm(9, 0, 0.4), unique(st))[st]
  data.frame(
    y = 5 + c(0, 3, 1)[a] + c(0, 1, 2)[match(b, c("x", "y", "z"))] +
      st_eff + stats::rnorm(n, 0, 1),
    a = a, b = b,
    g = sample(c("g1", "g2"), n, TRUE),
    age = stats::rnorm(n),
    stringsAsFactors = FALSE
  )
}

test_that("maihda_transformed_dimension_terms() flags only transformed dimensions", {
  sv <- c("a", "b")
  expect_equal(maihda_transformed_dimension_terms(y ~ factor(a) + b, sv), "factor(a)")
  expect_equal(maihda_transformed_dimension_terms(y ~ scale(a) + b, sv), "scale(a)")
  expect_equal(maihda_transformed_dimension_terms(y ~ I(a > 2) + b, sv), "I(a > 2)")
  expect_equal(maihda_transformed_dimension_terms(y ~ poly(a, 2) + b, sv), "poly(a, 2)")
  expect_equal(maihda_transformed_dimension_terms(y ~ splines::ns(a, 2) + b, sv),
               "splines::ns(a, 2)")
  # picked up inside an interaction, and reported once
  expect_equal(maihda_transformed_dimension_terms(y ~ factor(a) * b, sv), "factor(a)")
  # bar terms are not part of the fixed part
  expect_equal(maihda_transformed_dimension_terms(y ~ a + b + (1 | a:b), sv),
               character(0))

  # negative controls: bare dimensions, transformed COVARIATES, offsets, the response,
  # and a backtick-quoted non-syntactic dimension name (parses to a name, not a call)
  expect_equal(maihda_transformed_dimension_terms(y ~ a + b, sv), character(0))
  expect_equal(maihda_transformed_dimension_terms(y ~ a + b + log(age), sv),
               character(0))
  expect_equal(maihda_transformed_dimension_terms(y ~ a + b + a:log(age), sv),
               character(0))
  expect_equal(
    maihda_transformed_dimension_terms(y ~ a + offset(log(expo)), c("a", "expo")),
    character(0))
  expect_equal(maihda_transformed_dimension_terms(cbind(s, f) ~ a + b, c("s", "a")),
               character(0))
  expect_equal(
    maihda_transformed_dimension_terms(stats::as.formula("y ~ `a var` + b"),
                                       c("a var", "b")),
    character(0))
  expect_equal(maihda_transformed_dimension_terms(y ~ 1, sv), character(0))
  # an auto-binned dimension is matched through its .maihda_dim_* reconstruction too
  expect_equal(
    maihda_transformed_dimension_terms(y ~ factor(.maihda_dim_a) + b, sv,
                                       dim_terms = c(".maihda_dim_a", "b")),
    "factor(.maihda_dim_a)")
})

test_that("the message names the transformed dimension, not the whole stratum set", {
  expect_equal(
    maihda_transformed_dimension_vars("factor(a)", c("a", "b", "c")), "a")
  expect_error(
    maihda_check_no_transformed_dimension(y ~ factor(a) + b + c, c("a", "b", "c")),
    "dimension(s) a.", fixed = TRUE)
  expect_error(
    maihda_check_no_transformed_dimension(y ~ factor(a) + b + c, c("a", "b", "c")),
    "data$a <- factor(data$a)", fixed = TRUE)
  expect_silent(maihda_check_no_transformed_dimension(y ~ a + b, c("a", "b")))
})

test_that("the fixed-cell-interaction guard sees through a transformed dimension", {
  # `factor(a):b` involves only stratum dimensions, so it is a fixed cell-means
  # interaction and must be flagged even though `factor(a)` is not a bare dimension name.
  expect_equal(maihda_dimension_interaction_terms(y ~ factor(a) * b, c("a", "b")),
               "factor(a):b")
  expect_error(
    maihda_check_no_dimension_interaction(y ~ factor(a) * b, c("a", "b")),
    "interaction term(s) factor(a):b", fixed = TRUE)
  # a dimension-by-COVARIATE interaction stays legitimate in either spelling
  expect_equal(maihda_dimension_interaction_terms(y ~ a * age, c("a", "b")),
               character(0))
  expect_equal(maihda_dimension_interaction_terms(y ~ factor(a) * age, c("a", "b")),
               character(0))
})

test_that("maihda() rejects a transformed stratum dimension in every decomposition mode", {
  skip_on_cran()
  d <- maihda_tf_data()

  expect_error(
    suppressMessages(maihda(y ~ factor(a) + b + (1 | a:b), data = d, autobin = FALSE)),
    "transformed appearance")
  expect_error(
    suppressMessages(maihda(y ~ scale(a) + b + (1 | a:b), data = d, autobin = FALSE)),
    "transformed appearance")
  expect_error(
    suppressMessages(maihda(y ~ factor(a) + b + (1 | a:b), data = d, autobin = FALSE,
                            decomposition = "crossed-dimensions")),
    "transformed appearance")
  # `factor(a) * b` previously evaded the interaction guard; the transform check now
  # catches it first, with the more specific message.
  expect_error(
    suppressMessages(maihda(y ~ factor(a) * b + (1 | a:b), data = d, autobin = FALSE)),
    "transformed appearance")
  expect_error(
    suppressMessages(compare_maihda_groups(y ~ factor(a) + b + (1 | a:b), data = d,
                                           group = "g", autobin = FALSE)),
    "transformed appearance")
})

test_that("maihda() still accepts bare dimensions and transformed covariates", {
  skip_on_cran()
  d <- maihda_tf_data()

  a_bare <- suppressWarnings(suppressMessages(
    maihda(y ~ a + b + (1 | a:b), data = d, autobin = FALSE, interactions = FALSE)))
  # the null model carries NO dimension term -- the property the bug violated
  null_terms <- attr(stats::terms(maihda_nobars(a_bare$formula)), "term.labels")
  expect_false(any(c("a", "b") %in% null_terms))
  expect_true(is.finite(a_bare$pcv$pcv))

  # a transformed COVARIATE, and a dimension-by-covariate interaction, are untouched
  expect_no_error(suppressWarnings(suppressMessages(
    maihda(y ~ scale(age) + a + b + (1 | a:b), data = d, autobin = FALSE,
           interactions = FALSE))))
  expect_no_error(suppressWarnings(suppressMessages(
    maihda(y ~ a + b + a:age + (1 | a:b), data = d, autobin = FALSE,
           interactions = FALSE))))
  # and the documented remedy -- transform the COLUMN, write the bare name -- works
  d2 <- d
  d2$a <- factor(d2$a)
  a_fac <- suppressWarnings(suppressMessages(
    maihda(y ~ a + b + (1 | a:b), data = d2, autobin = FALSE, interactions = FALSE)))
  expect_true(a_fac$pcv$pcv > 0.9)
})

test_that("maihda_interactions() reports the transform rather than 'looks like a null model'", {
  skip_on_cran()
  d <- maihda_tf_data()
  dd <- make_strata(d, vars = c("a", "b"), autobin = FALSE)$data
  m_tf <- suppressMessages(fit_maihda(y ~ factor(a) + b + (1 | stratum), data = dd))
  expect_warning(maihda_interactions(m_tf), "transformed appearance")
  m_ok <- suppressMessages(fit_maihda(y ~ a + b + (1 | stratum), data = dd))
  expect_no_warning(maihda_interactions(m_ok))
})
