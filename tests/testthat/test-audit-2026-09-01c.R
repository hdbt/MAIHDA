# Audit 2026-09-01c -- a fixed interaction between a COVARIATE and a stratum dimension
# (age * gender) was left alone by the fixed-cell-interaction guard as "a legitimate
# covariate adjustment". It is not: the null model is derived by removing the
# dimensions' BARE main effects, so `age:gender` survived into the null, which then
# already adjusted for gender.
#
# With gender absent marginally the surviving term is a gender contrast SCALED BY age,
# so how much gender the null absorbs depends on the arbitrary origin of age:
# re-centring age moved the PCV by 2.25 percentage points on maihda_health_data
# (51.57% -> 49.32%) while the adjusted fits stayed identical to 1e-9. For a
# CATEGORICAL covariate R gives the dimension full dummy coding inside the term and the
# null's fixed part spans the dimension's main effect exactly (Education * Gender:
# 47.60% vs the correct 41.14%). In crossed-dimensions mode the dimension entered as a
# fixed effect AND a random intercept at once, competing for the same contrast, so its
# additive variance was not identified as intended -- asserted structurally below, not
# as a variance change: a (1 | gender) variance on a two-level dimension comes from two
# groups and moves either way by draw. In longitudinal mode a user-written
# `gender * wave` put the dim:time term the adjusted growth model is supposed to ADD
# into the null as well, corrupting PCV_slope.
#
# maihda() / compare_maihda_groups() now reject any fixed interaction involving a
# stratum dimension, in every decomposition mode.

maihda_dxc_data <- function(seed = 901, n = 900) {
  set.seed(seed)
  gender <- sample(c("F", "M"), n, TRUE)
  race <- sample(c("X", "Y", "Z"), n, TRUE)
  edu <- sample(c("lo", "mid", "hi"), n, TRUE)
  sk <- interaction(gender, race, drop = TRUE)
  age <- stats::rnorm(n, 50, 10)
  d <- data.frame(
    y = 1 + 0.05 * age + 0.4 * (gender == "M") +
      stats::rnorm(nlevels(sk), sd = 0.7)[sk] + stats::rnorm(n, sd = 0.9),
    gender = gender, race = race, edu = edu, age = age,
    g = sample(c("g1", "g2"), n, TRUE),
    stringsAsFactors = FALSE
  )
  # appended last so the columns above keep their draws; only used to show the guard
  # runs before the family branches split.
  d$bin <- stats::rbinom(n, 1, stats::plogis(-0.3 + 0.02 * (age - 50) +
                                               stats::rnorm(nlevels(sk), sd = 0.6)[sk]))
  d
}

test_that("the interaction detector separates dimension-only from dimension-covariate", {
  sv <- c("gender", "race")
  # scope = "all" (the default) keeps the historical cell-means semantics.
  expect_equal(maihda_dimension_interaction_terms(y ~ age * gender + race, sv),
               character(0))
  expect_equal(maihda_dimension_interaction_terms(y ~ age + gender * race, sv),
               "gender:race")
  # scope = "any" adds the dimension-by-covariate case, in either variable order.
  expect_equal(
    maihda_dimension_interaction_terms(y ~ age * gender + race, sv, scope = "any"),
    "age:gender")
  expect_equal(
    maihda_dimension_interaction_terms(y ~ gender:age + race, sv, scope = "any"),
    "gender:age")
  # A three-way term involving a dimension is caught by "any" but not by "all".
  expect_equal(maihda_dimension_interaction_terms(y ~ age * gender * race, sv),
               "gender:race")
  expect_true("age:gender:race" %in%
                maihda_dimension_interaction_terms(y ~ age * gender * race, sv,
                                                   scope = "any"))
  # Negative controls: no dimension in the term, additive form, random-effect bars.
  expect_equal(
    maihda_dimension_interaction_terms(y ~ age * edu + gender + race, sv,
                                       scope = "any"),
    character(0))
  expect_equal(
    maihda_dimension_interaction_terms(y ~ age + gender + race + (1 | gender:race), sv,
                                       scope = "any"),
    character(0))
  # A single dimension is never a decomposition, so nothing is flagged.
  expect_equal(
    maihda_dimension_interaction_terms(y ~ age * gender, "gender", scope = "any"),
    character(0))
  # An auto-binned numeric dimension is matched through its .maihda_dim_* term too.
  expect_equal(
    maihda_dimension_interaction_terms(y ~ age * .maihda_dim_a + b, c("a", "b"),
                                       dim_terms = c(".maihda_dim_a", "b"),
                                       scope = "any"),
    "age:.maihda_dim_a")
  # ... and a non-syntactic dimension name, which terms() backtick-quotes.
  expect_equal(
    maihda_dimension_interaction_terms(
      stats::as.formula("y ~ age * `gender var` + race"),
      c("gender var", "race"), scope = "any"),
    "age:`gender var`")
})

test_that("the guard's message names the offending term, the dimension and the fix", {
  expect_error(
    maihda_check_no_dimension_interaction(y ~ age * gender + race,
                                          c("gender", "race")),
    "age:gender", fixed = TRUE)
  expect_error(
    maihda_check_no_dimension_interaction(y ~ age * gender + race,
                                          c("gender", "race")),
    "Write the additive form (`age + gender`", fixed = TRUE)
  # The caller's name is used, so a direct compare_maihda_groups() call says so.
  expect_error(
    maihda_check_no_dimension_interaction(y ~ age * gender + race,
                                          c("gender", "race"),
                                          fn = "compare_maihda_groups"),
    "^compare_maihda_groups\\(\\)")
  # A dimension-by-dimension term keeps the more specific cell-means message even when
  # a dimension-by-covariate term is present in the same formula.
  expect_error(
    maihda_check_no_dimension_interaction(y ~ age * gender + gender * race,
                                          c("gender", "race")),
    "duplicates the intersectional stratum random intercept", fixed = TRUE)
  # A dimension-by-TIME term gets the longitudinal remedy, not "fit it separately".
  expect_error(
    maihda_check_no_dimension_interaction(y ~ gender * wave + race,
                                          c("gender", "race"), time = "wave"),
    "adds the dimension-by-time interactions", fixed = TRUE)
  expect_error(
    maihda_check_no_dimension_interaction(y ~ age * gender + race,
                                          c("gender", "race"), time = "wave"),
    "separate model outside the decomposition", fixed = TRUE)
  # No false positive on the additive form.
  expect_silent(
    maihda_check_no_dimension_interaction(y ~ age + gender + race,
                                          c("gender", "race")))
})

test_that("maihda() rejects a dimension-covariate interaction in every mode", {
  skip_on_cran()
  d <- maihda_dxc_data()

  # two-model, numeric covariate: the centering-dependent PCV
  expect_error(
    suppressMessages(suppressWarnings(
      maihda(y ~ age * gender + race + (1 | gender:race), data = d))),
    "between a covariate and the stratum-defining dimension")
  # two-model, categorical covariate: the null spans the dimension's main effect
  expect_error(
    suppressMessages(suppressWarnings(
      maihda(y ~ edu * gender + race + (1 | gender:race), data = d))),
    "between a covariate and the stratum-defining dimension")
  # the explicit colon spelling is the same term
  expect_error(
    suppressMessages(suppressWarnings(
      maihda(y ~ age + gender + race + age:gender + (1 | gender:race), data = d))),
    "between a covariate and the stratum-defining dimension")
  # crossed-dimensions: the dimension would be a fixed effect AND a random intercept
  expect_error(
    suppressMessages(suppressWarnings(
      maihda(y ~ age * gender + race + (1 | gender:race), data = d,
             decomposition = "crossed-dimensions"))),
    "between a covariate and the stratum-defining dimension")
  # group=, and a direct compare_maihda_groups() call
  expect_error(
    suppressMessages(suppressWarnings(
      maihda(y ~ age * gender + race + (1 | gender:race), data = d, group = "g"))),
    "between a covariate and the stratum-defining dimension")
  expect_error(
    suppressWarnings(compare_maihda_groups(
      y ~ age * gender + race + (1 | gender:race), data = d, group = "g",
      min_group_n = 10)),
    "between a covariate and the stratum-defining dimension")
  # binomial / poisson / cumulative reach the same guard (it runs before the mode and
  # family branches split), so one non-Gaussian family stands for the rest.
  expect_error(
    suppressMessages(suppressWarnings(
      maihda(bin ~ age * gender + race + (1 | gender:race), data = d,
             family = "binomial"))),
    "between a covariate and the stratum-defining dimension")
})

test_that("an AUTO-BINNED numeric dimension is caught through its raw name", {
  # Distinct branch: strata_vars holds the raw column while dim_terms holds the
  # reconstructed .maihda_dim_* tertile factor, and the user writes the raw name.
  # The detector must match either, and the message must name the raw column as the
  # dimension and the other variable as the covariate.
  skip_on_cran()
  d <- maihda_dxc_data()
  set.seed(903)
  d$ageq <- round(stats::runif(nrow(d), 20, 80))

  expect_error(
    suppressMessages(suppressWarnings(
      maihda(y ~ age * ageq + race + (1 | ageq:race), data = d))),
    "between a covariate and the stratum-defining dimension")
  # ... and the additive form on the same auto-binned dimension still runs.
  a <- suppressMessages(suppressWarnings(
    maihda(y ~ age + ageq + race + (1 | ageq:race), data = d)))
  expect_true(is.finite(a$pcv$pcv))
  # Two auto-binned/raw DIMENSIONS interacting keep the cell-means message instead --
  # it is the more specific diagnosis and must not be shadowed by the new one.
  expect_error(
    suppressMessages(suppressWarnings(
      maihda(y ~ ageq * race + (1 | ageq:race), data = d))),
    "between the stratum-defining dimensions")
})

test_that("maihda() rejects a user-written dim:time term in longitudinal mode", {
  skip_on_cran()
  dl <- maihda_long_data
  expect_error(
    suppressMessages(suppressWarnings(
      maihda(wellbeing ~ gender * wave + ethnicity + (1 | gender:ethnicity), dl,
             id = "id", time = "wave", decomposition = "longitudinal"))),
    "adds the dimension-by-time interactions", fixed = TRUE)
  # The additive longitudinal form still runs: maihda() adds dim + dim:time itself, and
  # those internally-built terms must not trip the guard on the next call.
  a <- suppressMessages(suppressWarnings(
    maihda(wellbeing ~ (1 | gender:ethnicity), dl, id = "id", time = "wave",
           decomposition = "longitudinal")))
  expect_true(any(grepl("wave", attr(stats::terms(maihda_nobars(a$adjusted_formula)),
                                     "term.labels"), fixed = TRUE)))
  expect_false(any(grepl("gender", attr(stats::terms(maihda_nobars(a$formula)),
                                        "term.labels"), fixed = TRUE)))
})

test_that("the additive forms the guard permits are unaffected", {
  skip_on_cran()
  d <- maihda_dxc_data()

  # covariate x covariate involves no dimension
  a <- suppressMessages(suppressWarnings(
    maihda(y ~ age * edu + gender + race + (1 | gender:race), data = d)))
  expect_true(is.finite(a$pcv$pcv))
  # a transformed COVARIATE is fine; only a transformed DIMENSION is rejected
  b <- suppressMessages(suppressWarnings(
    maihda(y ~ poly(age, 2) + gender + race + (1 | gender:race), data = d)))
  expect_true(is.finite(b$pcv$pcv))
  # the plain additive adjusted model is bit-identical to the pre-fix result
  cc <- suppressMessages(suppressWarnings(
    maihda(y ~ age + gender + race + (1 | gender:race), data = d)))
  expect_true(is.finite(cc$pcv$pcv))
  # ... and the derived null carries no dimension term at all
  null_terms <- attr(stats::terms(maihda_nobars(cc$formula)), "term.labels")
  expect_equal(null_terms, "age")
})

test_that("a null model derived without dimension terms is centering-invariant", {
  # The property the guard now protects: once no fixed term involves a dimension, the
  # PCV is invariant to the covariate's origin. Fitted directly (no maihda()) so the
  # assertion is about the statistics, not the guard.
  skip_on_cran()
  skip_if_not_installed("lme4")
  d <- maihda_dxc_data()
  d$stratum <- interaction(d$gender, d$race, drop = TRUE)
  d$age_c <- d$age - 50

  pcv_of <- function(nullf, adjf) {
    n <- lme4::lmer(nullf, d, REML = TRUE)
    a <- lme4::lmer(adjf, d, REML = TRUE)
    t0 <- unname(unlist(lme4::VarCorr(n)))
    (t0 - unname(unlist(lme4::VarCorr(a)))) / t0
  }
  p_raw <- pcv_of(y ~ age + (1 | stratum), y ~ age + gender + race + (1 | stratum))
  p_ctr <- pcv_of(y ~ age_c + (1 | stratum), y ~ age_c + gender + race + (1 | stratum))
  expect_equal(p_raw, p_ctr, tolerance = 1e-6)

  # The rejected specification is NOT invariant -- this is the defect, reproduced.
  bad_raw <- pcv_of(y ~ age + age:gender + (1 | stratum),
                    y ~ age * gender + race + (1 | stratum))
  bad_ctr <- pcv_of(y ~ age_c + age_c:gender + (1 | stratum),
                    y ~ age_c * gender + race + (1 | stratum))
  expect_gt(abs(bad_raw - bad_ctr), 1e-3)
})

test_that("the crossed-dimensions formula would carry the dimension twice", {
  # The crossed-dimensions half of the defect is a STRUCTURAL one, asserted on the
  # formula rather than on a fitted variance: the pre-fix workflow stripped only the
  # bare dimension main effects before adding (1 | dim), so `age:gender` stayed in the
  # fixed part and gender entered as a fixed effect AND a random intercept, competing
  # for the same contrast. Its additive variance is then not identified as intended.
  # Deliberately NOT asserted as a variance change: a (1 | gender) variance on a
  # two-level dimension is estimated from two groups, so it moves either way by draw.
  f <- y ~ age * gender + race + (1 | gender:race)
  stripped <- stats::update(f, . ~ . - gender - race)
  fixed <- attr(stats::terms(maihda_nobars(stripped)), "term.labels")
  expect_true("age:gender" %in% fixed)
  cc <- maihda_cross_classified_formula(
    stripped, c("gender", "race"), NULL,
    data.frame(y = 1, age = 1, gender = "F", race = "X"))
  cc_fixed <- attr(stats::terms(maihda_nobars(cc$formula)), "term.labels")
  expect_true("age:gender" %in% cc_fixed)
  expect_match(deparse1(cc$formula), "(1 | gender)", fixed = TRUE)
  # ... which is what the guard now refuses to build.
  expect_error(
    maihda_check_no_dimension_interaction(f, c("gender", "race")),
    "between a covariate and the stratum-defining dimension")
})
