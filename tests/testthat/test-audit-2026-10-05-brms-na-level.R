# Audit 2026-10-05: predict_maihda(type = "individual", allow_new_levels = TRUE) on
# an engine = "brms" fit predicted a row with a MISSING grouping value with a SAMPLED
# random effect instead of the documented zero effect. Under allow_new_levels = TRUE,
# maihda_prepare_prediction_data() keeps the stratum NA for a row missing a
# stratum-defining dimension (and passes a supplied NA stratum through), and
# maihda_brms_individual_prediction() kept every term whose row level was NA ("An NA
# grouping value is likewise left to the normal path"). brms then received that term
# with allow_new_levels = TRUE, read the NA as a new level, and its default
# sample_new_levels = "uncertainty" drew a new effect -- a different one on every call.
# On the reported model (y ~ x + (1 | a), 20 strata of 15, one chain of 800) the
# zero-effect value is 2.170913 and one row with a = NA returned 2.220580, 2.193850 and
# 2.139049 on three calls; seen and unseen rows were exact. lme4's allow.new.levels
# gives an NA level a zero effect, and so do maihda_wemix_linpred() and
# maihda_clmm_linpred(); the error message that sends a user here promises "the stratum
# random effect set to zero". (lme4 errors instead when an all-NA newdata column
# reaches it with no factor level -- the package's own stratum column, or a logical
# NA -- a separate defect.) The exemption came with the per-term rewrite of the
# 2026-07-25 audit (e4ba090): the stratum-only helper it replaced zeroed an NA stratum.
# The same held for an NA context level and an NA longitudinal id. Fixed by dropping a
# term for a row whose level is NA, exactly as for an unseen level; the caller's
# re_formula is still only ever narrowed. With allow_new_levels = FALSE brms refuses an
# NA level itself, as lme4 does, and that path is unchanged.
# The independent check found the same sampling for a grouping column ABSENT from
# newdata, which the code had kept on the false premise that brms would raise an
# error; under allow_new_levels = TRUE brms fills the column with NA and draws. Per the
# user it is now refused, as lme4 refuses it, unless the caller's re_formula leaves
# that term out.

# A brms-shaped model stub. maihda_brms_individual_prediction() reads only $formula
# and $data, and maihda_prepare_prediction_data() the strata table.
na05_stub <- function(formula = y ~ x + (1 | stratum)) {
  list(
    formula = formula,
    data = data.frame(y = c(0, 1, 0, 1), x = c(0, 1, 0, 1), wave = c(0, 1, 0, 1),
                      stratum = c("1", "2", "1", "2"), a = c("A", "B", "A", "B"),
                      school = c("s1", "s2", "s1", "s2"),
                      id = c("p1", "p2", "p3", "p4"), stringsAsFactors = FALSE),
    strata_info = data.frame(stratum = c("1", "2"), label = c("A", "B"),
                             a = c("A", "B"), stringsAsFactors = FALSE),
    strata_vars = "a", strata_sep = " x ")
}

# "FULL" when no scope is written (every term), "NONE" for re_formula = NA (fixed
# effects only), otherwise the grouping variables of the terms kept.
na05_scope_label <- function(dots) {
  nm <- intersect(c("re.form", "re_formula"), names(dots))
  if (length(nm) == 0L || is.null(dots[[nm[1L]]])) {
    return("FULL")
  }
  v <- dots[[nm[1L]]]
  if (is.logical(v) && length(v) == 1L && is.na(v)) {
    return("NONE")
  }
  if (!inherits(v, "formula")) {
    return(paste("UNREAD:", format(v)))
  }
  paste(sort(unlist(lapply(reformulas::findbars(v), function(b) all.vars(b[[3]])))),
        collapse = "+")
}

# The scope maihda_brms_individual_prediction() hands brms for each row of `nd` (by
# position), plus the dots of every block, with the block predictor mocked -- no Stan.
na05_scope <- function(object, nd, dots = list(allow_new_levels = TRUE),
                       allow = TRUE) {
  nd$.row <- seq_len(nrow(nd))
  scope <- rep(NA_character_, nrow(nd))
  blocks <- list()
  testthat::local_mocked_bindings(
    maihda_brms_predict_rows = function(object, nd, scale, dots) {
      scope[nd$.row] <<- na05_scope_label(dots)
      blocks[[length(blocks) + 1L]] <<- dots
      rep(0, nrow(nd))
    })
  MAIHDA:::maihda_brms_individual_prediction(object, nd, "link", allow, dots)
  list(scope = scope, dots = blocks)
}

na05_prepare <- function(object, nd) {
  MAIHDA:::maihda_prepare_prediction_data(object, nd, allow_new_levels = TRUE)
}

test_that("a row missing its stratum dimension gets no stratum term (brms, no Stan)", {
  obj <- na05_stub()
  prep <- na05_prepare(obj, data.frame(x = 0, a = c(NA, "A", "new"),
                                       stringsAsFactors = FALSE))
  # The premise: allow_new_levels = TRUE leaves the stratum of that row missing,
  # while an unseen combination gets an id no fitted stratum carries.
  expect_true(is.na(prep$stratum[1]))
  expect_identical(prep$stratum[2], "1")
  expect_false(prep$stratum[3] %in% c("1", "2", NA))

  # Missing, seen, unseen: the missing row is predicted like the unseen one.
  expect_identical(na05_scope(obj, prep)$scope, c("NONE", "FULL", "NONE"))

  # The reported call: one row whose only dimension is missing.
  one <- na05_prepare(obj, data.frame(x = 0, a = NA))
  expect_identical(na05_scope(obj, one)$scope, "NONE")
})

test_that("a supplied missing stratum gets no stratum term (brms, no Stan)", {
  obj <- na05_stub()
  prep <- na05_prepare(obj, data.frame(x = 0, stratum = c(NA, "1"),
                                       stringsAsFactors = FALSE))
  expect_true(is.na(prep$stratum[1]))
  expect_identical(na05_scope(obj, prep)$scope, c("NONE", "FULL"))
})

test_that("missing context and longitudinal levels are dropped term by term (no Stan)", {
  obj <- na05_stub(y ~ x + (1 | stratum) + (1 | school))
  prep <- na05_prepare(obj, data.frame(
    x = 0, a = c("A", NA, NA, "A", "A"), school = c(NA, "s1", NA, "s1", "zz"),
    stringsAsFactors = FALSE))
  # (seen, NA) keeps the stratum; (NA, seen) keeps the school; (NA, NA) keeps
  # nothing; (seen, seen) keeps both; (seen, unseen) keeps the stratum, as before.
  expect_identical(na05_scope(obj, prep)$scope,
                   c("stratum", "school", "NONE", "FULL", "stratum"))

  lng <- na05_stub(y ~ wave + (wave | id) + (wave | stratum))
  nd <- data.frame(wave = 1, id = c("p1", NA, NA), stratum = c("1", "1", NA),
                   stringsAsFactors = FALSE)
  expect_identical(na05_scope(lng, nd)$scope, c("FULL", "stratum", "NONE"))
})

test_that("a level is still kept when the training levels cannot be read (no Stan)", {
  # The conservative rule stands for a level that is present: with no training
  # levels to compare against, "zz" might be a fitted one, so its term is kept. A
  # missing level is not a fitted one whatever the training data hold.
  obj <- na05_stub(y ~ x + (1 | stratum) + (1 | school))
  obj$data$school <- NULL
  nd <- data.frame(x = 0, stratum = "1", school = c("zz", NA),
                   stringsAsFactors = FALSE)
  expect_identical(na05_scope(obj, nd)$scope, c("FULL", "stratum"))
})

test_that("a grouping column absent from newdata is refused, not sampled (no Stan)", {
  # brms fills an absent grouping column with NA and, under allow_new_levels = TRUE,
  # draws an effect for every row -- on a real context fit, 1.80, 1.75, 1.69 and 1.76
  # over four calls against a zero-school value of 1.77 -- where lme4 refuses
  # ("object 'school' not found"), as brms itself does without allow_new_levels. Per
  # the user it is refused as lme4 does.
  obj <- na05_stub(y ~ x + (1 | stratum) + (1 | school))
  nd <- data.frame(x = 0, stratum = c("1", NA), stringsAsFactors = FALSE)
  expect_error(na05_scope(obj, nd),
               "missing the grouping variable(s) of the model's random effects: school",
               fixed = TRUE)
  expect_error(na05_scope(obj, nd, dots = list(allow_new_levels = TRUE,
                                               re.form = ~ (1 | school))),
               "random effects: school", fixed = TRUE)
  lng <- na05_stub(y ~ wave + (wave | id) + (wave | stratum))
  expect_error(na05_scope(lng, data.frame(wave = 1, stratum = "1")),
               "random effects: id", fixed = TRUE)

  # A caller's re_formula that leaves the term out needs no column.
  expect_identical(na05_scope(obj, nd, dots = list(allow_new_levels = TRUE,
                                                   re_formula = ~ (1 | stratum)))$scope,
                   c("stratum", "NONE"))
  expect_identical(na05_scope(obj, nd, dots = list(allow_new_levels = TRUE,
                                                   re_formula = NA))$scope,
                   c("NONE", "NONE"))
  # Without allow_new_levels the rows reach brms untouched, and brms refuses them.
  expect_identical(na05_scope(obj, nd, dots = list(allow_new_levels = FALSE),
                              allow = FALSE)$scope, c("FULL", "FULL"))
})

test_that("without allow_new_levels an NA level still reaches brms's refusal (no Stan)", {
  # brms rejects an NA grouping level itself when allow_new_levels = FALSE ("Levels
  # 'NA' of grouping factor ... cannot be found"), as lme4 does; the fix must not
  # turn that refusal into a silent zero-effect prediction.
  obj <- na05_stub(y ~ x + (1 | stratum) + (1 | school))
  nd <- data.frame(x = 0, stratum = "1", school = c(NA, "s1"), stringsAsFactors = FALSE)
  r <- na05_scope(obj, nd, dots = list(allow_new_levels = FALSE), allow = FALSE)
  expect_identical(r$scope, c("FULL", "FULL"))
  expect_length(r$dots, 1L)
  expect_false(r$dots[[1]]$allow_new_levels)
})

test_that("a missing level only narrows the caller's re_formula (no Stan)", {
  obj <- na05_stub(y ~ x + (1 | stratum) + (1 | school))
  nd <- data.frame(x = 0, stratum = NA_character_, school = "s1",
                   stringsAsFactors = FALSE)
  sc <- function(...) na05_scope(obj, nd, dots = list(allow_new_levels = TRUE, ...))

  # Fixed effects only stays fixed effects only; the seen school is not added back.
  expect_identical(sc(re_formula = NA)$scope, "NONE")
  # A request for the seen school alone reaches brms exactly as written.
  school_only <- ~ (1 | school)
  r <- sc(re_formula = school_only)
  expect_identical(r$dots[[1]]$re_formula, school_only)
  # The missing stratum leaves a request that names it ...
  expect_identical(sc(re_formula = ~ (1 | stratum))$scope, "NONE")
  expect_identical(sc(re_formula = ~ (1 | stratum) + (1 | school))$scope, "school")
  # ... written back under the caller's own spelling, never both.
  r <- sc(re.form = ~ (1 | stratum) + (1 | school))
  expect_false("re_formula" %in% names(r$dots[[1]]))
  expect_identical(r$scope, "school")
  # A value outside brms's {NULL, NA, formula} contract stays the caller's.
  expect_identical(sc(re_formula = "nonsense")$dots[[1]]$re_formula, "nonsense")
})

test_that("lme4, WeMix and ordinal predict a missing dimension at a zero effect", {
  # The behaviour brms is brought into line with: on the likelihood engines a row
  # missing its stratum dimension takes the fixed-effects-only value, while a seen
  # row in the same newdata keeps its stratum effect.
  skip_on_cran()
  set.seed(1)
  J <- 20; n_per <- 15
  a <- rep(LETTERS[1:J], each = n_per)
  u <- stats::rnorm(J)
  x <- stats::rnorm(J * n_per)
  y <- 2 + 0.5 * x + u[match(a, LETTERS[1:J])] + stats::rnorm(J * n_per)
  d <- data.frame(y = y, x = x, a = a, w = stats::runif(J * n_per, 0.5, 2),
                  stringsAsFactors = FALSE)
  d$yo <- cut(d$y, stats::quantile(d$y, 0:4 / 4), include.lowest = TRUE,
              labels = c("lo", "mid", "hi", "top"), ordered_result = TRUE)
  nd <- data.frame(x = 1, a = c(NA, "A"), stringsAsFactors = FALSE)

  m <- fit_maihda(y ~ x + (1 | a), data = d)
  p <- predict_maihda(m, newdata = nd, scale = "link", allow_new_levels = TRUE)
  fe <- lme4::fixef(m$model)
  re <- lme4::ranef(m$model)$stratum
  s_a <- m$strata_info$stratum[m$strata_info$a == "A"]
  expect_equal(unname(p), unname(sum(fe) + c(0, re[s_a, 1])), tolerance = 1e-10)

  mw <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ x + (1 | a), data = d, sampling_weights = "w")))
  pw <- predict_maihda(mw, newdata = nd, scale = "link", allow_new_levels = TRUE)
  prep_w <- MAIHDA:::maihda_prepare_prediction_data(mw, nd, allow_new_levels = TRUE)
  fixed_w <- MAIHDA:::maihda_wemix_linpred(mw, newdata = prep_w, include_re = FALSE)
  expect_equal(unname(pw[1]), unname(fixed_w[1]), tolerance = 1e-10)
  expect_gt(abs(pw[2] - fixed_w[2]), 0.1)

  mo <- suppressWarnings(suppressMessages(
    fit_maihda(yo ~ x + (1 | a), data = d, engine = "ordinal")))
  po <- predict_maihda(mo, newdata = nd, scale = "link", allow_new_levels = TRUE)
  prep_o <- MAIHDA:::maihda_prepare_prediction_data(mo, nd, allow_new_levels = TRUE)
  fixed_o <- MAIHDA:::maihda_clmm_linpred(mo, newdata = prep_o, include_re = FALSE)
  expect_equal(unname(po[1]), unname(fixed_o[1]), tolerance = 1e-10)
  expect_gt(abs(po[2] - fixed_o[2]), 0.1)
})

test_that("brms predicts missing stratum and context levels at a zero effect (Stan)", {
  skip_on_cran()
  skip_if(Sys.getenv("MAIHDA_TEST_BRMS") != "true",
          "brms Stan tests are opt-in; set MAIHDA_TEST_BRMS=true to run them")
  skip_if_not_installed("brms")
  set.seed(1)
  J <- 20; n_per <- 15
  a <- rep(LETTERS[1:J], each = n_per)
  u <- stats::rnorm(J)
  x <- stats::rnorm(J * n_per)
  school <- sample(paste0("s", 1:8), J * n_per, replace = TRUE)
  v <- stats::setNames(stats::rnorm(8, sd = 1.2), paste0("s", 1:8))
  y <- 2 + 0.5 * x + u[match(a, LETTERS[1:J])] + v[school] + stats::rnorm(J * n_per)
  d <- data.frame(y = y, x = x, a = a, school = school, stringsAsFactors = FALSE)
  m <- suppressWarnings(suppressMessages(fit_maihda(
    y ~ x + (1 | a), data = d, engine = "brms", context = "school",
    chains = 1, iter = 800, refresh = 0, seed = 1)))

  # Rows: both missing, stratum dimension missing, school missing, both seen.
  nd <- data.frame(x = 0, a = c(NA, NA, "A", "A"), school = c(NA, "s1", NA, "s1"),
                   stringsAsFactors = FALSE)
  p1 <- predict_maihda(m, newdata = nd, scale = "link", allow_new_levels = TRUE)
  p2 <- predict_maihda(m, newdata = nd, scale = "link", allow_new_levels = TRUE)
  # A sampled effect differed from call to call; a zero one cannot.
  expect_identical(p1, p2)

  # Both missing: the posterior mean of the intercept, at x = 0.
  expect_equal(unname(p1[1]), brms::fixef(m$model)["Intercept", "Estimate"],
               tolerance = 1e-8)
  # The others: brms's own prediction with exactly the seen terms, read through
  # brms rather than the package's helpers.
  s_a <- m$strata_info$stratum[m$strata_info$a == "A"]
  at <- data.frame(x = 0, stratum = s_a, school = "s1", stringsAsFactors = FALSE)
  lp <- function(re) mean(brms::posterior_linpred(m$model, newdata = at,
                                                  re_formula = re))
  expect_equal(unname(p1[2]), lp(~ (1 | school)), tolerance = 1e-8)
  expect_equal(unname(p1[3]), lp(~ (1 | stratum)), tolerance = 1e-8)
  expect_equal(unname(p1[4]), lp(NULL), tolerance = 1e-8)
  # The kept effects are real ones, so the oracles above are not all one number.
  expect_gt(abs(p1[2] - p1[1]), 0.1)
  expect_gt(abs(p1[3] - p1[1]), 0.1)

  # The response scale follows (Gaussian: the same values, through fitted()).
  r1 <- predict_maihda(m, newdata = nd[1:3, ], scale = "response",
                       allow_new_levels = TRUE)
  ft <- function(re) stats::fitted(m$model, newdata = at, re_formula = re,
                                   summary = TRUE)[, "Estimate"]
  expect_equal(unname(r1), unname(c(ft(NA), ft(~ (1 | school)), ft(~ (1 | stratum)))),
               tolerance = 1e-8)

  # Without the school column the prediction is refused rather than sampled, unless
  # the caller's re_formula leaves the school out.
  no_school <- data.frame(x = 0, a = "A", stringsAsFactors = FALSE)
  expect_error(predict_maihda(m, newdata = no_school, allow_new_levels = TRUE),
               "random effects: school", fixed = TRUE)
  expect_equal(unname(predict_maihda(m, newdata = no_school, scale = "link",
                                     allow_new_levels = TRUE,
                                     re_formula = ~ (1 | stratum))),
               lp(~ (1 | stratum)), tolerance = 1e-8)
})
