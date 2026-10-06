# Audit 2026-10-05 (lme4): predict_maihda(newdata =) on an engine = "lme4" fit, two
# defects in how the grouping columns reach predict.merMod().
#
# 1. As reported: under allow_new_levels = TRUE a grouping column that is NA in every
#    newdata row ERRORED, "Invalid grouping factor specification, stratum" -- the
#    one-row prediction for a row missing its stratum dimension, which the package's
#    own message invites ("Pass allow_new_levels = TRUE to predict those rows with the
#    stratum random effect set to zero"), a supplied NA stratum, and an all-NA context
#    or longitudinal id. lme4's mkBlist() stops on a grouping factor with no value at
#    all; an NA escapes only as a factor level, which lme4 makes only for a factor
#    training column met by character or factor newdata -- hence the type boundary
#    measured (an NA_character_ context worked, a logical NA or any stratum NA did not).
#
# 2. Found while fixing it, and worse: lme4 2.0.1's levelfun() returns the fitted
#    random effects in the FITTED level order when newdata holds no new level, while
#    the random-effects design follows the newdata factor's levels. The stratum this
#    package rebuilds from the dimension columns is character, so it sorts as text
#    ("10" before "2"), and rows were SILENTLY given other strata's random effects --
#    at the default allow_new_levels = FALSE, 0.846 off fitted() on a 20-stratum fit
#    for rows of strata 2, 8, 10 and 15. Plain lme4 shows it with no package code (an
#    integer-trained grouping column predicted through a character copy). With a new
#    level present lme4 reorders by the newdata levels, which is why unseen-stratum
#    rows were right. The prediction-deviation panel's lme4 route had the same defect
#    for supplied data carrying a character stratum (2.21 off fitted()); per the user
#    it is fixed here too.
#
# Fixed by handing lme4 each pure grouping column as a factor in the fitted level
# order, with an NA under allow_new_levels = TRUE made a level no fitted level carries
# (maihda_lme4_grouping_newdata()). WeMix and ordinal already returned the zero effect
# for a missing stratum; brms does once its own NA-level fix is merged. An independent
# check found the guard missing an NA held as an explicit factor level (addNA()),
# which then became the new level even without allow_new_levels; it is pinned below.

l405_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

# The prediction, or NA when it errors -- so a defect fails each assertion on its own
# rather than aborting the block.
l405_pred <- function(...) tryCatch(unname(predict_maihda(...)), error = function(e) NA_real_)

l405_data <- function() {
  set.seed(1)
  J <- 20; n_per <- 15
  a <- rep(LETTERS[1:J], each = n_per)
  u <- stats::rnorm(J)
  x <- stats::rnorm(J * n_per)
  y <- 2 + 0.5 * x + u[match(a, LETTERS[1:J])] + stats::rnorm(J * n_per)
  data.frame(y = y, x = x, a = a, stringsAsFactors = FALSE)
}

test_that("a row missing its stratum dimension predicts at a zero stratum effect (lme4)", {
  skip_on_cran()
  m <- fit_maihda(y ~ x + (1 | a), data = l405_data())
  fe <- unname(lme4::fixef(m$model))
  p <- function(nd, ...) l405_pred(m, newdata = nd, scale = "link",
                                   allow_new_levels = TRUE, ...)

  # The reported call, and its spellings: every row NA.
  expect_equal(p(data.frame(x = 0, a = NA)), fe[1], tolerance = 1e-10)
  expect_equal(p(data.frame(x = 0, a = NA_character_)), fe[1], tolerance = 1e-10)
  expect_equal(p(data.frame(x = c(0, 1), a = c(NA, NA))), c(fe[1], sum(fe)),
               tolerance = 1e-10)
  # A supplied missing stratum, whatever its type.
  for (s in list(NA, NA_character_, NA_integer_)) {
    expect_equal(p(data.frame(x = 1, stratum = s)), sum(fe), tolerance = 1e-10)
  }
  # The response scale follows (Gaussian identity: the same value).
  expect_equal(l405_pred(m, newdata = data.frame(x = 0, a = NA), allow_new_levels = TRUE),
               fe[1], tolerance = 1e-10)

  # Unchanged: a missing dimension is still refused without allow_new_levels, with the
  # package's own message.
  expect_error(predict_maihda(m, newdata = data.frame(x = 0, a = NA)),
               "allow_new_levels = TRUE", fixed = TRUE)
})

test_that("seen strata keep their own random effects in rebuilt newdata (lme4)", {
  skip_on_cran()
  d <- l405_data()
  m <- fit_maihda(y ~ x + (1 | a), data = d)
  rows <- match(c("B", "H", "J", "O"), d$a)          # strata 2, 8, 10 and 15
  stratum <- m$model@frame$stratum[rows]
  expect_identical(as.integer(stratum), c(2L, 8L, 10L, 15L))
  ref <- unname(stats::fitted(m$model)[rows])
  nd <- d[rows, c("x", "a")]

  # The default call: the stratum is rebuilt from the dimension column, as character.
  prep <- MAIHDA:::maihda_prepare_prediction_data(m, nd)
  expect_type(prep$stratum, "character")
  expect_equal(l405_pred(m, newdata = nd, scale = "link"), ref, tolerance = 1e-10)
  expect_equal(l405_pred(m, newdata = nd, scale = "link", allow_new_levels = TRUE), ref,
               tolerance = 1e-10)
  # Beside a row missing its dimension (no new level then, so lme4 kept the fitted
  # order -- the case that went wrong) and beside an unseen combination (a new level,
  # which lme4 reorders by, so it was right before).
  with_na <- rbind(nd, data.frame(x = 0, a = NA))
  expect_equal(l405_pred(m, newdata = with_na, scale = "link", allow_new_levels = TRUE),
               c(ref, unname(lme4::fixef(m$model))[1]), tolerance = 1e-10)
  with_new <- rbind(nd, data.frame(x = 0, a = "new"))
  expect_equal(l405_pred(m, newdata = with_new, scale = "link",
                         allow_new_levels = TRUE)[1:4], ref, tolerance = 1e-10)
  # A supplied stratum: character went wrong, integer was always right.
  expect_equal(l405_pred(m, newdata = data.frame(x = nd$x, stratum = as.character(stratum)),
                         scale = "link"), ref, tolerance = 1e-10)
  expect_equal(l405_pred(m, newdata = data.frame(x = nd$x, stratum = stratum),
                         scale = "link"), ref, tolerance = 1e-10)
})

test_that("a binomial fit keeps its strata on the response scale (lme4)", {
  skip_on_cran()
  set.seed(5)
  n <- 720
  d <- data.frame(g1 = sample(c("a", "b", "c", "d"), n, TRUE),
                  g2 = sample(c("u", "v", "w"), n, TRUE), x = stats::rnorm(n),
                  stringsAsFactors = FALSE)
  eff <- stats::rnorm(12, sd = 1)
  d$yb <- stats::rbinom(n, 1, stats::plogis(0.2 + 0.4 * d$x +
                                            eff[as.integer(factor(paste(d$g1, d$g2)))]))
  m <- l405_quiet(fit_maihda(yb ~ x + (1 | g1:g2), data = d, family = "binomial"))
  fr <- m$model@frame
  rows <- c(match(2L, fr$stratum), match(11L, fr$stratum), match(9L, fr$stratum),
            match(12L, fr$stratum))
  ref <- unname(stats::fitted(m$model)[rows])
  expect_equal(l405_pred(m, newdata = d[rows, c("x", "g1", "g2")]), ref, tolerance = 1e-10)
})

test_that("missing context and longitudinal levels predict at a zero effect (lme4)", {
  skip_on_cran()
  d <- l405_data()
  set.seed(2)
  d$school <- sample(paste0("s", 1:8), nrow(d), replace = TRUE)
  d$y <- d$y + stats::setNames(stats::rnorm(8, sd = 1.2), paste0("s", 1:8))[d$school]
  m <- fit_maihda(y ~ x + (1 | a), data = d, context = "school")
  fe <- unname(lme4::fixef(m$model))[1]
  re <- lme4::ranef(m$model)
  u_a <- re$stratum[as.character(m$strata_info$stratum[m$strata_info$a == "A"]), 1]
  v_1 <- re$school["s1", 1]
  p <- function(nd) l405_pred(m, newdata = nd, scale = "link", allow_new_levels = TRUE)
  expect_equal(p(data.frame(x = 0, a = "A", school = NA)), fe + u_a, tolerance = 1e-10)
  expect_equal(p(data.frame(x = 0, a = "A", school = NA_character_)), fe + u_a,
               tolerance = 1e-10)
  expect_equal(p(data.frame(x = 0, a = NA, school = "s1")), fe + v_1, tolerance = 1e-10)
  expect_equal(p(data.frame(x = 0, a = NA, school = NA)), fe, tolerance = 1e-10)
  expect_equal(p(data.frame(x = 0, a = "A", school = c(NA, "s1"))),
               c(fe + u_a, fe + u_a + v_1), tolerance = 1e-10)
  # Unchanged: without allow_new_levels an NA context is refused by lme4, and a
  # context column absent from newdata is an error either way.
  expect_error(predict_maihda(m, newdata = data.frame(x = 0, a = "A", school = NA)))
  expect_error(predict_maihda(m, newdata = data.frame(x = 0, a = "A"),
                              allow_new_levels = TRUE), "school")
  # An NA held as a factor level is refused by lme4 as before, not turned into the
  # internal new level.
  na_level <- data.frame(x = 0, a = "A", school = addNA(factor(NA_character_)))
  err <- tryCatch(predict_maihda(m, newdata = na_level), error = conditionMessage)
  expect_type(err, "character")
  expect_false(grepl(".maihda_missing_level", err, fixed = TRUE))

  ids <- unique(maihda_long_data$id)[1:120]
  dl <- maihda_long_data[maihda_long_data$id %in% ids, ]
  ml <- l405_quiet(fit_maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
                              data = dl, id = "id", time = "wave"))
  fl <- lme4::fixef(ml$model); rl <- lme4::ranef(ml$model)
  tcol <- names(fl)[2]
  row <- dl[2, ]
  prep <- MAIHDA:::maihda_prepare_prediction_data(ml, row)
  s <- as.character(prep$stratum)
  tv <- prep[[tcol]]
  no_id <- unname(fl[1] + fl[2] * tv + rl$stratum[s, 1] + rl$stratum[s, 2] * tv)
  with_id <- no_id + rl$id[as.character(row$id), 1] + rl$id[as.character(row$id), 2] * tv
  one <- row; one$id <- NA
  expect_equal(l405_pred(ml, newdata = one, scale = "link", allow_new_levels = TRUE),
               no_id, tolerance = 1e-10)
  two <- rbind(row, row); two$id[2] <- NA
  expect_equal(l405_pred(ml, newdata = two, scale = "link", allow_new_levels = TRUE),
               c(with_id, no_id), tolerance = 1e-10)
})

test_that("the prediction-deviation panel keeps the strata of supplied data (lme4)", {
  skip_on_cran()
  m <- fit_maihda(y ~ x + (1 | a), data = l405_data())
  fr <- m$data
  chr <- fr
  chr$stratum <- as.character(chr$stratum)
  ref <- unname(stats::fitted(m$model))
  panel_fit <- function(data) tryCatch(
    MAIHDA:::maihda_prediction_panel_fitted(m$model, data, "gaussian",
                                            fitted_data = FALSE)$fit,
    error = function(e) NA_real_)
  expect_equal(panel_fit(fr), ref, tolerance = 1e-10)       # integer: always right
  expect_equal(panel_fit(chr), ref, tolerance = 1e-10)
  # Through the public function: the stratum summaries do not depend on the type.
  panel_data <- function(data) {
    p <- plot_prediction_deviation_panels(m, data = data, type = "gaussian")
    g <- if (inherits(p, "ggplot")) p else p[[1]]
    out <- g$data[, c("stratum", "fitted", "mean_fitted")]
    out$stratum <- as.integer(as.character(out$stratum))
    out[order(out$stratum), ]
  }
  expect_equal(panel_data(chr), panel_data(fr), tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("maihda_lme4_grouping_newdata hands lme4 factors in the fitted level order", {
  skip_on_cran()
  set.seed(9)
  d <- data.frame(g = sample(1:12, 240, TRUE), h = sample(c("p", "q"), 240, TRUE),
                  x = stats::rnorm(240), stringsAsFactors = FALSE)
  d$y <- d$x + stats::rnorm(12)[d$g] + stats::rnorm(240)
  fit <- l405_quiet(lme4::lmer(y ~ x + h + (1 | g) + (0 + x | h), data = d))
  f <- MAIHDA:::maihda_lme4_grouping_newdata
  fitted_levels <- levels(lme4::getME(fit, "flist")$g)

  out <- f(fit, data.frame(x = 0, h = "p", g = c("10", "2")))
  expect_s3_class(out$g, "factor")
  expect_identical(levels(out$g), fitted_levels)
  expect_identical(as.character(out$g), c("10", "2"))
  # h is a covariate as well as a grouping factor, so it is left alone.
  expect_identical(out$h, c("p", "p"))
  # An NA: kept for lme4 to refuse without allow_new, a new level with it -- one no
  # fitted level carries, even when a fitted level holds the default name.
  expect_true(anyNA(f(fit, data.frame(x = 0, h = "p", g = c(NA, 3L)))$g))
  # ... and so is an NA held as an explicit factor level, which anyNA() on the
  # factor itself does not see.
  held <- f(fit, data.frame(x = 0, h = "p", g = addNA(factor(c(NA, 3L)))))$g
  expect_true(anyNA(as.character(held)))
  na_new <- f(fit, data.frame(x = 0, h = "p", g = c(NA, 3L)), allow_new = TRUE)$g
  expect_false(anyNA(na_new))
  expect_false(as.character(na_new[1]) %in% fitted_levels)
  expect_identical(as.character(na_new[2]), "3")
  # A column absent from newdata stays absent.
  expect_false("g" %in% names(f(fit, data.frame(x = 0, h = "p"), allow_new = TRUE)))

  d2 <- d
  d2$g <- ifelse(d2$g == 1L, ".maihda_missing_level", as.character(d2$g))
  fit2 <- l405_quiet(lme4::lmer(y ~ x + (1 | g), data = d2))
  taken <- f(fit2, data.frame(x = 0, g = NA), allow_new = TRUE)$g
  expect_false(as.character(taken) %in% levels(lme4::getME(fit2, "flist")$g))
  expect_equal(unname(stats::predict(fit2, newdata = data.frame(x = 0, g = taken),
                                     allow.new.levels = TRUE)),
               unname(lme4::fixef(fit2)[1]), tolerance = 1e-10)
})
