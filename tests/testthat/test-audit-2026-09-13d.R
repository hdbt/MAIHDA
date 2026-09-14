# Audit 2026-09-13 (fourth pass): "A6. Ordinal predictions discard the fitted
# contrast matrix" -- CONFIRMED, and wider than reported.
#
# clmm has no predict() method and WeMix's needs the grouping re-resolved, so both
# engines rebuild the fixed-effect design by hand (maihda_clmm_linpred(),
# maihda_wemix_linpred(), and the clmm branch of the prediction-deviation panel). The
# rebuild called model.matrix() with neither the contrasts the fit used nor the
# levels it fitted. On this file's fixture, against the unmodified code:
#   * contrasts = list(f = "contr.sum") through fit_maihda()'s ... fitted, then every
#     prediction failed with "missing column(s): f1". A custom matrix whose column is
#     named like treatment coding's ("fhi") passed that name check and was silently
#     re-coded 0/1 instead of -1/+1: the fitted probability of the observed category
#     off by up to 0.180 against clmm's own.
#   * the same happened with no contrasts argument at all -- options(contrasts = )
#     in force at fit time (an error once restored), or a contrasts attribute on the
#     factor (silent, 0.054 in probability, under a misleading "contrasts dropped"
#     warning) -- and on the WeMix engine (0.293 on the link scale for the attribute),
#     where a binomial summary() itself failed once the options were restored and the
#     attribute route moved the design-weighted AUC silently (0.683 against 0.690).
#   * under the DEFAULT options, an ordered covariate with a level declared on the
#     data but absent from the fitted rows was mis-coded: the engines fit with unused
#     levels dropped (contr.poly on 3 levels), the rebuild kept them (contr.poly on
#     4). Silent: 0.145 in probability for clmm, 0.680 on the link scale for WeMix.
#   * clmm itself loses $contrasts when it drops an aliased column, so its own record
#     is not always there to reuse.
#
# FIX: the fit records its coding -- the fitted levels and every contrast as a
# numeric matrix at those levels -- while the options that shaped it are in force
# (maihda_engine_fixed_coding(), stored as $fixed_coding), and maihda_fixed_design()
# rebuilds every prediction design with it. Default-coded fits are bit-identical to
# the previous code (checked against the unmodified tree for ordinal and WeMix
# predictions -- unordered, ordered, character and logical covariates -- stratum
# tables, panel probabilities, the seeded proportional-odds p-value and maihda()'s
# PCV). A model saved without the record (every 0.2.1 fit) falls back to clmm's own
# record and the contrasts argument in its call; it cannot recover a coding set by
# options(contrasts = ) at fit time on an aliased clmm design or a WeMix fit, which
# errors when the column names differ and is silently wrong when they coincide. The
# oracles below are independent of the rebuild: clmm's own fitted probabilities and
# WeMix's own fitted values both come from the design the ENGINE built.

a6_data <- function(n = 600, seed = 2026) {
  set.seed(seed)
  d <- data.frame(gender = factor(sample(c("m", "f"), n, TRUE)),
                  race = factor(sample(c("a", "b", "c"), n, TRUE)),
                  f = factor(sample(c("lo", "hi"), n, TRUE), levels = c("lo", "hi")),
                  g = factor(sample(c("p", "q", "r"), n, TRUE)),
                  x = stats::rnorm(n))
  # An ordered covariate with a level ("top") declared but never observed.
  d$ed <- factor(sample(c("lo", "mid", "hi"), n, TRUE),
                 levels = c("lo", "mid", "hi", "top"), ordered = TRUE)
  cell <- interaction(d$gender, d$race, drop = TRUE)
  u <- stats::rnorm(nlevels(cell), sd = 0.6)[as.integer(cell)]
  eta <- 1.5 * (d$f == "hi") + 0.5 * d$x + 0.7 * (d$g == "q") - 0.4 * (d$g == "r") +
    0.8 * as.integer(d$ed) + u
  d$y <- factor(cut(eta + stats::rlogis(n), breaks = c(-Inf, 1, 2.3, 3.7, Inf),
                    labels = 1:4), ordered = TRUE)
  d$yc <- eta + stats::rnorm(n)
  # A stronger stratum effect for the binary outcome: with the one above, WeMix put
  # the stratum variance near 0 and its own summary() warned from cov2cor().
  d$yb <- as.integer(eta + 2 * u + stats::rlogis(n) > 2.3)
  d$w <- stats::runif(n, 0.5, 2)
  d
}

a6_fit <- function(...) suppressWarnings(suppressMessages(fit_maihda(...)))

# clmm's own fitted probability of each observed category (computed from the design
# clmm built) against the probability implied by the package's rebuilt location.
a6_clmm_gap <- function(fit, eta) {
  P <- maihda_ordinal_category_probs(eta, maihda_clmm_cutpoints(fit$model),
                                     fit$family$link)
  yi <- as.integer(fit$data$y)
  max(abs(P[cbind(seq_along(yi), yi)] - fit$model$fitted.values))
}

# WeMix's own fitted values, y - resid in row order (from the design WeMix built),
# against the package's rebuilt linear predictor.
a6_wemix_gap <- function(fit, y) {
  max(abs(maihda_wemix_linpred(fit, include_re = TRUE) -
            (fit$data[[y]] - unname(fit$model$resid))))
}

test_that("ordinal predictions reuse a sum contrast passed through fit_maihda()", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  d <- a6_data()
  fit <- a6_fit(y ~ f + g + x + (1 | gender:race), data = d, engine = "ordinal",
                family = "ordinal", contrasts = list(f = "contr.sum", g = "contr.sum"))
  # The fit really is sum-coded, so treatment columns cannot serve it.
  expect_setequal(names(fit$model$beta), c("f1", "g1", "g2", "x"))

  eta <- predict_maihda(fit, scale = "link")
  expect_lt(a6_clmm_gap(fit, eta), 1e-12)

  # The same location from the fitted design written out by hand.
  X <- stats::model.matrix(~ f + g + x, fit$data,
                           contrasts.arg = list(f = "contr.sum", g = "contr.sum"))
  re <- ordinal::ranef(fit$model)$stratum
  u <- stats::setNames(re[[1]], rownames(re))[as.character(fit$data$stratum)]
  expect_equal(eta, drop(X[, names(fit$model$beta)] %*% fit$model$beta) + unname(u),
               tolerance = 1e-12, ignore_attr = TRUE)

  # Response scale and a newdata batch agree with the in-sample prediction.
  sc <- predict_maihda(fit, scale = "response")
  expect_equal(sc, maihda_ordinal_eta_to_score(eta, maihda_clmm_cutpoints(fit$model),
                                               "logit"),
               tolerance = 1e-12, ignore_attr = TRUE)
  sub <- d$x > 0.5
  expect_equal(predict_maihda(fit, newdata = d[sub, ], scale = "link"), eta[sub],
               tolerance = 1e-12, ignore_attr = TRUE)
})

test_that("a custom contrast named like treatment coding is not silently re-coded", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  d <- a6_data()
  cm <- matrix(c(-1, 1), ncol = 1, dimnames = list(c("lo", "hi"), "hi"))
  fit <- a6_fit(y ~ f + g + x + (1 | gender:race), data = d, engine = "ordinal",
                family = "ordinal", contrasts = list(f = cm))

  # A column-name check cannot tell the two codings apart ...
  expect_true("fhi" %in% names(fit$model$beta))
  expect_true("fhi" %in% colnames(stats::model.matrix(~ f + g + x, d)))
  # ... and nothing errors: the predictions were simply wrong.
  eta <- predict_maihda(fit, scale = "link")
  expect_lt(a6_clmm_gap(fit, eta), 1e-12)
  expect_lt(max(abs(predict_maihda(fit, scale = "response") -
                      maihda_ordinal_eta_to_score(eta, maihda_clmm_cutpoints(fit$model),
                                                  "logit"))), 1e-12)

  # Check the VALUES of the rebuilt design: -1/+1, not treatment's 0/1.
  tt <- maihda_fitted_predict_terms(fit$formula, fit$data)$terms
  X <- maihda_fixed_design(tt, fit$data, maihda_object_fixed_coding(fit))$X
  expect_equal(unname(X[, "fhi"]), ifelse(fit$data$f == "hi", 1, -1))
  expect_equal(unname(X[, c("gq", "gr")]),
               unname(stats::model.matrix(~ g, fit$data)[, c("gq", "gr")]))
})

test_that("predictions do not depend on options(contrasts) at prediction time", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  d <- a6_data()
  op <- options(contrasts = c("contr.sum", "contr.poly"))
  on.exit(options(op), add = TRUE)
  fit_sum <- a6_fit(y ~ f + g + x + (1 | gender:race), data = d, engine = "ordinal",
                    family = "ordinal")
  options(op)
  expect_setequal(names(fit_sum$model$beta), c("f1", "g1", "g2", "x"))
  eta_sum <- predict_maihda(fit_sum, scale = "link")
  expect_lt(a6_clmm_gap(fit_sum, eta_sum), 1e-12)

  fit_trt <- a6_fit(y ~ f + g + x + (1 | gender:race), data = d, engine = "ordinal",
                    family = "ordinal")
  eta_trt <- predict_maihda(fit_trt, scale = "link")
  options(contrasts = c("contr.helmert", "contr.poly"))
  expect_identical(predict_maihda(fit_trt, scale = "link"), eta_trt)
  expect_identical(predict_maihda(fit_sum, scale = "link"), eta_sum)
  options(op)
})

test_that("a contrasts attribute on the factor is honoured, without a warning", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  d <- a6_data()
  cm3 <- matrix(c(-1, 1, 0, -1, -1, 2), ncol = 2,
                dimnames = list(c("p", "q", "r"), c("q", "r")))
  contrasts(d$g) <- cm3
  fit <- a6_fit(y ~ f + g + x + (1 | gender:race), data = d, engine = "ordinal",
                family = "ordinal")
  # Named like treatment coding again, so only the values give it away.
  expect_setequal(names(fit$model$beta), c("fhi", "gq", "gr", "x"))

  # The stored data keeps the attribute; re-levelling it used to warn "contrasts
  # dropped from factor g" while the treatment coding silently replaced it.
  expect_no_warning(eta <- predict_maihda(fit, scale = "link"))
  expect_lt(a6_clmm_gap(fit, eta), 1e-12)
  # newdata with and without the attribute predict alike.
  plain <- d
  attr(plain$g, "contrasts") <- NULL
  expect_equal(predict_maihda(fit, newdata = plain, scale = "link"), eta,
               tolerance = 1e-12, ignore_attr = TRUE)
})

test_that("a C() term keeps its fitted coding", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  d <- a6_data()
  fit <- a6_fit(y ~ C(g, contr.sum) + x + (1 | gender:race), data = d,
                engine = "ordinal", family = "ordinal")
  expect_no_warning(eta <- predict_maihda(fit, scale = "link"))
  expect_lt(a6_clmm_gap(fit, eta), 1e-12)
})

test_that("an ordered covariate with an unobserved declared level is coded as fitted", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  d <- a6_data()
  fit <- a6_fit(y ~ ed + x + (1 | gender:race), data = d, engine = "ordinal",
                family = "ordinal")
  # Default options: contr.poly on the three FITTED levels, while the stored data
  # still declares four.
  expect_identical(fit$model$xlevels$ed, c("lo", "mid", "hi"))
  expect_identical(levels(fit$data$ed), c("lo", "mid", "hi", "top"))
  expect_lt(a6_clmm_gap(fit, predict_maihda(fit, scale = "link")), 1e-12)

  # A newdata row AT that level is refused (lme4 and predict.lm refuse it too)
  # instead of being predicted from a coding the fit never had.
  nd <- d[1:3, ]
  nd$ed[1] <- "top"
  expect_error(predict_maihda(fit, newdata = nd, scale = "link"), "top")
})

test_that("clmm's lost $contrasts on an aliased design is recovered at fit time", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  d <- a6_data()
  d$x2 <- 2 * d$x
  fit <- a6_fit(y ~ g + x + x2 + (1 | gender:race), data = d, engine = "ordinal",
                family = "ordinal", contrasts = list(g = "contr.sum"))
  # clmm dropped x2 and, with it, the contrasts attribute of its design.
  expect_false("x2" %in% names(fit$model$beta))
  expect_null(fit$model$contrasts)
  expect_lt(a6_clmm_gap(fit, predict_maihda(fit, scale = "link")), 1e-12)
})

test_that("the other clmm prediction consumers use the fitted coding", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  d <- a6_data()
  fit <- a6_fit(y ~ f + g + x + (1 | gender:race), data = d, engine = "ordinal",
                family = "ordinal", contrasts = list(f = "contr.sum", g = "contr.sum"))
  yi <- as.integer(fit$data$y)
  pick <- cbind(seq_along(yi), yi)

  # Stratum table: its fixed part is the prediction-weighted mean, per stratum, of
  # the location built from the fitted design.
  sp <- maihda_stratum_predictions_ordinal(fit, summary(fit), scale = "link")
  X <- stats::model.matrix(~ f + g + x, fit$data,
                           contrasts.arg = list(f = "contr.sum", g = "contr.sum"))
  eta_fixed <- drop(X[, names(fit$model$beta)] %*% fit$model$beta)
  pw <- maihda_prediction_weights(fit)
  key <- as.character(fit$data$stratum)
  by_stratum <- vapply(split(seq_along(key), key),
                       function(i) stats::weighted.mean(eta_fixed[i], pw[i]), numeric(1))
  expect_identical(sp$stratum, names(by_stratum))
  expect_equal(sp$fixed_row, unname(by_stratum), tolerance = 1e-12)

  # Prediction-deviation panel: from the maihda_model's record, and from the bare
  # clmm's own $xlevels / $contrasts.
  pp <- maihda_prediction_panel_ordinal_probs(fit$model, fit$data,
                                              coding = maihda_object_fixed_coding(fit))
  expect_lt(max(abs(as.matrix(pp)[pick] - fit$model$fitted.values)), 1e-12)
  pp_bare <- maihda_prediction_panel_ordinal_probs(fit$model, fit$data)
  expect_lt(max(abs(as.matrix(pp_bare)[pick] - fit$model$fitted.values)), 1e-12)

  po <- maihda_proportional_odds_test(fit, n_sim = 9, seed = 7)
  expect_true(is.finite(po$p_value))
})

test_that("wemix predictions reuse contrasts that options() set at fit time", {
  skip_on_cran()
  skip_if_not_installed("WeMix")
  d <- a6_data()
  # WeMix's y - resid reproduces a default-coded rebuild, so it is a valid oracle.
  w_trt <- a6_fit(yc ~ g + x + (1 | gender:race), data = d, engine = "wemix",
                  sampling_weights = "w")
  expect_lt(a6_wemix_gap(w_trt, "yc"), 1e-10)

  op <- options(contrasts = c("contr.sum", "contr.poly"))
  on.exit(options(op), add = TRUE)
  w_opt <- a6_fit(yc ~ g + x + (1 | gender:race), data = d, engine = "wemix",
                  sampling_weights = "w")
  w_bin <- a6_fit(yb ~ g + x + (1 | gender:race), data = d, engine = "wemix",
                  sampling_weights = "w", family = "binomial")
  s_fit_time <- summary(w_bin)
  options(op)
  expect_setequal(names(w_opt$model$coef), c("(Intercept)", "g1", "g2", "x"))
  expect_lt(a6_wemix_gap(w_opt, "yc"), 1e-10)
  # A binomial fit's stratum standard errors and its AUC come from the same rebuild:
  # summary() failed once the fit-time options were restored. summary() drops an AUC
  # it cannot compute, so check there is one before comparing.
  s_now <- summary(w_bin)
  expect_identical(s_now$stratum_estimates$se, s_fit_time$stratum_estimates$se)
  expect_true(is.finite(s_fit_time$discriminatory_accuracy$auc))
  expect_identical(s_now$discriminatory_accuracy$auc, s_fit_time$discriminatory_accuracy$auc)
})

test_that("wemix predictions reuse the fitted levels and a factor's contrasts attribute", {
  skip_on_cran()
  skip_if_not_installed("WeMix")
  d <- a6_data()
  w_ord <- a6_fit(yc ~ ed + x + (1 | gender:race), data = d, engine = "wemix",
                  sampling_weights = "w")
  expect_setequal(names(w_ord$model$coef), c("(Intercept)", "ed.L", "ed.Q", "x"))
  expect_lt(a6_wemix_gap(w_ord, "yc"), 1e-10)

  contrasts(d$g) <- matrix(c(-1, 1, 0, -1, -1, 2), ncol = 2,
                           dimnames = list(c("p", "q", "r"), c("q", "r")))
  w_attr <- a6_fit(yc ~ g + x + (1 | gender:race), data = d, engine = "wemix",
                   sampling_weights = "w")
  expect_setequal(names(w_attr$model$coef), c("(Intercept)", "gq", "gr", "x"))
  expect_no_warning(gap <- a6_wemix_gap(w_attr, "yc"))
  expect_lt(gap, 1e-10)
})

test_that("the recorded coding holds numeric contrast matrices at the fitted levels", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  d <- a6_data()
  d$flag <- d$x > 0
  fit <- a6_fit(y ~ f + g + ed + flag + x + (1 | gender:race), data = d,
                engine = "ordinal", family = "ordinal", contrasts = list(g = "contr.sum"))
  rec <- fit$fixed_coding
  expect_identical(rec$xlev, fit$model$xlevels)
  expect_true(all(vapply(rec$contrasts, is.matrix, logical(1))))
  expect_setequal(names(rec$contrasts), c("f", "g", "ed", "flag"))
  expect_equal(unname(rec$contrasts$g), unname(stats::contr.sum(3)))
  expect_equal(unname(rec$contrasts$ed), unname(stats::contr.poly(3)))
  expect_identical(rownames(rec$contrasts$ed), c("lo", "mid", "hi"))
  # A logical is coded as factor(x, levels = c(FALSE, TRUE)).
  expect_identical(rownames(rec$contrasts$flag), c("FALSE", "TRUE"))
  expect_lt(a6_clmm_gap(fit, predict_maihda(fit, scale = "link")), 1e-12)

  # A model object without the record (fitted before it was kept) still predicts
  # from clmm's own record.
  legacy <- fit
  legacy$fixed_coding <- NULL
  expect_identical(predict_maihda(legacy, scale = "link"),
                   predict_maihda(fit, scale = "link"))
})

test_that("a saved model without the record recovers clmm's contrasts argument", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  # Every 0.2.1 fit lacks $fixed_coding. On an aliased design clmm has no $contrasts
  # either, but its call still holds the evaluated contrasts argument; without it this
  # same-named custom matrix would be re-coded silently (0.18 in probability).
  d <- a6_data()
  d$x2 <- 2 * d$x
  cm <- matrix(c(-1, 1), ncol = 1, dimnames = list(c("lo", "hi"), "hi"))
  fit <- a6_fit(y ~ f + g + x + x2 + (1 | gender:race), data = d, engine = "ordinal",
                family = "ordinal", contrasts = list(f = cm))
  expect_null(fit$model$contrasts)
  legacy <- fit
  legacy$fixed_coding <- NULL
  expect_lt(a6_clmm_gap(legacy, predict_maihda(legacy, scale = "link")), 1e-12)
})

test_that("a contrast function gone by prediction time does not matter", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  # The record holds the evaluated matrix, not the function's name, so a contrast
  # defined in one session (or by a package no longer attached) still predicts.
  skip_if(exists("a6_gone_contr", envir = globalenv(), inherits = FALSE),
          "the global environment already has an a6_gone_contr")
  assign("a6_gone_contr", function(n, contrasts = TRUE) {
    stats::contr.sum(n, contrasts) * 2
  }, envir = globalenv())
  on.exit(if (exists("a6_gone_contr", envir = globalenv(), inherits = FALSE)) {
    rm("a6_gone_contr", envir = globalenv())
  }, add = TRUE)

  d <- a6_data()
  fit <- a6_fit(y ~ g + x + (1 | gender:race), data = d, engine = "ordinal",
                family = "ordinal", contrasts = list(g = "a6_gone_contr"))
  rm("a6_gone_contr", envir = globalenv())
  expect_lt(a6_clmm_gap(fit, predict_maihda(fit, scale = "link")), 1e-12)
})

test_that("a contrast named by string is resolved as model.matrix() resolved it", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  # model.matrix() looks a contrast name up from the stats namespace outwards, so a
  # same-named object in the global environment -- where a user's own contrast
  # function lives -- does not replace stats' contr.sum in the fit. The record must
  # resolve the name the same way, or it stores a matrix the fit never used.
  skip_if(exists("contr.sum", envir = globalenv(), inherits = FALSE),
          "the global environment already has a contr.sum")
  assign("contr.sum", function(n, contrasts = TRUE, sparse = FALSE) {
    stats::contr.sum(n, contrasts) * 10
  }, envir = globalenv())
  on.exit(rm("contr.sum", envir = globalenv()), add = TRUE)

  d <- a6_data()
  fit <- a6_fit(y ~ g + x + (1 | gender:race), data = d, engine = "ordinal",
                family = "ordinal", contrasts = list(g = "contr.sum"))
  expect_equal(unname(fit$fixed_coding$contrasts$g), unname(stats::contr.sum(3)))
  expect_lt(a6_clmm_gap(fit, predict_maihda(fit, scale = "link")), 1e-12)
})

test_that("maihda_fixed_design() muffles only the superseded-contrasts warning", {
  skip_on_cran()
  d <- data.frame(g = factor(c("p", "q", "r", "q")), x = c(0.1, 0.4, -0.2, 1))
  coding <- list(xlev = list(g = c("p", "q", "r")),
                 contrasts = list(g = stats::contr.sum(3)))
  tt <- stats::delete.response(stats::terms(y ~ g + x))
  contrasts(d$g) <- stats::contr.helmert(3)
  expect_no_warning(des <- maihda_fixed_design(tt, d, coding))
  expect_equal(unname(des$X[, c("g1", "g2")]),
               unname(stats::contr.sum(3)[c(1, 2, 3, 2), ]))
  # Any other model.frame() warning still reaches the caller (here "variable 'g' is
  # not a factor", after which model.matrix() refuses to apply a contrast to it).
  num <- data.frame(g = c(1, 2, 3, 2), x = d$x)
  expect_warning(tryCatch(maihda_fixed_design(tt, num, coding),
                          error = function(e) NULL))
})

test_that("the superseded-contrasts warning is muffled in a translated locale too", {
  skip_on_cran()
  skip_if_not(exists("Sys.setLanguage", envir = baseenv()), "R < 4.2")
  skip_if_not(capabilities("NLS"))
  skip_if(Sys.getlocale("LC_CTYPE") %in% c("C", "POSIX"), "C locale")
  old <- Sys.getenv("LANGUAGE", unset = NA)
  on.exit({
    if (is.na(old)) Sys.unsetenv("LANGUAGE") else Sys.setenv(LANGUAGE = old)
    invisible(bindtextdomain(NULL))
  }, add = TRUE)
  suppressWarnings(Sys.setLanguage("it"))
  # model.frame.default() translates this warning; an English literal would miss it.
  skip_if(identical(gettextf("contrasts dropped from factor %s", "g", domain = "R-stats"),
                    "contrasts dropped from factor g"),
          "no Italian stats message catalog")

  d <- data.frame(g = factor(c("p", "q", "r", "q")), x = c(0.1, 0.4, -0.2, 1))
  contrasts(d$g) <- stats::contr.helmert(3)
  coding <- list(xlev = list(g = c("p", "q", "r")),
                 contrasts = list(g = stats::contr.sum(3)))
  tt <- stats::delete.response(stats::terms(y ~ g + x))
  expect_no_warning(maihda_fixed_design(tt, d, coding))
})
