# Audit 2026-10-14 (single finding 2026-10-14-lme4-covariate-scale): fit_maihda() on
# the lme4 engine recorded sound glmer fits as not converged when a covariate is far
# from zero or in large units -- CONFIRMED, wider than reported.
#
# lme4 takes the gradient and Hessian it checks a glmer fit with by central differences
# with a step of 1e-4 in the coefficients as the formula has them, and draws absolute
# lines on them (0.002 on a scaled gradient, 1e6 on an eigenvalue of the Hessian). On
# the finding's designs (2 x 3 to 2 x 3 x 4 strata, 180 to 1,200 rows, one covariate in
# four exact affine spellings) 21 of 24 birth-year and 21 of 24 currency-unit fits were
# recorded as not converged (the other three of each were singular, where lme4 runs no
# derivative check), every one within 7.1e-5 log-likelihood units of the best spelling;
# carried exactly to the maximum, the gradient check still fired in 13 of 21 and 20 of
# 21 of them and the eigenvalue check in 21 of 21 and 20 of 21. Not in the finding: age
# in plain years is flagged by the eigenvalue check in Poisson fits (18 of 18 at 9,600
# rows, 3 of 6 at 1,200) and in binomial fits on a few thousand rows, a count model is
# flagged even with a standardised covariate once the intercept's curvature, about
# 2 x rows x mean count, passes 1e6 (a mean count of about 50 on 9,600 rows), and on
# maihda_health_data
# Obese ~ Gender + Race + Education + Age + (1 | Gender:Race:Education) was recorded as
# not converged 0.0002 standard errors from its maximum. lmer fits are not affected (the
# fixed effects are not among the optimised parameters), and lme4 runs no derivative
# check from 10,000 rows on.
#
# Per the user, when lme4's checks fire on a glmer fit whose optimiser reports success,
# the fit is measured with the fixed effects as coefficients of centred and scaled
# columns: within 0.1 standard errors of a stationary point (sqrt(g' H^-1 g / 2) of the
# deviance, the line ordinal fits are judged by) it counts as converged, lme4's
# messages are kept in `lme4_overruled`, and print() says what lme4 warned of and how
# far the fit is. Measured on 53 bobyqa stops with optimiser code 0, the distance
# agrees with the log-likelihood shortfall (sqrt(2 * shortfall); quartiles of the ratio
# 1.00 / 1.01 / 1.08), and all 29 stops 0.1 standard errors or more short stay flagged;
# an independent check on 427 such stops agreed with an exact-maximum oracle on every
# verdict. The deviance function is rebuilt from a copy of the fit: lme4 writes into
# vectors it shares with the fit it was rebuilt from (the covariance factor, an offset
# given as an argument).
#
# No block draws a random number.

lc14_observe <- function(expr) {
  warns <- character(0)
  res <- withCallingHandlers(
    suppressMessages(expr),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  list(res = res, warns = warns)
}

# gender x ses strata of `rows` rows: a stratum's latent outcome is its shift plus the
# same fixed spread plus 0.4 z, with the same ages in every stratum, and the age in
# exact affine spellings. yb is the latent outcome above 0, yc a count.
lc14_years <- c(53, 49, 55, 51, 55, 43, 46, 55, 71, 56, 19, 53, 23, 41, 63, 45, 53, 54,
                54, 55, 38, 61, 23, 36, 55, 47, 51, 53, 41, 66)
lc14_strata <- function(shift = c(-0.9, 0.2, 0.5, -0.4, 0.8, -0.2), rows = 30L) {
  g <- expand.grid(gender = c("F", "M"), ses = c("low", "mid", "high"),
                   stringsAsFactors = FALSE)
  d <- g[rep(seq_len(nrow(g)), each = rows), , drop = FALSE]
  spread <- rep(seq(-1.45, 1.45, length.out = rows), times = nrow(g))
  d$age <- rep(rep_len(lc14_years, rows), times = nrow(g))
  d$z <- (d$age - 45) / 12
  d$birth <- 2020 - d$age
  d$income <- 40000 + 1250 * (d$age - 45)
  latent <- rep(shift, each = rows) + spread + 0.4 * d$z
  d$yb <- as.integer(latent > 0)
  d$yc <- as.integer(floor(exp(1 + 0.5 * latent)))
  rownames(d) <- NULL
  d
}
lc14_fit <- function(x, family = "binomial", data = lc14_strata(), ...) {
  response <- if (family == "binomial") "yb" else "yc"
  formula <- stats::as.formula(paste(response, "~", x, "+ (1 | gender:ses)"))
  suppressMessages(suppressWarnings(
    fit_maihda(formula, data = data, engine = "lme4", family = family, ...)))
}
lc14_variance <- function(m) {
  vc <- as.data.frame(lme4::VarCorr(m$model))
  vc$vcov[vc$grp == "stratum"]
}
lc14_loglik <- function(m) as.numeric(stats::logLik(m$model))
# bobyqa in both stages with a final trust radius of 0.3: it stops short of the
# maximum and reports success.
lc14_coarse <- lme4::glmerControl(optimizer = c("bobyqa", "bobyqa"),
                                  optCtrl = list(rhobeg = 0.5, rhoend = 0.3))

test_that("the centring and scaling map gives the coefficients of centred and scaled columns", {
  X <- cbind("(Intercept)" = 1, male = c(0, 1, 1, 0, 1, 0),
             birth = c(1960, 1971, 1955, 1999, 1983, 1966))
  A <- maihda_centre_scale_map(X)
  scaled <- cbind(1, (X[, 2] - mean(X[, 2])) / sd(X[, 2]),
                  (X[, 3] - mean(X[, 3])) / sd(X[, 3]))
  expect_equal(unname(X %*% A), unname(scaled))
  # So the same linear predictor: beta = A gamma.
  gamma <- c(0.3, -0.8, 0.5)
  expect_equal(drop(X %*% (A %*% gamma)), drop(scaled %*% gamma))
  # The intercept is found where it stands, not assumed to be the first column.
  moved <- X[, c(2, 1, 3)]
  expect_equal(unname(moved %*% maihda_centre_scale_map(moved)), unname(scaled[, c(2, 1, 3)]))
  # Without an intercept the columns are only scaled, and a constant column is left
  # as it is.
  no_intercept <- cbind(two = 2, birth = X[, 3])
  B <- maihda_centre_scale_map(no_intercept)
  expect_equal(B, diag(c(1, 1 / sd(X[, 3]))))
})

test_that("central differences give a quadratic's gradient and Hessian and stop at a bound", {
  Q <- matrix(c(4, 1, 0.5, 1, 3, -0.2, 0.5, -0.2, 2), 3, 3)
  b <- c(-1, 0.5, 2)
  f <- function(x) 0.5 * sum(x * drop(Q %*% x)) + sum(b * x) + 7
  x <- c(0.3, -1.2, 2.5)
  d <- maihda_central_derivatives(f, x)
  expect_equal(d$gradient, drop(Q %*% x) + b, tolerance = 1e-6)
  expect_equal(d$hessian, Q, tolerance = 1e-6)
  # A variance parameter within a step of its bound of 0 cannot be differenced both ways.
  expect_null(maihda_central_derivatives(f, c(5e-5, 1, 1), lower = c(0, -Inf, -Inf)))
  expect_false(is.null(maihda_central_derivatives(f, c(2e-4, 1, 1), lower = c(0, -Inf, -Inf))))
})

test_that("the distance from a stationary point does not depend on the covariate's units", {
  fits <- lapply(c(z = "z", age = "age", birth = "birth", income = "income"), lc14_fit)
  distance <- vapply(fits, function(m) maihda_lme4_stationary_distance(m$model), numeric(1))
  expect_true(all(distance < 1e-3))
  # The four spellings reach the same maximum.
  expect_lt(diff(range(vapply(fits, lc14_loglik, numeric(1)))), 1e-4)
  expect_lt(diff(range(vapply(fits, lc14_variance, numeric(1)))), 1e-4)

  # A fit that stopped short: the distance is the log-likelihood shortfall it implies.
  short <- lc14_fit("z", control = lc14_coarse)
  expect_identical(as.integer(short$model@optinfo$conv$opt), 0L)
  shortfall <- lc14_loglik(fits$z) - lc14_loglik(short)
  expect_gt(shortfall, 0.005)
  # A ratio, not expect_equal(tolerance =): under this suite's testthat edition a
  # scalar passes that whenever the absolute difference is below the tolerance.
  expect_lt(abs(maihda_lme4_stationary_distance(short$model) / sqrt(2 * shortfall) - 1), 0.1)
})

test_that("the distance is NA where the fit cannot say", {
  m <- lc14_fit("birth")$model
  # An lmer fit, whose fixed effects are not among the optimised parameters.
  g <- suppressMessages(suppressWarnings(
    fit_maihda(yc ~ birth + (1 | gender:ses), data = lc14_strata())))
  expect_true(is.na(maihda_lme4_stationary_distance(g$model)))
  # A glmer fit at nAGQ = 0, whose fixed effects were not optimised with the variance:
  # its deviance function takes the variance alone and is not called.
  m0 <- lc14_fit("birth", nAGQ = 0)
  expect_identical(utils::capture.output(r0 <- maihda_lme4_stationary_distance(m0$model)),
                   character(0))
  expect_true(is.na(r0))
  # A deviance that does not reproduce the fit's log-likelihood.
  deviance <- lme4::getME(m, "devfun")
  own <- c(lme4::getME(m, "theta"), lme4::fixef(m))
  expect_true(is.finite(maihda_lme4_distance_with(m, deviance)))
  expect_true(is.na(maihda_lme4_distance_with(m, function(p) deviance(p) + 1)))
  # A Hessian that is not positive definite: curved up along the variance parameter,
  # down along the fixed effects, with a gradient along the first only.
  target <- -2 * as.numeric(stats::logLik(m))
  saddle <- function(p) {
    target + (p[1] - own[1])^2 + 0.1 * (p[1] - own[1]) - sum((p[-1] - own[-1])^2)
  }
  expect_true(is.na(maihda_lme4_distance_with(m, saddle)))
  invisible(deviance(own))
  # Not a glmer fit at all.
  expect_true(is.na(maihda_lme4_stationary_distance(list())))
})

test_that("a glmer fit at its maximum is converged whatever the covariate's units", {
  for (family in c("binomial", "poisson")) {
    for (x in c("birth", "income")) {
      m <- lc14_fit(x, family)
      said <- m$model@optinfo$conv$lme4$messages
      # The premise: lme4's checks fire on this fit, and its optimiser reports success.
      expect_gt(length(said), 0L)
      expect_identical(as.integer(m$model@optinfo$conv$opt), 0L)
      dg <- m$diagnostics
      expect_true(dg$converged)
      expect_identical(dg$messages, character(0))
      expect_identical(dg$lme4_overruled, as.character(said))
      expect_lt(dg$stationary_distance, 1e-3)
    }
    for (x in c("z", "age")) {
      m <- lc14_fit(x, family)
      expect_length(m$model@optinfo$conv$lme4$messages, 0L)
      expect_true(m$diagnostics$converged)
      expect_null(m$diagnostics$stationary_distance)
      expect_null(m$diagnostics$lme4_overruled)
    }
  }
})

test_that("a fit short of a stationary point stays not converged, and says how far", {
  short <- lc14_fit("z", control = lc14_coarse)
  dg <- short$diagnostics
  expect_false(dg$converged)
  expect_identical(dg$messages, as.character(short$model@optinfo$conv$lme4$messages))
  expect_gte(dg$stationary_distance, 0.1)
  expect_null(dg$lme4_overruled)
  lines <- maihda_format_fit_diagnostics(dg)
  expect_identical(lines[1], "Convergence warnings reported by lme4:")
  expect_identical(utils::tail(lines, 2L), c(
    sprintf("In centred and scaled units the fit's gradient and Hessian put it %s standard",
            format(round(dg$stationary_distance, 2), nsmall = 2)),
    "  errors from a stationary point."))
})

test_that("the line is a tenth of a standard error, from both sides", {
  at <- function(distance) {
    local_mocked_bindings(maihda_lme4_stationary_distance = function(model) distance)
    lc14_fit("birth")$diagnostics
  }
  below <- at(0.0999)
  expect_true(below$converged)
  expect_identical(maihda_format_fit_diagnostics(below),
                   c("lme4 warned that the model failed to converge; in centred and scaled units the",
                     "  fit is 0.099 standard errors from the maximum."))
  on_line <- at(0.1)
  expect_false(on_line$converged)
  expect_null(on_line$lme4_overruled)
  expect_identical(utils::tail(maihda_format_fit_diagnostics(on_line), 2L),
                   c("In centred and scaled units the fit's gradient and Hessian put it 0.10 standard",
                     "  errors from a stationary point."))
})

test_that("a fit that cannot be measured keeps lme4's verdict", {
  local_mocked_bindings(maihda_lme4_stationary_distance = function(model) NA_real_)
  m <- lc14_fit("birth")
  said <- as.character(m$model@optinfo$conv$lme4$messages)
  expect_gt(length(said), 0L)
  expect_false(m$diagnostics$converged)
  expect_identical(m$diagnostics$messages, said)
  expect_null(m$diagnostics$stationary_distance)
  expect_null(m$diagnostics$lme4_overruled)
  expect_identical(maihda_format_fit_diagnostics(m$diagnostics)[1],
                   "Convergence warnings reported by lme4:")
})

test_that("a fit whose optimiser reports failure is not measured again", {
  capped <- lc14_fit("birth", control = lme4::glmerControl(
    optimizer = c("bobyqa", "bobyqa"), optCtrl = list(maxfun = 10)))
  expect_identical(as.integer(capped$model@optinfo$conv$opt), 1L)
  dg <- capped$diagnostics
  expect_false(dg$converged)
  expect_true(any(grepl("maximum number of function evaluations", dg$messages)))
  expect_null(dg$stationary_distance)
  expect_null(dg$lme4_overruled)
})

test_that("measuring the fit leaves it as lme4 made it", {
  d <- lc14_strata()
  d$stratum <- interaction(d$gender, d$ses, drop = TRUE)
  m <- suppressWarnings(lme4::glmer(yb ~ birth + (1 | stratum), data = d, family = binomial))
  # Copies, not references: lme4 writes into some of the fit's vectors in place.
  state <- function(m) unserialize(serialize(list(
    ranef = lme4::ranef(m), predict = stats::predict(m), b = lme4::getME(m, "b"),
    u = lme4::getME(m, "u"), theta = lme4::getME(m, "theta"), fixef = lme4::fixef(m),
    lambda = lme4::getME(m, "Lambdat")@x + 0, loglik = stats::logLik(m),
    vcov = as.matrix(stats::vcov(m)), frame = lapply(m@frame, function(column) column),
    offset = m@resp$offset + 0, mu = m@resp$mu + 0), NULL))
  # An offset given as an argument is a column of the model frame, and lme4 writes the
  # linear predictor into it.
  exposed <- transform(d, exposure = log(1 + (seq_len(nrow(d)) %% 5)))
  mo <- suppressWarnings(lme4::glmer(yc ~ birth + (1 | stratum), data = exposed,
                                     family = poisson, offset = exposure))
  before <- state(mo)
  expect_true(is.finite(maihda_lme4_stationary_distance(mo)))
  expect_identical(state(mo), before)
  expect_identical(mo@frame[["(offset)"]], exposed$exposure)
  # Through fit_maihda(), where lme4's checks fire on such a fit and it is measured.
  exposed$stratum <- NULL
  fo <- suppressMessages(suppressWarnings(fit_maihda(
    yc ~ birth + (1 | gender:ses), data = exposed, family = "poisson", offset = exposure)))
  expect_true(is.finite(fo$diagnostics$stationary_distance))
  expect_identical(fo$model@frame[["(offset)"]], exposed$exposure)
  expect_identical(fo$model@resp$offset, exposed$exposure)

  before <- state(m)
  dg <- maihda_fit_diagnostics(m)
  expect_true(dg$converged)
  expect_true(is.finite(dg$stationary_distance))
  expect_identical(state(m), before)
  # And it measures the same again.
  expect_identical(maihda_lme4_stationary_distance(m), dg$stationary_distance)

  # With the intercept alone, the last difference is taken away from the fit's variance.
  m1 <- suppressWarnings(lme4::glmer(yb ~ 1 + (1 | stratum), data = d, family = binomial))
  before <- state(m1)
  expect_true(is.finite(maihda_lme4_stationary_distance(m1)))
  expect_identical(state(m1), before)
  # A measurement that stops part way leaves the fit as it was too.
  local_mocked_bindings(maihda_central_derivatives = function(fun, x, ...) {
    fun(x + c(0.05, rep(0, length(x) - 1L)))
    stop("stopped part way")
  })
  before <- state(m)
  expect_true(is.na(maihda_lme4_stationary_distance(m)))
  expect_identical(state(m), before)
})

test_that("print and the app say what lme4 warned of and how far the fit is", {
  record <- function(said, distance, converged = distance < 0.1, engine = "lme4") {
    r <- list(engine = engine, singular = FALSE, singular_terms = character(0),
              stratum_singular = NA, converged = converged,
              messages = if (converged) character(0) else said,
              stationary_distance = distance)
    if (converged) r$lme4_overruled <- said
    r
  }
  gradient <- "Model failed to converge with max|grad| = 0.0048 (tol = 0.002, component 1)"
  eigen_said <- c("Model is nearly unidentifiable: very large eigenvalue\n - Rescale variables?",
                  "Model is nearly unidentifiable: large eigenvalue ratio\n - Rescale variables?")
  expect_identical(
    maihda_format_fit_diagnostics(record(c(gradient, eigen_said), 1.44e-4)),
    c("lme4 warned that the model failed to converge; in centred and scaled units the",
      "  fit is 0.00014 standard errors from the maximum."))
  expect_identical(
    maihda_format_fit_diagnostics(record(eigen_said, 3e-5)),
    c("lme4 warned that the model is nearly unidentifiable; in centred and scaled",
      "  units the fit is less than 0.0001 standard errors from the maximum."))
  expect_identical(
    maihda_format_fit_diagnostics(record(
      "Hessian is numerically singular: parameters are not uniquely determined", 0.0288)),
    c("lme4's convergence checks raised a warning; in centred and scaled units the fit",
      "  is 0.028 standard errors from the maximum."))
  expect_identical(
    maihda_format_fit_diagnostics(record(gradient, 0.6241699)),
    c("Convergence warnings reported by lme4:", paste0("  - ", gradient),
      "In centred and scaled units the fit's gradient and Hessian put it 0.62 standard",
      "  errors from a stationary point."))
  # Distances as printed: under the line never as the line.
  expect_identical(
    vapply(c(0.0999, 0.0288, 0.0024, 1.44e-4, 3e-5, 0.1, 0.6241699, 1.473, 12.345),
           maihda_format_stationary_distance, character(1)),
    c("0.099", "0.028", "0.0024", "0.00014", "less than 0.0001", "0.10", "0.62", "1.47",
      "12.35"))
  # Nothing for a fit that was not measured, nor for another engine's record.
  plain <- record(gradient, 0.01)
  plain$stationary_distance <- NULL
  plain$lme4_overruled <- NULL
  expect_identical(maihda_format_fit_diagnostics(plain), character(0))
  expect_identical(maihda_format_fit_diagnostics(record(gradient, 0.01, engine = "ordinal")),
                   character(0))

  # The app: each sentence is one item, the convergence messages sub-items.
  items <- maihda_app_diagnostics_items(maihda_format_fit_diagnostics(
    record(gradient, 1.44e-4)))
  expect_length(items, 1L)
  expect_identical(items[[1]]$text, paste0(
    "lme4 warned that the model failed to converge; in centred and scaled units the ",
    "fit is 0.00014 standard errors from the maximum."))
  items <- maihda_app_diagnostics_items(maihda_format_fit_diagnostics(
    record(gradient, 0.6241699)))
  expect_length(items, 2L)
  expect_identical(items[[1]]$sub, gradient)
  expect_identical(items[[2]]$text, paste0(
    "In centred and scaled units the fit's gradient and Hessian put it 0.62 standard ",
    "errors from a stationary point."))

  # And a live fit prints it.
  out <- utils::capture.output(print(lc14_fit("birth")))
  expect_true(any(grepl("lme4 warned that the model failed to converge; in centred and scaled",
                        out, fixed = TRUE)))
  expect_false(any(grepl("Convergence warnings reported by lme4", out, fixed = TRUE)))
})

test_that("compare_maihda_groups() does not name groups whose fits are at their maximum", {
  d <- rbind(cbind(country = "A", lc14_strata()),
             cbind(country = "B", lc14_strata(shift = c(0.6, -0.5, 0.1, 0.9, -0.7, 0.3))))
  obs <- lc14_observe(compare_maihda_groups(yb ~ birth + (1 | gender:ses), data = d,
                                            group = "country", family = "binomial"))
  # lme4 warns about every group's fit, as before ...
  expect_true(any(grepl("failed to converge", obs$warns, fixed = TRUE)))
  # ... but no group is named as not converged.
  expect_false(any(grepl("did not converge", obs$warns, fixed = TRUE)))
  expect_identical(obs$res$status, c("ok", "ok"))
})

test_that("the package's own data: an adjusted logistic MAIHDA with age in years", {
  skip_on_cran()
  data("maihda_health_data", package = "MAIHDA", envir = environment())
  hc <- maihda_health_data[stats::complete.cases(
    maihda_health_data[, c("Obese", "Age", "Gender", "Race", "Education")]), ]
  m <- suppressMessages(suppressWarnings(fit_maihda(
    Obese ~ Gender + Race + Education + Age + (1 | Gender:Race:Education),
    data = hc, family = "binomial")))
  said <- m$model@optinfo$conv$lme4$messages
  # The premise: lme4 calls it nearly unidentifiable for the units of Age alone.
  expect_true(any(grepl("very large eigenvalue", said, fixed = TRUE)))
  expect_false(any(grepl("failed to converge", said, fixed = TRUE)))
  expect_true(m$diagnostics$converged)
  expect_lt(m$diagnostics$stationary_distance, 0.01)
  out <- utils::capture.output(print(m))
  expect_true(any(grepl("lme4 warned that the model is nearly unidentifiable", out,
                        fixed = TRUE)))
})
