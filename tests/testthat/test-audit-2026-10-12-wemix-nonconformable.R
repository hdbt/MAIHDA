# Audit 2026-10-12 (packet 2026-10-12-wemix-nonconformable): fit_maihda(engine = "wemix",
# family = "binomial") stopped with
#   c("non-conformable matrix dimensions in Matrix::solve(R_l2, l2_points[rep(...
# when every stratum held exactly as many successes as failures -- CONFIRMED, and
# WeMix's own: WeMix::mix() called directly on the same rows raises it.
#
# WeMix 4.0.3 evaluates a binomial likelihood by adaptive quadrature centred on the
# conditional modes of the stratum effects. makeAdaptedPts() wraps the modes with
# Matrix() and adds that single column to a dense strata x points matrix. Matrix()
# returns a sparse matrix when more than half of its entries are exact zeros, and
# Matrix 1.7-3 adds a dense and a sparse matrix only when their dimensions agree. So
# mix() stops at an evaluation where more than half of the modes are exactly 0.0; in
# every case measured that was the first one, at the fixed effects of its unweighted
# starting fit. A stratum's mode is an exact zero there when those fixed effects put a
# probability of exactly one half on its rows and its residuals cancel without rounding.
#
# Measured on 92 designs, each fitted as WeMix ships and with that one line changed to
# keep the column dense: 20 stopped this way. With the same number of successes s of 40
# in each of 8 strata, s = 1..39, only s = 20 did; the other 38 stopped at a zero
# variance (the error of packet 2026-10-07-wemix-boundary). "Every stratum at one half"
# is neither necessary nor sufficient. 200 of 400 with the rows sorted does not stop
# this way: rounding leaves the modes at 1e-19. Counts 20, 20, 20, 20, 20, 20, 12, 28
# of 40 with the rows sorted do, although a fit exists: the same rows in the order used
# below return a variance of 0.0638, as lme4 with 13 quadrature points does. Of the 20,
# 18 end in the zero-variance error once the line is changed and 2 return a fit. The
# linear path and nQuad = 1 never reach that code.
#
# Per the user the fit stays an error, one that says what happened. The condition
# carries no call and Matrix translates its text ("nicht konforme Matrixdimensionen in
# ..."), so maihda_fit_wemix() tells the case by where mix() was when it stopped: a
# makeAdaptedPts() frame, the `+` it had called, a sparse `muHat2` in that frame and
# more than half of its `muHat` at zero. Every other WeMix error and every returned fit
# is as WeMix produced it.
#
# No block draws a random number.

nc12_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

# The error a call stops with, or NULL when it returns -- so a call that stops (or does
# not) fails each assertion it feeds rather than aborting the block.
nc12_error <- function(expr) tryCatch({ expr; NULL }, error = function(e) e)
nc12_text <- function(e) if (inherits(e, "error")) conditionMessage(e) else ""

# What maihda_wemix_zero_modes() reads off the stack at the error `expr` stops with,
# that error, and every warning raised on the way.
nc12_at_error <- function(expr) {
  got <- "no error"
  warns <- character(0)
  err <- tryCatch(
    withCallingHandlers(
      { expr; NULL },
      error = function(e) got <<- maihda_wemix_zero_modes(sys.calls(), sys.frames()),
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        invokeRestart("muffleWarning")
      }),
    error = function(e) e)
  list(got = got, error = err, warns = warns)
}

# The stack an error handler sees when mix() stops in makeAdaptedPts(): the calls, and
# one environment per call. `column` stands in, by its class, for the matrix Matrix()
# made of the modes.
nc12_pts_call <- quote(makeAdaptedPts(muHat, R/sqrt(PsiVec), nQuad, levels, q2, q3,
                                      n_l2, n_top))
nc12_add_call <- quote(solve(R_l2, l2_points[rep(1:nrow(l2_points), times = n_l2), ]) +
                         muHat2)
nc12_sparse <- structure(list(), class = c("dgCMatrix", "sparseMatrix"))
nc12_dense <- structure(list(), class = c("dgeMatrix", "denseMatrix"))

nc12_stack <- function(modes, column = nc12_sparse, pts = nc12_pts_call,
                       inner = nc12_add_call) {
  calls <- list(quote(mix(formula, data = data, weights = weights)),
                quote(getGradHess(func = main_lnl, x = est)),
                quote(func(x)),
                pts,
                inner,
                quote(.Ops.checkDim(d1, d2)),
                quote(stop(msg, call. = FALSE)),
                quote(h(simpleError(msg, call))))
  frames <- lapply(calls, function(cl) new.env())
  assign("muHat", modes, envir = frames[[4L]])
  assign("muHat2", column, envir = frames[[4L]])
  list(calls = calls, frames = frames)
}
nc12_read <- function(stack) maihda_wemix_zero_modes(stack$calls, stack$frames)

# Six modes at zero and two that are not: the reading wanted of it.
nc12_modes <- c(0, 0, 0, 0, 0, 0, -1.17, 1.17)
nc12_six <- c(zero = 6L, strata = 8L)

# Stand-ins for mix() and makeAdaptedPts() on a live stack. The addition is refused
# inside `+` and without a call, as Matrix refuses it; with `early` the same words are
# raised before the addition instead.
nc12_live <- function(modes, column, early = FALSE) {
  `+` <- function(e1, e2) stop("any words", call. = FALSE)
  nc12_solve <- function() {
    if (early) stop("any words", call. = FALSE)
    "the points"
  }
  makeAdaptedPts <- function(muHat) {
    muHat2 <- column
    points <- nc12_solve()
    points + muHat2
  }
  mix <- function() makeAdaptedPts(modes)
  mix()
}

# Eight strata of 40 rows with the given numbers of successes and unit weights. Within
# a stratum the rows alternate 1, 0 for as long as both remain, so that the residuals
# of a stratum at one half cancel without rounding in any order of summation.
nc12_data <- function(successes, n = 40L) {
  g <- expand.grid(d1 = 0:1, d2 = 0:1, d3 = 0:1)
  d <- g[rep(1:8, each = n), ]
  d$stratum <- rep(1:8, each = n)
  d$y <- unlist(lapply(successes, function(s) {
    pairs <- min(s, n - s)
    c(rep(c(1L, 0L), pairs), rep(as.integer(s > n - s), n - 2L * pairs))
  }))
  for (v in c("d1", "d2", "d3")) d[[v]] <- factor(d[[v]])
  d$w <- 1
  rownames(d) <- NULL
  d
}
nc12_all_half <- function() nc12_data(rep(20L, 8L))
# Six strata at one half and two that are not: a positive between-stratum variance.
nc12_six_half <- function() nc12_data(c(20L, 20L, 20L, 20L, 20L, 20L, 12L, 28L))

nc12_null <- y ~ 1 + (1 | stratum)
# Starting values with the intercept at exactly zero. mix() takes its starting fixed
# effects from them, so its first likelihood evaluation puts a probability of exactly
# one half on every row, whatever the unweighted fit of the platform returns (on the
# six-at-half rows above that fit's intercept is 1.3e-16 here, which is not enough).
nc12_start <- c(0, 0.1)

nc12_fit <- function(data, ...) {
  fit_maihda(nc12_null, data = data, engine = "wemix", family = "binomial",
             sampling_weights = "w", ...)
}
nc12_mix <- function(data, ...) {
  data$l2 <- 1
  WeMix::mix(nc12_null, data = data, weights = c("w", "l2"),
             family = stats::binomial(link = "logit"), ...)
}

nc12_diagnosis <- "WeMix::mix() stopped while evaluating the likelihood of this binomial model"
nc12_zero_variance <- "estimated the between-stratum variance of this binomial model at zero"

test_that("the stop in makeAdaptedPts() is told by the stack and not by the error's text", {
  # As it arrives: the column sparse, six of the eight modes (or all of them) at zero.
  expect_identical(nc12_read(nc12_stack(nc12_modes)), nc12_six)
  expect_identical(nc12_read(nc12_stack(rep(0, 8))), c(zero = 8L, strata = 8L))
  # A namespace-qualified call names the same function.
  qualified <- nc12_pts_call
  qualified[[1L]] <- call(":::", as.name("WeMix"), as.name("makeAdaptedPts"))
  expect_identical(nc12_read(nc12_stack(nc12_modes, pts = qualified)), nc12_six)
  # A mode that cannot be read does not count as a zero.
  expect_identical(nc12_read(nc12_stack(replace(nc12_modes, 8L, NA))), nc12_six)

  # Not it: a dense column; another call made from makeAdaptedPts(); another function.
  expect_null(nc12_read(nc12_stack(nc12_modes, column = nc12_dense)))
  expect_null(nc12_read(nc12_stack(nc12_modes, inner = quote(solve(R_l2, B)))))
  expect_null(nc12_read(nc12_stack(nc12_modes, pts = quote(makeAdaptedPts2(muHat)))))
  expect_null(nc12_read(nc12_stack(nc12_modes, pts = quote(obj$makeAdaptedPts(muHat)))))
  # Nor a sparse column beside modes of which no more than half are exact zeros.
  expect_null(nc12_read(nc12_stack(c(0, 0, 0, 0, 1, 1, 1, 1))))
  expect_null(nc12_read(nc12_stack(c(1e-19, -1e-19, 4e-21, 1e-300, 5e-324, 0, 1, 1))))
  # Nor a frame without the modes or without the column.
  expect_null(nc12_read(nc12_stack("six of them")))
  bare <- nc12_stack(nc12_modes)
  rm("muHat", envir = bare$frames[[4L]])
  expect_null(nc12_read(bare))
  bare <- nc12_stack(nc12_modes)
  rm("muHat2", envir = bare$frames[[4L]])
  expect_null(nc12_read(bare))
  # Nor makeAdaptedPts() as the last call: nothing was being added.
  whole <- nc12_stack(nc12_modes)
  expect_null(maihda_wemix_zero_modes(whole$calls[1:4], whole$frames[1:4]))

  # The zero-variance error of the same engine has no makeAdaptedPts() frame at all.
  near_pd <- list(quote(mix(formula, data = data, weights = weights)),
                  quote(getGradHess(main_lnl, est)), quote(func(xm)),
                  quote(nearPD2(Matrix(constructSigma(vcov, nRE, n_l2, n_top)))),
                  quote(nearPD(X, posd.tol = tol * sqrt(.Machine$double.eps))),
                  quote(stop(msg)))
  expect_null(maihda_wemix_zero_modes(near_pd, lapply(near_pd, function(cl) new.env())))

  # With two makeAdaptedPts() frames the innermost is the one read.
  outer_calls <- list(quote(makeAdaptedPts(muHat)), quote(func(x)))
  outer_frames <- list(new.env(), new.env())
  assign("muHat", c(1, 2, 3, 4), envir = outer_frames[[1L]])
  assign("muHat2", nc12_dense, envir = outer_frames[[1L]])
  expect_identical(
    maihda_wemix_zero_modes(c(outer_calls, whole$calls), c(outer_frames, whole$frames)),
    nc12_six)
  expect_null(
    maihda_wemix_zero_modes(c(whole$calls[1:5], outer_calls),
                            c(whole$frames[1:5], outer_frames)))
})

test_that("reading the stack can neither fail nor warn", {
  nothing <- function(calls, frames) {
    warned <- 0L
    got <- withCallingHandlers(
      maihda_wemix_zero_modes(calls, frames),
      warning = function(w) {
        warned <<- warned + 1L
        invokeRestart("muffleWarning")
      })
    is.null(got) && warned == 0L
  }
  expect_true(nothing(NULL, NULL))
  expect_true(nothing(list(), list()))
  expect_true(nothing(1:3, "a"))
  whole <- nc12_stack(nc12_modes)
  # Fewer frames than calls; something that is not an environment where the frame was.
  expect_true(nothing(whole$calls, whole$frames[1:3]))
  expect_true(nothing(whole$calls, replace(whole$frames, 4L, list("not a frame"))))
  # A binding that fails when read, and one that warns.
  failing <- nc12_stack(nc12_modes)
  rm("muHat2", envir = failing$frames[[4L]])
  makeActiveBinding("muHat2", function() stop("unreadable"), failing$frames[[4L]])
  expect_true(nothing(failing$calls, failing$frames))
  warning_one <- nc12_stack(nc12_modes)
  rm("muHat", envir = warning_one$frames[[4L]])
  makeActiveBinding("muHat", function() {
    warning("read with a warning")
    nc12_modes
  }, warning_one$frames[[4L]])
  expect_true(nothing(warning_one$calls, warning_one$frames))
})

test_that("an error handler reads the case off a live stack", {
  refused <- nc12_at_error(nc12_live(nc12_modes, nc12_sparse))
  expect_s3_class(refused$error, "error")
  expect_null(conditionCall(refused$error))
  expect_identical(refused$got, nc12_six)
  expect_identical(refused$warns, character(0))

  # The same refusal with a dense column, and the same words raised before the addition.
  expect_null(nc12_at_error(nc12_live(nc12_modes, nc12_dense))$got)
  before <- nc12_at_error(nc12_live(nc12_modes, nc12_sparse, early = TRUE))
  expect_s3_class(before$error, "error")
  expect_null(before$got)
})

test_that("the message gives the count, the cause, the way round and WeMix's words", {
  msg <- maihda_wemix_zero_modes_message(simpleError("WeMix's own words"), nc12_six)
  lines <- strsplit(msg, "\n", fixed = TRUE)[[1L]]
  expect_length(lines, 4L)
  expect_identical(
    lines[1L],
    paste0("engine = \"wemix\": ", nc12_diagnosis, ": the estimated effect of 6 of the ",
           "8 strata is exactly zero there, and mix() cannot go on when more than half ",
           "of them are."))
  expect_match(lines[2L], "exactly as many (weighted) successes as failures", fixed = TRUE)
  expect_match(lines[2L], "a limitation of WeMix", fixed = TRUE)
  expect_match(lines[3L], "nQuad = 1 (the Laplace approximation)", fixed = TRUE)
  expect_identical(lines[4L], "  WeMix's own error: WeMix's own words")
  # The first line is what a compare_maihda_groups() status cell holds.
  expect_identical(maihda_first_line(msg), lines[1L])
  # Counts are written out in full.
  expect_match(
    maihda_wemix_zero_modes_message(simpleError("x"), c(zero = 100000L, strata = 200000L)),
    "of 100000 of the 200000 strata", fixed = TRUE)
})

test_that("mix() stops in makeAdaptedPts() and fit_maihda() says what happened", {
  skip_on_cran()
  skip_if_not_installed("WeMix")
  doubled <- nc12_all_half()
  doubled$w <- 2
  cases <- list(list(data = nc12_all_half(), zero = 8L),
                list(data = nc12_six_half(), zero = 6L),
                list(data = doubled, zero = 8L))
  for (case in cases) {
    info <- paste(case$zero, "strata at one half, weight", case$data$w[1L])
    # WeMix on these rows with no package code in between, read by a handler of the
    # test's own. If this stops holding, WeMix or Matrix has changed and the branch in
    # maihda_fit_wemix() is no longer reached.
    own <- nc12_at_error(nc12_quiet(nc12_mix(case$data, start = nc12_start)))
    expect_s3_class(own$error, "error")
    expect_identical(own$got, c(zero = case$zero, strata = 8L), info = info)

    err <- nc12_error(nc12_fit(case$data, start = nc12_start))
    expect_s3_class(err, "error")
    expect_match(nc12_text(err), nc12_diagnosis, fixed = TRUE, info = info)
    expect_match(nc12_text(err),
                 paste0("of ", case$zero, " of the 8 strata is exactly zero"),
                 fixed = TRUE, info = info)
    expect_match(nc12_text(err), "nQuad = 1", fixed = TRUE, info = info)
    # WeMix's words are passed on, in whatever language Matrix put them.
    expect_match(nc12_text(err), paste0("WeMix's own error: ", nc12_text(own$error)),
                 fixed = TRUE, info = info)
  }

  # The reported call, with mix()'s own starting values. All strata are alike, so
  # mix() has no fit to return; which of its two errors comes first depends on whether
  # the unweighted fit's intercept is an exact zero. Either way the package's words.
  reported <- nc12_error(nc12_fit(nc12_all_half()))
  expect_s3_class(reported, "error")
  expect_true(grepl(nc12_diagnosis, nc12_text(reported), fixed = TRUE) ||
                grepl(nc12_zero_variance, nc12_text(reported), fixed = TRUE))
})

test_that("maihda() and the group comparison pass the diagnosis on", {
  skip_on_cran()
  skip_if_not_installed("WeMix")
  d <- nc12_all_half()
  d$stratum <- NULL
  err <- nc12_error(nc12_quiet(maihda(
    y ~ 1 + (1 | d1:d2:d3), data = d, engine = "wemix", family = "binomial",
    sampling_weights = "w", start = nc12_start)))
  expect_match(nc12_text(err), nc12_diagnosis, fixed = TRUE)

  both <- rbind(cbind(nc12_all_half(), grp = "g1"), cbind(nc12_six_half(), grp = "g2"))
  warns <- character(0)
  tab <- withCallingHandlers(
    suppressMessages(compare_maihda_groups(
      nc12_null, data = both, group = "grp", engine = "wemix", family = "binomial",
      sampling_weights = "w", start = nc12_start)),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  # The status cell holds the first line: the diagnosis with the group's own count.
  expect_identical(
    tab$status,
    paste0("fit failed: engine = \"wemix\": ", nc12_diagnosis, ": the estimated effect ",
           "of ", c(8, 6), " of the 8 strata is exactly zero there, and mix() cannot go ",
           "on when more than half of them are."))
  expect_length(grep("the model fit failed for group(s) g1, g2", warns, fixed = TRUE), 1L)
})

test_that("one quadrature point, or a start off the exact point, returns the fit", {
  skip_on_cran()
  skip_if_not_installed("WeMix")
  d <- nc12_six_half()
  lme4_variance <- function(points) {
    ref <- nc12_quiet(lme4::glmer(nc12_null, data = d, family = stats::binomial(),
                                  nAGQ = points))
    as.numeric(lme4::VarCorr(ref)$stratum)
  }

  # The way round the message names: the Laplace approximation, as lme4 computes it.
  laplace <- tryCatch(nc12_quiet(nc12_fit(d, start = nc12_start, nQuad = 1)),
                      error = function(e) e)
  expect_s3_class(laplace, "maihda_model")
  expect_equal(as.numeric(laplace$model$vars[1L]), lme4_variance(1L), tolerance = 1e-4)

  # The fit the error stood in front of: mix()'s 13 points from a start a millionth
  # away, beside lme4's 13 points.
  moved <- tryCatch(nc12_quiet(nc12_fit(d, start = c(1e-6, 0.1))), error = function(e) e)
  expect_s3_class(moved, "maihda_model")
  expect_equal(as.numeric(moved$model$vars[1L]), lme4_variance(13L), tolerance = 1e-3)
  expect_gt(as.numeric(moved$model$vars[1L]), 0.05)

  # With every stratum at one half there is no such fit: one quadrature point leads on
  # to the zero-variance diagnosis.
  alike <- nc12_error(nc12_fit(nc12_all_half(), start = nc12_start, nQuad = 1))
  expect_match(nc12_text(alike), nc12_zero_variance, fixed = TRUE)
})

test_that("other WeMix errors and returned fits are left as WeMix made them", {
  skip_on_cran()
  skip_if_not_installed("WeMix")

  # A zero variance without exact zeros among the modes (12 of 40 in every stratum)
  # keeps its own diagnosis.
  zero_variance <- nc12_error(nc12_fit(nc12_data(rep(12L, 8L))))
  expect_match(nc12_text(zero_variance), nc12_zero_variance, fixed = TRUE)
  expect_false(grepl(nc12_diagnosis, nc12_text(zero_variance), fixed = TRUE))

  # An error of WeMix's that is neither leaves in WeMix's words, with its call.
  other <- nc12_error(nc12_fit(nc12_all_half(), start = nc12_start, nQuad = 0))
  expect_match(nc12_text(other), "nQuad", fixed = TRUE)
  expect_false(grepl(nc12_diagnosis, nc12_text(other), fixed = TRUE))
  expect_true(is.call(conditionCall(other)))

  # A fit mix() returns is the fit mix() made.
  d <- nc12_six_half()
  fit <- nc12_quiet(nc12_fit(d, start = c(1e-6, 0.1)))
  own <- nc12_quiet(nc12_mix(d, start = c(1e-6, 0.1)))
  for (field in c("coef", "vars", "SE", "lnl", "theta", "invHessian", "cov_mat")) {
    expect_equal(fit$model[[field]], own[[field]], tolerance = 1e-8, info = field)
  }
})
