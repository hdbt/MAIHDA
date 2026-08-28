# Audit pass 2026-08-25: longitudinal fit-quality reporting (external report --
# singular fits and "puzzling" slope PCVs on irregular, unevenly-spaced data).
#
# The estimation was never at fault: fit_maihda(id=, time=) reproduces a hand-coded
# lmer of the same formula to the last digit. Four REPORTING defects were confirmed:
#
#   1 lme4 files "boundary (singular) fit" in optinfo$conv$lme4$messages, and
#     maihda_fit_diagnostics() read that list as the convergence verdict -- so a fit
#     whose optimizer returned code 0 was reported converged = FALSE, and the one
#     condition printed twice (once as "Singular fit", once under a "Convergence
#     warnings" header). It also double-counted such groups in
#     compare_maihda_groups()'s singular AND non-converged warnings.
#   2 The singular-fit report never said WHICH random-effects block was degenerate,
#     and asserted the between-stratum variance "may be unreliable" either way --
#     false when the boundary sits in a non-stratum block, which is routine for a
#     longitudinal fit whose (time | id) block has no person-level slope variation.
#   3 A RANK-DEFICIENT stratum growth block (a perfect intercept-slope correlation)
#     went undetected: maihda_stratum_growth_at_boundary_lme4() requires EVERY
#     variance to be at the boundary and maihda_variance_at_boundary() tests the
#     magnitude of a(t)'Sigma a(t), so both pass while the slope variance is pinned
#     to the intercept variance and PCV_slope is a ratio of boundary artefacts. The
#     longitudinal print() also had no adjusted-at-boundary note at all, though
#     calculate_pcv() has carried one for the cross-sectional path all along.
#   4 print() on a longitudinal maihda() showed the NULL model's diagnostics only, so
#     a non-converged ADJUSTED growth fit delivered the headline "N% additive"
#     silently.
#
# Plus the missing escape hatch the report needed: stratum_slope = FALSE fits
# (time | id) + (1 | stratum), giving a time-constant between-stratum variance.

# ---- Finding 1: a singular fit is not a convergence failure ------------------

test_that("maihda_fit_diagnostics separates singularity from non-convergence", {
  skip_on_cran()
  set.seed(101)
  dd <- data.frame(g = rep(seq_len(20), each = 5), x = rnorm(100))
  dd$y <- dd$x + rnorm(100)          # no group-level variance to find
  fm <- suppressMessages(lme4::lmer(y ~ x + (x | g), data = dd))
  skip_if_not(isTRUE(lme4::isSingular(fm)), "fixture did not land on the boundary")

  # lme4 files the singular note among its own messages, but the OPTIMIZER converged.
  expect_true(any(grepl("boundary (singular)", fm@optinfo$conv$lme4$messages,
                        fixed = TRUE)))
  expect_equal(as.integer(fm@optinfo$conv$opt), 0L)

  diag <- maihda_fit_diagnostics(fm)
  expect_true(diag$singular)
  expect_true(diag$converged)                       # was FALSE before this pass
  expect_length(diag$messages, 0)                   # not re-reported as convergence

  # The report says "Singular fit" once and never under a convergence header.
  lines <- maihda_format_fit_diagnostics(diag)
  expect_true(any(grepl("Singular fit", lines)))
  expect_false(any(grepl("Convergence warnings", lines)))

  # A genuine optimizer failure is still reported, singular or not.
  bad <- fm
  bad@optinfo$conv$opt <- 1L
  bad@optinfo$message <- "bobyqa -- maximum number of function evaluations exceeded"
  diag_bad <- maihda_fit_diagnostics(bad)
  expect_false(diag_bad$converged)
  expect_true(any(grepl("function evaluations", diag_bad$messages)))
})

# ---- Finding 2: name the degenerate block ------------------------------------

test_that("maihda_singular_terms_lme4 identifies the boundary block per term", {
  skip_on_cran()
  set.seed(202)
  n <- 40
  dd <- data.frame(g = rep(seq_len(n), each = 5))
  dd$x <- rnorm(nrow(dd))
  dd$h <- rep(seq_len(8), length.out = nrow(dd))
  b0 <- rnorm(n, 0, 1.5); b1 <- rnorm(n, 0, 0.8)
  dd$y <- b0[dd$g] + b1[dd$g] * dd$x + rnorm(nrow(dd), 0, 0.4)
  # A healthy two-block fit: (x | g) contributes 3 theta entries (Cholesky diagonal
  # at 1 and 3, off-diagonal at 2), (1 | h) contributes the 4th.
  fm <- suppressMessages(lme4::lmer(y ~ x + (x | g) + (1 | h), data = dd))
  expect_false(lme4::isSingular(fm))
  expect_length(maihda_singular_terms_lme4(fm), 0)
  expect_equal(names(lme4::getME(fm, "cnms")), c("g", "h"))

  # Drive a chosen Cholesky diagonal to the boundary directly, rather than hoping an
  # optimiser lands there: this pins the slice arithmetic that maps theta entries to
  # grouping factors, which is the part that can silently mis-attribute a block.
  zero_theta <- function(model, i) { model@theta[i] <- 0; model }
  expect_equal(maihda_singular_terms_lme4(zero_theta(fm, 4L)), "h")
  expect_equal(maihda_singular_terms_lme4(zero_theta(fm, 3L)), "g")
  expect_equal(maihda_singular_terms_lme4(zero_theta(fm, 1L)), "g")
  expect_setequal(maihda_singular_terms_lme4(zero_theta(zero_theta(fm, 3L), 4L)),
                  c("g", "h"))
  # The off-diagonal entry is unconstrained (lower = -Inf): a zero there is an
  # uncorrelated block, not a degenerate one.
  expect_length(maihda_singular_terms_lme4(zero_theta(fm, 2L)), 0)

  # Unreadable model -> character(0), so callers fall back to the global flag.
  expect_length(maihda_singular_terms_lme4(list()), 0)
})

test_that("the singular-fit report scopes its VPC caveat to the stratum block", {
  # Pure formatting: no fit needed, so the assertion does not ride on an optimiser.
  non_stratum <- list(engine = "lme4", singular = TRUE,
                      singular_terms = "id", stratum_singular = FALSE,
                      converged = TRUE, messages = character(0))
  lines <- maihda_format_fit_diagnostics(non_stratum)
  expect_true(any(grepl("block for 'id'", lines, fixed = TRUE)))
  expect_true(any(grepl("NOT at the boundary", lines, fixed = TRUE)))
  expect_false(any(grepl("may be unreliable", lines)))

  strat <- list(engine = "lme4", singular = TRUE,
                singular_terms = c("id", "stratum"), stratum_singular = TRUE,
                converged = TRUE, messages = character(0))
  lines2 <- maihda_format_fit_diagnostics(strat)
  expect_true(any(grepl("blocks for 'id' and 'stratum'", lines2, fixed = TRUE)))
  expect_true(any(grepl("may be unreliable", lines2)))

  # No per-term information (non-lme4 engine, unreadable theta): keep the old,
  # conservative wording rather than silently dropping the caveat.
  unknown <- list(engine = "wemix", singular = TRUE,
                  singular_terms = character(0), stratum_singular = NA,
                  converged = TRUE, messages = character(0))
  lines3 <- maihda_format_fit_diagnostics(unknown)
  expect_true(any(grepl("at least one variance component", lines3)))
  expect_true(any(grepl("may be unreliable", lines3)))
})

# ---- Finding 3: rank-deficient stratum growth block --------------------------

test_that("a rank-deficient stratum block is detected and flagged on the PCV", {
  skip_on_cran()
  set.seed(303)
  # Irregular, unevenly spaced occasions -- the reported design. The stratum-level
  # slope signal is negligible, so lme4 pins the intercept-slope correlation at 1.
  np <- 240
  pat <- data.frame(id = sprintf("p%03d", seq_len(np)),
                    a = sample(c("1", "2"), np, TRUE),
                    b = sample(c("1", "2"), np, TRUE),
                    cc = sample(c("1", "2", "3"), np, TRUE),
                    stringsAsFactors = FALSE)
  k <- sample(1:5, np, TRUE)
  rows <- do.call(rbind, lapply(seq_len(np), function(i)
    data.frame(id = pat$id[i], time = sort(round(runif(k[i], 0, 12), 2)),
               stringsAsFactors = FALSE)))
  dz <- merge(rows, pat, by = "id")
  st <- interaction(dz$a, dz$b, dz$cc, drop = TRUE)
  u0 <- stats::setNames(rnorm(nlevels(st), 0, 0.45), levels(st))
  v0 <- stats::setNames(rnorm(np, 0, 0.8), pat$id)
  dz$y <- 2 + u0[as.character(st)] + v0[dz$id] + 0.1 * dz$time +
    rnorm(nrow(dz), 0, 0.7)

  m <- suppressMessages(suppressWarnings(
    fit_maihda(y ~ time + (1 | a:b:cc), data = dz, id = "id", time = "time")))
  blk <- maihda_re_block(m, "stratum")
  rank_def <- maihda_stratum_growth_rank_deficient(m)
  skip_if_not(rank_def, "fixture stratum block was not rank-deficient")

  # The point of the finding: every variance is comfortably positive and
  # a(t)'Sigma a(t) is healthy, so BOTH pre-existing boundary tests say "fine".
  expect_true(all(diag(blk) > 0))
  expect_false(maihda_stratum_growth_at_boundary_lme4(m$model))
  expect_false(maihda_variance_at_boundary(
    maihda_var_at_time(blk, 0), maihda_residual_variance_lme4(m$model)))
  # ... while the correlation is pinned at the boundary.
  expect_gt(abs(stats::cov2cor(blk)[1, 2]), 0.999)

  a <- suppressMessages(suppressWarnings(
    maihda(y ~ time + (1 | a:b:cc), data = dz, id = "id", time = "time",
           decomposition = "longitudinal")))
  expect_true(isTRUE(a$pcv$null_rank_deficient))
  out <- paste(utils::capture.output(print(a$pcv)), collapse = "\n")
  expect_match(out, "RANK-DEFICIENT")
  expect_match(out, "stratum_slope = FALSE")

  # PCV_slope is NA, not a number: under a pinned correlation the slope variance is
  # not a free parameter, and the null and adjusted fits collapse onto different
  # directions, so the ratio does not compare the same quantity before and after
  # adjustment. Left ungated it is explosive rather than merely uncertain -- values
  # near -200 across replicates of this design.
  expect_true(is.na(a$pcv$pcv_slope))
  # The raw variances are still reported: seeing a denominator of ~1e-6 is exactly
  # what tells the analyst why the ratio is unavailable.
  expect_true(is.finite(a$pcv$var_slope_null))
  expect_true(is.finite(a$pcv$var_slope_adjusted))
  # The intercept side survives -- it rides on the intercept variance, which stays
  # estimable -- so gating it too would discard a usable result.
  expect_true(is.finite(a$pcv$pcv_intercept))
  expect_true(all(is.finite(a$pcv$pcv_t$pcv)))
})

test_that("PCV_slope is gated on identifiability, not on magnitude", {
  # The pre-existing denominator guard scales sqrt(var) against the RESIDUAL sd, but
  # a slope variance is in outcome-units^2 per time-unit^2, so its 1e-4 relative
  # threshold has no fixed meaning for that cell and cannot catch this. Pin the gap
  # with the observed numbers: a denominator small enough to send the ratio past
  # -190 still tests an order of magnitude above the threshold.
  expect_false(maihda_variance_at_boundary(1.317e-06, 0.49))
  expect_gt(abs((1.317e-06 - 2.602e-04) / 1.317e-06), 190)
  # Rescaling the slope variance into outcome units over the time span does not
  # rescue a magnitude test either -- identifiability is the right question.
  expect_false(maihda_variance_at_boundary(1.317e-06 * 12^2, 0.49))
})

test_that("the slope gate is the slope's own component, not the whole block", {
  skip_on_cran()
  # A QUADRATIC growth block can lose its I(time^2) dimension while estimating the
  # linear slope variance perfectly well. Gating PCV_slope on a block-level
  # "is anything at the boundary" verdict would discard a sound proportion, so the
  # gate reads the slope's own Cholesky diagonal.
  data(maihda_long_data, envir = environment())
  a <- suppressMessages(suppressWarnings(
    maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
           data = maihda_long_data, id = "id", time = "wave",
           time_degree = 2, decomposition = "longitudinal")))
  comps <- maihda_re_boundary_components_lme4(a$model$model)[["stratum"]]
  expect_length(comps, 3L)
  expect_false(unname(comps[1]))     # intercept fine
  expect_false(unname(comps[2]))     # LINEAR slope fine  -> PCV_slope survives
  expect_true(unname(comps[3]))      # quadratic term at the boundary

  expect_true(isTRUE(a$pcv$null_rank_deficient))        # block IS degenerate ...
  expect_false(isTRUE(a$pcv$null_slope_at_boundary))    # ... but not in the slope
  expect_true(is.finite(a$pcv$pcv_slope))               # so the proportion stands
  expect_gt(a$pcv$pcv_slope, 0)
  expect_lt(a$pcv$pcv_slope, 1)

  # The note still fires (the block is degenerate) but must not claim PCV_slope is NA.
  out <- paste(utils::capture.output(print(a$pcv)), collapse = "
")
  expect_match(out, "RANK-DEFICIENT")
  expect_match(out, "slope component itself is still estimable")
  expect_false(grepl("reported as NA", out, fixed = TRUE))
})

test_that("print.maihda_long_pcv carries an adjusted-at-boundary note", {
  # Pure print: build the object directly so the note does not depend on an
  # optimiser landing on the boundary.
  S <- matrix(c(0.2, 0.01, 0.01, 0.002), 2, 2,
              dimnames = list(c("(Intercept)", "t"), c("(Intercept)", "t")))
  x <- structure(
    list(pcv_intercept = 0.99, pcv_slope = 0.9,
         var_baseline_null = 0.2, var_baseline_adjusted = 0.002,
         var_slope_null = 0.002, var_slope_adjusted = 2e-4,
         ref_time = 0, time_center = 0,
         pcv_t = data.frame(time = 0:2, var_null = 0.2, var_adjusted = 0.002,
                            pcv = 0.99),
         Sigma_stratum_null = S, Sigma_stratum_adjusted = S,
         time = "t", time_degree = 1L, ml_refit = FALSE, estimation = "fitted",
         estimation_used = "fitted", null_at_boundary = FALSE,
         adjusted_at_boundary = TRUE, null_rank_deficient = FALSE,
         adjusted_rank_deficient = FALSE, stratum_slope = TRUE),
    class = "maihda_long_pcv")
  out <- paste(utils::capture.output(print(x)), collapse = "\n")
  expect_match(out, "ADJUSTED model's between-stratum")
  expect_match(out, "pinned near 100")

  # ... and stays silent when neither fit is at the boundary.
  x$adjusted_at_boundary <- FALSE
  out2 <- paste(utils::capture.output(print(x)), collapse = "\n")
  expect_false(grepl("ADJUSTED model's between-stratum", out2))
  expect_false(grepl("RANK-DEFICIENT", out2))
})

# ---- Finding 4: the adjusted growth fit's diagnostics are printed ------------

test_that("a longitudinal maihda() prints BOTH fits' diagnostics", {
  skip_on_cran()
  data(maihda_long_data, envir = environment())
  a <- suppressMessages(suppressWarnings(
    maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
           data = maihda_long_data, id = "id", time = "wave",
           decomposition = "longitudinal")))
  # Inject a diagnostic on the ADJUSTED fit only: before this pass print() read
  # x$model$diagnostics alone, so nothing about the adjusted fit could ever appear.
  a$model$diagnostics$singular <- FALSE
  a$model$diagnostics$converged <- TRUE
  a$model$diagnostics$messages <- character(0)
  a$model_adjusted$diagnostics$converged <- FALSE
  a$model_adjusted$diagnostics$messages <- "Model failed to converge with max|grad|"
  out <- paste(utils::capture.output(print(a)), collapse = "\n")
  expect_match(out, "Fit diagnostics (adjusted model)", fixed = TRUE)
  expect_match(out, "max|grad|", fixed = TRUE)

  # The label is opt-in, so single-fit callers keep the unlabelled header.
  plain <- paste(utils::capture.output(maihda_print_fit_diagnostics(
    list(engine = "lme4", singular = FALSE, singular_terms = character(0),
         stratum_singular = NA, converged = FALSE, messages = "boom"))),
    collapse = "\n")
  expect_match(plain, "Fit diagnostics:", fixed = TRUE)
})

# ---- stratum_slope = FALSE: a time-constant between-stratum variance ---------

test_that("stratum_slope = FALSE fits (time | id) + (1 | stratum)", {
  skip_on_cran()
  data(maihda_long_data, envir = environment())
  m <- suppressMessages(suppressWarnings(
    fit_maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
               data = maihda_long_data, id = "id", time = "wave",
               stratum_slope = FALSE)))
  bars <- vapply(reformulas::findbars(stats::formula(m$model)), deparse1, character(1))
  expect_true(any(grepl("^wave \\| id$", bars)))       # person growth block kept
  expect_true(any(grepl("^1 \\| stratum$", bars)))      # stratum intercept only
  expect_false(isTRUE(m$longitudinal_info$stratum_slope))

  # The stratum block is 1x1, so the between-stratum variance is constant in time.
  blk <- maihda_re_block(m, "stratum")
  expect_equal(dim(blk), c(1L, 1L))
  s <- suppressWarnings(summary(m))
  expect_equal(length(unique(round(s$longitudinal$var_stratum_t, 12))), 1L)
  # The VPC is NOT constant: the person-level slope variance is still in the
  # denominator. This is the distinction the argument's documentation turns on.
  expect_gt(length(unique(round(s$longitudinal$vpc_t$estimate, 8))), 1L)

  # Downstream scalar-per-stratum reporting stays well defined (no NA baseline
  # from padding intercept-only coefficients out to a slope column).
  ps <- predict(m, type = "strata")
  expect_false("slope" %in% names(ps))
  expect_false(anyNA(ps$baseline))
  expect_equal(ps$baseline, ps$intercept)
})

test_that("a stratum_slope = FALSE decomposition reports PCV_intercept only", {
  skip_on_cran()
  data(maihda_long_data, envir = environment())
  a <- suppressMessages(suppressWarnings(
    maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
           data = maihda_long_data, id = "id", time = "wave",
           decomposition = "longitudinal", stratum_slope = FALSE)))
  # BOTH growth models inherit the intercept-only stratum block.
  for (fit in list(a$model, a$model_adjusted)) {
    bars <- vapply(reformulas::findbars(stats::formula(fit$model)), deparse1,
                   character(1))
    expect_true(any(grepl("^1 \\| stratum$", bars)))
  }
  expect_true(is.na(a$pcv$pcv_slope))
  expect_true(is.finite(a$pcv$pcv_intercept))
  expect_false(isTRUE(a$pcv$stratum_slope))
  # PCV(t) is flat, because both blocks are time-constant.
  expect_equal(length(unique(round(a$pcv$pcv_t$pcv, 12))), 1L)

  out <- paste(utils::capture.output(print(a$pcv)), collapse = "\n")
  expect_match(out, "random intercept only")
  # No reported slope decomposition: neither the slope-variance line nor a
  # PCV_slope figure. (The explanatory note names "PCV_slope" while saying there
  # is none, so match the reported VALUE, not the bare word.)
  expect_false(grepl("variance at baseline", out, fixed = TRUE))
  expect_false(grepl("PCV_slope: +[0-9]", out))
  # The closing explanation must not advertise a PCV_slope that is not on the page.
  expect_false(grepl("a high PCV_slope", out, fixed = TRUE))
})

test_that("stratum_slope is validated and its refusal messages stay truthful", {
  skip_on_cran()
  data(maihda_long_data, envir = environment())
  expect_error(
    fit_maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
               data = maihda_long_data, id = "id", time = "wave",
               stratum_slope = NA),
    "single TRUE or FALSE")

  m <- suppressMessages(suppressWarnings(
    fit_maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
               data = maihda_long_data, id = "id", time = "wave",
               stratum_slope = FALSE)))
  # Still refused (the cross-sectional ranking path is not wired for a growth fit),
  # but no longer on the false ground that each stratum is a trajectory.
  expect_error(extract_between_variance(m), "time-constant")
  err <- tryCatch(plot(m, type = "predicted"), error = conditionMessage)
  expect_false(grepl("not a single value", err, fixed = TRUE))
  expect_match(err, "predict\\(type = \"strata\"\\)")

  # The slope-carrying default keeps the trajectory wording.
  m2 <- suppressMessages(suppressWarnings(
    fit_maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
               data = maihda_long_data, id = "id", time = "wave")))
  expect_error(extract_between_variance(m2), "time-varying")
})

# ---- the estimation itself was never in question ----------------------------

test_that("a longitudinal fit_maihda matches the hand-coded lmer exactly", {
  skip_on_cran()
  data(maihda_long_data, envir = environment())
  m <- suppressMessages(suppressWarnings(
    fit_maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
               data = maihda_long_data, id = "id", time = "wave")))
  hand_data <- m$data           # carries the resolved 'stratum' column
  hand <- suppressMessages(suppressWarnings(lme4::lmer(
    wellbeing ~ wave + (wave | id) + (wave | stratum), data = hand_data)))
  expect_equal(as.numeric(stats::logLik(m$model)),
               as.numeric(stats::logLik(hand)), tolerance = 1e-8)
  expect_equal(unname(lme4::getME(m$model, "theta")),
               unname(lme4::getME(hand, "theta")), tolerance = 1e-6)
})

# ---- Finding 4b: the cross-sectional two-model branch --------------------------
# Same defect as finding 4, in the branch that printed NEITHER model's diagnostics.
# The caveat existed all along on the PCV object, but print.maihda_analysis()
# formats its PCV inline instead of calling print(x$pcv), so it never reached the
# default output. The rule agreed with the maintainer: suppress ONLY the adjusted
# model's stratum-block singularity (the expected outcome of an additive
# decomposition, and what the boundary note is for); print everything else.

test_that("maihda_format_fit_diagnostics can suppress the expected stratum singularity", {
  strat <- list(engine = "lme4", singular = TRUE, singular_terms = "stratum",
                stratum_singular = TRUE, converged = TRUE, messages = character(0),
                adequacy = NULL)
  expect_true(any(grepl("Singular fit", maihda_format_fit_diagnostics(strat))))
  expect_length(maihda_format_fit_diagnostics(strat, suppress_stratum_singular = TRUE), 0)

  # A NON-stratum boundary is never expected and survives the suppression.
  ctx <- strat
  ctx$singular_terms <- "country"; ctx$stratum_singular <- FALSE
  lines <- maihda_format_fit_diagnostics(ctx, suppress_stratum_singular = TRUE)
  expect_true(any(grepl("block for 'country'", lines, fixed = TRUE)))

  # Convergence failures are never suppressed either.
  nc <- strat
  nc$converged <- FALSE; nc$messages <- "Model failed to converge with max|grad|"
  lines2 <- maihda_format_fit_diagnostics(nc, suppress_stratum_singular = TRUE)
  expect_false(any(grepl("Singular fit", lines2)))
  expect_true(any(grepl("max|grad|", lines2, fixed = TRUE)))

  # include_adequacy = FALSE drops only the adequacy caveats.
  ad <- list(engine = "lme4", singular = FALSE, singular_terms = character(0),
             stratum_singular = NA, converged = TRUE, messages = character(0),
             adequacy = list(overdispersion = list(flag = TRUE, ratio = 4.2,
                                                   chisq = 420, rdf = 100, p = 1e-8)))
  expect_gt(length(maihda_format_fit_diagnostics(ad)), 0)
  expect_length(maihda_format_fit_diagnostics(ad, include_adequacy = FALSE), 0)
})

test_that("a cross-sectional maihda() print surfaces the adjusted fit's boundary", {
  skip_on_cran()
  set.seed(9)
  n <- 1800
  d <- data.frame(a = sample(c("1","2"), n, TRUE), b = sample(c("1","2"), n, TRUE),
                  cc = sample(c("1","2","3"), n, TRUE), stringsAsFactors = FALSE)
  # Purely additive strata: the adjusted fit is expected to hit the boundary.
  d$y <- 2 + 0.8*(d$a=="2") + 0.5*(d$b=="2") + 0.6*(d$cc=="3") + rnorm(n, 0, 1)
  a <- suppressMessages(suppressWarnings(
    maihda(y ~ (1 | a:b:cc), data = d, interactions = FALSE)))
  skip_if_not(isTRUE(a$pcv$adjusted_at_boundary), "fixture adjusted fit not at boundary")

  out <- paste(utils::capture.output(print(a)), collapse = "\n")
  expect_match(out, "the adjusted model's between-stratum variance is at the singularity",
               fixed = TRUE)
  expect_match(out, "pinned near 100", fixed = TRUE)
  # ... and does NOT cry wolf with a diagnostics banner for the expected case.
  expect_false(grepl("Fit diagnostics (adjusted model)", out, fixed = TRUE))

  # The note is shared with print.pcv_result(), which keeps its own subject wording.
  pcv_out <- paste(utils::capture.output(print(a$pcv)), collapse = "\n")
  expect_match(pcv_out, "Model 2's between-stratum variance", fixed = TRUE)
  expect_match(pcv_out, "pinned near 100", fixed = TRUE)
})

test_that("a cross-sectional maihda() print surfaces non-converged and non-stratum fits", {
  skip_on_cran()
  set.seed(9)
  n <- 1800
  d <- data.frame(a = sample(c("1","2"), n, TRUE), b = sample(c("1","2"), n, TRUE),
                  cc = sample(c("1","2","3"), n, TRUE), stringsAsFactors = FALSE)
  st <- interaction(d$a, d$b, d$cc, drop = TRUE)
  iv <- stats::setNames(rnorm(nlevels(st), 0, 0.35), levels(st))
  d$y <- 2 + 0.8*(d$a=="2") + 0.5*(d$b=="2") + 0.6*(d$cc=="3") +
    iv[as.character(st)] + rnorm(n, 0, 1)
  a <- suppressMessages(suppressWarnings(
    maihda(y ~ (1 | a:b:cc), data = d, interactions = FALSE)))
  expect_false(grepl("Fit diagnostics", paste(utils::capture.output(print(a)),
                                              collapse = "\n")))   # healthy: silent

  # A non-converged adjusted fit must not stay silent (the original defect).
  nc <- a
  nc$model_adjusted$diagnostics$converged <- FALSE
  nc$model_adjusted$diagnostics$messages <- "Model failed to converge with max|grad| = 0.0071"
  out_nc <- paste(utils::capture.output(print(nc)), collapse = "\n")
  expect_match(out_nc, "Fit diagnostics (adjusted model)", fixed = TRUE)
  expect_match(out_nc, "max|grad| = 0.0071", fixed = TRUE)

  # A singular adjusted fit in a NON-stratum block is not the expected case, so the
  # suppression must not swallow it.
  ctx <- a
  ctx$model_adjusted$diagnostics$singular <- TRUE
  ctx$model_adjusted$diagnostics$singular_terms <- "country"
  ctx$model_adjusted$diagnostics$stratum_singular <- FALSE
  out_ctx <- paste(utils::capture.output(print(ctx)), collapse = "\n")
  expect_match(out_ctx, "Fit diagnostics (adjusted model)", fixed = TRUE)
  expect_match(out_ctx, "block for 'country'", fixed = TRUE)
})

# ---- Bell et al. (2024) eq. (5) trajectory VPCs ------------------------------
#
# Added in the same pass, from the same correspondence: the package reported only
# the discriminatory-accuracy VPC(t), which keeps the occasion-level residual in the
# denominator, while the paper it cites defines intercept and slope VPCs that
# exclude it. Both are now reported, with the distinction documented.

test_that("the trajectory VPCs are Bell eq. (5), excluding level-1 variance", {
  skip_on_cran()
  data(maihda_long_data, envir = environment())
  m <- suppressMessages(suppressWarnings(
    fit_maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
               data = maihda_long_data, id = "id", time = "wave")))
  s <- suppressWarnings(summary(m))
  Ss <- s$longitudinal$Sigma_stratum
  Si <- s$longitudinal$Sigma_id

  # Exactly the paper's ratios, computed from the raw covariance-block cells.
  expect_equal(s$longitudinal$vpc_intercept, Ss[1, 1] / (Ss[1, 1] + Si[1, 1]))
  expect_equal(s$longitudinal$vpc_slope, Ss[2, 2] / (Ss[2, 2] + Si[2, 2]))

  # The point of the distinction: the headline VPC keeps sigma^2_e, so it is
  # strictly smaller than the intercept VPC on the same fit.
  vpc_headline <- s$vpc$estimate
  expect_gt(s$longitudinal$vpc_intercept, vpc_headline)
  expect_equal(vpc_headline,
               Ss[1, 1] / (Ss[1, 1] + Si[1, 1] + s$longitudinal$var_resid))

  out <- paste(utils::capture.output(print(s)), collapse = "\n")
  expect_match(out, "Trajectory VPCs (Bell et al. 2024, eq. 5", fixed = TRUE)
  expect_match(out, "occasion-level variance excluded", fixed = TRUE)
  # The one sentence that says which question each VPC answers, so the larger
  # number is not reported as "the VPC".
  expect_match(out, "how intersectionally patterned trajectories are", fixed = TRUE)
  expect_match(out, "?summary.maihda_model", fixed = TRUE)
})

test_that("vpc_slope is NA when there is no stratum slope to take a share of", {
  skip_on_cran()
  data(maihda_long_data, envir = environment())
  m <- suppressMessages(suppressWarnings(
    fit_maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
               data = maihda_long_data, id = "id", time = "wave",
               stratum_slope = FALSE)))
  s <- suppressWarnings(summary(m))
  # A 1x1 stratum block would make the arithmetic return 0/(0 + slope_var_id) = 0,
  # which reads as "no intersectional patterning of trajectories" rather than
  # "not estimated". It must be NA.
  expect_true(is.na(s$longitudinal$vpc_slope))
  expect_true(is.finite(s$longitudinal$vpc_intercept))
  out <- paste(utils::capture.output(print(s)), collapse = "\n")
  expect_match(out, "Slope: NA", fixed = TRUE)
})

test_that("maihda_longitudinal_trajectory_vpc guards degenerate denominators", {
  S2 <- matrix(c(0.3, 0.02, 0.02, 0.004), 2, 2)
  I2 <- matrix(c(0.9, 0.01, 0.01, 0.006), 2, 2)
  v <- maihda_longitudinal_trajectory_vpc(S2, I2, 0)
  expect_equal(v$vpc_intercept, 0.3 / (0.3 + 0.9))
  expect_equal(v$vpc_slope, 0.004 / (0.004 + 0.006))

  # Intercept-only stratum block -> no slope share.
  S1 <- matrix(0.3, 1, 1)
  expect_true(is.na(maihda_longitudinal_trajectory_vpc(S1, I2, 0)$vpc_slope))

  # An all-zero pair is 0/0, not 0.
  Z <- matrix(0, 2, 2)
  z <- maihda_longitudinal_trajectory_vpc(Z, Z, 0)
  expect_true(is.na(z$vpc_intercept))
  expect_true(is.na(z$vpc_slope))
})
