# Audit pass 2026-10-13: the lme4 longitudinal PCV with few, weakly differing strata.
#
# Two things were measured on known-truth growth designs (12 strata, REML through
# maihda()), and neither changes a returned value:
#
#   1 Both PCVs are understated on average when the strata differ little after
#     adjustment next to the precision of one stratum's own estimate. That is said in
#     the help and the vignette; nothing here pins prose.
#   2 A rank-deficient LINEAR stratum block is rank 1, so its between-stratum variance
#     over time is the perfect square (s1 + s2 t)^2 with one zero. Near the NULL
#     model's zero the PCV divides by what the collapse leaves over (a PCV_intercept
#     of -2298% was reproduced, printed beside "carried by the intercept variance,
#     which stays estimable"); at the ADJUSTED model's zero it is exactly 1. print()
#     now names both times instead of asserting that the PCV is safe.
#
# No random numbers are drawn, so the suite's random stream is left as it was found.

# A rank 1 linear block with stratum deviations s[1] + s[2] * t.
rank1_block <- function(s, time = "wave") {
  matrix(c(s[1]^2, s[1] * s[2], s[1] * s[2], s[2]^2), 2, 2,
         dimnames = list(c("(Intercept)", time), c("(Intercept)", time)))
}

# A maihda_long_pcv as maihda_longitudinal_pcv() builds it, from two stratum blocks.
long_pcv_stub <- function(Sn, Sa, times, center = 0, time = "wave",
                          null_rank_deficient = TRUE, adjusted_rank_deficient = TRUE,
                          null_at_boundary = FALSE, adjusted_at_boundary = FALSE,
                          null_slope_at_boundary = null_rank_deficient) {
  vn <- maihda_var_at_time(Sn, times - center)
  va <- maihda_var_at_time(Sa, times - center)
  structure(
    list(pcv_intercept = (vn[1] - va[1]) / vn[1],
         pcv_slope = if (null_slope_at_boundary) NA_real_ else 1 - Sa[2, 2] / Sn[2, 2],
         var_baseline_null = vn[1], var_baseline_adjusted = va[1],
         var_slope_null = Sn[2, 2], var_slope_adjusted = Sa[2, 2],
         ref_time = times[1], time_center = center,
         pcv_t = data.frame(time = times, var_null = vn, var_adjusted = va,
                            pcv = (vn - va) / vn),
         Sigma_stratum_null = Sn, Sigma_stratum_adjusted = Sa,
         time = time, time_degree = nrow(Sn) - 1L, ml_refit = FALSE,
         estimation = "fitted", estimation_used = "fitted",
         null_at_boundary = null_at_boundary,
         adjusted_at_boundary = adjusted_at_boundary,
         null_rank_deficient = null_rank_deficient,
         null_slope_at_boundary = null_slope_at_boundary,
         adjusted_rank_deficient = adjusted_rank_deficient, stratum_slope = TRUE,
         variance_summary = "estimate", Sigma_summary = "estimate"),
    class = "maihda_long_pcv")
}

# The printed lines without colour codes, and the same text as one line.
printed_lines <- function(x) cli::ansi_strip(utils::capture.output(print(x)))
printed <- function(x) gsub("\\s+", " ", paste(printed_lines(x), collapse = " "))

old_claim <- "are still reported: they are carried by the intercept variance"

# ---- the time at which a rank-deficient linear block has no variance ---------

test_that("the zero time is where the block's own variance function vanishes", {
  s <- c(0.3, -0.1)                       # deviations 0.3 - 0.1 t: zero at t = 3
  S <- rank1_block(s)
  expect_equal(maihda_long_pcv_zero_time(S), 3)
  expect_equal(maihda_var_at_time(S, 3), 0)
  # The package's own a(t)' Sigma a(t) is the perfect square, on either side of it.
  tt <- c(-2, 0, 1.5, 3, 4, 10)
  expect_equal(maihda_var_at_time(S, tt), (s[1] + s[2] * tt)^2)

  # A positive correlation puts the zero before the origin.
  expect_equal(maihda_long_pcv_zero_time(rank1_block(c(0.3, 0.1))), -3)
  # An intercept variance of exactly zero puts it AT the origin.
  expect_equal(maihda_long_pcv_zero_time(rank1_block(c(0, 0.2))), 0)

  # The block is in the fit's centered coordinates; the time is on the user's scale.
  expect_equal(maihda_long_pcv_zero_time(S, center = 2000), 2003)
  expect_equal(maihda_var_at_time(S, 2003 - 2000), 0)
  expect_equal(maihda_long_pcv_zero_time(S, center = -5), -2)

  # Dimnames do not matter, and the result is a bare number.
  expect_equal(maihda_long_pcv_zero_time(unname(S), center = 1), 4)
  expect_null(attributes(maihda_long_pcv_zero_time(S)))
  expect_length(maihda_long_pcv_zero_time(S), 1L)
})

test_that("no time is named for a block that has none to name", {
  S <- rank1_block(c(0.3, -0.1))
  # Not a linear block: an intercept-only level, a quadratic block.
  expect_identical(maihda_long_pcv_zero_time(matrix(0.09, 1, 1)), NA_real_)
  S3 <- diag(c(0.09, 0.01, 0.001)); S3[1, 2] <- S3[2, 1] <- -0.03
  expect_identical(maihda_long_pcv_zero_time(S3), NA_real_)
  # No slope variance: the between-stratum variance is the same at every time.
  expect_identical(maihda_long_pcv_zero_time(rank1_block(c(0.3, 0))), NA_real_)
  S_neg <- S; S_neg[2, 2] <- -1e-12
  expect_identical(maihda_long_pcv_zero_time(S_neg), NA_real_)
  # Unreadable cells, or a centering that is not known.
  S_na <- S; S_na[1, 2] <- S_na[2, 1] <- NA_real_
  expect_identical(maihda_long_pcv_zero_time(S_na), NA_real_)
  S_na2 <- S; S_na2[2, 2] <- NA_real_
  expect_identical(maihda_long_pcv_zero_time(S_na2), NA_real_)
  expect_identical(maihda_long_pcv_zero_time(S, center = NULL), NA_real_)
  expect_identical(maihda_long_pcv_zero_time(S, center = NA_real_), NA_real_)
  expect_identical(maihda_long_pcv_zero_time(S, center = Inf), NA_real_)
  expect_identical(maihda_long_pcv_zero_time(NULL), NA_real_)
  expect_identical(maihda_long_pcv_zero_time(list()), NA_real_)
})

# ---- what the note says -------------------------------------------------------

test_that("the note names each rank-deficient block's time and what the PCV does there", {
  tt <- 0:4
  both <- maihda_long_pcv_zero_note("wave", 0.2619, 1.6373, tt)
  expect_identical(both, paste0(
    "PCV_intercept and PCV(wave) are still reported. A rank-deficient linear block ",
    "has no between-stratum variance at one time -- wave = 0.26 for the null model, ",
    "wave = 1.64 for the adjusted model -- and the PCV is not interpretable near ",
    "either: it divides by almost nothing near the null model's, and is 100% at the ",
    "adjusted model's. The further both lie from the reported times (0 to 4), the ",
    "more the PCV is carried by the intercept variance, which stays estimable."))

  null_only <- maihda_long_pcv_zero_note("wave", 0.2619, NA_real_, tt)
  expect_identical(null_only, paste0(
    "PCV_intercept and PCV(wave) are still reported. A rank-deficient linear block ",
    "has no between-stratum variance at one time, wave = 0.26 for the null model, ",
    "and the PCV is not interpretable near it: it divides by almost nothing there. ",
    "The further that time lies from the reported times (0 to 4), the more the PCV ",
    "is carried by the intercept variance, which stays estimable."))

  adj_only <- maihda_long_pcv_zero_note("wave", NA_real_, 1.6373, tt)
  expect_identical(adj_only, paste0(
    "PCV_intercept and PCV(wave) are still reported. A rank-deficient linear block ",
    "has no between-stratum variance at one time, wave = 1.64 for the adjusted model, ",
    "and the PCV is not interpretable near it: it is 100% there. The further that ",
    "time lies from the reported times (0 to 4), the more the PCV is carried by the ",
    "intercept variance, which stays estimable."))

  # Neither block has a time to name: the sentence for every other rank-deficient block.
  none <- paste0("PCV_intercept and PCV(wave) are still reported: they are carried by ",
                 "the intercept variance, which stays estimable.")
  expect_identical(maihda_long_pcv_zero_note("wave", NA_real_, NA_real_, tt), none)
  expect_identical(maihda_long_pcv_zero_note("wave", Inf, NaN, tt), none)
  expect_identical(maihda_long_pcv_zero_note("wave", numeric(0), NULL, tt), none)
  # The time variable is named as the user named it.
  expect_match(maihda_long_pcv_zero_note("age", 41.3, NA_real_, 40:44),
               "PCV(age) are still reported", fixed = TRUE)
})

test_that("the named time is rounded to the span of the reporting times", {
  note <- function(z, times) maihda_long_pcv_zero_note("t", z, NA_real_, times)
  z <- 12.34567
  expect_match(note(z, 0:4), "t = 12.35 for", fixed = TRUE)            # span 4
  expect_match(note(z, c(0, 9.9)), "t = 12.35 for", fixed = TRUE)      # span below 10
  expect_match(note(z, c(0, 10)), "t = 12.3 for", fixed = TRUE)        # span 10
  expect_match(note(z, c(40, 80)), "t = 12.3 for", fixed = TRUE)       # span 40
  expect_match(note(z, c(0, 400)), "t = 12 for", fixed = TRUE)         # span 400
  expect_match(note(z, c(0, 4000)), "t = 12 for", fixed = TRUE)        # never negative digits
  expect_match(note(z, c(0, 0.5)), "t = 12.346 for", fixed = TRUE)     # span 0.5
  expect_match(note(-7.5, 0:4), "t = -7.50 for", fixed = TRUE)
  # The reporting times are quoted as their range, whatever their order.
  expect_match(note(z, c(3, 0, 4, 1, 2)), "the reported times (0 to 4),", fixed = TRUE)
  expect_match(note(z, c(2000, 2004)), "the reported times (2000 to 2004),", fixed = TRUE)
  expect_match(note(z, c(0.5, NA, 2.5)), "the reported times (0.5 to 2.5),", fixed = TRUE)
  # One reporting time, or none: two decimals, and no range where there is none.
  expect_match(note(z, 3), "t = 12.35 for", fixed = TRUE)
  expect_match(note(z, numeric(0)), "t = 12.35 for", fixed = TRUE)
  expect_match(note(z, numeric(0)), "from the reported times, the more", fixed = TRUE)
  expect_match(note(z, NULL), "from the reported times, the more", fixed = TRUE)
})

# ---- print(): the stored blocks, the centering, the flags ----------------------

test_that("print() names the times of the flagged blocks on the user's time scale", {
  # Null deviations 0.03 - 0.10 t (zero at 0.3), adjusted 0.15 - 0.10 t (zero at 1.5),
  # on a calendar axis the fit centered at 2000.
  Sn <- rank1_block(c(0.03, -0.10), "year")
  Sa <- rank1_block(c(0.15, -0.10), "year")
  x <- long_pcv_stub(Sn, Sa, times = 2000:2004, center = 2000, time = "year")
  # The object is the case the note is about: a baseline PCV far below zero ...
  expect_equal(x$pcv_intercept, 1 - 0.0225 / 0.0009)
  expect_lt(x$pcv_intercept, -20)
  out <- printed(x)
  expect_match(out, "RANK-DEFICIENT", fixed = TRUE)
  expect_match(out, "year = 2000.30 for the null model, year = 2001.50 for the adjusted model",
               fixed = TRUE)
  expect_match(out, "it divides by almost nothing near the null model's, and is 100% at the adjusted model's",
               fixed = TRUE)
  expect_match(out, "from the reported times (2000 to 2004),", fixed = TRUE)
  expect_false(grepl(old_claim, out, fixed = TRUE))
  # ... and the rest of the note is still there.
  expect_match(out, "reported as NA", fixed = TRUE)
  expect_match(out, "stratum_slope = FALSE", fixed = TRUE)
  expect_match(out, "says so explicitly.", fixed = TRUE)
  # What the two sentences say is true of the object's own variances.
  at <- function(S, t) maihda_var_at_time(S, t - x$time_center)
  expect_equal(at(Sn, 2000.3), 0)
  expect_equal(1 - at(Sa, 2001.5) / at(Sn, 2001.5), 1)
  expect_lt(1 - at(Sa, 2000.31) / at(Sn, 2000.31), -1e4)
  # The paragraph is wrapped like the lines around it, and print() returns its input.
  utils::capture.output(res <- withVisible(print(x)))
  expect_identical(res$value, x)
  expect_false(res$visible)
  lines <- printed_lines(x)
  from <- grep("are still reported", lines, fixed = TRUE)
  expect_length(from, 1L)
  expect_true(all(nchar(lines[from:length(lines)]) <= 80L))
})

test_that("print() names a time only for a block that is flagged, linear and not all zero", {
  Sn <- rank1_block(c(0.03, -0.10))
  Sa <- rank1_block(c(0.15, -0.10))
  full <- matrix(c(0.2, 0.01, 0.01, 0.002), 2, 2,
                 dimnames = list(c("(Intercept)", "wave"), c("(Intercept)", "wave")))

  # Only the null block is flagged. The adjusted block would have a zero too, but a
  # block that is not flagged gets no time, whatever its cells look like.
  null_only <- printed(long_pcv_stub(Sn, Sa, 0:4, adjusted_rank_deficient = FALSE))
  expect_match(null_only, "The null model's between-stratum block is RANK-DEFICIENT",
               fixed = TRUE)
  expect_match(null_only, "at one time, wave = 0.30 for the null model, and the PCV",
               fixed = TRUE)
  expect_match(null_only, "it divides by almost nothing there", fixed = TRUE)
  expect_false(grepl("for the adjusted model", null_only, fixed = TRUE))
  expect_false(grepl("1.50", null_only, fixed = TRUE))

  # Only the adjusted block is flagged: the null slope is estimable and PCV_slope stays.
  adj_only <- printed(long_pcv_stub(full, Sa, 0:4, null_rank_deficient = FALSE))
  expect_match(adj_only, "The adjusted model's between-stratum block is RANK-DEFICIENT",
               fixed = TRUE)
  expect_match(adj_only, "at one time, wave = 1.50 for the adjusted model, and the PCV",
               fixed = TRUE)
  expect_match(adj_only, "it is 100% there", fixed = TRUE)
  expect_match(adj_only, "slope component itself is still estimable", fixed = TRUE)
  expect_false(grepl("for the null model", adj_only, fixed = TRUE))
  expect_false(grepl(old_claim, adj_only, fixed = TRUE))

  # A block with EVERY variance at the boundary has no variance at any time: it is
  # given no time, and its own note speaks for it.
  dust <- rank1_block(c(3e-7, -1e-7))
  all_zero <- printed(long_pcv_stub(Sn, dust, 0:4, adjusted_at_boundary = TRUE))
  expect_match(all_zero, "wave = 0.30 for the null model, and the PCV", fixed = TRUE)
  expect_false(grepl("for the adjusted model", all_zero, fixed = TRUE))
  expect_match(all_zero, "ADJUSTED model's between-stratum", fixed = TRUE)
  null_zero <- printed(long_pcv_stub(dust, Sa, 0:4, null_at_boundary = TRUE))
  expect_false(grepl("for the null model", null_zero, fixed = TRUE))
  expect_match(null_zero, "wave = 1.50 for the adjusted model, and the PCV", fixed = TRUE)

  # Neither block is flagged: no rank-deficiency note, and no time.
  healthy <- printed(long_pcv_stub(full, full / 4, 0:4, null_rank_deficient = FALSE,
                                   adjusted_rank_deficient = FALSE))
  expect_false(grepl("RANK-DEFICIENT", healthy, fixed = TRUE))
  expect_false(grepl("no between-stratum variance at one time", healthy, fixed = TRUE))
  expect_false(grepl("are still reported", healthy, fixed = TRUE))
})

test_that("print() keeps its sentence where a block has no single zero to name", {
  # A quadratic block that lost a dimension is rank 2: its variance over time is not a
  # perfect square, so no time is named.
  S3 <- diag(c(0.09, 0.01, 0.001))
  dimnames(S3) <- list(c("(Intercept)", "wave", "I(wave^2)"),
                       c("(Intercept)", "wave", "I(wave^2)"))
  x3 <- long_pcv_stub(S3, S3 / 4, 0:4, null_slope_at_boundary = FALSE)
  out3 <- printed(x3)
  expect_match(out3, "RANK-DEFICIENT", fixed = TRUE)
  expect_match(out3, paste0("PCV_intercept and PCV(wave) are still reported: they are ",
                            "carried by the intercept variance, which stays estimable. ",
                            "A trajectory decomposition needs either more measurement ",
                            "occasions per stratum, or stratum_slope = FALSE"),
               fixed = TRUE)
  expect_false(grepl("at one time", out3, fixed = TRUE))

  # A linear block whose slope variance is zero is flat over time: nothing to name.
  flat <- rank1_block(c(0.3, 0))
  out_flat <- printed(long_pcv_stub(flat, flat / 4, 0:4))
  expect_match(out_flat, old_claim, fixed = TRUE)
  expect_false(grepl("at one time", out_flat, fixed = TRUE))

  # An object stored without its time centering cannot place the time on the user's
  # scale, so none is printed.
  Sn <- rank1_block(c(0.03, -0.10))
  x_old <- long_pcv_stub(Sn, rank1_block(c(0.15, -0.10)), 0:4)
  x_old$time_center <- NULL
  out_old <- printed(x_old)
  expect_match(out_old, old_claim, fixed = TRUE)
  expect_false(grepl("at one time", out_old, fixed = TRUE))
})

# ---- a fitted pair --------------------------------------------------------------

test_that("the time print() names is where the fitted pair's own variance is smallest", {
  skip_on_cran()
  # The bundled panel on a time axis that does not start at zero, so the fit is
  # centered internally; under estimation = "ML" its adjusted stratum block collapses.
  data(maihda_long_data, envir = environment())
  d <- maihda_long_data
  d$wave <- d$wave + 10
  a <- suppressMessages(suppressWarnings(
    maihda(wellbeing ~ wave + (1 | gender:ethnicity:education), data = d, id = "id",
           time = "wave", decomposition = "longitudinal", estimation = "ML")))
  p <- a$pcv
  skip_if_not(isTRUE(p$adjusted_rank_deficient) && !isTRUE(p$adjusted_at_boundary),
              "the adjusted block did not land on the boundary")
  expect_equal(p$time_center, 10)
  # Independent of the stored block and of the centering: the between-stratum variance
  # the PCV object reports at each time is a quadratic in time; find its lowest point.
  mid <- mean(p$pcv_t$time)
  tc <- p$pcv_t$time - mid
  co <- unname(stats::coef(stats::lm(p$pcv_t$var_adjusted ~ tc + I(tc^2))))
  lowest <- mid - co[2] / (2 * co[3])
  zero <- maihda_long_pcv_zero_time(p$Sigma_stratum_adjusted, p$time_center)
  expect_equal(zero, lowest, tolerance = 1e-6)
  # The variance there is nothing next to the variance at the reporting times.
  expect_lt(maihda_var_at_time(p$Sigma_stratum_adjusted, zero - p$time_center),
            1e-4 * min(p$pcv_t$var_adjusted))
  out <- printed(p)
  expect_match(out, sprintf("wave = %s for the adjusted model", formatC(zero, digits = 2, format = "f")),
               fixed = TRUE)
  expect_match(out, "from the reported times (10 to 14),", fixed = TRUE)
  if (!isTRUE(p$null_rank_deficient)) {
    expect_false(grepl("for the null model", out, fixed = TRUE))
  }
})
