# Audit 2026-10-28 (single finding, series-pcv-basis): the variance basis named by
# print() of stepwise_pcv(), pcv_importance() and compare_maihda_groups().
#
# F01  For maximum-likelihood variances (glmer, wemix or ordinal fits, or lmer fits
#      made by maximum likelihood) the three prints called the basis "as fitted (REML
#      for Gaussian lmer, matching summary())" and gave no note, where calculate_pcv()
#      on the same pair names maximum likelihood and adds its note. The results record
#      ml_variances, and the prints use calculate_pcv()'s label and note without its
#      estimation = "REML" pointer; under maihda(group =) the note is printed once
#      (the user's decisions, 2026-10-28). A column subset of a comparison keeps the
#      recorded basis.

# Leaves R's random-number state as it found it (see setup-seed.R).
s28_keep_rng <- function(expr) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    if (had) {
      assign(".Random.seed", old, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  })
  force(expr)
}

s28_quiet <- function(expr) s28_keep_rng(suppressWarnings(suppressMessages(expr)))

# print() output without colour codes or trailing blanks.
s28_lines <- function(x) {
  sub("\\s+$", "", cli::ansi_strip(utils::capture.output(print(x))))
}

# The lines of the maximum-likelihood note in printed output: from its first line to
# the next blank line or the end.
s28_note_block <- function(out) {
  i <- grep("^Note: these are maximum-likelihood variances", out)
  if (length(i) == 0) return(character(0))
  j <- i[1]
  while (j < length(out) && nzchar(out[j + 1])) j <- j + 1
  out[i[1]:j]
}

# calculate_pcv()'s note where estimation = "REML" does not apply, as its print shows
# it for maximum-likelihood variances.
s28_pcv_note <- function() {
  p <- structure(list(pcv = 0.5, pvc = 0.5, var_model1 = 1, var_model2 = 0.5,
                      estimation = "fitted", estimation_used = "fitted",
                      adjusted_at_boundary = FALSE, ml_refit_failed = FALSE,
                      ml_variances = TRUE, reml_available = FALSE, bootstrap = FALSE),
                 class = c("pcv_result", "pvc_result"))
  s28_note_block(s28_lines(p))
}

# No two blank lines in a row.
s28_expect_single_blanks <- function(out, info = NULL) {
  blank <- out == ""
  expect_false(any(blank[-1] & blank[-length(blank)]), info = info)
}

# Two groups of 3 x 4 strata of 30, built without random numbers: a binary outcome y,
# a Gaussian z, an ordered o and design weights w.
s28_data <- function() {
  g <- expand.grid(a = factor(1:3), b = factor(1:4), grp = c("g1", "g2"), rep = 1:30,
                   stringsAsFactors = FALSE)
  ia <- as.integer(g$a)
  ib <- as.integer(g$b)
  ig <- match(g$grp, c("g1", "g2"))
  eta <- -0.4 + 0.5 * (ia - 2) + 0.3 * (ib - 2.5) + 0.6 * sin(3.1 * ia * ib)
  g$y <- as.integer(((g$rep * 37 + 11 * ig) %% 30) / 30 < stats::plogis(eta))
  g$z <- 1 + 0.6 * ia - 0.4 * ib + 0.5 * sin(2.3 * ia * ib) + sin(1.7 * seq_len(nrow(g)))
  g$o <- cut(g$z, stats::quantile(g$z, c(0, 1, 2, 3) / 3), include.lowest = TRUE,
             labels = c("lo", "mid", "hi"), ordered_result = TRUE)
  g$w <- 1 + (seq_len(nrow(g)) %% 3) / 2
  g
}

s28_strata <- function() make_strata(s28_data(), c("a", "b"))$data

ml_line <- "Variance basis: as fitted (maximum likelihood)"
reml_line <- "Variance basis: as fitted (REML for Gaussian lmer, matching summary())"

test_that("a glmer stepwise series names maximum-likelihood variances as calculate_pcv() does", {
  d <- s28_strata()
  sw <- s28_quiet(stepwise_pcv(d, "y", c("a", "b"), family = "binomial"))
  expect_true(attr(sw, "ml_variances"))

  # calculate_pcv() on the null model and the last step.
  null <- s28_quiet(fit_maihda(y ~ 1 + (1 | stratum), data = d, family = "binomial"))
  full <- s28_quiet(fit_maihda(y ~ a + b + (1 | stratum), data = d, family = "binomial"))
  p <- s28_quiet(calculate_pcv(null, full))
  expect_lt(abs(sw$Total_PCV[3] - p$pcv), 1e-10)
  expect_true(p$ml_variances)
  p$reml_available <- FALSE
  want <- s28_lines(p)

  out <- s28_lines(sw)
  basis <- grep("^Variance basis", out, value = TRUE)
  expect_identical(basis, ml_line)
  expect_identical(basis, grep("^Variance basis", want, value = TRUE))
  note <- s28_note_block(out)
  expect_identical(note, s28_note_block(want))
  expect_identical(note, s28_pcv_note())
  expect_length(note, 3L)
  # The note follows the basis line, one blank line apart, and names no REML option.
  i <- match(ml_line, out)
  expect_identical(out[i + 1], "")
  expect_identical(out[i + 2], note[1])
  expect_false(any(grepl("REML", out, fixed = TRUE)))
  s28_expect_single_blanks(out)

  # Before the boundary notes when a step is at the boundary.
  b <- sw
  attr(b, "boundary_steps") <- "b"
  out_b <- s28_lines(b)
  expect_lt(match(note[1], out_b), grep("^Note: the between-stratum variance is at", out_b))
  s28_expect_single_blanks(out_b)

  # A table saved without the record prints as it did.
  old <- sw
  attr(old, "ml_variances") <- NULL
  out_old <- s28_lines(old)
  expect_identical(grep("^Variance basis", out_old, value = TRUE), reml_line)
  expect_length(s28_note_block(out_old), 0L)
})

test_that("a glmer attribution names maximum-likelihood variances and notes them before the latent-scale caveat", {
  d <- s28_strata()
  imp <- s28_quiet(pcv_importance(d, "y", c("a", "b"), family = "binomial"))
  expect_true(imp$ml_variances)
  out <- s28_lines(imp)
  expect_identical(grep("^Variance basis", out, value = TRUE), ml_line)
  note <- s28_note_block(out)
  expect_identical(note, s28_pcv_note())
  i <- match(note[1], out)
  expect_identical(out[i - 1], "")
  expect_lt(grep("(efficiency).", out, fixed = TRUE), i)
  j <- i + length(note)
  expect_identical(out[j], "")
  expect_match(out[j + 1], "^Latent-scale caveat")
  expect_false(any(grepl("REML", out, fixed = TRUE)))
  s28_expect_single_blanks(out)

  # With no latent-scale caveat after it, the dominance block opens with its own blank
  # line.
  dom <- imp
  dom$method <- "dominance"
  dom$family_name <- "gaussian"
  dom$conditional <- matrix(c(0.1, 0.2, 0.3, 0.4), 2,
                            dimnames = list(c("a", "b"), c("size_0", "size_1")))
  dom$complete_dominance <- matrix(c(NA, FALSE, TRUE, NA), 2,
                                   dimnames = list(c("a", "b"), c("a", "b")))
  out_dom <- s28_lines(dom)
  k <- match(note[1], out_dom) + length(note)
  expect_identical(out_dom[k], "")
  expect_match(out_dom[k + 1], "^Conditional dominance")
  s28_expect_single_blanks(out_dom)

  # An object saved without the field prints as it did.
  old <- imp
  old$ml_variances <- NULL
  out_old <- s28_lines(old)
  expect_identical(grep("^Variance basis", out_old, value = TRUE), reml_line)
  expect_length(s28_note_block(out_old), 0L)
})

test_that("REML lmer series keep their label and no note, and an ML refit is labelled as one", {
  d <- s28_strata()
  sw <- s28_quiet(stepwise_pcv(d, "z", c("a", "b")))
  imp <- s28_quiet(pcv_importance(d, "z", c("a", "b")))
  expect_false(attr(sw, "ml_variances"))
  expect_false(imp$ml_variances)
  for (out in list(s28_lines(sw), s28_lines(imp))) {
    expect_identical(grep("^Variance basis", out, value = TRUE), reml_line)
    expect_length(s28_note_block(out), 0L)
  }

  # estimation = "ML" compares maximum-likelihood refits, which the basis line names;
  # the note is for the default basis only, as in calculate_pcv().
  sw_ml <- s28_quiet(stepwise_pcv(d, "z", c("a", "b"), estimation = "ML"))
  imp_ml <- s28_quiet(pcv_importance(d, "z", c("a", "b"), estimation = "ML"))
  expect_true(attr(sw_ml, "ml_variances"))
  expect_true(imp_ml$ml_variances)
  for (out in list(s28_lines(sw_ml), s28_lines(imp_ml))) {
    expect_identical(grep("^Variance basis", out, value = TRUE),
                     paste("Variance basis:", maihda_pcv_basis_label("ML")))
    expect_length(s28_note_block(out), 0L)
  }

  # One model with maximum-likelihood variances is enough, in either position: a REML
  # null against an adjusted model fitted by maximum likelihood.
  null_reml <- s28_quiet(fit_maihda(z ~ 1 + (1 | stratum), data = d))
  adj_ml <- s28_quiet(fit_maihda(z ~ a + b + (1 | stratum), data = d, REML = FALSE))
  expect_false(maihda_pcv_ml_variances(null_reml))
  expect_true(maihda_pcv_ml_variances(null_reml, adj_ml))
  expect_true(maihda_pcv_ml_variances(adj_ml, null_reml))
  expect_true(s28_quiet(calculate_pcv(null_reml, adj_ml))$ml_variances)
})

test_that("a group comparison names the basis and notes it below a table with PCVs", {
  d <- s28_data()
  cmp <- s28_quiet(compare_maihda_groups(y ~ 1 + (1 | a:b), data = d, group = "grp",
                                         family = "binomial"))
  expect_true(attr(cmp, "ml_variances"))
  out <- s28_lines(cmp)
  expect_identical(grep("^Variance basis", out, value = TRUE), ml_line)
  note <- s28_note_block(out)
  expect_identical(note, s28_pcv_note())
  expect_identical(utils::tail(out, length(note)), note)
  expect_identical(out[length(out) - length(note)], "")
  expect_false(any(grepl("REML", out, fixed = TRUE)))
  s28_expect_single_blanks(out)

  # calculate_pcv() on one group's pair names the same basis.
  rows <- d[d$grp == "g1", ]
  n1 <- s28_quiet(fit_maihda(y ~ 1 + (1 | a:b), data = rows, family = "binomial"))
  a1 <- s28_quiet(fit_maihda(y ~ a + b + (1 | a:b), data = rows, family = "binomial"))
  p1 <- s28_quiet(calculate_pcv(n1, a1))
  expect_lt(abs(cmp$pcv[cmp$group == "g1"] - p1$pcv), 1e-6)
  expect_identical(paste("Variance basis:", maihda_pcv_result_basis_label(p1)), ml_line)

  # Rows keep everything; without the pcv column the line stays and the note goes.
  expect_identical(s28_note_block(s28_lines(cmp[1, ])), note)
  cols <- cmp[, c("group", "n", "vpc")]
  expect_identical(attr(cols, "ml_variances"), TRUE)
  expect_identical(attr(cols, "estimation_used"), "fitted")
  out_cols <- s28_lines(cols)
  expect_identical(grep("^Variance basis", out_cols, value = TRUE), ml_line)
  expect_length(s28_note_block(out_cols), 0L)
  # So does a pcv column without values.
  no_values <- cmp
  no_values$pcv <- NA_real_
  out_na <- s28_lines(no_values)
  expect_identical(grep("^Variance basis", out_na, value = TRUE), ml_line)
  expect_length(s28_note_block(out_na), 0L)

  # A comparison saved without the record prints as it did.
  old <- cmp
  attr(old, "ml_variances") <- NULL
  out_old <- s28_lines(old)
  expect_identical(grep("^Variance basis", out_old, value = TRUE), reml_line)
  expect_length(s28_note_block(out_old), 0L)
})

test_that("crossed-dimensions, maximum-likelihood lmer and REML lmer comparisons", {
  d <- s28_data()
  # No PCV in a crossed-dimensions table, so the line alone.
  cc <- s28_quiet(compare_maihda_groups(y ~ 1 + (1 | a:b), data = d, group = "grp",
                                        family = "binomial",
                                        decomposition = "crossed-dimensions"))
  expect_true(attr(cc, "ml_variances"))
  out_cc <- s28_lines(cc)
  expect_identical(grep("^Variance basis", out_cc, value = TRUE), ml_line)
  expect_length(s28_note_block(out_cc), 0L)

  # lmer fits made by maximum likelihood through `...`.
  ml <- s28_quiet(compare_maihda_groups(z ~ 1 + (1 | a:b), data = d, group = "grp",
                                        REML = FALSE))
  expect_true(attr(ml, "ml_variances"))
  out_ml <- s28_lines(ml)
  expect_identical(grep("^Variance basis", out_ml, value = TRUE), ml_line)
  expect_identical(s28_note_block(out_ml), s28_pcv_note())
  # calculate_pcv() on one group's pair, fitted the same way, has the table's PCV and
  # prints the same basis line and note.
  rows <- d[d$grp == "g1", ]
  n1 <- s28_quiet(fit_maihda(z ~ 1 + (1 | a:b), data = rows, REML = FALSE))
  a1 <- s28_quiet(fit_maihda(z ~ a + b + (1 | a:b), data = rows, REML = FALSE))
  p1 <- s28_quiet(calculate_pcv(n1, a1))
  expect_lt(ml$pcv[ml$group == "g1"], 0.99)
  expect_lt(abs(ml$pcv[ml$group == "g1"] - p1$pcv), 1e-6)
  out_p1 <- s28_lines(p1)
  expect_identical(grep("^Variance basis", out_p1, value = TRUE), ml_line)
  expect_identical(s28_note_block(out_p1), s28_note_block(out_ml))

  reml <- s28_quiet(compare_maihda_groups(z ~ 1 + (1 | a:b), data = d, group = "grp"))
  expect_false(attr(reml, "ml_variances"))
  out_reml <- s28_lines(reml)
  expect_identical(grep("^Variance basis", out_reml, value = TRUE), reml_line)
  expect_length(s28_note_block(out_reml), 0L)

  # estimation = "ML": the PCVs compare maximum-likelihood refits of REML fits, as
  # calculate_pcv() records it; the basis line names the refit and there is no note.
  refit <- s28_quiet(compare_maihda_groups(z ~ 1 + (1 | a:b), data = d, group = "grp",
                                           estimation = "ML"))
  expect_true(attr(refit, "ml_variances"))
  used <- attr(refit, "estimation_used")
  expect_true(used %in% c("ML", "mixed"))
  out_refit <- s28_lines(refit)
  expect_identical(grep("^Variance basis", out_refit, value = TRUE),
                   paste("Variance basis:", maihda_pcv_basis_label(used)))
  expect_length(s28_note_block(out_refit), 0L)

  # A column subset keeps the basis actually used: a brms comparison's posterior
  # medians, a mixed ML comparison.
  keep <- c("group", "n", "vpc", "var_between", "pcv")
  post <- reml
  attr(post, "engine") <- "brms"
  attr(post, "estimation_used") <- "posterior"
  mixed <- reml
  attr(mixed, "estimation") <- "ML"
  attr(mixed, "estimation_used") <- "mixed"
  for (case in list(list(x = post, basis = "posterior"), list(x = mixed, basis = "mixed"))) {
    want <- paste("Variance basis:", maihda_pcv_basis_label(case$basis))
    expect_identical(grep("^Variance basis", s28_lines(case$x[, keep]), value = TRUE),
                     want, info = case$basis)
    expect_identical(grep("^Variance basis", s28_lines(case$x[keep]), value = TRUE),
                     want, info = case$basis)
    expect_identical(attr(case$x[1, keep], "estimation_used"), case$basis)
  }
})

test_that("maihda(group =) prints the maximum-likelihood note once", {
  d <- s28_data()
  an <- s28_quiet(maihda(y ~ a + b + (1 | a:b), data = d, family = "binomial",
                         group = "grp"))
  expect_true(an$pcv$ml_variances)
  expect_true(attr(an$groups, "ml_variances"))
  out <- s28_lines(an)
  starts <- grep("^Note: these are maximum-likelihood variances", out)
  expect_length(starts, 1L)
  g <- grep("^Group comparison by", out)
  expect_lt(starts, g)
  expect_identical(s28_note_block(out), s28_pcv_note())
  # The group table keeps its basis line; the PCV block names the basis by its note.
  expect_identical(sum(out == ml_line), 1L)
  expect_gt(match(ml_line, out), g)
  s28_expect_single_blanks(out)

  # Without a pooled PCV the group table carries the note.
  no_pcv <- an
  no_pcv$pcv <- NULL
  out2 <- s28_lines(no_pcv)
  starts2 <- grep("^Note: these are maximum-likelihood variances", out2)
  expect_length(starts2, 1L)
  expect_gt(starts2, grep("^Group comparison by", out2))
  expect_identical(s28_note_block(out2), s28_pcv_note())
  s28_expect_single_blanks(out2)
})

test_that("wemix and ordinal series name maximum-likelihood variances", {
  skip_on_cran()
  skip_if_not_installed("WeMix")
  skip_if_not_installed("ordinal")
  skip_slow_test()
  d <- s28_strata()
  raw <- s28_data()
  series <- list(
    wemix_stepwise = s28_quiet(stepwise_pcv(d, "z", c("a", "b"), sampling_weights = "w")),
    wemix_groups = s28_quiet(compare_maihda_groups(z ~ 1 + (1 | a:b), data = raw,
                                                   group = "grp", sampling_weights = "w")),
    ordinal_stepwise = s28_quiet(stepwise_pcv(d, "o", c("a", "b"))),
    ordinal_groups = s28_quiet(compare_maihda_groups(o ~ 1 + (1 | a:b), data = raw,
                                                     group = "grp")))
  expect_identical(attr(series$wemix_stepwise, "estimation_used"), "fitted")
  expect_identical(attr(series$wemix_groups, "engine"), "wemix")
  expect_identical(attr(series$ordinal_groups, "engine"), "ordinal")
  for (nm in names(series)) {
    x <- series[[nm]]
    expect_true(attr(x, "ml_variances"), info = nm)
    out <- s28_lines(x)
    expect_identical(grep("^Variance basis", out, value = TRUE), ml_line, info = nm)
    expect_identical(s28_note_block(out), s28_pcv_note(), info = nm)
  }
})
