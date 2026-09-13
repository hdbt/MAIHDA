# Audit 2026-09-13 (second pass): the low-replicate warning described the wrong
# interval, and fired on the very n_boot the help page recommends.
#
# maihda_validate_bootstrap_args() is shared by two bootstraps that fail for
# different reasons. The VPC, PCV and importance bootstraps form a percentile
# interval from two tail quantiles, where a flat "below ~200 replications is
# unstable" is fair. summary(df_method = "bootstrap") forms a different thing:
# an interval symmetric about the estimate at a critical value taken from ONE
# |t*| order statistic, at the rank maihda_boot_critical_rank() returns. The
# shared warning called that a percentile interval and fired at n_boot = 199 --
# which the same help page recommends as "the usual compromise for a GLMM".
#
# What actually makes that cut-off unsteady is how few draws lie at or beyond it,
# and that depends on the level, not on n_boot alone: 199 draws put 10 beyond a
# 95% cut-off but only 2 beyond a 99% one, while 99 draws already put 10 beyond a
# 90% one. The check is now level-aware, with its own wording, and the threshold
# of ten draws is exactly what makes the conventional 99, 199 and 999 the
# smallest counts that clear the 10%, 5% and 1% levels. The percentile branch,
# and therefore every other caller, is untouched.

w_of <- function(expr) {
  w <- character(0)
  withCallingHandlers(expr, warning = function(cond) {
    w <<- c(w, conditionMessage(cond))
    invokeRestart("muffleWarning")
  })
  w
}

test_that("the percentile branch is untouched and is still the default", {
  # Every other caller (calculate_pcv, compare_maihda_groups, pcv_importance,
  # the app) takes the default, so the default must stay the percentile rule.
  expect_match(w_of(maihda_validate_bootstrap_args(199, 0.95)), "percentile interval")
  expect_match(w_of(maihda_validate_bootstrap_args(199, 0.95)), "below ~200 replications")
  expect_length(w_of(maihda_validate_bootstrap_args(200, 0.95)), 0L)
  expect_length(w_of(maihda_validate_bootstrap_args(1000, 0.99)), 0L)
  # ... and naming it explicitly gives the same warning as the default does.
  expect_identical(w_of(maihda_validate_bootstrap_args(50, 0.95, "percentile")),
                   w_of(maihda_validate_bootstrap_args(50, 0.95)))

  # The hard floor and the conf_level check are unchanged on both branches.
  expect_error(maihda_validate_bootstrap_args(9, 0.95), "whole number")
  expect_error(maihda_validate_bootstrap_args(9, 0.95, "studentized"), "whole number")
  expect_error(maihda_validate_bootstrap_args(199, 1.4, "studentized"), "conf_level")
  expect_error(maihda_validate_bootstrap_args(199, 0.95, "nonsense"))
})

test_that("the studentized branch counts draws beyond the cut-off, not replicates", {
  # The documented recommendation passes clean at the level it is documented for.
  expect_length(w_of(maihda_validate_bootstrap_args(199, 0.95, "studentized")), 0L)
  # ... and so do the other conventional counts at their own levels.
  expect_length(w_of(maihda_validate_bootstrap_args(99, 0.90, "studentized")), 0L)
  expect_length(w_of(maihda_validate_bootstrap_args(999, 0.99, "studentized")), 0L)

  # Ten draws beyond the cut-off is the threshold, and 99 / 199 / 999 are exactly
  # the smallest counts that reach it at the 10% / 5% / 1% levels: one draw fewer
  # warns. That correspondence is the whole justification for the conventional
  # counts, so pin it rather than the constants.
  for (cl in c(0.90, 0.95, 0.99)) {
    need <- c("0.9" = 99L, "0.95" = 199L, "0.99" = 999L)[[as.character(cl)]]
    expect_identical(maihda_boot_critical_rank(need, cl), 10L, info = as.character(cl))
    expect_lt(maihda_boot_critical_rank(need - 1L, cl), 10L)
    expect_length(w_of(maihda_validate_bootstrap_args(need, cl, "studentized")), 0L)
    expect_match(w_of(maihda_validate_bootstrap_args(need - 1L, cl, "studentized")),
                 sprintf("at least %d draws", need))
  }

  # The wording names the rank and the level instead of calling it a percentile
  # interval: what the reader needs in order to judge the number.
  msg <- w_of(maihda_validate_bootstrap_args(99, 0.95, "studentized"))
  expect_match(msg, "rank 5 of the 99 |t*| draws", fixed = TRUE)
  expect_match(msg, "95% fixed-effect interval", fixed = TRUE)
  expect_false(grepl("percentile", msg, fixed = TRUE))
  # 199 draws are ample at 95% and thin at 99%: the level, not the count, decides.
  expect_match(w_of(maihda_validate_bootstrap_args(199, 0.99, "studentized")),
               "rank 2 of the 199", fixed = TRUE)

  # Too few draws to reach the level at all leaves an unbounded interval, which
  # says so rather than reporting a rank of zero.
  unb <- w_of(maihda_validate_bootstrap_args(18, 0.95, "studentized"))
  expect_match(unb, "unbounded", fixed = TRUE)
  expect_match(unb, "at least 199 draws", fixed = TRUE)
  expect_identical(maihda_boot_critical_rank(18L, 0.95), 0L)

  # A caller forming both intervals gets both checks, not just the first.
  expect_length(w_of(maihda_validate_bootstrap_args(50, 0.99,
                                                    c("percentile", "studentized"))), 2L)

  # A level within a whisker of 1 needs more draws than an integer can hold. The
  # advice must still be a number, and the count must not be coerced -- doing that
  # raised R's own "NAs introduced by coercion" beside a message reading
  # "Use at least NA draws".
  extreme <- w_of(maihda_validate_bootstrap_args(50, 1 - 1e-10, "studentized"))
  expect_length(extreme, 1L)
  expect_match(extreme, "Use at least [0-9]+ draws")
  expect_false(grepl("NA", extreme, fixed = TRUE))
})

test_that("summary() applies the rule that matches the interval it forms", {
  skip_on_cran()
  set.seed(3)
  g <- expand.grid(d1 = factor(1:2), d2 = factor(1:2), d3 = factor(1:2))
  d <- g[rep(seq_len(8), each = 40), , drop = FALSE]
  d$stratum <- factor(rep(seq_len(8), each = 40))
  u <- stats::rnorm(8, 0, 0.5)
  d$y <- u[d$stratum] + stats::rnorm(nrow(d))
  rownames(d) <- NULL
  fit <- fit_maihda(y ~ d1 + (1 | stratum), data = d)
  boot_warn <- function(w) w[grepl("n_boot", w, fixed = TRUE)]

  # The finding: following the help page's own recommendation drew a warning.
  set.seed(1)
  expect_length(boot_warn(w_of(
    summary(fit, df_method = "bootstrap", n_boot = 199))), 0L)

  # A genuinely thin count still warns, in the studentized wording.
  set.seed(1)
  thin <- boot_warn(w_of(summary(fit, df_method = "bootstrap", n_boot = 12)))
  expect_match(thin, "fixed-effect interval", fixed = TRUE)
  expect_false(any(grepl("percentile", thin, fixed = TRUE)))

  # The VPC interval really is a percentile one, and keeps the old warning.
  set.seed(1)
  expect_match(boot_warn(w_of(summary(fit, bootstrap = TRUE, n_boot = 199))),
               "percentile interval", fixed = TRUE)

  # Asking for both at once applies both rules, and each only when it is due:
  # 199 draws are thin for the percentile interval and ample for the studentized
  # one at 95%, so exactly one warning should come back, the percentile one.
  set.seed(1)
  both <- boot_warn(w_of(summary(fit, bootstrap = TRUE, df_method = "bootstrap",
                                 n_boot = 199)))
  expect_length(both, 1L)
  expect_match(both, "percentile interval", fixed = TRUE)
})

test_that("the help page's stated rule matches the constant behind it", {
  skip_on_cran()
  man <- testthat::test_path("..", "..", "man")
  skip_if_not(dir.exists(man), "man/ is not available from this test run")
  rd <- gsub("[[:space:]]+", " ",
             paste(readLines(file.path(man, "summary.maihda_model.Rd"), warn = FALSE),
                   collapse = " "))
  # The page names a threshold and three counts that follow from it. The counts
  # are pinned by behaviour above; this ties the prose to the same constant, so a
  # change to one cannot leave the other standing.
  expect_true(grepl("ten draws beyond the cut-off", rd, fixed = TRUE))
  expect_true(grepl("level-aware", rd, fixed = TRUE))
  expect_identical(maihda_boot_critical_rank(99L, 0.90), 10L)
  expect_identical(maihda_boot_critical_rank(199L, 0.95), 10L)
  expect_identical(maihda_boot_critical_rank(999L, 0.99), 10L)
})
