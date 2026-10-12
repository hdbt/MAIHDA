# Audit 2026-10-12 (packet 2026-10-12-compare-failure-wording): what
# compare_maihda_groups() says when a group's fit or its PCV decomposition fails --
# CONFIRMED. Wording only: every number, row and column is as it was.
#
# Three things were wrong with it.
#   1. The aggregated PCV-failure warning ended "A stratum dimension that is constant
#      within a group is the usual cause." after EVERY failure. That is the cause when
#      a dimension has one level among a group's rows, and not otherwise: a binomial
#      wemix fit stopped at a zero between-stratum variance got the sentence with two
#      levels of every dimension in every group, and so did an lme4 group whose null
#      between-stratum variance is zero.
#   2. Each group carried its own copy of the error message. The wemix zero-variance
#      message is three lines, so two groups failing alike gave a warning of 1,260
#      characters with the message twice; failed null fits raised one warning a group.
#   3. A failed null fit put the whole message in the `status` cell: 498 characters
#      over three lines there, and a 501-character line in the printed table.
#
# Per the user: the cause sentence is given only to a group whose error is a one-level
# factor, naming the group and the dimension; a reason several groups share is given
# once, with those groups, in the PCV-failure warning and in ONE warning for failed
# fits; and the status cell holds the first line of the message.
#
# The one-level case is told by its erroring call and not by its text. Every engine
# builds its design with model.matrix(), which stops in `contrasts<-`; the text is
# translated ("Kontraste koennen nur auf Faktoren ...", with an o-umlaut, on a German
# system). Which dimension it is comes from the strata the group's null fit populated.
#
# No block draws a random number: the fixtures are built from fixed sequences.

cf12_observe <- function(expr) {
  warns <- character(0)
  res <- withCallingHandlers(
    suppressMessages(expr),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  list(res = res, warns = warns)
}

# How often `what` occurs in `x`.
cf12_count <- function(what, x) {
  hits <- gregexpr(what, x, fixed = TRUE)[[1L]]
  if (identical(as.integer(hits[1L]), -1L)) 0L else length(hits)
}

cf12_spread <- seq(-1.45, 1.45, length.out = 30)
cf12_old_sentence <- "is the usual cause"
cf12_cause <- "has a single level"

# Groups g1 and g2 hold women only, g3 both genders; three ses levels; 30 rows a
# stratum, each stratum's outcomes its own mean plus the same fixed spread.
cf12_one_level_data <- function(gender_factor = FALSE) {
  mk <- function(country, genders) {
    g <- expand.grid(gender = genders, ses = c("low", "mid", "high"),
                     stringsAsFactors = FALSE)
    d <- g[rep(seq_len(nrow(g)), each = 30), ]
    d$country <- country
    d
  }
  d <- rbind(mk("g1", "F"), mk("g2", "F"), mk("g3", c("F", "M")))
  eff <- c(F.low = -0.5, F.mid = 0.1, F.high = 0.6, M.low = 0.3, M.mid = -0.4,
           M.high = 0.2)
  d$y <- 1 + unname(eff[paste(d$gender, d$ses, sep = ".")]) +
    rep(cf12_spread, times = nrow(d) / 30)
  if (gender_factor) d$gender <- factor(d$gender, levels = c("F", "M"))
  rownames(d) <- NULL
  d
}

# Group g1: every stratum holds the same 30 outcomes, so its null between-stratum
# variance is zero. Group g2: the strata differ. Both dimensions vary in both groups.
cf12_zero_null_data <- function() {
  g <- expand.grid(gender = c("F", "M"), ses = c("low", "mid", "high"),
                   stringsAsFactors = FALSE)
  d1 <- g[rep(seq_len(nrow(g)), each = 30), ]
  d1$y <- rep(cf12_spread, times = nrow(g))
  d1$country <- "g1"
  d2 <- d1
  d2$y <- d1$y + rep(c(-0.9, 0.2, 0.5, -0.4, 0.8, -0.2), each = 30)
  d2$country <- "g2"
  d <- rbind(d1, d2)
  rownames(d) <- NULL
  d
}

# Two groups of eight strata x `n` rows whose success counts are additive on the
# probability scale (12, 16, 16, 20, 16, 20, 20, 24 of 40): with n = 40 the adjusted
# model's variance is zero and the null model's is not; with n = 20 the null model's
# is zero too. A binomial wemix fit stops at a zero variance.
cf12_wemix_data <- function(n) {
  g <- expand.grid(d1 = 0:1, d2 = 0:1, d3 = 0:1)
  s <- round(n * (0.3 + 0.1 * (g$d1 + g$d2 + g$d3)))
  one <- g[rep(1:8, each = n), ]
  one$y <- unlist(lapply(1:8, function(j) rep(c(1L, 0L), c(s[j], n - s[j]))))
  for (v in c("d1", "d2", "d3")) one[[v]] <- factor(one[[v]])
  one$w <- 1
  d <- rbind(cbind(one, grp = "g1"), cbind(one, grp = "g2"))
  rownames(d) <- NULL
  d
}
cf12_diagnosis <- "estimated the between-stratum variance of this binomial model at zero"

test_that("a status cell keeps the first line of an error message", {
  three <- paste0("engine = \"wemix\": the diagnosis.\n",
                  "  What a zero variance means.\n",
                  "  engine = \"brms\" returns a fit.")
  expect_identical(maihda_first_line(three), "engine = \"wemix\": the diagnosis.")
  # A one-line message is the cell as it always was.
  expect_identical(maihda_first_line("Matrix seems negative semi-definite"),
                   "Matrix seems negative semi-definite")
  # Blank leading lines are passed over, and the line comes without its indentation.
  expect_identical(maihda_first_line("\n\n   indented first line  \nsecond"),
                   "indented first line")
  expect_identical(maihda_first_line("first\r\nsecond"), "first")
  expect_identical(maihda_first_line(""), "")
  # Bytes that are not text in their declared encoding cannot be split: the message
  # is kept as it came, without a warning from the splitting.
  bad <- "first\xff\nsecond"
  Encoding(bad) <- "UTF-8"
  expect_silent(kept <- maihda_first_line(bad))
  expect_identical(kept, bad)
})

test_that("a one-level factor is told by its erroring call, not its text", {
  words <- "contrasts can be applied only to factors with 2 or more levels"
  from <- quote(`contrasts<-`(`*tmp*`, value = contr.funs[1 + isOF[nn]]))

  # R's condition as it arrives, in English and as a German system prints it.
  expect_true(maihda_one_level_factor_error(simpleError(words, from)))
  expect_true(maihda_one_level_factor_error(simpleError(
    "Kontraste koennen nur auf Faktoren mit 2 oder mehr Stufen angewendet werden",
    from)))
  # A namespace-qualified call names the same function.
  qualified <- as.call(list(call("::", as.name("stats"), as.name("contrasts<-")),
                            as.name("x")))
  expect_true(maihda_one_level_factor_error(simpleError("any text", qualified)))

  # The same words from anywhere else are not it, nor is an error without a call.
  expect_false(maihda_one_level_factor_error(
    simpleError(words, quote(model.matrix.default(mt, mf)))))
  expect_false(maihda_one_level_factor_error(simpleError(words, quote(contrasts(x)))))
  expect_false(maihda_one_level_factor_error(simpleError(words)))

  # What R itself stops with for a one-level factor. If this stops holding, the
  # design-matrix code raises that error elsewhere and the cause is no longer told.
  own <- tryCatch(
    stats::model.matrix(~ f, data.frame(f = factor(rep("a", 3)))),
    error = function(e) e)
  expect_s3_class(own, "error")
  expect_true(maihda_one_level_factor_error(own))
})

test_that("the single-level dimensions are read from the strata a fit populated", {
  info <- data.frame(stratum = 1:4, label = c("a", "b", "c", "d"),
                     n = c(30L, 30L, 30L, 0L),
                     gender = c("F", "F", "F", "M"),
                     ses = c("low", "mid", "high", "low"),
                     stringsAsFactors = FALSE)
  model <- list(strata_info = info, strata_vars = c("gender", "ses"))
  # The stratum with the second gender holds no analytic row, so it does not count.
  expect_identical(maihda_single_level_dims(model), "gender")

  info$n[4] <- 12L
  expect_identical(
    maihda_single_level_dims(list(strata_info = info,
                                  strata_vars = c("gender", "ses"))),
    character(0))

  # Both dimensions constant; an auto-binned one is held as its bins.
  one <- info[1, ]
  one$age <- factor("Low", levels = c("Low", "Mid", "High"))
  expect_identical(
    maihda_single_level_dims(list(strata_info = one,
                                  strata_vars = c("gender", "ses", "age"))),
    c("gender", "ses", "age"))

  # Nothing to read from: no cause is named.
  expect_identical(maihda_single_level_dims(list(strata_info = NULL,
                                                 strata_vars = c("gender", "ses"))),
                   character(0))
  expect_identical(maihda_single_level_dims(list(strata_info = info[, c("stratum",
                                                                        "gender")],
                                                 strata_vars = "gender")),
                   character(0))
})

test_that("a reason several groups share is listed once, with those groups", {
  reasons <- list(g1 = "same words", g2 = "other words", g3 = "same words", g4 = "")
  expect_identical(maihda_failure_entries(reasons),
                   c("g1, g3: same words", "g2: other words", "g4: unknown error"))

  # A message of several lines stays whole, once.
  long <- "first line.\n  second line."
  expect_identical(maihda_failure_entries(list(a = long, b = long)),
                   "a, b: first line.\n  second line.")

  # The cause sentence goes to the groups it was found for, and names them and the
  # dimension. Groups that differ in the dimension are listed apart.
  one_level <- list(g1 = "words", g2 = "words", g3 = "words", g4 = "words")
  entries <- maihda_failure_entries(
    one_level, list(g1 = "gender", g2 = "gender", g3 = c("gender", "ses"),
                    g4 = character(0)))
  expect_identical(entries, c(
    paste0("g1, g2: words\n  In g1, g2 the stratum dimension 'gender' has a single ",
           "level, so its main effect cannot be estimated there."),
    paste0("g3: words\n  In g3 the stratum dimensions 'gender', 'ses' have a single ",
           "level each, so their main effects cannot be estimated there."),
    "g4: words"))
})

test_that("failed fits: the first line in status and one warning for all groups", {
  d <- cf12_zero_null_data()
  three <- paste0("the engine's diagnosis.\n",
                  "  What it means.\n",
                  "  What to use instead.")

  # Every group fails alike.
  local_mocked_bindings(fit_maihda = function(...) stop(three, call. = FALSE))
  obs <- cf12_observe(compare_maihda_groups(y ~ (1 | gender:ses), data = d,
                                            group = "country"))
  expect_identical(obs$res$status, rep("fit failed: the engine's diagnosis.", 2))
  expect_false(any(grepl("\n", obs$res$status, fixed = TRUE)))
  expect_true(all(is.na(obs$res$vpc)))

  failed <- grep("the model fit failed for group(s)", obs$warns, fixed = TRUE,
                 value = TRUE)
  expect_length(failed, 1L)
  expect_match(failed, "for group(s) g1, g2.", fixed = TRUE)
  expect_match(failed, paste0("\ng1, g2: ", three), fixed = TRUE)
  expect_identical(cf12_count("What it means.", failed), 1L)
  # One warning in all: none per group.
  expect_false(any(grepl("Group 'g", obs$warns, fixed = TRUE)))

  # Groups failing for different reasons each keep their own.
  local_mocked_bindings(fit_maihda = function(formula, data, ...) {
    stop("reason of ", data$country[1], call. = FALSE)
  })
  two <- cf12_observe(compare_maihda_groups(y ~ (1 | gender:ses), data = d,
                                            group = "country"))
  expect_identical(two$res$status, c("fit failed: reason of g1",
                                     "fit failed: reason of g2"))
  failed2 <- grep("the model fit failed for group(s)", two$warns, fixed = TRUE,
                  value = TRUE)
  expect_length(failed2, 1L)
  expect_match(failed2, "\ng1: reason of g1\ng2: reason of g2", fixed = TRUE)
})

test_that("a one-level dimension is named as the cause, for its groups only", {
  for (as_factor in c(FALSE, TRUE)) {
    d <- cf12_one_level_data(gender_factor = as_factor)
    obs <- cf12_observe(compare_maihda_groups(y ~ (1 | gender:ses), data = d,
                                              group = "country"))
    res <- obs$res
    info <- paste("gender as factor:", as_factor)
    expect_identical(res$status, rep("ok", 3), info = info)
    expect_identical(res$pcv_status[res$group %in% c("g1", "g2")],
                     rep("failed", 2), info = info)
    expect_true(res$pcv_status[res$group == "g3"] %in% c("ok", "singular"), info = info)
    expect_true(is.finite(res$pcv[res$group == "g3"]), info = info)

    w <- grep("per-group PCV decomposition failed", obs$warns, fixed = TRUE,
              value = TRUE)
    expect_length(w, 1L)
    expect_match(w, "failed for group(s) g1, g2.", fixed = TRUE, info = info)
    # The two groups share R's message, whatever language it is in: one entry.
    expect_identical(cf12_count("\ng1, g2: ", w), 1L, info = info)
    expect_match(w, paste0("\n  In g1, g2 the stratum dimension 'gender' has a single ",
                           "level, so its main effect cannot be estimated there."),
                 fixed = TRUE, info = info)
    expect_false(grepl("'ses'", w, fixed = TRUE), info = info)
    expect_false(grepl("g3", w, fixed = TRUE), info = info)
    expect_false(any(grepl(cf12_old_sentence, obs$warns, fixed = TRUE)), info = info)
  }
})

test_that("a zero null variance is not blamed on a constant dimension", {
  d <- cf12_zero_null_data()
  obs <- cf12_observe(compare_maihda_groups(y ~ (1 | gender:ses), data = d,
                                            group = "country"))
  res <- obs$res
  expect_identical(res$status, c("ok", "ok"))
  expect_identical(res$pcv_status, c("failed", "ok"))
  expect_lt(res$var_between[1], 1e-8)

  w <- grep("per-group PCV decomposition failed", obs$warns, fixed = TRUE, value = TRUE)
  expect_length(w, 1L)
  # The reason is calculate_pcv()'s own, and nothing is said about dimensions.
  own <- tryCatch(
    suppressMessages(suppressWarnings(calculate_pcv(
      fit_maihda(y ~ (1 | gender:ses), data = d[d$country == "g1", ]),
      fit_maihda(y ~ gender + ses + (1 | gender:ses), data = d[d$country == "g1", ])))),
    error = function(e) conditionMessage(e))
  expect_match(w, paste0("\ng1: ", own), fixed = TRUE)
  expect_false(grepl(cf12_cause, w, fixed = TRUE))
  expect_false(any(grepl(cf12_old_sentence, obs$warns, fixed = TRUE)))
})

test_that("wemix groups stopped at a zero variance share one copy of the message", {
  skip_on_cran()
  skip_if_not_installed("WeMix")
  f <- y ~ 1 + (1 | d1:d2:d3)

  # Adjusted fits at zero: the PCV fails in both groups, for the same reason.
  d <- cf12_wemix_data(40)
  obs <- cf12_observe(compare_maihda_groups(f, data = d, group = "grp",
                                            engine = "wemix", family = "binomial",
                                            sampling_weights = "w"))
  res <- obs$res
  expect_identical(res$status, c("ok", "ok"))
  expect_identical(res$pcv_status, c("failed", "failed"))
  expect_true(all(is.na(res$pcv)))
  # The null fits are untouched: each group's VPC is the one its own fit gives.
  one <- suppressMessages(suppressWarnings(summary(fit_maihda(
    f, data = d[d$grp == "g1", ], engine = "wemix", family = "binomial",
    sampling_weights = "w"))))
  expect_equal(res$vpc, rep(one$vpc$estimate, 2), tolerance = 1e-10)

  w <- grep("per-group PCV decomposition failed", obs$warns, fixed = TRUE, value = TRUE)
  expect_length(w, 1L)
  expect_match(w, "failed for group(s) g1, g2.", fixed = TRUE)
  expect_identical(cf12_count(cf12_diagnosis, w), 1L)
  expect_identical(cf12_count("engine = \"brms\" with 'sampling_weights'", w), 1L)
  # Every dimension has two levels in both groups: no cause is named.
  expect_false(grepl(cf12_cause, w, fixed = TRUE))
  expect_false(any(grepl(cf12_old_sentence, obs$warns, fixed = TRUE)))

  # Null fits at zero: the cell holds the diagnosis line, the warning all three.
  d20 <- cf12_wemix_data(20)
  obs20 <- cf12_observe(compare_maihda_groups(f, data = d20, group = "grp",
                                              engine = "wemix", family = "binomial",
                                              sampling_weights = "w"))
  err <- tryCatch(
    suppressMessages(suppressWarnings(fit_maihda(
      f, data = d20[d20$grp == "g1", ], engine = "wemix", family = "binomial",
      sampling_weights = "w"))),
    error = function(e) conditionMessage(e))
  expect_match(err, cf12_diagnosis, fixed = TRUE)
  lines <- strsplit(err, "\n", fixed = TRUE)[[1]]
  expect_gt(length(lines), 1L)
  expect_identical(obs20$res$status, rep(paste0("fit failed: ", lines[1]), 2))
  failed <- grep("the model fit failed for group(s)", obs20$warns, fixed = TRUE,
                 value = TRUE)
  expect_length(failed, 1L)
  expect_match(failed, paste0("\ng1, g2: ", err), fixed = TRUE)
  expect_identical(cf12_count(cf12_diagnosis, failed), 1L)
})
