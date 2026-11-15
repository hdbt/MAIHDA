# Audit 2026-11-15 (rdevel-one-level-factor): under R-devel, compare_maihda_groups() no
# longer said that a stratum dimension with a single level was why a group's PCV
# decomposition failed -- CONFIRMED. Wording only: every number, row and column is as
# it was. Reported by CRAN's checks: three failures of
# test-audit-2026-10-12-compare-failure-wording.R on the r-devel flavours.
#
# maihda_one_level_factor_error() told that error by its erroring call, `contrasts<-`,
# which model.matrix.default() calls to set a factor's contrasts. From R-devel r90878
# (2026-11-12) model.matrix.default() stops itself before that call:
#   cannot set contrasts for variable 'gender':
#     contrasts can be applied only to factors with 2 or more levels
# so the erroring call is model.matrix.default(...). Its first line is new and its
# second is the sentence `contrasts<-` gives; R translates both.
#
# An error raised by model.matrix.default() is that error when its message ends with
# the sentence as R words it in the session's language. An error raised by
# `contrasts<-` is one as before, whatever it says. The first and the last block build
# R-devel's condition by hand, so they test it on every R version.
#
# No block draws a random number: the fixtures are built from fixed sequences, and
# lme4's lFormula() and lmer() draw nothing.

rd15_sentence <- "contrasts can be applied only to factors with 2 or more levels"

# The sentence in the words R has for it in this session.
rd15_r_sentence <- function() gettext(rd15_sentence, domain = "R-stats")

# What model.matrix.default() raises for a one-level factor `variable` when it stops
# itself, with the call lme4 makes for its fixed-effect design.
rd15_devel_error <- function(variable, first = "cannot set contrasts for variable") {
  simpleError(paste0(first, " '", variable, "':\n  ", rd15_r_sentence()),
              quote(model.matrix.default(fixedform, fr, contrasts)))
}

# Groups g1 and g2 hold women only, g3 both genders; three ses levels; 30 rows a
# stratum, each stratum's outcomes its own mean plus the same fixed spread.
rd15_one_level_data <- function() {
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
    rep(seq(-1.45, 1.45, length.out = 30), times = nrow(d) / 30)
  rownames(d) <- NULL
  d
}

# compare_maihda_groups() on that data, with the adjusted model of a one-gender group
# stopping with `error` and every other fit left to fit_maihda().
rd15_compare <- function(error) {
  fit <- fit_maihda
  local_mocked_bindings(fit_maihda = function(formula, data, ...) {
    fixed <- attr(stats::terms(maihda_nobars(formula)), "term.labels")
    if ("gender" %in% fixed && length(unique(data$gender)) == 1L) {
      stop(error)
    }
    fit(formula, data, ...)
  })
  warns <- character(0)
  res <- withCallingHandlers(
    suppressMessages(compare_maihda_groups(y ~ (1 | gender:ses),
                                           data = rd15_one_level_data(),
                                           group = "country")),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  list(res = res,
       warning = grep("per-group PCV decomposition failed", warns, fixed = TRUE,
                      value = TRUE))
}

test_that("the one-level-factor error model.matrix.default() raises itself is told", {
  devel <- rd15_devel_error("gender")
  from <- conditionCall(devel)
  expect_true(maihda_one_level_factor_error(devel))
  # Whatever variable it names and however its first line is worded: only the
  # sentence is read.
  expect_true(maihda_one_level_factor_error(rd15_devel_error("a b:c")))
  expect_true(maihda_one_level_factor_error(
    rd15_devel_error("gender", first = "Kontraste lassen sich nicht setzen fuer")))
  # A namespace-qualified call names the same function.
  qualified <- as.call(list(call(":::", as.name("stats"),
                                 as.name("model.matrix.default")),
                            as.name("mt"), as.name("mf")))
  expect_true(maihda_one_level_factor_error(
    simpleError(conditionMessage(devel), qualified)))

  # model.matrix.default()'s other errors are not it.
  expect_false(maihda_one_level_factor_error(
    simpleError("model frame and formula mismatch in model.matrix()", from)))
  expect_false(maihda_one_level_factor_error(
    simpleError("'contrasts.arg' argument must be named", from)))
  # The sentence ends the message, in R's words for it.
  expect_false(maihda_one_level_factor_error(
    simpleError(paste0(conditionMessage(devel), "\n  and more"), from)))
  expect_false(maihda_one_level_factor_error(simpleError(
    "cannot set contrasts for variable 'gender':\n  a factor needs two levels", from)))
  # The same message from another function, or without a call, is not it.
  expect_false(maihda_one_level_factor_error(simpleError(
    conditionMessage(devel), quote(model.frame.default(formula = mt, data = d)))))
  expect_false(maihda_one_level_factor_error(simpleError(conditionMessage(devel))))
  # A message that is not one string is not it, and does not stop the check.
  for (msg in list(NULL, NA_character_, c(conditionMessage(devel), "second"), 1)) {
    odd <- structure(class = c("simpleError", "error", "condition"),
                     list(message = msg, call = from))
    expect_false(maihda_one_level_factor_error(odd))
  }

  # `contrasts<-` raising it is told by the call alone.
  earlier <- quote(`contrasts<-`(`*tmp*`, value = contr.funs[1 + isOF[nn]]))
  expect_true(maihda_one_level_factor_error(simpleError(rd15_r_sentence(), earlier)))
  expect_true(maihda_one_level_factor_error(simpleError("in any words", earlier)))
})

test_that("what R and lme4 stop with for a one-level factor is told", {
  skip_on_cran()
  caught <- function(expr) tryCatch(expr, error = function(e) e)
  # A factor, and a character column, which model.matrix() makes a factor.
  for (column in list(factor(rep("a", 4)), rep("a", 4))) {
    own <- caught(stats::model.matrix(~ f, data.frame(f = column)))
    expect_s3_class(own, "error")
    expect_true(maihda_one_level_factor_error(own))
  }
  # lme4 building the fixed-effect design of a model with such a main effect.
  d <- data.frame(y = seq(1, 12), f = "a", g = rep(c("u", "v", "w"), 4),
                  stringsAsFactors = FALSE)
  via_lme4 <- caught(lme4::lFormula(y ~ f + (1 | g), data = d))
  expect_s3_class(via_lme4, "error")
  expect_true(maihda_one_level_factor_error(via_lme4))

  # Other errors a design matrix stops with are not it.
  two <- data.frame(f = factor(c("a", "b", "a", "b")))
  unnamed <- caught(stats::model.matrix(~ f, two, contrasts.arg = list("contr.sum")))
  expect_s3_class(unnamed, "error")
  expect_false(maihda_one_level_factor_error(unnamed))
  absent <- caught(stats::model.matrix(~ nope, two))
  expect_s3_class(absent, "error")
  expect_false(maihda_one_level_factor_error(absent))
})

test_that("a group's single-level dimension is named when model.matrix.default() stops", {
  skip_on_cran()
  devel <- rd15_devel_error("gender")
  obs <- rd15_compare(devel)
  res <- obs$res
  expect_identical(res$status, rep("ok", 3))
  expect_identical(res$pcv_status[res$group %in% c("g1", "g2")], rep("failed", 2))
  expect_true(all(is.na(res$pcv[res$group %in% c("g1", "g2")])))
  expect_true(res$pcv_status[res$group == "g3"] %in% c("ok", "singular"))
  expect_true(is.finite(res$pcv[res$group == "g3"]))

  # R's message once for the two groups, then the dimension with one level there.
  expect_length(obs$warning, 1L)
  expect_match(obs$warning, paste0(
    "\ng1, g2: ", conditionMessage(devel),
    "\n  In g1, g2 the stratum dimension 'gender' has a single level, so its main ",
    "effect cannot be estimated there."), fixed = TRUE)
  expect_false(grepl("'ses'", obs$warning, fixed = TRUE))
  expect_false(grepl("g3", obs$warning, fixed = TRUE))

  # Another error of model.matrix.default() in the same groups is listed without it.
  other <- simpleError("model frame and formula mismatch in model.matrix()",
                       conditionCall(devel))
  obs_other <- rd15_compare(other)
  expect_identical(obs_other$res$pcv_status[obs_other$res$group %in% c("g1", "g2")],
                   rep("failed", 2))
  expect_length(obs_other$warning, 1L)
  expect_match(obs_other$warning, paste0("\ng1, g2: ", conditionMessage(other)),
               fixed = TRUE)
  expect_false(grepl("has a single level", obs_other$warning, fixed = TRUE))
})
