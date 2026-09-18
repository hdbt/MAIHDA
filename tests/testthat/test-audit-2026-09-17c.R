# Audit 2026-09-17 (c): compare_maihda_groups() turned each of fit_maihda()'s
# up-front refusals into a failed fit, and a warning, for EVERY group.
#
# fit_maihda() refuses the lme4-only fitting arguments (weights=, subset=, offset=)
# for the wemix, brms and ordinal engines before it fits anything, and maihda() stops
# on the same message through its pooled fit. compare_maihda_groups() wrapped each
# per-group fit in tryCatch, so the identical refusal came back as a
# "fit failed: ..." status row plus one warning per group, and the call RETURNED a
# comparison table. Measured over 3 engines x 3 arguments (plus a sampling_weights +
# weights conflict): every combination returned with one warning per group where
# fit_maihda() and maihda() stopped once. Four sibling refusals behaved the same way:
# engine = "wemix" without sampling_weights, engine = "ordinal" with
# sampling_weights, engine = "ordinal" on a non-ordinal outcome, and an unsupported
# wemix family.
#
# The partial spellings were a second, wider defect. maihda_resolve_engine_dots()
# renames a partial spelling to the argument the ENGINE would bind it to, and the
# non-lme4 engines have no such argument to bind: WeMix::mix() has no subset/offset
# (and its weights is supplied by the package), brms::brm() has none of the three,
# and clmm() has no offset. So the refusals never saw them, and the engines then
# differed wildly: `offs = off` on the ordinal engine FITTED, silently without the
# offset, warning only that an "unknown control element" was ignored; WeMix stopped
# with "unused argument" and the whole vector deparsed into the message; brms stopped
# only after compiling its Stan model, 96 s in, with "passing unknown arguments".
#
# FIX (both scope choices the user's): the refusals move into shared helpers --
# maihda_refuse_engine_dots() and maihda_check_weights_conflict() -- which
# fit_maihda() and compare_maihda_groups() both call, so the messages cannot drift;
# compare_maihda_groups() calls them once, before any group is fitted, and also
# refuses the four sibling incompatibilities there. maihda_resolve_unbound_lme4_dots()
# renames a spelling the engine cannot bind but that abbreviates one of the three,
# so every engine refuses it with the normal message (brms `w = 500` still binds to
# warmup, `o` to opencl, and the lme4 path is untouched). Confinement: 56 resolved
# name sets and fits with legitimate engine arguments are identical() to the previous
# code.

c17_data <- function() {
  set.seed(20260917)
  n <- 360
  d <- data.frame(a = factor(sample(c("p", "q"), n, TRUE)),
                  b = factor(sample(c("r", "s", "t"), n, TRUE)),
                  grp = sample(c("g1", "g2"), n, TRUE))
  d$x <- stats::rnorm(n)
  d$w <- stats::runif(n, 0.5, 2)
  d$pw <- stats::runif(n, 0.5, 2)
  d$keep <- d$x > -1.6
  d$off <- stats::rnorm(n)
  sid <- as.integer(interaction(d$a, d$b))
  d$y <- 1 + 0.4 * d$x + stats::rnorm(6, sd = 0.6)[sid] + stats::rnorm(n)
  d$cnt <- stats::rpois(n, 2)
  d$yo <- factor(sample(c("lo", "mid", "hi"), n, TRUE),
                 levels = c("lo", "mid", "hi"), ordered = TRUE)
  d
}

# The error message, plus every warning raised on the way to it. Warnings are
# COUNTED, not suppressed: one per group is the symptom this pass is about.
c17_observe <- function(expr) {
  warns <- character(0)
  msg <- tryCatch(
    withCallingHandlers(
      suppressMessages({ force(expr); "returned without an error" }),
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        invokeRestart("muffleWarning")
      }),
    error = conditionMessage)
  list(msg = msg, warns = warns)
}

c17_gaussian <- y ~ x + (1 | a:b)
c17_ordered <- yo ~ x + (1 | a:b)

# engine settings x one forwarded argument, for the cases that need no engine package.
c17_cases <- function(d) {
  list(
    list(f = c17_gaussian, extra = list(engine = "wemix", sampling_weights = "w"),
         dots = list(subset = d$keep)),
    list(f = c17_gaussian, extra = list(engine = "wemix", sampling_weights = "w"),
         dots = list(offset = d$off)),
    list(f = c17_gaussian, extra = list(engine = "brms"), dots = list(weights = d$pw)),
    list(f = c17_gaussian, extra = list(engine = "brms"), dots = list(subset = d$keep)),
    list(f = c17_gaussian, extra = list(engine = "brms"), dots = list(offset = d$off)),
    list(f = c17_gaussian, extra = list(engine = "brms", sampling_weights = "w"),
         dots = list(subset = d$keep)),
    list(f = c17_ordered, extra = list(engine = "ordinal", family = "ordinal"),
         dots = list(weights = d$pw)),
    list(f = c17_ordered, extra = list(engine = "ordinal", family = "ordinal"),
         dots = list(subset = d$keep)),
    list(f = c17_ordered, extra = list(engine = "ordinal", family = "ordinal"),
         dots = list(offset = d$off)),
    # Design weights and precision weights together: refused before the engine guard.
    list(f = c17_gaussian, extra = list(engine = "wemix", sampling_weights = "w"),
         dots = list(weights = d$pw)),
    list(f = c17_gaussian, extra = list(engine = "brms", sampling_weights = "w"),
         dots = list(weights = d$pw))
  )
}

c17_fit_msg <- function(cs, d) {
  c17_observe(do.call(fit_maihda, c(list(cs$f, d), cs$extra, cs$dots)))$msg
}
c17_compare <- function(cs, d) {
  c17_observe(do.call(compare_maihda_groups,
                      c(list(cs$f, d, group = "grp"), cs$extra, cs$dots)))
}

test_that("compare_maihda_groups refuses lme4-only arguments once, as fit_maihda does", {
  d <- c17_data()
  cases <- c17_cases(d)
  # fit_maihda()'s own message for each case, taken before the mock below.
  expect_length(cases, 11)
  fit_msgs <- vapply(cases, c17_fit_msg, character(1), d = d)
  expect_true(all(grepl("not supported by engine|Supply either", fit_msgs)))

  # Reaching a per-group fit is the failure this pass fixes: make it visible.
  local_mocked_bindings(fit_maihda = function(...) stop("a fit was attempted"))
  for (i in seq_along(cases)) {
    obs <- c17_compare(cases[[i]], d)
    expect_identical(obs$msg, fit_msgs[[i]])
    expect_length(obs$warns, 0)
  }
})

test_that("a partial spelling the engine cannot bind is refused like the full name", {
  skip_on_cran()
  d <- c17_data()
  pkg_of <- c(wemix = "WeMix", brms = "brms", ordinal = "ordinal")
  engines <- names(pkg_of)[vapply(pkg_of, requireNamespace, logical(1), quietly = TRUE)]
  skip_if(length(engines) == 0, "no engine package installed")
  spec_of <- list(
    wemix = list(f = c17_gaussian, extra = list(engine = "wemix", sampling_weights = "w")),
    brms = list(f = c17_gaussian, extra = list(engine = "brms")),
    ordinal = list(f = c17_ordered, extra = list(engine = "ordinal", family = "ordinal"))
  )
  # partial spelling -> the argument it abbreviates, which the message must name.
  spellings <- list(subs = "subset", offs = "offset", weig = "weights")
  cases <- list()
  for (eng in engines) {
    for (tag in names(spellings)) {
      cases[[length(cases) + 1L]] <- list(
        f = spec_of[[eng]]$f, extra = spec_of[[eng]]$extra,
        dots = stats::setNames(list(switch(tag, subs = d$keep, offs = d$off, weig = d$pw)),
                               tag),
        named = spellings[[tag]])
    }
  }
  # Every fit_maihda() message first: local_mocked_bindings() below lasts for the
  # whole test, so a message taken after it would be the mock's.
  expect_length(cases, 3L * length(engines))
  fit_msgs <- vapply(cases, c17_fit_msg, character(1), d = d)
  for (i in seq_along(cases)) {
    expect_match(fit_msgs[[i]], "not supported by engine|Supply either")
    expect_match(fit_msgs[[i]], cases[[i]]$named, fixed = TRUE)
  }

  local_mocked_bindings(fit_maihda = function(...) stop("a fit was attempted"))
  for (i in seq_along(cases)) {
    obs <- c17_compare(cases[[i]], d)
    expect_identical(obs$msg, fit_msgs[[i]])
    expect_length(obs$warns, 0)
  }
})

test_that("the ordinal engine no longer fits an offs = spelling without the offset", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  d <- c17_data()
  # It used to fit every row with no offset at all, warning only about an unknown
  # control element.
  expect_error(
    suppressMessages(fit_maihda(c17_ordered, d, engine = "ordinal", family = "ordinal",
                                offs = d$off)),
    "Argument(s) not supported by engine = \"ordinal\": offset", fixed = TRUE)
  # The same fit without that argument is untouched.
  m <- suppressWarnings(suppressMessages(
    fit_maihda(c17_ordered, d, engine = "ordinal", family = "ordinal")))
  expect_identical(nrow(m$data), nrow(d))
})

test_that("compare_maihda_groups refuses the sibling engine incompatibilities up front", {
  d <- c17_data()
  siblings <- list(
    list(f = c17_gaussian, extra = list(engine = "wemix")),
    list(f = c17_ordered, extra = list(engine = "ordinal", family = "ordinal",
                                       sampling_weights = "w")),
    list(f = c17_gaussian, extra = list(engine = "ordinal")),
    list(f = c17_gaussian, extra = list(engine = "wemix", family = "poisson",
                                        sampling_weights = "w")),
    # A family spec neither engine can resolve fails here too, not per group.
    list(f = c17_gaussian, extra = list(family = "typo"))
  )
  expect_length(siblings, 5)
  fit_msgs <- vapply(siblings, function(cs) c17_fit_msg(c(cs, list(dots = NULL)), d),
                     character(1))
  expect_true(all(grepl("wemix|ordinal|Unsupported family", fit_msgs)))
  local_mocked_bindings(fit_maihda = function(...) stop("a fit was attempted"))
  for (i in seq_along(siblings)) {
    obs <- c17_compare(c(siblings[[i]], list(dots = NULL)), d)
    expect_identical(obs$msg, fit_msgs[[i]])
    expect_length(obs$warns, 0)
  }
})

test_that("several arguments, two spellings of one, and a look-alike name", {
  d <- c17_data()
  both <- c17_observe(fit_maihda(c17_gaussian, d, engine = "brms",
                                 subset = d$keep, offset = d$off))$msg
  # Both are named, in the canonical order, in one message.
  expect_match(both, "engine = \"brms\": subset, offset.", fixed = TRUE)
  local_mocked_bindings(fit_maihda = function(...) stop("a fit was attempted"))
  obs <- c17_compare(list(f = c17_gaussian, extra = list(engine = "brms"),
                          dots = list(subset = d$keep, offset = d$off)), d)
  expect_identical(obs$msg, both)
  expect_length(obs$warns, 0)

  # Renaming must not merge two spellings of one argument silently: both are named.
  dup <- c17_compare(list(f = c17_gaussian, extra = list(engine = "brms"),
                          dots = list(subset = d$keep, subs = d$keep)), d)
  expect_match(dup$msg, "supplied more than once", fixed = TRUE)
  expect_match(dup$msg, "'subset' and 'subs'", fixed = TRUE)

  # A name that merely starts like one of them is not an abbreviation of it.
  expect_identical(names(maihda_resolve_engine_dots(list(offset_col = 1), "brms",
                                                    "gaussian")), "offset_col")
})

test_that("design weights beside precision weights are refused before normalizing", {
  d <- c17_data()
  d$pw[1:3] <- 0
  # fit_maihda() refuses the conflict before it warns about the zero weights; so
  # must the wrapper, or a call it will never fit warns about its sample first.
  fit <- c17_observe(fit_maihda(c17_gaussian, d, engine = "wemix",
                                sampling_weights = "w", weights = d$pw))
  expect_match(fit$msg, "Supply either 'sampling_weights'", fixed = TRUE)
  expect_length(fit$warns, 0)
  local_mocked_bindings(fit_maihda = function(...) stop("a fit was attempted"))
  obs <- c17_compare(list(f = c17_gaussian,
                          extra = list(engine = "wemix", sampling_weights = "w"),
                          dots = list(weights = d$pw)), d)
  expect_identical(obs$msg, fit$msg)
  expect_length(obs$warns, 0)
})

test_that("the lme4 engine still takes weights, subset and offset per group", {
  d <- c17_data()
  # The refusal must not fire for lme4: the groups are still fitted.
  reached <- 0L
  local_mocked_bindings(fit_maihda = function(...) {
    reached <<- reached + 1L
    stop("a fit was attempted")
  })
  obs <- c17_compare(list(f = c17_gaussian, extra = list(),
                          dots = list(weights = d$pw, subset = d$keep, offset = d$off)), d)
  expect_identical(obs$msg, "returned without an error")
  expect_gt(reached, 0L)
  expect_true(any(grepl("fit failed", obs$warns, fixed = TRUE)))
})

test_that("a real lme4 group comparison with those arguments is unaffected", {
  skip_on_cran()
  d <- c17_data()
  cg <- suppressWarnings(suppressMessages(
    compare_maihda_groups(c17_gaussian, d, group = "grp",
                          weights = pw, subset = keep, offset = off)))
  expect_s3_class(cg, "data.frame")
  expect_setequal(as.character(cg$group), c("g1", "g2"))
  expect_true(all(cg$status == "ok"))
  expect_true(all(is.finite(cg$vpc)))
})

test_that("only a spelling the engine cannot bind is renamed", {
  skip_on_cran()
  one <- function(tag, value, engine, family = "gaussian") {
    names(maihda_resolve_engine_dots(stats::setNames(list(value), tag), engine, family))
  }
  skip_if_not_installed("brms")
  # Bound by brms: warmup and opencl, not weights/offset.
  expect_identical(one("w", 500, "brms"), "w")
  expect_identical(one("o", NULL, "brms"), "o")
  expect_identical(one("subs", TRUE, "brms"), "subset")
  expect_identical(one("weig", 1, "brms"), "weights")
  expect_identical(one("offs", 1, "brms"), "offset")
  expect_identical(one("iter", 100, "brms"), "iter")
  skip_if_not_installed("ordinal")
  expect_identical(one("offs", 1, "ordinal"), "offset")
  expect_identical(one("nAGQ", 5, "ordinal"), "nAGQ")
  expect_identical(one("thresh", "flexible", "ordinal"), "thresh")
  skip_if_not_installed("WeMix")
  expect_identical(one("weig", 1, "wemix"), "weights")
  expect_identical(one("nQuad", 13, "wemix"), "nQuad")
  # lme4 binds all three itself, so nothing here changes it.
  expect_identical(one("subs", TRUE, "lme4"), "subset")
  expect_identical(one("verb", 1, "lme4"), "verb")
})
