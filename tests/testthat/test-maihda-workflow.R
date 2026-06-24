make_workflow_data <- function(seed = 4001, n = 360) {
  set.seed(seed)
  d <- data.frame(
    country = rep(c("A", "B", "C"), length.out = n),
    gender = sample(c("F", "M"), n, replace = TRUE),
    race = sample(c("X", "Y"), n, replace = TRUE),
    age = rnorm(n)
  )
  sk <- interaction(d$gender, d$race, drop = TRUE)
  d$y <- 1 + 0.3 * d$age + rnorm(nlevels(sk), sd = 0.7)[sk] + rnorm(n, sd = 0.4)
  d$bin <- rbinom(n, 1, plogis(-0.2 + 0.4 * d$age + rnorm(nlevels(sk), sd = 0.6)[sk]))
  d
}

test_that("maihda() returns a consistent bundle; groups NULL without group", {
  d <- make_workflow_data()
  a <- suppressMessages(maihda(y ~ age + (1 | gender:race), data = d))

  expect_s3_class(a, "maihda_analysis")
  expect_s3_class(a$model, "maihda_model")
  expect_s3_class(a$summary, "maihda_summary")
  expect_null(a$groups)

  # Overall VPC equals fitting + summarising directly
  m_direct <- fit_maihda(y ~ age + (1 | gender:race), data = d)
  expect_equal(a$summary$vpc$estimate, summary(m_direct)$vpc$estimate, tolerance = 1e-8)
})

test_that("the formula builders quote non-syntactic dimension names safely", {
  # Regression: maihda_adjusted_formula() (two-model path) and remove_terms()
  # (crossed-dimensions path) built RHS terms with sprintf("`%s`", .), which breaks
  # a legal name containing a backtick. maihda_quote_name() escapes it correctly.
  expect_identical(maihda_quote_name("a`b"), "`a\\`b`")
  expect_identical(maihda_quote_name("gender var"), "`gender var`")

  set.seed(4101)
  n <- 900
  gv <- sample(c("F", "M"), n, TRUE)
  rc <- sample(c("A", "B", "C"), n, TRUE)
  sk <- interaction(gv, rc, drop = TRUE)
  d <- data.frame(
    y = 2 + 0.8 * (gv == "M") + rnorm(nlevels(sk), sd = 1.2)[sk] + rnorm(n, sd = 1),
    `gender var` = gv, race = rc, check.names = FALSE
  )

  # Two-model path: the adjusted formula must add the backtick-quoted main effects.
  a2 <- suppressMessages(suppressWarnings(maihda(y ~ (1 | `gender var`:race), data = d)))
  expect_match(paste(deparse(a2$adjusted_formula), collapse = " "), "`gender var`",
               fixed = TRUE)
  expect_true(is.finite(a2$pcv$pvc))

  # Crossed-dimensions path strips the explicit main effects via remove_terms().
  a3 <- suppressMessages(suppressWarnings(
    maihda(y ~ `gender var` + race + (1 | `gender var`:race), data = d,
           decomposition = "crossed-dimensions")))
  expect_true(is.finite(a3$summary$vpc$estimate))
})

test_that("maihda() with group attaches a comparison equal to compare_maihda_groups()", {
  d <- make_workflow_data(4002)
  a <- suppressMessages(maihda(y ~ age + (1 | gender:race), data = d, group = "country"))

  expect_s3_class(a$groups, "maihda_group_comparison")
  expect_setequal(a$groups$group, c("A", "B", "C"))

  cmp_direct <- compare_maihda_groups(y ~ age + (1 | gender:race), d, group = "country")
  expect_equal(a$groups$vpc, cmp_direct$vpc, tolerance = 1e-8)
  expect_equal(a$groups$var_between, cmp_direct$var_between, tolerance = 1e-8)
})

test_that("maihda() forwards comparison-only args without leaking them into lmer", {
  d <- make_workflow_data(4007)
  # min_group_n and shared_strata are compare_maihda_groups args, not lmer args;
  # they must not be passed through to the model fitter.
  expect_no_error(
    a <- suppressMessages(maihda(y ~ age + (1 | gender:race), data = d, group = "country",
                min_group_n = 10, shared_strata = TRUE))
  )
  expect_s3_class(a$groups, "maihda_group_comparison")
})

test_that("maihda() auto-detects a binary outcome consistently for model and groups", {
  d <- make_workflow_data(4003)
  expect_warning(
    a <- suppressMessages(maihda(bin ~ age + (1 | gender:race), data = d, group = "country")),
    "binary", ignore.case = TRUE
  )
  expect_equal(a$model$family$family, "binomial")
  # The group comparison used the same resolved family (no silent gaussian fallback)
  expect_equal(attr(a$groups, "family"), "binomial")
})

test_that("print and summary methods work for maihda_analysis", {
  d <- make_workflow_data(4004)
  a <- suppressMessages(maihda(y ~ age + (1 | gender:race), data = d, group = "country"))

  expect_output(print(a), "MAIHDA Analysis")
  expect_output(print(a), "VPC/ICC")
  expect_output(print(a), "Group comparison by 'country'")

  s <- summary(a)
  expect_s3_class(s, "maihda_summary")
  expect_s3_class(attr(s, "groups"), "maihda_group_comparison")
})

test_that("plot.maihda_analysis dispatches to model and group plots", {
  d <- make_workflow_data(4005)
  a <- suppressMessages(maihda(y ~ age + (1 | gender:race), data = d, group = "country"))

  expect_s3_class(plot(a, type = "vpc"), "ggplot")
  expect_s3_class(plot(a, type = "predicted"), "ggplot")
  expect_s3_class(plot(a, type = "group_vpc"), "ggplot")
  expect_s3_class(plot(a, type = "group_components"), "ggplot")

  # group plots require a group argument
  a_nogroup <- suppressMessages(maihda(y ~ age + (1 | gender:race), data = d))
  expect_error(plot(a_nogroup, type = "group_vpc"), "group", ignore.case = TRUE)
})

test_that("compute_maihda_ternary_data is classed and plots via plot()", {
  skip_if_not_installed("ggtern")
  d <- make_workflow_data(4006)
  model <- fit_maihda(y ~ age + (1 | gender:race), data = d)

  td <- compute_maihda_ternary_data(model, verbose = FALSE)
  expect_s3_class(td, "maihda_ternary")
  expect_s3_class(plot(td), "ggplot")

  expect_warning(plot_maihda_ternary(td), "deprecated")
})

test_that("maihda() fits the adjusted model and reports a PCV", {
  d <- make_workflow_data(4101)
  a <- suppressMessages(maihda(y ~ age + (1 | gender:race), data = d))

  expect_s3_class(a$model_adjusted, "maihda_model")
  expect_s3_class(a$summary_adjusted, "maihda_summary")
  expect_s3_class(a$pcv, "pvc_result")
  expect_true(is.finite(a$pcv$pvc))

  # Adjusted formula = null formula + the dimensions' additive main effects;
  # the null model carries NO dimension main effects.
  adj_terms <- attr(stats::terms(reformulas::nobars(a$adjusted_formula)), "term.labels")
  null_terms <- attr(stats::terms(reformulas::nobars(a$formula)), "term.labels")
  expect_true(all(c("gender", "race") %in% adj_terms))
  expect_false(any(c("gender", "race") %in% null_terms))

  # PCV matches a manual calculate_pvc() on the two fitted models
  manual <- calculate_pvc(a$model, a$model_adjusted)
  expect_equal(a$pcv$pvc, manual$pvc, tolerance = 1e-8)
})

test_that("maihda() print and summary surface the PCV and adjusted model", {
  d <- make_workflow_data(4102)
  a <- suppressMessages(maihda(y ~ age + (1 | gender:race), data = d))

  expect_output(print(a), "PCV \\(null -> adjusted\\)")
  expect_output(print(a), "Adjusted formula")

  s <- summary(a)
  expect_s3_class(attr(s, "pcv"), "pvc_result")
  expect_s3_class(attr(s, "adjusted"), "maihda_summary")
})

test_that("maihda() enters an auto-binned numeric dimension as the tertile factor", {
  set.seed(4103)
  n <- 600
  d <- data.frame(gender = sample(c("F", "M"), n, TRUE), ses = rnorm(n))
  skb <- cut(d$ses, stats::quantile(d$ses, c(0, 1/3, 2/3, 1)), include.lowest = TRUE)
  ig <- interaction(d$gender, skb, drop = TRUE)
  d$math <- 1 + rnorm(nlevels(ig), sd = 0.8)[ig] + rnorm(n, sd = 0.5)

  a <- suppressWarnings(suppressMessages(maihda(math ~ 1 + (1 | gender:ses), data = d)))
  adj_terms <- attr(stats::terms(reformulas::nobars(a$adjusted_formula)), "term.labels")

  # make_strata() leaves the original `ses` numeric, so the additive main effect must
  # be the reconstructed tertile FACTOR (matching the strata), not the linear column.
  expect_true(".maihda_dim_ses" %in% adj_terms)
  expect_false("ses" %in% adj_terms)
  binned <- a$model_adjusted$original_data[[".maihda_dim_ses"]]
  expect_s3_class(binned, "factor")
  expect_equal(nlevels(binned), 3L)
  expect_true(is.finite(a$pcv$pvc))
})

test_that("maihda() errors with a single stratum dimension (no intersection)", {
  d <- make_workflow_data(4104)
  expect_error(
    maihda(y ~ age + (1 | gender), data = d),
    "at least two stratum dimensions"
  )
})

test_that("maihda() errors when the stratum dimensions are unidentifiable", {
  set.seed(4106)
  n <- 200
  d <- data.frame(stratum = factor(rep(seq_len(8), each = 25)), x = rnorm(n))
  d$y <- 1 + 0.4 * d$x + rnorm(8, sd = 0.6)[d$stratum] + rnorm(n, sd = 0.3)
  # A hand-built 'stratum' column records no dimensions, so the adjusted model cannot
  # be built; maihda() errors rather than returning a null-only result.
  expect_error(
    maihda(y ~ x + (1 | stratum), data = d),
    "could not be identified"
  )
})

test_that("plot.maihda_analysis routes decomposition types to the adjusted model", {
  d <- make_workflow_data(4105)
  a <- suppressMessages(maihda(y ~ age + (1 | gender:race), data = d))

  p_analysis <- plot(a, type = "effect_decomp")
  p_adjusted <- plot(a$model_adjusted, type = "effect_decomp",
                     summary_obj = a$summary_adjusted)
  expect_s3_class(p_analysis, "ggplot")
  # Same plot as drawing effect_decomp on the adjusted model directly
  expect_equal(p_analysis$labels$title, p_adjusted$labels$title)
  expect_equal(nrow(p_analysis$data), nrow(p_adjusted$data))
})

# Capture messages from an expression without letting them print.
maihda_capture_messages <- function(expr) {
  msgs <- character(0)
  val <- withCallingHandlers(
    expr,
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  list(value = val, messages = msgs)
}

test_that("maihda() accepts the dimensions' main effects in the formula (no PCV-of-0 bug)", {
  d <- make_workflow_data(4201)
  # Fully-specified adjusted form: the dimensions are written as additive fixed
  # effects. This must NOT collapse the null and adjusted into the same model (the old
  # behaviour gave PCV = 0 and a mislabelled VPC); maihda() derives the null by
  # dropping the dimension main effects.
  res <- maihda_capture_messages(
    maihda(y ~ age + gender + race + (1 | gender:race), data = d))
  b <- res$value
  expect_false(any(grepl("added the additive main effect", res$messages)))  # nothing auto-added

  null_terms <- attr(stats::terms(reformulas::nobars(b$formula)), "term.labels")
  adj_terms  <- attr(stats::terms(reformulas::nobars(b$adjusted_formula)), "term.labels")
  expect_false(any(c("gender", "race") %in% null_terms))   # null excludes the dimensions
  expect_true(all(c("gender", "race") %in% adj_terms))      # adjusted includes them
  expect_true(is.finite(b$pcv$pvc) && b$pcv$pvc != 0)

  # Equivalent to omitting the main effects (the legacy form): identical null VPC + PCV.
  a <- suppressMessages(maihda(y ~ age + (1 | gender:race), data = d))
  expect_equal(b$pcv$pvc, a$pcv$pvc, tolerance = 1e-8)
  expect_equal(b$summary$vpc$estimate, a$summary$vpc$estimate, tolerance = 1e-8)
})

test_that("maihda() adds the dimensions' main effects with a message when omitted", {
  d <- make_workflow_data(4202)
  expect_message(
    maihda(y ~ age + (1 | gender:race), data = d),
    "added the additive main effect"
  )
})

test_that("maihda() handles a partially-specified formula (some dimensions present)", {
  d <- make_workflow_data(4203)
  res <- maihda_capture_messages(
    maihda(y ~ age + gender + (1 | gender:race), data = d))
  cc <- res$value
  # The message names only the MISSING dimension (race), not the supplied one (gender).
  expect_true(any(grepl("race", res$messages)))
  expect_false(any(grepl("gender", res$messages)))

  null_terms <- attr(stats::terms(reformulas::nobars(cc$formula)), "term.labels")
  adj_terms  <- attr(stats::terms(reformulas::nobars(cc$adjusted_formula)), "term.labels")
  expect_false(any(c("gender", "race") %in% null_terms))   # null excludes ALL dimensions
  expect_true(all(c("gender", "race") %in% adj_terms))      # adjusted includes ALL of them
})

test_that("maihda() group decomposition is the same with or without dims in the formula", {
  d <- make_workflow_data(4204)
  g_dims   <- suppressMessages(maihda(y ~ age + gender + race + (1 | gender:race), data = d,
                                      group = "country", min_group_n = 10))
  g_nodims <- suppressMessages(maihda(y ~ age + (1 | gender:race), data = d,
                                      group = "country", min_group_n = 10))
  expect_equal(g_dims$groups$vpc, g_nodims$groups$vpc, tolerance = 1e-8)
  expect_equal(g_dims$groups$pcv, g_nodims$groups$pcv, tolerance = 1e-8)
})
