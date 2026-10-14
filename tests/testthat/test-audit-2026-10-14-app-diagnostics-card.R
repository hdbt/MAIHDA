# 2026-10-14 (packet 2026-10-14-app-diagnostics-card): the "Fit diagnostics" card of
# the Shiny app. Presentation only: no number, no wording and no print() output
# changes.
#
# The card turned every line of maihda_format_fit_diagnostics() into a bullet. Those
# lines are broken for a console, so one sentence was spread over two or three
# bullets ("shares with caution." stood alone on the app's default NHANES setup with
# decomposition = "crossed-dimensions"), and a convergence message showed as a bullet
# starting with "- ". Per the user: one bullet per diagnostic, its sentences as
# running text, the convergence messages listed under it.
#
# No block draws a random number: the fixtures are built from fixed sequences.

dc14_keep_seed <- function(expr) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    if (had) {
      assign(".Random.seed", old, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  }, add = TRUE)
  force(expr)
}

dc14_quiet <- function(expr) suppressMessages(suppressWarnings(expr))

# A diagnostics record with these parts.
dc14_record <- function(terms = character(0), messages = character(0),
                        adequacy = NULL) {
  list(engine = "lme4", singular = length(terms) > 0, singular_terms = terms,
       stratum_singular = if (length(terms) == 0L) NA else "stratum" %in% terms,
       converged = length(messages) == 0L, messages = messages, adequacy = adequacy)
}
dc14_cc <- list(dim_groups = c(gender = "gender", ses = "ses"),
                interaction_group = "stratum")
dc14_adequacy <- list(
  overdispersion = list(flag = TRUE, ratio = 4.2, p = 1e-8, family = "poisson"),
  zeroinflation = list(flag = TRUE, observed = 412, expected = 301.4, ratio = 1.37),
  re_normality = list(flag = TRUE, group = "stratum", term = "(Intercept)",
                      skew = 1.23, excess_kurtosis = 2.1, shapiro_p = 0.004),
  autocorrelation = list(flag = TRUE, acf1 = 0.41, n_pairs = 2400, gap = 1))

# Every word of the lines, in order: what the items must carry, nothing lost or added.
dc14_words <- function(x) unlist(strsplit(trimws(x), "\\s+"))
dc14_item_words <- function(items) {
  dc14_words(unlist(lapply(items, function(it) c(it$text, it$sub))))
}
dc14_line_words <- function(lines) dc14_words(sub("^\\s*- ", "", lines))

test_that("each diagnostic becomes one item with its sentences joined", {
  # A crossed-dimensions fit with a dimension at zero: three wrapped sentences.
  lines <- maihda_format_fit_diagnostics(dc14_record(c("ses", "gender")),
                                         cc_info = dc14_cc)
  items <- maihda_app_diagnostics_items(lines)
  expect_length(items, 1L)
  expect_identical(
    items[[1L]]$text,
    paste("Singular fit: the random-effects blocks for 'ses' and 'gender' are at (or near)",
          "the boundary (a variance estimated at ~0, or a perfect correlation between an",
          "intercept and a slope). 'ses' and 'gender' are stratum dimensions of this",
          "crossed-dimensions model: the additive variance has nothing from them, so read",
          "the additive and interaction shares with caution. The between-stratum variance",
          "and the VPC come from the terms not at the boundary."))
  expect_length(items[[1L]]$sub, 0L)

  # The sentence a plain fit breaks over two lines.
  lines <- maihda_format_fit_diagnostics(dc14_record("school"))
  items <- maihda_app_diagnostics_items(lines)
  expect_length(items, 1L)
  expect_match(items[[1L]]$text,
               paste("The 'stratum' block is NOT at the boundary, so the between-stratum",
                     "variance and the VPC/PCV are unaffected by this."), fixed = TRUE)

  # Convergence messages are the header's sub-items, without their "- ".
  msgs <- c("Model failed to converge with max|grad| = 0.0071 (tol = 0.002, component 1)",
            "Model is nearly unidentifiable: very large eigenvalue")
  items <- maihda_app_diagnostics_items(
    maihda_format_fit_diagnostics(dc14_record(messages = msgs)))
  expect_length(items, 1L)
  expect_identical(items[[1L]]$text, "Convergence warnings reported by lme4:")
  expect_identical(items[[1L]]$sub, msgs)

  # Each adequacy caveat is one item of two sentences.
  items <- maihda_app_diagnostics_items(
    maihda_format_fit_diagnostics(dc14_record(adequacy = dc14_adequacy)))
  expect_length(items, 4L)
  expect_identical(
    items[[2L]]$text,
    paste("Zero inflation: 412 observed zeros vs 301 expected (ratio 1.37). The",
          "likelihood underpredicts zeros; the VPC/PCV may be misleading."))
  expect_true(all(vapply(items, function(it) length(it$sub) == 0L, logical(1))))

  # All of them at once, in the formatter's order, and no word lost.
  lines <- maihda_format_fit_diagnostics(
    dc14_record("stratum", messages = msgs, adequacy = dc14_adequacy), cc_info = dc14_cc)
  items <- maihda_app_diagnostics_items(lines)
  expect_length(items, 6L)
  expect_match(items[[1L]]$text, "^Singular fit:")
  expect_identical(items[[2L]]$sub, msgs)
  expect_match(items[[3L]]$text, "^Overdispersion:")
  expect_identical(dc14_item_words(items), dc14_line_words(lines))

  # Nothing to show, nothing returned; an indented first line still starts an item.
  expect_length(maihda_app_diagnostics_items(character(0)), 0L)
  odd <- maihda_app_diagnostics_items(c("  first", "  second", "Next: x", "  - m"))
  expect_identical(odd[[1L]]$text, "first second")
  expect_identical(odd[[2L]]$text, "Next: x")
  expect_identical(odd[[2L]]$sub, "m")

  # The console print is untouched: still one line per line.
  out <- utils::capture.output(maihda_print_fit_diagnostics(
    dc14_record(messages = msgs)))
  expect_identical(out[1:4], c("Fit diagnostics:", "  Convergence warnings reported by lme4:",
                               paste0("    - ", msgs)))
})

# The app's objects, sourced as test-run_app.R does it.
dc14_source_app <- function() {
  app_file <- system.file("shiny", "app.R", package = "MAIHDA")
  if (app_file == "") {
    app_file <- file.path("inst", "shiny", "app.R")
  }
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  app_env <- new.env(parent = globalenv())
  suppressPackageStartupMessages(suppressWarnings(sys.source(app_file, envir = app_env)))
  app_env
}

# The card's HTML for a fitted model: the alert box and the bullets in it.
dc14_card <- function(app_env, res, model = res$model) {
  html <- NULL
  shiny::testServer(app_env$server, {
    session$setInputs(dataset = "pisa")
    model_results(model)
    summary_results(res$summary_obj)
    null_summary_results(res$summary_obj)
    decomposition_results(res$decomposition)
    decomposition_mode(res$decomposition_mode)
    fitted_family("gaussian")
    pcv_results(NULL)
    stepwise_results(NULL)
    session$flushReact()
    html <<- as.character(output$model_summary_ui$html)
  })
  start <- regexpr("<div class=\"alert alert-warning\">", html, fixed = TRUE)
  if (start < 0) {
    return(list(html = html, alert = NA_character_, li = NA_integer_))
  }
  alert <- substring(html, start)
  alert <- substring(alert, 1, regexpr("</div>", alert, fixed = TRUE) - 1L)
  li <- gregexpr("<li>", alert, fixed = TRUE)[[1]]
  list(html = html, alert = alert, li = if (li[1] < 0) 0L else length(li))
}

test_that("the app's card shows one bullet per diagnostic", {
  skip_on_cran()
  for (pkg in maihda_app_required_packages()) skip_if_not_installed(pkg)
  g <- expand.grid(d1 = paste0("a", 1:4), d2 = paste0("b", 1:5), stringsAsFactors = FALSE)
  mu <- c(-0.9, -0.2, 0.3, 0.8)[match(g$d1, paste0("a", 1:4))] +
    c(-0.7, -0.3, 0, 0.4, 0.6)[match(g$d2, paste0("b", 1:5))]
  d <- g[rep(seq_len(nrow(g)), each = 30), ]
  d$y <- rep(mu, each = 30) + rep(seq(-1.45, 1.45, length.out = 30), times = nrow(g))
  res <- dc14_quiet(maihda_app_fit_models(d, "y", c("d1", "d2"),
                                          decomposition = "crossed-dimensions"))
  expect_identical(res$model$diagnostics$singular_terms, "stratum")

  dc14_keep_seed({
    app_env <- dc14_source_app()

    # The singular fit: one bullet, its three sentences together.
    card <- dc14_card(app_env, res)
    expect_false(is.na(card$alert))
    expect_equal(card$li, 1L)
    expect_match(card$alert, paste(
      "is at (or near) the boundary (a variance estimated at ~0, or a perfect",
      "correlation between an intercept and a slope). 'stratum' is the",
      "interaction term of this crossed-dimensions model: the strata are estimated as",
      "additive (interaction share 0%). The between-stratum variance and the VPC come",
      "from the terms not at the boundary."), fixed = TRUE)

    # With two convergence messages: one more bullet, and the messages listed under
    # it without a leading "- ".
    nc <- res$model
    nc$diagnostics$converged <- FALSE
    nc$diagnostics$messages <- c("first message", "second message")
    card <- dc14_card(app_env, res, model = nc)
    expect_equal(card$li, 4L)
    expect_match(card$alert, paste0(
      "<li>\\s*Convergence warnings reported by lme4:\\s*<ul>\\s*",
      "<li>first message</li>\\s*<li>second message</li>\\s*</ul>\\s*</li>"))
    expect_false(grepl("<li>- ", card$alert, fixed = TRUE))

    # A fit without diagnostics: no box, the green line instead.
    clean <- res$model
    clean$diagnostics$singular <- FALSE
    clean$diagnostics$singular_terms <- character(0)
    card <- dc14_card(app_env, res, model = clean)
    expect_true(is.na(card$alert))
    expect_match(card$html, "Model converged with no singularity warnings.", fixed = TRUE)
  })
})
