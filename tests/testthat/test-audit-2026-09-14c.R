# Audit 2026-09-14c: plot(type = "effect_decomp") decomposed a prediction that still
# carried non-intersectional random effects -- CONFIRMED, and wider than the
# crossed-dimensions case it was reported against.
#
# plot_effect_decomposition() predicted EVERY random effect for its total and took the
# global mean from that prediction. A contextual (1 | site) from context = (or any
# grouping other than the stratum and, crossed-dimensions, the dimension effects)
# therefore:
#   * crossed-dimensions: stayed in the component labelled "Additive (dimension random
#     effects)", as each stratum's site composition. On the auditor's design (16 strata
#     crossed with 16 sites, the site mix following dimension a) that component spanned
#     -1.965..2.318 against -0.590..0.556 from the fit's own dimension effects, up to
#     1.810 of it site composition, although the same fit put the dimension variances at
#     0.063 and 0.078 against 2.703 for site; the strata were also ranked differently
#     from the stratum predictions behind plot(type = "predicted"), which already
#     excluded the context.
#   * two-model and crossed-dimensions alike: entered the global mean as the
#     row-weighted mean of the site effects, shifting every bar by that constant, so the
#     plotted deviations stopped averaging to zero -- 0.966 on this file's unbalanced
#     two-model design, 0.003 to 0.275 on six designs with lognormal site sizes and
#     site effects independent of size.
#
# FIX: the total prediction, and the global mean taken from it, is scoped to the
# intersectional random effects with the helpers the stratum predictions and the
# intersectional-scope AUC already use (maihda_da_re_scopes(),
# maihda_re_form_for_groups()), on the lme4 and brms paths. A fit without another
# grouping still predicts every random effect, and its plot data were identical() to
# the old ones on lme4 null / adjusted / crossed / binomial / weighted fits and on
# WeMix and clmm fits.

c14_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

c14_decomp <- function(p) {
  for (ly in p$layers) {
    if (!is.null(ly$data) && "additive_dev" %in% names(ly$data)) {
      return(as.data.frame(ly$data))
    }
  }
  stop("no stratum-level layer in the effect decomposition plot")
}

# Per-row pieces read straight off the lme4 fit -- fixef(), ranef() and the model
# frame -- so no expectation below runs through the package's own scoping helpers.
c14_rows <- function(fit) {
  fr <- fit@frame
  re <- lme4::ranef(fit)
  list(xb = as.numeric(lme4::getME(fit, "X") %*% lme4::fixef(fit)),
       u = lapply(stats::setNames(names(re), names(re)),
                  function(g) re[[g]][as.character(fr[[g]]), "(Intercept)"]),
       stratum = as.character(fr$stratum))
}

# The auditor's design: 16 strata (a x b) crossed with 16 sites; 80% of the rows sit
# in one of four sites tied to their level of a, and the site effects climb with a.
c14_crossed_data <- function() {
  set.seed(762)
  dg <- expand.grid(a = factor(1:4), b = factor(1:4))
  d <- dg[rep(seq_len(nrow(dg)), each = 60), ]
  primary <- (as.integer(d$a) - 1L) * 4L + sample(1:4, nrow(d), TRUE)
  d$site <- factor(ifelse(stats::runif(nrow(d)) < 0.8, primary,
                          sample(1:16, nrow(d), TRUE)))
  u_site <- rep(c(-2, -0.6, 0.6, 2), each = 4) + stats::rnorm(16, sd = 0.5)
  d$y <- u_site[d$site] + 0.15 * (as.integer(d$a) - 2.5) +
    0.2 * (as.integer(d$b) - 2.5) + stats::rnorm(nrow(d))
  d
}

# Twelve sites of 15 to 400 rows, the largest carrying the largest effects, so the
# row-weighted mean site effect is far from the (zero) mean over sites.
c14_unbalanced_data <- function() {
  set.seed(11)
  sizes <- round(exp(seq(log(15), log(400), length.out = 12)))
  site_eff <- sort(stats::rnorm(12, sd = 1.2))
  site <- factor(rep(seq_len(12), sizes))
  n <- length(site)
  d <- data.frame(site = site, a = sample(c("a1", "a2", "a3"), n, TRUE),
                  b = sample(c("b1", "b2"), n, TRUE), stringsAsFactors = FALSE)
  st <- interaction(d$a, d$b, drop = TRUE)
  d$y <- 1 + 0.4 * (d$a == "a2") - 0.3 * (d$b == "b2") +
    stats::rnorm(nlevels(st), sd = 0.4)[st] + site_eff[d$site] + stats::rnorm(n)
  d
}

test_that("the crossed-dimensions decomposition keeps a context effect out of the additive component", {
  skip_on_cran()
  a <- c14_quiet(maihda(y ~ 1 + (1 | a:b), data = c14_crossed_data(),
                        decomposition = "crossed-dimensions", context = "site"))
  m <- a$model
  expect_true(all(c("a", "b", "stratum", "site") %in% names(lme4::ranef(m$model))))

  pd <- c14_decomp(plot(a, type = "effect_decomp"))
  k <- as.character(pd$stratum)
  o <- c14_rows(m$model)
  scoped <- o$xb + o$u$a + o$u$b + o$u$stratum
  total <- tapply(scoped, o$stratum, mean) - mean(scoped)
  additive <- tapply(o$xb + o$u$a + o$u$b, o$stratum, mean) - mean(scoped)
  # (tapply() returns 1-d arrays; as.numeric() drops their dim as well as the names.)
  expect_equal(pd$additive_dev, as.numeric(additive[k]), tolerance = 1e-8)
  expect_equal(pd$total_dev, as.numeric(total[k]), tolerance = 1e-8)
  # The interaction component is still the stratum random effect itself.
  expect_equal(pd$intersectional_dev, as.numeric(tapply(o$u$stratum, o$stratum, mean)[k]),
               tolerance = 1e-8)

  # The fixture has teeth: the site composition the old plot folded into the additive
  # component exceeds 1 on the link scale for some stratum.
  site_mix <- tapply(o$u$site, o$stratum, mean) - mean(o$u$site)
  expect_gt(max(abs(site_mix)), 1)

  # The view now ranks the strata as the stratum predictions do.
  pr <- maihda_stratum_predictions_lme4(m, a$summary, scale = "link")
  expect_identical(k, as.character(pr$stratum)[order(pr$predicted_row)])
})

test_that("the two-model decomposition centres on the intersectional global mean", {
  skip_on_cran()
  a <- c14_quiet(maihda(y ~ 1 + (1 | a:b), data = c14_unbalanced_data(),
                        context = "site"))
  m <- a$model_adjusted
  pd <- c14_decomp(plot(a, type = "effect_decomp"))
  k <- as.character(pd$stratum)
  o <- c14_rows(m$model)

  g <- mean(o$xb + o$u$stratum)
  expect_equal(pd$additive_dev, as.numeric((tapply(o$xb, o$stratum, mean) - g)[k]),
               tolerance = 1e-8)
  n_j <- as.numeric(table(o$stratum)[k])
  expect_equal(sum(n_j * pd$total_dev) / sum(n_j), 0, tolerance = 1e-10)
  u <- stats::setNames(a$summary_adjusted$stratum_estimates$random_effect,
                       as.character(a$summary_adjusted$stratum_estimates$stratum))
  expect_equal(pd$intersectional_dev, unname(u[k]), tolerance = 1e-12)

  # Teeth: the row-weighted mean site effect the old global mean carried.
  expect_gt(abs(mean(o$u$site)), 0.5)
})

test_that("without another grouping the decomposition still uses every random effect", {
  # Negative control: passes before and after the fix.
  skip_on_cran()
  d <- c14_crossed_data()
  a <- c14_quiet(maihda(y ~ 1 + (1 | a:b), data = d[, c("a", "b", "y")],
                        decomposition = "crossed-dimensions"))
  pd <- c14_decomp(plot(a, type = "effect_decomp"))
  k <- as.character(pd$stratum)
  eta <- as.numeric(stats::predict(a$model$model, type = "link"))
  s <- as.character(a$model$data$stratum)
  u <- stats::setNames(a$summary$stratum_estimates$random_effect,
                       as.character(a$summary$stratum_estimates$stratum))
  total <- tapply(eta, s, mean) - mean(eta)
  expect_equal(pd$total_dev, as.numeric(total[k]), tolerance = 1e-12)
  expect_equal(pd$additive_dev, as.numeric(total[k] - u[k]), tolerance = 1e-12)
})

test_that("the brms decomposition requests the intersectional scope (Stan-free)", {
  skip_if_not_installed("brms")

  # posterior_linpred() stand-in: the posterior-mean linear predictor is the fixed part
  # plus the random effects of the groupings re_formula keeps (NULL = every grouping of
  # the model, NA = none), returned as two draws averaging to it.
  effects <- list(a = c(a1 = 0.3, a2 = -0.3), b = c(b1 = 0.1, b2 = -0.1),
                  stratum = c("1" = 0.05, "2" = -0.02, "3" = 0.04),
                  site = c(s1 = 1.5, s2 = -1.5))
  requested <- list()
  fake <- NULL
  local_mocked_bindings(
    posterior_linpred = function(object, re_formula = NULL, ...) {
      # list() wrapping, so a NULL re_formula is recorded rather than dropped.
      requested[length(requested) + 1L] <<- list(re_formula)
      groups <- if (is.null(re_formula)) {
        vapply(reformulas::findbars(fake$formula), function(b) deparse(b[[3]]), "")
      } else if (!inherits(re_formula, "formula") && is.na(re_formula)) {
        character(0)
      } else {
        vapply(reformulas::findbars(re_formula), function(b) deparse(b[[3]]), "")
      }
      eta <- rep(0.2, nrow(fake$data))
      for (g in groups) eta <- eta + unname(effects[[g]][fake$data[[g]]])
      rbind(eta + 0.01, eta - 0.01)
    },
    .package = "brms"
  )
  make <- function(formula, data, cc_info = NULL) {
    structure(list(
      model = structure(list(family = list(family = "gaussian", link = "identity")),
                        class = "brmsfit"),
      engine = "brms", formula = formula, data = data,
      family = list(family = "gaussian", link = "identity"),
      sampling_weights = NULL, cc_info = cc_info
    ), class = "maihda_model")
  }
  summ <- list(stratum_estimates = data.frame(
    stratum = names(effects$stratum), random_effect = unname(effects$stratum),
    stringsAsFactors = FALSE))
  groups_of <- function(f) vapply(reformulas::findbars(f), function(b) deparse(b[[3]]), "")

  # Crossed-dimensions with a context: stratum 1 sits in s1, stratum 3 in s2.
  d <- data.frame(stratum = c("1", "1", "2", "2", "3", "3"),
                  a = c("a1", "a1", "a1", "a1", "a2", "a2"),
                  b = c("b1", "b1", "b2", "b2", "b1", "b1"),
                  site = c("s1", "s1", "s1", "s2", "s2", "s2"),
                  stringsAsFactors = FALSE)
  fake <- make(y ~ 1 + (1 | a) + (1 | b) + (1 | stratum) + (1 | site), d,
               cc_info = list(dim_groups = c(a = "a", b = "b"),
                              interaction_group = "stratum", dim_labels = c("a", "b")))
  pd <- c14_decomp(plot_effect_decomposition(fake, summ))
  k <- as.character(pd$stratum)
  scoped <- c("1" = 0.65, "2" = 0.38, "3" = 0.04)          # 0.2 + u_a + u_b + u_stratum
  total <- scoped - mean(scoped)                             # equal stratum sizes
  expect_equal(pd$total_dev, unname(total[k]), tolerance = 1e-12)
  expect_equal(pd$additive_dev, unname(total[k] - effects$stratum[k]), tolerance = 1e-12)
  expect_setequal(groups_of(requested[[1]]), c("a", "b", "stratum"))
  expect_true(is.na(requested[[2]]))

  # Two-model with a context: five of the six rows sit in s1, so the row mean of the
  # site effects is (5 * 1.5 - 1.5) / 6 = 1, which must not enter the global mean.
  requested <- list()
  d2 <- d
  d2$site <- c("s1", "s1", "s1", "s1", "s1", "s2")
  fake <- make(y ~ 1 + (1 | stratum) + (1 | site), d2)
  pd <- c14_decomp(plot_effect_decomposition(fake, summ))
  k <- as.character(pd$stratum)
  g <- 0.2 + mean(effects$stratum)
  expect_equal(pd$additive_dev, rep(0.2 - g, 3), tolerance = 1e-12)
  expect_equal(mean(pd$total_dev), 0, tolerance = 1e-12)
  expect_setequal(groups_of(requested[[1]]), "stratum")

  # No other grouping: the total asks for every random effect, as before.
  requested <- list()
  fake <- make(y ~ 1 + (1 | stratum), d2)
  plot_effect_decomposition(fake, summ)
  expect_null(requested[[1]])
})
