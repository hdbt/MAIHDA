# Audit 2026-09-26 -- a design-weighted binomial fit must not be plotted on the
# link scale under Gaussian labels.
#
# FINDING. plot_prediction_deviation_panels(fit, type = "auto") -- the default, and
# what plot(fit, type = "all") uses -- routed a wemix (WeMix::mix, sampling-weighted)
# BINOMIAL fit through the GAUSSIAN branch and plotted per-stratum LOG-ODDS labelled
# "Fitted Value", where every other engine plotted probabilities under "Predicted
# Probability". Measured on the fixture below:
#
#   wemix  type = "auto"        range [-0.414975, 0.787588]   <- log-odds, negative
#   wemix  type = "binomial"    range [ 0.405296, 0.669213]   <- probabilities
#   lme4   type = "auto"        range [ 0.406591, 0.619677]   <- probabilities
#
# Nothing warned. The same model on the same data gave two different quantities on
# an axis whose label also differed, depending only on the engine.
#
# ROOT CAUSE. maihda_prediction_panel_auto_type() reads maihda_family(), which is
# stats::family() with a brmsfit fallback -- and stats::family() has no method for a
# WeMixResults, so it returned NULL and auto_type() fell through to its "gaussian"
# default. The maihda_model wrapper knew the family all along
# (maihda_model_family(fb)$family == "binomial"), but the panel unwraps the fit
# before asking.
#
# FIX, and why it sits at maihda_family() rather than at the one call site: a
# WeMixResults CARRIES its family object at attr(, "resp")$family -- a real stats
# family for both the LM (WeMixLMResp) and GLM (WeMixGLMResp) paths, and the same
# object WeMix's own predict.WeMixResults() reads to invert the link. Resolving it
# there fixes every caller at once, including a bare WeMixResults handed in without
# its wrapper, which maihda_model_family() cannot help.
#
# SCOPE (measured by tracing maihda_family() over the public API on wemix, ordinal
# and lme4 fits, not inferred). Exactly two call sites ever saw a NULL:
#   * maihda_prediction_panel_auto_type()  -- WeMixResults. THE DEFECT.
#   * maihda_count_vpc_note()              -- WeMixResults and clmm, but BENIGN:
#     it returns NULL for every non-count family anyway, and neither engine can fit
#     one (maihda_wemix_check_family() allows gaussian-identity and binomial-logit
#     only; the ordinal engine is cumulative). Pinned below so it stays benign.
# The longitudinal maihda_family() sites are unreachable here -- a longitudinal
# MAIHDA refuses wemix and ordinal outright -- and the discriminatory-accuracy
# sites read the wrapper's family first.
#
# THE ORDINAL ENGINE IS NOT AFFECTED: stats::family() is undefined for a clmm too,
# but auto_type() tests for the clmm CLASS before it ever asks for a family. A clmm
# carries a link and no family object, so maihda_family() is still NULL there and is
# left that way rather than synthesising one.
#
# CONCURRENT PIN, ELSEWHERE -- ALREADY UPDATED. A parallel pass measured this same
# misroute and pinned the WRONG behaviour on purpose, so that fixing it would be
# visible rather than accidental:
#
#   expect_identical(maihda_prediction_panel_auto_type(fb$model), "gaussian")
#   expect_true(is.null(maihda_family(fb$model)))
#   expect_identical(maihda_model_family(fb)$family, "binomial")
#   expect_true(any(panel_fitted(plot_prediction_deviation_panels(fb, type = "auto")) < 0))
#
# Those four have been flipped to the corrected routing ("binomial"; not NULL; still
# "binomial"; nothing < 0, and "auto" identical to the type = "binomial" panel). They
# live in test-audit-2026-09-26e.R on the branch that carries the wemix
# transformed-response panel fix -- renamed from -26c.R, which collided with a third
# pass's file of that name. That file NEEDS this one's R/utils_maihda.R: on its own
# branch those expectations fail 7/24, and only the two changes together are green.
# The rest of it, including its own sum(pb) == 2.974711829268, is unaffected -- the
# response-scale values did not move.
#
# BEHAVIOUR CHANGE. This changes what the panel plots, and its labels, for every
# wemix binomial fit. Of 180 measurements taken across wemix/lme4/ordinal fits at
# HEAD, 176 were numerically identical after the fix; the four that moved are the
# family lookup itself (both wemix families), the routing, and the panel pinned here.

audit_0926d_data <- function() {
  # The RNG stream matters: y/ly are drawn BEFORE b, so those lines stay even though
  # the binomial fit does not use them -- drop them and the fit differs.
  set.seed(1); n <- 600
  d <- data.frame(g1 = factor(sample(c("a", "b"), n, TRUE)),
                  g2 = factor(sample(c("p", "q", "r"), n, TRUE)),
                  x  = rnorm(n))
  d <- make_strata(d, c("g1", "g2"))$data
  st <- as.integer(factor(d$stratum))
  d$y <- exp(0.4 * d$x + rnorm(length(unique(st)), 0, 0.5)[st] +
               rnorm(n, 0, 0.4) + 1)
  d$w <- runif(n, 0.5, 2)
  d$ly <- log(d$y)
  d$b <- rbinom(n, 1, plogis(0.5 * d$x + rnorm(length(unique(st)), 0, 0.4)[st]))
  d
}

# The stratum panel is the SECOND panel; its $data$fitted is what the y axis shows.
audit_0926d_panel <- function(p) {
  stats::setNames(as.numeric(p[[2]]$data$fitted), as.character(p[[2]]$data$stratum))
}

audit_0926d_q <- function(e) suppressMessages(suppressWarnings(e))


test_that("a WeMixResults reports the family it was fitted with", {
  skip_on_cran()
  skip_if_not_installed("WeMix")
  d <- audit_0926d_data()
  q <- audit_0926d_q

  fb <- q(fit_maihda(b ~ x + (1 | stratum), data = d, sampling_weights = "w",
                     family = binomial()))
  fg <- q(fit_maihda(ly ~ x + (1 | stratum), data = d, sampling_weights = "w"))
  expect_identical(fb$engine, "wemix")
  expect_s3_class(fb$model, "WeMixResults")

  # The premise of the defect: stats::family() really has no method here, so a
  # brmsfit-style fallback is the only route to the family.
  expect_error(stats::family(fb$model))

  # ... and the fit really does carry one, on BOTH wemix paths.
  for (f in list(fb, fg)) {
    fam <- maihda_family(f$model)
    expect_s3_class(fam, "family")
    # it is the fitted object's own family, not the wrapper's ...
    expect_identical(fam$family, attr(f$model, "resp")$family$family)
    # ... and the two agree, so no comparability key moved.
    expect_identical(fam$family, maihda_model_family(f)$family)
    expect_identical(fam$link, maihda_model_family(f)$link)
  }
  expect_identical(maihda_family(fb$model)$family, "binomial")
  expect_identical(maihda_family(fb$model)$link, "logit")
  expect_identical(maihda_family(fg$model)$family, "gaussian")
  expect_identical(maihda_family(fg$model)$link, "identity")
  expect_identical(maihda_model_family_key(fb), "binomial(logit)")
  expect_identical(maihda_model_family_key(fg), "gaussian(identity)")

  # Not a blanket "anything without a family() method gets one": a clmm still has no
  # family object, and none is invented for it.
  skip_if_not_installed("ordinal")
  labs <- c("low", "mid", "high")
  d$yo <- factor(labs[1 + (d$x > -0.4) + (d$x > 0.7)], levels = labs, ordered = TRUE)
  fo <- q(fit_maihda(yo ~ x + (1 | stratum), data = d, family = "ordinal"))
  expect_s3_class(fo$model, "clmm")
  expect_null(maihda_family(fo$model))
  # and it never needed one here: the class test runs first
  expect_identical(maihda_prediction_panel_auto_type(fo$model), "ordinal")
})


test_that("the wemix binomial panel plots probabilities, not log-odds", {
  skip_on_cran()
  skip_if_not_installed("WeMix")
  d <- audit_0926d_data()
  q <- audit_0926d_q
  fb <- q(fit_maihda(b ~ x + (1 | stratum), data = d, sampling_weights = "w",
                     family = binomial()))

  # THE DEFECT: this returned "gaussian".
  expect_identical(maihda_prediction_panel_auto_type(fb$model), "binomial")

  pa <- audit_0926d_panel(q(plot_prediction_deviation_panels(fb, type = "auto")))
  pb <- audit_0926d_panel(q(plot_prediction_deviation_panels(fb, type = "binomial")))
  o <- names(pb)

  # "auto" now IS the binomial branch -- same values, same order.
  expect_identical(pa[o], pb[o])
  # ... and those values did NOT move: this is the response-scale sum measured at
  # HEAD from type = "binomial", which the fix must leave alone.
  expect_equal(sum(pb), 2.974711829268, tolerance = 1e-10)
  expect_length(pb, 6L)

  # Probabilities, not log-odds. At HEAD the "auto" panel went negative.
  expect_true(all(pa > 0 & pa < 1))
  expect_equal(unname(range(pb[o])), c(0.405296105, 0.669212562), tolerance = 1e-8)

  # The HEAD values are the weighted-mean LOG-ODDS. Recorded so the regression stays
  # recognisable, and to show the branches are genuinely different quantities: the
  # mean of plogis(eta) is not plogis(mean of eta), so the panel could not have been
  # rescued by transforming what it already plotted.
  head_logodds <- c(-0.414975106, -0.259623016, -0.208948971,
                    -0.033344581, 0.063654847, 0.787587946)
  expect_false(isTRUE(all.equal(sort(unname(pb)), sort(plogis(head_logodds)),
                                tolerance = 1e-12)))

  # An OUTSIDE oracle, not just the sibling branch: the stratum value is the
  # sampling-weighted mean of plogis(eta) over the stratum's rows, with eta built
  # from the wrapper's own linear predictor.
  eta <- maihda_wemix_linpred(fb, include_re = TRUE)
  w <- fb$data$w
  wm <- function(v) vapply(split(seq_along(v), as.character(fb$data$stratum)),
                           function(i) sum(w[i] * v[i]) / sum(w[i]), numeric(1))
  expect_equal(unname(pb[o]), unname(wm(plogis(eta))[o]), tolerance = 1e-8)
  # and it is NOT the weighted mean of eta, which is what HEAD plotted
  expect_equal(sort(unname(wm(eta))), sort(head_logodds), tolerance = 1e-8)

  # The labels say so too. The distribution title sits on the FIRST panel, the axis
  # label on the second, and even the size legend changes meaning with the branch.
  p_auto <- q(plot_prediction_deviation_panels(fb, type = "auto"))
  expect_identical(p_auto[[1]]$labels$title, "Distribution of Predicted Probabilities")
  expect_identical(p_auto[[2]]$labels$y, "Predicted Probability")
  expect_identical(p_auto[[2]]$labels$size, "|Deviance\nResidual|")

  # plot(type = "all") is the route most users reach this through: a named list of
  # ggplots, whose "prediction_deviation" entry is this same panel.
  all_p <- q(plot(fb, type = "all"))
  expect_true("prediction_deviation" %in% names(all_p))
  expect_identical(audit_0926d_panel(all_p$prediction_deviation)[o], pb[o])
  expect_identical(all_p$prediction_deviation[[2]]$labels$y, "Predicted Probability")
})


test_that("a bare WeMixResults, with no wrapper, routes on scale too", {
  skip_on_cran()
  skip_if_not_installed("WeMix")
  d <- audit_0926d_data()
  q <- audit_0926d_q
  fb <- q(fit_maihda(b ~ x + (1 | stratum), data = d, sampling_weights = "w",
                     family = binomial()))

  # Why this case matters: the fix reads the FITTED object's family, so it works
  # where maihda_model_family() cannot -- a WeMixResults handed in on its own. At HEAD
  # this panel was log-odds too.
  bare <- audit_0926d_panel(
    q(plot_prediction_deviation_panels(fb$model, data = fb$data, type = "auto")))
  expect_true(all(bare > 0 & bare < 1))
  expect_identical(
    q(plot_prediction_deviation_panels(fb$model, data = fb$data,
                                       type = "auto"))[[2]]$labels$y,
    "Predicted Probability")

  # ONLY the family is recovered, not the design. The sampling weights live on the
  # wrapper, so a bare fit aggregates its strata with UNIT weights and its panel is the
  # UNWEIGHTED mean of the same row probabilities -- pre-existing behaviour of
  # maihda_prediction_panel_prior_weights(), unchanged here, and pinned so the two
  # routes are not mistaken for interchangeable.
  expect_true(all(maihda_prediction_panel_prior_weights(NULL, fb$model, fb$data) == 1))
  p_row <- as.numeric(stats::predict(fb$model, newdata = fb$data, type = "response"))
  grp <- split(seq_len(nrow(fb$data)), as.character(fb$data$stratum))
  unwt <- vapply(grp, function(i) mean(p_row[i]), numeric(1))
  expect_equal(unname(bare[names(unwt)]), unname(unwt), tolerance = 1e-10)

  # ... and it therefore does NOT equal the wrapped, design-weighted panel.
  wrapped <- audit_0926d_panel(q(plot_prediction_deviation_panels(fb, type = "auto")))
  expect_false(isTRUE(all.equal(unname(bare[names(unwt)]),
                                unname(wrapped[names(unwt)]), tolerance = 1e-6)))
})


test_that("wemix and lme4 now label and scale the same panel the same way", {
  skip_on_cran()
  skip_if_not_installed("WeMix")
  d <- audit_0926d_data()
  q <- audit_0926d_q
  fb <- q(fit_maihda(b ~ x + (1 | stratum), data = d, sampling_weights = "w",
                     family = binomial()))
  fl <- q(fit_maihda(b ~ x + (1 | stratum), data = d, engine = "lme4",
                     family = binomial()))

  pw <- q(plot_prediction_deviation_panels(fb, type = "auto"))
  pl <- q(plot_prediction_deviation_panels(fl, type = "auto"))
  expect_identical(maihda_prediction_panel_auto_type(fl$model), "binomial")
  expect_identical(pw[[2]]$labels$y, pl[[2]]$labels$y)
  expect_identical(pw[[1]]$labels$title, pl[[1]]$labels$title)
  expect_identical(pw[[2]]$labels$size, pl[[2]]$labels$size)

  # Both on (0, 1) now; the lme4 range is its own measured value (a different
  # estimator on the same data, so the numbers differ -- only the SCALE must agree).
  lv <- audit_0926d_panel(pl)
  expect_equal(unname(range(lv)), c(0.406591, 0.619677), tolerance = 1e-5)
  expect_true(all(lv > 0 & lv < 1))
})


test_that("the wemix gaussian panel is untouched", {
  skip_on_cran()
  skip_if_not_installed("WeMix")
  d <- audit_0926d_data()
  q <- audit_0926d_q
  fg <- q(fit_maihda(ly ~ x + (1 | stratum), data = d, sampling_weights = "w"))

  # Resolving the family changed nothing here: the fit was routed to "gaussian" by
  # the no-family DEFAULT before, and is routed there by its actual family now.
  expect_identical(maihda_prediction_panel_auto_type(fg$model), "gaussian")
  p <- q(plot_prediction_deviation_panels(fg, type = "auto"))
  expect_identical(p[[2]]$labels$y, "Fitted Value")
  expect_identical(p[[1]]$labels$title, "Distribution of Fitted Values")
  expect_identical(p[[2]]$labels$size, "Deviation\nMagnitude")

  v <- audit_0926d_panel(p)
  eta <- maihda_wemix_linpred(fg, include_re = TRUE)
  w <- fg$data$w
  wm <- vapply(split(seq_along(eta), as.character(fg$data$stratum)),
               function(i) sum(w[i] * eta[i]) / sum(w[i]), numeric(1))
  expect_equal(unname(v[names(wm)]), unname(wm), tolerance = 1e-8)
  # an identity link, so "auto" and an explicit "gaussian" are the same panel
  expect_identical(v, audit_0926d_panel(
    q(plot_prediction_deviation_panels(fg, type = "gaussian"))))
})


test_that("reading the family neither mutates the fit nor accepts a non-family", {
  skip_on_cran()
  skip_if_not_installed("WeMix")
  d <- audit_0926d_data()
  q <- audit_0926d_q
  fb <- q(fit_maihda(b ~ x + (1 | stratum), data = d, sampling_weights = "w",
                     family = binomial()))

  # The family is read straight off the fit, so a caller that edits what it gets back
  # must not reach through into the stored object.
  before <- utils::capture.output(print(attr(fb$model, "resp")$family))
  got <- maihda_family(fb$model)
  got$family <- "TAMPERED"
  expect_identical(utils::capture.output(print(attr(fb$model, "resp")$family)), before)
  expect_identical(maihda_family(fb$model)$family, "binomial")
  expect_equal(maihda_linkinv(maihda_family(fb$model))(0), 0.5)

  # The guard is inherits(, "family"), NOT a bare non-NULL test: a malformed response
  # component falls back to NULL -- the old behaviour, which maihda_model_family() then
  # covers from the wrapper -- rather than yielding a half-family that would route the
  # panel on a family name with no linkinv behind it.
  fake <- fb$model
  attr(fake, "resp") <- list(family = list(family = "bogus"))
  expect_null(maihda_family(fake))
  attr(fake, "resp") <- NULL
  expect_null(maihda_family(fake))

  # Canonicalisation is a no-op for the only two families wemix can fit, so nothing is
  # renamed on the way through.
  expect_identical(maihda_normalize_family_name("binomial"), "binomial")
  expect_identical(maihda_normalize_family_name("gaussian"), "gaussian")
})


test_that("the other maihda_family() NULL site stays benign", {
  skip_on_cran()
  skip_if_not_installed("WeMix")
  d <- audit_0926d_data()
  q <- audit_0926d_q

  # maihda_count_vpc_note() was the only OTHER call site a wemix fit reached with a
  # NULL family. It returns NULL for every non-count family, so resolving the family
  # cannot change its answer -- but only because wemix refuses count families at all.
  # Pin both halves: if that gate is ever relaxed, this breaks rather than silently
  # dropping the approximation note.
  fb <- q(fit_maihda(b ~ x + (1 | stratum), data = d, sampling_weights = "w",
                     family = binomial()))
  fg <- q(fit_maihda(ly ~ x + (1 | stratum), data = d, sampling_weights = "w"))
  expect_null(maihda_count_vpc_note(fb))
  expect_null(maihda_count_vpc_note(fg))
  expect_error(
    q(fit_maihda(b ~ x + (1 | stratum), data = d, sampling_weights = "w",
                 family = poisson())),
    "gaussian\\(identity\\) and binomial\\(logit\\)")
})
