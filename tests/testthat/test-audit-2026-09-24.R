# Regression tests for the 2026-09-24 audit finding.
#
#   F01 [P1] maihda_prepare_prediction_data() named an unseen stratum combination by
#            its display LABEL, in the same namespace as the internal stratum ids
#            ("1", "2", ...) every engine looks the random effect up by. A
#            one-dimension label such as "1" -- or "12" under sep = "" -- therefore
#            IS a fitted id, so allow_new_levels = TRUE silently returned that
#            stratum's random effect instead of the documented zero-effect
#            (fixed-effects-only) prediction: measured on one 20-stratum fixture the
#            lme4 value was 1.5080 against a correct 2.1651, and WeMix, ordinal and
#            brms were wrong by their own stratum effects the same way. Unseen
#            combinations now get generated ids checked against every id the model
#            carries. The supplied-'stratum' contradiction check was defeated by the
#            same collision (it compared a label against ids), including with the
#            default allow_new_levels = FALSE, and now rejects it. Classifying a row
#            as unseen also moved from the strata TABLE to the strata the FIT used,
#            so a combination whose rows all left the analytic sample -- for which
#            WeMix and ordinal silently gave the zero-effect prediction by default --
#            is refused like any other unseen stratum.

# 20 strata labelled A..T, so the ids 1..20 and the labels occupy disjoint value
# sets EXCEPT where a prediction row names a numeric-looking level.
maihda_f01_data <- function() {
  set.seed(1)
  J <- 20; n_per <- 15
  a <- rep(LETTERS[1:J], each = n_per)
  u <- stats::rnorm(J, sd = 1)
  x <- stats::rnorm(J * n_per)
  y <- 2 + 0.5 * x + u[match(a, LETTERS[1:J])] + stats::rnorm(J * n_per, sd = 1)
  d <- data.frame(y = y, x = x, a = a, w = stats::runif(J * n_per, 0.5, 2),
                  stringsAsFactors = FALSE)
  d$yo <- cut(d$y, stats::quantile(d$y, 0:4 / 4), include.lowest = TRUE,
              labels = c("lo", "mid", "hi", "top"), ordered_result = TRUE)
  d
}

# ---- the id generator and its inputs (no fit) --------------------------------

test_that("generated unseen-stratum ids avoid every id the model carries", {
  f <- MAIHDA:::maihda_unseen_stratum_ids

  # One id per distinct key, the same id wherever a key repeats.
  ids <- f(c("k1", "k2", "k1"), taken = as.character(1:20))
  expect_length(ids, 3L)
  expect_identical(ids[1], ids[3])
  expect_false(ids[1] == ids[2])
  expect_false(any(ids %in% as.character(1:20)))

  # A prefix alone guarantees nothing: when the drawn ids are themselves taken the
  # generator must draw again rather than hand back a colliding value.
  first <- f("k1", taken = character())
  again <- f(c("k1", "k2"), taken = c(first, paste0(".", first)))
  expect_false(any(again %in% c(first, paste0(".", first))))
  expect_identical(anyDuplicated(again), 0L)

  # Still collision-free when the whole first batch is taken.
  blocked <- f(c("k1", "k2", "k3"), taken = f(c("k1", "k2", "k3"), taken = character()))
  expect_false(any(blocked %in% f(c("k1", "k2", "k3"), taken = character())))

  # paste0() recycles a zero-length argument to "", so the empty case is explicit.
  expect_identical(f(character(), taken = character()), character())
})

test_that("maihda_taken_stratum_ids covers fitted, declared and table ids", {
  obj <- list(
    data = data.frame(stratum = factor(c("1", "2"), levels = c("1", "2", "3")),
                      stringsAsFactors = FALSE),
    strata_info = data.frame(stratum = c(1L, 2L, 3L, 4L),
                             label = c("A", "B", "C", "D"),
                             stringsAsFactors = FALSE))
  taken <- MAIHDA:::maihda_taken_stratum_ids(obj)
  expect_true(all(c("1", "2", "3", "4") %in% taken))   # fitted, declared, table-only
  expect_identical(anyDuplicated(taken), 0L)
  # A model with no stratum column at all yields the table ids alone, not an error.
  expect_identical(MAIHDA:::maihda_taken_stratum_ids(list(strata_info = obj$strata_info)),
                   c("1", "2", "3", "4"))
})

test_that("combination keys separate combinations that share a display label", {
  f <- MAIHDA:::maihda_stratum_combination_keys
  sep <- " x "
  # Two different combinations whose values contain the separator paste to the same
  # label ("a x b x c"), so the label cannot identify them; the keys must.
  d <- data.frame(d1 = c("a x b", "a"), d2 = c("c", "b x c"), stringsAsFactors = FALSE)
  labs <- MAIHDA:::maihda_stratum_labels(d, c("d1", "d2"), sep)
  expect_identical(labs[1], labs[2])
  keys <- f(d, c("d1", "d2"))
  expect_false(keys[1] == keys[2])
  # Identical rows share a key.
  same <- f(data.frame(d1 = c("a", "a"), d2 = c("b", "b"), stringsAsFactors = FALSE),
            c("d1", "d2"))
  expect_identical(same[1], same[2])
  # A numeric dimension is auto-binned first, so two values in one bin share a key.
  info <- list(v = list(breaks = c(0, 50, 100), labels = c("lo", "hi")))
  binned <- f(data.frame(g = c("F", "F", "F"), v = c(10, 20, 80),
                         stringsAsFactors = FALSE), c("g", "v"), info)
  expect_identical(binned[1], binned[2])
  expect_false(binned[1] == binned[3])
})

# ---- lme4 end to end ---------------------------------------------------------

test_that("an unseen combination labelled like a stratum id predicts at a zero effect", {
  skip_on_cran()
  d <- maihda_f01_data()
  m <- fit_maihda(y ~ x + (1 | a), data = d)

  # "1" and "7" are fitted stratum ids (labels A and G); as dimension VALUES they
  # are combinations the model never saw.
  expect_true(all(c("1", "7") %in% maihda_known_strata(m)))
  expect_false(any(c("1", "7") %in% d$a))

  nd <- data.frame(x = 0, a = c("1", "new", "7"), stringsAsFactors = FALSE)
  p <- predict_maihda(m, newdata = nd, allow_new_levels = TRUE)

  # Oracle: lme4's own prediction with the random effects dropped.
  fixed_only <- unname(stats::predict(m$model, newdata = data.frame(x = 0),
                                      re.form = NA))
  expect_equal(unname(p), rep(fixed_only, 3), tolerance = 1e-12)

  # What the collision was worth: the borrowed effects are materially non-zero, so
  # the assertion above could not pass by accident.
  re <- lme4::ranef(m$model)$stratum
  expect_gt(abs(re["1", 1]), 0.3)
  expect_gt(abs(re["7", 1]), 0.3)
  # The value returned is not the one the collision produced (fixed + that stratum's
  # effect), and the two are far apart.
  expect_gt(abs(unname(p[1]) - unname(fixed_only + re["1", 1])), 0.3)
  expect_gt(abs(unname(p[3]) - unname(fixed_only + re["7", 1])), 0.3)

  # A seen stratum still carries its own effect (the fix must not zero everything).
  seen <- predict_maihda(m, newdata = data.frame(x = 0, a = "A", stringsAsFactors = FALSE))
  expect_equal(unname(seen), unname(fixed_only + re["1", 1]), tolerance = 1e-12)
  expect_gt(abs(seen - fixed_only), 0.3)

  # Rejected by default, and the message names the LABELS the user supplied, not the
  # internal ids they are given.
  expect_error(predict_maihda(m, newdata = nd), "not present when the model was fit")
  expect_error(predict_maihda(m, newdata = nd), "fit: 1, new, 7.", fixed = TRUE)
})

test_that("prepared unseen rows carry ids no fitted stratum uses, one per combination", {
  skip_on_cran()
  d <- maihda_f01_data()
  m <- fit_maihda(y ~ x + (1 | a), data = d)

  nd <- data.frame(x = 0, a = c("1", "2", "new", "1", "A"), stringsAsFactors = FALSE)
  prep <- MAIHDA:::maihda_prepare_prediction_data(m, nd, allow_new_levels = TRUE)
  new_ids <- prep$stratum[1:4]

  # The invariant: no unseen row may be given an id the model carries.
  expect_false(any(new_ids %in% MAIHDA:::maihda_taken_stratum_ids(m)))
  # One id per distinct combination: rows 1 and 4 are the same combination.
  expect_identical(new_ids[1], new_ids[4])
  expect_identical(anyDuplicated(new_ids[1:3]), 0L)
  # The seen row keeps its training id.
  expect_identical(prep$stratum[5], "1")

  # Several new combinations at once each get the fallback.
  p <- predict_maihda(m, newdata = nd, allow_new_levels = TRUE)
  fixed_only <- unname(stats::predict(m$model, newdata = data.frame(x = 0), re.form = NA))
  expect_equal(unname(p[1:4]), rep(fixed_only, 4), tolerance = 1e-12)
  expect_false(isTRUE(all.equal(unname(p[5]), fixed_only)))
})

test_that("a custom separator cannot make a label collide with a stratum id", {
  skip_on_cran()
  set.seed(2)
  d2 <- expand.grid(p = c("x", "y", "z"), q = c("r", "s", "t", "v"), rep = 1:12,
                    stringsAsFactors = FALSE)
  d2$xx <- stats::rnorm(nrow(d2))
  u <- stats::rnorm(12)[as.integer(factor(paste(d2$p, d2$q)))]
  d2$y <- 1 + 0.3 * d2$xx + u + stats::rnorm(nrow(d2))
  st <- make_strata(d2, vars = c("p", "q"), sep = "")
  m <- fit_maihda(y ~ xx + (1 | stratum), data = st$data)

  # With sep = "" the unseen combinations ("1", "2") and ("1", "1") label as "12"
  # and "11", both fitted stratum ids in a 12-stratum fit.
  expect_true(all(c("11", "12") %in% maihda_known_strata(m)))
  nd <- data.frame(xx = 0, p = c("1", "1"), q = c("2", "1"), stringsAsFactors = FALSE)
  expect_identical(
    MAIHDA:::maihda_stratum_labels(nd, m$strata_vars, m$strata_sep), c("12", "11"))

  p <- predict_maihda(m, newdata = nd, allow_new_levels = TRUE)
  fixed_only <- unname(stats::predict(m$model, newdata = data.frame(xx = 0), re.form = NA))
  expect_equal(unname(p), rep(fixed_only, 2), tolerance = 1e-12)
  re <- lme4::ranef(m$model)$stratum
  expect_gt(max(abs(re[c("11", "12"), 1])), 0.1)
})

# ---- the supplied-'stratum' sibling -------------------------------------------

test_that("a supplied stratum that is also an unseen combination's label is refused", {
  skip_on_cran()
  d <- maihda_f01_data()
  m <- fit_maihda(y ~ x + (1 | a), data = d)
  bad <- data.frame(x = 0, a = "1", stratum = "1", stringsAsFactors = FALSE)

  # Refused with AND without the opt-in: allow_new_levels chooses how an unseen
  # stratum is handled, it does not license lending a fitted stratum's effect.
  expect_error(predict_maihda(m, newdata = bad), "would take that stratum's random effect")
  expect_error(predict_maihda(m, newdata = bad, allow_new_levels = TRUE),
               "would take that stratum's random effect")
  expect_error(predict_maihda(m, newdata = bad), "not present when the model was fit")

  # The first offending row is the one reported, whichever kind it is.
  mixed <- data.frame(x = 0, a = c("B", "1"), stratum = c("1", "1"),
                      stringsAsFactors = FALSE)
  first <- expect_error(predict_maihda(m, newdata = mixed))
  expect_match(conditionMessage(first), "row 1", fixed = TRUE)
  expect_false(grepl("would take that stratum's random effect",
                     conditionMessage(first), fixed = TRUE))
  expect_error(predict_maihda(m, newdata = mixed[2:1, , drop = FALSE]),
               "would take that stratum's random effect")

  # Controls, all unchanged: a consistent supplied id, a contradicting id for a SEEN
  # combination (the original message), and a supplied label that is not an id.
  expect_silent(predict_maihda(
    m, newdata = data.frame(x = 0, a = "A", stratum = "1", stringsAsFactors = FALSE)))
  expect_error(predict_maihda(
    m, newdata = data.frame(x = 0, a = "B", stratum = "1", stringsAsFactors = FALSE)),
    "identify stratum '2'")
  round_trip <- predict_maihda(
    m, newdata = data.frame(x = 0, a = "new", stratum = "new", stringsAsFactors = FALSE),
    allow_new_levels = TRUE)
  expect_equal(unname(round_trip),
               unname(stats::predict(m$model, newdata = data.frame(x = 0), re.form = NA)),
               tolerance = 1e-12)
})

# ---- a stratum the table holds but the fit never used -------------------------

test_that("a stratum with no analytic rows is refused like any unseen stratum", {
  skip_on_cran()
  d <- maihda_f01_data()
  d$x[d$a == "T"] <- NA          # every row of stratum "T" leaves the analytic sample
  m <- suppressWarnings(fit_maihda(y ~ x + (1 | a), data = d))

  # The strata table still holds it, with zero analytic rows, and the fit does not.
  expect_true("T" %in% m$strata_info$label)
  expect_identical(m$strata_info$n[m$strata_info$label == "T"], 0L)
  expect_false(as.character(m$strata_info$stratum[m$strata_info$label == "T"]) %in%
                 maihda_known_strata(m))

  nd <- data.frame(x = 0, a = "T", stringsAsFactors = FALSE)
  expect_error(predict_maihda(m, newdata = nd),
               "not present when the model was fit")
  p <- predict_maihda(m, newdata = nd, allow_new_levels = TRUE)
  expect_equal(unname(p),
               unname(stats::predict(m$model, newdata = data.frame(x = 0), re.form = NA)),
               tolerance = 1e-12)

  # The ordinal engine silently returned the zero-effect value here by default,
  # since its linpred helper maps an unknown stratum to 0.
  skip_if_not_installed("ordinal")
  mo <- suppressWarnings(suppressMessages(
    fit_maihda(yo ~ x + (1 | a), data = d, engine = "ordinal")))
  expect_error(predict_maihda(mo, newdata = nd, scale = "link"),
               "not present when the model was fit")
  expect_equal(
    unname(predict_maihda(mo, newdata = nd, scale = "link", allow_new_levels = TRUE)),
    unname(MAIHDA:::maihda_clmm_linpred(mo, newdata = nd, include_re = FALSE)),
    tolerance = 1e-12)
})

# ---- the promise is lme4's: zero the unseen effect, keep the others ------------

test_that("a contextual fit keeps the seen school effect for a colliding stratum", {
  skip_on_cran()
  set.seed(11)
  n <- 1200
  vals <- c("1", "2", "3", "4", "5", "7")        # "6" is never observed
  d <- data.frame(a = sample(vals, n, TRUE),
                  school = sample(paste0("s", 1:10), n, TRUE),
                  x = stats::rnorm(n), stringsAsFactors = FALSE)
  us <- stats::rnorm(10, sd = 0.8)[as.integer(factor(d$school))]
  ua <- stats::rnorm(6, sd = 1.2)[as.integer(factor(d$a))]
  d$y <- 1 + 0.4 * d$x + ua + us + stats::rnorm(n, sd = 0.7)
  m <- suppressMessages(fit_maihda(y ~ x + (1 | a), data = d, context = "school"))

  # Six strata, so the id "6" exists while the dimension VALUE "6" was never seen.
  expect_true("6" %in% maihda_known_strata(m))
  expect_false("6" %in% d$a)

  nd <- data.frame(x = 0, a = "6", school = "s1", stringsAsFactors = FALSE)
  p <- predict_maihda(m, newdata = nd, allow_new_levels = TRUE)
  fixed_only <- unname(stats::predict(m$model, newdata = data.frame(x = 0),
                                      re.form = NA))
  re_school <- lme4::ranef(m$model)$school["s1", 1]
  re_stratum6 <- lme4::ranef(m$model)$stratum["6", 1]

  # Only the unseen stratum effect is dropped; the school the row DID appear in
  # keeps its effect (lme4's allow.new.levels contract).
  expect_equal(unname(p), unname(fixed_only + re_school), tolerance = 1e-12)
  expect_gt(abs(re_school), 0.5)                       # the kept effect matters
  expect_gt(abs(re_stratum6), 0.5)                     # so does the borrowed one
  expect_gt(abs(unname(p) - unname(fixed_only + re_school + re_stratum6)), 0.5)
})

test_that("a crossed-dimensions fit gives an unseen intersection the additive value", {
  skip_on_cran()
  set.seed(31)
  d <- expand.grid(g = c("1", "2", "3"), r = c("1", "2", "3", "4", "5"), rep = 1:20,
                   stringsAsFactors = FALSE)
  d <- d[!(d$g == "1" & d$r == "2"), ]          # this intersection is never observed
  d$x <- stats::rnorm(nrow(d))
  ug <- stats::rnorm(3, sd = 0.9)[as.integer(factor(d$g))]
  ur <- stats::rnorm(5, sd = 0.7)[as.integer(factor(d$r))]
  ui <- stats::rnorm(14, sd = 1.1)[as.integer(factor(paste(d$g, d$r)))]
  d$y <- 1 + 0.3 * d$x + ug + ur + ui + stats::rnorm(nrow(d), sd = 0.5)
  # sep = "" makes the unseen intersection ("1", "2") label as "12", a fitted id.
  st <- make_strata(d, vars = c("g", "r"), sep = "")
  an <- suppressMessages(suppressWarnings(
    maihda(y ~ x + (1 | stratum), data = st$data, decomposition = "crossed-dimensions")))
  m <- if (!is.null(an$adjusted_model)) an$adjusted_model else an$model
  expect_s3_class(m, "maihda_model")

  nd <- data.frame(x = 0, g = "1", r = "2", stringsAsFactors = FALSE)
  lab <- MAIHDA:::maihda_stratum_labels(nd, m$strata_vars, m$strata_sep)
  expect_identical(lab, "12")
  expect_true(lab %in% maihda_known_strata(m))

  # The estimand for an unseen intersection here is the ADDITIVE prediction: fixed
  # effects plus both dimension random effects, with no interaction effect. It is not
  # the fixed-effects-only value, and it is not the additive value plus some other
  # intersection's interaction effect, which is what the collision produced.
  p <- predict_maihda(m, newdata = nd, allow_new_levels = TRUE)
  additive <- unname(stats::predict(m$model, newdata = nd,
                                    re.form = ~ (1 | g) + (1 | r),
                                    allow.new.levels = TRUE))
  fixed_only <- unname(stats::predict(m$model, newdata = nd, re.form = NA))
  re_borrowed <- lme4::ranef(m$model)$stratum[lab, 1]
  expect_equal(unname(p), additive, tolerance = 1e-12)
  expect_gt(abs(additive - fixed_only), 0.5)          # the dimension effects are kept
  expect_gt(abs(re_borrowed), 0.2)                    # the interaction effect matters
  expect_gt(abs(unname(p) - (additive + re_borrowed)), 0.2)
})

test_that("a longitudinal fit drops the unseen stratum's whole trajectory", {
  skip_on_cran()
  set.seed(32)
  np <- 160
  dl <- expand.grid(pid = 1:np, t = 0:3)
  vals <- c("1", "2", "3", "4", "5", "7", "8", "9")   # 8 strata -> ids 1..8; "6" absent
  dl$a <- vals[(dl$pid %% 8) + 1]
  us <- stats::rnorm(8, sd = 0.9)[as.integer(factor(dl$a))]
  sl <- stats::rnorm(8, sd = 0.25)[as.integer(factor(dl$a))]
  up <- stats::rnorm(np, sd = 0.5)[dl$pid]
  dl$y <- 1 + 0.4 * dl$t + us + sl * dl$t + up + stats::rnorm(nrow(dl), sd = 0.4)
  dl$pid <- factor(dl$pid)
  m <- suppressMessages(suppressWarnings(
    fit_maihda(y ~ (1 | a), data = dl, id = "pid", time = "t")))

  expect_true("6" %in% maihda_known_strata(m))        # a fitted id ...
  expect_false("6" %in% dl$a)                         # ... and never a dimension value

  nd <- data.frame(t = c(0, 3), a = "6",
                   pid = factor(c("1", "1"), levels = levels(dl$pid)),
                   stringsAsFactors = FALSE)
  p <- predict_maihda(m, newdata = nd, allow_new_levels = TRUE)
  prep <- MAIHDA:::maihda_prepare_prediction_data(m, nd, allow_new_levels = TRUE)

  # Oracle: lme4 on the same prepared rows with a stratum value that cannot be a
  # fitted level, so the person effect survives and the stratum trajectory does not.
  orc <- prep
  orc$stratum <- "zz_never_a_level"
  oracle <- unname(stats::predict(m$model, newdata = orc, allow.new.levels = TRUE))
  expect_equal(unname(p), oracle, tolerance = 1e-12)

  # The collision borrowed stratum 6's intercept AND slope, so its error grew with
  # time -- the tell that a whole trajectory, not just a level, was lent.
  bad <- prep
  bad$stratum <- "6"
  borrowed <- unname(stats::predict(m$model, newdata = bad, allow.new.levels = TRUE))
  gap <- abs(borrowed - oracle)
  expect_gt(gap[1], 0.02)
  expect_gt(gap[2], gap[1] * 2)
})

# ---- the other engines --------------------------------------------------------

test_that("wemix and ordinal give the zero-effect value for a colliding label", {
  skip_on_cran()
  skip_if_not_installed("WeMix")
  skip_if_not_installed("ordinal")
  d <- maihda_f01_data()
  nd <- data.frame(x = 0, a = c("1", "7"), stringsAsFactors = FALSE)

  mw <- suppressMessages(fit_maihda(y ~ x + (1 | a), data = d, engine = "wemix",
                                    sampling_weights = "w"))
  pw <- predict_maihda(mw, newdata = nd, scale = "link", allow_new_levels = TRUE)
  expect_equal(unname(pw),
               unname(MAIHDA:::maihda_wemix_linpred(mw, newdata = nd, include_re = FALSE)),
               tolerance = 1e-12)
  rw <- MAIHDA:::maihda_wemix_ranef_vector(mw)
  expect_gt(min(abs(rw[c("1", "7")])), 0.2)     # the effects it used to borrow

  mo <- suppressWarnings(suppressMessages(
    fit_maihda(yo ~ x + (1 | a), data = d, engine = "ordinal")))
  po <- predict_maihda(mo, newdata = nd, scale = "link", allow_new_levels = TRUE)
  expect_equal(unname(po),
               unname(MAIHDA:::maihda_clmm_linpred(mo, newdata = nd, include_re = FALSE)),
               tolerance = 1e-12)
  ro <- MAIHDA:::maihda_clmm_stratum_ranef(mo)
  expect_gt(min(abs(ro$random_effect[ro$stratum %in% c("1", "7")])), 0.2)
  # The response scale follows: the expected category score at a zero effect.
  ps <- predict_maihda(mo, newdata = nd, scale = "response", allow_new_levels = TRUE)
  expect_equal(unname(ps),
               unname(MAIHDA:::maihda_ordinal_eta_to_score(
                 MAIHDA:::maihda_clmm_linpred(mo, newdata = nd, include_re = FALSE),
                 MAIHDA:::maihda_clmm_cutpoints(mo$model), mo$family$link)),
               tolerance = 1e-12)
})

# The brms path zeroes an unseen level by dropping its grouping term from the
# re_formula, so what it hands brms is the whole story -- captured without Stan.
maihda_f01_brms_scope <- function(nd) {
  obj <- list(
    formula = y ~ x + (1 | stratum),
    data = data.frame(y = c(0, 1, 0, 1), x = c(0, 1, 0, 1),
                      stratum = c("1", "2", "1", "2"),
                      a = c("A", "B", "A", "B"), stringsAsFactors = FALSE),
    strata_info = data.frame(stratum = c("1", "2"), label = c("A", "B"),
                             a = c("A", "B"), stringsAsFactors = FALSE),
    strata_vars = "a", strata_sep = " x ")
  prep <- MAIHDA:::maihda_prepare_prediction_data(obj, nd, allow_new_levels = TRUE)
  seen <- NULL
  testthat::local_mocked_bindings(
    maihda_brms_predict_rows = function(object, nd, scale, dots) {
      seen <<- c(seen, list(list(rows = nrow(nd), dots = dots)))
      rep(0, nrow(nd))
    })
  MAIHDA:::maihda_brms_individual_prediction(obj, prep, "link", TRUE,
                                             list(allow_new_levels = TRUE))
  list(prep = prep, seen = seen)
}

test_that("brms drops the stratum term for a combination labelled like an id", {
  # Row 1's dimension value "1" is also the fitted stratum id "1": brms used to keep
  # the stratum term and predict it WITH that stratum's random effect.
  got <- maihda_f01_brms_scope(data.frame(x = 0, a = "1", stringsAsFactors = FALSE))
  expect_false(got$prep$stratum %in% c("1", "2"))
  expect_length(got$seen, 1L)
  re_form <- got$seen[[1]]$dots$re_formula
  expect_true(is.logical(re_form) && is.na(re_form))   # fixed effects only

  # A seen combination is untouched: brms keeps every term (no re_formula injected).
  seen_row <- maihda_f01_brms_scope(data.frame(x = 0, a = "A", stringsAsFactors = FALSE))
  expect_identical(seen_row$prep$stratum, "1")
  expect_false("re_formula" %in% names(seen_row$seen[[1]]$dots))
})

test_that("brms predicts a colliding label at a zero stratum effect (Stan)", {
  skip_on_cran()
  skip_if(Sys.getenv("MAIHDA_TEST_BRMS") != "true",
          "brms Stan tests are opt-in; set MAIHDA_TEST_BRMS=true to run them")
  skip_if_not_installed("brms")
  d <- maihda_f01_data()
  m <- suppressWarnings(suppressMessages(fit_maihda(
    y ~ x + (1 | a), data = d, engine = "brms",
    chains = 1, iter = 800, refresh = 0, seed = 1)))

  nd <- data.frame(x = 0, a = c("1", "new", "A"), stringsAsFactors = FALSE)
  p <- predict_maihda(m, newdata = nd, scale = "link", allow_new_levels = TRUE)
  fixed_only <- MAIHDA:::maihda_brms_linpred_mean(
    m$model, newdata = data.frame(x = 0, stratum = "1"), re_formula = NA)

  expect_equal(unname(p[1]), unname(fixed_only), tolerance = 1e-8)
  expect_equal(unname(p[2]), unname(fixed_only), tolerance = 1e-8)
  # The seen row keeps stratum A's effect, which is what row 1 used to borrow.
  re <- brms::ranef(m$model)$stratum[, "Estimate", "Intercept"]
  expect_equal(unname(p[3]), unname(fixed_only + re["1"]), tolerance = 1e-8)
  expect_gt(abs(re["1"]), 0.3)
})
