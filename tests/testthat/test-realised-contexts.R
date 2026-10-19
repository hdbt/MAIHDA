# The realised partition of a brms fit with a context: finite-population summaries
# (Gelman 2005, eq. 8 with only the mean removed; Section 3.5) of the effects of the
# strata and context levels in the data. In each posterior draw, a grouping factor
# with C levels in the data and effects u_k has the realised variance
#   sum_k (u_k - mean(u))^2 / (C - 1),
# the strata's (the stratum intercepts, or the zeta_j of a crossed-dimensions fit) is
# var_stratum, each context's is per_context and their sum context_var_total, and any
# other random effect's is other_var. The shares divide by the sum of these and the
# within-stratum variance: vpc_stratum, vpc_context and vpc_context_total. Each is
# summarised by the posterior median and the equal-tailed interval. The shares are
# given where the within-stratum variance stands apart from the variance components:
# a Gaussian identity fit and a binomial, Bernoulli or cumulative fit on its latent
# scale. The oracle numbers below come from a separate implementation of these
# definitions that does not call the package.

# Evaluate `expr` under a fixed seed and put the caller's RNG state back, so this file
# leaves the random stream the later test files read as it found it (see setup-seed.R).
rcx_with_seed <- function(seed, expr) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    if (had) {
      assign(".Random.seed", old, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  }, add = TRUE)
  set.seed(seed)
  force(expr)
}

# Printed lines without colour codes.
rcx_print <- function(expr) cli::ansi_strip(utils::capture.output(expr))

rcx_q <- function(x, conf_level = 0.95) {
  a <- 1 - conf_level
  c(stats::median(x), stats::quantile(x, c(a / 2, 1 - a / 2), names = FALSE))
}

# Five draws built without random numbers: 6 strata, a 3-level and a 2-level context,
# a 4-level further random effect and a within-stratum variance per draw.
rcx_draws <- function() {
  nd <- 5L
  list(V = matrix(sin(seq_len(30) * 1.7) / 2, nd, 6,
                  dimnames = list(NULL, paste0("s", 1:6))),
       A = matrix(cos(seq_len(15) * 0.9) / 3, nd, 3,
                  dimnames = list(NULL, c("a1", "a2", "a3"))),
       B = matrix(sin(seq_len(10) * 2.1) / 4, nd, 2, dimnames = list(NULL, c("b1", "b2"))),
       G = matrix(cos(seq_len(20) * 1.3) / 5, nd, 4,
                  dimnames = list(NULL, paste0("g", 1:4))),
       w = c(1.1, 0.9, 1.3, 1.0, 1.2))
}

# A contextual maihda_model shell around a brmsfit stub that carries its own posterior
# draws (an `r_` column marks that level effects were kept) and level-effect draws for
# the mocked reader. Every stratum sits in every level of `site`; `region` and `school`
# are spread over the rows.
rcx_stub <- function(family = "gaussian", link = "identity",
                     contexts = c("site", "region"), other = FALSE, level_draws = TRUE) {
  dr <- rcx_draws()
  dat <- expand.grid(stratum = paste0("s", 1:6), site = c("a1", "a2", "a3"),
                     stringsAsFactors = FALSE)
  dat$region <- rep(c("b1", "b2"), length.out = nrow(dat))
  dat$school <- rep(paste0("g", 1:4), length.out = nrow(dat))
  draws <- data.frame(b_Intercept = rep(0, 5), sd_stratum__Intercept = seq(0.3, 0.7, 0.1),
                      sd_site__Intercept = seq(0.5, 0.9, 0.1),
                      sd_region__Intercept = seq(0.2, 0.6, 0.1), sigma = sqrt(dr$w))
  if (!"region" %in% contexts) {
    draws$sd_region__Intercept <- NULL
  }
  if (other) {
    draws$sd_school__Intercept <- seq(0.1, 0.5, 0.1)
  }
  if (level_draws) {
    draws[["r_stratum[s1,Intercept]"]] <- dr$V[, 1]
  }
  fam <- list(family = family, link = link)
  structure(
    list(model = structure(list(formula = y ~ 1, family = fam, draws = draws,
                                levels = list(stratum = dr$V, site = dr$A, region = dr$B,
                                              school = dr$G)),
                           class = "brmsfit"),
         engine = "brms", formula = y ~ 1, family = fam, data = dat,
         context_info = list(context_vars = contexts)),
    class = "maihda_model")
}

# A crossed-dimensions maihda_model shell: a 2 x 3 grid of strata, each in both levels
# of `site`.
rcx_cc_stub <- function(family = "bernoulli", link = "logit", site_draws = TRUE) {
  cells <- expand.grid(ga = c("p", "q"), gb = c("x", "y", "z"), stringsAsFactors = FALSE)
  nd <- 5L
  Ua <- matrix(sin(1:10), nd, 2, dimnames = list(NULL, c("p", "q")))
  Ub <- matrix(cos(1:15) / 2, nd, 3, dimnames = list(NULL, c("x", "y", "z")))
  Vs <- matrix(sin(seq_len(30) * 2.3) / 3, nd, 6, dimnames = list(NULL, paste0("s", 1:6)))
  S <- matrix(c(0.3, -0.2, 0.1, 0.4, -0.1, -0.3, 0.2, 0, -0.4, 0.1), nd, 2,
              dimnames = list(NULL, c("u", "v")))
  dat <- cells
  dat$stratum <- paste0("s", 1:6)
  dat <- rbind(dat, dat)
  dat$site <- rep(c("u", "v"), each = 6)
  draws <- data.frame(b_Intercept = rep(0, 5), sd_stratum__Intercept = seq(0.1, 0.5, 0.1),
                      sd_ga__Intercept = seq(0.2, 1, 0.2),
                      sd_gb__Intercept = seq(0.3, 0.7, 0.1),
                      sd_site__Intercept = seq(0.4, 0.8, 0.1), sigma = seq(0.9, 1.3, 0.1))
  draws[["r_stratum[s1,Intercept]"]] <- Vs[, 1]
  levels <- list(ga = Ua, gb = Ub, stratum = Vs)
  if (site_draws) {
    levels$site <- S
  }
  fam <- list(family = family, link = link)
  z <- Vs + Ua[, match(cells$ga, colnames(Ua))] + Ub[, match(cells$gb, colnames(Ub))]
  structure(
    list(model = structure(list(formula = y ~ 1, family = fam, draws = draws,
                                levels = levels), class = "brmsfit"),
         engine = "brms", formula = y ~ 1, family = fam, data = dat,
         cc_info = list(dim_groups = c(sex = "ga", region = "gb"),
                        interaction_group = "stratum"),
         context_info = list(context_vars = "site"),
         # for the tests: the strata's effects and the context's level effects
         truth = list(z = z, site = S)),
    class = "maihda_model")
}

rcx_mock_brms <- function(env = parent.frame()) {
  local_mocked_bindings(
    maihda_posterior_draws_brms = function(model) model$draws,
    maihda_level_draws_brms = function(model, groups) model$levels[groups],
    .env = env
  )
}

test_that("the realised partition follows its definitions", {
  dr <- rcx_draws()
  sv <- apply(dr$V, 1, stats::var)
  sa <- apply(dr$A, 1, stats::var)
  sb <- apply(dr$B, 1, stats::var)
  sg <- apply(dr$G, 1, stats::var)
  r <- maihda_context_realised_from_draws(sv, list(site = sa, region = sb), other = sg,
                                          within = dr$w, n_strata = 6L)

  # The separate implementation's values, at full precision.
  expect_equal(c(r$var_stratum, r$var_stratum_ci),
               c(0.15623583482948078, 0.13170061536591882, 0.16882619005419094),
               tolerance = 1e-12)
  expect_equal(c(r$per_context[["site"]], unname(r$per_context_ci["site", ])),
               c(0.085276180035156732, 0.05519155706689486, 0.10148469672432121),
               tolerance = 1e-12)
  expect_equal(c(r$per_context[["region"]], unname(r$per_context_ci["region", ])),
               c(0.025559527309432145, 0.020339087502429114, 0.092151982641671154),
               tolerance = 1e-12)
  expect_equal(c(r$context_var_total, r$context_var_total_ci),
               c(0.123306261951887, 0.10735799976717769, 0.15429708630083455),
               tolerance = 1e-12)
  expect_equal(c(r$other_var, r$other_var_ci),
               c(0.0014493411711179851, 0.00023356651709177025, 0.0029607589390600632),
               tolerance = 1e-12)
  expect_equal(c(r$vpc_stratum, r$vpc_stratum_ci),
               c(0.11098123221225603, 0.10431078123332715, 0.11612481265709071),
               tolerance = 1e-12)
  expect_equal(c(r$vpc_context[["site"]], unname(r$vpc_context_ci["site", ])),
               c(0.055058380643555262, 0.041823045871093679, 0.080130570922880576),
               tolerance = 1e-12)
  expect_equal(c(r$vpc_context[["region"]], unname(r$vpc_context_ci["region", ])),
               c(0.016180157432112922, 0.015778685539080232, 0.076067026429232576),
               tolerance = 1e-12)
  expect_equal(c(r$vpc_context_total, r$vpc_context_total_ci),
               c(0.098190834390964904, 0.071931622642365842, 0.12147457168045027),
               tolerance = 1e-12)

  # The same, from the definitions written out with base R.
  total <- sv + sa + sb + sg + dr$w
  expect_equal(c(r$var_stratum, r$var_stratum_ci), rcx_q(sv))
  expect_equal(c(r$context_var_total, r$context_var_total_ci), rcx_q(sa + sb))
  expect_equal(c(r$vpc_stratum, r$vpc_stratum_ci), rcx_q(sv / total))
  expect_equal(c(r$vpc_context_total, r$vpc_context_total_ci), rcx_q((sa + sb) / total))
  expect_equal(unname(r$vpc_context_ci["region", ]), rcx_q(sb / total)[2:3])
  # With two levels the realised variance is half the squared difference.
  expect_equal(sb, (dr$B[, 1] - dr$B[, 2])^2 / 2)
  expect_identical(names(r), c(
    "var_stratum", "var_stratum_ci", "per_context", "per_context_ci",
    "context_var_total", "context_var_total_ci", "other_var", "other_var_ci",
    "vpc_stratum", "vpc_stratum_ci", "vpc_context", "vpc_context_ci",
    "vpc_context_total", "vpc_context_total_ci", "n_strata", "conf_level", "method"))
  expect_identical(dimnames(r$per_context_ci), list(c("site", "region"), c("lower", "upper")))
  expect_identical(dimnames(r$vpc_context_ci), dimnames(r$per_context_ci))
  expect_identical(r$n_strata, 6L)
  expect_identical(r$conf_level, 0.95)
  expect_identical(r$method, "posterior")

  # Without another random effect the total is the strata, the contexts and the
  # within-stratum variance; without a within-stratum variance there are no shares.
  r0 <- maihda_context_realised_from_draws(sv, list(site = sa, region = sb), within = dr$w)
  expect_false(any(c("other_var", "other_var_ci") %in% names(r0)))
  expect_equal(c(r0$vpc_stratum, r0$vpc_stratum_ci),
               c(0.11104157141383512, 0.1044594946184087, 0.1161650022546373),
               tolerance = 1e-12)
  expect_equal(c(r0$vpc_context_total, r0$vpc_context_total_ci),
               c(0.098304290770759209, 0.072045553478672045, 0.12149583830893217),
               tolerance = 1e-12)
  rn <- maihda_context_realised_from_draws(sv, list(site = sa))
  expect_false(any(grepl("^vpc", names(rn))))
  expect_identical(rn$per_context_ci["site", ], c(lower = rn$per_context_ci[1, 1],
                                                 upper = rn$per_context_ci[1, 2]))

  # The interval follows conf_level, and the point estimate can be the posterior mean.
  r90 <- maihda_context_realised_from_draws(sv, list(site = sa, region = sb), other = sg,
                                            within = dr$w, conf_level = 0.9)
  expect_equal(c(r90$vpc_stratum, r90$vpc_stratum_ci),
               c(0.11098123221225603, 0.10431384341913071, 0.11592994199247345),
               tolerance = 1e-12)
  expect_identical(r90$conf_level, 0.9)
  rm <- maihda_context_realised_from_draws(sv, list(site = sa, region = sb), other = sg,
                                           within = dr$w, point = "mean")
  expect_equal(rm$vpc_stratum, 0.11006359043251643, tolerance = 1e-12)
  expect_equal(rm$var_stratum, 0.15181643668832137, tolerance = 1e-12)
  expect_identical(rm$vpc_stratum_ci, r$vpc_stratum_ci)

  # A context with one level in the data has no realised variance, and so no shares.
  r1 <- maihda_context_realised_from_draws(sv, list(site = rep(NA_real_, 5)), within = dr$w)
  expect_true(is.na(r1$per_context[["site"]]) && is.na(r1$vpc_stratum))
})

test_that("each grouping factor's levels in the fitted data are read once, by name", {
  dr <- rcx_draws()
  # stats::var() of each draw's level effects, up to rounding: stats::var() rounds its
  # mean and products differently, so the last bit can differ (it does on arm64).
  fp <- maihda_fp_var_draws(dr$A)
  expect_null(attributes(fp))
  expect_length(fp, 5L)
  expect_lt(max(abs(fp / apply(dr$A, 1, stats::var) - 1)), 1e-12)
  expect_identical(maihda_fp_var_draws(dr$A[, 1, drop = FALSE]), rep(NA_real_, 5))
  d <- data.frame(site = c("a2", NA, "a1", "a2", "a3"), stringsAsFactors = FALSE)
  expect_identical(maihda_levels_in_data(d, "site"), c("a2", "a1", "a3"))
  expect_error(maihda_levels_in_data(d, "region"), "do not hold the 'region' column")

  # Levels in another order, and a level the posterior holds that the data do not (as
  # brms keeps with drop_unused_levels = FALSE): only the data's levels count.
  m <- rcx_stub()
  A4 <- cbind(dr$A[, c(3, 1, 2)], a9 = 100)
  m$model$levels$site <- A4
  local_mocked_bindings(maihda_level_draws_brms = function(model, groups) model$levels[groups])
  v <- maihda_realised_group_vars_brms(m, c("site", "stratum"))
  expect_identical(names(v), c("site", "stratum"))
  expect_equal(v$site, apply(dr$A, 1, stats::var))
  expect_equal(v$stratum, apply(dr$V, 1, stats::var))
  # A level of the data the posterior lacks is an error.
  m$model$levels$site <- dr$A[, 1:2]
  expect_error(maihda_realised_group_vars_brms(m, "site"),
               "A level of 'site' in the fitted data has no effect in the brms posterior")
})

test_that("a brms contextual summary carries the realised partition of its draws", {
  skip_if_not_installed("brms")
  rcx_mock_brms()
  dr <- rcx_draws()
  m <- rcx_stub(other = TRUE)
  expect_no_warning(res <- maihda_context_summary_brms(m, m$context_info, conf_level = 0.95))
  ctx <- res$context
  # Placed last; everything else is what a fit without level effects gives.
  m0 <- rcx_stub(other = TRUE, level_draws = FALSE)
  res0 <- maihda_context_summary_brms(m0, m0$context_info, conf_level = 0.95)
  expect_identical(names(ctx), c(names(res0$context), "realised"))
  expect_identical(ctx[names(ctx) != "realised"], res0$context)
  expect_identical(res[c("variance_components", "vpc_result")],
                   res0[c("variance_components", "vpc_result")])

  # The realised partition is the pure helper's on the level effects in the data.
  sigma2 <- m$model$draws$sigma^2
  expect_equal(ctx$realised, maihda_context_realised_from_draws(
    apply(dr$V, 1, stats::var),
    list(site = apply(dr$A, 1, stats::var), region = apply(dr$B, 1, stats::var)),
    other = apply(dr$G, 1, stats::var), within = sigma2, n_strata = 6L))
  expect_identical(ctx$realised$n_strata, 6L)
  # A credible level and the posterior mean reach it.
  r90 <- maihda_context_summary_brms(m, m$context_info, conf_level = 0.9)$context$realised
  expect_identical(r90$conf_level, 0.9)
  expect_equal(r90$vpc_context_total_ci, rcx_q(
    (apply(dr$A, 1, stats::var) + apply(dr$B, 1, stats::var)) /
      (apply(dr$V, 1, stats::var) + apply(dr$A, 1, stats::var) +
         apply(dr$B, 1, stats::var) + apply(dr$G, 1, stats::var) + sigma2), 0.9)[2:3])
  rmean <- maihda_context_summary_brms(m, m$context_info, 0.95, point = "mean")$context
  expect_equal(rmean$realised$var_stratum, mean(apply(dr$V, 1, stats::var)))

  # One context, no other random effect.
  m1 <- rcx_stub(contexts = "site")
  r1 <- maihda_context_summary_brms(m1, m1$context_info, 0.95)$context$realised
  expect_false("other_var" %in% names(r1))
  expect_identical(names(r1$per_context), "site")
  expect_equal(c(r1$vpc_stratum, r1$vpc_stratum_ci), rcx_q(
    apply(dr$V, 1, stats::var) /
      (apply(dr$V, 1, stats::var) + apply(dr$A, 1, stats::var) + sigma2)))
})

test_that("the shares follow the family", {
  skip_if_not_installed("brms")
  rcx_mock_brms()
  dr <- rcx_draws()
  sv <- apply(dr$V, 1, stats::var)
  sa <- apply(dr$A, 1, stats::var)
  give <- function(family, link, within = 1) {
    m <- rcx_stub(family, link, contexts = "site")
    maihda_context_realised_brms(m, "site", character(0), rep(within, 5), 0.95)
  }
  # Latent scales: the within-stratum variance is a constant of the link.
  for (fl in list(c("binomial", "logit"), c("bernoulli", "logit"), c("cumulative", "logit"),
                  c("binomial", "probit"), c("bernoulli", "cloglog"),
                  c("gaussian", "identity"))) {
    r <- give(fl[1], fl[2], within = 2.5)
    expect_equal(r$vpc_stratum, stats::median(sv / (sv + sa + 2.5)))
    expect_equal(r$vpc_context_total, stats::median(sa / (sv + sa + 2.5)))
  }
  # Count families: their within-stratum variance is built from the variance components,
  # so they get the variances only.
  for (fam in c("poisson", "negbinomial")) {
    r <- give(fam, "log")
    expect_false(any(grepl("^vpc", names(r))))
    expect_equal(r$var_stratum, stats::median(sv))
    expect_equal(r$per_context[["site"]], stats::median(sa))
  }
  # Through the summary: a binomial fit's latent within-stratum variance.
  m <- rcx_stub("binomial", "logit", contexts = "site")
  r <- maihda_context_summary_brms(m, m$context_info, 0.95)$context$realised
  expect_equal(c(r$vpc_stratum, r$vpc_stratum_ci), rcx_q(sv / (sv + sa + pi^2 / 3)))
})

test_that("a failure to read the level effects warns and leaves the rest as it was", {
  skip_if_not_installed("brms")
  rcx_mock_brms()
  m <- rcx_stub()
  m0 <- rcx_stub(level_draws = FALSE)
  res0 <- maihda_context_summary_brms(m0, m0$context_info, 0.95)
  local_mocked_bindings(maihda_level_draws_brms = function(model, groups) {
    stop("no level draws here")
  })
  expect_warning(res <- maihda_context_summary_brms(m, m$context_info, 0.95),
                 "realised contextual partition.*no level draws here")
  expect_identical(res, res0)

  # A context level of the data the posterior lacks: the same warning, naming it.
  local_mocked_bindings(maihda_level_draws_brms = function(model, groups) {
    lv <- model$levels[groups]
    lv$region <- lv$region[, 1, drop = FALSE]
    lv
  })
  expect_warning(res <- maihda_context_summary_brms(m, m$context_info, 0.95),
                 "A level of 'region' in the fitted data has no effect")
  expect_null(res$context$realised)
})

test_that("a crossed-dimensions summary with a context: one set of realised strata", {
  skip_if_not_installed("brms")
  rcx_mock_brms()
  m <- rcx_cc_stub()
  z <- m$truth$z
  between <- apply(z, 1, stats::var)
  ctx <- apply(m$truth$site, 1, stats::var)
  expect_no_warning(res <- maihda_cc_summary_brms(m, m$cc_info, conf_level = 0.95))
  d <- res$decomposition$realised
  r <- res$context$realised
  # The realised strata's VPC has the context's realised variance in its denominator,
  # and the realised partition's stratum share is that VPC.
  expect_equal(c(d$vpc, d$vpc_ci), rcx_q(between / (between + ctx + pi^2 / 3)))
  expect_identical(r$vpc_stratum, d$vpc)
  expect_identical(r$vpc_stratum_ci, d$vpc_ci)
  expect_identical(r$var_stratum, d$between_var)
  expect_identical(r$var_stratum_ci, d$between_var_ci)
  expect_equal(c(r$per_context[["site"]], unname(r$per_context_ci["site", ])), rcx_q(ctx))
  expect_equal(c(r$vpc_context_total, r$vpc_context_total_ci),
               rcx_q(ctx / (between + ctx + pi^2 / 3)))
  expect_identical(r$n_strata, 6L)
  expect_false("other_var" %in% names(r))
  # A credible level and the posterior mean reach it.
  r90 <- maihda_cc_summary_brms(m, m$cc_info, conf_level = 0.9)$context$realised
  expect_identical(r90$conf_level, 0.9)
  expect_equal(c(r90$vpc_context_total, r90$vpc_context_total_ci),
               rcx_q(ctx / (between + ctx + pi^2 / 3), 0.9))
  rmean <- maihda_cc_summary_brms(m, m$cc_info, 0.95, point = "mean")$context$realised
  expect_equal(rmean$per_context[["site"]], mean(ctx))
  expect_equal(rmean$vpc_stratum, mean(between / (between + ctx + pi^2 / 3)))
  # The parameter partition is what it was.
  m0 <- m
  m0$model$draws[["r_stratum[s1,Intercept]"]] <- NULL
  res0 <- suppressWarnings(maihda_cc_summary_brms(m0, m0$cc_info, conf_level = 0.95))
  expect_identical(res$context[names(res$context) != "realised"], res0$context)
  expect_identical(res$vpc_result, res0$vpc_result)

  # A count family: variances only, and no VPC in the realised strata.
  mp <- rcx_cc_stub(family = "poisson", link = "log")
  local_mocked_bindings(maihda_residual_variance_draws_brms = function(model, draws, ...) {
    rep(0.4, nrow(draws))
  })
  resp <- maihda_cc_summary_brms(mp, mp$cc_info, conf_level = 0.95)
  expect_null(resp$decomposition$realised$vpc)
  expect_false(any(grepl("^vpc", names(resp$context$realised))))
  expect_equal(resp$context$realised$per_context[["site"]], stats::median(ctx))
})

test_that("the note points to the realised partition only when it is printed", {
  base <- maihda_context_few_levels_note(c(ctx = 2L))
  expect_identical(maihda_context_few_levels_note(c(ctx = 2L), pointer = "none"), base)
  block <- maihda_context_few_levels_note(c(ctx = 2L), pointer = "block")
  expect_identical(paste(trimws(block), collapse = " "), paste(
    "Note: context 'ctx' has 2 levels. With engine = \"brms\" the variance of such a",
    "context depends on the prior on its standard deviation (with two levels roughly in",
    "proportion to that prior's scale), and the context share and the VPC depend on it;",
    "the realised partition below depends little on it when each context has enough",
    "data. See ?maihda."))
  lines <- maihda_context_few_levels_note(c(school = 4L), pointer = "lines")
  expect_identical(paste(trimws(lines), collapse = " "), paste(
    "Note: context 'school' has 4 levels. With engine = \"brms\" the variance of such a",
    "context depends on the prior on its standard deviation, and the context share and",
    "the VPC depend on it; the realised VPC and context share below depend little on it",
    "when each context has enough data. See ?maihda."))
  for (note in list(block, lines)) {
    expect_true(all(startsWith(note, "  ")))
    expect_true(all(nchar(note) <= 77L))
  }
  expect_identical(maihda_context_few_levels_note(integer(0), pointer = "block"),
                   character(0))
  expect_error(maihda_context_few_levels_note(c(ctx = 2L), pointer = "elsewhere"))
})

test_that("print() shows the realised block below the contextual partition", {
  dr <- rcx_draws()
  sv <- apply(dr$V, 1, stats::var)
  r <- maihda_context_realised_from_draws(
    sv, list(site = apply(dr$A, 1, stats::var), region = apply(dr$B, 1, stats::var)),
    other = apply(dr$G, 1, stats::var), within = dr$w, n_strata = 6L)
  f4 <- function(est, ci) sprintf("%.4f [%.4f, %.4f]", est, ci[1], ci[2])
  pc <- function(est, ci) sprintf("%.1f%% [%.1f%%, %.1f%%]", 100 * est, 100 * ci[1],
                                  100 * ci[2])
  block <- c(
    "Realised strata and contexts (6 strata; 3 levels of 'site', 2 levels of 'region'):",
    "  (Posterior 95% credible interval)",
    sprintf("  %-26s %s  share %s", "Between-stratum variance:",
            f4(r$var_stratum, r$var_stratum_ci), pc(r$vpc_stratum, r$vpc_stratum_ci)),
    sprintf("  %-26s %s  share %s", "Context 'site' variance:",
            f4(r$per_context[["site"]], r$per_context_ci["site", ]),
            pc(r$vpc_context[["site"]], r$vpc_context_ci["site", ])),
    sprintf("  %-26s %s  share %s", "Context 'region' variance:",
            f4(r$per_context[["region"]], r$per_context_ci["region", ]),
            pc(r$vpc_context[["region"]], r$vpc_context_ci["region", ])),
    sprintf("  %-26s %s  share %s", "All contexts combined:",
            f4(r$context_var_total, r$context_var_total_ci),
            pc(r$vpc_context_total, r$vpc_context_total_ci)),
    sprintf("  %-26s %s", "Other random effects:", f4(r$other_var, r$other_var_ci)),
    "  Note: finite-population summaries of the effects of these strata and",
    "  contexts (Gelman 2005); they depend little on the prior when each has",
    "  enough data. See ?maihda.",
    "")
  expect_identical(rcx_print(maihda_print_context_realised(r, c(site = 3L, region = 2L))),
                   block)

  # Without shares (a count family) the lines carry the variances only.
  rn <- maihda_context_realised_from_draws(sv, list(site = apply(dr$A, 1, stats::var)),
                                           n_strata = 6L)
  outn <- rcx_print(maihda_print_context_realised(rn, c(site = 3L)))
  expect_identical(outn[3:4], c(
    sprintf("  %-25s %s", "Between-stratum variance:", f4(rn$var_stratum, rn$var_stratum_ci)),
    sprintf("  %-25s %s", "Context 'site' variance:",
            f4(rn$per_context[["site"]], rn$per_context_ci["site", ]))))
  expect_false(any(grepl("share|All contexts|Other random", outn)))

  # Inside the contextual partition: below its notes, the few-levels note pointing to it.
  ctx <- list(context_vars = "site", var_stratum = 0.2, per_context = c(site = 0.6),
              context_var_total = 0.6, within_var = 3.29, other_var = 0,
              vpc_stratum = 0.05, vpc_context = c(site = 0.15), vpc_context_total = 0.15,
              vpc_context_total_ci = c(0.01, 0.7), bootstrap = FALSE, method = "posterior",
              conf_level = 0.95, context_levels = c(site = 3L))
  out0 <- rcx_print(maihda_print_context_partition(ctx))
  ctx$realised <- rn
  out <- rcx_print(maihda_print_context_partition(ctx))
  note0 <- maihda_context_few_levels_note(c(site = 3L))
  note <- maihda_context_few_levels_note(c(site = 3L), pointer = "block")
  i <- match(note[1], out)
  expect_identical(out[seq_len(i - 1L)], out0[seq_len(i - 1L)])
  expect_identical(out[i + seq_along(note) - 1L], note)
  expect_identical(out[i + length(note)], "")
  expect_identical(out[-seq_len(i + length(note))], outn)
  expect_identical(out0[i + seq_along(note0) - 1L], note0)
  expect_length(out0, i + length(note0))
})

test_that("print() of a two-model analysis shows the realised VPC and context share", {
  g <- expand.grid(d1 = c("a", "b"), d2 = c("p", "q"), site = c("s1", "s2", "s3"),
                   stringsAsFactors = FALSE)
  d <- g[rep(seq_len(nrow(g)), each = 15), ]
  rownames(d) <- NULL
  d$y <- 0.6 * (d$d1 == "b") + 0.4 * (d$d2 == "q") +
    c(s1 = -0.3, s2 = 0, s3 = 0.3)[d$site] + sin(seq_len(nrow(d)) * 12.9898)
  a <- rcx_with_seed(1, suppressWarnings(suppressMessages(
    maihda(y ~ 1 + (1 | d1:d2), data = d, context = "site"))))
  # An lme4 partition has no realised element, and its prints have no realised lines.
  expect_null(a$summary$context$realised)
  out_a <- rcx_print(print(a))
  out_s <- rcx_print(print(a$summary))
  expect_false(any(grepl("Realised", c(out_a, out_s), fixed = TRUE)))

  # The same analysis with the partition of a brms fit that carries one.
  dr <- rcx_draws()
  r <- maihda_context_realised_from_draws(apply(dr$V, 1, stats::var),
                                          list(site = apply(dr$A, 1, stats::var)),
                                          within = dr$w, n_strata = 4L)
  b <- a
  b$summary$context$method <- "posterior"
  b$summary$context$context_levels <- c(site = 3L)
  b$summary$context$realised <- r
  out_b <- rcx_print(print(b))
  note <- maihda_context_few_levels_note(c(site = 3L), pointer = "lines")
  i <- grep("^Context share \\(null\\):", out_b)
  expect_length(i, 1L)
  expect_identical(out_b[i + seq_along(note)], note)
  k <- i + length(note)
  expect_identical(out_b[k + 1:2], c(
    sprintf("Realised VPC (null):           %.4f [%.4f, %.4f]", r$vpc_stratum,
            r$vpc_stratum_ci[1], r$vpc_stratum_ci[2]),
    sprintf("Realised context share (null): %.4f [%.4f, %.4f]", r$vpc_context_total,
            r$vpc_context_total_ci[1], r$vpc_context_total_ci[2])))
  expect_identical(out_b[-(i + seq_len(length(note) + 2L))], out_a)

  # A partition without shares (a count family): no lines, and the note has no pointer.
  b$summary$context$realised <- maihda_context_realised_from_draws(
    apply(dr$V, 1, stats::var), list(site = apply(dr$A, 1, stats::var)), n_strata = 4L)
  out_c <- rcx_print(print(b))
  note_plain <- maihda_context_few_levels_note(c(site = 3L))
  expect_identical(out_c[i + seq_along(note_plain)], note_plain)
  expect_identical(out_c[-(i + seq_along(note_plain))], out_a)
  expect_false(any(grepl("^Realised", out_c)))

  # The summary prints the block below the partition; nothing else moves.
  b$summary$context$realised <- r
  out_bs <- rcx_print(print(b$summary))
  block <- rcx_print(maihda_print_context_realised(r, c(site = 3L)))
  j <- match(block[1], out_bs)
  expect_identical(out_bs[j + seq_along(block) - 1L], block)
  # The note's lines end one blank line above the block.
  note_b <- maihda_context_few_levels_note(c(site = 3L), pointer = "block")
  at_note <- j - length(note_b) - 2L + seq_along(note_b)
  expect_identical(out_bs[at_note], note_b)
  expect_identical(out_bs[j - 1L], "")
  expect_identical(out_bs[-c(at_note, j + seq_along(block) - 1L)], out_s)
})

test_that("the help page defines the realised partition and cites the source", {
  skip_on_cran()
  man <- testthat::test_path("..", "..", "man")
  # Present when the suite runs from the package tree; R CMD check copies only tests/.
  skip_if_not(dir.exists(man), "man/ is not available from this test run")
  flat <- function(f) {
    gsub("[[:space:]]+", " ", paste(readLines(file.path(man, f), warn = FALSE),
                                    collapse = " "))
  }
  rd <- flat("maihda.Rd")
  expect_true(grepl("\\section{Realised strata and contexts}{", rd, fixed = TRUE))
  expect_true(grepl("10.1214/009053604000001048", rd, fixed = TRUE))
  expect_true(grepl("(C - 1)", rd, fixed = TRUE))
  expect_true(grepl("depend little on the prior on the standard deviations when each",
                    rd, fixed = TRUE))
  for (field in c("var_stratum", "per_context", "context_var_total", "other_var",
                  "vpc_stratum", "vpc_context", "vpc_context_total")) {
    expect_true(grepl(sprintf("\\item{\\code{%s}}", field), rd, fixed = TRUE) ||
                  grepl(sprintf("\\code{%s}", field), rd, fixed = TRUE))
  }
  expect_true(grepl("the realised partition in \\code{realised}",
                    flat("summary.maihda_model.Rd"), fixed = TRUE))
})

test_that("a real brms fit: the realised partition from the fit's own level effects", {
  # Compiles a Stan model, so it is OPT-IN: set MAIHDA_TEST_BRMS=true to run it.
  skip_on_cran()
  skip_if(Sys.getenv("MAIHDA_TEST_BRMS") != "true",
          "brms Stan tests are opt-in; set MAIHDA_TEST_BRMS=true to run them")
  skip_if_not_installed("brms")

  # 2 x 2 strata crossed with three contexts, 60 binomial trials a cell. brms draws its
  # sampler seed from R's stream, so the fit sits inside the seeded block too.
  m <- rcx_with_seed(20261018, {
    g <- expand.grid(d1 = c("a", "b"), d2 = c("p", "q"), ctx = c("k1", "k2", "k3"),
                     stringsAsFactors = FALSE)
    eta <- -0.6 + 0.5 * (g$d1 == "b") + 0.4 * (g$d2 == "q") +
      c(k1 = -0.3, k2 = 0, k3 = 0.3)[g$ctx] + stats::rnorm(nrow(g), 0, 0.1)
    g$n <- 60L
    g$s <- stats::rbinom(nrow(g), 60L, stats::plogis(eta))
    suppressWarnings(suppressMessages(
      fit_maihda(s | trials(n) ~ 1 + (1 | d1:d2), data = g, context = "ctx",
                 engine = "brms", family = "binomial", chains = 2, iter = 1000,
                 refresh = 0)))
  })
  # summary() passes on brms's own warnings about the sampler (divergent transitions).
  s <- suppressWarnings(summary(m))
  r <- s$context$realised
  expect_false(is.null(r))
  expect_identical(r$n_strata, 4L)

  # The partition assembled here from brms's own level-effect draws.
  re <- brms::ranef(m$model, summary = FALSE)
  strata <- unique(as.character(m$data$stratum))
  sv <- apply(re$stratum[, strata, "Intercept"], 1, stats::var)
  cv <- apply(re$ctx[, c("k1", "k2", "k3"), "Intercept"], 1, stats::var)
  total <- sv + cv + pi^2 / 3
  expect_equal(c(r$var_stratum, r$var_stratum_ci), rcx_q(sv), tolerance = 1e-10)
  expect_equal(c(r$per_context[["ctx"]], unname(r$per_context_ci["ctx", ])), rcx_q(cv),
               tolerance = 1e-10)
  expect_equal(c(r$vpc_stratum, r$vpc_stratum_ci), rcx_q(sv / total), tolerance = 1e-10)
  expect_equal(c(r$vpc_context_total, r$vpc_context_total_ci), rcx_q(cv / total),
               tolerance = 1e-10)

  # print() shows the block below the partition, under a note that points to it.
  out <- rcx_print(suppressWarnings(print(s)))
  note <- maihda_context_few_levels_note(c(ctx = 3L), pointer = "block")
  i <- match(note[1], out)
  expect_identical(out[i + seq_along(note) - 1L], note)
  expect_identical(out[i + length(note) + 1L],
                   "Realised strata and contexts (4 strata; 3 levels of 'ctx'):")
})
