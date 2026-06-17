# data-raw/precompute_youth_vignette.R
# Precompute everything the ESS youth-unemployment case-study vignette needs, and
# cache *lightweight* summaries to vignettes/youth_precomputed.rds plus static
# figures to vignettes/figures/. The vignette ships no microdata (ESS is license-
# restricted), so -- unlike the other vignettes -- it cannot fit anything live;
# it renders entirely from this cache.
#
# Prereq: data-raw/fetch_ess_youth.R has been run (creates the gitignored frames
# data-raw/ess/ess_youth_analysis.rds and ess_youth_neet.rds).
#
# Re-run from the package root:  Rscript data-raw/precompute_youth_vignette.R
# (brms is slow; set SKIP_BRMS=1 to refresh only the fast lme4/figure parts.)

suppressMessages(devtools::load_all(quiet = TRUE))
suppressMessages({library(dplyr); library(ggplot2)})
RUN_BRMS <- !identical(Sys.getenv("SKIP_BRMS"), "1")

figdir <- "vignettes/figures"; dir.create(figdir, showWarnings = FALSE, recursive = TRUE)
du <- readRDS("data-raw/ess/ess_youth_analysis.rds")
dn <- readRDS("data-raw/ess/ess_youth_neet.rds")
sv <- function(p, f, w = 7, h = 4.5) { ggsave(file.path(figdir, f), p, width = w, height = h, dpi = 150); cat("fig:", f, "\n") }

## ===================== UNEMPLOYMENT (cross-classified) =====================
set.seed(2026)
au <- maihda(unemployed ~ gender + migration + education + (1 | gender:migration:education),
             data = du, context = "country", family = "binomial",
             response_vpc = TRUE, seed = 2026, interactions = FALSE)
ctx <- au$summary$context
da  <- au$summary$discriminatory_accuracy

# additivity: additive vs saturated fixed-effects fit, and obs vs predicted by stratum
m_add <- glm(unemployed ~ gender + migration + education, binomial, du)
m_sat <- glm(unemployed ~ gender * migration * education, binomial, du)
lrt   <- anova(m_add, m_sat, test = "LRT")
du$padd <- predict(m_add, type = "response")
strata <- du |>
  group_by(gender, migration, education) |>
  summarise(n = n(), obs = mean(unemployed), pred = mean(padd), .groups = "drop") |>
  mutate(label = paste(gender, migration, education, sep = "."))

unemp <- list(
  varcomp = data.frame(
    component = c("Between-stratum (intersectional)", "Between-country (contextual)", "Within (residual)"),
    variance = c(ctx$var_stratum, ctx$context_var_total, ctx$within_var),
    vpc      = c(ctx$vpc_stratum, ctx$vpc_context_total, NA)),
  vpc_stratum = ctx$vpc_stratum, vpc_country = ctx$vpc_context_total,
  vpc_response = tryCatch(au$summary$vpc_response$vpc, error = function(e) NA_real_),
  auc = da$auc, mor = da$mor,
  cases = sum(du$unemployed), controls = sum(du$unemployed == 0),
  pcv = list(value = au$pcv$pvc, var_null = au$pcv$var_model1,
             var_adj = au$pcv$var_model2,
             singular = au$pcv$var_model2 < 1e-6),
  fe = list(chisq = lrt$Deviance[2], df = lrt$Df[2], p = lrt$`Pr(>Chi)`[2],
            aic_add = AIC(m_add), aic_sat = AIC(m_sat)),
  maxresid = max(abs(strata$obs - strata$pred)),
  strata = strata)

# figures
sv(plot(au, type = "context_vpc") + theme(plot.subtitle = element_text(size = 9)),
   "ess_partition.png", 9.5, 5.5)
sv(plot(au, type = "effect_decomp") +
     labs(subtitle = "Every stratum deviation is the grey additive (fixed-effect) component; the orange interaction (random) component is ~0") +
     theme(plot.subtitle = element_text(size = 10)),
   "ess_effect_decomp.png", 10, 6)
sv(ggplot(strata, aes(pred, obs)) +
     geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
     geom_point(aes(size = n, colour = education), alpha = 0.85) +
     scale_x_continuous(labels = scales::percent, limits = c(0, 0.26)) +
     scale_y_continuous(labels = scales::percent, limits = c(0, 0.26)) +
     labs(title = "Youth unemployment is additive across intersections",
          subtitle = sprintf("18 strata; additive model reproduces every observed rate to within %.1f pp", 100 * unemp$maxresid),
          x = "Additive-model predicted rate", y = "Observed rate", size = "n", colour = "Education") +
     theme_minimal(base_size = 12), "ess_additivity.png", 8, 5)
ctry <- du |> group_by(country) |> summarise(n = n(), rate = mean(unemployed), .groups = "drop")
sv(ggplot(ctry, aes(reorder(country, rate), rate)) +
     geom_hline(yintercept = mean(du$unemployed), linetype = 2, colour = "grey50") +
     geom_point(aes(size = n), colour = "#2c7fb8") + coord_flip() +
     scale_y_continuous(labels = scales::percent) +
     labs(title = "Country context dwarfs intersection",
          subtitle = sprintf("By country (dashed = pooled %.1f%%). Country VPC %.1f%% > stratum VPC %.1f%%",
                             100 * mean(du$unemployed), 100 * unemp$vpc_country, 100 * unemp$vpc_stratum),
          x = NULL, y = "Youth unemployment rate", size = "n") +
     theme_minimal(base_size = 11), "ess_country.png", 8, 6.5)

## ===================== NEET (the interaction probe) =====================
na2 <- glm(neet ~ gender + migration + education, binomial, dn)
n2w <- glm(neet ~ (gender + migration + education)^2, binomial, dn)
n3w <- glm(neet ~ gender * migration * education, binomial, dn)
lrt2 <- anova(na2, n2w, test = "LRT"); lrt3 <- anova(na2, n3w, test = "LRT")
an <- maihda(neet ~ gender + migration + education + (1 | gender:migration:education),
             data = dn, context = "country", family = "binomial", interactions = FALSE)
cell <- dn |> group_by(gender, migration) |> summarise(rate = mean(neet), n = n(), .groups = "drop")

neet <- list(
  rate = mean(dn$neet), cell = cell,
  fe2 = list(chisq = lrt2$Deviance[2], df = lrt2$Df[2], p = lrt2$`Pr(>Chi)`[2]),
  fe3 = list(chisq = lrt3$Deviance[2], df = lrt3$Df[2], p = lrt3$`Pr(>Chi)`[2]),
  lme4 = list(vpc_stratum = an$summary$context$vpc_stratum, pcv = an$pcv$pvc,
              var_null = an$pcv$var_model1, var_adj = an$pcv$var_model2,
              singular = an$pcv$var_model2 < 1e-6))

precomp <- list(
  meta = list(rounds = "9-11 (2018-2022)", n_unemp = nrow(du), n_neet = nrow(dn),
              n_countries = n_distinct(du$country), n_strata = nrow(strata),
              unemp_rate = mean(du$unemployed), neet_rate = mean(dn$neet),
              min_cell_unemp = min(strata$n)),
  unemp = unemp, neet = neet)
saveRDS(precomp, "vignettes/youth_precomputed.rds", version = 2)
cat(sprintf("\nFAST done. unemp VPC stratum=%.3f country=%.3f AUC=%.3f PCV=%.2f | NEET 2way p=%.3f\n",
            unemp$vpc_stratum, unemp$vpc_country, unemp$auc, unemp$pcv$value, neet$fe2$p))

## ===================== brms (slow; the singular-vs-interval lesson) =====================
if (RUN_BRMS) {
  suppressMessages({library(brms); library(posterior)}); options(mc.cores = 4)
  dn <- dn |> mutate(stratum = interaction(gender, migration, education, drop = TRUE, sep = "."))

  cat("\n[brms fit1] explicit gender x migration interaction ...\n")
  f1 <- brm(neet ~ gender * migration + education + (1 | country), family = bernoulli(),
            data = dn, prior = prior(normal(0, 1), class = "b"),
            chains = 4, iter = 4000, warmup = 1500, seed = 2026, refresh = 0, silent = 2)
  dr <- as_draws_df(f1)
  ints <- lapply(c("b_genderfemale:migrationsecond_gen", "b_genderfemale:migrationfirst_gen"),
                 function(v) { s <- dr[[v]]; data.frame(
                   term = sub("^b_genderfemale:migration", "female x ", v),
                   or_med = exp(median(s)), or_lo = exp(quantile(s, .025)),
                   or_hi = exp(quantile(s, .975)), p_gt0 = mean(s > 0)) })
  interaction <- do.call(rbind, ints); rownames(interaction) <- NULL

  cat("[brms fit2] MAIHDA adjusted (stratum + country REs), mild sd prior ...\n")
  f2 <- brm(neet ~ gender + migration + education + (1 | stratum) + (1 | country),
            family = bernoulli(), data = dn,
            prior = c(prior(normal(0, 1), class = "b"), prior(normal(0, 0.5), class = "sd")),
            chains = 4, iter = 4000, warmup = 1500, seed = 2026,
            control = list(adapt_delta = 0.95), refresh = 0, silent = 2)
  sds <- posterior_summary(f2, variable = "sd_stratum__Intercept")
  np <- nuts_params(f2); ndraws <- length(unique(np$Iteration)) * length(unique(np$Chain))
  diag <- list(max_rhat = max(rhat(f2), na.rm = TRUE),
               min_ess = min(neff_ratio(f2), na.rm = TRUE) * ndraws,
               divergences = sum(np$Parameter == "divergent__" & np$Value == 1))

  precomp$neet$brms <- list(
    interaction = interaction,
    sd_stratum = list(med = sds[1, "Estimate"], lo = sds[1, "Q2.5"], hi = sds[1, "Q97.5"]),
    sd_lme4 = 0,
    diag = diag, prior = "normal(0,0.5) on class='sd'; normal(0,1) on class='b'",
    settings = "4 chains x 4000 iter (1500 warmup), adapt_delta = 0.95")
  saveRDS(precomp, "vignettes/youth_precomputed.rds", version = 2)

  # the climax figure: singular ML point vs calibrated brms interval for the
  # interaction (adjusted between-stratum) standard deviation
  figdat <- data.frame(
    method = factor(c("lme4 (ML)", "brms (Bayesian)"), levels = c("lme4 (ML)", "brms (Bayesian)")),
    sd = c(0, precomp$neet$brms$sd_stratum$med),
    lo = c(NA, precomp$neet$brms$sd_stratum$lo),
    hi = c(NA, precomp$neet$brms$sd_stratum$hi))
  sv(ggplot(figdat, aes(method, sd, colour = method)) +
       geom_pointrange(aes(ymin = lo, ymax = hi), na.rm = TRUE, linewidth = 0.9) +
       geom_point(size = 3) +
       labs(title = "NEET interaction: a singular ML point vs. a calibrated brms interval",
            subtitle = "Adjusted-model between-stratum SD: lme4 collapses to exactly 0; brms gives an interval",
            x = NULL, y = "Between-stratum SD (logit scale)") +
       theme_minimal(base_size = 12) + theme(legend.position = "none"), "ess_neet_brms.png", 9, 4.5)

  cat(sprintf("\nbrms done. migrant-women OR=%.2f [%.2f, %.2f] P(>0)=%.3f | stratum SD=%.3f [%.3f, %.3f] rhat=%.3f div=%d\n",
              interaction$or_med[2], interaction$or_lo[2], interaction$or_hi[2], interaction$p_gt0[2],
              precomp$neet$brms$sd_stratum$med, precomp$neet$brms$sd_stratum$lo,
              precomp$neet$brms$sd_stratum$hi, diag$max_rhat, diag$divergences))
}
cat("saved vignettes/youth_precomputed.rds\n")
