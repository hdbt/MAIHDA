# data-raw/maihda_youth_data.R
# Simulate `maihda_youth_data`: a cross-national, cross-sectional youth dataset for
# demonstrating MAIHDA on the questions youth-research institutes care about --
# the school-to-work transition (NEET status) and youth subjective well-being --
# at the intersection of gender, migration background, and social origin.
#
# The dataset is SYNTHETIC. It is modelled on the design of established youth
# surveys (the Youth Survey Luxembourg, the German Youth Institute's AID:A, and the
# youth sub-sample of the European Social Survey) and CALIBRATED to published
# Eurostat-style figures for 16-29-year-olds (NEET prevalence and life
# satisfaction), but it is not real microdata -- the real surveys require data-use
# agreements and cannot be redistributed inside a package. It follows the same
# convention as the package's other simulated datasets (`maihda_sim_data`,
# `maihda_long_data`, `maihda_sparse_data`): a KNOWN generative structure with a
# baked-in intersectional interaction, recorded as attr(maihda_youth_data, "truth")
# and documented in R/maihda_youth_data.R, so a vignette can claim *recovery* of a
# real interaction rather than merely report numbers.
#
# Design:
#   * 3 stratum dimensions -> 2 x 2 x 3 = 12 intersectional strata
#     (gender x migration x parental_edu)
#   * a higher-level grouping variable `country` (Luxembourg, Germany, France,
#     Italy) for the maihda(group = ) / compare_maihda_groups() workflow
#   * two outcomes built from one shared latent disadvantage so they correlate, as
#     in reality: `neet` (binary) and `wellbeing` (continuous)
#   * a true interaction share of ~33% of the between-stratum variance, on both
#     outcomes, concentrated in named compounding-(dis)advantage strata
#   * intersectional inequality AMPLIFIED in the higher-NEET countries, so the
#     VPC/ICC genuinely differs across the grouping variable
#
# Base R only (stats::uniroot calibrates the per-country NEET prevalence).
# Run from the package root:  Rscript data-raw/maihda_youth_data.R

set.seed(2026)

## ---- intersectional dimensions: 2 x 2 x 3 = 12 strata ----
gender_lv    <- c("Women", "Men")
migration_lv <- c("None", "Migration background")
parental_lv  <- c("Low", "Medium", "High")

grid <- expand.grid(gender = gender_lv, migration = migration_lv,
                    parental_edu = parental_lv, stringsAsFactors = FALSE)
S <- nrow(grid)                                            # 12
stratum_key <- apply(grid, 1, paste, collapse = ":")

## ---- additive main effects on a latent DISADVANTAGE scale ----
# Higher value = more disadvantaged (higher NEET risk, lower well-being). Social
# origin (parental education) carries the strongest, best-established gradient.
g_fx <- c(Women = 0.30, Men = -0.30)
m_fx <- c("None" = -0.50, "Migration background" = 0.50)
e_fx <- c(Low = 0.80, Medium = 0.00, High = -0.80)
raw_add <- g_fx[grid$gender] + m_fx[grid$migration] + e_fx[grid$parental_edu]

## ---- genuine intersectional departures beyond additivity ----
# The intersectional interaction is a genuine departure from the additive
# prediction, so by construction it must be ORTHOGONAL to the dimensions' main
# effects (otherwise the adjusted model simply absorbs it and there is no
# intersectionality to recover). We build two interpretable, textbook interaction
# contrasts on centred ("effect") codes -- which are exactly orthogonal to the main
# effects in this balanced 12-stratum grid -- and concentrate the signal in the
# policy-relevant corner: the migration penalty COMPOUNDS with low social origin,
# and again for young women.
gc <- c(Women = 0.5, Men = -0.5)[grid$gender]                       # gender, centred
mc <- c("None" = -0.5, "Migration background" = 0.5)[grid$migration] # migration, centred
pl <- c(Low = 1, Medium = 0, High = -1)[grid$parental_edu]          # social origin, centred
# migration x social-origin: the migration penalty is larger when parental
# education is low; gender x migration: larger again for women.
raw_int <- 1.00 * (mc * pl) + 0.60 * (gc * mc)

std <- function(x) (x - mean(x)) / stats::sd(x)
add0 <- std(raw_add)   # unit variance across the 12 strata
int0 <- std(raw_int)   # unit variance, orthogonal to the additive main effects
names(add0) <- names(int0) <- stratum_key

## ---- place the components on each outcome's scale to hit the target share ----
target_share <- 0.33   # interaction's share of the between-stratum variance

# Gaussian well-being (a life-satisfaction index, higher = better): between-stratum
# variance 0.50 at amplitude 1; residual SD 1.8 -> VPC ~ 0.13 at amplitude 1.
between_w <- 0.50
s_add_w <- sqrt(between_w * (1 - target_share))
s_int_w <- sqrt(between_w * target_share)
sigma_w <- 1.8

# Binary NEET (latent logit scale): between-stratum variance 1.00 at amplitude 1
# -> latent VPC ~ 0.23 at amplitude 1.
between_b <- 1.00
a_add <- sqrt(between_b * (1 - target_share))
a_int <- sqrt(between_b * target_share)

# Shared stratum-level signal (same shape on both outcomes -> NEET and low
# well-being correlate across strata).
disadv_w <- s_add_w * add0 + s_int_w * int0   # lowers well-being
disadv_b <- a_add  * add0 + a_int  * int0      # raises NEET log-odds
names(disadv_w) <- names(disadv_b) <- stratum_key

## ---- countries: the grouping variable, with calibrated prevalence + amplitude ----
country_lv    <- c("Luxembourg", "Germany", "France", "Italy")
n_per_country <- 750L

# Amplitude multiplier on the between-stratum signal: intersectional inequality is
# sharper in the higher-NEET countries, so the VPC genuinely differs across groups.
amp <- c(Luxembourg = 0.80, Germany = 0.90, France = 1.10, Italy = 1.40)

# Marginal NEET prevalence targets (Eurostat-style, ages 15-29); the per-country
# intercept is calibrated to these below.
neet_target <- c(Luxembourg = 0.07, Germany = 0.09, France = 0.13, Italy = 0.19)

# Mean well-being targets per country (life-satisfaction scale, higher = better).
wb_mean <- c(Luxembourg = 7.4, Germany = 7.3, France = 7.0, Italy = 6.8)

# Stratum-dimension composition is held CONSTANT across countries (same share with a
# migration background, same parental-education mix), so the cross-country VPC
# differences are driven by the intended `amp` (the between-stratum signal) rather
# than by confounding between composition and a country's baseline level.
mig_prob <- c(Luxembourg = 0.35, Germany = 0.35, France = 0.35, Italy = 0.35)

beta_age_neet <- -0.04   # older youth a little less likely to be NEET
beta_age_wb   <-  0.02   # and report marginally higher well-being

make_country <- function(cty) {
  n <- n_per_country
  gender    <- sample(gender_lv, n, replace = TRUE, prob = c(0.50, 0.50))
  migration <- sample(migration_lv, n, replace = TRUE,
                      prob = c(1 - mig_prob[[cty]], mig_prob[[cty]]))
  parental_edu <- sample(parental_lv, n, replace = TRUE, prob = c(0.30, 0.40, 0.30))
  age <- sample(16:29, n, replace = TRUE)
  key <- paste(gender, migration, parental_edu, sep = ":")
  age_c <- age - 22

  # Well-being: disadvantage lowers it; intersectional inequality scaled by amp.
  wb_eta <- wb_mean[[cty]] - amp[[cty]] * unname(disadv_w[key]) + beta_age_wb * age_c
  wellbeing <- wb_eta + stats::rnorm(n, 0, sigma_w)

  # NEET: calibrate the intercept so the marginal prevalence hits the target.
  lp_noint <- amp[[cty]] * unname(disadv_b[key]) + beta_age_neet * age_c
  b0 <- stats::uniroot(
    function(b) mean(stats::plogis(b + lp_noint)) - neet_target[[cty]],
    lower = -12, upper = 6)$root
  neet <- ifelse(stats::runif(n) < stats::plogis(b0 + lp_noint), "Yes", "No")

  data.frame(country = cty, gender = gender, migration = migration,
             parental_edu = parental_edu, age = age, neet = neet,
             wellbeing = round(wellbeing, 2), stringsAsFactors = FALSE)
}

maihda_youth_data <- do.call(rbind, lapply(country_lv, make_country))

## ---- assemble, identify, and factor-code with explicit level order ----
maihda_youth_data$id <- sprintf("Y%04d", seq_len(nrow(maihda_youth_data)))
maihda_youth_data <- maihda_youth_data[, c("id", "country", "gender", "migration",
                                           "parental_edu", "age", "neet", "wellbeing")]
maihda_youth_data$country      <- factor(maihda_youth_data$country, levels = country_lv)
maihda_youth_data$gender       <- factor(maihda_youth_data$gender, levels = gender_lv)
maihda_youth_data$migration    <- factor(maihda_youth_data$migration, levels = migration_lv)
maihda_youth_data$parental_edu <- factor(maihda_youth_data$parental_edu, levels = parental_lv)
maihda_youth_data$neet         <- factor(maihda_youth_data$neet, levels = c("No", "Yes"))
rownames(maihda_youth_data) <- NULL

## ---- record the ground truth (exact, because we built it) ----
truth <- list(
  target_interaction_share = target_share,
  dimensions = c("gender", "migration", "parental_edu"),
  n_strata = S,
  additive_main_effects = list(gender = g_fx, migration = m_fx, parental_edu = e_fx),
  interaction_form = paste("migration x social-origin + gender x migration",
                           "(orthogonal to the main effects); compounding",
                           "disadvantage concentrated in Women x Migration",
                           "background x Low parental education"),
  interaction_by_stratum = round(int0, 3),   # standardised interaction per stratum
  gaussian = list(between_var = between_w, residual_var = sigma_w^2,
                  interaction_share = target_share,
                  vpc = between_w / (between_w + sigma_w^2)),
  binary_latent = list(between_var = between_b, residual_var = pi^2 / 3,
                       interaction_share = target_share,
                       vpc = between_b / (between_b + pi^2 / 3)),
  country_amplitude = amp,         # multiplies the between-stratum variance by amp^2
  neet_prevalence_target = neet_target,
  wellbeing_mean_target = wb_mean,
  note = paste("VPCs above are at amplitude 1; the per-country between-stratum",
               "variance is amp^2 x the value shown, so the VPC differs by country.",
               "target_interaction_share is the construction target at the stratum",
               "level; the additive/interaction split a fitted model RECOVERS",
               "differs by outcome (the logit link attenuates it for the binary",
               "outcome, and pooling countries of different amplitude inflates it",
               "for the Gaussian one), so read the PCV the model reports, not 0.33.")
)
attr(maihda_youth_data, "truth") <- truth

## ---- report (paste these numbers into the roxygen doc) ----
cells <- table(apply(maihda_youth_data[, c("gender", "migration", "parental_edu")],
                     1, paste, collapse = ":"))
cat(sprintf("N=%d  strata=%d  median cell=%g  min cell=%g\n",
            nrow(maihda_youth_data), S,
            stats::median(as.numeric(cells)), min(as.numeric(cells))))
cat("NEET prevalence by country (target vs realised):\n")
print(rbind(target  = neet_target,
            realised = round(tapply(maihda_youth_data$neet == "Yes",
                                    maihda_youth_data$country, mean), 3)))
cat("Mean well-being by country (target vs realised):\n")
print(rbind(target  = wb_mean,
            realised = round(tapply(maihda_youth_data$wellbeing,
                                    maihda_youth_data$country, mean), 2)))

stopifnot(
  nrow(maihda_youth_data) == n_per_country * length(country_lv),
  !anyNA(maihda_youth_data[, c("country", "gender", "migration",
                               "parental_edu", "neet", "wellbeing")]),
  min(as.numeric(cells)) > 0   # every intersectional stratum is populated
)

usethis::use_data(maihda_youth_data, overwrite = TRUE, compress = "xz")
