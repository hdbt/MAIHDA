# Simulate the bundled longitudinal MAIHDA dataset `maihda_long_data`.
#
# A long-format (one row per person-occasion) panel for demonstrating
# longitudinal / growth-curve MAIHDA: repeated measurements of a wellbeing score
# over five waves for individuals in 12 intersectional strata
# (gender x ethnicity x education). The between-stratum trajectory differences are
# constructed to be MOSTLY ADDITIVE (gender/ethnicity/education main effects on the
# intercept and slope) with one genuine multiplicative interaction, so the
# longitudinal PCV (high PCV_slope, but < 1) is demonstrable. No external packages
# are required -- this uses base R only.
#
# Run from the package root:  Rscript data-raw/maihda_long_data.R

set.seed(2024)

n_persons <- 600L
waves <- 0:4

gender_lv <- c("Women", "Men")
ethnicity_lv <- c("EthA", "EthB", "EthC")
education_lv <- c("Low", "High")

# Per-person (time-invariant) characteristics.
person <- data.frame(
  id = sprintf("P%04d", seq_len(n_persons)),
  gender = sample(gender_lv, n_persons, replace = TRUE, prob = c(0.52, 0.48)),
  ethnicity = sample(ethnicity_lv, n_persons, replace = TRUE,
                     prob = c(0.5, 0.3, 0.2)),
  education = sample(education_lv, n_persons, replace = TRUE, prob = c(0.55, 0.45)),
  age = round(stats::runif(n_persons, 30, 60)),
  stringsAsFactors = FALSE
)

# Additive stratum effects on the baseline level (intercept) and the rate of
# change (slope on wave): each dimension contributes independently.
g_int <- c(Women = 0.20, Men = -0.20)
e_int <- c(EthA = 0.30, EthB = 0.00, EthC = -0.30)
d_int <- c(Low = -0.40, High = 0.40)
g_slp <- c(Women = -0.05, Men = 0.05)
e_slp <- c(EthA = 0.05, EthB = 0.00, EthC = -0.05)
d_slp <- c(Low = -0.08, High = 0.08)

add_int <- g_int[person$gender] + e_int[person$ethnicity] + d_int[person$education]
add_slp <- g_slp[person$gender] + e_slp[person$ethnicity] + d_slp[person$education]

# Genuine intersectional interactions beyond additive: a few specific strata
# deviate from the additive prediction on BOTH the baseline level and the slope.
# This is the multiplicative part the adjusted model cannot absorb, so the PCV is
# clearly below 1 (and the adjusted fit retains some between-stratum variance --
# it is not a degenerate, singular boundary fit). Trajectory differences remain
# mostly additive, so PCV_slope stays high.
cell <- paste(person$gender, person$ethnicity, person$education, sep = "_")
int_int <- c("Women_EthC_Low" = -0.35, "Men_EthA_High" = 0.30,
             "Women_EthB_High" = -0.25)
int_slp <- c("Women_EthC_Low" = -0.30, "Men_EthA_High" = 0.18,
             "Women_EthB_High" = 0.12)
add_int <- add_int + ifelse(cell %in% names(int_int), int_int[cell], 0)
add_slp <- add_slp + ifelse(cell %in% names(int_slp), int_slp[cell], 0)

# Person-level random intercepts and slopes (within-stratum heterogeneity).
v0 <- stats::rnorm(n_persons, 0, 0.50)
v1 <- stats::rnorm(n_persons, 0, 0.15)

grand_mean <- 5.0
beta_age <- 0.010   # per year, centred below
age_c <- person$age - mean(person$age)

rows <- vector("list", length(waves))
for (k in seq_along(waves)) {
  w <- waves[k]
  eta <- grand_mean + beta_age * age_c +
    (add_int + v0) + (add_slp + v1) * w
  y <- eta + stats::rnorm(n_persons, 0, 0.60)
  rows[[k]] <- data.frame(
    id = person$id,
    wave = w,
    gender = person$gender,
    ethnicity = person$ethnicity,
    education = person$education,
    age = person$age,
    wellbeing = round(y, 3),
    stringsAsFactors = FALSE
  )
}

maihda_long_data <- do.call(rbind, rows)
# A binary companion outcome (low wellbeing) for the logistic longitudinal path.
maihda_long_data$low_wellbeing <- as.integer(
  maihda_long_data$wellbeing < stats::quantile(maihda_long_data$wellbeing, 0.40))

# Order by person then wave (natural long format) and reset row names.
maihda_long_data <- maihda_long_data[order(maihda_long_data$id, maihda_long_data$wave), ]
rownames(maihda_long_data) <- NULL

# Factor coding with explicit, reproducible level order.
maihda_long_data$gender <- factor(maihda_long_data$gender, levels = gender_lv)
maihda_long_data$ethnicity <- factor(maihda_long_data$ethnicity, levels = ethnicity_lv)
maihda_long_data$education <- factor(maihda_long_data$education, levels = education_lv)

stopifnot(
  nrow(maihda_long_data) == n_persons * length(waves),
  anyDuplicated(maihda_long_data$id) > 0  # genuinely repeated measures
)

save(maihda_long_data, file = "data/maihda_long_data.rda", compress = "xz")
message("Wrote data/maihda_long_data.rda (", nrow(maihda_long_data), " rows, ",
        n_persons, " persons x ", length(waves), " waves).")
