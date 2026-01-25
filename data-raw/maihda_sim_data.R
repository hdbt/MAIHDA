# Script to generate example dataset for MAIHDA package
set.seed(123)

# Create a simulated health dataset with intersectional categories
n <- 500

# Generate categorical variables
gender <- sample(c("Male", "Female"), n, replace = TRUE)
race <- sample(c("White", "Black", "Hispanic", "Asian"), n, replace = TRUE, 
               prob = c(0.6, 0.2, 0.15, 0.05))
education <- sample(c("High School", "Some College", "College", "Graduate"), n, 
                   replace = TRUE, prob = c(0.3, 0.25, 0.3, 0.15))

# Generate continuous variables
age <- round(rnorm(n, mean = 45, sd = 12))
age <- pmax(18, pmin(80, age))  # Constrain between 18 and 80

# Create intersectional effects
# Base health outcome
base_health <- 70

# Main effects
gender_effect <- ifelse(gender == "Female", 2, 0)

race_effect <- numeric(n)
race_effect[race == "White"] <- 5
race_effect[race == "Asian"] <- 3
race_effect[race == "Hispanic"] <- -2
race_effect[race == "Black"] <- -3

education_effect <- numeric(n)
education_effect[education == "Graduate"] <- 8
education_effect[education == "College"] <- 5
education_effect[education == "Some College"] <- 2
education_effect[education == "High School"] <- 0

age_effect <- -0.2 * (age - 45)

# Add some intersectional effects (interaction terms)
intersect_effect <- rep(0, n)
intersect_effect[gender == "Female" & race == "Black"] <- -5
intersect_effect[gender == "Male" & education == "High School"] <- -3
intersect_effect[race == "Hispanic" & education %in% c("College", "Graduate")] <- 4

# Random noise
random_noise <- rnorm(n, mean = 0, sd = 10)

# Combine effects
health_outcome <- base_health + gender_effect + race_effect + education_effect + 
                 age_effect + intersect_effect + random_noise

# Round to 1 decimal place
health_outcome <- round(health_outcome, 1)

# Create the dataset
maihda_sim_data <- data.frame(
  id = 1:n,
  gender = gender,
  race = race,
  education = education,
  age = age,
  health_outcome = health_outcome,
  stringsAsFactors = FALSE
)

# Save the dataset (would use usethis::use_data in actual development)
# For now, save as RData file manually
save(maihda_sim_data, file = "data/maihda_sim_data.rda", compress = "xz")
