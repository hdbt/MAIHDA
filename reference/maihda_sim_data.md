# Simulated Health Data for MAIHDA Examples

A simulated dataset containing health outcomes and demographic variables
for 500 individuals. This dataset is designed to demonstrate
intersectional health inequalities suitable for MAIHDA analysis. The
data includes main effects and intersectional effects between gender,
race, and education.

## Usage

``` r
data("maihda_sim_data")
```

## Format

A data frame with 500 observations on the following 6 variables:

- `id`:

  Integer identifier for each individual

- `gender`:

  Character variable indicating gender ("Male" or "Female")

- `race`:

  Character variable indicating race/ethnicity ("White", "Black",
  "Hispanic", or "Asian")

- `education`:

  Character variable indicating education level ("High School", "Some
  College", "College", or "Graduate")

- `age`:

  Numeric variable for age in years (range: 18-80)

- `health_outcome`:

  Numeric variable representing a health score (higher is better)

## Details

The health outcome was simulated with:

- Main effects of gender, race, education, and age

- Intersectional effects (e.g., Black women, men with high school
  education)

- Random noise with standard deviation of 10

The data demonstrates typical patterns in health inequalities research
where outcomes vary both by individual characteristics and their
intersections.

## Examples

``` r
data(maihda_sim_data)

# View structure
str(maihda_sim_data)
#> 'data.frame':    500 obs. of  6 variables:
#>  $ id            : int  1 2 3 4 5 6 7 8 9 10 ...
#>  $ gender        : chr  "Male" "Male" "Male" "Female" ...
#>  $ race          : chr  "White" "White" "White" "White" ...
#>  $ education     : chr  "High School" "College" "High School" "Graduate" ...
#>  $ age           : num  63 44 51 48 43 44 57 43 21 43 ...
#>  $ health_outcome: num  78.9 86.4 75.1 88.3 90.3 72.2 65.5 83.3 68.6 81.1 ...

# Summary statistics
summary(maihda_sim_data)
#>        id           gender              race            education        
#>  Min.   :  1.0   Length:500         Length:500         Length:500        
#>  1st Qu.:125.8   Class :character   Class :character   Class :character  
#>  Median :250.5   Mode  :character   Mode  :character   Mode  :character  
#>  Mean   :250.5                                                           
#>  3rd Qu.:375.2                                                           
#>  Max.   :500.0                                                           
#>       age        health_outcome  
#>  Min.   :18.00   Min.   : 43.60  
#>  1st Qu.:37.00   1st Qu.: 67.80  
#>  Median :45.00   Median : 76.40  
#>  Mean   :45.31   Mean   : 75.95  
#>  3rd Qu.:53.00   3rd Qu.: 84.42  
#>  Max.   :80.00   Max.   :109.70  

if (FALSE) { # \dontrun{
# Example MAIHDA analysis
library(MAIHDA)

# Create strata
strata_result <- make_strata(maihda_sim_data, vars = c("gender", "race"))

# Fit model
model <- fit_maihda(health_outcome ~ age + (1 | stratum), 
                   data = strata_result$data)

# Summarize
summary_maihda(model)

# Visualize
plot_maihda(model, type = "caterpillar")
} # }
```
