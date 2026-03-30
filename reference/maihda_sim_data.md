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

# \donttest{
# Example MAIHDA analysis
library(MAIHDA)

# Create strata
strata_result <- make_strata(maihda_sim_data, vars = c("gender", "race"))

# Fit model
model <- fit_maihda(health_outcome ~ age + (1 | stratum), 
                   data = strata_result$data)

# Summarize
summary_maihda(model)
#> MAIHDA Model Summary
#> ====================
#> 
#> Variance Partition Coefficient (VPC/ICC):
#>   Estimate: 0.1907
#> 
#> Variance Components:
#>                  component variance     sd proportion
#>   Between-stratum (random)    26.27  5.125     0.1907
#>  Within-stratum (residual)   111.49 10.559     0.8093
#>                      Total   137.76 11.737     1.0000
#> 
#> Fixed Effects:
#>         term estimate
#>  (Intercept)  84.9021
#>          age  -0.2647
#> 
#> Stratum Estimates (first 10):
#>  stratum stratum_id           label random_effect     se lower_95 upper_95
#>        1          1    Female_Asian        -1.310 3.1488   -7.482   4.8613
#>        2          2    Female_Black        -6.536 1.4336   -9.346  -3.7263
#>        3          3 Female_Hispanic         1.357 1.7539   -2.081   4.7944
#>        4          4    Female_White         7.942 0.8502    6.275   9.6080
#>        5          5      Male_Asian        -2.978 2.5427   -7.961   2.0059
#>        6          6      Male_Black        -3.114 1.5202   -6.094  -0.1345
#>        7          7   Male_Hispanic        -0.916 1.4079   -3.676   1.8435
#>        8          8      Male_White         5.556 0.8447    3.901   7.2119

# Visualize
plot_maihda(model, type = "caterpillar")

# }
```
