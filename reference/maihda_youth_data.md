# Simulated Cross-National Youth Data for MAIHDA

A simulated, cross-sectional sample of 16–29-year-olds in four European
countries, for demonstrating MAIHDA on the questions youth-research
institutes study: the school-to-work transition (**NEET** status – Not
in Employment, Education or Training) and youth **subjective
well-being**, at the intersection of **gender, migration background, and
social origin** (parental education). It supports the two-model
decomposition
([`maihda`](https://hdbt.github.io/MAIHDA/reference/maihda.md)), the
binary discriminatory-accuracy workflow (`family = "binomial"`), and the
cross-group comparison (`maihda(group = "country")`).

## Usage

``` r
maihda_youth_data
```

## Format

A data frame with 3000 rows (750 per country) and 8 variables:

- id:

  Person identifier.

- country:

  Country (`Luxembourg`/`Germany`/`France`/`Italy`); the higher-level
  grouping variable for `maihda(group = "country")`.

- gender:

  Gender (`Women`/`Men`); a stratum dimension.

- migration:

  Migration background (`None`/`Migration background`); a stratum
  dimension.

- parental_edu:

  Parental education / social origin (`Low`/`Medium`/`High`); a stratum
  dimension.

- age:

  Age in years (16–29); an individual-level covariate.

- neet:

  NEET status (`No`/`Yes`); the binary outcome.

- wellbeing:

  Subjective well-being (a life-satisfaction index, higher is better);
  the continuous outcome.

## Source

Simulated for the MAIHDA package; **not real microdata**. The design
(age band, dimensions, and outcomes) follows established youth surveys –
the Youth Survey Luxembourg
(<https://www.youth-in-luxembourg.lu/project/youth-survey-luxembourg-ysl/>),
the German Youth Institute's AID:A
(<https://www.dji.de/en/about-us/projects/projekte/aida.html>), and the
youth sub-sample of the European Social Survey – and the NEET
prevalences and life-satisfaction levels are calibrated to Eurostat
figures for 15–29-year-olds
(<https://ec.europa.eu/eurostat/web/education-and-training/data>). The
real surveys require data-use agreements and cannot be redistributed in
a package, so the bundled data are synthetic, following the same
convention as the package's other simulated datasets
([`maihda_sim_data`](https://hdbt.github.io/MAIHDA/reference/maihda_sim_data.md),
[`maihda_long_data`](https://hdbt.github.io/MAIHDA/reference/maihda_long_data.md),
[`maihda_sparse_data`](https://hdbt.github.io/MAIHDA/reference/maihda_sparse_data.md)).

## Details

The 3 stratum dimensions form 2 x 2 x 3 = 12 intersectional strata,
evaluated within each of 4 countries (the higher-level grouping
variable). Two outcomes are generated from one shared latent
disadvantage, so they correlate as in reality. The between-stratum
differences are mostly additive (the dimensions' main effects) with a
genuine intersectional interaction – the migration penalty is steeper at
low social origin, and again for young women – built orthogonal to the
main effects so it survives adjustment and is recoverable (a sub-1 PCV).
The intensity of intersectional inequality is amplified in the
higher-NEET countries, so the VPC/ICC genuinely differs across the
grouping variable.

The exact generative parameters are stored as
`attr(maihda_youth_data, "truth")`: the additive main effects, the
standardised interaction per stratum (`interaction_by_stratum`), the
between-stratum variances and VPCs (at amplitude 1), the per-country
amplitude multipliers, and the calibration targets. The split a fitted
model *recovers* (the PCV) differs from the construction target by
outcome, because the logit link attenuates the interaction for the
binary outcome and pooling countries of different amplitude inflates it
for the Gaussian one – read the PCV the model reports rather than the
nominal target.

## See also

[`maihda`](https://hdbt.github.io/MAIHDA/reference/maihda.md),
[`maihda_discriminatory_accuracy`](https://hdbt.github.io/MAIHDA/reference/maihda_discriminatory_accuracy.md),
[`maihda_interactions`](https://hdbt.github.io/MAIHDA/reference/maihda_interactions.md)

## Examples

``` r
data(maihda_youth_data)
# \donttest{
# Well-being: the two-model VPC/ICC + PCV (additive vs intersectional):
wb <- maihda(wellbeing ~ gender + migration + parental_edu +
               (1 | gender:migration:parental_edu), data = maihda_youth_data)
wb
#> MAIHDA Analysis
#> ===============
#> 
#> Null formula:    wellbeing ~ (1 | stratum)
#> Adjusted formula:wellbeing ~ gender + migration + parental_edu + (1 | stratum)
#> Engine: lme4 | Family: gaussian
#> VPC/ICC (null): 0.1394
#> PCV (null -> adjusted): 0.7192
#> Between-stratum variance: 0.4576 (null) -> 0.1285 (adjusted)
#>   ~71.9% of the between-stratum variance is additive (the dimensions' main
#>   effects); the remainder is the between-stratum variance remaining after the
#>   additive main effects -- a model-dependent quantity
#> Strata: 12
#> 
#> Use summary() for variance components and plot(type = ...) for figures.

# NEET: the binary discriminatory-accuracy workflow (AUC / MOR ride along):
neet <- maihda(neet ~ gender + migration + parental_edu +
                 (1 | gender:migration:parental_edu),
               data = maihda_youth_data, family = "binomial")
#> Binary outcome 'neet' recoded to 0/1: 'No' = 0 (reference), 'Yes' = 1 (modeled event). Set the factor levels (or supply a 0/1 outcome) to control which level is the event.
neet
#> MAIHDA Analysis
#> ===============
#> 
#> Null formula:    neet ~ (1 | stratum)
#> Adjusted formula:neet ~ gender + migration + parental_edu + (1 | stratum)
#> Engine: lme4 | Family: binomial
#> VPC/ICC (null): 0.2033
#> PCV (null -> adjusted): 0.7464
#> Between-stratum variance: 0.8396 (null) -> 0.2129 (adjusted)
#>   ~74.6% of the between-stratum variance is additive (the dimensions' main
#>   effects); the remainder is the between-stratum variance remaining after the
#>   additive main effects -- a model-dependent quantity
#> 
#> Discriminatory accuracy (null model):
#>   AUC: 0.727 | MOR: 2.397 | cases/controls: 382/2618
#>   Adjusted-model AUC: 0.723
#> Strata: 12
#> 
#> Use summary() for variance components and plot(type = ...) for figures.

# Compare intersectional inequality across countries:
by_country <- maihda(neet ~ gender + migration + parental_edu +
                       (1 | gender:migration:parental_edu),
                     data = maihda_youth_data, group = "country",
                     family = "binomial")
#> Binary outcome 'neet' recoded to 0/1: 'No' = 0 (reference), 'Yes' = 1 (modeled event). Set the factor levels (or supply a 0/1 outcome) to control which level is the event.
#> Binary outcome 'neet' recoded to 0/1: 'No' = 0 (reference), 'Yes' = 1 (modeled event). Set the factor levels (or supply a 0/1 outcome) to control which level is the event.
#> Binary outcome 'neet' recoded to 0/1: 'No' = 0 (reference), 'Yes' = 1 (modeled event). Set the factor levels (or supply a 0/1 outcome) to control which level is the event.
#> Binary outcome 'neet' recoded to 0/1: 'No' = 0 (reference), 'Yes' = 1 (modeled event). Set the factor levels (or supply a 0/1 outcome) to control which level is the event.
#> Binary outcome 'neet' recoded to 0/1: 'No' = 0 (reference), 'Yes' = 1 (modeled event). Set the factor levels (or supply a 0/1 outcome) to control which level is the event.
by_country
#> MAIHDA Analysis
#> ===============
#> 
#> Null formula:    neet ~ (1 | stratum)
#> Adjusted formula:neet ~ gender + migration + parental_edu + (1 | stratum)
#> Engine: lme4 | Family: binomial
#> VPC/ICC (null): 0.2033
#> PCV (null -> adjusted): 0.7464
#> Between-stratum variance: 0.8396 (null) -> 0.2129 (adjusted)
#>   ~74.6% of the between-stratum variance is additive (the dimensions' main
#>   effects); the remainder is the between-stratum variance remaining after the
#>   additive main effects -- a model-dependent quantity
#> 
#> Discriminatory accuracy (null model):
#>   AUC: 0.727 | MOR: 2.397 | cases/controls: 382/2618
#>   Adjusted-model AUC: 0.723
#> Strata: 12
#> 
#> Group comparison by 'country':
#> MAIHDA Group Comparison
#> =======================
#> 
#> Group variable: country 
#> Engine: lme4  | Family: binomial  | Strata: shared/global 
#> 
#>       group   n n_strata    vpc var_between var_other var_residual    pcv
#>      France 750       12 0.2535      1.1169         0         3.29 0.9057
#>     Germany 750       12 0.1310      0.4958         0         3.29 0.6344
#>       Italy 750       12 0.3570      1.8265         0         3.29 0.8311
#>  Luxembourg 750       12 0.1568      0.6119         0         3.29 0.9338
#>  var_between_adjusted status
#>               0.10532     ok
#>               0.18128     ok
#>               0.30859     ok
#>               0.04051     ok
#> 
#> Use summary() for variance components and plot(type = ...) for figures.
# }
```
