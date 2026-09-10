# Warn when a MAIHDA fit's fixed part cannot represent the outcome's mean

A MAIHDA reads the stratum random intercept's variance as
between-stratum inequality. When the fixed design does not span the
intercept – a `0 +` or `- 1` formula whose surviving terms are all
numeric, so no factor is coded as cell means – there is nothing in the
fixed part to carry the outcome's mean and the random intercept absorbs
it instead. The between-stratum variance, and with it the VPC, the MOR
and any PCV computed from the fit, then measure the location of the
outcome rather than the strata.

## Usage

``` r
maihda_warn_no_grand_mean(formula, data_list)
```

## Arguments

- formula:

  The fitted (resolved) model formula.

- data_list:

  List of candidate data frames; the first on which the fixed design can
  be built decides. The analytic model frame may name only derived
  columns (`scale(x)`), the pre-fit frame only raw ones, so both are
  offered.

## Value

Invisibly, `TRUE` if a warning was emitted.

## Details

Silent unless the span test can be run and provably fails, so a formula
whose design could not be built is never warned about on the strength of
a written `0 +`.
