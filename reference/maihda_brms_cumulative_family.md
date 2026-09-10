# Cumulative brmsfamily for the brms ordinal engine

Turns whatever the caller supplied for an ordinal `engine = "brms"` fit
into a proper
[`brms::cumulative()`](https://paulbuerkner.com/brms/reference/brmsfamily.html)
family, preserving the cumulative-specific options a user set.

## Usage

``` r
maihda_brms_cumulative_family(family)
```

## Arguments

- family:

  The family carried into the brms branch of
  [`fit_maihda()`](https://hdbt.github.io/MAIHDA/reference/fit_maihda.md):
  a `brmsfamily`, or the plain cumulative marker list.

## Value

A
[`brms::cumulative()`](https://paulbuerkner.com/brms/reference/brmsfamily.html)
family object.

## Details

A rebuild is needed because `family` may still be the plain
`list(family = "cumulative", link = )` marker that the family-string
path and
[`maihda_cumulative`](https://hdbt.github.io/MAIHDA/reference/maihda_cumulative.md)
produce, which `brm()` cannot use. It must be a rebuild rather than a
pass-through for the same reason. Rebuilding from the link ALONE,
however, silently reset every other option to brms's defaults:
`brms::cumulative(threshold = "equidistant")` fitted the flexible model,
with no warning that a different model had been fitted from the one
asked for. `threshold` and `link_disc` are therefore carried across when
the supplied family actually has them; the marker lists carry neither,
so those paths still get brms's defaults.

Note the contrast with the `clmm` engine, where a structured threshold
also has to be handled on the way OUT
([`maihda_clmm_cutpoints`](https://hdbt.github.io/MAIHDA/reference/maihda_clmm_cutpoints.md)).
brms needs no such treatment: its generated quantities expand the
constrained thresholds into the full `b_Intercept[1..K-1]` vector, so
[`brms::fixef()`](https://rdrr.io/pkg/nlme/man/fixed.effects.html)
reports all \\K-1\\ cut points under every threshold structure and
[`maihda_brms_ordinal_thresholds`](https://hdbt.github.io/MAIHDA/reference/maihda_brms_ordinal_thresholds.md)
reads them unchanged.
