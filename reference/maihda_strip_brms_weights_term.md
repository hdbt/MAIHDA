# Remove the package-internal weights(.maihda_sw) addition term from a formula

The inverse of the injection
[`maihda_brms_weights_formula()`](https://hdbt.github.io/MAIHDA/reference/maihda_brms_weights_formula.md)
performs for the reserved `.maihda_sw` column.
[`maihda()`](https://hdbt.github.io/MAIHDA/reference/maihda.md) and
[`compare_maihda_groups()`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
derive the null / adjusted models from a sampling-weighted fit's
*stored* formula, which already carries the injected
`weights(.maihda_sw)` term; that term must be stripped before the
derived model is re-prepared, or the reserved-column guard (`.maihda_sw`
both in the formula and in the carried-over data) and the
weights-formula rewrite (which rejects an existing
[`weights()`](https://rdrr.io/r/stats/weights.html) term) abort the
package's own refit. Only the internal `weights(.maihda_sw)` call is
removed – a user's own
[`weights()`](https://rdrr.io/r/stats/weights.html) term (a genuine
conflict) is left in place so it is still rejected – and any other
addition term (e.g. `trials(n)`) is preserved.

## Usage

``` r
maihda_strip_brms_weights_term(formula)
```

## Arguments

- formula:

  A model formula.

## Value

The formula with the internal `weights(.maihda_sw)` term removed, or the
input unchanged when it carries no such term.
