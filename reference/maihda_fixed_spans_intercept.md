# Does a formula's fixed part carry the grand mean?

`TRUE` when the fixed-effect design the formula implies – built on
`data` the way the engines build it,
[`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) on the
bar-free formula – carries the all-ones (intercept) vector in its COLUMN
SPACE, whether or not the formula writes `1`. Writing `0 +` does not
settle it: with no intercept R codes the FIRST factor term as cell
means, so `y ~ 0 + f + x` still spans the grand mean, while `y ~ 0 + x`
with a numeric `x` does not. Membership is tested by QR rank – the
all-ones vector lies in colspace(X) exactly when rank(\[X, 1\]) ==
rank(X) – mirroring `maihda_re_lhs_spans_intercept()` for the random
part – but only for the `0 +` / `- 1` spellings: a formula that writes
an intercept returns `TRUE` without building a design, because
[`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) always
emits the all-ones column for one. Returns `NA` – never a guess – when
the design cannot be built or ranked (an absent or all-NA column, a
parse error), so neither caller rewrites a formula or warns about one on
the strength of the written `0 +` alone.

## Usage

``` r
maihda_fixed_spans_intercept(formula, data)
```

## Arguments

- formula:

  A model formula, with or without random-effect bars.

- data:

  Data frame the fixed design is built on.

## Value

`TRUE` if the fixed design spans the intercept, `FALSE` if it provably
does not, `NA` if the design could not be built.
