# Subset a MAIHDA group comparison

Indexing method for
[`compare_maihda_groups`](https://hdbt.github.io/MAIHDA/reference/compare_maihda_groups.md)
results that preserves the comparison's metadata attributes (group
variable, engine, family, ...) so the
[`print()`](https://rdrr.io/r/base/print.html) header stays populated
after a row/column subset. Plain `[.data.frame` keeps the S3 class but
drops those attributes, which would otherwise blank the printed header.

## Usage

``` r
# S3 method for class 'maihda_group_comparison'
x[...]
```

## Arguments

- x:

  A `maihda_group_comparison` object.

- ...:

  Row/column indices forwarded to `[.data.frame`.

## Value

A `maihda_group_comparison` with the metadata attributes carried over,
or a bare vector when the selection drops to a single column.
