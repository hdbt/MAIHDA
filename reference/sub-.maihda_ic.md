# Subset MAIHDA information criteria

Indexing method for
[`maihda_ic`](https://hdbt.github.io/MAIHDA/reference/maihda_ic.md)
results that preserves the `ic_primary` metadata attribute (which names
the criterion the `delta` column is computed on) after a row/column
subset. Plain `[.data.frame` keeps the `maihda_ic` class but drops that
attribute, which would leave
[`print.maihda_ic()`](https://hdbt.github.io/MAIHDA/reference/print.maihda_ic.md)
testing a dropped (`NULL`) attribute as a scalar condition and erroring
after the table prints.

## Usage

``` r
# S3 method for class 'maihda_ic'
x[...]
```

## Arguments

- x:

  A `maihda_ic` object.

- ...:

  Row/column indices forwarded to `[.data.frame`.

## Value

A `maihda_ic` with `ic_primary` carried over, or a bare vector when the
selection drops to a single column.
