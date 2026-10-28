# Read an element of a PCV result

`$` and `[[` read a
[`calculate_pcv`](https://hdbt.github.io/MAIHDA/reference/calculate_pcv.md)
result as a list. Reading its deprecated `pvc` element, a copy of `pcv`,
warns; read `pcv` instead.

## Usage

``` r
# S3 method for class 'pcv_result'
x$name

# S3 method for class 'pcv_result'
x[[i, ...]]
```

## Arguments

- x:

  A `pcv_result` object from
  [`calculate_pcv`](https://hdbt.github.io/MAIHDA/reference/calculate_pcv.md).

- name, i:

  The element's name (or, for `[[`, its position).

- ...:

  Passed on to the list method.

## Value

The element.
