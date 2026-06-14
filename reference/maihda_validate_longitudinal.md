# Validate a longitudinal (id / time) specification

Validate a longitudinal (id / time) specification

## Usage

``` r
maihda_validate_longitudinal(
  id,
  time,
  time_degree,
  data,
  engine = "lme4",
  sampling_weights = NULL,
  context = NULL
)
```

## Arguments

- id:

  Single column name: the person/unit identifier (level 2).

- time:

  Single column name: a numeric measurement-time variable (level 1).

- time_degree:

  Integer \>= 1: polynomial degree of the growth curve (1 = linear).
  brms supports degree 1 only.

- data:

  The model data.

- engine:

  The fitting engine; only "lme4"/"brms" support the 3-level growth
  structure.

- sampling_weights, context:

  Must be NULL – design-weighted and contextual longitudinal models are
  out of scope.

## Value

A list `list(id, time, time_degree)`.
