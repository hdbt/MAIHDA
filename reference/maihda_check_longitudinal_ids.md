# Guard against longitudinal ids reused across strata

`(time | id)` treats every row sharing an id value as the SAME person,
so ids that are only unique within a site or group (person "1" in every
site) silently merge different people's trajectories into one level-2
unit. An id appearing in more than one stratum is the observable symptom
of that (a person's intersectional stratum is a person-level attribute,
so a genuinely unique id maps to exactly one stratum); reject it with
guidance rather than guessing a nesting. Called after strata resolution,
when the `stratum` column exists. Rows missing id or stratum are ignored
(the engines drop them).

## Usage

``` r
maihda_check_longitudinal_ids(data, id)
```

## Arguments

- data:

  The model data carrying the resolved `stratum` column.

- id:

  Name of the person/unit identifier column.

## Value

`NULL`, invisibly; stops on ambiguous ids.
