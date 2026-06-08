# Join non-empty caption lines for a plot

Pastes its arguments into a single newline-separated caption, dropping
any that are NULL, NA, or empty. Returns NULL when nothing remains, so a
plot with no caption content keeps a clean (caption-free) look.

## Usage

``` r
maihda_compose_caption(...)
```

## Arguments

- ...:

  Character scalars (or NULL).

## Value

A single newline-joined string, or NULL.
