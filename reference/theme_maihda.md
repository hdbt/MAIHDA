# MAIHDA plot theme

A ggplot2 theme that applies the MAIHDA brand identity to plot *chrome*
only – navy titles and axis labels and soft slate gridlines. It encodes
no data: the colourblind-safe data palettes used by the package's
plotting functions are left untouched.

## Usage

``` r
theme_maihda(base_size = 11, base_family = NULL)
```

## Arguments

- base_size:

  Base font size in points. Default \`11\`.

- base_family:

  Base font family. Defaults to \`getOption("maihda.font", "")\`, i.e.
  the graphics device default unless you have opted into a brand font
  globally. Pass a string to force a family for a single plot.

## Value

A ggplot2 theme object that can be added to a plot with \`+\`.

## Details

Built on \[ggplot2::theme_minimal()\], so it composes with additional
\`+ theme()\` calls in the usual way (later settings win).

## Fonts

By default the theme uses the graphics device's default font, which is
safe on every device (including the PostScript/PDF devices used by \`R
CMD check\` and many rendering back-ends). To render figures in the
brand font (Montserrat) – matching the hex logo – set
\`options(maihda.font = "Montserrat")\` or pass \`base_family =
"Montserrat"\`, and use a graphics device that can resolve that family
(e.g. ragg or showtext). Forcing a font that the active device cannot
find produces "invalid font type" errors, which is why it is opt-in
rather than automatic.

## Examples

``` r
library(ggplot2)
ggplot(mtcars, aes(mpg, wt)) +
  geom_point(colour = "#0072B2") +
  labs(title = "MAIHDA brand theme", x = "MPG", y = "Weight") +
  theme_maihda()
```
