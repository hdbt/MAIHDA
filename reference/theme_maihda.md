# MAIHDA plot theme

A ggplot2 theme that applies the MAIHDA brand identity to plot *chrome*
only – navy titles and axis labels, soft slate gridlines, and the brand
font (Montserrat) when it is available, with a safe fallback to the
device default font otherwise. It encodes no data: the colourblind-safe
data palettes used by the package's plotting functions are left
untouched.

## Usage

``` r
theme_maihda(base_size = 11, base_family = NULL)
```

## Arguments

- base_size:

  Base font size in points. Default \`11\`.

- base_family:

  Base font family. When \`NULL\` (the default), the brand font
  "Montserrat" is used if it is registered with the graphics device,
  otherwise the device default. Pass a string to force a family, or set
  \`options(maihda.font = "...")\` to override globally.

## Value

A ggplot2 theme object that can be added to a plot with \`+\`.

## Details

Built on \[ggplot2::theme_minimal()\], so it composes with additional
\`+ theme()\` calls in the usual way (later settings win).

## Examples

``` r
library(ggplot2)
#> 
#> Attaching package: ‘ggplot2’
#> The following objects are masked from ‘package:ggtern’:
#> 
#>     aes, annotate, ggplot, ggplotGrob, ggplot_build, ggplot_gtable,
#>     ggsave, layer_data, theme_bw, theme_classic, theme_dark,
#>     theme_gray, theme_light, theme_linedraw, theme_minimal, theme_void
ggplot(mtcars, aes(mpg, wt)) +
  geom_point(colour = "#0072B2") +
  labs(title = "MAIHDA brand theme", x = "MPG", y = "Weight") +
  theme_maihda()
```
