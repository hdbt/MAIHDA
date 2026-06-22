library(hexSticker)
library(ggplot2)
library(showtext)
library(sf)

font_add_google("Montserrat", "montserrat")
showtext_auto()

# Define circle centers and radius
centers <- data.frame(
  x = c(-0.5, 0.5, 0),
  y = c(0.3, 0.3, -0.4),
  r = c(1, 1, 1)
)

# Create sf polygons for the circles
circles_sf <- lapply(1:3, function(i) {
  st_buffer(st_point(as.numeric(centers[i, 1:2])), dist = centers$r[i])
})
circles_sfc <- st_sfc(circles_sf)

# Find the intersection of all three circles
intersection_all <- st_intersection(st_intersection(circles_sfc[[1]], circles_sfc[[2]]), circles_sfc[[3]])

# (Optional) Find pairwise intersections, subtracting the all-intersection center
int_12 <- st_difference(st_intersection(circles_sfc[[1]], circles_sfc[[2]]), intersection_all)
int_13 <- st_difference(st_intersection(circles_sfc[[1]], circles_sfc[[3]]), intersection_all)
int_23 <- st_difference(st_intersection(circles_sfc[[2]], circles_sfc[[3]]), intersection_all)

# Base colors: bright CMY-style primaries, one per social dimension.
col1 <- "#E64585" # Pink / magenta
col2 <- "#45E6E6" # Cyan
col3 <- "#E6C545" # Yellow

# Intersectionality is the *non-additive* part, so overlaps are MULTIPLICATIVE
# (a subtractive product), not additive. Each pairwise overlap is the product of
# its two parents -> an emergent secondary (pink x cyan = blue, pink x yellow =
# red, cyan x yellow = green) that is not a sum of the parts.
multiply_colors <- function(...) {
  hexes <- list(...)
  rgbs <- lapply(hexes, function(x) col2rgb(x) / 255)
  blended <- Reduce(`*`, rgbs)
  rgb(blended[1, 1], blended[2, 1], blended[3, 1])
}

# Lift a colour's value (HSV) to a legibility floor, preserving hue/saturation,
# so dark products still read against the navy background.
lift_value <- function(hex, min_v = 0.62) {
  h <- rgb2hsv(col2rgb(hex))
  hsv(h[1, 1], h[2, 1], max(h[3, 1], min_v))
}

col_12 <- lift_value(multiply_colors(col1, col2)) # pink x cyan   -> blue
col_13 <- lift_value(multiply_colors(col1, col3)) # pink x yellow -> red
col_23 <- lift_value(multiply_colors(col2, col3)) # cyan x yellow -> green

# The full intersection is the irreducible intersectional effect: an emergent
# quantity NOT predictable from the parts (a literal 3-way product collapses to
# a muddy near-black), so it is set deliberately as a luminous hero colour.
col_all <- "#7E33A8" # Emergent violet

# Plot
p <- ggplot() +
  # Base circles
  geom_sf(data = circles_sfc[[1]], fill = col1, color = "white", linewidth = 1.2) +
  geom_sf(data = circles_sfc[[2]], fill = col2, color = "white", linewidth = 1.2) +
  geom_sf(data = circles_sfc[[3]], fill = col3, color = "white", linewidth = 1.2) +
  # Overwrite pairwise intersections with distinct emergent colors
  geom_sf(data = int_12, fill = col_12, color = "white", linewidth = 1) +
  geom_sf(data = int_13, fill = col_13, color = "white", linewidth = 1) +
  geom_sf(data = int_23, fill = col_23, color = "white", linewidth = 1) +
  # Overwrite the central intersection with an emergent color (e.g., bright pink/magenta)
  geom_sf(data = intersection_all, fill = col_all, color = "white", linewidth = 1.2) +
  theme_void() +
  theme(legend.position = "none",
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA)) +
  coord_sf(xlim = c(-2, 2), ylim = c(-1.8, 1.8), expand = FALSE)

sticker(p,
        package = "MAIHDA",
        p_size = 30,
        p_color = "#FFFFFF",            # White text
        p_family = "montserrat",
        p_fontface = "bold",
        p_y = 1.5,
        s_x = 1,
        s_y = 0.82,
        s_width = 1.5,                  # enlarge Venn to fill the body
        s_height = 1.2,
        url = "hdbt.github.io/MAIHDA",  # fills the empty lower point
        u_size = 5.5,
        u_color = "#AEB6BF",            # soft slate-grey, doesn't fight the title
        u_family = "montserrat",
        u_x = 1,
        u_y = 0.08,
        h_fill = "#2C3E50",             # Dark blue background
        h_color = "#FFFFFF",            # White border
        h_size = 1.6,                   # thicker so the silhouette survives small
        filename = "man/figures/logo.png",
        dpi = 600)

cat("Blue background hex sticker saved\n")
