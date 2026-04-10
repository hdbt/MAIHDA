library(hexSticker)
library(ggplot2)
library(showtext)
library(sf)

font_add_google("Montserrat", "montserrat")
showtext_auto()

set.seed(404)

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

# Base colors (approx. 70% saturation, high value)
col1 <- "#E64585" # Pink
col2 <- "#45E6E6" # Cyan
col3 <- "#E6C545" # Yellow

# Function to multiply hex colors (simulates interaction effects)
multiply_colors <- function(...) {
  hexes <- list(...)
  rgbs <- lapply(hexes, function(x) col2rgb(x) / 255)
  blended <- Reduce(`*`, rgbs)
  rgb(blended[1,1], blended[2,1], blended[3,1])
}

col_12 <- multiply_colors(col1, col2)
col_13 <- multiply_colors(col1, col3)
col_23 <- multiply_colors(col2, col3)
col_all <- multiply_colors(col1, col2, col3)

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
        p_y = 1.45,
        s_x = 1,
        s_y = 0.85,
        s_width = 1.4,
        s_height = 1.1,
        h_fill = "#2C3E50",             # Dark blue background
        h_color = "#FFFFFF",            # White border
        h_size = 1.2,
        filename = "man/figures/logo.png",
        dpi = 600)

cat("Blue background hex sticker saved\n")
