library(hexSticker)
library(ggplot2)
library(showtext)

font_add_google("Montserrat", "montserrat")
showtext_auto()

set.seed(404)

circles <- data.frame(
  x = c(-0.5, 0.5, 0),
  y = c(0.3, 0.3, -0.4),
  r = rep(1, 3),
  group = factor(1:3)
)

circle_points <- function(center_x, center_y, r, n_points = 100) {
  theta <- seq(0, 2*pi, length.out = n_points)
  data.frame(
    x = center_x + r * cos(theta),
    y = center_y + r * sin(theta)
  )
}

circle_data <- do.call(rbind, lapply(1:3, function(i) {
  pts <- circle_points(circles$x[i], circles$y[i], circles$r[i])
  pts$group <- factor(i)
  pts
}))

p <- ggplot(circle_data, aes(x = x, y = y, group = group, fill = group)) +
  geom_polygon(alpha = 0.6, color = "white", size = 1.2) +
  scale_fill_manual(values = c("#FFD93D", "#6BCF7F", "#4D9DE0")) +
  theme_void() +
  theme(legend.position = "none",
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA)) +
  coord_fixed(xlim = c(-2, 2), ylim = c(-1.8, 1.8))

sticker(p,
        package = "MAIHDA",
        p_size = 22,
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
