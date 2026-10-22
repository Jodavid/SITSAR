## Octagonal logo of the SITSAR package (run from the package root)
library(grid)
load("data/SanFrancisco150.rda")

# Pauli RGB of the San Francisco image
x <- SanFrancisco150
C11 <- Re(x[, , 1]); C22 <- Re(x[, , 2]); C33 <- Re(x[, , 3]); C13 <- Re(x[, , 5])
sc <- function(v) { v <- v / quantile(v, 0.985); v[v > 1] <- 1; v^0.55 }
pauli <- array(c(sc((C11 + C33 - 2 * C13) / 2), sc(2 * C22), sc((C11 + C33 + 2 * C13) / 2)),
               c(150, 150, 3))

octagon <- function(r, cx = 0.5, cy = 0.5) {
  a <- pi / 8 + (0:7) * pi / 4
  list(x = cx + r * cos(a), y = cy + r * sin(a))
}

draw_logo <- function() {
  grid.newpage()
  o_out <- octagon(0.495); o_mid <- octagon(0.462); o_in <- octagon(0.445)
  # outer ring (teal -> blue gradient)
  grid.polygon(o_out$x, o_out$y, gp = gpar(
    col = NA, fill = linearGradient(c("#2DD4BF", "#0EA5E9", "#6366F1"),
                                    x1 = 0, y1 = 1, x2 = 1, y2 = 0)))
  grid.polygon(o_mid$x, o_mid$y, gp = gpar(col = NA, fill = "#F8FAFC"))
  # background inside the octagon
  inner <- polygonGrob(o_in$x, o_in$y)
  pushViewport(viewport(clip = inner))
  bg <- "#0E2242"
  grid.rect(gp = gpar(col = NA, fill = bg))
  # PolSAR image at the bottom, blended into the background
  grid.raster(pauli, x = 0.5, y = 0.02, width = 0.98, height = 0.98 * 0.62,
              just = c("centre", "bottom"), interpolate = TRUE)
  grid.rect(x = 0.5, y = 0.02, width = 1, height = 0.62, just = c("centre", "bottom"),
            gp = gpar(col = NA, fill = linearGradient(
              paste0(bg, c("FF", "D9", "40", "00")),
              stops = c(0, 0.25, 0.7, 1), x1 = 0.5, y1 = 1, x2 = 0.5, y2 = 0)))
  # radar footprint lines
  for (k in 1:3) {
    grid.lines(x = c(0.665, 0.18 + 0.10 * k), y = c(0.735, 0.08),
               gp = gpar(col = "#5EEAD4", alpha = 0.18, lwd = 1.2))
  }
  # radar waves from the satellite
  th <- seq(200, 250, length.out = 60) * pi / 180
  for (k in 1:4) {
    r <- 0.06 + 0.055 * k
    grid.lines(0.665 + r * cos(th), 0.735 + r * sin(th),
               gp = gpar(col = "#5EEAD4", alpha = 1 - 0.18 * k, lwd = 7 - k, lineend = "round"))
  }
  # satellite
  pushViewport(viewport(x = 0.67, y = 0.74, width = 0.22, height = 0.22, angle = 35))
  grid.rect(0.5, 0.5, 0.22, 0.16, gp = gpar(fill = "#E2E8F0", col = "#0B1A33", lwd = 2))
  for (s in c(-1, 1)) {
    grid.rect(0.5 + s * 0.31, 0.5, 0.34, 0.13,
              gp = gpar(fill = linearGradient(c("#38BDF8", "#1D4ED8")), col = "#0B1A33", lwd = 2))
    for (q in 1:3) grid.lines(0.5 + s * (0.14 + 0.085 * q) * c(1, 1), c(0.435, 0.565),
                              gp = gpar(col = "#0B1A33", lwd = 1))
    grid.lines(0.5 + s * c(0.11, 0.14), c(0.5, 0.5), gp = gpar(col = "#0B1A33", lwd = 2))
  }
  grid.circle(0.5, 0.36, 0.045, gp = gpar(fill = "#F8FAFC", col = "#0B1A33", lwd = 2))
  popViewport()
  # name
  grid.text("SITSAR", x = 0.503, y = 0.497, gp = gpar(fontface = "bold", fontsize = 92,
            col = "#00000066", fontfamily = "sans"))
  grid.text("SITSAR", x = 0.5, y = 0.505, gp = gpar(fontface = "bold", fontsize = 92,
            col = "#FFFFFF", fontfamily = "sans"))
  grid.lines(c(0.25, 0.75), c(0.425, 0.425), gp = gpar(col = "#2DD4BF", lwd = 3))
  grid.text("SAR · PolSAR · Information Theory", x = 0.5, y = 0.385,
            gp = gpar(fontsize = 19, col = "#CBD5E1", fontfamily = "sans"))
  popViewport()
  # thin inner border
  grid.polygon(o_in$x, o_in$y, gp = gpar(col = "#0B1A33", lwd = 2, fill = NA))
}

ragg::agg_png("man/figures/logo.png", width = 600, height = 600, background = "transparent", res = 72)
draw_logo()
dev.off()
svglite::svglite("logo/logo.svg", width = 600 / 72, height = 600 / 72, bg = "transparent")
draw_logo()
dev.off()
