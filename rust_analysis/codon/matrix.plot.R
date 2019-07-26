matrix.plot <- function(weights, type = "tri", col_weights = NULL,
                        pal_cols = c("lightblue", "blue"), pal_size = 255,
                        xlab = NULL, ylab = NULL) {
  area_from_height <- function(h, type) {
    if (type == "tri") {
      area <- h**2 / sqrt(3)
    }
    if (type == "circle") {
      area <- pi * (h / 2) ** 2
    }
    area
  }
  area_from_radius <- function(r, type) {
    if (type == "tri") {
      area <- area_from_height(3 * r / 2, type = type)
    }
    if (type == "circle") {
      area <- pi * r ** 2
    }
    area
  }
  radius_from_area <- function(area, type) {
    if (type == "tri") {
      r <- 2 * sqrt(area) / sqrt(3 * sqrt(3))
    }
    if (type == "circle") {
      r <- sqrt(area / pi)
    }
    r
  }
  draw_grid <- function(coords, fg) {
    symbols(coords, add = TRUE, inches = FALSE, fg = fg, bg = NA,
            rectangles = matrix(1, nrow = length(coords$x), ncol = 2))
  }
  draw_triangle <- function(x, y, radius, col = NA) {
    polygon(x = c(x, x - sqrt(3) * radius / 2, x + sqrt(3) * radius / 2),
            y = c(y - radius, y + radius / 2, y + radius / 2),
            border = NA, col = col)
  }
  draw_circle <- function(x, y, radius, angles, col = NA) {
    circ_x <- cos(angles) * radius + x
    circ_y <- sin(angles) * radius + y
    polygon(x = circ_x, y = circ_y, border = NA, col = col)
  }
  draw_circles <- function(coords, radii, col = NA) {
    symbols(coords, add = TRUE, inches = FALSE, circles = radii,
            fg = col, bg = col)
  }
  if (type == "tri") {
    # calculation of radii and direction of triangles
    areas <- abs(weights)
    norm_factor <- area_from_radius(0.5, type) / max(areas)
    norm_area <- areas * norm_factor
    radii <- sign(weights) * as.numeric(apply(norm_area, c(1,2),
                                              radius_from_area, type))
  }
  if (type == "circle") {
    stopifnot(all(weights >= 0))
    min_weight <- min(weights)
    max_weight <- max(weights)
    max_area <- area_from_radius(0.5, type)
    min_area <- area_from_radius(0.025, type)
    if (min_area < min_weight) {
      min_area <- 0
      min_weight <- 0
    }
    norm_slope <- (max_area - min_area) / (max_weight - min_weight)
    norm_intercept <- min_area - norm_slope * min_weight
    norm_area <- norm_intercept + weights * norm_slope
    radii <- as.numeric(apply(norm_area, c(1,2), radius_from_area, type))
  }
  # calculation of color gradient
  if (is.null(col_weights)) {
    col_weights <- weights
  }
  pal = colorRampPalette(pal_cols)
  cols = pal(pal_size)
  col_step = (pal_size - 1) / (max(col_weights) - min(col_weights))
  col_i = round((col_weights - min(col_weights)) * col_step) + 1
  leg_p = c(0, 0.25, 0.5, 0.75, 1)
  leg_steps = quantile(1:pal_size, probs = leg_p)
  leg_labels = round(min(col_weights) + (leg_steps - 1) / col_step, digits = 2)
  print(leg_labels)
  # calculation of coordinates and plotting
  coords <- xy.coords(data.frame(x = rep(1:nrow(weights), each = ncol(weights)),
                                 y = rep(1:ncol(weights), nrow(weights))),
                      xlab = xlab, ylab = ylab, log = "", recycle = FALSE)
  par(mar = c(1, 1, 2, 2))
  plot.new()
  plot.window(xlim = c(0.5, max(coords$x) + 0.5),
              ylim = c(max(coords$y) + 0.5, 0.5),
              log = "", asp = 1, xaxs = "i", yaxs = "i")
  if (type == "tri") {
    for (i in 1:length(coords$x)) {
      draw_triangle(coords$x[i], coords$y[i], radii[i], cols[col_i[i]])
    }
  }
  if (type == "circle") {
#    draw_circles(coords, radii, cols[col_i])
    increment <- 2 * pi / 100
    angles <- seq(0, 2 * pi - increment, by = increment)
    for (i in 1:length(coords$x)) {
      draw_circle(coords$x[i], coords$y[i], radii[i], angles, cols[col_i[i]])
    }
  }
  axis(3, 1:max(coords$x), colnames(weights), tick = FALSE)
  text(x = -0.25, y = 1:max(coords$y), labels = rownames(weights),
       srt = 0, adj = 0, xpd = TRUE)
  x1 <- max(coords$x) + 1
  x2 <- x1 + 1
  y1 <- 0
  y2 <- max(coords$y) + 1
  rasterImage(as.raster(matrix(cols, ncol = 1)), x1, y1, x2, y2)
  text(x = x2, y = y1 + 0.5 + (y2 - y1 - 1) * rev(leg_p),
       labels = leg_labels, srt = 0, adj = c(0, 0.5), xpd = TRUE)
  print(y1 + 0.5 + (y2 - y1 - 1) * rev(leg_p))
  draw_grid(coords, "lightgray")
}
