#!/usr/bin/env Rscript

require("viridis")
library("geneticcode")
library("RColorBrewer")
source("rustr.R")

pal <- viridis_pal(option = "A")
# pal <- function(x) brewer.pal(n = x, "Dark2")

plot_params <- list(
    "topology" = list(legend_pos = "left",
                      colors = c("#a6611a", "#018571"),
                      filter = FALSE, normalized = TRUE, log = "y"))
read_boundaries <- list(
    "di" = list(c(-16.5, -14.5), c(3.5, 4.5)),
    "rp" = list(c(-5.5, -4.5), c(3.5, 4.5)),
    "tr_all" = list(c(-6.5, -3.5), c(3.5, 6.5)))
x_shifts <- list(
  "topology" = 0)
plot_func <- "topology"
for (read_type in names(read_boundaries)) {
  rust_dir <- paste(read_type, plot_func, sep = "_")
  rust_files <- dir(path = rust_dir)
  rust_file <- rust_files[grep("meta", rust_files)]
  rust_base <- substr(rust_file, 1, nchar(rust_file) - 9)
  rust_data <- read.rust(paste(rust_dir, rust_base, sep = "/"))
  read_boundary <- read_boundaries[[read_type]]
  params <- plot_params[[plot_func]]
  file_dir <- "topology"
  file_name <- paste(read_type, "reduced", "raw",
                     "rust_norm.pdf", sep = "_")
  words <- c("Intramembrane", "Transmembrane")
  pdf(file.path(file_dir, file_name))
  colors <- if (params[["colors"]][1] == "auto") pal(length(words)) else params[["colors"]]
  plot.rust(rust_data, what = c("observed", "present"),words_subset = words,
            normalized = params[["normalized"]], col = colors,
            legend_pos = params[["legend_pos"]], log = params[["log"]],
            shift = x_shifts[[plot_func]],
            ap_site = TRUE, boxes = read_boundary, box_col = "gray95")
  dev.off()
}
quit(save = "no", status = 0)

