#!/usr/bin/env Rscript

source("rustr.R")

for (dtype in c("rp", "tr", "di")) {
  for (signalp in c("include", "exclude")) {
    dir_name <- paste(dtype, "aminoacid", "signalp", signalp, sep = "_")
    rust_file <- file.path(dir_name, "RUST_aminoacid_02-May-19_09.20")
    rust_data <- read.rust(rust_file)
    fig <- logo.rust(rust_data,
                     threshold = FALSE, use_expected = TRUE)
    fig <- fig + ylim(c(-6, 3.5))
    ggsave(paste("Fig", "aminoacid", dtype, "signalp", signalp, "q1_unfiltered__full_logo.pdf", sep = "_"), plot = fig)
  }
}
