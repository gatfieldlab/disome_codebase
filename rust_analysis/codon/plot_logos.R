#!/usr/bin/env Rscript

source("rustr.R")

for (dtype in c("rp", "tr", "di")) {
  rust_file <- paste(dtype, "aminoacid/RUST_aminoacid_12-Jun-18_11.37", sep = "_")
  rust_data <- read.rust(rust_file)
  fig <- logo.rust(rust_data,
                   threshold = FALSE, use_expected = TRUE)
  fig <- fig + ylim(c(-6, 3.5))
  ggsave(paste("Fig", "aminoacid", dtype, "unfiltered__full_logo.pdf", sep = "_"), plot = fig)
}
