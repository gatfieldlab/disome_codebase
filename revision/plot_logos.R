#!/usr/bin/env Rscript

source("rustr.R")

rust_60 <- read.rust("di_aminoacid_len/RUST_aminoacid_27-Feb-20_15.13")
rust_63 <- read.rust("di_aminoacid_len/RUST_aminoacid_27-Feb-20_15.14")
fig60 <- logo.rust(rust_60, threshold = FALSE, use_expected = TRUE)
fig60 <- fig60 + ylim(c(-8, 4.5))
fig63 <- logo.rust(rust_63, threshold = FALSE, use_expected = TRUE)
fig63 <- fig63 + ylim(c(-8, 4.5))
ggsave(paste("Fig", "aminoacid", "di", "unfiltered_59-60_full_logo.pdf",
             sep = "_"), plot = fig60)
ggsave(paste("Fig", "aminoacid", "di", "unfiltered_62-63_full_logo.pdf",
             sep = "_"), plot = fig63)
