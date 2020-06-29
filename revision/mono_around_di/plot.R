#!/usr/bin/env Rscript

library("RColorBrewer")
col.pal <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))
r_data_file <- "R_cache_data.RData" 
trimean <- function(x) mean(quantile(x, c(0.25, 0.5, 0.5, 0.75)))

data_file <- list(
  "mono" = list("file"= "mono_di_te_by_mono_matrix.tsv",
                "peakpos" = "V26"
                ),
  "di" = list("file" = "mono_di_te_by_di_matrix.tsv",
              "peakpos" = "V66"
              )
)

if (file.exists(r_data_file)) {
  load(r_data_file)
  writeLines(paste("loaded data with timestamp:", psmd.timestamp))
} else {
  psmd <- list()
  quant <- list()
  groups <- list()
#  groups <- list()
  for (data_type in names(data_file)) {
    psmd[[data_type]] <- list()
    quant[[data_type]] <- list()
    groups[[data_type]] <- list()
    writeLines(paste("processing", data_type))
    dat <- read.table(data_file[[data_type]][["file"]])
    # sort_pos <- data_file[[data_type]][["peakpos"]]
    for (sort_pos in c("V26", "V66")) {
      writeLines(paste("I got", sort_pos, "as sorting column"))
      dat_sorted <- dat[order(dat[sort_pos]), 2:81]
      writeLines("data sorted")
      breaks <- unique(quantile(dat[[sort_pos]], seq(0, 1, 1/100)))
      writeLines(paste("I got", length(breaks), "many breaks"))
      p_groups <- cut(dat_sorted[[sort_pos]],
                      breaks = breaks, include.lowest = T)
      psmd[[data_type]][[sort_pos]] <- as.matrix(
        apply(dat_sorted, 2,
              function(x) aggregate(x ~ p_groups, FUN=trimean)[,2])
      )
      quant[[data_type]][[sort_pos]] <- as.data.frame(
        apply(dat_sorted, 2, quantile, c(0.25, 0.5, 0.75))
      )
      groups[[data_type]][[sort_pos]] <- as.data.frame(table(p_groups))
    }
  }
  psmd.timestamp <- Sys.time()
  save(psmd, psmd.timestamp, quant, groups, file = r_data_file)
  writeLines(paste("saved data with timestamp:", psmd.timestamp))
}

layout.matrix  <- matrix(c(0,  2,  2,  0,  5,  0,
                           16, 3,  1,  0,  4,  6,
                           0,  3,  15, 0,  17, 6,
                           0,  8,  8,  0,  11, 0,
                           18, 9,  7,  0,  10, 12,
                           0,  9,  13, 0,  14, 12),
                         nrow = 6, ncol = 6, byrow = TRUE)
layout.heights <- c(2, 6, 1, 2, 6, 1)
layout.widths <- c(0.5, 2, 10, 0.2, 10, 2)

## layout.matrix  <- matrix(c(2, 1, 4, 3, 5), nrow = 5, ncol = 1)
## layout.heights <- c(2, 6, 2, 6, 1)
## layout.widths <- 1

pdf("Fig_mono_di_peak_matrix.pdf", width = 21, height = 8)
layout(mat = layout.matrix, heights = layout.heights, widths = layout.widths)
for (dat_type in names(data_file)) {
  for (sort_pos in c("V26", "V66")) {

    # Plot 1: Heatmap of matrix
    par(mai = c(0, 0, 0, 0))
    image(t(psmd[[dat_type]][[sort_pos]]), col = col.pal(n = 40),
          breaks = quantile(c(psmd[[dat_type]][[sort_pos]]),
                            seq(0, 1, 1/40)), axes = F)

    # Plot 2: Quantiles of matrix (disome sorted)
    if (sort_pos == "V26") {
      par(mar = c(0, 13.3, 0, 0))
    } else {
      par(mai = c(0, 0, 0, 0))
    }
#plot(1:80, apply(quant[["di"]], 2, function(x) (x[1] + 2*x[2] + x[3])/4),
#     type="n", xaxs="i", axes=F)
#lines(1:40, apply(quant[["di"]][1:40], 2, function(x) (x[1] + 2*x[2] + x[3])/4))
#lines(41:80, apply(quant[["di"]][41:80], 2, function(x) (x[1] + 2*x[2] + x[3])/4))
    plot(1:80, quant[[dat_type]][[sort_pos]][3,], ylab = "density x 10-3",
         type="n", xaxs="i", axes = F,
         ylim = c(0, max(quant[[dat_type]][[sort_pos]][3, ])),
         xlim = c(0.5, 80.5))
    polygon(c(1:40, 40:1),
            c(quant[[dat_type]][[sort_pos]][1, 1:40],
              quant[[dat_type]][[sort_pos]][3, 40:1]),
            border = "gray", col = "gray")
    polygon(c(41:80, 80:41),
            c(quant[[dat_type]][[sort_pos]][1, 41:80],
              quant[[dat_type]][[sort_pos]][3, 80:41]),
            border = "gray", col = "gray")
    abline(v=seq(0, 80, 5), lty = 2, col = "lightgray")
    lines(1:40, quant[[dat_type]][[sort_pos]][2, 1:40])
    lines(41:80, quant[[dat_type]][[sort_pos]][2, 41:80])
    if (sort_pos == "V26") {
      at = axTicks(2)
      #      labels = formatC(at, format = "e", digits = 0)
      labels = at / 1e-03
      axis(2, at = at, labels = labels, las = 1)
    }

    # Plot 3: Histogram
    par(mar = c(2.625, 0, 0, 0))
    direction <- if (sort_pos == "V26") -1 else 1
    bardata <- groups[[dat_type]][[sort_pos]]$Freq * direction
    barplot(bardata, 1.038, space = 0, names = FALSE, horiz = TRUE, xaxt = 'n',
            ylim = c(0, length(bardata)), border = "grey30")
    at = axTicks(1)
    labels = abs(at) / 1e03
    axis(1, at = at, labels = labels)
  }
}
# Bottom x labels - positions
for (i in 1:2) {
  par(mar = c(3.3, 0.5, 0, 0.5))
  plot(1:80, rep(1, 80), xlim = c(1, 80), xaxs = "i", axes = F, type = "n")
  at = c(5, 10, 15, 20, 23, 24, 25, 30, 35, 40)
  labels = c("-20", "-15", "-10", "-5", "E", "P", "A", "+5", "+10", "+15")
  axis(1, at = at, labels = labels)
  axis(1, at = at + 40, labels = labels)
}
# sorted/pulled by monosome / disome labels
for (dtype in c("monosomes", "disomes")) {
  par(mai = c(0, 0, 0, 0))
  plot(1, 1, type = "n", axes = F)
  legend("center", paste("Sorted by", dtype), bty = "n")

  par(mai = c(0, 0, 0, 0))
  plot(1, 1, type = "n", axes = F)
  text(c(1, 1), paste("Pulled by", dtype), pos = 1, srt = 90)
}
dev.off()
      

