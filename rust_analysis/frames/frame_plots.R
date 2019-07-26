#!/usr/bin/env Rscript

source("../rustr.R")
require("ggplot2")
require("RColorBrewer")

mlist <- dir(pattern = "*meta.txt")
blist <- unlist(lapply(mlist, substr, 1, 26))
i <- 0
frame_data <- data.frame()
for (f in 0:2) {
    for (l in 55:64) {
       i <- i + 1
        tmp <- read.rust(blist[i])
        dfn <- paste("l", l, "_f", f, sep = "")
        print(paste(dfn, tmp$meta_data$profile$Offset,
        	        tmp$meta_data$profile$Frames, sep=" - "))
        frame_data <- rbind(frame_data, data.frame(length = l, frame = f,
        	                                       pos = -15:74, kl = tmp$kl))
    }
}


plot_base <- ggplot(data = frame_data,
	                aes(x = pos, y = kl, group = frame, color = factor(frame)))
colors <- rev(RColorBrewer::brewer.pal(n = 4, name = "BuPu"))[1:3]
fig <- plot_base +
       facet_grid(length ~ .) +
       geom_line() +
       scale_color_manual(values = colors) +
       geom_vline(xintercept = c(45, 48),
       	          linetype = "dashed", color = "gray") +
       theme_light()
ggsave("read_len_frame_KL_plots.pdf", plot = fig)