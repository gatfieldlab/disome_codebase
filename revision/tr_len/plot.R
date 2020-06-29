#!/usr/bin/env Rscript

r_data_file <- "R_cache_data.RData" 

data_file <- list(
  "all" = list(
    "file" = "di_matrix_cutoff_10_flag.tsv",
    "label" = "all_transcripts",
    "title" = "All transcripts"),
  "signalp" = list(
    "file" = "di_matrix_cutoff_10_flag_signalp.tsv",
    "label" = "signalp_transcripts",
    "title" = "Signal-peptide transcripts"))


if (file.exists(r_data_file)) {
  load(r_data_file)
  writeLines(paste("loaded data with timestamp:", psdm.timestamp))
} else {
  psdm <- list()
  for (data_type in names(data_file)) {
    psdm[[data_type]] <- read.table(
      data_file[[data_type]][["file"]],
      colClasses = c(
        "integer", "integer", "numeric", "numeric", "numeric", "numeric"),
      col.names = c("rank", "len", "start", "end", "norm", "density"))
  }
  psdm.timestamp <- Sys.time()
  save(psdm, psdm.timestamp, file = r_data_file)
  writeLines(paste("saved data with timestamp:", psdm.timestamp))
}

# Blue palette adapted from Color Brewer
pal.blue <- c("#c6dbef66", "#9ecae166", "#6baed666",
              "#4292c666", "#2171b566", "#08459466")
col.orange <- "#f9982066"
limits <- list("start" = c(-500, 6000),
               "end" = c(-6000, 500))
for (dat_type in names(data_file)) {
  for (sort_pos in names(limits)) {
    png(paste("Fig", data_file[[dat_type]][["label"]],
              sort_pos, "metagene_profile.png", sep = "_"),
        width = 800, height = 800)
    plot(psdm[[dat_type]][[sort_pos]], psdm[[dat_type]][["rank"]],
         col = col.orange, pch=16, cex=0.2, xlim=limits[[sort_pos]],
         xlab = paste("Position relative to", sort_pos, "(nt)"),
         ylab = "Transcript rank by CDS length",
         main = data_file[[dat_type]][["title"]])
    points(psdm[[dat_type]][[sort_pos]], psdm[[dat_type]][["rank"]],
    col = pal.blue[as.numeric(
      cut(psdm[[dat_type]][["density"]],
          breaks=c(0, 0.0008, 0.0027, 0.0086, 0.01, 0.05, 1)))],
    pch = 16, cex = 0.2)
    dev.off()
  }
}

      

