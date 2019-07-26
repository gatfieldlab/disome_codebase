#!/usr/bin/env Rscript

# Analysis code for disome project
#
VERSION <- "3.3.2"
args <- commandArgs(trailingOnly = TRUE)
run_env_name <- environmentName(environment())

# ====  Source in functions, settings and third-party code ====

source("functions.R")
source("settings.R")

# ==== Seet up versions ========================================================

code_version <- package_version(
  c(VERSION, SETTINGS_VERSION, FUNCTIONS_VERSION))
if (length(args) < 2) {
  run_version <- random_str()
  ver_msg <- "randomly generated"
} else {
  run_version <- args[2]
  ver_msg <- "user supplied"
}
logger(paste("Running environment:", run_env_name))
logger(paste("Using", ver_msg, "run version:", run_version))
logger(paste("Settings code version:", SETTINGS_VERSION))
logger(paste("Functions code version:", FUNCTIONS_VERSION))
logger(paste("Analysis code version", VERSION))
rm(ver_msg)

# ==== Set up data environment =================================================

data_env <- new.env()
prev_env <- new.env()
res_folder <- "."
if (length(args) > 0) {
  res_folder <- args[1]
  dir.create(file.path(res_folder), showWarnings = FALSE)
}
logger(paste("Using output folder:", res_folder))

# =======  Set up data file  ===================================================

if (length(args) > 2) {
  project_data_file <- args[3]
}
logger(paste("Using data file:", project_data_file))

# =======  Read previous results  ==============================================

logger("Loading previously saved data")
load_prev_data(registered_env = data_env, data_file = project_data_file)
attach(data_env)
copy_registrar(from_env = data_env, to_env = prev_env)

# =======  Local functions  ====================================================

log_n_register <- function(res_file, descr, lvl = 2, res_type = "Figure") {
  logger(paste(res_type, " '", res_file, "' was created", sep = ""),
         level = lvl)
  register_result_file(res_file, descr,
                       registered_env = data_env, run_ver = run_version)
}
get_res_file <- function(file_name) {
  dest_folder <- get("res_folder", envir = .GlobalEnv)
  file.path(dest_folder, file_name)
}

pass_dependencies <- function(module_name, dependencies) {
  dependency_met <- check_dependency(dependencies, data_env)
  if (! dependency_met) {
    stop(paste("[", module_name, "] Depends on:", dependencies))
  }
}

are_you_sure <- function() {
  if (interactive()) {
    answer <- readline("Are you sure ([Yy]es)? ")
  } else {
    cat("Are you sure ([Yy]es? ")
    answer <- readLines("stdin", n = 1)
  }
  substr(answer, 1, 1) %in% c("y", "Y")
}

# ======  Sample setup  ========================================================

if (do_module("Setup")) {

  logger("Initial setup")
  logger("creating samples", level = 2)
  for (sample in c("tr", "mono", "di")) {
    cur_samples <- read.table(
      paste(sample, get_option("Setup", "count_files_ext"), sep = ""),
      col.names = get_option("Setup", "count_cols"), stringsAsFactors = FALSE,
      header = FALSE, sep = "\t")
    cur_samples <- cur_samples[order(factor(cur_samples$treatment,
                                            levels = treatment_order),
                                     cur_samples$rep), ]
    cur_samples$treatment <- factor(cur_samples$treatment,
                                    levels = treatment_order)
    data_obj_name <- paste(sample, "samples", sep = "_")
    assign(data_obj_name, cur_samples)
    # ---- Data registration ----
    register_data_obj(data_obj_name, registered_env = data_env,
                      run_ver = run_version, by_name = T)
  }

  # ==== Gene annotations ====

  logger("creating annotation tables", level = 2)
  annot <- read.table(get_option("Setup", "cds_models"), header = FALSE,
                      sep = "\t", stringsAsFactors = FALSE,
                      col.names = get_option("Setup", "cds_cols"))
  composite <- annot[annot$tr_type == "composite", ]
  composite_sizes <- with(composite, cbind(x5utr = cds_start,
                                           cds = cds_end - cds_start,
                                           x3utr = size - cds_end, utr = NA,
                                           not_in_db = NA))
  row.names(composite_sizes) <- composite$gene_id
  annot$status <- factor(annot$status)
  annot$tr_type <- factor(annot$tr_type)
  gene_annotation <- read.table(get_option("Setup", "annotation"), header = T,
                                stringsAsFactors = F, sep = "\t",
                                col.names = get_option("Setup", "annot_cols"))
  single_cds_models <- names(which(table(annot$gene_id) == 2))
  # Size of raw count files (# genes)
  gene_ids_used <- read.table(get_option("Setup", "gids"))
  tsize <- nrow(gene_ids_used)
  # ---- Data registration ----
  register_data_obj(tsize, annot, composite_sizes,
                    single_cds_models, gene_annotation,
                    registered_env = data_env, run_ver = run_version)

  # ==== Experimental design ====

  logger("creating experimental design", level = 2)
  treatments <- factor(treatment_order, levels = treatment_order)
  master_design <- data.frame()
  for (dset in c("tr", "mono", "di")) {
    samples <- subset(get(paste(dset, "samples", sep = "_")), use == "Y")
    valid <- unique(samples[c("treatment", "rep")])
    design_df <- data.frame(
      treatment = valid$treatment,
      row.names = paste(valid$treatment, valid$rep, sep = "."))
    if (length(master_design) == 0) {
      master_design <- design_df
    } else {
      if (! identical(design_df, master_design)) {
        stop(paste("Experimental design in", dset, "doesn't match other",
                   "datasets' design"))
      }
    }
  }
  design_df <- master_design
  design_flat <- master_design$treatment
  design_reps <- factor(row.names(master_design))
  reps <- table(design_flat)
  design_exp <- formula(~ treatment)

  # ==== Read-length restrictions ====

  read_len_mono <- c(28, 32)
  read_len_di <- c(56, 64)

  # ---- Data registration ----
  register_data_obj(design_df, design_flat, design_reps, design_exp,
                    reps, treatments, read_len_mono, read_len_di,
                    registered_env = data_env, run_ver = run_version)
  # ---- Garbage collection ----
  rm(gene_ids_used, composite)
  register_module("Setup", registered_env = data_env, run_ver = run_version)
}

# ======  Mapping analysis  ====================================================

if (do_module("Mapping")) {
  pass_dependencies("Mapping", c("Setup"))
  logger("Reading mapping log data and generating plots")
  meta_map_data <- read.table(
    get_option("Mapping", "map_data_file"), sep = "\t",
    col.names = get_option("Mapping", "map_data_cols"))
  map_type_order <- c("mouse_mt-tRNA", "mouse_tRNA", "unmapped", "human_rRNA",
                      "mouse_rRNA", "mouse_cDNA", "genomic")
  map_labels <- c("mt-tRNA", "tRNA", "unmapped", "human rRNA", 
                  "rRNA", "cDNA", "genome")
  map_colors <- brewer.pal(length(map_labels), "BrBG")
  blank_theme <- theme_minimal() + theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 14, face = "bold"))
  for (dset in c("mono", "di", "tr")) {
    logger(paste("processing", dset, "mapping outcomes"), level = 2)
    samples <- subset(get(paste(dset, "samples", sep = "_")), use == "Y")
    cur_map <- merge(samples, meta_map_data, by.x = "name", by.y = "sample")
    cur_map$map_type <- factor(cur_map$map_type, levels = map_type_order,
                               ordered = TRUE)
    aggr_maptype <- aggregate(count ~ map_type, cur_map, sum)
    aggr_maptype <- cbind(aggr_maptype,
                      pc = aggr_maptype$count / sum(aggr_maptype$count) * 100)
    aggr_treat <- aggregate(count ~ treatment + rep + map_type, cur_map, sum)
    # Plotting single aggregate picharts
    fig_res <- get_res_file(
      paste("Fig_mapping_piechart_", dset, ".pdf", sep = ""))
    fig_base <- ggplot(
      aggr_maptype, aes(x = factor(1), y = count, fill = map_type,
      label = ifelse(pc >= 1.5, paste(round(pc), "%", sep = ""), "")))
    fig <- fig_base +
      geom_bar(stat = "identity", width = 1) +
      coord_polar("y", start = 0) +
      geom_text(size = 3, position = position_stack(vjust = 0.5)) +
      scale_fill_manual("Map types", values = map_colors,
                        labels = map_labels) +
      ggtitle(dset) + blank_theme
    ggsave(fig_res, width = 10, height = 7, plot = fig, useDingbats = FALSE)
    log_n_register(fig_res, paste(
      "Piechart of mapping outcome for", dset), lvl = 3)
    # Plotting sample level piecharts
    fig_res <- get_res_file(
      paste("Fig_mapping_ind_samples_", dset, ".pdf", sep = ""))
    fig_base <- ggplot(
      aggr_treat, aes(x = factor(1), y = count, fill = map_type))
    fig <- fig_base +
      geom_bar(stat = "identity", width = 1, position = position_fill()) +
      coord_polar("y", start = 0) +
      facet_grid(treatment ~ rep) +
      scale_fill_manual("Map types", values = map_colors,
                        labels = map_labels) +
      ggtitle(dset) + blank_theme
    ggsave(fig_res, width = 10, heigh = 10, plot = fig, useDingbats = FALSE)
    log_n_register(fig_res, paste(
      "Piechart of mapping outcome for individual", dset, "samples"), lvl = 3)
    # Plotting sequence run level piecharts
    fig_res <- get_res_file(
      paste("Fig_mapping_ind_runs_", dset, ".pdf", sep = ""))
    fig_base <- ggplot(
      data.frame(cur_map, run = substr(cur_map$name, 6, 10)),
      aes(x = factor(1), y = count, fill = map_type))
    fig <- fig_base +
      geom_bar(stat = "identity", width = 1, position = position_fill()) +
      coord_polar("y", start = 0) +
      facet_grid(treatment + rep ~ run) +
      scale_fill_manual("Map types", values = map_colors,
                        labels = map_labels) +
      ggtitle(dset) + blank_theme
    ggsave(fig_res, width = 10, heigh = 10, plot = fig, useDingbats = FALSE)
    log_n_register(fig_res, paste(
      "Piechart of mapping outcome for individual", dset, "seq-runs"), lvl = 3)
    # Writing summary tables
    table_res <- get_res_file(
      paste("Table_mapping_outcome_", dset, ".csv", sep = ""))
    write.csv(aggr_treat, file = table_res)
    log_n_register(table_res,
      paste("Table for mapping outcome of", dset), lvl = 3, res_type = "Table")
  }

  # ---- Data registration ----
  register_data_obj(meta_map_data,
                    registered_env = data_env, run_ver = run_version)

  # Garbage collection
  # rm(map_data, nohits, hitonly, di_aggr, di_sums, di_simple, rrna, trna,
  #    di_verysimple, mono_map_data, mono_simple, m_trna, m_genome,
  #    mono_verysimple, blank_theme, fig_base, fig_res, fig)
  register_module("Mapping", registered_env = data_env, run_ver = run_version)
}

# ==== Insert size analysis ====================================================

if (do_module("Insert_size")) {
  pass_dependencies("Insert_size", c("Setup"))
  logger("Reading insert size data and generating plots")
  for (dset in c("mono", "di")) {
    insert_size <- data.frame()
    samples <- subset(get(paste(dset, "samples", sep = "_")), use == "Y")
    for (i in 1:nrow(samples)) {
      sample <- samples$treatment[i]
      rep <- samples$rep[i]
      file_base <- samples$name[i]
      insert_file <- paste("../insert_lengths/", file_base,
                           "_insert_sizes.tsv", sep = "")
      insert_len <- cbind(sample = sample, rep = rep,
                          read.table(insert_file,
                                     col.names = c("length", "count")))
      insert_size <- rbind(insert_size, insert_len)
      assign(paste(dset, "insert_size", sep = "_"), insert_size)
    }
  }
  # Meta analysis
  meta_insert_size <- rbind(cbind(type = "Di",
                                  aggregate(count ~ length,
                                            di_insert_size, sum)),
                            cbind(type = "Mono",
                                  aggregate(count ~ length,
                                            mono_insert_size, sum)))
  res_file <- get_res_file("Table_read_size_counts.tsv")
  write.table(meta_insert_size, file = res_file, quote = FALSE,
              col.names = TRUE, row.names = FALSE, sep = "\t")
  log_n_register(res_file, "Table of read counts by read size")
  # Cumulative analysis
  di_cum_insert_size <- subset(meta_insert_size, type == "Di")
  di_cum_insert_size$cumsum <- cumsum(di_cum_insert_size$count)
  di_cum_insert_size$cumpr <- di_cum_insert_size$cumsum /
                              sum(di_cum_insert_size$count)
  # ---- Plots ----
  # Di vs monosome
  res_fig <- get_res_file("Fig_mapped_insert_len_distribution.pdf")
  figbase <- ggplot(meta_insert_size, aes(x = length, y = count, fill = type))
  fig <- figbase +
          geom_bar(stat = "identity") +
          scale_x_continuous("Mapped insert length (nt)", limits = c(20, 79)) +
          scale_fill_manual(values = c(di_col, mono_col),
                            labels = c("Disome", "Monosome")) +
          theme_grey() +
          ylab("Count")
  ggsave(res_fig, width = 4, height = 4, plot = fig, useDingbats = F)
  log_n_register(res_fig,
    "Barplot of mapped read length distribution - mono vs disome")

  # Disome cumulative
  res_fig <- get_res_file("Fig_mapped_insert_len_cum_dist_disome.pdf")
  pdf(res_fig, width = 6, height = 7)
  with(di_cum_insert_size,
    plot(length, cumpr, pch = 16, type = "b",
         xlab = "Mapped insert length (nt)", ylab = "Cumulative percentage"))
  plot_min <- par("usr")[3]
  y_data <- subset(di_cum_insert_size, length %in% read_len_di)$cumpr
  segments(c(56, 64), c(plot_min, plot_min),
           c(56, 64), y_data, lty = 2, col = "grey")
  dev.off()
  log_n_register(res_fig,
    "Cumulative distribution of mapped insert length for disomes")

  # ---- Data registration ----
  register_data_obj(meta_insert_size, di_cum_insert_size,
                    registered_env = data_env, run_ver = run_version)
  # Garbage collection
  rm(insert_size, sample, rep, file_base, insert_file, insert_len, res_fig,
     figbase, fig, y_data, plot_min)
  register_module("Insert_size", registered_env = data_env,
                  run_ver = run_version)
}

# ==== Read in the counts ======================================================

if (do_module("Count")) {
  pass_dependencies("Count", c("Setup"))
  generic_count <- function(path = NULL, label = "", tag = "") {
    if (is.null(path)) {
      label <- "original"
      tag <- ""
    }
    logger(paste("Reading in the", label, "raw count data"))
    for (dset in c("mono", "di", "tr")) {
      logger(paste("processing", dset, "tables"), level = 2)
      # Read counts
      samples <- subset(get(paste(dset, "samples", sep = "_")), use == "Y")
      if (!is.null(path)) {
        samples$path <- path
      }
      counts <- read_counts(samples, "all", tsize)
      # Extract region info
      regions <- unlist(
        lapply(strsplit(colnames(counts), split = ".", fixed = T),
               function(x) paste(x[2:length(x)], collapse = "_")))
      # Merge seperate runs of same sample
      ttr <- t(counts)
      ttr_agr <- aggregate(
        x = ttr, by = list(rep(samples$rep, each = 5),
                           rep(samples$treatment, each = 5), regions), sum)
      ttr_names <- apply(ttr_agr, 1, function(x) paste(x[1:3], collapse = "_"))
      counts <- t(ttr_agr[, 4:ncol(ttr_agr)])
      colnames(counts) <- ttr_names
      if (tag != "") {
        count_obj_name <- paste(dset, tag, sep = "_")
      } else {
        count_obj_name <- dset
      }
      assign(count_obj_name, counts, envir = parent.frame())
      # Split data according to regions
      for (region in c("cds", "x5utr", "x3utr")) {
        data_obj_name <- paste(count_obj_name, region, sep = "_")
        assign(data_obj_name, grep(region, colnames(counts)),
               envir = parent.frame())
        register_data_obj(data_obj_name, registered_env = data_env,
                          run_ver = run_version, by_name = T)
      }
      d_5utr <- get(paste(count_obj_name, "x5utr", sep = "_"))
      d_cds <- get(paste(count_obj_name, "cds", sep = "_"))
      d_3utr <- get(paste(count_obj_name, "x3utr", sep = "_"))
      gene_wide_count <- counts[, d_5utr] + counts[, d_cds] + counts[, d_3utr]
      colnames(gene_wide_count) <- design_reps
      gene_count_name <- paste(count_obj_name, "gene", "count", sep = "_")
      assign(gene_count_name, gene_wide_count, envir = parent.frame())
      # ---- Data registration ----
      register_data_obj(count_obj_name, gene_count_name,
                        registered_env = data_env, run_ver = run_version,
                        by_name = T)
    }
    # Some garbage collection
    rm(ttr, counts, samples, ttr_agr, regions, gene_wide_count,
       gene_count_name, count_obj_name)
    list(label = label, tag = tag, path = path)
  }
  if (!exists("count_tables")) {
    count_tables <- list()
  }
  tables_to_count <- get_option("Count", "tables")
  for (table_props in tables_to_count) {
    count_attrs <- generic_count(
      path = table_props[["path"]],
      label = table_props[["label"]],
      tag = table_props[["tag"]])
    count_tables[[count_attrs[["label"]]]] <- count_attrs
  }
  register_data_obj(count_tables,
                    registered_env = data_env, run_ver = run_version)
  register_module("Count", registered_env = data_env, run_ver = run_version)
}

# ==== Filtering ===============================================================

if (do_module("Filter")) {
  pass_dependencies("Filter", c("Setup", "Count"))
  logger("Generating basic gene filters")
  for (label in names(count_tables)) {
    tag <- count_tables[[label]][["tag"]]
    idex <- if (tag != "") paste("_", tag, sep = "") else ""
    logger(paste("for", label, "counts"), level = 2)
    # Low-level filtering based on raw counts
    logger("filtering gene lists on raw counts", level = 3)
    for (dset in c("mono", "di", "tr")) {
      cset <- if (tag != "") paste(dset, "_", tag, sep = "") else dset
      cur_gene <- get(paste(cset, "gene", "count", sep = "_"))
      cur_cds <- get(cset)[, get(paste(cset, "cds", sep = "_"))]
      # filters based on CDS counts
      cur_filter <- get_count_filter(cur_cds, min.fraction = 1 / 3)
      cur_filter_ids <- names(cur_filter)
      # filters based on gene-wide counts
      cur_filter_gene <- get_count_filter(
        cur_gene, min.count = 30, min.fraction = 1 / 2)
      cur_filter_gene_ids <- names(cur_filter_gene)
      # filters excluding non-DB CDS counts
      cur_filter_db_ids <- intersect(cur_filter_ids, row.names(composite_sizes))
      # ---- Data registration ----
      obj_names <- c(paste(cset, "filter", sep = "_"),
                     paste(cset, "filter_ids", sep = "_"),
                     paste(cset, "filter_gene", sep = "_"),
                     paste(cset, "filter_gene_ids", sep = "_"),
                     paste(cset, "filter_db_ids", sep = "_"))
      obj_sources <- c("cur_filter", "cur_filter_ids", "cur_filter_gene",
                       "cur_filter_gene_ids", "cur_filter_db_ids")
      len_obj <- length(obj_names)
      if (len_obj !=  length(obj_sources))
        stop("Number of assignees and data objects do not match in [Filter]")
      for (oid in 1:len_obj) {
        assign(obj_names[oid], get(obj_sources[oid]))
        register_data_obj(obj_names[oid], registered_env = data_env,
                          run_ver = run_version, by_name = T)
      }
    }
    logger("filtering gene lists across seq types", level = 3)
    expressed_ids <- intersect(
      di_filter_ids, intersect(mono_filter_ids, tr_filter_ids))
    expressed_single_cds <- intersect(single_cds_models, expressed_ids)
    expressed_gene_ids <- intersect(
      di_filter_gene_ids, intersect(mono_filter_gene_ids, tr_filter_gene_ids))
    # ---- Data registration ----
    obj_names <- c(paste("expressed_ids", idex, sep = ""),
                   paste("expressed_single_cds", idex, sep = ""),
                   paste("expressed_gene_ids", idex, sep = ""))
    obj_sources <- c("expressed_ids", "expressed_single_cds",
                     "expressed_gene_ids")
    len_obj <- length(obj_names)
    if (len_obj !=  length(obj_sources))
      stop("Number of assignees and data objects do not match in [Filter]")
    for (oid in 1:len_obj) {
      assign(obj_names[oid], get(obj_sources[oid]))
      register_data_obj(obj_names[oid], registered_env = data_env,
                        run_ver = run_version, by_name = T)
    }
  }
  register_module("Filter", registered_env = data_env, run_ver = run_version)
}

# ==== Normalization of count data =============================================

if (do_module("Normalize")) {
  pass_dependencies("Normalize", c("Setup", "Count", "Filter"))
  logger("Data normalization")

  get_treatment_means <- function(dataset) {
    colnames(dataset) <- design_flat
    dataset <- melt(dataset, value.name = "count")
    dataset <- dcast(aggregate(count ~ Var2 + Var1, dataset, mean),
                     Var1 ~ Var2, value.var = "count")
    rownames(dataset) <- dataset[, 1]
    dataset[, 2:ncol(dataset)]
  }

  count_labels <- get_option("Normalize", "count-labels")
  if (is.null(count_labels)) count_labels <- names(count_tables)
  for (label in count_labels) {
    tag <- count_tables[[label]][["tag"]]
    idex <- if (tag != "") paste("_", tag, sep = "") else ""
    logger(paste("for", label), level = 2)
    # Quartile based library size normalization
    logger("quartile based library size normalization", level = 3)
    sample_norm_df <- data.frame()
    for (dset in c("mono", "di", "tr")) {
      logger(paste("normalizing", dset, "count data-sets"), level = 4)
      cset <- if (tag != "") paste(dset, "_", tag, sep = "") else dset
      cur_cds <- get(cset)[, get(paste(dset, "cds", sep = "_"))]
      cur_filter <- get(paste(cset, "filter", sep = "_"))
      cur_filter_ids <- get(paste(cset, "filter_ids", sep = "_"))
      cur_filter_db_ids <- get(paste(cset, "filter_db_ids", sep = "_"))
      norm_factors <- get_norm_factors(cur_cds[cur_filter, ], design_flat)
      if (is_option_set("Normalize", "print-density-db")) {
        dtype <- paste(toupper(substring(dset, 1, 1)), substring(dset, 2),
                       sep = "", collapse = " ")
        samples <- subset(get(paste(dset, "samples", sep = "_")), use == "Y")
        samples$path <- get_option("Normalize", "density_folder")
        samples$ext <- get_option("Normalize", "density_extension")
        samples$treatment <- paste(dtype, samples$treatment, sep = "")
        design_index <- as.numeric(factor(samples$sample,
                                          levels = unique(samples$sample)))
        sample_norm_df <- rbind(sample_norm_df,
                                cbind(samples, norm_factors[design_index]))
      }
      norm_cds <- t(apply(cur_cds[cur_filter, ], 1, "/", norm_factors))
      colnames(norm_cds) <- design_reps
      norm_cds_mean <- get_treatment_means(norm_cds)
      logger(paste("calculating RPKM values for", dset), level = 4)
      libsize <- get_norm_factors(cur_cds[cur_filter, ],
                                  design_flat, "efflibsize")
      rpkm_cds <- norm_cds[cur_filter_db_ids, ] /
        composite_sizes[cur_filter_db_ids, "cds"] / libsize * 1e+09
      rpkm_cds_mean <- get_treatment_means(rpkm_cds)
      # ---- Data registration ----
      obj_names <- c(paste("norm_factors", cset, sep = "_"),
                     paste("norm_cds", cset, sep = "_"),
                     paste("norm_cds", cset, "mean", sep = "_"),
                     paste("libsize", cset, sep = "_"),
                     paste("rpkm_cds", cset, sep = "_"),
                     paste("rpkm_cds", cset, "mean", sep = "_"))
      obj_sources <- c("norm_factors", "norm_cds", "norm_cds_mean", "libsize",
                       "rpkm_cds", "rpkm_cds_mean")
      len_obj <- length(obj_names)
      if (len_obj !=  length(obj_sources))
        stop("Number of assignees and data objects do not match in [Normalize]")
      for (oid in 1:len_obj) {
        assign(obj_names[oid], get(obj_sources[oid]))
        register_data_obj(obj_names[oid], registered_env = data_env,
                          run_ver = run_version, by_name = T)
      }
      # Garbage collection
      rm(list = obj_sources)
      rm(obj_names, obj_sources)
    }
    if (is_option_set("Normalize", "print-density-db")) {
      filename <- paste("density_files_db_autogen", idex, ".tsv", sep = "")
      write.table(sample_norm_df, file = filename,
                  sep = "\t", quote = F, col.names = F, row.names = F)
      rm(sample_norm_df, samples, design_index)
    }

    # Ribosome occupancy calculations

    logger("calculating normalized ribosome occupancies", level = 3)
    for (dset in c("mono", "di")) {
      logger(paste("based on", dset), level = 4)
      cset <- if (tag != "") paste(dset, "_", tag, sep = "") else dset
      rpkm_cds <- get(paste("rpkm_cds", cset, sep = "_"))
      expressed_ids <- get(paste("expressed_ids", idex, sep = ""))
      eff <- rpkm_cds[expressed_ids, ] / rpkm_cds_tr[expressed_ids, ]
      eff_noinf <- as.data.frame(no.inf(eff))
      eff_log <- log2(eff)
      eff_log_noinf <- no.inf(eff_log)
      colnames(eff_log_noinf) <- design_reps
      eff_log_noinf_mean <- get_treatment_means(eff_log_noinf)
      # ---- Data registration ----
      obj_names <- c(paste("eff", cset, sep = "_"),
                     paste("eff", cset, "noinf", sep = "_"),
                     paste("eff", cset, "log", sep = "_"),
                     paste("eff", cset, "log_noinf", sep = "_"),
                     paste("eff", cset, "log_noinf_mean", sep = "_"))
      obj_sources <- c("eff", "eff_noinf", "eff_log",
                       "eff_log_noinf", "eff_log_noinf_mean")
      len_obj <- length(obj_names)
      if (len_obj !=  length(obj_sources))
        stop("Number of assignees and data objects do not match in [Normalize]")
      for (oid in 1:len_obj) {
        assign(obj_names[oid], get(obj_sources[oid]))
        register_data_obj(obj_names[oid], registered_env = data_env,
                          run_ver = run_version, by_name = T)
      }
    }
    # Garbage collection
    rm(list = obj_sources)
    rm(rpkm_cds, obj_names, obj_sources)
  }
  register_module("Normalize", registered_env = data_env,
                  run_ver = run_version)
}

# ==== Signal peptide count distribution =======================================

if (do_module("Signalp_distribution")) {
  pass_dependencies("Signalp_distribution", c("Setup", "Count", "Gene_lists"))
  logger("Performing analysis of distribution of reads to signal region")
  robust <- gene_list[["robust_expressed_single_ids"]]
  robuste <- names(which(composite_sizes[robust, "cds"] > 75 * 3 * 2))
  rsizes <- composite_sizes[robuste, "cds"] - 75 * 3
  sig_dist <- data.frame()
  for (dset in c("mono", "di", "tr")) {
    logger(paste("analysing", dset), level = 2)
    cds <- get(paste(dset, "cds", sep = "_"))
    ori_count <- get(dset)[robuste, cds]
    f75_count <- get(paste(dset, "f75", sep = "_"))[robuste, cds]
    r75_count <- ori_count - f75_count
    norm_f75 <- f75_count / (75 * 3) / get(
      paste("norm_factors", dset, sep = "_"))
    norm_r75 <- r75_count / rsizes / get(
      paste("norm_factors", dset, sep = "_"))
    sum_f75 <- rowSums(norm_f75)
    sum_r75 <- rowSums(norm_r75)
    tempdf <- data.frame(read_type = dset,
                         signalp = robuste %in% gene_list[["signalp_ids"]])
    sig_dist <- rbind(sig_dist,
                      cbind(tempdf, region = "first_75",
                            proportion = sum_f75 / (sum_f75 + sum_r75)),
                      cbind(tempdf, region = "rest",
                            proportion = sum_r75 / (sum_f75 + sum_r75)))
  }
  sig_dist$read_type <- factor(sig_dist$read_type,
                               levels = c("di", "mono", "tr"))
  res_fig <- get_res_file("Fig_distribution_signalP_first75_rest.pdf")
  signalp_labeller <- c("TRUE" = "coding for signal peptide",
                        "FALSE" = "not coding for signal peptide")
  readtype_labeller <- c("di" = "Disome", "mono" = "Monosome")
  fig <- ggplot(data = subset(sig_dist, read_type != "tr"),
                aes(y=proportion, x = region, fill = signalp)) +
    geom_violin() +
    facet_grid(read_type ~ signalp,
               labeller = labeller(read_type = readtype_labeller,
                                   signalp = signalp_labeller)) +
    stat_summary(fun.y="median", geom="point", size=3, show.legend = FALSE) +
    theme_minimal() +
    scale_fill_manual(values = c("gray", "brown2"), labels = c("No", "Yes"),
                      name = "Coding for\nsignal peptide?") +
    scale_x_discrete(labels = c("1 to 75 codon", "76 - end"))
  ggsave(res_fig, plot = fig)
  log_n_register(res_fig, paste(
    "Violin plot for distribution of reads to the first 75 codons & the rest"))
  register_data_obj(sig_dist, registered_env = data_env,
                    run_ver = run_version)
  register_module("Signalp_distribution", registered_env = data_env,
                  run_ver = run_version)
}


# ==== Differential expression analysis ========================================

if (do_module("Diff_exp")) {
  pass_dependencies("Diff_exp", c("Setup", "Count", "Filter", "Normalize"))
  logger("Performing analysis of differential expression using DESeq2")
  for (dset in c("mono", "di", "tr")) {
    logger(paste("analysing", dset), level = 2)
    logger("creating DESeq2 count and dispersion tables", level = 3)
    count <- get(dset)[get(paste(dset, "filter", sep = "_")),
                       get(paste(dset, "cds", sep = "_"))]
    colnames(count) <- design_reps
    cd <- get_cds(count, colData = design_df, design = design_exp,
                  get(paste("norm_factors", dset, sep = "_")))
    theta <- get_dispersions(cd)
    logger("performing negative binomial GLM with Wald statistics", level = 3)
    cd <- nbinomWaldTest(cd)
    res <- results(cd, contrast = c("treatment", "ZT12", "ZT00"))
    res <- res[order(res$padj), ]
    siglist <- rownames(res)[which(res$padj < 0.05)]
    # ---- Data registration ----
    obj_names <- c(paste("cd", dset, sep = "_"),
                   paste("theta", dset, sep = "_"),
                   paste("res", dset, sep = "_"),
                   paste("siglist", dset, sep = "_"))
    obj_sources <- c("cd", "theta", "res", "siglist")
    len_obj <- length(obj_names)
    if (len_obj !=  length(obj_sources))
      stop("Number of assignees and data objects do not match in [Diff_exp]")
    for (oid in 1:len_obj) {
      assign(obj_names[oid], get(obj_sources[oid]))
      register_data_obj(obj_names[oid], registered_env = data_env,
                        run_ver = run_version, by_name = T)
    }
  }
  # Garbage collection
  rm(list = obj_sources)
  rm(obj_names, obj_sources)
  register_module("Diff_exp", registered_env = data_env, run_ver = run_version)

}

# ==== Principal component analysis ============================================

if (do_module("PCA")) {
  pass_dependencies("PCA", c("Diff_exp"))
  logger("Performing principal component analysis (PCA)")
  pca_cols <- brewer.pal(n = 3, name = "Dark2")
  for (dset in c("mono", "di", "tr")) {
    logger(paste("analysing", dset), level = 2)
    cd <- get(paste("cd", dset, sep = "_"))
    vsd <- varianceStabilizingTransformation(cd, blind = TRUE)
    pca <- getPCA(vsd, intgroup = c("treatment"), ntop = 400)
    pca_p <- list(x = as.numeric(pca$x[, 1]), y = as.numeric(pca$x[, 2]))
    pca_labels <- as.character(unlist(lapply(reps, function(x) 1:x)))
    res_fig <- get_res_file(paste("Fig", dset, "PCA.pdf", sep = "_"))
    pdf(res_fig, useDingbats = FALSE)
    eqscplot(pca_p, col = "white", xlab = "PC.1", ylab = "PC.2", main = dset)
    pretty_grid(5, col = "grey")
    points(pca_p, col = rep(pca_cols, reps),
           pch = 19, cex = 1.3)
    text(pca_p, pca_labels, pos = 1)
    legend("bottomleft", legend = as.character(treatments),
           col = pca_cols, pch = 16, pt.cex = 1.3, ncol = 1,
           box.lwd = 0.5, inset = c(0.01, 0.02))
    dev.off()
    log_n_register(res_fig, paste("PCA for", dset))
  }
  # Garbage collection
  rm(cd, vsd, pca_p, pca_labels, pca_cols, res_fig)
  register_module("PCA", registered_env = data_env, run_ver = run_version)

}

# ==== TE change analysis ======================================================

if (do_module("Diff_TE")) {
  pass_dependencies("Diff_TE", c("Setup", "Count", "Filter", "Normalize"))
  logger("Performing differential TE analysis")
  if (is_option_set("Diff_TE", "xtail")) {
    require(xtail)
    xtail_bins <- get_option("Diff_TE", "xtail_bins")
    # This assumes that first condition in treatments is the Control
    ctrl <- as.character(treatments[1])
    for (treatment in treatments[-1]) {
      comparison <- c(ctrl, treatment)
      select <- design_flat %in% comparison
      mrna <- norm_cds_tr[expressed_ids, tr_cds][, select]
      condition <- as.character(design_flat[select])
      for (dset in c("mono", "di")) {
        logger(paste("comparing", treatment, "to Ctrl in", dset), level = 2)
        rpf <- get(paste("norm_cds", dset, sep = "_"))
        rpf <- rpf[expressed_ids, get(paste(dset, "cds", sep = "_"))][, select]
        res <- xtail(mrna = mrna, rpf = rpf, condition = condition,
                     normalize = FALSE, bins = xtail_bins)
        res_table <- resultsTable(res, sort.by = "pvalue.adjust",
                                  log2FCs = TRUE, log2Rs = TRUE)
        res_file <- get_res_file(paste(
          "Table_xtail", dset, treatment, "vs", ctrl, "cds.tsv", sep = "_"))
        write.table(res_table, res_file,
                    quote = FALSE, sep = "\t", col.names = NA)
        log_n_register(res_file,
          paste("xtail results table for", treatment,
                "vs", ctl, "using", dset), lvl = 3, res_type = "Table")
        # ---- Data registration ----
        obj_names <- c(paste("xtail_res", dset, treatment, sep = "_"),
                       paste("xtail_table", dset, treatment, sep = "_"))
        obj_sources <- c("res", "res_table")
        for (oid in 1:length(obj_names)) {
          assign(obj_names[oid], get(obj_sources[oid]))
          register_data_obj(obj_names[oid], registered_env = data_env,
                            run_ver = run_version, by_name = T)
        }
      }
    }
    # Garbage collection
    rm(list = obj_sources)
    rm(treatment, comparison, select, condition, dset, mrna, rpf, res_file,
       obj_sources, obj_names)
  }
  if (is_option_set("Diff_TE", "riborex")) {
    # require(riborex)
  }
  register_module("Diff_TE", registered_env = data_env, run_ver = run_version)

}

# ==== Summary tables for general counts etc ===================================

if (do_module("Summary")) {

  logger("Preparing and printing summary tables")
  if (check_dependency(c("Setup", "Count", "Filter"), data_env)) {
    count_summary <- data.frame(
      aggregate(sample ~ rep + treatment,
                subset(tr_samples, use == "Y"), unique),
      tr_cds_sums = colSums(tr[, tr_cds]),
      rp_cds_sums = colSums(mono[, mono_cds]),
      di_cds_sums = colSums(di[, di_cds]),
      tr_sums = colSums(tr_gene_count),
      rp_sums = colSums(mono_gene_count),
      di_sums = colSums(di_gene_count),
      tr_usable = colSums(tr[tr_filter, tr_cds]),
      rp_usable = colSums(mono[mono_filter, mono_cds]),
      di_usable = colSums(di[di_filter, di_cds]))
    res_file <- get_res_file("Table_count_summary.tsv")
    write.table(count_summary, res_file, sep = "\t",
                col.names = NA, quote = FALSE)
    log_n_register(res_file, "Summary table for counts",
                   lvl = 2, res_type = "Table")
    # Garbage collection
    rm(count_summary, res_file)
  }
  if (check_dependency(c("Diff_exp"), data_env)) {
    count_summary <- NA
    # Garbage collection
    rm(count_summary)
  }
  register_module("Summary", registered_env = data_env, run_ver = run_version)

}

# ==== Custom gene lists =======================================================

if (do_module("Gene_lists")) {
  pass_dependencies("Gene_lists", c("Setup", "Filter"))
  make_list <- function(gids, name) {
    g_list <- list()
    g_name <- paste(name, "ids", sep = "_")
    t_name <- paste(name, "trs", sep = "_")
    g_list[[g_name]] <- intersect(unique(gids), expressed_ids)
    g_list[[t_name]] <- subset(annot, gene_id %in% g_list[[g_name]] &
                                      tr_type != "composite", select = tr_id)
    g_list
  }
  logger("Creating custom gene lists")
  gene_list <- list()
  data_folder <- get_option("Gene_lists", "folder")
  # Signal peptide predictions
  signalp <- read.table(
    file.path(data_folder, "signal_peptide_prots2.txt"),
    sep = "\t", stringsAsFactors = FALSE)
  signalp <- signalp[-4]
  signalp <- signalp[-1, ]
  colnames(signalp) <- c("gene_id", "tr_id", "prot_id")
  gene_list <- append(gene_list, make_list(gids = signalp$gene_id,
                                           name = "signalp"))
 # 7-transmembrane proteins
  tm7 <- read.table(
    file.path(data_folder, "mouse_7mem_gids.txt"),
    stringsAsFactors = FALSE, col.names = c("gene_id"))
  gene_list <- append(gene_list, make_list(gids = tm7$gene_id,
                                           name = "7tm"))
  # Ribosomal proteins
  for (rp_type in c("Mrpl", "Mrps", "Rpl", "Rps")) {
    gene_list <- append(gene_list, make_list(
      gids = subset(gene_annotation, grepl(rp_type, gene_name),
                    gene_id)$gene_id,
      name = rp_type))
  }
  # robustly expressed single isoform coding genes
  robust_expressed_single <- expressed_single_cds[
    which(rowMeans(rpkm_cds_tr_mean[expressed_single_cds, ]) > 5)]
  gene_list <- append(gene_list, make_list(
    gids = robust_expressed_single, name = "robust_expressed_single"))

  # Phylop scores.. low ..medium .. high
  phylop_table <- read.table(
    paste(data_folder, "tr_cds_mean_phylop_scores.v2.tsv", sep = "/"),
    sep = "\t", stringsAsFactor = FALSE, header = FALSE,
    col.names = c("gene_id", "tr_id", "con_score"))
  phylop_genes <- subset(phylop_table, gene_id %in% single_cds_models &
                                         tr_id %in% annot$tr_id)
  phylop_lims <- quantile(phylop_genes$con_score,
                          probs = c(0, 0.25, 0.5, 0.75, 1))
  phylop_groups <- cut(phylop_genes$con_score, phylop_lims,
                       labels = c("Q1", "Q2", "Q3", "Q4"))
  for (level in levels(phylop_groups)) {
    gene_list <- append(gene_list, make_list(
      gids = phylop_genes[phylop_groups == level, ]$gene_id,
      name = paste("phylop", level, sep = "_")))
  }

  # Export
  if (get_option("Gene_lists", "export")) {
    for (list_name in names(gene_list)) {
      res_file <- get_res_file(paste(list_name, "txt", sep ="."))
      write.table(gene_list[[list_name]],
                  file = res_file,
                  col.names = FALSE, row.names = FALSE, quote = FALSE)
      log_n_register(res_file, paste("Gene list", list_name),
                     lvl = 2, res_type = "Table")
    }
  }
  # Entrez ID lists
  if (get_option("Gene_lists", "entrez")) {
    expressed_entrez <- get_id_from_biomart(
      gid = expressed_ids, what_to_get = c("entrezgene", "external_gene_name"),
      version = 91)
  } else {
    expressed_entrez <- c()
  }

  # ---- Data registration ----
  register_data_obj(signalp, expressed_entrez, gene_list, phylop_table,
                    registered_env = data_env, run_ver = run_version)
  # garbage collection
  rm(robust_expressed_single)
  register_module("Gene_lists", registered_env = data_env,
                  run_ver = run_version)

}

# ==== Percentages along features ==============================================

if (do_module("Percent")) {
  pass_dependencies("Percent", c("Setup", "Count", "Filter"))
  logger("Calculating read percentages per transcript / per region")
  for (dset in c("mono", "di", "tr")) {
    cdf <- get(dset)
    cdf_x3utr <- get(paste(dset, "x3utr", sep = "_"))
    cdf_cds <- get(paste(dset, "cds", sep = "_"))
    cdf_x5utr <- get(paste(dset, "x5utr", sep = "_"))
    curtype_pc <- NULL
    for (cursam in 1: length(cdf_cds)) {
        cur_cdf <- cdf[, c(cdf_x5utr[cursam], cdf_cds[cursam],
                           cdf_x3utr[cursam])]
        curtype_pc <- cbind(curtype_pc, cur_cdf / rowSums(cur_cdf))
    }
    pc_cols <- ncol(curtype_pc)
    tpc <- data.frame(row.names = expressed_ids)
    for (feature in c("x5utr", "cds", "x3utr")) {
      select <- grep(feature, colnames(curtype_pc))
      tpc[[feature]] <- rowMeans(curtype_pc[expressed_ids, select])
    }
    assign(paste(dset, "pc", "mean", sep = "_"), tpc)
    assign(paste(dset, "pc", sep = "_"), curtype_pc)
  }
  composite_pc <- composite_sizes[, 1:3] / rowSums(composite_sizes[, 1:3])

  # selection of subset of genes
  # -- at least 100 nts in each feature
  filter_feature_size <- names(which(composite_sizes[, "x5utr"] > 100 &
                                     composite_sizes[, "cds"] > 100 &
                                     composite_sizes[, "x3utr"] > 100))
  # -- at least 30 counts within gene for mono, di and tr
  genes_for_pc <- intersect(expressed_gene_ids, filter_feature_size)

  selections <- list("all" =  expressed_ids,
                     "non-short" = filter_feature_size,
                     "non-short-exp" = genes_for_pc,
                     "single" = expressed_single_cds)
  # Common settings for plots
  fig_aes <- aes(Feature, Mean, fill = Type)
  fig_cols <- c(di_col, mono_col, tr_col, "#ff7f00")
  pos_dodge <- position_dodge(width = 0.85)

  for (selection in names(selections)) {
    sel_ids <- selections[[selection]]
    sel_pc <- data.frame(
      rbind(data.frame(type = "Disome", di_pc_mean[sel_ids, ]),
            data.frame(type = "Monosome", mono_pc_mean[sel_ids, ]),
            data.frame(type = "Total RNA", tr_pc_mean[sel_ids, ]),
            data.frame(type = "Random", composite_pc[sel_ids, ])))
    reshaped_pc <- melt(sel_pc, measure.vars = c(2, 3, 4),
                        variable.name = "feature")
    sel_pc_df <- data.frame(
      Type = rep(c("Disome", "Monosome", "totalRNA", "size"), each = 3),
      Feature = rep(c("5'UTR", "CDS", "3'UTR"), times = 4),
      Mean = aggregate(value ~ feature + type,
                     reshaped_pc, mean, na.rm = T)$value,
      SD = aggregate(value ~ feature + type,
                     reshaped_pc, sd, na.rm = T)$value,
      SE = aggregate(value ~ feature + type,
                     reshaped_pc, sd, na.rm = T)$value / sqrt(length(sel_ids)))
    sel_pc_df$Feature <- ordered(sel_pc_df$Feature,
                                levels(sel_pc_df$Feature)[c(2, 3, 1)])
    sel_pc_df$Type <- ordered(sel_pc_df$Type,
                              levels(sel_pc_df$Type)[c(1, 2, 4, 3)])


    # Plot
    figbase <- ggplot(sel_pc_df, fig_aes)
    res_fig <- get_res_file(paste("Fig_read_dist_features_",
                                   selection, ".pdf", sep = ""))
    fig <- figbase +
             geom_bar(stat = "identity", position = pos_dodge, width = 0.8) +
             geom_errorbar(aes(ymax = Mean + 2 * SE, ymin = Mean - 2 * SE),
                               position = pos_dodge, width = 0.2) +
             scale_fill_manual(values = fig_cols) +
             theme_grey(base_size = 18)
    ggsave(res_fig, fig, useDingbats = F)
    log_n_register(res_fig,
      paste("Barplot of read distribution to transcript features for",
            selection))
    obj_name <- paste(selection, "pc_df", sep = "_")
    assign(obj_name, sel_pc_df)
    register_data_obj(obj_name, registered_env = data_env,
                      run_ver = run_version, by_name = T)
  }

  # ---- Data registration ----
  register_data_obj(di_pc_mean, mono_pc_mean, tr_pc_mean, composite_pc,
                    registered_env = data_env, run_ver = run_version)
  # Garbage collection
  rm(cdf, cdf_x3utr, cdf_x5utr, cdf_cds, curtype_pc, tpc, pc_cols,
     filter_feature_size, genes_for_pc, selections, fig_aes, fig_cols,
     pos_dodge, sel_ids, sel_pc, reshaped_pc, sel_pc_df, figbase, fig, res_fig)
  register_module("Percent", registered_env = data_env, run_ver = run_version)

}

# ==== Trinucleotide plots =====================================================

if (do_module("Trinucleotide")) {
  pass_dependencies("Trinucleotide", c("Setup", "Gene_lists"))
  logger("Calculating and plotting trinucleotide periodicity")
  # Hardcoded variables are: 21nt margin, 400nt window, x/ylim ranges
  # These need to be all changed by hand in different sections below if
  # one is modified in one place...
  if (is_option_set("Trinucleotide", "merge-treatments")){
    merge_treatments <- TRUE
    local_treatments <- c("All")
  } else {
    merge_treatments <- FALSE
    local_treatments <- treatments
  }
  if (is_option_set("Trinucleotide", "run-scripts")) {
    logger("running scripts to generate raw data for trinucleotide analysis",
           level = 2)
    norm_density_prog <- "norm_density.py"
    norm_density_ver <- system(paste(norm_density_prog, "--version"),
                               intern = T)
    if (package_version(strsplit(norm_density_ver, " ")[[1]][2]) >= 1) {
      write.table(gene_list[["robust_expressed_single_trs"]],
                  file = "robust.single.trs.txt",
                  quote = FALSE, row.names = FALSE, col.names = FALSE)
      common_args <- paste(
        "--density-db-file density_files_db_autogen.tsv",
        "--tr-list robust.single.trs.txt",
        "--specs-file", get_option("Trinucleotide", "specs"),
        "--cds", get_option("Setup", "cds_models"),
        "--region '*cds-21:cds*+21' --rep 1,2")
      base_struc <- list(
        c("monosome", "mono", "std", "--type Mono"),
        c("disome", "disome", "std", "--type Di"),
        c("disome", "disome", "15",
          "--type Di --offset 15 --length 59,60,62,63"))
      struc <- list()
      for (treatment in local_treatments) {
        for (d_args in base_struc) {
          d_args[1] <- paste(treatment, d_args[1], sep = "/")
          d_args[2] <- paste(treatment, d_args[2], sep = "_")
          if (merge_treatments) {
            d_args[4] <- paste(d_args[4], "--treatments",
                               paste(treatments, collapse = ","))
          } else {
            d_args[4] <- paste(d_args[4], "--treatments", treatment)
          }
          struc[[length(struc) + 1]] <- d_args
        }
      }
      for (d_args in struc) {
        logger(paste("generating normalized", d_args[1],
                     "densities, Asite =", d_args[3]), level = 3)
        res_file <- get_res_file(paste("Data", d_args[2], d_args[3],
                                       "all_norm_density.txt", sep = "_"))
        logger(paste("Running:", common_args, d_args[4]), level = 3)
        args <- paste(common_args, d_args[4], ">", res_file)
        if (system(paste(norm_density_prog, args)) == 0) {
          log_n_register(res_file, paste(
            "Normalized densities/position for", d_args[1], d_args[3],
            "A-site corr reads using", norm_density_ver),
                         lvl = 3, res_type = "Data")
        } else {
          logger(paste("Could not successfully create", res_file),
                 level = 3, is.message = T)
        }
      }
    # Garbage collection
    rm(norm_density_prog, norm_density_ver, common_args, struc,
       d_args, res_file, args)
    } else {
      logger(paste("version found:", norm_density_ver, "is not sufficient"),
             is.message = T, level = 3)
    }
  }
  logger("Reading normalized density data", level = 2)
  tri_struc <- list("mono" = list(name = "mono", corr = c("std")),
                    "di" = list(name = "disome", corr = c("15", "std")))
  dirs <- list("5p" = T, "3p" = F)
  y_max <- 0
  signalp_trs <- gene_list[["signalp_trs"]]$tr_id
  for (treatment in local_treatments) {
    for (dset in names(tri_struc)) {
      dname <- tri_struc[[dset]]$name
      corrections <- tri_struc[[dset]]$corr
      for (corr in corrections) {
        res_file <- get_res_file(
          paste("Data", treatment, dname, corr,
                "all_norm_density.txt", sep = "_"))
        nd <- read.table(res_file)
        for (direction in names(dirs)) {
          logger(paste("Extracting 400nt from", direction, "end of", dname,
                       "/", treatment, "with", corr, "correction"), level = 3)
          nd_400 <- extract_regions(nd, from_the_beginning = dirs[[direction]],
                                    size = 423)
          nd_400_sig <- subset(nd_400, V1 %in% signalp_trs)
          nd_400_nosig <- subset(nd_400, V1 %ni% signalp_trs)
          nd_400_mean <- aggregate(V3 ~ V2, nd_400, mean)
          nd_400_sig_mean <- aggregate(V3 ~ V2, nd_400_sig, mean)
          nd_400_nosig_mean <- aggregate(V3 ~ V2, nd_400_nosig, mean)
          y_max <- max(y_max, nd_400_mean$V3,
                       nd_400_sig_mean$V3, nd_400_nosig_mean$V3)
          assign(paste(dset, "density", treatment, corr, direction,
                       "400", sep = "_"), nd_400_mean)
          assign(paste(dset, "density", treatment, corr, direction,
                       "400", "sig", sep = "_"), nd_400_sig_mean)
          assign(paste(dset, "density", treatment, corr, direction,
                       "400", "nosig", sep = "_"), nd_400_nosig_mean)
        }
      }
    }
  }
  logger("Plotting disome vs monosome trinucleotide plots", level = 2)
  plot_args <- list("5p" = list("xlim" = c(-20, 200),
                                "ylim" = c(0, y_max), "corr" = -22),
                    "3p" = list("xlim" = c(-200, 20),
                                "ylim" = c(0, y_max), "corr" = -403))
  for (treatment in local_treatments) {
    for (direction in c("5p", "3p")) {
      for (corr in c("std", "15")) {
        cur_args <- plot_args[[direction]]
        res_fig <- get_res_file(paste("Fig_trinuc_mono_vs_di", treatment,
                                      corr, direction, "200nt.pdf", sep = "_"))
        pdf(res_fig)
        with(get(paste("mono", "density", treatment, "std",
                       direction, "400", sep = "_")),
          plot(V2 + cur_args$corr, V3, type = "l",
               xlab = "Position relative to start site (nt)",
               ylab = "Mean normalized read density",
               xlim = cur_args$xlim, ylim = cur_args$ylim))
        with(get(paste("di", "density", treatment, corr,
                       direction, "400", sep = "_")),
          lines(V2 + cur_args$corr, V3, col = "red"))
        dev.off()
        log_n_register(res_fig, paste(
          "Trinucleotide periodicity of mono vs disome for ", treatment,
          "(", corr, ") within 200nt at ", direction, sep = ""), lvl = 3)
      }
    }
  }
  logger("Plotting signal vs non-signal-peptide transcripts", level = 2)
  obj_names <- c()
  plot_args <- list("5p" = list("xlim" = c(-20, 400),
                                "ylim" = c(0, y_max), "corr" = -22),
                    "3p" = list("xlim" = c(-400, 20),
                                "ylim" = c(0, y_max), "corr" = -403))
  for (treatment in local_treatments) {
    for (direction in c("5p", "3p")) {
      cur_args <- plot_args[[direction]]
      for (corr in c("std", "15")) {
        for (dset in c("mono", "di")) {
          if (dset == "mono" & corr != "std") next
          res_fig <- get_res_file(
            paste("Fig_trinuc_signalp", dset, treatment,
                  corr, direction, "400nt.pdf", sep = "_"))
          pdf(res_fig)
          dset_nosig <- paste(dset, "density", treatment, corr, direction,
                              "400", "nosig", sep = "_")
          dset_sig <- paste(dset, "density", treatment, corr, direction,
                            "400", "sig", sep = "_")
          obj_names <- c(obj_names, dset_nosig, dset_sig)
          with(get(dset_nosig),
            plot(V2 + cur_args$corr, V3, type = "l",
                 xlab = "Position relative to start site (nt)",
                 ylab = "Mean normalized read density",
                 xlim = cur_args$xlim, ylim = cur_args$ylim))
          with(get(dset_sig),
            lines(V2 + cur_args$corr, V3, col = "red"))
          dev.off()
          log_n_register(res_fig, paste(
            "Trinucleotide periodicity of signal peptide - mono vs disome(",
            corr, ") within 400nt at ", direction, sep = ""), lvl = 3)
        }
      }
    }
  }
  # ---- Data registration ----
  for (obj_name in obj_names) {
    register_data_obj(obj_name, registered_env = data_env,
                      run_ver =  run_version, by_name = TRUE)
  }
  # Garbage collection
  rm(tri_struc, dirs, dset, dname, corrections, corr, nd, nd_400, nd_400_mean,
     nd_400_sig, nd_400_sig_mean, plot_args, direction, cur_args, dset_nosig,
     dset_sig, local_treatments, merge_treatments, obj_names, obj_name, res_fig)
  register_module("Trinucleotide", registered_env = data_env,
                  run_ver = run_version)

}

# ==== Mono vs Di global comparisons ===========================================

if (do_module("Di_vs_mono")) {
  pass_dependencies("Di_vs_mono", c("Setup", "Filter",
                                    "Normalize", "Gene_lists"))
  count_label <- get_option("Di_vs_mono", "count-label")
  if (is.null(count_label)) count_label <- "original"
  tag <- count_tables[[count_label]][["tag"]]
  idex <- if (tag != "") paste("_", tag, sep = "") else ""
  logger(paste("Performing global monosome vs disome comparisons for",
               count_label))
  struc <- list(
    list(name = "SignalP", short = "sig", gids = gene_list[["signalp_ids"]]),
    list(name = "conservation_Q1", short = "phylop_q1",
         gids = gene_list[["phylop_Q1_ids"]]),
    list(name = "conservation_Q2", short = "phylop_q2",
         gids = gene_list[["phylop_Q2_ids"]]),
    list(name = "conservation_Q3", short = "phylop_q3",
         gids = gene_list[["phylop_Q3_ids"]]),
    list(name = "conservation_Q4", short = "phylop_q4",
         gids = gene_list[["phylop_Q4_ids"]]))
  if (get_option("Di_vs_mono", "merge-treatments")) {
    local_treatments <- "All"
    merge_treatments <- TRUE
  } else {
    local_treatments <- treatments
    merge_treatments <- FALSE
  }
  obj_names <- c()
  for (treatment in local_treatments) {
    logger(paste("analysing", treatment), level = 2)
    if (merge_treatments) {
      treat_cols <- design_flat == design_flat
    } else {
      treat_cols <- design_flat == treatment
    }
    cur_eff_mono <- get(paste("eff_mono", idex, "_log_noinf", sep = ""))
    cur_eff_di <- get(paste("eff_di", idex, "_log_noinf", sep = ""))
    eff_mono_mean <- rowMeans(cur_eff_mono[, treat_cols], na.rm = T)
    eff_mono_sd <- apply(cur_eff_mono[, treat_cols], 1, sd, na.rm = T)
    eff_di_mean <- rowMeans(cur_eff_di[, treat_cols], na.rm = T)
    eff_di_sd <- apply(cur_eff_di[, treat_cols], 1, sd, na.rm = T)
    eff_df <- data.frame("mono_mean" = eff_mono_mean, "di_mean" = eff_di_mean,
                         "mono_sd" = eff_mono_sd, "di_sd" = eff_di_sd)
    logger("fitting Deming regression line for all genes", level = 3)
    eff_deming_all_std <- deming(di_mean ~ mono_mean, data = eff_df,
                                 xstd = mono_sd, ystd = di_sd)
    # Objects to register and keep
    count_treat <- paste(treatment, idex, sep = "")
    df_name <- paste("eff_df", count_treat, sep = "_")
    assign(df_name, eff_df)
    deming_name <- paste("eff_deming_all", count_treat, sep = "_")
    assign(deming_name, eff_deming_all_std)
    obj_names <- c(obj_names, c(df_name, deming_name))
    for (genelist in struc) {
      # Scatterplot with deming lines
      res_fig <- get_res_file(
        paste("Fig", count_treat, "eff_mono_vs_di", genelist$name,
              "deming_lines.pdf", sep = "_"))
      logger(
        paste("fitting Deming regression line for", genelist$name), level = 3)
      sub_deming <- deming(di_mean ~ mono_mean, data = eff_df,
                           xstd = mono_sd, ystd = di_sd,
                           subset = rownames(eff_df) %in% genelist$gids)
      deming_name <- paste("eff_deming", genelist$short, count_treat, sep = "_")
      assign(deming_name, sub_deming)
      obj_names <- c(obj_names, deming_name)
      pdf(res_fig, width = 6, height = 7)
      eqscplot(eff_mono_mean, eff_di_mean, pch = 16, col = rgb(0, 0, 0, 0.2),
               cex = 0.6, xlab = paste("Mean monosome TE"),
               ylab = paste("Mean disome TE"), main = treatment, type = "n")
      pretty_grid(1)
      points(eff_mono_mean[expressed_ids %ni% genelist$gids],
             eff_di_mean[expressed_ids %ni% genelist$gids],
             pch = 16, col = rgb(0, 0, 0, 0.2), cex = 0.6)
      points(eff_mono_mean[genelist$gids],
             eff_di_mean[genelist$gids],
             pch = 16, col = rgb(1, 0, 0, 0.4), cex = 0.6)
      abline(0, 1, col = "lightgrey", lt = 2)
      abline(coef(eff_deming_all_std), col = "black")
      abline(coef(sub_deming), col = "red")
      legend("topleft",
             c(paste("All:", round(coef(eff_deming_all_std)[2], 3),
                     ", 95% CI:", paste(round(eff_deming_all_std$ci[2, ], 3),
                                        collapse = "-")),
               paste(paste(genelist$name, ":", sep = ""),
                     round(coef(sub_deming)[2], 3),
                     ", 95% CI:", paste(round(sub_deming$ci[2, ], 3),
                                        collapse = "-"))),
             col = c("black", "red"), pch = 16, lty = 1, bty = "n")
      dev.off()
      log_n_register(res_fig, paste(
        "Scatterplot of mono vs disome mean TEs with Deming trend lines :",
        count_label, treatment, genelist$name), lvl = 3)
      # Same plot with marginal densities
      res_fig <- get_res_file(
        paste("Fig", count_treat, "eff_mono_vs_di", genelist$name,
              "marginal_density.pdf", sep = "_"))
      tdf <- data.frame(x = eff_mono_mean,
                        y = eff_di_mean,
                        k = "All", stringsAsFactors = F)
      tdf[genelist$gids, ]$k <- genelist$short
      mcols <- c("#00000033", "#FF000066")
      cairo_pdf(res_fig, width = 6, height = 7)
      scatter_with_marginal_densities(
        tdf, c("All", genelist$short), my.cols = mcols, margin.xy = "xy",
        ylim = c(-8, 6), log.xy = "", main = treatment, multiplier.x = 1.6,
        multiplier.y = 1.7, lwd.den = 2, den.solid = F, cex = 0.8,
        xlab = "Normalized monosome density, TE (log2)",
        ylab = "Normalized disome density (log2)")
      abline(0, 1, col = "grey", lt = 2)
      abline(coef(eff_deming_all_std), col = "black")
      abline(coef(sub_deming), col = "red")
      legend(-7, 5.8, c("", ""), col = c("red", "red"),
             pch = 16, lty = 1, bty = "n", cex = 1.2)
      legend(-7, 5.8,
             c(paste("All (N = ", length(expressed_ids), ")", sep = ""),
               paste(genelist$name, " (N = ", length(genelist$gids), ")",
                     sep = "")),
             col = c("black", "red"),
             pch = c("\u25D0", NA), lty = 1, bty = "n", cex = 1.2, pt.cex = 0.8)
      dev.off()
      log_n_register(res_fig, paste(
        "Scatterplot of mono vs disome mean TEs with marginal densities :",
        count_label, treatment, genelist$name), lvl = 3)
    }
    # Mono vs di TE with conservation coloring
    if (is_option_set("Di_vs_mono", "conservation")) {
      res_fig <- get_res_file(
        paste("Fig", count_treat, "eff_mono_vs_di-conservation.pdf", sep = "_"))
      expressed_cons <- subset(phylop_table, gene_id %in% expressed_ids)
      col_pal <- viridis(nrow(expressed_cons))
      colors <- col_pal[rank(expressed_cons$con_score)]
      pdf(res_fig)
      eqscplot(eff_mono_mean[expressed_cons$gene_id],
               eff_di_mean[expressed_cons$gene_id],
               pch = 16, col = colors,
               cex = 0.6, xlab = paste("Mean monosome TE"),
               ylab = paste("Mean disome TE"), main = treatment, type = "n")
      pretty_grid(1)
      points(eff_mono_mean[expressed_cons$gene_id],
             eff_di_mean[expressed_cons$gene_id],
             pch = 16, col = colors, cex = 0.6)
      abline(0, 1, col = "lightgrey", lt = 2)
      dev.off()
      log_n_register(res_fig,
                      paste("Scatterplot of mono vs disome mean TE colored by ",
                            "PhyloP scores :", count_label, treatment), lvl = 3)
    }
    # Mono TE vs di/mono ratio
    res_fig <- get_res_file(
      paste("Fig", count_treat, "eff_mono_vs_di-mono_ratio.pdf", sep = "_"))
    cur_norm_mono <- get(paste("norm_cds_mono", idex, sep = ""))
    cur_norm_di <- get(paste("norm_cds_di", idex, sep = ""))
    di_mono_ratio <- 
      rowMeans(no.inf(log(cur_norm_di[expressed_ids, treat_cols] /
                          cur_norm_mono[expressed_ids, treat_cols])), na.rm = T)
    di_mono_lm <- lm(di_mono_ratio ~ eff_mono_mean)
    pdf(res_fig)
    plot(eff_mono_mean, di_mono_ratio, pch = 16, col = rgb(0, 0, 0, 0.2),
         cex = 0.6, xlab = paste("Mean monosome TE"),
         ylab = paste("Mean log(disome/monosome)"), type = "n")
    pretty_grid(1)
    points(eff_mono_mean, di_mono_ratio,
           pch = 16, col = rgb(0, 0, 0, 0.2), cex = 0.6)
    abline(h = 0, col = "lightgrey", lt = 2)
    abline(di_mono_lm)
    legend("topleft", lwd = 1, bty = "n", cex = 0.6,
           legend = c(as.expression(lmp(di_mono_lm, "linear regression"))))
    dev.off()
    log_n_register(res_fig, paste("Scatterplot of mean disome/monosome ratio",
      "vs mono TEs :", treatment), lvl = 3)
    di_mono_name <- paste("di_mono_ratio", count_treat, sep = "_")
    assign(di_mono_name, di_mono_ratio)
    obj_names <- c(obj_names, di_mono_name)
  }
  # ZT12 vs ZT2 density plots
  if (is_option_set("Di_vs_mono", "zt-density")) {
    res_fig <- get_res_file(paste("Fig_di_mono_ratio_density_btw_ZT12_2",
                                  idex, ".pdf", sep = ""))
    pdf(res_fig)
    zt_df <- na.omit(get(paste("di_mono_ratio_ZT12", idex, sep = "")) -
                       get(paste("di_mono_ratio_ZT02", idex, sep = "")))
    zt_grp <- names(zt_df) %in% union(gene_list[["Rpl_ids"]],
                                         gene_list[["Rps_ids"]])
    sm.density.compare(as.numeric(zt_df), zt_grp, model="equal", xlim=c(-2, 2),
                       col=c("black","red"), col.band="lightgray",
                       xlab = "log(delta(disome/monosome ZT12 - ZT2))",
                       nboot = 10000)
    abline(v=0, lty=2, col="gray")
    legend("topleft", c("Rp(l/s) N=57", "others N=8558"), col=c("red", "black"),
           lty=c(2,1), bty="n")
    dev.off()
    log_n_register(
      res_fig, paste(
        "Comparison of kernel densities of delta di/mono ZT12 vs ZT2"), lvl = 2)
  }
  # ---- Data registration ----
  for (obj_name in obj_names) {
    register_data_obj(obj_name, registered_env = data_env,
                      run_ver =  run_version, by_name = TRUE)
  }
  # Garbage collection
  rm(eff_mono_mean, eff_di_mean, eff_mono_sd, eff_di_sd, eff_df, treat_cols,
     df_name, deming_name, eff_deming_all_std, sub_deming, tdf, mcols, res_fig,
     di_mono_ratio, di_mono_lm, di_mono_name, obj_names)
  register_module("Di_vs_mono", registered_env = data_env,
                  run_ver = run_version)

}

# ==== Treatment vs treatment scatterplots =====================================

if (do_module("Treatment_plots")) {
  pass_dependencies("Treatment_plots", c("Setup", "Filter",
                                         "Normalize", "Gene_lists"))
  logger("Plotting treatment scatters with various highlights")
  datasets <- list(
    "rpkm" = list(varname = "rpkm_cds_{}_mean",
                  log_it = TRUE, subs = c("tr", "mono", "di")),
    "norm" = list(varname = "norm_cds_{}_mean",
                  log_it = TRUE, subs = c("tr", "mono", "di")),
    "eff" = list(varname = "eff_{}_log_noinf_mean",
                 log_it = FALSE, subs = c("mono", "di")))
  comparisons <- list(
    list(xaxis = "ZT00", yaxis = "ZT12", datasets = c("rpkm", "norm", "eff"),
         highlight = c("all", "signalp_ids"),
         abline = TRUE),
    list(xaxis = "ZT00", yaxis = "ZT02", datasets = c("rpkm", "norm", "eff"),
         highlight = c("all", "signalp_ids"),
         abline = TRUE),
    list(xaxis = "ZT02", yaxis = "ZT12", datasets = c("rpkm", "norm", "eff"),
         highlight = c("all", "signalp_ids"),
         abline = TRUE))
  for (item in comparisons) {
    x_treat <- item$xaxis
    y_treat <- item$yaxis
    for (dataset in item$datasets) {
      cur_dataset <- datasets[[dataset]]
      for (dtype in cur_dataset$subs) {
        cur_data <- get(gsub("\\{\\}", dtype, cur_dataset$varname))
        for (highlight in item$highlight) {
          res_fig <- get_res_file(
            paste("Fig", dataset, dtype, x_treat, "vs", y_treat, highlight,
                  "scatterplot.pdf", sep = "_"))
          if (highlight == "all") {
            all_ids <- expressed_ids
          } else {
            highlight_ids <- gene_list[[highlight]]
            all_ids <- setdiff(expressed_ids, highlight_ids)
          }
          if (cur_dataset$log_it) {
            log_arg <- "xy"
          } else {
            log_arg <- ""
          }
          pdf(res_fig, width = 6, height = 7)
          plot(cur_data[all_ids, x_treat], cur_data[all_ids, y_treat],
               xlab = paste(x_treat, dtype, dataset),
               ylab = paste(y_treat, dtype, dataset),
               type = "n", log = log_arg)
          grid()
          points(cur_data[all_ids, x_treat], cur_data[all_ids, y_treat],
                 pch = 16, col = rgb(0, 0, 0, 0.4), cex = 0.8)
          if (highlight != "all") {
            points(cur_data[highlight_ids, x_treat],
                   cur_data[highlight_ids, y_treat],
                   pch = 16, col = rgb(1, 0, 0, 0.4), cex = 0.8)
          }
          if (item$abline) {
            abline(0, 1, col = "lightgrey", lt = 2)
          }
          dev.off()
          log_n_register(res_fig,
            paste("Scatterplot of", x_treat, "vs", y_treat, dtype, dataset,
                  "with", highlight, "genes highlighted"), lvl = 3)
        }
      }
    }
  }
  register_module("Treatment_plots", registered_env = data_env,
                  run_ver = run_version)

}

# ==== GSEA and pathway analysis ===============================================

if (do_module("GSEA-gage")) {
  locally_required_packages = c("gage", "GO.db", "org.Mm.eg.db")
  local_dependency_met <- unlist(lapply(locally_required_packages, require,
                                        quietly = T, character.only = T))
  count_label <- get_option("GSEA-gage", "count-label")
  if (is.null(count_label)) count_label <- "original"
  tag <- count_tables[[count_label]][["tag"]]
  idex <- if (tag != "") paste("_", tag, sep = "") else ""
  if (get_option("GSEA-gage", "merge-treatments")) {
    local_treatments <- "All"
    merge_treatments <- TRUE
  } else {
    local_treatments <- treatments
    merge_treatments <- FALSE
  }
  if(all(local_dependency_met)) {
    logger("Performing Generally Applicable Geneset/Pathway (gage) Analysis")
    logger(paste("using", count_label))
    gs_struc <- list()
    # KEGG Database
    if (is_option_set("GSEA-gage", "KEGG")) {
      logger("extracting genesets from KEGG DB", level = 2)
      kegg_mmu <- kegg.gsets("mmu")
      kegg_gs <- kegg_mmu$kg.sets[kegg_mmu$sigmet.idx]
      gs_struc <- c(gs_struc, list("kegg" = list(name = "KEGG", gs = kegg_gs))
      )
      rm(kegg_mmu, kegg_gs)
    }
    # GO Database
    if (is_option_set("GSEA-gage", "GO")) {
      logger("extracting genesets from GO DB", level = 2)
      go_mmu <- go.gsets("mouse")
      go_bp_gs <- go_mmu$go.sets[go_mmu$go.subs$BP]
      go_cc_gs <- go_mmu$go.sets[go_mmu$go.subs$CC]
      go_mf_gs <- go_mmu$go.sets[go_mmu$go.subs$MF]
      gs_struc <- c(gs_struc, list(
        "go_bp" = list(name = "GO.BioProcess", gs = go_bp_gs),
        "go_cc" = list(name = "GO.CellComponent", gs = go_cc_gs), 
        "go_mf" = list(name = "GO.MolFunction", gs = go_mf_gs)))
      rm(go_mmu, go_bp_gs, go_cc_gs, go.mf.gs)
    }
    # Molecular Signatures Database
    if (is_option_set("GSEA-gage", "MSigDB")) {
      logger("extracting genesets from MSigDB", level = 2)
      msigdb_file <- file.path(
        get_option("GSEA-gage", "folder"),
        get_option("GSEA-gage", "msig-db-file"))
      logger(paste("using", msigdb_file), level = 3)
      load(msigdb_file)
      gs_struc <- c(gs_struc, list(
        "msig_h" = list(name = "MSig.Hallmark", gs = Mm.H.gs), 
        "msig_c2" = list(name = "MSig.Curated", gs = Mm.C2.gs), 
        "msig_c3" = list(name = "MSig.Motif", gs = Mm.C3.gs), 
        "msig_c4" = list(name = "MSig.Computation", gs = Mm.C4.gs), 
        "msig_c5" = list(name = "MSig.GO", gs = Mm.C5.gs), 
        "msig_c6" = list(name = "MSig.Oncogenic", gs = Mm.C6.gs), 
        "msig_c7" = list(name = "MSig.Immunologic", gs = Mm.C7.gs)))
      rm(Mm.H.gs, Mm.C2.gs, Mm.C3.gs, Mm.C4.gs, Mm.C5.gs, Mm.C6.gs, Mm.C7.gs)
    }
    gs_ids <- c()
    gs_all <- list()
    for (ds in names(gs_struc)) {
      gs_all <- c(gs_all, gs_struc[[ds]]$gs)
      gs_ids <- union(gs_ids, unlist(gs_struc[[ds]]$gs))
    }
    gs_struc[["all"]] <- list(name = "All.Pathways", gs = gs_all)
    cur_exp_ids <- get(paste("expressed_ids", idex, sep = ""))
    expressed_gage_ids <- unlist(lapply(cur_exp_ids, subs_ids_geneset,
                                        expressed_entrez, gs_ids))
    print(length(expressed_gage_ids))
    cur_eff_mono <- get(paste("eff_mono", idex, "_log_noinf", sep = ""))
    cur_eff_di <- get(paste("eff_di", idex, "_log_noinf", sep = ""))
    for (treatment in local_treatments) {
      logger(paste("analysing", treatment), level = 2)
      if (merge_treatments) {
        treat_cols <- design_flat == design_flat
      } else {
        treat_cols <- design_flat == treatment
      }
      count_treat <- paste(treatment, idex, sep = "")
      eff_df <- get(paste("eff_df", count_treat, sep = "_"))
      eff_deming <- get(paste("eff_deming_all", count_treat, sep = "_"))
      if (is_option_set("GSEA-gage", "rotate")) {
        df_rotator <- function (x, y, ...) {
          for (i in 1:ncol(x)) {
            r <- rotate_points(x[, i], y[, i], ...)
            x[, i] <- r[, 1]
            y[, i] <- r[, 2]
          }
          cbind(x, y)
        }
        logger("rotating efficiency data around Deming regression",
               level = 2)
        deming_coef <- coef(eff_deming)
        p_x <- deming_coef[1] / (1 - deming_coef[2])
        theta <- tanh((deming_coef[2] - 1) / (1 + deming_coef[2]))
        eff_gage_df <- df_rotator(cur_eff_mono[, treat_cols],
                                  cur_eff_di[, treat_cols],
                                  point = c(p_x, p_x), theta = theta)
      } else {
        eff_gage_df <- cbind(cur_eff_mono[, treat_cols],
                             cur_eff_di[, treat_cols])
      }
      rownames(eff_gage_df) <- expressed_gage_ids
      col_width <- sum(treat_cols)
      gage_mono <- 1:col_width
      gage_di <- (col_width + 1):(2 * col_width)
      gage_sig_gs <- list()
      for (ds_n in names(gs_struc)) {
        ds <- gs_struc[[ds_n]]
        dir.create(file.path(res_folder, ds$name), showWarnings = F)
        logger(paste("generating gene lists/plots for significant",
                     ds$name, "pathways"), level = 2)
        eff_gage <- gage(eff_gage_df, gsets = ds$gs,
                         ref = gage_mono, samp = gage_di, compare = "as.group")
        ds_sig_list <- list()
        for (trend in c("greater", "less")) {
          fig_base <- paste("Fig", ds$name, trend, "", sep = "_")
          ds_sig_list[[trend]] <- list()
          pathways <- names(which(eff_gage[[trend]][, 4] < 1e-04))
          for (pathway in pathways) {
            p_ids <- cur_exp_ids[which(
              rownames(eff_gage_df) %in% ds$gs[[pathway]])]
            ds_sig_list[[trend]][[pathway]]$gid <- p_ids
            p_deming <- deming(di_mean ~ mono_mean, data = eff_df,
                               xstd = mono_sd, ystd = di_sd,
                               subset = rownames(eff_df) %in% p_ids)
            ds_sig_list[[trend]][[pathway]]$deming <- p_deming
            res_fig <- paste(fig_base,
                             make_file_name(pathway), ".pdf", sep = "")
            res_fig <- file.path(ds$name, res_fig)
            res_fig <- get_res_file(res_fig)
            pdf(res_fig, width = 6, height = 7)
            eqscplot(eff_df$mono_mean, eff_df$di_mean, pch = 16,
                     col = rgb(0, 0, 0, 0.2), cex = 0.6,
                     xlab = paste("Mean monosome TE"),
                     ylab = paste("Mean disome TE"), type = "n")
            pretty_grid(1)
            points(eff_df[cur_exp_ids %ni% p_ids, ]$mono_mean,
                   eff_df[cur_exp_ids %ni% p_ids, ]$di_mean,
                   pch = 16, col = rgb(0, 0, 0, 0.2), cex = 0.6)
            points(eff_df[p_ids, ]$mono_mean, eff_df[p_ids, ]$di_mean,
                   pch = 16, col = rgb(1, 0, 0, 0.7), cex = 0.6)
            abline(0, 1, col = "lightgrey", lt = 2)
            abline(coef(eff_deming), col = "black")
            abline(coef(p_deming), col = "red")
            legend("topleft",
                   c(paste("All (", length(cur_exp_ids), "):",
                           round(coef(eff_deming)[2], 3), ", 95% CI:",
                           paste(round(eff_deming$ci[2, ], 3),
                                 collapse = "-")),
                     paste(pathway, "(", length(p_ids), ")",
                           round(coef(p_deming)[2], 3), ", 95% CI:",
                           paste(round(p_deming$ci[2, ], 3),
                                 collapse = "-"))),
                   col = c("black", "red"), pch = 16,
                   lty = 1, bty = "n", cex = 0.6)
            dev.off()
            log_n_register(res_fig,
                           "Scatterplot of mono vs disome mean TEs with Deming trend lines",
                           lvl = 3)
          }
        }
        gage_sig_gs[[ds_n]] <- ds_sig_list
        # ---- Data registration of current gage results ----
        eff_obj_name <- paste("eff_gage", ds_n, count_treat, sep = "_")
        assign(eff_obj_name, eff_gage)
        register_data_obj(eff_obj_name, registered_env = data_env,
                          run_ver = run_version, by_name = T)
      }
      # eff_gage_df is not registered
      # gage_sig_gs is not registered
    }
    # ---- Data registration ----
    register_data_obj(gs_struc,
                      registered_env = data_env, run_ver =  run_version)
    # Garbage collection
    suppressWarnings(
      rm(gs_ids, gs_all, expressed_gage_ids, gage_mono, gage_di, pathways,
         p_ids, locally_required_packages, local_dependency_met, pathway,
         p_deming, res_fig, fig_base, ds_sig_list, trend, ds, ds_n,
         eff_obj_name, khier, korg)
    )
    register_module("GSEA-gage", registered_env = data_env,
                    run_ver = run_version)
  } else {
    logger("Can NOT perform gage analysis - required libs could not be loaded",
           is.message = TRUE)
  }
}

# ==== Distance analysis =======================================================

if (do_module("Distance_analysis")) {
  pass_dependencies("Distance_analysis", c("Setup", "Gene_lists"))
  logger("Calculation of complete distance analysis of features")
  pcf <- function(x) {
    x$V1 / (x$V1 - x$V2)
  }
  obj_names <- c()
  structs <- c("3,30,Structured,6,30,Unstructured,3,30,Structured",
               "3,30,Unstructured,6,30,Structured,3,30,Unstructured")
  struct_types <- c("std", "rev")
  bootstrap_n <- get_option("Distance_analysis", "bootstrap-N")
  what_to_do <- get_option("Distance_analysis", "what-to-do")
  if ("run-scripts" %in% what_to_do) {
    logger("running scripts to generate raw data for distance analysis",
           level = 2)
    dist_prog <- "distance_analyzer.py"
    dist_ver <- "version 1.1.1" # system(paste(dist_prog, "--version"),
                  #        intern = T)
    if (package_version(strsplit(dist_ver, " ")[[1]][2]) >= 1) {
      common_args <- c(
        "--mode", "normalized", "--density-db-file",
        "density_files_db_autogen.tsv", "--specs-file",
        get_option("Distance_analysis", "specs"),
        "--treatments", paste(treatments, collapse = ","),
        "--cds", get_option("Setup", "cds_models"),
        "--fasta", get_option("Distance_analysis", "fasta"),        
        "secondary_struct_simple")
      for (read_type in c("Tr", "Mono", "Di")) {
        for (struct in 1:length(structs)) {
          logger(paste("launching", dist_prog, "for", read_type, "using",
                       structs[struct]), level = 3)
          args = c("--bootn", bootstrap_n, "--type", read_type,
                   "--mask", structs[struct], common_args)
          if ("randomized" %in% what_to_do) {
            system2(dist_prog, args = c("--randomized", args), wait = TRUE)
            raw <- list()
            file_base <- paste(
              "raw", paste(strsplit(structs[struct], split = ",")[[1]],
                           collapse = "_"), "rand", sep = "_")
            for (i in 1:bootstrap_n - 1) {
              raw[[i+1]] <- read.table(file = paste(file_base, i, sep = "_"))
            }
            obj_name <- paste("raw", read_type, struct_types[struct], "rand",
                              sep = "_")
            assign(obj_name, raw)
            obj_names <- c(obj_names, obj_name)
            system2("rm", paste(file_base, "*", sep = ""))
          }
          if ("data" %in% what_to_do) {
            system2(dist_prog, args = args, wait = TRUE)
            file_base <- paste(
              "raw", paste(strsplit(structs[struct], split = ",")[[1]],
                           collapse = "_"), sep = "_")
            raw <- read.table(file = paste(file_base, 0, sep = "_"))
            obj_name <- paste("raw", read_type, struct_types[struct], sep = "_")
            assign(obj_name, raw)
            obj_names <- c(obj_names, obj_name)
            system2("rm", paste(file_base, "*", sep = ""))
          }
        }
      }
    }
  }
  if ("calc" %in% what_to_do) {
    logger("Analysing raw distance data", level = 2)
    for (read_type in c("Tr", "Mono", "Di")) {
      for (struct in 1:length(structs)) {
        logger(paste("processing raw data for", read_type, "using",
                     structs[struct]), level = 3)
        logger("calculating joined edge percentages", level = 4)
        raw_obj <- paste("raw", read_type, struct_types[struct], sep = "_")
        raw_dat <- get(raw_obj)
        raw_rand <- get(paste(raw_obj, "rand", sep = "_"))
        den_rand <- sapply(raw_rand, function(x) {
          density(pcf(x), from = -2, to = 3, weights = x$V3 / sum(x$V3))$y })
        den <- density(pcf(raw_dat), from = -2, to = 3,
                       weights = raw_dat$V3 / sum(raw_dat$V3))
        den_quantiles <- apply(den_rand, 1, quantile,
                               c(0.01, 0.025, 0.05, 0.125, 0.25, 0.5,
                                 0.75, 0.875, 0.95, 0.975, 0.99))
        den_obj <- paste("density", read_type, struct_types[struct],
                         sep = "_")
        den_rand_obj <- paste(den_obj, "quan", "rand", sep = "_")
        logger("calculating left edge probabilities", level = 4)
        left_sums <- lapply(raw_rand, function(x) {
          aggregate(V3 ~ V1, subset(x, V2 < 0), sum)
        })
        minmax <- sapply(left_sums, function(x) range(x$V1))
        sum_range <- min(minmax):max(minmax)
        left_res <- data.frame(matrix(, nrow = length(sum_range), ncol = 0),
                               row.names = sum_range)
        for (i in 1:length(left_sums)) {
          left_res[as.character(left_sums[[i]]$V1), i] <- left_sums[[i]]$V3
        }
#        left_res <- t(left_res[order(as.numeric(rownames(left_res))),])
        left_res <- t(left_res)
        left_obj <- paste("left_junc", read_type, struct_types[struct],
                          sep = "_")
        logger("calculating right edge probabilities", level = 4)
        right_sums <- lapply(raw_rand, function(x) {
          aggregate(V3 ~ V2, subset(x, V1 >= 0), sum)
        })
        minmax <- sapply(right_sums, function(x) range(x$V2))
        sum_range <- min(minmax):max(minmax)
        right_res <- data.frame(matrix(, nrow = length(sum_range), ncol = 0),
                                row.names = sum_range)
        for (i in 1:length(right_sums)) {
          right_res[as.character(right_sums[[i]]$V2), i] <- right_sums[[i]]$V3
        }
#        right_res <- t(right_res[order(as.numeric(rownames(right_res))),])
        right_res <- t(right_res)
        right_obj <- paste("right_junc", read_type, struct_types[struct],
                           sep = "_")
        assign(den_obj, den)
        assign(den_rand_obj, den_quantiles)
        assign(left_obj, left_res)
        assign(right_obj, right_res)
        obj_names <- c(obj_names, den_obj, den_rand_obj, left_obj, right_obj)
      }
    } 
  }
  if ("clean-rand" %in% what_to_do) {
    logger("cleaning-up large randomization-raw data objects", level = 2)
    big_len <- 3 * length(struct_types)
    big_objs <- paste(rep("raw", big_len), c("Tr", "Mono", "Di"), struct_types,
                      rep("rand", big_len), sep = "_")
    writeLines("Object to be deleted:")
    print(big_objs)
    if (are_you_sure()) {
      for (big_obj in big_objs) {
        rm_data_obj(big_obj, registered_env = data_env, by_name = TRUE)
      }
      logger(
        "large data objects are removed from registered project data",
        level = 3)
    } else {
      logger(
        "skipped by user's choice", level = 3)
    }
  }
  if ("plot" %in% what_to_do) {
    logger("plotting distance densities", level = 2)
    quantiles <- list(
      "2%" = c("1%", "99%", "gray90"),
      "5%" = c("2.5%", "97.5%", "gray80"),
      "10%" = c("5%", "95%", "gray70"),
      "25%" = c("12.5%", "87.5%", "gray60"),
      "50%" = c("25%", "75%", "gray50"),
      "median" = c("50%", "50%", "gray40"))
    for (read_type in c("Tr", "Mono", "Di")) {
      for (struct in 1:length(structs)) {
        logger(paste("reading now data for", read_type, "using",
                     structs[struct]), level = 3)
        # Density plot
        den_obj <- paste("density", read_type, struct_types[struct],
                         sep = "_")
        den_rand_obj <- paste(den_obj, "quan", "rand", sep = "_")
        den <- get(den_obj)
        den_quantiles <- get(den_rand_obj)
        res_fig <- get_res_file(
          paste("Fig", "distance", read_type, struct_types[struct],
                "density.pdf", sep = "_"))
        pdf(res_fig, height = 7, width = 7)
        plot(den, ylim = c(0, max(den$y, den_quantiles)),
             xlab = "Normalized distance from domain boundaries",
             main = paste(read_type, struct_types[struct], "( N =", den$n,
                          "Bandwidth =", den$bw, ")"))
        for (quan_area in names(quantiles)) {
          quan_data <- quantiles[[quan_area]]
          if (quan_area == "median") {
            col <- NA
            border <- quan_data[3]
          } else {
            col <- quan_data[3]
            border <- NA
          }
          polygon(c(den$x, rev(den$x)),
                  c(den_quantiles[quan_data[1],],
                    rev(den_quantiles[quan_data[2],])),
                  col = col, border = border, lwd = 2)
        }
        lines(den, lwd = 2)
        abline(v = 0, lty = 2)
        abline(v = 1, lty = 2)
        dev.off()
        log_n_register(res_fig, paste(
          "Kernel density plots of secondary structure features for", read_type,
          struct_types[struct]), lvl = 4)
        # edge plots
        raw_dat <- get(paste("raw", read_type, struct_types[struct], sep = "_"))
        left_prob <- aggregate(V3 ~ V1, subset(raw_dat, V2 < 0), sum)
        right_prob <- aggregate(V3 ~ V2, subset(raw_dat, V1 >= 0), sum)
        for (edge in c("left", "right")) {
          probs <- get(paste(edge, "prob", sep = "_"))
          junc <- get(paste(edge, "junc", read_type, struct_types[struct],
                            sep = "_"))
          offset <- abs(as.integer(colnames(junc)[1])) + 1
          res_fig <- get_res_file(
          paste("Fig", "distance", read_type, struct_types[struct],
                edge, "edge.pdf", sep = "_"))
          pdf(res_fig, height = 7, width = 7)
          boxplot(junc, outline = FALSE, col = "lightgray",
                  ylim = c(0, max(junc, probs$V3, na.rm = TRUE)))
          points(probs[, 1] + offset, probs$V3, lwd = 1.2, cex = 1.4)
          abline(v = offset, lty = 2)
          dev.off()
          log_n_register(res_fig, paste(
            edge, "densities of secondary structure features for", read_type,
            struct_types[struct]), lvl = 4)
        }
      }
    }
  }
  # ---- Data registration ----
  for (obj_name in obj_names) {
    register_data_obj(obj_name, registered_env = data_env,
                      run_ver =  run_version, by_name = TRUE)
  }
  # Garbage collection

  register_module("Distance_analysis", registered_env = data_env,
                  run_ver = run_version)  
}

# ==== Domain density distribution analysis ====================================

if (do_module("Domain_distribution")) {
  pass_dependencies("Domain_distribution", c("Setup"))
  logger("Analysis of density distribution across domains")
  obj_names <- c()
  structs <- list(
    "coiled-coil region" = list(
      list(name = "SignalP", short = "sig",
           trids = gene_list[["signalp_trs"]]$tr_id)))
    ## "Helical; Signal-anchor for type II membrane protein" = list(
    ##   list(name = "SignalP", short = "sig",
    ##        trids = gene_list[["signalp_trs"]]$tr_id)),
    ## "signal peptide" = list(
    ##   list(name = "SignalP", short = "sig",
    ##        trids = gene_list[["signalp_trs"]]$tr_id)),
    ## "Ig-like C2-type" = list(
    ##   list(name = "SignalP", short = "sig",
    ##        trids = gene_list[["signalp_trs"]]$tr_id)),
    ## "EF-hand 1" = list(
    ##   list(name = "SignalP", short = "sig",
    ##        trids = gene_list[["signalp_trs"]]$tr_id)),
    ## "Fibronectin type-III" = list(
    ##   list(name = "SignalP", short = "sig",
    ##        trids = gene_list[["signalp_trs"]]$tr_id)),
    ## "Helical; Name=1" = list(
    ##   list(name = "SignalP", short = "sig",
    ##        trids = gene_list[["signalp_trs"]]$tr_id)),
    ## "Helical; Name=2" = list(
    ##   list(name = "SignalP", short = "sig",
    ##        trids = gene_list[["signalp_trs"]]$tr_id)),
    ## "Helical; Name=3" = list(
    ##   list(name = "SignalP", short = "sig",
    ##        trids = gene_list[["signalp_trs"]]$tr_id)))
    ## "PH" = list(
    ##   list(name = "SignalP", short = "sig",
    ##        trids = gene_list[["signalp_trs"]]$tr_id)))
  bootstrap_n <- get_option("Domain_distribution", "bootstrap-N")
  what_to_do <- get_option("Domain_distribution", "what-to-do")
  up_bound <- get_option("Domain_distribution", "upstream-boundary")
  down_bound <- get_option("Domain_distribution", "downstream-boundary")
  if ("run-scripts" %in% what_to_do) {
    logger("running scripts to generate raw data for domain analysis",
           level = 2)
    dist_prog <- "distribution_analyzer.py"
    dist_ver <-  system(paste(dist_prog, "--version"), intern = T)
    if (package_version(strsplit(dist_ver, " ")[[1]][2]) >= 1) {
      if (all(reps == reps[1])) {
        replicates <- 1:reps[1]
      } else {
        logger("replicate numbers are not consistent, using rep=1",
               level = 3)
        replicates <- 1
      }
      common_args <- c(
        "--mode", "normalized", "--density-db-file",
        "density_files_db_autogen.tsv", "--specs-file",
        get_option("Domain_distribution", "specs"),
        "--treatments", paste(treatments, collapse = ","),
        "--rep", paste(replicates, collapse = ","), 
        "--cds", get_option("Setup", "cds_models"),
        "--keyword", shQuote(paste(names(structs), collapse = "|"), type="sh"),
        "--region", shQuote(paste("custom*-", up_bound, ":custom*+",
                                  down_bound, sep = ""), type="sh"),
        "--uniprot", get_option("Domain_distribution", "uniprot"),
        "--randomized", "--include-tr", "--match-type", "starts",
        "--bootn", bootstrap_n)
      for (read_type in c("Tr", "Mono", "Di")) {
        logger(paste("launching", dist_prog, "for", read_type),
               level = 3)
        args = c("--type", read_type, common_args)
        print(args)
        system2(dist_prog, args = args, wait = TRUE)
        raw <- list(
          "density" = list(),
          "structs" = structs,
          "args" = args
        )
        for (domain_i in 1:length(structs)) {
          file_base <- paste("domain", domain_i - 1, sep = "_")
          domain_name <- names(structs)[domain_i]
          raw$density[[domain_name]] <- list()
          for (i in 0:bootstrap_n) {
            raw$density[[domain_name]][[i+1]] <- read.table(
              file = paste(file_base, i, sep = "_"),
              col.names=c("tr_id", "pos", "density"))
          }
          system2("rm", args = c(paste(file_base, "*", sep = "")), wait = TRUE)
        }
        obj_name <- paste("raw", "domain", read_type, sep = "_")
        assign(obj_name, raw)
        obj_names <- c(obj_names, obj_name)
      }
    } else {
      logger(paste(dist_prog, "version", dist_ver, "is not sufficient"),
             level = 3)
    }
  }
  if ("calc" %in% what_to_do) {
    logger("Analysing raw domain data", level = 2)
    for (read_type in c("Tr", "Mono", "Di")) {
      raw_obj <- paste("raw", "domain", read_type, sep = "_")
      aggr_obj <- paste("aggr", "domain", read_type, sep = "_")
      raw_dat <- get(raw_obj)
      if (exists(aggr_obj)) {
        aggr_dat <- get(aggr_obj)
      } else {
        aggr_dat <- list()
      }
      for (domain_name in names(raw_dat$structs)) {
        logger(paste("processing raw data for", read_type, "using",
                     domain_name), level = 3)
        if (domain_name %in% names(aggr_dat)) {
          logger(paste(domain_name, "is already there for", read_type,
                       ": overwriting it!"), level = 4)
        }
        aggr_dat[[domain_name]] <- list(
          "density" = list(),
          "tr_lists" = list(),
          "args" = raw_dat$args)
        aggr_dat[[domain_name]]$density[["all"]] <- sapply(
          raw_dat$density[[domain_name]], function(x) {
            aggregate(density ~ pos, x, mean, trim=0.02)$density
          })
        raw_trids <- unique(raw_dat$density[[domain_name]][[1]]$tr_id)
        aggr_dat[[domain_name]]$tr_lists[["all"]] <- raw_trids
        for (tr_list in raw_dat$structs[[domain_name]]) {
          tr_list_name <- tr_list[["short"]]
          tr_non_name <- paste("non", tr_list_name, sep = "_")
          tr_ids <- tr_list[["trids"]]
          aggr_dat[[domain_name]]$density[[tr_list_name]] <- sapply(
            raw_dat$density[[domain_name]], function(x) {
              subdata <- subset(x, tr_id %in% tr_ids)
              if (nrow(subdata)) {
                aggregate(density ~ pos, subdata,
                          mean, trim=0.02)$density
              }
            })
          aggr_dat[[domain_name]]$density[[tr_non_name]] <- sapply(
            raw_dat$density[[domain_name]], function(x) {
              subdata <- subset(x, tr_id %ni% tr_ids)
              if (nrow(subdata)) {
                aggregate(density ~ pos, subdata,
                          mean, trim=0.02)$density
              }
            })
          aggr_dat[[domain_name]]$tr_lists[[tr_list_name]] <- raw_trids[
            raw_trids %in% tr_ids] 
          aggr_dat[[domain_name]]$tr_lists[[tr_non_name]] <- raw_trids[
            raw_trids %ni% tr_ids]
        }
      }
      assign(aggr_obj, aggr_dat)
      obj_names <- c(obj_names, aggr_obj)
    }
  }
  if ("clean-rand" %in% what_to_do) {
    logger("cleaning-up large randomization-raw data objects", level = 2)
    big_objs <- paste("raw", "domain", c("Tr", "Mono", "Di"), sep = "_")
    writeLines("Object to be deleted:")
    print(big_objs)
    if (are_you_sure()) {
      for (big_obj in big_objs) {
        rm_data_obj(big_obj, registered_env = data_env,
                    run_ver = run_version, by_name = TRUE)
      }
      logger(
        "large data objects are removed from registered project data",
        level = 3)
    } else {
      logger(
        "skipped by user's choice", level = 3)
    }
  }
  if ("plot" %in% what_to_do) {
    logger("plotting distance densities", level = 2)
    for (read_type in c("Tr", "Mono", "Di")) {
      aggr_dat <- get(paste("aggr", "domain", read_type, sep = "_"))
      for (domain_name in names(aggr_dat)) {
        domain_data <- aggr_dat[[domain_name]]$density
        for (tr_list_name in names(domain_data)) {
          plot_data <- t(domain_data[[tr_list_name]])
          if (nrow(plot_data) == bootstrap_n + 1) {
            res_fig <- get_res_file(
              paste("Fig", "domain", read_type, domain_name,
                    tr_list_name, "density", "distrib.pdf", sep = "_"))
            pdf(res_fig, height = 7, width = 7)
            x_lab <- c(-up_bound:-1, 1:down_bound)
            x_val <- 1:length(x_lab)
            rdata <- data.frame(plot_data[2:(bootstrap_n + 1),])
            boxdata <- boxplot(rdata, names = x_lab, plot = FALSE)
            boxdata$stats <- mod.boxplot.stats(rdata, boxplot.percentiles)$stats
            bxp(boxdata, outline = FALSE, ylim = c(0, 0.01),
                boxfill = "lightgray", lwd = 0.5, xaxt = "n")
            points(plot_data[1,], type = "p", pch = 20, col = "red",
                   lwd = 0.5, cex = 1.2)
            fit <- loess(plot_data[1,] ~ x_val, span = 0.2)
            plx <- predict(fit, se = TRUE) 
            lines(x_val, plx$fit, col = "red")
            lines(x_val, plx$fit - qt(0.975, plx$df) * plx$se,
                  lty = 2, col = "red")
            lines(x_val, plx$fit + qt(0.975, plx$df) * plx$se,
                  lty = 2, col = "red")
            axis(1, at = x_val, labels = FALSE, tcl = -0.25)
            tick_pos <- seq(1, length(x_lab), 5)
            axis(1, at = tick_pos, labels = x_lab[tick_pos])
            tr_n <- length(aggr_dat[[domain_name]]$tr_lists[[tr_list_name]])
            legend("topright",
                   legend=c(paste(domain_name, ", N=", tr_n, sep = "")),
                   pch = 20, col = "red", bty = "n")
            dev.off()
            log_n_register(res_fig, paste(
              "Density distribution plots of", read_type, "for",
              domain_name, "domain using", tr_list_name), lvl = 4)
          }
        }
      }
    }
  }
  # ---- Data registration ----
  for (obj_name in obj_names) {
    register_data_obj(obj_name, registered_env = data_env,
                      run_ver =  run_version, by_name = TRUE)
  }
  # Garbage collection

  register_module("Domain_distribution", registered_env = data_env,
                  run_ver = run_version)  
}


# ==== Save current results ====================================================

logger("Saving results/data of current analysis")
save_cur_data(registered_env = data_env, data_file = project_data_file)
registry_history_file <- get_res_file("Registry_history.tsv")
logger(paste("Updating registry history flat file", registry_history_file))
update_registry_history(registry_history_file,
                        new_env = data_env, old_env = prev_env)
if (run_env_name != "R_GlobalEnv") quit(save = "no", status = 0)
