FUNCTIONS_VERSION <- "2.1.0"

## Required packages -- Please all packages here
my_packages <- c("edgeR", "DESeq2", "ggplot2", "deming", "MASS",
                 "reshape2", "RColorBrewer", "crayon", "sm",
                 "viridis")
dependency <- unlist(
  lapply(my_packages, require, quietly = TRUE, character.only = TRUE))
dependency <- data.frame(is.ok = dependency, row.names = my_packages)

## Functions with no dependency

# --- Object registration, saving and loading ---
load_prev_data <- function ( registered_env, data_file ) {
  if (file.exists(data_file)) {
    load(data_file, envir = registered_env, verbose = FALSE)
  } else {
    registered_env$registrar <- list()
  }
}

save_cur_data <- function ( registered_env, data_file ) {
  save(list = ls(registered_env), file = data_file, envir = registered_env)
}

register_data_obj <- function ( ..., registered_env, run_ver,
                                code_ver = code_version, by_name = FALSE) {
  t <- date()
  if (by_name) {
    obj_names <- list(...)
    obj_list <- Map(get, obj_names)
  } else {
    obj_list <- list(...)
    obj_names <- as.list(as.character(match.call()))[-1]
    obj_names <- obj_names[1:length(obj_list)]
  }
  Map(function(...) assign(..., envir = registered_env), obj_names, obj_list)
  Map(function(...) registered_env$registrar[[...]] <- list(
    "time" = t, "code_ver" = code_ver, "run_ver" = run_ver), obj_names)
  invisible(NULL)
}

rm_data_obj <- function ( ..., registered_env, run_ver,
                         code_ver = code_version, by_name = FALSE) {
  t <- date()
  if (by_name) {
    obj_names <- list(...)
  } else {
    obj_names <- as.list(as.character(match.call()))[-1]
    obj_names <- obj_names[1:length(list(...))]
  }
  rm(list = as.character(obj_names), envir = registered_env)
  Map(function(...) registered_env$registrar[[...]] <- list(
    "time" = t, "code_ver" = code_ver, "run_ver" = run_ver,
    "descr" = "removed"), obj_names)
  invisible(NULL)  
}

register_result_file <- function(file_name, descr, registered_env, run_ver,
                                 code_ver=code_version) {
  t <- date()
  registered_env$registrar$result_files[[file_name]] <- list(
    "time" = t, "code_ver" = code_ver, "run_ver" = run_ver, "descr" = descr)
  invisible(NULL)
}

register_module <- function(module_name, registered_env, run_ver,
                            code_ver = code_version) {
  t <- timestamp(quiet = TRUE, prefix = "", suffix = "")
  registered_env$registrar$modules[[module_name]] <- list(
    "time" = t, "code_ver" = code_ver, "run_ver" = run_ver, "descr" = "OK")
  invisible(NULL)
}

get_registrar <- function(registered_env) {
  registrar_table <- data.frame()
  for (data_name in names(registered_env$registrar)) {
    if (data_name %in% c("result_files", "modules")) next
    data_obj <- registered_env$registrar[[data_name]]
    registrar_table <- rbind(registrar_table,
      cbind("name" = data_name, "timestamp" = data_obj$time,
            "description" = if ("descr" %in% names(data_obj)) data_obj$descr
              else "-",
            "code_version" = paste(data_obj$code_ver, collapse = "_"),
            "run_version" = data_obj$run_ver))
  }
  registrar_table[] <- lapply(registrar_table, as.character)
  registrar_table
}

copy_registrar <- function(from_env, to_env) {
  to_env$registrar <- from_env$registrar
}

get_registered_obj <- function(registered_env, collection) {
  obj_table <- data.frame()
  obj_collection <- registered_env$registrar[[collection]]
  for (obj in names(obj_collection)) {
    data_obj <- obj_collection[[obj]]
    descr <- if (! is.null(data_obj$descr)) data_obj$descr else "-"
    code_ver <- paste(data_obj$code_ver, collapse = "_")
    obj_table <- rbind(obj_table,
      cbind("name" = obj, "timestamp" = data_obj$time,
            "description" = descr, "code_version" = code_ver,
            "run_version" = data_obj$run_ver))
  }
  obj_table[] <- lapply(obj_table, as.character)
  obj_table
}

get_registered_modules <- function(registered_env) {
  get_registered_obj(registered_env, collection = "modules")
}

get_registered_files <- function(registered_env) {
  get_registered_obj(registered_env, collection = "result_files")
}

update_registry_history <- function(history_file, new_env, old_env) {
  get_diff <- function(new_table, old_table) {
    diff_table <- data.frame()
    for (new_index in 1:nrow(new_table)) {
      obj <- new_table$name[new_index]
      old_index <- which(old_table$name == obj)
      if (length(old_index) == 0 |
          ! all(as.character(old_table[old_index, ]) ==
            as.character(new_table[new_index, ]))) {
        diff_table <- rbind(diff_table, new_table[new_index, ])
      }
    }
    diff_table
  }
  update_struct <- list(
    list(command = get_registered_files, label = "Result_file"),
    list(command = get_registrar, label = "Data_object"),
    list(command = get_registered_modules, label = "Module_execution"))
  updates <- data.frame()
  for (section in update_struct) {
    res <- get_diff(section$command(new_env), section$command(old_env))
    if (length(res)) updates <- rbind(updates, cbind(type = section$label, res))
  }
  date_str <- as.POSIXct(updates$timestamp,
                         format = "%a %b %d %H:%M:%S %Y", tz = "UTC")
  updates <- updates[order(date_str), ]
  write.table(updates, file = history_file, append = TRUE, sep = "\t",
              quote = FALSE, col.names = FALSE, row.names = FALSE)
}

check_dependency <- function(module_list, registered_env) {
  run_modules <- get_registered_modules(registered_env)$name
  all(module_list %in% run_modules)
}

# --- Count filtering ---
filter_by_count_number <- function (one.row, min.count=10, min.fraction=0.25) {
  length(which(one.row > min.count)) / length(one.row) > min.fraction
}

get_count_filter <- function (countTable, min.count = 10, min.fraction = 1 / 6,
                              logical = FALSE, reverse = FALSE) {
  f <- apply(countTable, 1, filter_by_count_number, min.count, min.fraction)
  if (reverse) f <- !f
  if (logical) return(f) else return(which(f))
}

# --- Utility functions ---
no.inf <- function(d.f, repl = NaN) {
  apply(d.f, c(1, 2), function(x) if (is.infinite(x)) repl else x )
}

"%ni%" <- Negate("%in%")
"%ini%" <- function(x, y) !is.na(x)
"]" <- function(x, y) x[length(x) + 1 - y]

make_file_name <- function(s) {
  t <- gsub(" ", "_", s)
  t <- gsub("/", "-", t)
  t
}

random_str <- function ( length=12 ) {
  paste(sample(c(0:9, letters, LETTERS), length, replace = TRUE),
        collapse = "")
}

get_par_limits <- function(p=4) {
  d <- (100 + 2 * p) / p
  usr <- par("usr")
  xr <- (usr[2] - usr[1]) / d
  yr <- (usr[4] - usr[3]) / d
  xlim <- c(usr[1] + xr, usr[2] - xr)
  ylim <- c(usr[3] + yr, usr[4] - yr)
  list(xlim = xlim, ylim = ylim)
}

log_scatter <- function (x, y, steps = c(1, 1), log.ax = "xy", ...) {
  needx <- grepl("x", log.ax)
  needy <- grepl("y", log.ax)
  if (needx) xaxt <- "n" else xaxt <- "s"
  if (needy) yaxt <- "n" else yaxt <- "s"
  plot(x, y, log = log.ax, xaxt = xaxt, yaxt = yaxt, ...)
  usr <- par("usr")
  if (needx) {
    xr <- (usr[2] - usr[1]) / 27 # 27 = (100 + 2*4) / 4
    xlim <- c(usr[1] + xr, usr[2] - xr)
    atx <- log10(axTicks(1))
    atx <- seq(min(atx), floor(xlim[2]), steps[1])
    labelsx <- sapply(atx, function(x) as.expression(bquote(10 ^ .(x))))
    axis(1, at = 10 ^ atx, labels = labelsx)
  }
  if (needy) {
    yr <- (usr[4] - usr[3]) / 27
    ylim <- c(usr[3] + yr, usr[4] - yr)
    aty <- log10(axTicks(2))
    aty <- seq(min(aty), floor(ylim[2]), steps[2])
    labelsy <- sapply(aty, function(x) as.expression(bquote(10^ .(x))))
    axis(2, at = 10 ^ aty, labels = labelsy)
  }
}

lmp <- function (modelobject, label = "") {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  s <- summary(modelobject)
  f <- s$fstatistic
  f.p <- pf(f[1], f[2], f[3], lower.tail=FALSE)
  t <- summary(modelobject)$coefficients
  t.p <- t[2, 4]
  t.t <- t[2, 3]
  b <- t[2, 1]
  leg.text <- bquote(paste(
    .(label), ": ", italic(b), "=", .(round(b, 3)), ", ", italic(t), "(",
    .(f[3]), ")=", .(round(abs(t.t), 2)), ", ", italic(p), "=",
    .(signif(t.p, 3)), "; ", italic(R) ^ 2, "=", .(round(s$r.squared, 3)),
    ", ", italic(F), "(", .(f[2]), ",", .(f[3]), ")=", .(round(f[1], 2)),
    ", ", italic(p), "=", .(signif(f.p, 3))))
  return(leg.text)
}

pretty_grid <- function (distx, disty=distx, col="lightgray",
                         lty="dotted", lwd=par("lwd")) {
  usr <- floor(par("usr"))
  atx <- c(usr[1]:usr[2])[usr[1]:usr[2] %% distx == 0]
  aty <- c(usr[3]:usr[4])[usr[3]:usr[4] %% disty == 0]
  abline(v = atx, col = col, lty = lty, lwd = lwd)
  abline(h = aty, col = col, lty = lty, lwd = lwd)
}

logger <- function(msg, level=1, is.message=F) {
  lvl_pads <- c(" ", "   ... ", "     .. ", "       .. ")
  lvl_cols <- list(green, blue, cyan, yellow)
  level <- if (level < 1) 1
           else if (level > length(lvl_pads)) length(lvl_pads)
           else level
  msg <- paste(lvl_pads[level], msg, sep = "")
  log_color <- lvl_cols[[level]]
  msg_color <- red
  if (is.message) message(msg_color(msg)) else writeLines(log_color(msg))
}

# --- Read read-counts from text files ---
# Low-level function to read from count text file
read_countTable <- function (filename, nrows, col.names) {
    read.table(textConnection(readLines(filename)[1:nrows]),
               stringsAsFactors = FALSE, sep = "\t",
               col.names = col.names)
}

# High-level function to make a data frame from
# supplied samples and fields (cds, utr etc)
read_counts <- function (samples, type = c("cds"), nrows,
                         col.names = c("GeneID", "x5utr", "cds", "x3utr",
                                       "utr", "not_in_db")) {
    filenames <- with(samples, paste(path, name, ext, sep = ""))
    first_tb <- read_countTable(filenames[1], nrows, col.names)
    if (any(c("all", "All", "a", "A") %in% type)) {
        read_cols <- which(!colnames(first_tb) %in% "GeneID")
    } else {
        read_cols <- which(colnames(first_tb) %in% setdiff(type, "GeneID"))
    }
    datafr <- data.frame(first_tb[c(1, read_cols)])
    for (i in 2:nrow(samples)) {
        datafr <- cbind(datafr, read_countTable(filenames[i], nrows,
                                                col.names)[read_cols])
    }
    row.names(datafr) <- datafr$GeneID
    datafr <- datafr[-1]
    colnames(datafr) <- paste(rep(samples$sample, each = length(read_cols)),
                              colnames(first_tb)[read_cols], sep = ".")
    datafr
}

## Functions with dependencies

getPCA <- function(x, intgroup, ntop=500) {
  require("IRanges")
  require("genefilter")
  rv <- rowVars(assay(x))
  select <- order(rv, decreasing = TRUE)[seq_len(ntop)]
  prcomp(t(assay(x)[select, ]))
}

# --- Differential expression analysis ---
if (dependency["DESeq2", "is.ok"]) {
  get_cds <- function (counts, colData, design, normfactors, ...) {
      cds <- DESeqDataSetFromMatrix(counts, colData = colData, design = design)
      sizeFactors(cds) <- normfactors
      cds <- estimateDispersions(cds, ...)
      cds
  }
  get_dispersions <- function(cds) 1 / dispersions(cds)
} else {
  get_cds <- function (...) {
    stop("'get_cds' function depends on DESeq2 package")
  }
  get_dispersions <- function (...) {
    stop("'get_diepersions' function depends on DESeq2 package")
  }
}

# --- Normalization functions ---
if (dependency["edgeR", "is.ok"]) {
    get_norm_factors <- function(counts, group, what = "normfactor",
                                 method="upperquartile", ...) {
        dge_counts <- DGEList(counts = counts, group = group)
        dge_counts <- calcNormFactors(dge_counts, method = method, ...)
        eff.libsize <- dge_counts$samples$lib.size *
                       dge_counts$samples$norm.factors
        eff.libgeo <- exp(mean(log(eff.libsize)))
        if (what == "normfactor") return(eff.libsize / eff.libgeo)
        if (what == "efflibsize") return(eff.libgeo)
        if (what == "dgenorm") return(dge_counts$samples$norm.factors)
    }
} else {
    get_norm_factors <- function (...) {
        stop("'get_norm_factors' function depends on edgeR package")
    }
}

# --- High level plotting functions ---
scatter_with_marginal_densities <- function(
    d.f, ordered_keys, my.cols, margin.xy="xy", multiplier = 0.8,
    multiplier.x = multiplier, multiplier.y = multiplier, fg.cols = NA,
    den.cols = NA, lwd.den = NA, main = "", xlab="", ylab="",
    label.it = FALSE, xlim=c(), ylim=c(), log.xy = "xy", den.solid = TRUE,
    pch=21, verbose = FALSE, ...) {
  dots <- list(...)
  nmdots <- names(dots)
  for (ax in c("x", "y")) {
    if (grepl(ax, margin.xy)) {
      density.m <- list()
      i <- 1
      for (key in ordered_keys) {
        if (key == "All" | key == "all") {
          if (verbose) writeLines(
            "Special key: *all* is used instead of key matching")
          density.m[[i]] <- density(get(ax, d.f), na.rm = TRUE)
        } else {
          density.m[[i]] <- density(subset(d.f, k == key, get(ax))[, 1],
                                    na.rm = TRUE)
        }
        i <- i + 1
      }
      assign(paste("density", ax, sep = "."), density.m)
    }
  }
  if (length(xlim) == 0) {
    xlim <- c(min(d.f$x), max(d.f$x))
  }
  if (length(ylim) == 0) {
    ylim <- c(min(d.f$y), max(d.f$y))
  }
  # canvas
  if (grepl("x|y", log.xy)) {
    log_scatter(1, 1, log.ax = log.xy, ylim = ylim, xlim = xlim, type = "n",
               main = main, xlab = xlab, ylab = ylab, ...)
    grid("grey")
  } else {
    eqscplot(1, 1, ylim = ylim, xlim = xlim, type = "n", main = main,
             xlab = xlab, ylab = ylab, ...)
    pretty_grid(2, col = "grey")
  }
  real_limits <- get_par_limits(p = 2)
  xlim <- real_limits$xlim
  ylim <- real_limits$ylim
  # scatterplots
  i <- 1
  if (pch < 21 | pch > 25) fg.cols <- my.cols
  for (key in ordered_keys) {
    points(subset(d.f, k == key, 1:2), type = "p", col = fg.cols[i],
           bg = my.cols[i], pch = pch, ...)
    i <- i + 1
  }
  # marginal densities
  if (length(den.cols) == 1 && is.na(den.cols)) den.cols <- my.cols
  if (den.solid) {
    den.cols <- paste(substring(den.cols, 1, 7), "FF", sep = "") 
  } else {
    den.cols <- if (length(den.cols[0]) > 7) den.cols else {
      paste(substring(den.cols, 1, 7), "BF", sep="")
    }
  }
  if (is.na(lwd.den)) {
    lwd.den <- if ("lwd" %in% nmdots) dots$lwd
    else par("lwd")
  }
  if (grepl("x", margin.xy)) {
    divider <- max(sapply(density.x, function(x) max(x$y)))
    for (i in 1:length(ordered_keys)) {
      points(density.x[[i]]$x,
             density.x[[i]]$y / divider * multiplier.x + ylim[1],
             type = "l", col = den.cols[i], lwd = lwd.den, ...)
    }
  }
  if (grepl("y", margin.xy)) {
    divider <- max(sapply(density.y, function(x) max(x$y)))
    for (i in 1:length(ordered_keys)) {
      points(density.y[[i]]$y / divider * multiplier.y + xlim[1],
             density.y[[i]]$x,
             type = "l", col = den.cols[i], lwd = lwd.den, ...)
    }
  }
  if (label.it) {
    d.f.sub <- subset(d.f, name != "")
    text(d.f.sub$x, d.f.sub$y, d.f.sub$name,
         pos = d.f.sub$pos, offset = 0.25, cex = 0.6)
  }
}


## BiomaRt Functions
get_id_from_biomart <- function(gid, what_to_get = c("entrezgene",
                                                     "external_gene_id"),
                                version = 75) {
  if (require(biomaRt)) {
    ens_mart <- useEnsembl(biomart = "ensembl",
                           dataset = "mmusculus_gene_ensembl",
                           version = version)
    getBM(attributes = c("ensembl_gene_id", what_to_get),
          filters = "ensembl_gene_id", values = gid, mart = ens_mart)
  } else {
    print("BiomaRt package could not be loaded..")
    FALSE
  }
}

# Kegg etc pathway analysis related Functions
subs_ids_geneset <- function(gid, id_db, geneset_ids, prefix="",
                             ensembl_col = "ensembl_gene_id",
                             target_col = "entrezgene") {
  found <- which(id_db[ensembl_col] == gid)
  len_found <- length(found)
  if (len_found == 0) return(paste(prefix, "NA", sep = ""))
  best_case <- paste(prefix, id_db[found[1], target_col], sep = "")
  if (len_found == 1) return(best_case)
  for (i in found) {
    eid <- paste(prefix, id_db[i, target_col], sep = "")
    if (eid %in% geneset_ids) return(eid)
  }
  best_case
}

extract_regions <- function (df, from_the_beginning=T, size=222) {
  f <- df[order(df$V1), 3]
  sizes <- table(df$V1)
  slc <- c()
  gids <- c()
  cpos <- 0
  for (gid in levels(df$V1)) {
    csize <- sizes[gid]
    if (csize >= size) {
      gids <- c(gids, gid)
      cur <- if (from_the_beginning) cpos + 1:size else { 
        cpos + (csize - size + 1):csize
      }
      slc <- c(slc, cur)
    }
    cpos <- cpos + csize
  }
  return(cbind(data.frame(V1 = rep(gids, each = size),
                          V2 = rep(1:size, length(gids))), V3 = f[slc]))
}

mod.boxplot.stats <- function(x, func, ...) {
  args <- list(x)
  namedargs <- if (!is.null(attributes(args)$names))
    attributes(args)$names != ""
  else rep_len(FALSE, length(args))
  groups <- if (is.list(x))
    x
  else args[!namedargs]
  if (0L == (n <- length(groups)))
    stop("invalid first argument")
  if (length(class(groups)))
    groups <- unclass(groups)
  cls <- vapply(groups, function(x) class(x)[1L], "")
  cl <- if (all(cls == cls[1L]))
          cls[1L]
  for (i in 1L:n) groups[i] <- list(
    func(unclass(groups[[i]]), ...))
  stats <- matrix(0, nrow = 5L, ncol = n)
  ct <- 1
  for (i in groups) {
    stats[, ct] <- i$stats
    ct <- ct + 1
  }
  list(stats = stats)
}

boxplot.percentiles <- function (x, qprob=c(.05, .25, .5, .75, .95)) {
  list(stats = quantile(x, qprob, na.rm = TRUE))
} 

rotate_points <- function(x, y, point=c(0, 0), theta=pi/3) {
  v <- matrix(c(x, y, rep(1, length(x))), ncol=3)
  # shift to center
  locator1 <- matrix(c(1, 0, -point[1],
                       0, 1, -point[2],
                       0, 0, 1), ncol=3, byrow=T)
  # shift back to original
  locator2 <- matrix(c(1, 0, point[1],
                       0, 1, point[2],
                       0, 0, 1), ncol=3, byrow=T)
  # counter-clockwise
  rotator <- matrix(c(cos(theta), -sin(theta), 0,
                      sin(theta), cos(theta), 0,
                      0, 0, 1), ncol=3, byrow=T)

  transformer <- locator1 %*% rotator %*% locator2
  v_rotated <- v %*% transformer
  v_rotated[,1:2]
}
