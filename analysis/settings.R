SETTINGS_VERSION <- "3.2.2"

project_data_file <- "disome_registered_data_obj.RData"
project_spec_file <- "../config/disome_specs.yaml"

## MODULES

modules <- list(
    "Setup" = FALSE,  # OK
    "Mapping" = FALSE,  # OK
    "Insert_size" = FALSE, # OK
    "Gene_lists" = FALSE, # OK
    "Trinucleotide" = FALSE, # OK
    "Count" = FALSE, # OK
    "Filter" = FALSE, # OK
    "Normalize" = FALSE, # OK
    "Diff_exp" = FALSE, # OK
    "Diff_TE" = FALSE, # OK
    "Summary" = FALSE, # OK
    "Percent" = FALSE, # OK
    "Di_vs_mono" = FALSE, # OK
    "Treatment_plots" = FALSE, # OK
    "PCA" = FALSE, # OK
    "Distance_analysis" = FALSE, # OK
    "Domain_distribution" = FALSE, # OK
    "Signalp_distribution" = FALSE, # OK
    "GSEA-gage" = TRUE
)

# options
module_options <- list(

  "Setup" = list(
    "gids" = "../cds_preparation/disome_filtered_genes.tsv",
    "count_files_ext" = "_count_files.tsv",
    "count_cols" = c("sample", "treatment", "rep",
                     "path", "name", "ext", "use"),
    "cds_models" = "../cds_preparation/disome_prepared_cds.tsv",
    "cds_cols" = c("gene_id", "status", "tr_id", "tr_type", "prot_id",
                   "tsl", "start_ok", "end_ok", "size", "cds_start",
                   "cds_end", "flag"),
    "annotation" = "../cds_preparation/disome_gene_annotation.tsv",
    "annot_cols" = c("gene_id", "gene_name", "gene_type",
                     "tr_id", "tr_name", "tr_type")
  ),

  "Count" = list(
    "tables" = list(
      list(path = NULL, label = "", tag = ""),
      list(path = "../count_pri_data/", label = "is_primary", tag = "pri"),
      list(path = "../count_first75_data/", label = "first_75", tag = "f75")
    )
  ),

  "Normalize" = list(
    "count-labels" = c("original", "pri"),
    "print-density-db" = TRUE,
    "density_folder" = "../density_data",
    "density_extension" = "_5prime_count.bgz"
  ),

  "Gene_lists" = list(
    "folder" = "../extra_data",
    "export" = TRUE,
    "entrez" = TRUE
  ),

  "GSEA-gage" = list(
    "folder" = "../extra_data",
    "msig-db-file" = "MSigDB_ver_6.2.1_mmu.RData",
    "merge-treatments"= TRUE,
    "count-label" = "is_primary", 
    "KEGG" = TRUE,
    "GO" = FALSE,
    "MSigDB" = TRUE,
    "rotate" = TRUE
  ),

  "Mapping" = list(
    "map_data_file" = "disome_mapping_data.tsv",
    "map_data_cols" = c("sample", "map_type", "category", "count")
  ),

  "Trinucleotide" = list(
    "run-scripts" = TRUE,
    "specs" = project_spec_file,
    "merge-treatments" = TRUE
  ),

  "Distance_analysis" = list(
    "what-to-do" = c("plot"), #c("clean-rand", "run-scripts", "data", "randomized", "calc", "plot"),
    "specs" = project_spec_file,
    "fasta" = "../rust_analysis/structure/structure.fa",
    "bootstrap-N" = 10000
  ),

  "Domain_distribution" = list(
    "what-to-do" = c("plot", "calc", "run-scripts"), #"clean-rand"), # calc
    "upstream-boundary" = 15,
    "downstream-boundary" = 90,
    "specs" = project_spec_file,
    "uniprot" = "/local/databases/uniprot/uniprot_mouse_features_04_2019.v2.json",
    "bootstrap-N" = 1250
   ),

  "Diff_TE" = list(
    "xtail" = TRUE,
    "xtail_bins" = 1000,
    "riborex" = FALSE
  ),

  "Di_vs_mono" = list(
    "count-label" = "is_primary",
    "zt-density" = FALSE,
    "conservation" = TRUE,
    "merge-treatments" = TRUE
  )

)

# functions
do_module <- function (module_name) isTRUE(modules[[module_name]])

get_option <- function(module_name, option_name) {
   module_options[[module_name]][[option_name]]
}

is_option_set <- function(module_name, option_name) {
  isTRUE(get_option(module_name, option_name))
}

## DESIGN

require(yaml)
proj_specs <- yaml.load_file(project_spec_file)
treatment_order <- proj_specs$specs$treatment_levels
num_reps <- proj_specs$specs$num_reps

## COLORS

# Main colors
require(graphics)
get_analogous_colors <- function (color, d=15) {
  hsv_col <- rgb2hsv(col2rgb(color))
  rotations <- c(-d / 360, d / 360)
  h <- (hsv_col[1] + rotations) %% 1
  c(hsv(h[1], hsv_col[2], hsv_col[3]), hsv(h[2], hsv_col[2], hsv_col[3]))
}

require("RColorBrewer")
color_n <- 7 # Need odd number!
color_palette <- brewer.pal(n = color_n, "BrBG")

mono_col <- color_palette[color_n]
di_col <- color_palette[1]
tr_col <- color_palette[(1 + color_n) / 2]

# Analogs
mono_col_analogs <- get_analogous_colors(mono_col, 30)
di_col_analogs <- get_analogous_colors(di_col, 30)

# Variants
mono_col_var <- color_palette[color_n - 1]
di_col_var <- color_palette[2]

# Color functions
mono_col_func <- colorRampPalette(c(mono_col_analogs[1],
                                    mono_col_var, mono_col_analogs[2]))
di_col_func  <- colorRampPalette(c(di_col_analogs[1],
                                   di_col_var, di_col_analogs[2]))
