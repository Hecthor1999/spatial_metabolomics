# This R script calculates the wassersteindistance by using the mwass.cpp file in /storage/homefs/ha25g949/pannet_metabolism/mwass.cpp
# then it calcualtes the batch balanced k nearest neighbour graph using bbknnR, runs Leiden clustering on the graph and then visualizes both the full image and the PCOAs for all images in either negative or positive modes
# the input of this script is an imzml file with an ibd and a metadata file

# arguments
# -1 list of images to analyze
# -2 output directory
# -3 input directory

# Set library directory
.libPaths("~/R/library")
my_paths <- c(.libPaths()[1],"/storage/software/epyc2.9/software/R-bundle-CRAN/2024.11-foss-2024a",.libPaths()[2])
.libPaths(my_paths)

# Configure BioCManager to use Posit Package Manager:
options(BioC_mirror = "https://packagemanager.posit.co/bioconductor/latest")
options(BIOCONDUCTOR_CONFIG_FILE = "https://packagemanager.posit.co/bioconductor/latest/config.yaml")

# Set the Bioconductor version to prevent defaulting to a newer version:
Sys.setenv("R_BIOC_VERSION" = "3.20")

# Configure a CRAN snapshot compatible with Bioconductor 3.20:
options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/rhel9/latest"))

# install and load packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", lib="~/R/library")

library(BiocManager)

BiocManager::install("scran", type = "source", ask = FALSE, force = TRUE)
BiocManager::install("scater", type = "source")

if (!require("dplyr", quietly = TRUE))
    install.packages("dplyr", lib="~/R/library")

#if (!require("Cardinal", quietly = TRUE))
#BiocManager::install("Cardinal", lib="~/R/library")

if (!require("BiocParallel", quietly = TRUE))
BiocManager::install("BiocParallel", lib="~/R/library", force = TRUE)

BiocManager::install("leidenbase", type = "source", update = FALSE)
BiocManager::install("igraph", type = "source", update = FALSE)
BiocManager::install("beachmat",lib="~/R/library")
BiocManager::install("Cardinal", lib="~/R/library", update = FALSE)

install.packages("tidyverse")
install.packages("openxlsx")
install.packages("purr")
install.packages("ggplot2")
install.packages("ape")
install.packages("bbknnR")

library(tidyverse)
library(scran) 
library(scater)
library(openxlsx)
library(purrr)
library(ggplot2)
library(ape)
library(cluster)
library(Cardinal)
library(SingleCellExperiment)
library(Matrix)
library(dplyr)

# command line arguments
args <- commandArgs(trailingOnly = TRUE)

# data directory
data_dir <- args[3]

# the meta data is required to match the image files from the different modes to the correct sample
# This list is mainly hand curated 
DESI_meta_data <- read.xlsx(
  "/storage/homefs/ha25g949/pannet_metabolism/240208_sample_list_DESI_RNA_coregistration.xlsx",
  sheet = "sample_overview"
) %>% mutate(patient_uid = paste0(aPnumber, "_", B_extended))

# remove the first numbers to match with the IMZML files
DESI_meta_data$sample <- gsub("^[^_]*_", "", DESI_meta_data$sample)

# There is one sample that with a slightly different result that needs to be checked
#IMZML_1k_files = data.frame(top_1k = list.files(data_dir, pattern = "CLMC.txt", recursive = T)) %>% 
#  filter(!grepl("ambiguous", top_1k)) %>% 
#  mutate(sample = gsub(" Analyte.*$", "", basename(top_1k)))

IMZML_1k_files <- data.frame(
  top_1k = list.files(data_dir, pattern = "\\.imzML$", recursive = TRUE, full.names = FALSE)
) %>%
  filter(!grepl("ambiguous", top_1k)) %>%
  mutate(sample = tools::file_path_sans_ext(basename(top_1k))) %>%
  mutate(sample = gsub("^\\d+_", "", sample)) 

# remove bottom to match with DESI meta data
IMZML_1k_files$sample <- gsub("[[:alpha:]]_Analyte_1_CLMC_1", "N", IMZML_1k_files$sample)


DESI_meta_data_joined <- left_join(DESI_meta_data, IMZML_1k_files, by = "sample") %>%
filter(!is.na(top_1k))

rm(IMZML_1k_files)

processMetaData = function(filepath) {
  step_df = data.frame(axis = c("X", "Y"), 
                       step = c(NA, NA))
  
  con = file(filepath, "r")
  while ( TRUE ) {
    line = readLines(con, n = 1)
    
    if (length(line) != 0){
      if (grepl("[XY]Step", line)){
        axis = str_extract(line, pattern = "[XY]")
        value = str_extract(line, pattern = "[0-9\\.]+$")
        step_df$step[step_df$axis == axis] = value
      } else{
        next
      }
      
      sum(is.na(step_df$step))
      if (sum(is.na(step_df$step)) == 0){
        close(con)
        return(step_df)
      }
    }
  }
}

# --- IMPORTANT ---
# When loading the data edge pixels are always removed 
# These pixels frequently contain artifacts 
# -----------------
getHDIPeaks = function(image_name,
                       remove_edge = T){
  # This is the optimized data loading function
  # The old readr based function used in 230619_PanNET_DESI_coregistration is no longer supported
  # first the colnames are extracted to fix the missing values in the front and at the end
  # In between the colnames are converted to numeric and trailing zeroes are dropped - this should not be problematic
  raw_colnames = data.table::fread(file.path(data_dir, 
                                             image_name),
                                   sep = "\t",
                                   skip = 2, 
                                   nrows = 1) %>% 
    unlist()
  raw_colnames[1:3] = c("index", "x", "y")
  # The last two columns of the data have not names and contain no useful information
  raw_colnames = c(raw_colnames, c("to_remove_1", "to_remove_2"))
  n_cols = length(raw_colnames)
  # After loading the data the data frame is reordered and values are rounded to integer
  hdi_image = data.table::fread(file.path(data_dir, 
                                          image_name),
                                sep = "\t",
                                skip = 3,
                                col.names = raw_colnames) 
  # Selecting and sorting the peaks 
  data_mat = as.matrix(hdi_image[, 4:(n_cols - 2)])
  data_peaks = as.numeric(colnames(data_mat))
  data_mat = data_mat[, order(data_peaks)]
  data_peaks = sort(data_peaks)
  # Generic colnames should simplify addressing peaks
  colnames(data_mat) = paste0("p", seq_len(ncol(data_mat)))
  # The data is rounded to integer mainly to save on data size 
  data_mat = round(data_mat, 0)
  mode(data_mat) = "integer"
  
  # Making the sample meta data 
  colData = hdi_image[,c(1:3)] %>% 
    as.data.frame() %>% 
    mutate(x_orig = match(x, sort(unique(x))),
           x_new = x_orig - 1,
           y_orig = match(y, sort(unique(y))),
           y_new = y_orig -1,
           is_border = x_orig == 1 | x_orig == max(x_orig) | y_orig == 1 | y_orig == max(y_orig))
  
  # extracting general image information
  image_res = processMetaData(file.path(data_dir, 
                                        gsub("/imaging.*", "", image_name), 
                                        "_extern.inf")) %>% 
    mutate(step = as.integer(as.numeric(step) * 1e3))
  
  image_metadata = list(sample = gsub(" Analyte.*$", "", basename(image_name)),
                        x_res = image_res$step[1],
                        y_res = image_res$step[2])
  
  data_object = SingleCellExperiment(
    assays = list(counts = t(data_mat)),
    colData = DataFrame(colData),
    rowData = DataFrame(peak_index = seq_len(ncol(data_mat)),
                        mz = data_peaks),
    metadata = image_metadata)
  
  # Regions with very low signal are removed - these represent scan errors
  # The threshold is chosen empirically 
  failed_pixel = colSums(assay(data_object, "counts")) < 5e5
  metadata(data_object)$failed_pixel = sum(failed_pixel)
  if (sum(failed_pixel) > 0)
    data_object = data_object[, !failed_pixel]
  if (remove_edge)
    data_object = data_object[, !data_object$is_border]
  gc()
  return(data_object)
}

getImzMLPeaks <- function(imzml_path, remove_edge = TRUE, bin_size = 0.01, peak_not_aligned = TRUE) {

  # Load imzML data
  img <- readMSIData(imzml_path)

  # Bin spectra to fixed m/z intervals
  #img_binned <- bin(img_binned, resolution = bin_size, units = "mz")

  # Get intensity matrix (pixels × m/z)
  if (!peak_not_aligned ) {
  intensity_mat <- spectra(img, "intensity")  # Each row = pixel, columns = m/z bins
  mz_vals <- fData(img)$mz
  }

  # if using peak picking with a referece
  if ( peak_not_aligned ){ 
  intensity_mat <- do.call(rbind, as.list(spectra(img, "intensity")))
  mz_list <- do.call(rbind, as.list(spectra(img, "mz")))
  mz_vals <- mz_list[1, ]
  }

  # Convert to dense matrix if needed
  #if (inherits(intensity_mat, "sparseMatrix")) {
  #  intensity_mat <- as.matrix(intensity_mat)
  #}

  # Spatial coordinates
  coords <- coord(img)
  colnames(coords) <- c("x", "y")

  # Only keep intensity for these pixels
  # intensity_mat <- intensity_mat[seq_len(nrow(coords)), , drop = FALSE]

  # Remove edge pixels if requested
  if (! peak_not_aligned ) {
  if (remove_edge) {
  x_ord <- rank(coords[, "x"])
  y_ord <- rank(coords[, "y"])
  is_border <- x_ord == 1 | x_ord == max(x_ord) | y_ord == 1 | y_ord == max(y_ord)
  intensity_mat <- intensity_mat[, !is_border, drop = FALSE]  # columns = pixels
  coords <- coords[!is_border, , drop = FALSE]
    }
  }


  # Remove low-signal pixels
  #total_signal <- Matrix::rowSums(intensity_mat)   # works with sparse matrices
#
  ## Only keep pixels that exist in coords and pass signal threshold
  #keep <- which(total_signal > 5e5)                # indices, not logical
  #intensity_mat <- intensity_mat[keep, , drop = FALSE]
  #coords <- coords[keep, , drop = FALSE]

  # DONT USE THIS FEATURE WHEN USING PEAK PICKED DATA Transpose so features (peaks) are rows and pixels are columns 
  if ( peak_not_aligned ){ 
  intensity_mat <- t(intensity_mat)
  }

  # Create SingleCellExperiment object
  sce <- SingleCellExperiment(
  assays = list(counts = intensity_mat),
  colData = DataFrame(coords),
  rowData = DataFrame(mz = mz_vals),
  metadata = list(
    file = basename(imzml_path),
    x_res = if (nrow(coords) > 1) mean(diff(sort(unique(coords[, "x"])))) else NA,
    y_res = if (nrow(coords) > 1) mean(diff(sort(unique(coords[, "y"])))) else NA
  )
)
  return(sce)
}

# --- SAFE VERSION OF LOADER ---
#
#neg_rows <- DESI_meta_data %>% filter(DESI_mode == "negative")
#cat("Found", nrow(neg_rows), "negative mode files to load\n")
#HDI_neg_1k <- list()
#
#for (i in seq_len(nrow(neg_rows))) {
#  sample_name <- neg_rows$sample[i]
#  f <- neg_rows$top_1k[i]
#
#  cat("→ [", i, "/", nrow(neg_rows), "] Loading ", sample_name, " ...\n", sep = "")
#
#  res <- tryCatch({
#    getHDIPeaks(f)
#  }, error = function(e) {
#    cat("  ⚠️ Error: ", e$message, "\n")
#    NULL
#  })
#
#  if (!is.null(res)) {
#    HDI_neg_1k[[sample_name]] <- res
#    cat("  ✅ Success (", ncol(res), " pixels, ", nrow(res), " peaks)\n", sep = "")
#  }
#}

neg_rows <- DESI_meta_data_joined %>% filter(DESI_mode == "negative")
cat("Found", nrow(neg_rows), "negative mode files to load\n")
IMZML_neg_1k <- list()

for (i in seq_len(nrow(neg_rows))) {
  sample_name <- neg_rows$sample[i]
  f <- neg_rows$top_1k[i]

  if (is.na(f) || f == "") {
    cat("→ [", i, "/", nrow(neg_rows), "] Skipping ", sample_name, " (no file found)\n", sep = "")
    next
  }

  cat("→ [", i, "/", nrow(neg_rows), "] Loading ", sample_name, " ...\n", sep = "")
  cat("   Path: ", file.path(data_dir, f), "\n", sep = "")

  res <- tryCatch({
    getImzMLPeaks(file.path(data_dir, f), bin_size = 0.01,n_peaks = 1000)
  }, error = function(e) {
    cat("  ⚠️ Error: ", e$message, "\n")
    NULL
  })

  if (!is.null(res)) {
    IMZML_neg_1k[[sample_name]] <- res
    cat("  ✅ Success (", ncol(res), " pixels, ", nrow(res), " peaks)\n", sep = "")
  }
}

mwass_Rbase = function(n, nx, m, mx){
  i = 1
  j = 1
  ret = 0
  #gamma = matrix(0, nrow = length(n), ncol = length(m))
  while (i <= length(n) && j <= length(m)){
    d = min(n[i], m[j])
    ret = ret + d * abs(nx[i] - mx[j])
    n[i] = n[i] - d
    m[j] = m[j] - d
    #gamma[i, j] = d
    if (n[i] == 0)
      i = i + 1
    else
      j = j + 1
  }
  #list(ret = ret, gamma = gamma)
  ret
}

library(Rcpp)
sourceCpp("/storage/homefs/ha25g949/pannet_metabolism/mwass.cpp")

# -----------------------------
# LOAD PACKAGES
# -----------------------------
library(bbknnR)
library(leidenbase)
library(igraph)

# -----------------------------
# INPUT DATA
# -----------------------------
# Assume:
#   data_list = list of feature matrices, one per image
#   img_names = vector of image names, same length as data_list
#   IMZML_neg_1k = list of SingleCellExperiment or Spatial objects with coordinates

# img_indices <- eval(parse(text = args[1]))

# to get out the normal tissue images from inside the R script
#img_names <- names(IMZML_neg_1k)[!grepl("(NN|nN)$", names(IMZML_neg_1k))]

# to get only the normal tissue images from inside the R script
#img_names <- names(IMZML_neg_1k)[grepl("(NN|nN)$", names(IMZML_neg_1k))]

#img_names = names(IMZML_neg_1k)[img_indices]
img_names = names(IMZML_neg_1k)

#data_list = lapply(img_names, function(nm) {
#  tmp = IMZML_neg_1k[[nm]]
#  assay(tmp, "TIC") = t(t(assay(tmp, "counts")) / colSums(assay(tmp, "counts")))
#  assay(tmp, "TIC")
#})

data_list <- lapply(img_names, function(nm) {
  tmp <- IMZML_neg_1k[[nm]]

  # Normalize counts
  assay(tmp, "TIC") <- t(t(assay(tmp, "counts")) / colSums(assay(tmp, "counts")))

  # Extract matrix
  mat <- assay(tmp, "TIC")

  # Assign m/z values as rownames (rounded to avoid floating mismatches)
  if (!is.null(rowData(tmp)$mz)) {
    rownames(mat) <- round(rowData(tmp)$mz, 4)
  } else {
    warning(paste("No m/z values found for", nm))
  }

  mat
})

# Remove any NULL entries (failed loads)
IMZML_neg_1k <- Filter(Negate(is.null), IMZML_neg_1k)

# Safety check
if (length(IMZML_neg_1k) == 0) {
  stop("No valid IMZML_neg_1k objects were loaded. Check your input files.")
}

# Sanity check 
data_list <- Filter(Negate(is.null), data_list)
common_rows <- Reduce(intersect, lapply(data_list, rownames))
data_list <- lapply(data_list, function(x) x[common_rows, , drop = FALSE])

# Combine all data into one big matrix
all_data <- do.call(cbind, data_list)

# Create batch labels (one per cell)
batch_vector <- rep(seq_along(data_list), times = sapply(data_list, ncol))

# bbknnR expects cells (spectra) as rows → transpose
embedding <- t(all_data)

# -----------------------------
# RUN BBKNN 
# -----------------------------

bbknn_result <- RunBBKNN(
  object = embedding,
  batch_list = batch_vector,
  neighbors_within_batch = 3, # reduce from 10 if there are small regions which connect to others but shouldnt
  method = "nndescent",
  metric = "euclidean"
)

# Extract the connectivity graph
graph <- bbknn_result$connectivities
if (is.null(graph)) {
  stop("No connectivities returned by RunBBKNN(). Check bbknnR version or returned object structure.")
}

# Ensure sparse matrix format
if (!inherits(graph, "dgCMatrix")) {
  graph <- as(graph, "dgCMatrix")
}

# -----------------------------
# CONVERT TO iGRAPH
# -----------------------------
# Leidenbase requires an igraph object (weighted, undirected)
g <- graph_from_adjacency_matrix(graph, mode = "undirected", weighted = TRUE, diag = FALSE)

# -----------------------------
# RUN LEIDEN CLUSTERING
# -----------------------------
leiden_res <- leiden_find_partition(
  g,
  partition_type = "RBConfigurationVertexPartition",
  resolution_parameter = 1.6 #0.8 gives us 4 clusters for the tumors, make finer clusters by increasing the resolution paramenter
)


leiden_clusters <- as.integer(leiden_res$membership)

# -----------------------------
# SPLIT CLUSTERS BY IMAGE
# -----------------------------
n_cells_img <- sapply(data_list, ncol)
clust_split <- split(leiden_clusters, rep(seq_along(n_cells_img), times = n_cells_img))

# -----------------------------
# CONSISTENT CLUSTER COLORS
# -----------------------------

# Get all unique clusters across all images
all_clusters <- sort(unique(leiden_clusters))

# Define a fixed global color palette
#cluster_colors <- setNames(
#  RColorBrewer::brewer.pal(min(length(all_clusters), 12), "Paired"),
#  all_clusters
#)

# for a palette with more than 12 colors
library(colorspace)

cluster_colors <- setNames(
  qualitative_hcl(length(all_clusters), palette = "Dynamic"),
  all_clusters
)

# -----------------------------
# PLOT RESULTS
# -----------------------------
name_output_dir <- args[2]

if (!dir.exists(name_output_dir)) dir.create(name_output_dir)

for (idx in seq_along(img_names)) {
  nm <- img_names[idx]
  tmp <- IMZML_neg_1k[[nm]]

  plot_data <- as.data.frame(colData(tmp))
  plot_data$clust <- as.character(clust_split[[idx]])

  p <- ggplot(plot_data, aes(x = y, y = x, fill = clust)) +
    geom_tile() +
    coord_fixed() +
    scale_y_reverse() +
    scale_x_reverse() +
    theme_classic() +
    labs(title = paste0(nm, " (Leiden)")) +
    scale_fill_manual(
      values = cluster_colors,
      na.value = "grey80",
      drop = FALSE
    )

  print(p)
  ggsave(
    filename = file.path(name_output_dir, paste0(nm, "_bbknn_leiden_31_plot.png")),
    plot = p, width = 6, height = 5, dpi = 300
  )
}

message("Leiden Plots saved")

# -----------------------------
# PCOA VISUALIZATION
# -----------------------------
library(ape)
library(ggplot2)

# Compute distance matrix
# use the BBKNN embedding 
dist_mat <- dist(embedding)

# Perform Principal Coordinate Analysis (PCoA)
pcoa_result <- ape::pcoa(dist_mat)

# Extract first two PCoA components
pcoa_df <- as.data.frame(pcoa_result$vectors[, 1:2])
colnames(pcoa_df) <- c("PCoA1", "PCoA2")

#cluster + image labels
pcoa_df$cluster <- as.factor(leiden_clusters)
pcoa_df$image <- rep(img_names, times = n_cells_img)

# Use the same consistent global cluster colors
my_colors <- cluster_colors

# -----------------------------
# PLOT PER-IMAGE PCOA
# -----------------------------
if (!dir.exists(name_output_dir)) dir.create(name_output_dir)
for (nm in img_names) {
  tmp_df <- subset(pcoa_df, image == nm)

  p <- ggplot(tmp_df, aes(x = PCoA1, y = PCoA2, color = cluster)) +
    geom_point(size = 2.5, alpha = 0.8) +
    scale_color_manual(values = my_colors, drop = FALSE) +
    labs(
      title = paste("PCoA -", nm),
      x = "PCoA1",
      y = "PCoA2"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "right"
    )

  print(p)
  ggsave(
    filename = file.path(name_output_dir, paste0(nm, "bbknn_leiden_31_pcoa_plot.png")),
    plot = p, width = 6, height = 5, dpi = 300
  )
}

message("PCoA plots saved")

# -----------------------------
# GLOBAL PCOA FOR ALL IMAGES
# -----------------------------
library(ape)
library(ggplot2)

# Compute distance matrix on all spectra
dist_mat <- dist(embedding)

# Perform PCoA
pcoa_result <- ape::pcoa(dist_mat)

# Extract first two PCoA components
pcoa_df <- as.data.frame(pcoa_result$vectors[, 1:2])
colnames(pcoa_df) <- c("PCoA1", "PCoA2")

# Add Leiden cluster assignments and image labels
pcoa_df$cluster <- as.factor(leiden_clusters)
pcoa_df$image <- rep(img_names, times = n_cells_img)

my_colors <- cluster_colors

# -----------------------------
# GLOBAL PCOA PLOT
# -----------------------------
if (!dir.exists(name_output_dir)) dir.create(name_output_dir)
p_global <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = cluster)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = my_colors, drop = FALSE) +
  labs(
    x = "PCoA1",
    y = "PCoA2",
    title = "Global PCoA of All Images (Leiden Clusters)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid = element_line(color = "grey90"),
    plot.title = element_text(hjust = 0.5)
  )

print(p_global)

# Save figure
ggsave(
  filename = file.path(name_output_dir, "_bbknn_leiden_31_Global_PCoA_All_Images.png"),
  plot = p_global, width = 7, height = 6, dpi = 300
)

message("Global PCoA plot saved")