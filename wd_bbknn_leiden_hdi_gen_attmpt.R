# This R script calculates the wassersteindistance by using the mwass.cpp file in /storage/homefs/ha25g949/pannet_metabolism/mwass.cpp
# then it calcualtes the batch balanced k nearest neighbour graph using bbknnR, runs Leiden clustering on the graph and then visualizes both the full image and the PCOAs for all images in either negative or positive modes

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

data_dir = "/storage/research/igmp_grp_perren/raw_data_DESI_imaging/PanNET_Umara/lockmass_corrected"

# the meta data is required to match the image files from the different modes to the correct sample
# This list is mainly hand curated 
DESI_meta_data = read.xlsx("/storage/homefs/ha25g949/pannet_metabolism/240208_sample_list_DESI_RNA_coregistration.xlsx") %>% 
  mutate(patient_uid = paste0(aPnumber, "_", B_extended))

# There is one sample that with a slightly different result that needs to be checked
HDI_1k_files = data.frame(top_1k = list.files(data_dir, pattern = "CLMC.txt", recursive = T)) %>% 
  filter(!grepl("ambiguous", top_1k)) %>% 
  mutate(sample = gsub(" Analyte.*$", "", basename(top_1k)))


DESI_meta_data = left_join(DESI_meta_data,
                           HDI_1k_files,
                           by = "sample")

rm(HDI_1k_files)

DESI_meta_data = DESI_meta_data %>% 
  filter(!is.na(top_1k))

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

# --- SAFE VERSION OF LOADER ---

neg_rows <- DESI_meta_data %>% filter(DESI_mode == "negative")
cat("Found", nrow(neg_rows), "negative mode files to load\n")
HDI_neg_1k <- list()

for (i in seq_len(nrow(neg_rows))) {
  sample_name <- neg_rows$sample[i]
  f <- neg_rows$top_1k[i]

  cat("→ [", i, "/", nrow(neg_rows), "] Loading ", sample_name, " ...\n", sep = "")

  res <- tryCatch({
    getHDIPeaks(f)
  }, error = function(e) {
    cat("  ⚠️ Error: ", e$message, "\n")
    NULL
  })

  if (!is.null(res)) {
    HDI_neg_1k[[sample_name]] <- res
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
library(SingleCellExperiment)
library(ggplot2)
library(bbknnR)
library(leidenbase)
library(Matrix)
library(igraph)

# -----------------------------
# INPUT DATA
# -----------------------------
# Assume:
#   data_list = list of feature matrices, one per image
#   img_names = vector of image names, same length as data_list
#   HDI_neg_1k = list of SingleCellExperiment or Spatial objects with coordinates

img_names = names(HDI_neg_1k)[c(1:31)]

data_list = lapply(img_names, function(nm) {
  tmp = HDI_neg_1k[[nm]]
  assay(tmp, "TIC") = t(t(assay(tmp, "counts")) / colSums(assay(tmp, "counts")))
  assay(tmp, "TIC")
})

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
  neighbors_within_batch = 10,
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
  resolution_parameter = 0.8
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
cluster_colors <- setNames(
  RColorBrewer::brewer.pal(min(length(all_clusters), 12), "Paired"),
  all_clusters
)

#cluster_colors <- setNames(
#  scales::hue_pal()(length(all_clusters)),
#  all_clusters
#)

# -----------------------------
# PLOT RESULTS
# -----------------------------
if (!dir.exists("results")) dir.create("results")

for (idx in seq_along(img_names)) {
  nm <- img_names[idx]
  tmp <- HDI_neg_1k[[nm]]

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
    filename = file.path("results", paste0(nm, "_bbknn_leiden_31_plot.png")),
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
    filename = file.path("results", paste0(nm, "bbknn_leiden_31_pcoa_plot.png")),
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
  filename = file.path("results", "_bbknn_leiden_31_Global_PCoA_All_Images.png"),
  plot = p_global, width = 7, height = 6, dpi = 300
)

message("Global PCoA plot saved")
