# Load packages

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

if (!require("dplyr", quietly = TRUE))
    install.packages("dplyr", lib="~/R/library")

if (!require("Cardinal", quietly = TRUE))
BiocManager::install("Cardinal", lib="~/R/library")

if (!require("BiocParallel", quietly = TRUE))
BiocManager::install("BiocParallel", lib="~/R/library", force = TRUE)

# Load packages
library(Cardinal)
library(dplyr)


args <- commandArgs(trailingOnly=TRUE)

# load data
path_aligned_data <- args[1]
output_path <- args[2]
file_name <- args[3]
complete_path <- paste0(path_aligned_data,"/", file_name)
data <- readMSIData(complete_path)

# Calculate min and max for peak bins
mz_values <- mz(data)
mz_min <- mz_values - 0.5
mz_max <- mz_values + 0.5
pixel_max_hit <- rep(1, length(mz(data))) # dummy var
num_pixels <- nrow(data)*featureData(data)$freq
percent_1_hit <- rep(100.0, length(mz(data))) # dummy var
concentration <- rep(1,length(mz(data))) # dummy var
median_intensity <- apply(as.matrix(intensity(data)),MARGIN = 1,median)
mz_estimated <-  mz(data)
difference <- mz_max - mz_min

# Combine into a data frame
ranges_df <- data.frame(min=mz_min, max=mz_max, pixel_max_hits=pixel_max_hit,num_pixels=num_pixels, percent_1_hit=percent_1_hit, concentration = concentration,media_intensity = median_intensity, mz_estimated= mz_estimated, difference=difference)

# save
write.csv(ranges_df, file=file.path(output_path, paste0(file_name, ".csv")), row.names=TRUE)
