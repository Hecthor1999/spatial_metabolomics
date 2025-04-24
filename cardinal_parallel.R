# This Rscript aligns peaks using the cardinal package and it is called by the cardinal_align.sh bash script

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

args <- commandArgs(trailingOnly = TRUE)
# path to the directory where files are and saving directory
data_path <- args[1]
output_preprocessed_path <- args[2]
output_picked_path <- args[3]
output_aligned_path <- args[4]
# overwrite <- args[10]

# list with the names of the files we want to use
# get all imzml files if no specific ones specified
files <- if (length(args) > 1 && !is.na(args[10])){
    unlist(strsplit(args[10], ","))
} else {
    list.files(path=data_path, pattern="\\.imzML$", all.files=FALSE, full.names=FALSE)
}

#read data
mse <- list()
for (file in files) {
    file_path <- file.path(data_path, file)
    cat("Processing file:", file_path, "\n")
    
    if (!file.exists(file_path)) {
        cat("File not found:", file_path, "\n")
    } else {
        mse[[file]] <- readMSIData(file_path)  
    }
}

# preprocessing without recalibration
baseline_red_meth <- "snip"
extra_args <- "10"
# use as arguments "width" when performing snip and "smooth" when performing locmin.

mse <-  lapply(mse, function(data) {
    data %>%
    normalize(method = "rms") %>%
    smooth(method = "sgolay") %>%
    reduceBaseline(method = baseline_red_meth, width=extra_args)
  }) %>%
  lapply(process)


lapply(seq_along(mse), function(i) {
  file_name <- names(mse)[i]  
  output_file <- file.path(output_preprocessed_path, paste0(file_name, ".imzML"))  
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)  
  writeMSIData(mse[[i]], file = output_file, bundle = FALSE)
})

#thresholding
threshold <- args[5] # 1 

mse <- lapply(mse, function(data) {
    # modify spectra(data) directly to maintain structure
    spectra(data) <- lapply(spectra(data), function(spectrum) {
    spectrum[spectrum < threshold] <- 0
    return(spectrum)  
    })        
    return(data)  
})

# pick peaks
signal_noise_ratio <- args[6] # 10

picked_peaks <- mse %>%
    lapply(function(data) {
    data %>%
    peakPick(method="diff",SNR=signal_noise_ratio) 
  }) %>% lapply(process)

# do filtering method after first peak picking
picked_peaks <- picked_peaks %>%
  lapply(function(data) {
  data %>%
  peakPick(method="filter",SNR=signal_noise_ratio) #isPeaks = TRUE #try 
}) %>% lapply(process)

# write picked data
lapply(seq_along(picked_peaks), function(i) {
  file_name <- runNames(picked_peaks[[i]])  
  output_file <- file.path(output_picked_path, paste0(file_name, ".imzML"))  
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)  
  if (file.exists(output_file)) {
    file.remove(output_file)
  }
  writeMSIData(picked_peaks[[i]], file = output_file, bundle = FALSE)
})

align_freq_thres <- args[7] # 0.1

# align peaks
aligned <- lapply(picked_peaks, function(data) {
    aligned <- peakAlign(data, tolerance=100, units="ppm")  
    #subsetFeatures(aligned, freq > align_freq_thres)
  }) %>% lapply(process)

# write aligned data
lapply(seq_along(aligned), function(i) {
  file_name <- runNames(aligned[[i]])
  output_file <- file.path(output_aligned_path, paste0(file_name, ".imzML"))  
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)  
  writeMSIData(aligned[[i]], file = output_file, bundle = FALSE)
})

# do peak alignment and freq filtering
aligned_mse <- lapply(picked_peaks, function(data) {
    aligned <- peakAlign(data, tolerance=100, units="ppm")  
    subsetFeatures(aligned, freq > align_freq_thres)
  }) %>% lapply(process)


output_freq_aligned_path <- args[8]
# write freq filtered aligned data
lapply(seq_along(aligned_mse), function(i) {
  file_name <- runNames(aligned_mse[[i]])
  output_file <- file.path(output_freq_aligned_path, paste0(file_name, ".imzML"))  
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)  
  writeMSIData(aligned_mse[[i]], file = output_file, bundle = FALSE)
})

# peak picking based on a reference
ref_mz <- unlist(lapply(aligned_mse, mz))
reference_based_peak_pick <- lapply(mse, function(data) {
    peakPick(data, ref= ref_mz, type="area",
    tolerance=80, units="ppm")
  }) %>% lapply(process)

# write aligned data
output_reference_based__path <- args[9]

lapply(seq_along(reference_based_peak_pick), function(i) {
  file_name <- runNames(reference_based_peak_pick[[i]])
  output_file <- file.path(output_reference_based_peak_pick_path, paste0(file_name, ".imzML"))  
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)  
  writeMSIData(reference_based_peak_pick[[i]], file = output_file, bundle = FALSE)
})

# create ranges file for umaia molecular matching
