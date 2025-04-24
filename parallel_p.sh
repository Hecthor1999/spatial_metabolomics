#!/bin/bash
#SBATCH --time=05:00:00
#SBATCH --ntasks=1
#SBATCH --job-name="cardinal_positive"
#SBATCH --mail-user=hector.arribasarias@unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=32G
#SBATCH --error=/storage/homefs/ha25g949/pannet_metabolism/scripts/log/%x_%A_%a.e
#SBATCH --output=/storage/homefs/ha25g949/pannet_metabolism/scripts/log/%x_%A_%a.o
#SBATCH --partition=epyc2 
#SBATCH --array=1-30 # placeholder, the total number of jobs is updated automatically

module load R/4.4.2-gfbf-2024a
module load Workspace_Home

export R_LIBS=$HOME/R/library
mkdir -p $R_LIBS

# positive 
mkdir -p /storage/homefs/ha25g949/pannet_metabolism/parallel
mkdir -p /storage/homefs/ha25g949/pannet_metabolism/parallel/unnormalized
mkdir -p /storage/homefs/ha25g949/pannet_metabolism/parallel/unnormalized/total
mkdir -p /storage/homefs/ha25g949/pannet_metabolism/parallel/unnormalized/total/positive
mkdir -p /storage/homefs/ha25g949/pannet_metabolism/parallel/unnormalized/total/positive/frequency_0.1
mkdir -p /storage/homefs/ha25g949/pannet_metabolism/parallel/unnormalized/total/positive/frequency_0.1/preprocessed
mkdir -p /storage/homefs/ha25g949/pannet_metabolism/parallel/unnormalized/total/positive/frequency_0.1/peak_picked
mkdir -p /storage/homefs/ha25g949/pannet_metabolism/parallel/unnormalized/total/positive/frequency_0.1/peak_aligned
mkdir -p /storage/homefs/ha25g949/pannet_metabolism/parallel/unnormalized/total/positive/frequency_0.1/peak_pick_ref
mkdir -p /storage/homefs/ha25g949/pannet_metabolism/parallel/unnormalized/total/positive/frequency_0.1/peak_aligned_non_filtered

# define the directory containing input files
INPUT_DIR="/storage/research/igmp_grp_perren/raw_data_DESI_imaging/PanNET_Umara/imzML/positive_mode"

# Create a file list dynamically
FILE_LIST=($(ls $INPUT_DIR/*.imzML | awk -F'/' '{print $NF}'))
TOTAL_FILES=${#FILE_LIST[@]}

# Check if SLURM_ARRAY_TASK_ID is within range
if [ "$SLURM_ARRAY_TASK_ID" -gt "$TOTAL_FILES" ]; then
    echo "Job index $SLURM_ARRAY_TASK_ID exceeds total file count ($TOTAL_FILES). Exiting."
    exit 1
fi

# Select the file corresponding to this job array index
INPUT_FILE=${FILE_LIST[$SLURM_ARRAY_TASK_ID-1]}

# Run R script with the selected file
Rscript /storage/homefs/ha25g949/pannet_metabolism/scripts/parallelized/cardinal_parallel.R \
 "$INPUT_DIR" \
 /storage/homefs/ha25g949/pannet_metabolism/parallel/unnormalized/total/positive/frequency_0.1/preprocessed \
 /storage/homefs/ha25g949/pannet_metabolism/parallel/unnormalized/total/positive/frequency_0.1/peak_picked \
 /storage/homefs/ha25g949/pannet_metabolism/parallel/unnormalized/total/positive/frequency_0.1/peak_aligned \
 1 \
 10 \
 0.1 \
 /storage/homefs/ha25g949/pannet_metabolism/parallel/unnormalized/total/positive/frequency_0.1/peak_aligned_non_filtered \
 /storage/homefs/ha25g949/pannet_metabolism/parallel/unnormalized/total/positive/frequency_0.1/peak_pick_ref \
 "$INPUT_FILE"