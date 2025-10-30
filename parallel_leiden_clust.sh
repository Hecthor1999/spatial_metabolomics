#!/bin/bash

#SBATCH --time=04:00:00
#SBATCH --mail-user=hector.arribasarias@unibe.ch
#SBATCH --cpus-per-task=8
#SBATCH --job-name=leiden_clustering
#SBATCH --partition=epyc2
#SBATCH --mail-type=end
#SBATCH --mem=8G
#SBATCH --output=/storage/homefs/ha25g949/pannet_metabolism/parallel/rms/total/negative/frequency_0.1/bg_rm/log/leiden_clustering_%j_%A.o
#SBATCH --error=/storage/homefs/ha25g949/pannet_metabolism/parallel/rms/total/negative/frequency_0.1/bg_rm/log/leiden_clustering_%j_%A.e
#SBATCH --array=1-16

module load Workspace_Home
module load Anaconda3/2024.02-1
eval "$(conda shell.bash hook)"

conda activate json_qpath_env

mkdir -p "/storage/homefs/ha25g949/pannet_metabolism/parallel/rms/total/negative/frequency_0.1/bg_rm/leiden_clustering/peak_pick_ref"
cd "/storage/homefs/ha25g949/pannet_metabolism/parallel/rms/total/negative/frequency_0.1/bg_rm/peak_pick_ref"

# define the directory containing input files
INPUT_DIR="/storage/homefs/ha25g949/pannet_metabolism/parallel/rms/total/negative/frequency_0.1/bg_rm/peak_pick_ref"

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
python /storage/homefs/ha25g949/pannet_metabolism/scripts/parallelized/spatial_metabolomics/parallelized_leiden_clust.py \
"/storage/homefs/ha25g949/pannet_metabolism/parallel/rms/total/negative/frequency_0.1/bg_rm/peak_pick_ref/$INPUT_FILE" \
"/storage/homefs/ha25g949/pannet_metabolism/parallel/rms/total/negative/frequency_0.1/bg_rm/leiden_clustering/peak_pick_ref" \
"$INPUT_FILE"


#mkdir -p "/storage/homefs/ha25g949/pannet_metabolism/parallel/rms/total/positive/frequency_0.1/leiden_clustering"
#cd "/storage/homefs/ha25g949/pannet_metabolism/parallel/rms/total/positive/frequency_0.1/peak_pick_ref"
#
# # define the directory containing input files
#INPUT_DIR="/storage/homefs/ha25g949/pannet_metabolism/parallel/rms/total/positive/frequency_0.1/peak_pick_ref"
# 
# # Create a file list dynamically
#FILE_LIST=($(ls $INPUT_DIR/*.imzML | awk -F'/' '{print $NF}'))
#TOTAL_FILES=${#FILE_LIST[@]}
#
## Check if SLURM_ARRAY_TASK_ID is within range
#if [ "$SLURM_ARRAY_TASK_ID" -gt "$TOTAL_FILES" ]; then
#    echo "Job index $SLURM_ARRAY_TASK_ID exceeds total file count ($TOTAL_FILES). Exiting."
#    exit 1
#fi
#
## Select the file corresponding to this job array index
#INPUT_FILE=${FILE_LIST[$SLURM_ARRAY_TASK_ID-1]}
#
## Run R script with the selected file
#python /storage/homefs/ha25g949/pannet_metabolism/scripts/parallelized/spatial_metabolomics/parallelized_leiden_clust.py \
#"/storage/homefs/ha25g949/pannet_metabolism/parallel/rms/total/positive/frequency_0.1/peak_pick_ref/$INPUT_FILE" \
#"/storage/homefs/ha25g949/pannet_metabolism/parallel/rms/total/positive/frequency_0.1/leiden_clustering" \
#"$INPUT_FILE"

