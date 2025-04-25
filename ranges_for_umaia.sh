#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=5G
#SBATCH --output=/storage/homefs/ha25g949/pannet_metabolism/scripts/parallelized/log/%x_%A_%a.o
#SBATCH --error=/storage/homefs/ha25g949/pannet_metabolism/scripts/parallelized/log/%x_%A_%a.e
#SBATCH --job-name=create_ranges
#SBATCH --mail-type=end
#SBATCH --mail-user=hector.arribasarias@unibe.ch
#SBATCH --partition=epyc2 
#SBATCH --array=1-30 # placeholder, the total number of jobs is updated automatically


module load R/4.4.2-gfbf-2024a
module load Workspace_Home

# for negative

#mkdir -p /storage/homefs/ha25g949/pannet_metabolism/parallel/total/negative/ranges/ranges_peakpick_reference
#
#DATA_PATH=/storage/homefs/ha25g949/pannet_metabolism/parallel/total/negative/peak_pick_ref
#OUTPUT_PATH=/storage/homefs/ha25g949/pannet_metabolism/parallel/total/negative/ranges/ranges_peakpick_reference
#
#FILE_LIST=($(ls $DATA_PATH/*.imzML | awk -F'/' '{print $NF}'))
#TOTAL_FILES=${#FILE_LIST[@]}
#
## Check if SLURM_ARRAY_TASK_ID is within range
#if [ "$SLURM_ARRAY_TASK_ID" -gt "$TOTAL_FILES" ]; then
#    echo "Job index $SLURM_ARRAY_TASK_ID exceeds total file count ($TOTAL_FILES). Exiting."
#    exit 1
#fi
#
#INPUT_FILE=${FILE_LIST[$SLURM_ARRAY_TASK_ID-1]}
#
#Rscript create_ranges_for_umaia.R $DATA_PATH $OUTPUT_PATH $INPUT_FILE
#
## for positive
#
#mkdir -p /storage/homefs/ha25g949/pannet_metabolism/parallel/total/positive/ranges/ranges_peakpick_reference
#
#DATA_PATH=/storage/homefs/ha25g949/pannet_metabolism/parallel/total/positive/peak_pick_ref
#OUTPUT_PATH=/storage/homefs/ha25g949/pannet_metabolism/parallel/total/positive/ranges/ranges_peakpick_reference
#
#FILE_LIST=($(ls $DATA_PATH/*imzML | awk -F '/' '{print $NF}'))
#TOTAL_FILES=${#FILE_LIST[@]}
#
## Check if SLURM_ARRAY_TASK_ID is within range
#if [ "$SLURM_ARRAY_TASK_ID" -gt "$TOTAL_FILES" ]; then
#    echo "Job index $SLURM_ARRAY_TASK_ID exceeds total file count ($TOTAL_FILES). Exiting."
#    exit 1
#fi
#
#INPUT_FILE=${FILE_LIST[$SLURM_ARRAY_TASK_ID-1]}
#
#Rscript create_ranges_for_umaia.R $DATA_PATH $OUTPUT_PATH $INPUT_FILE

# for negative

#mkdir -p /storage/homefs/ha25g949/pannet_metabolism/parallel/total/negative/ranges/ranges_peakpick_reference
#
#DATA_PATH=/storage/homefs/ha25g949/pannet_metabolism/parallel/total/negative/peak_pick_ref
#OUTPUT_PATH=/storage/homefs/ha25g949/pannet_metabolism/parallel/total/negative/ranges/ranges_peakpick_reference
#
#FILE_LIST=($(ls $DATA_PATH/*.imzML | awk -F'/' '{print $NF}'))
#TOTAL_FILES=${#FILE_LIST[@]}
#
## Check if SLURM_ARRAY_TASK_ID is within range
#if [ "$SLURM_ARRAY_TASK_ID" -gt "$TOTAL_FILES" ]; then
#    echo "Job index $SLURM_ARRAY_TASK_ID exceeds total file count ($TOTAL_FILES). Exiting."
#    exit 1
#fi
#
#INPUT_FILE=${FILE_LIST[$SLURM_ARRAY_TASK_ID-1]}
#
#Rscript /storage/homefs/ha25g949/pannet_metabolism/scripts/parallelized/create_ranges_for_umaia.R $DATA_PATH $OUTPUT_PATH $INPUT_FILE
#
## for positive
#
#mkdir -p /storage/homefs/ha25g949/pannet_metabolism/parallel/total/positive/ranges/ranges_peakpick_reference
#
#DATA_PATH=/storage/homefs/ha25g949/pannet_metabolism/parallel/total/positive/peak_pick_ref
#OUTPUT_PATH=/storage/homefs/ha25g949/pannet_metabolism/parallel/total/positive/ranges/ranges_peakpick_reference
#
#FILE_LIST=($(ls $DATA_PATH/*imzML | awk -F '/' '{print $NF}'))
#TOTAL_FILES=${#FILE_LIST[@]}
#
## Check if SLURM_ARRAY_TASK_ID is within range
#if [ "$SLURM_ARRAY_TASK_ID" -gt "$TOTAL_FILES" ]; then
#    echo "Job index $SLURM_ARRAY_TASK_ID exceeds total file count ($TOTAL_FILES). Exiting."
#    exit 1
#fi
#
#INPUT_FILE=${FILE_LIST[$SLURM_ARRAY_TASK_ID-1]}
#
#Rscript /storage/homefs/ha25g949/pannet_metabolism/scripts/parallelized/create_ranges_for_umaia.R $DATA_PATH $OUTPUT_PATH $INPUT_FILE


# for negative

mkdir -p /storage/homefs/ha25g949/pannet_metabolism/parallel/total/negative/frequency_0.01/ranges

DATA_PATH=/storage/homefs/ha25g949/pannet_metabolism/parallel/total/negative/frequency_0.01/peak_aligned
OUTPUT_PATH=/storage/homefs/ha25g949/pannet_metabolism/parallel/total/negative/frequency_0.01/ranges

FILE_LIST=($(ls $DATA_PATH/*.imzML | awk -F'/' '{print $NF}'))
TOTAL_FILES=${#FILE_LIST[@]}

# Check if SLURM_ARRAY_TASK_ID is within range
if [ "$SLURM_ARRAY_TASK_ID" -gt "$TOTAL_FILES" ]; then
    echo "Job index $SLURM_ARRAY_TASK_ID exceeds total file count ($TOTAL_FILES). Exiting."
    exit 1
fi

INPUT_FILE=${FILE_LIST[$SLURM_ARRAY_TASK_ID-1]}

Rscript /storage/homefs/ha25g949/pannet_metabolism/scripts/parallelized/create_ranges_for_umaia.R $DATA_PATH $OUTPUT_PATH $INPUT_FILE

# for positive

mkdir -p /storage/homefs/ha25g949/pannet_metabolism/parallel/total/positive/frequency_0.01/ranges

DATA_PATH=/storage/homefs/ha25g949/pannet_metabolism/parallel/total/positive/frequency_0.01/peak_aligned
OUTPUT_PATH=/storage/homefs/ha25g949/pannet_metabolism/parallel/total/positive/frequency_0.01/ranges

FILE_LIST=($(ls $DATA_PATH/*imzML | awk -F '/' '{print $NF}'))
TOTAL_FILES=${#FILE_LIST[@]}

# Check if SLURM_ARRAY_TASK_ID is within range
if [ "$SLURM_ARRAY_TASK_ID" -gt "$TOTAL_FILES" ]; then
    echo "Job index $SLURM_ARRAY_TASK_ID exceeds total file count ($TOTAL_FILES). Exiting."
    exit 1
fi

INPUT_FILE=${FILE_LIST[$SLURM_ARRAY_TASK_ID-1]}

Rscript /storage/homefs/ha25g949/pannet_metabolism/scripts/parallelized/create_ranges_for_umaia.R $DATA_PATH $OUTPUT_PATH $INPUT_FILE