#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --job-name="wd_hdi1000"
#SBATCH --mail-user=hector.arribasarias@unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --cpus-per-task=8
#SBATCH --mem=500G
#SBATCH --error=/storage/homefs/ha25g949/pannet_metabolism/scripts/log/%x_%j_bbknn.e
#SBATCH --output=/storage/homefs/ha25g949/pannet_metabolism/scripts/log/%x_%j_bbknn.o
#SBATCH --partition=epyc2

module load R/4.4.2-gfbf-2024a
module load Workspace_Home

echo "[$(date)] Starting job $SLURM_JOB_NAME ($SLURM_JOB_ID)"

# Define input/output directories
input_dir="/storage/research/igmp_grp_perren/raw_data_DESI_imaging/PanNET_Umara/lockmass_corrected"
output_base="/storage/homefs/ha25g949/pannet_metabolism/scripts/log/results"

# List all files
mapfile -t all_files < <(find "$input_dir" -maxdepth 1 -type f -printf "%f\n")

######################################################################
# 1. Run script WITHOUT normal samples (NN)
indices=()
for i in "${!all_files[@]}"; do
    file="${all_files[$i]}"
    if [[ ! "$file" =~ NN ]]; then
        indices+=($((i+1)))
    fi
done
r_vector="c($(IFS=,; echo "${indices[*]}"))"
Rscript /storage/homefs/ha25g949/pannet_metabolism/scripts/wd_bbknn_leiden_hdi_gen_attmpt2.R \
    "$r_vector" "$output_base/results_without_normals" "$input_dir"

######################################################################
# 2. Run script WITH ONLY normal samples (NN)
indices=()
for i in "${!all_files[@]}"; do
    file="${all_files[$i]}"
    if [[ "$file" =~ NN ]]; then
        indices+=($((i+1)))
    fi
done
r_vector="c($(IFS=,; echo "${indices[*]}"))"
Rscript /storage/homefs/ha25g949/pannet_metabolism/scripts/wd_bbknn_leiden_hdi_gen_attmpt2.R \
    "$r_vector" "$output_base/results_normals" "$input_dir"

######################################################################
# 3. Run script with ALL samples
indices=()
for i in "${!all_files[@]}"; do
    indices+=($((i+1)))
done
r_vector="c($(IFS=,; echo "${indices[*]}"))"
Rscript /storage/homefs/ha25g949/pannet_metabolism/scripts/wd_bbknn_leiden_hdi_gen_attmpt2.R \
    "$r_vector" "$output_base/results_all" "$input_dir"

