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

Rscript /storage/homefs/ha25g949/pannet_metabolism/scripts/wd_bbknn_leiden_hdi_gen_attmpt.R
