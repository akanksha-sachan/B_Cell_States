#!/bin/bash
#SBATCH --job-name="reconstruct"
#SBATCH --output="reconstruct_S1.log"

#SBATCH -p GPU-shared
#SBATCH -N 1 # by default for gpu-shared you can only request 1 node
#SBATCH -t 12:00:00 # wall-time
#SBATCH --gres=gpu:v100-16:1 # how many cards do you want?

# Set the working directory to where data and makefiles are located
WORK_DIR="/ocean/projects/cis240075p/asachan/datasets/B_Cell/multiome_1st_donor_UPMC_aggr/dictys_outs"
# Change directory to the working directory
cd $WORK_DIR || { echo "Error: Could not change directory to $WORK_DIR"; exit 1; }

########### run dictys network reconstruct ###########
module load anaconda3/2022.10
source activate dictys
# 40 cores on 1 node, try to give 20 next time for faster inference
dictys network reconstruct --device cuda:0 --nth 6 tmp_dynamic/Subset1/expression.tsv.gz tmp_dynamic/Subset1/binlinking.tsv.gz tmp_dynamic/Subset1/net_weight.tsv.gz tmp_dynamic/Subset1/net_meanvar.tsv.gz tmp_dynamic/Subset1/net_covfactor.tsv.gz tmp_dynamic/Subset1/net_loss.tsv.gz tmp_dynamic/Subset1/net_stats.tsv.gz