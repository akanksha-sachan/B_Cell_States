#!/bin/bash
#SBATCH --job-name="TF_TSS"
#SBATCH --output="binlinking_S1.log"

#SBATCH -p RM
#SBATCH -N 1 # by default for gpu-shared you can only request 1 node
#SBATCH -t 6:00:00 # wall-time
#SBATCH --ntasks-per-node=128 # number of cores per node; memory is 2gb per node for rm-shared

# Set the working directory to where data and makefiles are located
WORK_DIR="/ocean/projects/cis240075p/asachan/datasets/B_Cell/multiome_1st_donor_UPMC_aggr/dictys_outs"

# Change directory to the working directory
cd $WORK_DIR || { echo "Error: Could not change directory to $WORK_DIR"; exit 1; }

################ RUN dictys in dynamic mode ################
# bash runs all commands sequentially
module load anaconda3/2022.10
source activate dictys
dictys chromatin binding tmp_dynamic/Subset1/wellington.tsv.gz tmp_dynamic/Subset1/homer.tsv.gz tmp_dynamic/Subset1/binding.tsv.gz
dictys chromatin tssdist tmp_dynamic/Subset1/expression.tsv.gz tmp_dynamic/Subset1/wellington.tsv.gz data/gene.bed tmp_dynamic/Subset1/tssdist.tsv.gz
dictys chromatin linking tmp_dynamic/Subset1/binding.tsv.gz tmp_dynamic/Subset1/tssdist.tsv.gz tmp_dynamic/Subset1/linking.tsv.gz
# top 20 TF-Target links are selected to create adjacency matrix of TF-Target network
dictys chromatin binlinking tmp_dynamic/Subset1/linking.tsv.gz tmp_dynamic/Subset1/binlinking.tsv.gz 20
