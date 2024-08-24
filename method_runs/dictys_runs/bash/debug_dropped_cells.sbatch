#!/bin/bash
#SBATCH --job-name="debug" 
#SBATCH --output="debug.log" 
#SBATCH -p RM-shared
#SBATCH -N 1 # number of nodes
#SBATCH --ntasks-per-node=32 # number of cores per node; memory is 2gb per node for rm-shared
#SBATCH -t 1:00:00

module load anaconda3/2022.10
source activate dictys
# set path to store all output files
exp_mtx_path="/ocean/projects/cis240075p/asachan/datasets/B_Cell/multiome_1st_donor_UPMC_aggr/dictys_outs/data/expression.tsv.gz"



