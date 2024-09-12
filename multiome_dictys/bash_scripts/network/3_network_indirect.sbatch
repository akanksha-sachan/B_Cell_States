#!/bin/bash
#SBATCH --job-name="postprocess"
#SBATCH --output="links_postprocess_S1.log"

#SBATCH -p RM-512
#SBATCH -N 1 
#SBATCH -t 3:00:00 # wall-time
#SBATCH --ntasks-per-node=128 # number of cores per node; memory is 2gb per node for rm-shared

# Set the working directory to where data and makefiles are located
WORK_DIR="/ocean/projects/cis240075p/asachan/datasets/B_Cell/multiome_1st_donor_UPMC_aggr/dictys_outs"
# Change directory to the working directory
cd $WORK_DIR || { echo "Error: Could not change directory to $WORK_DIR"; exit 1; }

################ RUN dictys in dynamic mode ################
module load anaconda3/2022.10
source activate dictys

# get indirect effects of grn (feedback loops)
dictys network indirect --nth 64 --fi_meanvar tmp_dynamic/Subset1/net_meanvar.tsv.gz tmp_dynamic/Subset1/net_weight.tsv.gz tmp_dynamic/Subset1/net_covfactor.tsv.gz tmp_dynamic/Subset1/net_iweight.tsv.gz
# normalize the indirect effects
dictys network normalize --nth 64 tmp_dynamic/Subset1/net_iweight.tsv.gz tmp_dynamic/Subset1/net_meanvar.tsv.gz tmp_dynamic/Subset1/net_covfactor.tsv.gz tmp_dynamic/Subset1/net_inweight.tsv.gz
