#!/bin/bash
#SBATCH --job-name="homer_motifs" 
#SBATCH --output="homer_motifs.log" 
#SBATCH -p RM
#SBATCH -N 1 # number of nodes
#SBATCH --ntasks-per-node=128 # number of cores per node; memory is 2gb per node for rm-shared
#SBATCH -t 2:00:00


source activate dictys

dictys_data_path='/ocean/projects/cis240075p/asachan/datasets/B_Cell/multiome_1st_donor_UPMC_aggr/dictys_data'
# basedir is where homer is installed on PSC; -c = 2 for captialized gene names
dictys_helper motif_homer.sh -b basedir -o human -c 2 > $dictys_data_path/motifs.motif

