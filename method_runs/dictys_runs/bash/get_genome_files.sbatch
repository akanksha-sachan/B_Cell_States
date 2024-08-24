#!/bin/bash
#SBATCH --job-name="genome_files" 
#SBATCH --output="genome_files.log" 
#SBATCH -p RM-shared
#SBATCH -N 1 # number of nodes
#SBATCH --ntasks-per-node=64 # number of cores per node; memory is 2gb per node for rm-shared
#SBATCH -t 3:00:00


source activate dictys
# set path to store all output files
dictys_data_path="/ocean/projects/cis240075p/asachan/datasets/B_Cell/multiome_1st_donor_UPMC_aggr/dictys_data"

####### get motifs pwm file (.motif format) #######
# homer 4.11.0 doesn't have the directory knownTFs on the psc installation (using hocomoco instead)
# wget -q -o /dev/null -O $dictys_data_path/motifs.motif 'https://hocomoco11.autosome.org/final_bundle/hocomoco11/full/HUMAN/mono/HOCOMOCOv11_full_HUMAN_mono_homer_format_0.0001.motif'

####### get reference genome #######
#dictys_helper genome_homer.sh -b /opt/packages/homer/4.11.0/ hg38 $dictys_data_path/genome
#BUG: the script trys to chmod perl scripts in the basedir of homer > not allowed for cluster homer.
#DEBUG: directly copy the ref from the path 
# cp -R  /opt/packages/homer/4.11.0/data/genomes/hg38 $dictys_data_path/genome

####### extract genes in bed format from gtf file #######
dictys_helper gene_gtf.sh $dictys_data_path/gene.gtf $dictys_data_path/gene.bed
rm $dictys_data_path/gene.gtf

