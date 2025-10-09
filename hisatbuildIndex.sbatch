#!/usr/bin/env bash

#SBATCH --job-name=execute_hisat2-build #for use on Alpine supercomputer as of 20251009
#SBATCH --nodes=1
#SBATCH --ntasks=20 
#SBATCH --partition=amilan
#SBATCH --qos=normal    
#SBATCH --time=0:20:00   
#SBATCH --mail-type=END                     
#SBATCH --mail-user=lmalsick@colostate.edu 
#SBATCH --output=log_hisat2-build_%J.txt

# Build hisat2 indexes for a. jamacensis --make sure file is in your directory! use .fna

hisat2-build -p ${SLURM_NTASKS} GCF_021234435.1_CSHL_Jam_final_genomic.fna ajgene

# Check the build

echo -e "\n\nINDEX-BUILD: inspecting indexes:"

hisat2-inspect -s ajgene

# Capture version number

echo -e "\n\nINDEX-BUILD: version:"

hisat2-build --version
