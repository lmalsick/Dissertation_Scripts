#!/usr/bin/env bash

#SBATCH --job-name=execute_hisat2-build
#SBATCH --nodes=1
#SBATCH --ntasks=20 
#SBATCH --partition=amilan
#SBATCH --qos=normal    
#SBATCH --time=0:30:00   
#SBATCH --mail-type=END                     
#SBATCH --mail-user=lmalsick@colostate.edu 
#SBATCH --output=log_hisat2-build_%J.txt

# Use this step if you want exons/splice sites for rMATs analysis
echo "Extracting splice sites and exons from GTF..."

hisat2_extract_splice_sites.py GCF_021234435.1_CSHL_Jam_final_genomic.gtf > ajgene_splice_sites.txt
hisat2_extract_exons.py GCF_021234435.1_CSHL_Jam_final_genomic.gtf > ajgene_exons.txt

echo "Extraction complete."

# Build hisat2 indexes for a. jamaicensis with splice site and exon information -- make sure .fna is in your directory!
echo "Building HISAT2 index with annotation..."

hisat2-build \
  -p ${SLURM_NTASKS} \
  --ss ajgene_splice_sites.txt \
  --exon ajgene_exons.txt \
  GCF_021234435.1_CSHL_Jam_final_genomic.fna \
  ajgene

# Check the build
echo -e "\n\nINDEX-BUILD: inspecting indexes:"

hisat2-inspect -s ajgene

# Capture version number
echo -e "\n\nINDEX-BUILD: version:"
hisat2-build --version
