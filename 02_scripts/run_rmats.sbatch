#!/usr/bin/env bash
#SBATCH --job-name=rMATS
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=2-00:00:00
#SBATCH --partition=amilan
#SBATCH --qos=long
#SBATCH --mail-type=END
#SBATCH --mail-user=lmalsick@colostate.edu
#SBATCH --output=log_rMATS_%j.txt

# Load modules if needed:
# module load hisat2 rMATS

# Execute:
# Usage: bash rMATS_run.sh <group1_bam_list.txt> <group2_bam_list.txt> <annotation.gtf> <output_dir> [threads] [read_length]
bash rMATS_run.sh wt_bamfiles.txt v816_bamfiles.txt /scratch/alpine/lmalsick@colostate.edu/Genewiz_sequencing/02_scripts/GCF_021234435.1_CSHL_Jam_final_genomic.gtf 20 150
