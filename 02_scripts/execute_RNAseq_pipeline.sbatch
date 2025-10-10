#!/usr/bin/env bash

#SBATCH --job-name=RNAseq_pipeline 
#SBATCH --nodes=1                           # this script is designed to run on one node
#SBATCH --ntasks=16                        
#SBATCH --time=04:00:00                    
#SBATCH --partition=amem                    
#SBATCH --qos=mem                        
#SBATCH --mail-type=END                    
#SBATCH --mail-user=lmalsick@colostate.edu            
#SBATCH --output=log_RNAseq_pipe_%j.txt     
export TMPDIR=$SLURM_SCRATCH
export TMP=$TMPDIR

######### Instructions ###########

# Modify your SLURM entries above to fit your choices
# Below, modify the SECOND argument to point to YOUR metadata.file
# Below, you don't need to change $SLURM_NTASKS. It will automatically populate whatever you put in --ntasks=# above.
# Execute this script using $ sbatch execute_RNAseq_pipeline.sbatch

## Execute the RNA-seq_pipeline to run the pipeline
bash RNAseq_analyzer_251009.sh ../01_input/metadata_infected.txt  $SLURM_NTASKS


## Execute the cleanup script to zip .fastq files and delete extra files
#bash RNAseq_cleanup_251009.sh ../01_input/metadata.txt
#this script crashes, not sure why
