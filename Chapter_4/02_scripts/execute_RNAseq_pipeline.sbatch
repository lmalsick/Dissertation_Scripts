#!/usr/bin/env bash

#SBATCH --job-name=RNAseq_pipeline 
#SBATCH --nodes=1                          
#SBATCH --ntasks=20                       
#SBATCH --time=48:00:00                    
#SBATCH --partition=amilan                    
#SBATCH --qos=long                        
#SBATCH --mail-type=END                     
#SBATCH --mail-user=######@colostate.edu            # add your e-mail here
#SBATCH --output=log_RNAseq_pipe_%j.txt     
export TMPDIR=$SLURM_SCRATCH
export TMP=$TMPDIR

######### Instructions ###########

# Modify your SLURM entries above to fit your choices
# Below, modify the SECOND argument to point to YOUR metadata.file
# Below, you don't need to change $SLURM_NTASKS. It will automatically populate whatever you put in --ntasks=# above.
# Execute this script using $ sbatch execute_RNAseq_pipeline.sbatch



## Execute the RNA-seq_pipeline to run the pipeline
bash RNAseq_analyzer_250408.sh ../01_input/metadata.txt  $SLURM_NTASKS
