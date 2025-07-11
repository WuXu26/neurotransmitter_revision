#!/bin/bash
#SBATCH -p workq
#SBATCH -n 48
#SBATCH -t 10:00:00
#SBATCH -A loni_tsr_6
#SBATCH -N 1
#SBATCH -J Drug_Triplets_Keys_Cal
#SBATCH -o %j-output_Drug_Triplets_Keys_Cal.out
#SBATCH -e %j-error_Drug_Triplets_Keys_Cal.err

cd $SLURM_SUBMIT_DIR

module load python/3.7.6
python Drug_Triplets_Keys_Cal.py



