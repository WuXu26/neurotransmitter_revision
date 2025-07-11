#!/bin/bash
#SBATCH -p workq
#SBATCH -n 48
#SBATCH -t 50:00:00
#SBATCH -A loni_tsr_6
#SBATCH -N 1
#SBATCH -J usr_similarity_atomic_distance_asp_glu
#SBATCH -o %j-output_usr_similarity_atomic_distance_asp_glu.out
#SBATCH -e %j-error_usr_similarity_atomic_distance_asp_glu.err

cd $SLURM_SUBMIT_DIR

module load python/3.7.6
python usr_similarity_atomic_distance_asp_glu.py





