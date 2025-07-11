#!/bin/bash
#SBATCH -p workq
#SBATCH -n 48
#SBATCH -t 10:00:00
#SBATCH -A loni_tsr_6
#SBATCH -N 1
#SBATCH -J ari_generalised_asp_glu_usr_2g
#SBATCH -o %j-output_ari_generalised_asp_glu_usr_2g.out
#SBATCH -e %j-error_ari_generalised_asp_glu_usr_2g.err

cd $SLURM_SUBMIT_DIR

module load python/3.7.6
python ari_generalised_asp_glu_usr_2g.py





