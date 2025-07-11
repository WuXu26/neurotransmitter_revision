#!/bin/bash
#SBATCH -p workq
#SBATCH -n 64
#SBATCH -t 00:10:00
#SBATCH -A loni_tsr_7
#SBATCH -N 1
#SBATCH -J input_converter
#SBATCH -o output_input_converter.out
#SBATCH -e error_input_converter.err

module purge
module load python/3.11.5-anaconda

python3 input_converter.py
