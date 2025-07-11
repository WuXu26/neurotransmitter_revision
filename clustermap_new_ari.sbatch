#!/bin/bash
#SBATCH -p workq
#SBATCH -n 64
#SBATCH -t 05:00:00
#SBATCH -A loni_tsr_6
#SBATCH -N 1
#SBATCH -J clustermap
#SBATCH -o output_clustermap_new_ari.out
#SBATCH -e error_clustermap_new_ari.err

module load python/3.11.5-anaconda

python3 -m venv myenv
source myenv/bin/activate
pip3 install --upgrade pip
pip3 install -r /ddnB/work/wxx6941/TSR/code/code/neurotransmitter_revision/output_tsr/output_tsr.h5_res_csv/requirements.txt

python3 /ddnB/work/wxx6941/TSR/code/code/neurotransmitter_revision/output_tsr/output_tsr.h5_res_csv/clustermap_new_ari.py -p /ddnB/work/wxx6941/TSR/code/code/neurotransmitter_revision/output_tsr/output_tsr.h5_res_csv/generalised.csv
