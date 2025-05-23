#!/usr/bin/bash
#SBATCH -J IC50_ChEMBL_7
#SBATCH -c 1
#SBATCH -p ile
#SBATCH --mem=1GB
#SBATCH -o out/%j.out
#SBATCH -e out/e%j.out
source activate Process
python process_ic50_batch.py -nth 7 -total 1
