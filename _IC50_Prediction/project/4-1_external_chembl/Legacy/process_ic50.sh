#!/usr/bin/bash
#SBATCH -J ChEMBL_IC50
#SBATCH -c 1
#SBATCH -p ile
#SBATCH --mem=40GB
#SBATCH -o out/%j.out
#SBATCH -e out/e%j.out

source activate Process
python process_ic50.py -total 1
