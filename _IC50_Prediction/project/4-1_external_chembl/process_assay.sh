#!/usr/bin/bash
#SBATCH -J ChEMBL_Assay
#SBATCH -c 1
#SBATCH -p ile
#SBATCH --mem=10GB
#SBATCH -o out/%j.out
#SBATCH -e out/e%j.out

source activate Process
python process_assay.py
