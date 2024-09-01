#!/usr/bin/bash
#SBATCH -J Assay_7
#SBATCH -c 1
#SBATCH -p cath
#SBATCH -x lysine
#SBATCH --mem=1GB
#SBATCH -o out/%j.out
#SBATCH -e out/e%j.out
source activate Process
python process_assay.py -nth 7
