#!/usr/bin/bash
#SBATCH -J GSVA_PPI
#SBATCH -c 36
#SBATCH -p met
#SBATCH --mem=240GB
#SBATCH -o out/%j.out
#SBATCH -e out/e%j.out

source activate Process
Rscript gsva_ppa.R 5
