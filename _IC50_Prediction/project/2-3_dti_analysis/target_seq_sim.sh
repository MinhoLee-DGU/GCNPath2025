#!/usr/bin/bash
#SBATCH -n 1
#SBATCH -c 48
#SBATCH -x lysine
#SBATCH -o %j.out
#SBATCH -e e%j.out
#SBATCH -J MoA_Seq
Rscript target_seq_sim.R
