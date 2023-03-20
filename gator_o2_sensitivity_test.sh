#!/bin/bash
#SBATCH -p short
#SBATCH -J nextflow_O2
#SBATCH -o run.o
#SBATCH -e run.e
#SBATCH -t 0-12:00
#SBATCH --mem=8G
#SBATCH --mail-type=END
#SBATCH --mail-user=ajitj_nirmal@dfci.harvard.edu


conda activate scimap