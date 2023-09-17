#!/bin/bash
#SBATCH --job-name=MSR_project
#SBATCH --output=%x_%j.out
#SBATCH --mail-user=h.van.der.pol@umail.leidenuniv.nl
#SBATCH --mail-type=ALL
#SBATCH --partition=cpu-medium,cpu-long,mem
#SBATCH --time=23:00:00
#SBATCH --ntasks=1
#SBATCH --mem=10G
#SBATCH --nodes=1

module load R/4.0.5-foss-2020b
Rscript MSR.R