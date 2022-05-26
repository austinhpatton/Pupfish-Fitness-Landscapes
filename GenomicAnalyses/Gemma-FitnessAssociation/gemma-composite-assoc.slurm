#!/bin/bash
#SBATCH --output=GEMMA_%J.out
#SBATCH --mail-type=ALL
#SBATCH --time=72:00:00
#SBATCH -p savio2
#SBATCH --account=ac_fishes

# Run GEMMA analysis on survival
cd /global/scratch/austinhpatton/fitnessGenomics/results/GEMMA

PREFIX='/global/scratch/austinhpatton/fitnessGenomics/results/GEMMA'

# Conduct Gemma LMM to test for associations between genotypes and composite fitness. 
/global/home/users/austinhpatton/software/gemma-0.98.1 -bfile ${PREFIX}/Hybs-LinkageThinned -lmm 1 \
-k ./output/Hybs-LinkageThinned.cXX.txt \
-c /global/home/users/austinhpatton/fitnessGenomics/popFiles/GEMMA-Covariates.txt \
-o Hybs-LinkageThinned-Surv-LMM-Linear
