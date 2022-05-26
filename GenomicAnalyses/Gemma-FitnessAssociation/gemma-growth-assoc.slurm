#!/bin/bash
#SBATCH --output=traits_%J.out
#SBATCH --mail-type=ALL
#SBATCH --time=72:00:00
#SBATCH --mail-user=austinhpatton@berkeley.edu
#SBATCH -p savio2
#SBATCH --account=ac_fishes

# This analysis will run gemma to test for associations between genotypes and growth.

module load r
module load gnu-parallel

# Run GEMMA analysis on survival
cd /global/scratch/austinhpatton/fitnessGenomics/results/GEMMA/traits

####################################################################################
# First estimate the covariance matrix
cd /global/scratch/austinhpatton/fitnessGenomics/results/GEMMA

PREFIX='/global/scratch/austinhpatton/fitnessGenomics/results/GEMMA/Hybs-LinkageThinned-GrowthOnly'
# Okay, use the general linkage-thinned SNPs
# First we need to convert the SNPs to the necessary PED format
/global/home/users/austinhpatton/software/plink2 \
--vcf /global/scratch/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/Hybs-Miss.85-maf05-DP7-LinkageThinned-Hybrids-FINAL.vcf.gz --out $PREFIX \
--allow-extra-chr \
--make-bed \
--keep /global/home/users/austinhpatton/fitnessGenomics/popFiles/Growth-Only-popfile.txt

# And replace the .fam file produced that is missing phenotype. 
cp ~/fitnessGenomics/GemmaTraits/growth-only.fam /global/scratch/austinhpatton/fitnessGenomics/results/GEMMA/Hybs-LinkageThinned-GrowthOnly.fam

# Compute kinship matrix
/global/home/users/austinhpatton/software/gemma-0.98.1 -bfile ${PREFIX} -gk -o Hybs-LinkageThinned-GrowthOnly


# Alright, now run GEMMA
PREFIX='/global/scratch/austinhpatton/fitnessGenomics/results/GEMMA/traits'
TRAITPATH='/global/home/users/austinhpatton/fitnessGenomics/GemmaTraits'

cd /global/scratch/austinhpatton/fitnessGenomics/results/GEMMA/traits

cp /global/scratch/austinhpatton/fitnessGenomics/results/GEMMA/Hybs-LinkageThinned-GrowthOnly.bed growth-only.bed
cp /global/scratch/austinhpatton/fitnessGenomics/results/GEMMA/Hybs-LinkageThinned-GrowthOnly.bim growth-only.bim
cp ${TRAITPATH}/growth-only.fam  growth-only.fam

/global/home/users/austinhpatton/software/gemma-0.98.1 -bfile ${PREFIX}/growth-only -lmm 1 \
-k /global/scratch/austinhpatton/fitnessGenomics/results/GEMMA/output/Hybs-LinkageThinned-GrowthOnly.cXX.txt \
-c /global/home/users/austinhpatton/fitnessGenomics/popFiles/GEMMA-Covariates-GrowthOnly.txt \
-o Hybs-LinkageThinned-growth-only-LMM-Linear

rm growth-only.bed
rm growth-only.bim
rm growth-only.fam
    
####################################################################################
