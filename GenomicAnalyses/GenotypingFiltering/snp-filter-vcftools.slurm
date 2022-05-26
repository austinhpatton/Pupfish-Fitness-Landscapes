#!/bin/bash
#SBATCH --output=filterSnps_%J.out
#SBATCH --mail-type=ALL
#SBATCH --time=72:00:00
#SBATCH --mail-user=austinhpatton@berkeley.edu
#SBATCH -p savio2
#SBATCH --account=ac_fishes
module load bio/vcftools/0.1.15

VCF='/global/scratch/austinhpatton/fitnessGenomics/genomes/GATK-FINAL/pupfishFitness_RTM1-2_Carib-FinalSamps-FINAL-SNPS.vcf.gz'
maf5dp7='/global/scratch/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/Hybs-Miss.85-maf05-DP7'
dp7='/global/scratch/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/Hybs-Miss.85-DP7'

vcftools \
	--gzvcf $VCF \
	--out $maf5dp7 \
	--maf 0.05 \
	--min-alleles 2 \
	--max-alleles 2 \
	--min-meanDP 7 \
	--max-meanDP 100 \
	--max-missing 0.85 \
	--remove-indels \
  --remove-indv LLHF011.2 \
	--remove-indv LLLF002 \
	--remove-indv LLHF003 \
  --remove-indv LLHF026 \
  --remove-indv CP06H05 \
  --remove-indv CP04B10 \
  --remove-indv CP07B03 \
	--remove-indv CP11C07 \
  --remove-indv CP14D08 \
  --remove-indv CP17A08 \
  --remove-indv CP17H10 \
	--remove-indv CP20A04 \
  --remove-indv APHFAS \
  --remove-indv CUATESS1 \
  --remove-indv MEG-Q1 \
  --remove-indv NCCA10-2 \
  --remove-indv NCCA7-2 \
  --remove-indv NCCA8-2 \
	--remove-indv BAVA10-2 \
	--remove-indv VENA1-2 \
	--recode-INFO-all \
	--recode

/global/home/users/austinhpatton/software/htslib-1.11/bgzip /global/scratch/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/Hybs-Miss.85-maf05-DP7.recode.vcf
bcftools reheader -s SampRename.txt -o /global/scratch/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/Hybs-Miss.85-maf05-DP7.reheadered.vcf.gz /global/scratch/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/Hybs-Miss.85-maf05-DP7.recode.vcf.gz
