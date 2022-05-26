#!/bin/bash
#SBATCH --output=thinPCA_%J.out
#SBATCH --mail-type=ALL
#SBATCH --time=48:00:00
#SBATCH --mail-user=austinhpatton@berkeley.edu
#SBATCH -p savio2
#SBATCH --account=ac_fishes
module load bcftools

cd /global/scratch/austinhpatton/fitnessGenomics/results/PCA

VCF='/global/scratch/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/Hybs-Miss.85-maf05-DP7.reheadered.vcf.gz'

#First thin the sites so as to reduce linkage
/global/home/users/austinhpatton/software/plink-1.9 --noweb \
--vcf $VCF \
--double-id \
--allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 10['kb'] 50 0.5 \
--out /global/scratch/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/Hybs-Miss.85-maf05-DP7-LinkageThinned

/global/home/users/austinhpatton/software/plink-1.9 --noweb \
--vcf $VCF \
--recode vcf \
--double-id \
--extract /global/scratch/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/Hybs-Miss.85-maf05-DP7-LinkageThinned.prune.in \
--allow-extra-chr \
--set-missing-var-ids @:# \
--out /global/scratch/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/Hybs-Miss.85-maf05-DP7-LinkageThinned

/global/home/users/austinhpatton/software/htslib-1.11/bgzip \
/global/scratch/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/Hybs-Miss.85-maf05-DP7-LinkageThinned.vcf

# Then, perform a PCA
/global/home/users/austinhpatton/software/plink-1.9 \
--vcf /global/scratch/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/Hybs-Miss.85-maf05-DP7-LinkageThinned.vcf.gz \
--double-id \
--allow-extra-chr \ 
--set-missing-var-ids @:# \
--make-bed \
--pca \
--out /global/scratch/austinhpatton/fitnessGenomics/results/PCA/AllSamp-PCA

# Remove double IDs
bcftools query -l /global/scratch/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/Hybs-Miss.85-maf05-DP7-LinkageThinned.vcf.gz > tmp
cut -d'_' -f2 tmp > tmp2
paste tmp tmp2 > plink-reheader

# Fix the headers
bcftools reheader -s /global/home/users/austinhpatton/fitnessGenomics/scripts/plink-reheader.txt \
-o /global/scratch/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/Hybs-Parents-Miss.85-maf05-DP7-LinkageThinned.FINAL.vcf.gz \
/global/scratch/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/Hybs-Miss.85-maf05-DP7-LinkageThinned.vcf.gz

####

# Pull out only the hybrids from these thinned data
bcftools view -O z -o /global/scratch/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/Hybs-Miss.85-maf05-DP7-LinkageThinned-Hybrids-FINAL.vcf.gz \
--samples-file ~/fitnessGenomics/popFiles/Surv-NonSurv-pop_file.txt /global/scratch/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/Hybs-Parents-Miss.85-maf05-DP7-LinkageThinned.FINAL.vcf.gz 

#And again, perform a PCA
/global/home/users/austinhpatton/software/plink-1.9 \
--vcf /global/scratch/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/Hybs-Miss.85-maf05-DP7-LinkageThinned-Hybrids-FINAL.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--make-bed \
--pca \
--out /global/scratch/austinhpatton/fitnessGenomics/results/PCA/Hybrid-PCA

################################################################################
# Repeat the process using just hybrids and parents from CP and LL
# Pull out only the hybrids from these thinned data
bcftools view -O z -o /global/scratch/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/CP-LL-HybParents-Miss.85-maf05-DP7-LinkageThinned-Hybrids-FINAL.vcf.gz \
--samples-file ~/fitnessGenomics/popFiles/CP-LL-HybParents.txt /global/scratch/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/Hybs-Parents-Miss.85-maf05-DP7-LinkageThinned.FINAL.vcf.gz 

#And again, perform a PCA
/global/home/users/austinhpatton/software/plink-1.9 \
--vcf /global/scratch/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/CP-LL-HybParents-Miss.85-maf05-DP7-LinkageThinned-Hybrids-FINAL.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--make-bed \
--pca \
--out /global/scratch/austinhpatton/fitnessGenomics/results/PCA/CP-LL-HybParents-PCA


# Also just pull out hybrids and parents from CP and LL, but from un-thinned data
bcftools view -O z -o /global/scratch/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/CP-LL-HybParents-Miss.85-maf05-DP7-Unthinned-Hybrids-FINAL.vcf.gz \
--samples-file ~/fitnessGenomics/popFiles/CP-LL-HybParents.txt ${VCF}
