#!/bin/bash
#SBATCH --output=sup_admix_%J.out
#SBATCH --error=sup_admix_%J.err
#SBATCH --mail-type=ALL
#SBATCH --time=72:00:00
#SBATCH --mail-user=austinhpatton@berkeley.edu
#SBATCH -p savio
#SBATCH --account=fc_fishes

#################################################################################################
module load bcftools

#################################################################################################
# Pull out the fish to be used in LL analysis. Generalist and specialist parents from LL and 
# Osprey (the two are a single body of water), and LL hybrids. 

bcftools view -O z -o \
/global/scratch/users/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/LL-Admix.vcf.gz \ 
--samples-file ~/fitnessGenomics/popFiles/LL-Admixture-set.txt --force-samples \ 
/global/scratch/users/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/Hybs-Parents-Miss.85-maf05-DP7-LinkageThinned.FINAL.vcf.gz

FILE='/global/scratch/users/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/LL-Admix'

# Generate the input file in plink format
/global/scratch/users/austinhpatton/software/plink-1.9 \
--vcf ${FILE}.vcf.gz --make-bed --out $FILE --allow-extra-chr

# ADMIXTURE does not accept chromosome names that are not human chromosomes. We will thus just 
# exchange the first column by 0
awk '{$1=0;print $0}' $FILE.bim > $FILE.bim.tmp
mv $FILE.bim.tmp $FILE.bim

#################################################################################################
# Do the same for CP. Pull out the fish to be used in CP analysis. Generalist and specialist 
# parents from CP and CP hybrids. 

bcftools view -O z -o \ 
/global/scratch/users/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/CP-Admix.vcf.gz \
--samples-file ~/fitnessGenomics/popFiles/CP-Admixture-set.txt --force-samples \
/global/scratch/users/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/Hybs-Parents-Miss.85-maf05-DP7-LinkageThinned.FINAL.vcf.gz

FILE='/global/scratch/users/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/CP-Admix'

# Generate the input file in plink format
/global/scratch/users/austinhpatton/software/plink-1.9 \
--vcf ${FILE}.vcf.gz --make-bed --out $FILE --allow-extra-chr

# ADMIXTURE does not accept chromosome names that are not human chromosomes. We will thus just
# exchange the first column by 0
awk '{$1=0;print $0}' $FILE.bim > $FILE.bim.tmp
mv $FILE.bim.tmp $FILE.bim

#################################################################################################
# And now run the analyses, unsupervised, assuming K=3. 
LL='/global/scratch/users/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/LL-Admix.bed'
/global/scratch/users/austinhpatton/software/admixture_linux-1.3.0/admixture ${LL} 3 -j24 > admixture-LL-K3.log

CP='/global/scratch/users/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/CP-Admix.bed'
/global/scratch/users/austinhpatton/software/admixture_linux-1.3.0/admixture ${CP} 3 -j24 > admixture-CP-K3.log

grep "CV" admixture-K*.log | awk '{print $3,$4}' | cut -c 4,7-20 > Admixture-HybParents.cv.error

#################################################################################################
# Repeat the analyses, this time supervised, assuming K=3.
LL='/global/scratch/users/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/LL-Admix.bed'
/global/scratch/users/austinhpatton/software/admixture_linux-1.3.0/admixture ${LL} 3 -j30 --supervised > admixture-LL-Supervised-K3.log

CP='/global/scratch/users/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/CP-Admix.bed'
/global/scratch/users/austinhpatton/software/admixture_linux-1.3.0/admixture ${CP} 3 -j30 --supervised > admixture-CP-Supervised-K3.log


#################################################################################################
# Now do these again, but including only the samples from RTM2
#################################################################################################

# First for little lake/osprey lakes
bcftools view -O z -o \
/global/scratch/users/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/LL-RTM2-Admix.vcf.gz \
--samples-file ~/fitnessGenomics/popFiles/LL-RTM2-Admixture-set.txt --force-samples \
/global/scratch/users/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/Hybs-Parents-Miss.85-maf05-DP7-LinkageThinned.FINAL.vcf.gz

FILE='/global/scratch/users/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/LL-RTM2-Admix'

# Generate the input file in plink format
/global/scratch/users/austinhpatton/software/plink-1.9 \
--vcf ${FILE}.vcf.gz --make-bed --out $FILE --allow-extra-chr

# ADMIXTURE does not accept chromosome names that are not human chromosomes. We will thus just
# exchange the first column by 0
awk '{$1=0;print $0}' $FILE.bim > $FILE.bim.tmp
mv $FILE.bim.tmp $FILE.bim

###############################################################

# then for crescent pond
bcftools view -O z -o \
/global/scratch/users/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/CP-RTM2-Admix.vcf.gz \
--samples-file ~/fitnessGenomics/popFiles/CP-RTM2-Admixture-set.txt --force-samples \
/global/scratch/users/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/Hybs-Parents-Miss.85-maf05-DP7-LinkageThinned.FINAL.vcf.gz

FILE='/global/scratch/users/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/CP-RTM2-Admix'

# Generate the input file in plink format
/global/scratch/users/austinhpatton/software/plink-1.9 \
--vcf ${FILE}.vcf.gz --make-bed --out $FILE --allow-extra-chr

# ADMIXTURE does not accept chromosome names that are not human chromosomes. We will thus just
# exchange the first column by 0
awk '{$1=0;print $0}' $FILE.bim > $FILE.bim.tmp
mv $FILE.bim.tmp $FILE.bim

###############################################################
# First the unsupervised analyses
LL='/global/scratch/users/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/LL-RTM2-Admix.bed'
/global/scratch/users/austinhpatton/software/admixture_linux-1.3.0/admixture ${LL} 3 -j30 > admixture-LL-RTM2-K3.log

CP='/global/scratch/users/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/CP-RTM2-Admix.bed'
/global/scratch/users/austinhpatton/software/admixture_linux-1.3.0/admixture ${CP} 3 -j30 > admixture-CP-RTM2-K3.log

###############################################################
# Then the upervised analyses
LL='/global/scratch/users/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/LL-RTM2-Admix.bed'
/global/scratch/users/austinhpatton/software/admixture_linux-1.3.0/admixture ${LL} 3 -j30 --supervised > admixture-LL-RTM2-Supervised-K3.log

CP='/global/scratch/users/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/CP-RTM2-Admix.bed'
/global/scratch/users/austinhpatton/software/admixture_linux-1.3.0/admixture ${CP} 3 -j30 --supervised > admixture-CP-RTM2-Supervised-K3.log
