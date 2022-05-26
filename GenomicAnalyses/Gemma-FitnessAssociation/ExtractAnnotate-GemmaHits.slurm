#!/bin/bash
#SBATCH --output=get-fit-hits_%J.out
#SBATCH --mail-type=ALL
#SBATCH --time=72:00:00
#SBATCH --mail-user=austinhpatton@berkeley.edu
#SBATCH -p savio2
#SBATCH --account=ac_fishes
module load bio/vcftools/0.1.15
module load bcftools
module load java

## This script will:
#1) Pull out VCFs in hybrids that includes:
#    a) All significant GEMMA hits for growth (excluding non-survivors)
#    b) All significant GEMMA hits for composite fitness
#    c) Sites that are sourced from either standing genetic variation, de novo mutation, introgression or their combinations
#    (also recode these as 012 format for fitness landscapes)
#2) use snpEff to annotate these files, producing summary tables for each
#3) correspond the CBRO gene ID to gene card

# Get annotations
vcftools \
	--gzvcf /global/scratch/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/Hybs-Miss.85-maf05-DP7-LinkageThinned-Hybrids-FINAL.vcf.gz \
        --out /global/scratch/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/SiteAnnotations/Hybs-GrowthOnly-SigHits \
        --positions /global/home/users/austinhpatton/fitnessGenomics/SitesToAnnotate/GrowthOnly-AllSig-Sites.txt \
        --recode-INFO-all \
	--recode

vcftools \
	--gzvcf /global/scratch/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/Hybs-Miss.85-maf05-DP7-LinkageThinned-Hybrids-FINAL.vcf.gz \
        --out /global/scratch/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/SiteAnnotations/Hybs-Growth-SigHits \
        --positions /global/home/users/austinhpatton/fitnessGenomics/SitesToAnnotate/Growth-AllSig-Sites.txt \
        --recode-INFO-all \
        --recode

# Now, annotate each, and spit the output into their own summary csv
# Growth (without non-survivors)
java -Xmx8g -jar ~/software/snpEff/snpEff.jar Cbro -csvStats /global/scratch/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/SiteAnnotations/Hybs-GrowthOnly-SigHits \
/global/scratch/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/SiteAnnotations/Hybs-GrowthOnly-SigHits.recode.vcf > /global/scratch/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/SiteAnnotations/Hybs-GrowthOnly-SigHits.ann.vcf

# Composite fitness
java -Xmx8g -jar ~/software/snpEff/snpEff.jar Cbro -csvStats /global/scratch/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/SiteAnnotations/Hybs-Composite-SigHits \ 
/global/scratch/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/SiteAnnotations/Hybs-Growth-SigHits.recode.vcf > /global/scratch/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/SiteAnnotations/Hybs-Growth-SigHits.ann.vcf

# Pull out per-site summaries of effects, including the annotated genes
# First Growth (without non-survivors)
cat /global/scratch/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/SiteAnnotations/Hybs-GrowthOnly-SigHits.ann.vcf | ~/software/snpEff/scripts/vcfEffOnePerLine.pl | java -jar ~/software/snpEff/SnpSift.jar extractFields -s "," -e "." \
/global/scratch/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/SiteAnnotations/Hybs-GrowthOnly-SigHits.ann.vcf CHROM POS REF ALT "EFF[*].EFFECT" "EFF[*].AA" "EFF[*].IMPACT" "EFF[*].GENE" > \
/global/scratch/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/SiteAnnotations/Hybs-GrowthOnly-SigHits-SiteSumm.txt

# Then Composite
cat /global/scratch/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/SiteAnnotations/Hybs-Composite-SigHits.ann.vcf | ~/software/snpEff/scripts/vcfEffOnePerLine.pl | java -jar ~/software/snpEff/SnpSift.jar extractFields -s "," -e "." \
/global/scratch/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/SiteAnnotations/Hybs-Composite-SigHits.ann.vcf CHROM POS REF ALT "EFF[*].EFFECT" "EFF[*].AA" "EFF[*].IMPACT" "EFF[*].GENE" > \
/global/scratch/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/SiteAnnotations/Hybs-Composite-SigHits-SiteSumm.txt

# Now, duplicate and append the last column which is listed as 'CBRO' gene identifiers, but replace with gene cards
# Growth (excluding non-survivors)
cut -f8 /global/scratch/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/SiteAnnotations/Hybs-GrowthOnly-SigHits-SiteSumm.txt > tmp
grep -r -E -o ".{0,10}CBRO.{0,10}" tmp | sed 's/-//g' | sed 's/,//g' | sort | uniq > cbro.genes
while read gene; do grep $gene ~/software/snpEff/data/Cbro/putative-functions-final.txt; done < cbro.genes | \
sed 's/Protein of unknown function/UNKNOWN: Protein of unknown function/g' | sed "s/: /\t/g" | cut -f1,2 > tmp2 && mv tmp2 gene.cards
for i in $(cat cbro.genes); do card=$(grep $i gene.cards | cut -f2); sed -i "s/$i/$card/g" tmp; done
paste /global/scratch/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/SiteAnnotations/Hybs-GrowthOnly-SigHits-SiteSumm.txt tmp > tmp2 
mv tmp2 /global/scratch/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/SiteAnnotations/Hybs-GrowthOnly-SigHits-SiteSumm.txt
rm tmp && rm gene.cards && rm cbro.genes

# Composite
cut -f8 /global/scratch/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/SiteAnnotations/Hybs-Composite-SigHits-SiteSumm.txt > tmp
grep -r -E -o ".{0,10}CBRO.{0,10}" tmp | sed 's/-//g' | sed 's/,//g' | sort | uniq > cbro.genes
while read gene; do grep $gene ~/software/snpEff/data/Cbro/putative-functions-final.txt; done < cbro.genes | \
sed 's/Protein of unknown function/UNKNOWN: Protein of unknown function/g' | sed "s/: /\t/g" | cut -f1,2 > tmp2 && mv tmp2 gene.cards
for i in $(cat cbro.genes); do card=$(grep $i gene.cards | cut -f2); sed -i "s/$i/$card/g" tmp; done
paste /global/scratch/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/SiteAnnotations/Hybs-Composite-SigHits-SiteSumm.txt tmp > tmp2 
mv tmp2 /global/scratch/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/SiteAnnotations/Hybs-Composite-SigHits-SiteSumm.txt
rm tmp && rm gene.cards && rm cbro.genes
