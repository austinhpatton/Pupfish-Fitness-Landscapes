#!/bin/bash
#SBATCH --output=sweeping012_%J.out
#SBATCH --mail-type=ALL
#SBATCH --time=72:00:00
#SBATCH --mail-user=austinhpatton@berkeley.edu
#SBATCH -p savio3
#SBATCH --account=co_fishes
module load bio/vcftools/0.1.15
module load bcftools

vcftools \
	--gzvcf /global/scratch/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/Hybs-Parents-Miss.85-maf05-DP7-LinkageThinned.FINAL.vcf.gz \
  --out /global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/sgv-sites \
	--keep /global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/hyb-parents-keep.txt \
  --positions /global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/sgv-sites.txt \
  --012

sed -i 's/-1/NA/g' /global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/sgv-sites.012
sed 's/\t/:/g' /global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/sgv-sites.012.pos | awk 'BEGIN { ORS = " " } { print }' > tmp
echo ID $(cat tmp) > tmp2
awk 'FNR==NR{a[NR]=$1;next}{$1=a[FNR]}1' /global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/sgv-sites.012.indv /global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/sgv-sites.012 > tmp 
cat tmp2 tmp >/global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/sgv-sites.012-mat.txt
rm tmp && rm tmp2 

vcftools \
	--gzvcf /global/scratch/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/Hybs-Parents-Miss.85-maf05-DP7-LinkageThinned.FINAL.vcf.gz \
  --out /global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/sgv-intro-sites \
  --keep /global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/hyb-parents-keep.txt \
  --positions /global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/sgv-intro-sites.txt \
  --012

sed -i 's/-1/NA/g' /global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/sgv-intro-sites.012
sed 's/\t/:/g' /global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/sgv-intro-sites.012.pos | awk 'BEGIN { ORS = " " } { print }' > tmp
echo ID $(cat tmp) > tmp2
awk 'FNR==NR{a[NR]=$1;next}{$1=a[FNR]}1' /global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/sgv-intro-sites.012.indv /global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/sgv-intro-sites.012 > tmp
cat tmp2 tmp >/global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/sgv-intro-sites.012-mat.txt
rm tmp && rm tmp2 

vcftools \
	--gzvcf /global/scratch/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/Hybs-Parents-Miss.85-maf05-DP7-LinkageThinned.FINAL.vcf.gz \
  --out /global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/sgv-denovo-sites \
  --keep /global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/hyb-parents-keep.txt \
  --positions /global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/sgv-denovo-sites.txt \
  --012

sed -i 's/-1/NA/g' /global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/sgv-denovo-sites.012
sed 's/\t/:/g' /global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/sgv-denovo-sites.012.pos | awk 'BEGIN { ORS = " " } { print }' > tmp
echo ID $(cat tmp) > tmp2
awk 'FNR==NR{a[NR]=$1;next}{$1=a[FNR]}1' /global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/sgv-denovo-sites.012.indv /global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/sgv-denovo-sites.012 > tmp
cat tmp2 tmp >/global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/sgv-denovo-sites.012-mat.txt
rm tmp && rm tmp2 

vcftools \
	--gzvcf /global/scratch/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/Hybs-Parents-Miss.85-maf05-DP7-LinkageThinned.FINAL.vcf.gz \
  --out /global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/sgv-intro-denovo-sites \
  --keep /global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/hyb-parents-keep.txt \
  --positions /global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/sgv-intro-denovo-sites.txt \
  --012

sed -i 's/-1/NA/g' /global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/sgv-intro-denovo-sites.012
sed 's/\t/:/g' /global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/sgv-intro-denovo-sites.012.pos | awk 'BEGIN { ORS = " " } { print }' > tmp
echo ID $(cat tmp) > tmp2
awk 'FNR==NR{a[NR]=$1;next}{$1=a[FNR]}1' /global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/sgv-intro-denovo-sites.012.indv /global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/sgv-intro-denovo-sites.012 > tmp
cat tmp2 tmp >/global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/sgv-intro-denovo-sites.012-mat.txt#rm tmp && rm tmp2 

vcftools \
	--gzvcf /global/scratch/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/Hybs-Parents-Miss.85-maf05-DP7-LinkageThinned.FINAL.vcf.gz \
  --out /global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/denovo-sites \
  --keep /global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/hyb-parents-keep.txt \
  --positions /global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/denovo-sites.txt \
  --012

sed -i 's/-1/NA/g' /global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/denovo-sites.012
sed 's/\t/:/g' /global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/denovo-sites.012.pos | awk 'BEGIN { ORS = " " } { print }' > tmp
echo ID $(cat tmp) > tmp2
awk 'FNR==NR{a[NR]=$1;next}{$1=a[FNR]}1' /global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/denovo-sites.012.indv /global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/denovo-sites.012 > tmp
cat tmp2 tmp >/global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/denovo-sites.012-mat.txt
rm tmp && rm tmp2 

vcftools \
	--gzvcf /global/scratch/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/Hybs-Parents-Miss.85-maf05-DP7-LinkageThinned.FINAL.vcf.gz \
  --out /global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/intro-sites \
  --keep /global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/hyb-parents-keep.txt \
  --positions /global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/intro-sites.txt \
  --012

sed -i 's/-1/NA/g' /global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/intro-sites.012
sed 's/\t/:/g' /global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/intro-sites.012.pos | awk 'BEGIN { ORS = " " } { print }' > tmp
echo ID $(cat tmp) > tmp2
awk 'FNR==NR{a[NR]=$1;next}{$1=a[FNR]}1' /global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/intro-sites.012.indv /global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/intro-sites.012 > tmp
cat tmp2 tmp >/global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/intro-sites.012-mat.txt
rm tmp && rm tmp2 

vcftools \
	--gzvcf /global/scratch/austinhpatton/fitnessGenomics/genomes/FilteredSNPs/Hybs-Parents-Miss.85-maf05-DP7-LinkageThinned.FINAL.vcf.gz \
  --out /global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/denovo-intro-sites \
  --keep /global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/hyb-parents-keep.txt \
  --positions /global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/denovo-intro-sites.txt \
  --012

sed -i 's/-1/NA/g' /global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/denovo-intro-sites.012
sed 's/\t/:/g' /global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/denovo-intro-sites.012.pos | awk 'BEGIN { ORS = " " } { print }' > tmp
echo ID $(cat tmp) > tmp2
awk 'FNR==NR{a[NR]=$1;next}{$1=a[FNR]}1' /global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/denovo-intro-sites.012.indv /global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/denovo-intro-sites.012 > tmp
cat tmp2 tmp >/global/scratch/austinhpatton/fitnessGenomics/results/SweepingSNPs/denovo-intro-sites.012-mat.txt
rm tmp && rm tmp2 


