#!/bin/bash
#SBATCH --output=gatk-SNP_%J.out
#SBATCH --mail-type=ALL
#SBATCH --time=72:00:00
#SBATCH --mail-user=austinhpatton@berkeley.edu
#SBATCH -p savio2
#SBATCH --account=fc_fishes
module load java

# Pull out RAW indels 
java -jar ~/software/GenomeAnalysisTK.jar -T SelectVariants \
    -V /global/scratch/austinhpatton/fitnessGenomics/genomes/GATK-RAW/pupfishFitness_RTM1-2_Carib250_assembly-sorted.vcf.gz \
    -R /global/scratch/ejrichards/Cyprinodon_ref/jeremy_align/asm.racon.fasta \
    -o /global/scratch/austinhpatton/fitnessGenomics/genomes/GATK-RAW/pupfishFitness_RTM1-2_Carib-FinalSamps-RAW-SNPS.vcf.gz \
    -selectType SNP \
    --exclude_sample_name CPH125 \
    --exclude_sample_name CPH06 

java -jar ~/software/GenomeAnalysisTK.jar -T VariantFiltration \
    -V /global/scratch/austinhpatton/fitnessGenomics/genomes/GATK-RAW/pupfishFitness_RTM1-2_Carib-FinalSamps-RAW-SNPS.vcf.gz \
    -R /global/scratch/ejrichards/Cyprinodon_ref/jeremy_align/asm.racon.fasta \
    -o /global/scratch/austinhpatton/fitnessGenomics/genomes/GATK-FINAL/pupfishFitness_RTM1-2_Carib-FinalSamps-FINAL-SNPS.vcf.gz \
    -filter "QD < 2.0" -filterName "QD2" \
    -filter "QUAL < 20.0" -filterName "QUAL20" \
    -filter "FS > 60.0" -filterName "FS60" \
    -filter "MQRankSum < -12.5" -filterName "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" -filterName "ReadPosRankSum-8"
