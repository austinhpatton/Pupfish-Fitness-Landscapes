# Here we are relating morphological distance to genetic distance in sequenced hybrids
library(vcfR)
library(ggplot2)
library(reshape2)

setwd('~/Dropbox/Research/Martin-Berkeley-Postdoc/Pupfish-Fitness-Genomics/Results/')
samps <- read.table('./FixedSweeping/Hyb-Parent-Metadata.tsv', header = T)
samps <- samps[which(samps$Species == 'Hybrid'),]

# First using morphological data
morph <- read.table('../Data/Pupfish-final-trait-measurements-growth-admix.txt', header = T)
meta <- morph[,c(1:5,37)]
morph <- morph[,-c(1:5,37:39)]
rownames(morph) <- meta$ID

morph <- morph[which(rownames(morph) %in% samps$ID),]
meta <- meta[which(meta$ID %in% samps$ID),]
morph.dist <- dist(morph)

# then SNPs - this is easy
snps <- read.vcfR('./SNPs/CP-LL-HybParents-Miss.85-maf05-DP7-LinkageThinned-Hybrids-FINAL.vcf.gz')
snps <- vcfR2genlight(snps)

genet.dist <- dist(snps)
rm(snps)

genet.dist <- as.matrix(genet.dist)
genet.dist <- genet.dist[which(rownames(genet.dist) %in% as.character(samps$ID)), 
                         which(colnames(genet.dist) %in% as.character(samps$ID))]
genet.dist <- as.dist(genet.dist)

# Now reorder to make sure the two matrices match in order
genet.dist <- as.matrix(genet.dist)
morph.dist <- as.matrix(morph.dist)
row.idx <- match(rownames(morph.dist), rownames(genet.dist))
col.idx <- match(colnames(morph.dist), colnames(genet.dist))
genet.dist <- genet.dist[row.idx, col.idx]

genet.dist <- as.dist(genet.dist)
morph.dist <- as.dist(morph.dist)

      
dists <- data.frame('Morph' = c(morph.dist), 'Geno' = c(genet.dist))
dists <- melt(dist)

fit <- lm(Morph ~ Geno, data = dists)
dists.p <- 
  ggplot(dists, aes(x = Geno, y = Morph)) + 
  geom_point(alpha = 0.5, size = 1) +
  stat_smooth(method = lm) +
  ylab('Morphological Distance') +
  xlab('Genetic Distance') +
  annotate(geom = 'text', x = 450, y = 15, hjust = 0,
           label = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 3))) + 
  annotate(geom = 'text', x = 450, y = 14.5, hjust = 0,
           label = paste("Slope =",signif(fit$coef[[2]], 3))) +
  annotate(geom = 'text', x = 450, y = 14, hjust = 0,
           label = paste("P =",signif(summary(fit)$coef[2,4], 3))) +
  theme_classic(base_size = 18)

ggsave(dists.p, filename = 'Morphological-Vs-Genetic-Distance.pdf', height = 7, width = 7)


