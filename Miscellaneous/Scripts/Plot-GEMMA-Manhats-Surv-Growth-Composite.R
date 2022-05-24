# Analyze and present a simple LMM in GEMMA
library(tidyverse)
library(data.table)
library(ggplot2)
library(progressr)
library(foreach)
library(wesanderson)
library(cowplot)
library(ggrepel)

setwd('~/Dropbox/Research/Martin-Berkeley-Postdoc/Pupfish-Fitness-Genomics/Results/GEMMA/Survival/output/')

res <- read_tsv('./Hybs-LinkageThinned-Surv-LMM-Linear.assoc.txt')
res <- res[which(res$af >= 0.05),]

# Obtain p/q-values for plotting
res$log10_q <- -log10(p.adjust(res$p_wald, method = 'fdr'))
res$p <- -log10(res$p_wald)

res$RelPosition <- cumsum(res$ps)

cols <- c('grey0', 'grey44', '#2b83ba', '#d7191c')

# read in the file
scaffLen <- read.table('~/Dropbox/Research/Martin-Berkeley-Postdoc/Pupfish-Fitness-Genomics/Data/Cyprinodon_ref/Scaff-Lengths.txt')

colnames(scaffLen) <- c('chr', 'length')
scaffLen <- scaffLen[order(scaffLen$length, decreasing = T),]

scaffLen$CumPos <- c(0, cumsum(scaffLen$length)[-length(scaffLen$length)])

params <- data.frame('chr'= NA, 'rs' = NA, 'ps' = NA, 'n_miss' = NA,
                     'allele1' = NA, 'allele0' = NA, 'af' = NA, 
                     'beta' = NA, 'se' = NA, 'logl_H1' = NA,
                     'l_remle' = NA, 'p_wald' = NA, 'log10_q' = NA,
                     'RelPosition' = NA, 'ScaffCol' = NA)

get.relPos.1 <- 
  function(scaffLen, res, chr){
    tmp <- res[which(res$chr == chr),]
    tmp <- merge(scaffLen[,c(1,3)], tmp, by = 'chr')
    tmp$RelPosition <- tmp$ps + tmp$CumPos
    tmp$ScaffCol <- rep('grey0', nrow(tmp))
    return(tmp)
  }

get.relPos.2 <- 
  function(scaffLen, res, chr){
    tmp <- res[which(res$chr == chr),]
    tmp <- merge(scaffLen[,c(1,3)], tmp, by = 'chr')
    tmp$RelPosition <- tmp$ps + tmp$CumPos
    tmp$ScaffCol <- rep('grey44', nrow(tmp))
    return(tmp)
  }


with_progress({
  p <- progressor(along = 1:length(scaffLen$chr))
  lmm.res <- 
    foreach(i = 1:length(scaffLen$chr), .combine = 'rbind') %do% {
      p(sprintf("x=%g", i))
      if(i%%2 == 1){
        get.relPos.1(scaffLen, res, chr = scaffLen$chr[i])
      }else{
        get.relPos.2(scaffLen, res, chr = scaffLen$chr[i])
      }
    }
})

# cols <- RColorBrewer::brewer.pal(9, c('YlOrRd'))
# newcol <- colorRampPalette(cols[2:7])
# ncols <- 10
# cols <- newcol(ncols)

# write.table(lmm.res, file = 'GEMMA-LMM-1per-10Kb-RESULTS.txt', sep = '\t', quote = F, row.names = F)
lmm.res <- read_tsv('GEMMA-LMM-1per-10Kb-RESULTS.txt')

# Recolor outliers as blue (FDR) or red (bonferroni)
lmm.res[which(lmm.res$log10_q > 1.30103),17] <- 'FDR-sig'
lmm.res[which(lmm.res$p >= -log10(0.05 / nrow(lmm.res))),17] <- 'Bonf-sig'
bonf.p <- -log10(0.05 / nrow(lmm.res))

to.drop <- which(lmm.res$ScaffCol %in% c('grey0', 'grey44') & lmm.res$p < 3.5)
lmm.res <- lmm.res[-to.drop[round(seq(1, length(to.drop), 1.25))],]

lmm.res$ScaffCol <- factor(lmm.res$ScaffCol, 
                           levels = c('grey0', 'grey44',
                                      'FDR-sig', 'Bonf-sig'))

y.axis <- expression(paste("-Log10 (", italic(p), '-value)'))

lmm.res[-which(lmm.res$chr %in% unique(lmm.res$chr)[1:24] | lmm.res$ScaffCol %in% c('FDR-sig', 'Bonf-sig')),] 
lmm.res <- mutate(lmm.res, ScaffCol = ifelse(test = ScaffCol %in% c('grey0', 'grey44') & chr %in% unique(lmm.res$chr)[-c(1:24)], 
                                             yes = as.character('grey68'), no = as.character(ScaffCol))) 
lmm.res$ScaffCol <- as.factor(lmm.res$ScaffCol)



# #Now, thin non-outliers since it's crazy to plot otherwise
# to.drop <- which(lmm.res$ScaffCol %in% c('grey0', 'grey44') & lmm.res$p < 3.5)
# lmm.res <- lmm.res[-to.drop[round(seq(1, length(to.drop), 1.25))],]

surv.manhat <-
  ggplot(lmm.res, aes(RelPosition, p, colour = ScaffCol, size = ScaffCol)) + 
  geom_hline(yintercept = -log10(0.05 / nrow(lmm.res)), color = '#d7191c', 
             lty = 2, alpha = 0.5, size = 1.25) +
  geom_point(alpha = 0.75) +
  scale_color_manual(values = c('grey0', 'grey44', 'grey68')) +
  scale_size_manual(values = c(1.5,1.5,1.5)) +
  labs(title = 'Survival') +
  xlab('Relative Position') + 
  ylab(y.axis) +
  theme_light(base_size = 24) +
  scale_x_continuous(breaks = NULL) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = 'none')

png('Gemma-LMM-Bonf-LinkageThinned-survivorship.png', width = 12, height = 4, units = 'in', res = 600)
surv.manhat
dev.off()

pdf('Gemma-LMM-Bonf-LinkageThinned-survivorship.pdf', width = 12, height = 4)
surv.manhat
dev.off()

######################################################################################################
# Plot with data for composite fitness

res <- read_tsv('../../Traits/output/Hybs-LinkageThinned-growth-LMM-Linear.assoc.txt')
res <- res[which(res$af >= 0.05),]

# Obtain p/q-values for plotting
res$log10_q <- -log10(p.adjust(res$p_wald, method = 'fdr'))
res$p <- -log10(res$p_wald)

res$RelPosition <- cumsum(res$ps)

cols <- c('grey0', 'grey44', '#2b83ba', '#d7191c')

# read in the file
scaffLen <- read.table('~/Dropbox/Research/Martin-Berkeley-Postdoc/Pupfish-Fitness-Genomics/Data/Cyprinodon_ref/Scaff-Lengths.txt')

colnames(scaffLen) <- c('chr', 'length')
scaffLen <- scaffLen[order(scaffLen$length, decreasing = T),]

scaffLen$CumPos <- c(0, cumsum(scaffLen$length)[-length(scaffLen$length)])

# Rescale position by relative position; 
#res <- merge(res, refPos[c(3,16,17)], by = 'rs')
with_progress({
  p <- progressor(along = 1:length(scaffLen$chr))
  growth.res <- 
    foreach(i = 1:length(scaffLen$chr), .combine = 'rbind') %do% {
      p(sprintf("x=%g", i))
      if(i%%2 == 1){
        get.relPos.1(scaffLen, res, chr = scaffLen$chr[i])
      }else{
        get.relPos.2(scaffLen, res, chr = scaffLen$chr[i])
      }
    }
})

res <- growth.res

#write.table(res, file = '../../Traits/output/Composite-GEMMA-LMM-1per-10Kb-RESULTS.txt', sep = '\t', quote = F, row.names = F)
comp.res <- read_tsv('../../Traits/output/Composite-GEMMA-LMM-1per-10Kb-RESULTS.txt')

# Recolor outliers as blue (FDR) or red (bonferroni)
comp.res[which(comp.res$log10_q > 1.30103),17] <- 'FDR-sig'
comp.res[which(comp.res$p >= -log10(0.05 / nrow(comp.res))),17] <- 'Bonf-sig'
bonf.p <- -log10(0.05 / nrow(comp.res))

to.drop <- which(comp.res$ScaffCol %in% c('grey0', 'grey44') & comp.res$p < 3.5)
comp.res <- comp.res[-to.drop[round(seq(1, length(to.drop), 1.25))],]

comp.res$ScaffCol <- factor(comp.res$ScaffCol, 
                       levels = c('grey0', 'grey44',
                                  'FDR-sig', 'Bonf-sig'))

y.axis <- expression(paste("-Log10 (", italic(p), '-value)'))

# comp.annots <- 
#   data.frame('Gene' = c('CSAD', 'CSAD', 'CSAD', 'GLCCI1', 'INO80C', 
#                         'INO80C', 'INO80C', 'MAG', 'PIM2', 'MED25'),
#              'RelPosition' = c(677459848, 677460326, 677459697,
#                                298797665, 286230740, 286230787,
#                                286231015, 282928254, 561324138,
#                                337515010),
#              'p' = c(6.66, 6.82, 6.0, 5.18, 5.51, 
#                      5.66, 5.73, 5.18, 6.50, 5.35),
#              'ScaffCol' = rep('grey0', 10))

comp.annots <- 
  data.frame('Gene' = c('CSAD', 'GLCCI1', 'INO80C', 
                        'MAG', 'PIM2', 'MED25'),
             'RelPosition' = c(677460326, 298797665, 286231015, 
                               282928254, 561324138, 337515010),
             'p' = c(6.82, 5.18, 5.73, 
                     5.18, 6.50, 5.35),
             'ScaffCol' = rep('grey0', 6))

comp.res[-which(comp.res$chr %in% unique(comp.res$chr)[1:24] | comp.res$ScaffCol %in% c('FDR-sig', 'Bonf-sig')),] 
comp.res <- mutate(comp.res, ScaffCol = ifelse(test = ScaffCol %in% c('grey0', 'grey44') & chr %in% unique(comp.res$chr)[-c(1:24)], 
                                             yes = as.character('grey68'), no = as.character(ScaffCol))) 
comp.res$ScaffCol <- as.factor(comp.res$ScaffCol)

comp.manhat <-
  ggplot(comp.res, aes(RelPosition, p, colour = ScaffCol, size = ScaffCol)) + 
  geom_hline(yintercept = bonf.p, color = '#d7191c', 
             lty = 2, alpha = 0.5, size = 1.25) +
  geom_point(alpha = 0.75) +
  scale_color_manual(values = c('#d7191c', '#2b83ba', 'grey0', 'grey44', 'grey68')) +
  scale_size_manual(values = c(4.5,3,1.5,1.5,1.5)) +
  labs(title = 'Composite Fitness') +
  xlab('Relative Position') + 
  ylab(y.axis) +
  theme_light(base_size = 24) +
  scale_x_continuous(breaks = NULL) +
  geom_label_repel(data = comp.annots, aes(label = Gene), 
                   box.padding = 1, size = 4, 
                   fill = 'white') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = 'none')
comp.manhat

png('../../Traits/output/Gemma-LMM-Bonf-LinkageThinned-composite.png', width = 12, height = 4, units = 'in', res = 600)
comp.manhat
dev.off()

pdf('../../Traits/output/Gemma-LMM-Bonf-LinkageThinned-composite.pdf', width = 12, height = 4)
comp.manhat
dev.off()






######################################################################################################
# Plot with data for growth

res <- read_tsv('../../Traits/output/Hybs-LinkageThinned-growth-only-LMM-Linear.assoc.txt')
res <- res[which(res$af >= 0.05),]

# Obtain p/q-values for plotting
res$log10_q <- -log10(p.adjust(res$p_wald, method = 'fdr'))
res$p <- -log10(res$p_wald)

res$RelPosition <- cumsum(res$ps)

cols <- c('grey0', 'grey44', '#2b83ba', '#d7191c')

# read in the file
scaffLen <- read.table('~/Dropbox/Research/Martin-Berkeley-Postdoc/Pupfish-Fitness-Genomics/Data/Cyprinodon_ref/Scaff-Lengths.txt')

colnames(scaffLen) <- c('chr', 'length')
scaffLen <- scaffLen[order(scaffLen$length, decreasing = T),]

scaffLen$CumPos <- c(0, cumsum(scaffLen$length)[-length(scaffLen$length)])

# Rescale position by relative position; 
#res <- merge(res, refPos[c(3,16,17)], by = 'rs')
with_progress({
  p <- progressor(along = 1:length(scaffLen$chr))
  growth.res <- 
    foreach(i = 1:length(scaffLen$chr), .combine = 'rbind') %do% {
      p(sprintf("x=%g", i))
      if(i%%2 == 1){
        get.relPos.1(scaffLen, res, chr = scaffLen$chr[i])
      }else{
        get.relPos.2(scaffLen, res, chr = scaffLen$chr[i])
      }
    }
})

res <- growth.res

#write.table(res, file = '../../Traits/output/Growth-GEMMA-LMM-1per-10Kb-RESULTS.txt', sep = '\t', quote = F, row.names = F)
grow.res <- read_tsv('../../Traits/output/Growth-GEMMA-LMM-1per-10Kb-RESULTS.txt')

# Recolor outliers as blue (FDR) or red (bonferroni)
grow.res[which(grow.res$log10_q > 1.30103),17] <- 'FDR-sig'
grow.res[which(grow.res$p >= -log10(0.05 / nrow(grow.res))),17] <- 'Bonf-sig'
bonf.p <- -log10(0.05 / nrow(grow.res))

to.drop <- which(grow.res$ScaffCol %in% c('grey0', 'grey44') & grow.res$p < 3.5)
grow.res <- grow.res[-to.drop[round(seq(1, length(to.drop), 1.25))],]

grow.res$ScaffCol <- factor(grow.res$ScaffCol, 
                       levels = c('grey0', 'grey44',
                                  'FDR-sig', 'Bonf-sig'))

y.axis <- expression(paste("-Log10 (", italic(p), '-value)'))

grow.res[-which(grow.res$chr %in% unique(grow.res$chr)[1:24] | grow.res$ScaffCol %in% c('FDR-sig', 'Bonf-sig')),] 
grow.res <- mutate(grow.res, ScaffCol = ifelse(test = ScaffCol %in% c('grey0', 'grey44') & chr %in% unique(grow.res$chr)[-c(1:24)], 
                                             yes = as.character('grey68'), no = as.character(ScaffCol))) 
grow.res$ScaffCol <- as.factor(grow.res$ScaffCol)

grow.annots <- 
  data.frame('Gene' = c('CSAD', 'METTL21E'),
             'RelPosition' = c(677460174, 253184989),
             'p' = c(6.23, 6.14),
             'ScaffCol' = rep('grey0', 2))



growth.manhat <-
  ggplot(grow.res, aes(RelPosition, p, colour = ScaffCol, size = ScaffCol)) + 
  geom_hline(yintercept = bonf.p, color = '#d7191c', 
             lty = 2, alpha = 0.5, size = 1.25) +
  geom_point(alpha = 0.75) +
  scale_color_manual(values = c('#d7191c', '#2b83ba', 'grey0', 'grey44', 'grey68')) +
  scale_size_manual(values = c(4.5,3,1.5,1.5,1.5)) +
  labs(title = 'Growth') +
  xlab('Relative Position') + 
  ylab(y.axis) +
  theme_light(base_size = 24) +
  scale_x_continuous(breaks = NULL) +
  geom_label_repel(data = grow.annots, aes(label = Gene), 
                   box.padding = 1, size = 4, 
                   fill = 'white') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = 'none')

png('../../Traits/output/Gemma-LMM-Bonf-LinkageThinned-growth.png', width = 12, height = 4, units = 'in', res = 600)
growth.manhat
dev.off()

pdf('../../Traits/output/Gemma-LMM-Bonf-LinkageThinned-growth.pdf', width = 12, height = 4)
growth.manhat
dev.off()

library(cowplot)


comp.manhat <-
  ggplot(comp.res, aes(RelPosition, p, colour = ScaffCol, size = ScaffCol)) + 
  geom_hline(yintercept = bonf.p, color = '#d7191c', 
             lty = 2, alpha = 0.5, size = 1.25) +
  geom_point(alpha = 0.75) +
  scale_color_manual(values = c('#d7191c', '#2b83ba', 'grey0', 'grey44', 'grey68')) +
  scale_size_manual(values = c(4.5,3,1.5,1.5,1.5)) +
  labs(title = 'Composite Fitness') +
  xlab('Relative Position') + 
  ylab(y.axis) +
  theme_light(base_size = 24) +
  scale_x_continuous(breaks = NULL) +
  geom_label_repel(data = comp.annots, aes(label = Gene), 
                   box.padding = 1.5, size = 6, 
                   fill = 'white') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = 'none')

growth.manhat <-
  ggplot(grow.res, aes(RelPosition, p, colour = ScaffCol, size = ScaffCol)) + 
  geom_hline(yintercept = bonf.p, color = '#d7191c', 
             lty = 2, alpha = 0.5, size = 1.25) +
  geom_point(alpha = 0.75) +
  scale_color_manual(values = c('#d7191c', '#2b83ba', 'grey0', 'grey44', 'grey68')) +
  scale_size_manual(values = c(4.5,3,1.5,1.5,1.5)) +
  labs(title = 'Growth') +
  xlab('Relative Position') + 
  ylab(y.axis) +
  theme_light(base_size = 24) +
  scale_x_continuous(breaks = NULL) +
  geom_label_repel(data = grow.annots, aes(label = Gene), 
                   box.padding = 2, size = 6, 
                   fill = 'white') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = 'none')

all.manhats <- 
  plot_grid(surv.manhat, 
            growth.manhat,
            comp.manhat,
            ncol = 1, align = 'V')

ggsave(all.manhats, filename = '../../Traits/output/Gemma-LMM-Bonf-LinkageThinned-SurvCompGroth.pdf', 
       width = 22.75, height = 16, useDingbats = F)
ggsave(all.manhats, filename = '../../Traits/output/Gemma-LMM-Bonf-LinkageThinned-SurvCompGroth.png', 
       width = 22.75, height = 16, dpi = 600)
