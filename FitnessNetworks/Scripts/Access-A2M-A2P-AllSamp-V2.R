# This script will determine and compare whether genotypic 
# paths between Generalists -> Molluscivore are more or less 
# accessible than those between Generalists -> Scale-eaters

# THIS IS A RE-ANALYSES, USING ALL SAMPLES

# We will do so by:
# 1) Constructing networks from SGV + Intro (for molluscivores)
# 2) Constructing networks from SGV + Intro + De Novo 
# (for scale eaters)

# Then, we will calculate:
# 1) the scaled number of accessible paths and 
# 2) the length of those paths

# Lastly we will compare these two metric for 
# 1) Generalist -> Molluscivore and 
# 2) Generalists to Scale-Eaters


library(ggplot2)
library(tidyverse)
library(igraph)
library(maditr)
library(cowplot)
library(ggforce)
library(pracma)

setwd('~/Dropbox/Research/Martin-Berkeley-Postdoc/Pupfish-Fitness-Genomics/Results/FixedSweeping_AllRedo/')
source('../../Scripts/NetworkSummaryFunctions-V2.R')

# Read in metadata so we can associate genotypes with lake/species/hybrid
meta <- read.table('../FixedSweeping/Hyb-Parent-Metadata.tsv', header = T)

dat <- read.table('../Hybrid-GamDat.tsv', header = T)
dat <- merge(meta, dat[,-c(2,5,7:8)], by = 'ID', all = T)
base.dir <- '~/Dropbox/Research/Martin-Berkeley-Postdoc/Pupfish-Fitness-Genomics/Results/FixedSweeping_AllRedo/'
source.dir <- '~/Dropbox/Research/Martin-Berkeley-Postdoc/Pupfish-Fitness-Genomics/Results/FixedSweeping/'

#######################################################################
# Okay, begin with scale eaters

genos1 <- read.table(paste0(source.dir, 'sgv-sites-P.012-mat.txt'), header = T)
genos2 <- read.table(paste0(source.dir, 'intro-sites-P.012-mat.txt'), header = T)
genos3 <- read.table(paste0(source.dir, 'denovo-sites-P.012-mat.txt'), header = T)

genos <- cbind(genos1, genos2[,-1], genos3[,-1])

# Now, drop out RTM1
genos <- genos[-which(genos$ID %in% rtm1),]

sgv.sites <- which(colnames(genos[-1]) %in% colnames(genos1[-1]))
intro.sites <- which(colnames(genos[-1]) %in% colnames(genos2[-1]))
dnv.sites <- which(colnames(genos[-1]) %in% colnames(genos3[-1]))

write.table(genos, 'sgv-denovo-intro-sites-P.012-mat-AllRedo.txt', sep = '\t')
genos <- merge(dat, genos, by = 'ID', all = T)
# Now, we are going to generate equal numbers of random 
# permutations that have some number of snps sourced from sgv, dnv, or intro
# and so on. Then, randomly sample 1000 from these without replacement
perms <- matrix(nrow = 12000, ncol = 5)
for(i in 1:2000){
  perms[i,] <-
    c(sample(sgv.sites, size = 1),
      sample(intro.sites, size = 1),
      sample(dnv.sites, size = 3))
  perms[i,] <- perms[i,][order(perms[i,])]
} 
for(i in 2001:4000){
  perms[i,] <-
    c(sample(sgv.sites, size = 1),
      sample(intro.sites, size = 2),
      sample(dnv.sites, size = 2))
  perms[i,] <- perms[i,][order(perms[i,])]
}
for(i in 4001:6000){
  perms[i,] <-
    c(sample(sgv.sites, size = 1),
      sample(intro.sites, size = 3),
      sample(dnv.sites, size = 1))
  perms[i,] <- perms[i,][order(perms[i,])]
}
for(i in 6001:8000){
  perms[i,] <-
    c(sample(sgv.sites, size = 2),
      sample(intro.sites, size = 1),
      sample(dnv.sites, size = 2))
  perms[i,] <- perms[i,][order(perms[i,])]
}
for(i in 8001:10000){
  perms[i,] <-
    c(sample(sgv.sites, size = 2),
      sample(intro.sites, size = 2),
      sample(dnv.sites, size = 1))
  perms[i,] <- perms[i,][order(perms[i,])]
}
for(i in 10001:12000){
  perms[i,] <-
    c(sample(sgv.sites, size = 3),
      sample(intro.sites, size = 1),
      sample(dnv.sites, size = 1))
  perms[i,] <- perms[i,][order(perms[i,])]
}

perms <- unique(perms)
perms <- perms[sample(1:12000, size = 5000, replace = F),]
write.table(perms, paste0(base.dir, 'sgv-denovo-intro-permuations-P-AllRedo.txt'), row.names = F, quote = F, col.names = F)

# Now, construct networks and summarize
all.comp.p <- 
  summarize.nets(dat = dat, genet.source = 'sgv-denovo-intro', gam.res = NULL, fitness = 'GrowSurv', spp = 'P',
                 n = 5000, base = base.dir, dir = paste0(base.dir, 'NetworkSummaries-raw/GrowthSurvival/'), raw = 'growth',
                 return.simple = TRUE)
all.surv.p <- 
  summarize.nets(dat = dat, genet.source = 'sgv-denovo-intro', gam.res = NULL, fitness = 'survival', spp = 'P',
                 n = 5000, base = base.dir, dir = paste0(base.dir, 'NetworkSummaries-raw/Survival/'), raw = 'survival',
                 return.simple = TRUE)
all.grow.p <- 
  summarize.nets(dat = dat, genet.source = 'sgv-denovo-intro', gam.res = NULL, fitness = 'growth', spp = 'P',
                 n = 5000, base = base.dir, dir = paste0(base.dir, 'NetworkSummaries-raw/Growth/'), raw = 'growth',
                 return.simple = TRUE)

#################################################################
# Do the same for molluscivores, using sgv and introgressed loci
genos1 <- read.table(paste0(source.dir, 'sgv-sites-M.012-mat.txt'), header = T)
genos2 <- read.table(paste0(source.dir, 'intro-sites-M.012-mat.txt'), header = T)

genos <- cbind(genos1, genos2[,-1])

# Now, drop out RTM1
genos <- genos[-which(genos$ID %in% rtm1),]

sgv.sites <- which(colnames(genos[-1]) %in% colnames(genos1[-1]))
intro.sites <- which(colnames(genos[-1]) %in% colnames(genos2[-1]))

write.table(genos, 'sgv-intro-sites-M.012-mat-AllRedo.txt', sep = '\t')
genos <- merge(dat, genos, by = 'ID', all = T)
# Now, we are going to generate equal numbers of random 
# permutations that have some number of snps sourced from sgv, dnv, or intro
# and so on. Then, randomly sample 1000 from these without replacement
perms <- matrix(nrow = 12000, ncol = 5)
for(i in 1:3000){
  perms[i,] <-
    c(sample(sgv.sites, size = 1),
      sample(intro.sites, size = 4))
  perms[i,] <- perms[i,][order(perms[i,])]
}
for(i in 3001:6000){
  perms[i,] <-
    c(sample(sgv.sites, size = 2),
      sample(intro.sites, size = 3))
  perms[i,] <- perms[i,][order(perms[i,])]
}
for(i in 6001:9000){
  perms[i,] <-
    c(sample(sgv.sites, size = 3),
      sample(intro.sites, size = 2))
  perms[i,] <- perms[i,][order(perms[i,])]
}
for(i in 9001:12000){
  perms[i,] <-
    c(sample(sgv.sites, size = 4),
      sample(intro.sites, size = 1))
  perms[i,] <- perms[i,][order(perms[i,])]
}

perms <- unique(perms)
perms <- perms[sample(1:12000, size = 5000, replace = F),]
write.table(perms, paste0(base.dir, 'sgv-intro-permuations-M-AllRedo.txt'), row.names = F, quote = F, col.names = F)

# Now, construct networks and summarize
all.comp.m <- 
  summarize.nets(dat = dat, genet.source = 'sgv-intro', gam.res = NULL, fitness = 'GrowSurv', spp = 'M',
                 n = 5000, base = base.dir, dir = paste0(base.dir, 'NetworkSummaries-raw/GrowthSurvival/'), raw = 'growth',
                 return.simple = TRUE)
all.surv.m <- 
  summarize.nets(dat = dat, genet.source = 'sgv-intro', gam.res = NULL, fitness = 'survival', spp = 'M',
                 n = 5000, base = base.dir, dir = paste0(base.dir, 'NetworkSummaries-raw/Survival/'), raw = 'survival',
                 return.simple = TRUE)
all.grow.m <- 
  summarize.nets(dat = dat, genet.source = 'sgv-intro', gam.res = NULL, fitness = 'growth', spp = 'M',
                 n = 5000, base = base.dir, dir = paste0(base.dir, 'NetworkSummaries-raw/Growth/'), raw = 'growth',
                 return.simple = TRUE)


#######################################################################
# Okay, now we can ask whether genotypic paths separating generalists #
# from each specialist differs in their accessibility                 #
#######################################################################
# Assign colors as the specialist to which paths are being sampled
setwd(base.dir)
m.col <- '#709AE1FF'
p.col <- '#F05C3BFF'

base_breaks <- function(n = 5){
  function(x) {
    axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
  }
}

plotdir <- '~/Dropbox/Research/Martin-Berkeley-Postdoc/Pupfish-Fitness-Genomics/Results/FixedSweeping_AllRedo/NetworkSummaries-raw/'

for(fitness in c('Survival', 'Growth', 'GrowthSurvival')){
  setwd(paste0('~/Dropbox/Research/Martin-Berkeley-Postdoc/Pupfish-Fitness-Genomics/Results/FixedSweeping_AllRedo/NetworkSummaries-raw/', fitness))
  if(fitness == 'Survival'){
    fit <- 'survival'
  }
  if(fitness == 'Growth'){
    fit <- 'growth'
  }
  if(fitness == 'GrowthSurvival'){
    fit <- 'GrowSurv'
  }
  
  Results <- 
    data.frame('Trajectory' = c('AtoM', 'AtoP'),
               'Source' = c('SGV + Introgression', 'SGV + Introgression + De Novo'), 
               'NumPaths.Med' = NA, 'NumPaths.Mean' = NA, 'NumPaths.SE' = NA,
               'RelNumPaths.Med' = NA, 'RelNumPaths.Mean' = NA, 'RelNumPaths.SE' = NA,
               'MinDistance.Med' = NA, 'MinDistance.Mean' = NA, 'MinDistance.SE' = NA,
               'MeanDistance.Med' = NA, 'MeanDistance.Mean' = NA, 'MeanDistance.SE' = NA,
               'MedDistance.Med' = NA, 'MedDistance.Mean' = NA, 'MedDistance.SE' = NA,
               'MeanFit.Med' = NA, 'MeanFit.Mean' = NA, 'MeanFit.SE' = NA,
               'MeanFit.ByLen.Med' = NA, 'MeanFit.ByLen.Mean' = NA, 'MeanFit.ByLen.SE' = NA,
               'MedFit.Med' = NA, 'MedFit.Mean' = NA, 'MedFit.SE' = NA,
               'MedFit.ByLen.Med' = NA, 'MedFit.ByLen.Mean' = NA, 'MedFit.ByLen.SE' = NA,
               'NumNodes.Med' = NA, 'NumNodes.Mean' = NA, 'NumNodes.SE' = NA,
               'NumEdges.Med' = NA, 'NumEdges.Mean' = NA, 'NumEdges.SE' = NA)
  
  m.all.summs <- na.omit(read.table(paste0('./sgv-intro/AtoM-MedMeanFit-PerNet-sgv-intro-M-', fit, '.tsv'), header = T))
  m.all.summs <- get.rel(m.all.summs)
  
  p.all.summs <- na.omit(read.table(paste0('./sgv-denovo-intro/AtoP-MedMeanFit-PerNet-sgv-denovo-intro-P-', fit, '.tsv'), header = T))
  p.all.summs <- get.rel(p.all.summs)
  
  
  m.all.summs$Trajectory <- 'AtoM'
  p.all.summs$Trajectory <- 'AtoP'
  m.all.summs$Source <- 'SGV+Intro'
  p.all.summs$Source <- 'SGV+Intro+DNV'
  
  meds <- 
    rbind(m.all.summs, 
          p.all.summs)
  
  meds$Source <- factor(meds$Source, levels = c('SGV+Intro','SGV+Intro+DNV'))
  
  # Begin by getting means and SE for the permutations, for each measure of fitness:
  m.all <- c()
  p.all <- c()
  
  for(i in 2:12){
    m.all <- c(m.all, summary_func(m.all.summs, i)[c(1:2,4)])
    p.all <- c(p.all, summary_func(p.all.summs, i)[c(1:2,4)])
    
    names(m.all) <- NULL
    names(p.all) <- NULL
  }
  
  Results[1:2,-c(1:2)] <- rbind(m.all, p.all)
  
  write.table(Results, file = paste0('../', fitness, 'AllSites-Accessibility-MedMeanSE-AllRedo.txt'), sep = '\t', quote = F)
  
  
  # Now, look at the relative number of accessible paths
  plts <- prep.plot(meds, 'RelNumPaths', mode = 'param', multiplier = 1.5)
  
  meds$X <- 
    gsub(meds$Trajectory, pattern = "AtoM", 
         replacement = "Generalist to Molluscivore") %>%
    gsub(., pattern = "AtoP", replacement = "Generalist to Scale-eater")
  
  RelNumPaths.p <- 
    ggplot(meds, aes(x = X, y = RelNumPaths, fill = Trajectory)) + 
    scale_fill_manual(values = c(m.col, p.col)) +
    geom_violin(alpha = 0.5) +
    geom_sina(alpha = 0.2, aes(color = Trajectory), size = 1.5, stroke = 0.0, bins = 50) +
    scale_color_manual(values = c(m.col, p.col)) +
    stat_boxplot(geom ='errorbar', width = 0.3, size = 1) +
    geom_boxplot(width = 0.075, outlier.colour = NA, notch = T) +
    # geom_point(position = position_jitterdodge(jitter.height = 0.1, jitter.width = 3), alpha = 0.05) + 
    # geom_text(x=0.5, y=plts$coords[2], label = plts$p,
    #           hjust = 0, size = 5) +
    scale_y_continuous(trans = log10_trans(), breaks = base_breaks(), limits = plts$coords) + 
    annotation_logticks(sides = 'l') +
    ylab(expression(frac("# Accessible Paths", "# Nodes in Network"))) + 
    theme_classic(base_size = 18) + 
    theme(legend.position = 'none', 
          axis.title.x = element_blank()) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 13))
  RelNumPaths.p
  
  # Calculate odds ratios for these
  RelNumPaths.or <- pairwiseOR(dat = meds, variable = "RelNumPaths", by.source = F)

  # Also include results for raw number of paths, and number of nodes for completeness
  NumPaths.or <- pairwiseOR(dat = meds, variable = "NumPaths", by.source = F)
  NumNodes.or <- pairwiseOR(dat = meds, variable = "NumNodes.InNet", by.source = F)
  
  # And plot the relationship between number of accessible paths, and either number of 
  # nodes or edges in networks (to justify scaling the number of paths).
  
  ### First, how strong is the relationship between the number of nodes in a network and the 
  ### number of edges?
  edge.node.mod1 <- 
    glm(NumEdges.InNet ~ NumNodes.InNet, 
        data = meds, family = "poisson")
  edge.node.mod2 <- 
    glm(NumEdges.InNet ~ NumNodes.InNet + Trajectory, 
        data = meds, family = "poisson")
  edge.node.mod3 <- 
    glm(NumEdges.InNet ~ NumNodes.InNet * Trajectory, 
        data = meds, family = "poisson")
  
  aic <- c(AIC(edge.node.mod1), AIC(edge.node.mod2), AIC(edge.node.mod3))
  best.mod <- list(edge.node.mod1, edge.node.mod2, edge.node.mod3)[[which(aic == min(aic))]]
  # Come up with a neat formula to include in the plot
  if(which(aic == min(aic)) == 1){
    edge.node.form <- "glm(# Edges ~ # Nodes, family = 'poisson')"
  }else{
    if(which(aic == min(aic)) == 2){
      edge.node.form <- "glm(# Edges ~ # Nodes + Trajectory, family = 'poisson')"
    }else{
      edge.node.form <- "glm(# Edges ~ # Nodes * Trajectory, family = 'poisson')"
    }
  }
  
  best.mod.summ <- summary(best.mod)
  p.edge.node <- best.mod.summ$coefficients[2,4]
  edge.node.coef <- paste0('exp(Coef) = ', round(exp(best.mod.summ$coefficients[2,1]), 4))
  
  if(p.edge.node < 0.0001){p.edge.node <- "P < 0.0001"}else{p.edge.node <- paste0("P = ", round(p.edge.node, 4))}
  
  # Get predicted values
  pred <- best.mod$fitted.values
  meds$pred <- pred
  
  edges.nodes.p <- 
    ggplot(data = meds, aes(y = NumEdges.InNet, x = NumNodes.InNet)) +
    geom_point(alpha = 0.5, size = 3.5, pch  = 21, aes(fill = Trajectory),
               position = position_jitterdodge(jitter.height = 0.1, 
                                               jitter.width = 0.1, 
                                               dodge.width = 0)) + 
    stat_smooth(aes(y = pred), method = 'glm', color ='black', size =2) + 
    annotate(geom = 'text', label = edge.node.form, size = 5,
             y = max(meds$NumEdges.InNet), x = min(meds$NumNodes.InNet),
             hjust = 0, vjust = -1, color = 'black') +
    annotate(geom = 'text', label = paste0(p.edge.node, '; ', edge.node.coef), size = 5,
             y = max(meds$NumEdges.InNet), x = min(meds$NumNodes.InNet),
             hjust = 0, vjust = 1, color = 'black') +
    scale_fill_manual(values = c(m.col, p.col)) +
    scale_y_continuous(trans = log10_trans(), breaks = base_breaks()) + 
    scale_x_continuous(trans = log10_trans(), breaks = base_breaks()) + 
    annotation_logticks(sides = 'lb') +
    ylab('# of Edges in Network') + 
    xlab('# of Nodes in Network') +
    theme_classic(base_size = 18) +
    theme(legend.position = 'none')
  edges.nodes.p
  
  # And a little table of the models/their fits. 
  edge.node.mods <- 
    data.frame(Model = c(format(edge.node.mod1$formula), 
                         format(edge.node.mod2$formula), 
                         format(edge.node.mod3$formula)),
               Family = rep('Poisson', 3),
               AIC = c(edge.node.mod1$aic, edge.node.mod2$aic, edge.node.mod3$aic),
               Coefficient = c(edge.node.mod1[[1]][2][[1]],
                               edge.node.mod2[[1]][2][[1]],
                               edge.node.mod3[[1]][2][[1]]),
               expCoefficient = exp(c(edge.node.mod1[[1]][2][[1]],
                                      edge.node.mod2[[1]][2][[1]],
                                      edge.node.mod3[[1]][2][[1]])),
               `P-value` = c(summary(edge.node.mod1)[[12]][2,4],
                             summary(edge.node.mod2)[[12]][2,4],
                             summary(edge.node.mod3)[[12]][2,4]))
  
  ######
  # Now, how do the number of accessible paths scale with the number of nodes in the network?
  path.node.mod1 <- 
    glm(NumPaths ~ NumNodes.InNet, 
        data = meds, family = "poisson")
  path.node.mod2 <- 
    glm(NumPaths ~ NumNodes.InNet + Trajectory, 
        data = meds, family = "poisson")
  path.node.mod3 <- 
    glm(NumPaths ~ NumNodes.InNet * Trajectory, 
        data = meds, family = "poisson")
  
  aic <- c(AIC(path.node.mod1), AIC(path.node.mod2), AIC(path.node.mod3))
  best.mod <- list(path.node.mod1, path.node.mod2, path.node.mod3)[[which(aic == min(aic))]]
  # Come up with a neat formula to include in the plot
  if(which(aic == min(aic)) == 1){
    path.node.form <- "glm(# Paths ~ # Nodes, family = 'poisson')"
  }else{
    if(which(aic == min(aic)) == 2){
      path.node.form <- "glm(# Paths ~ # Nodes + Trajectory, family = 'poisson')"
    }else{
      path.node.form <- "glm(# Paths ~ # Nodes * Trajectory, family = 'poisson')"
    }
  }
  
  best.mod.summ <- summary(best.mod)
  p.path.node <- best.mod.summ$coefficients[2,4]
  path.node.coef <- paste0('exp(Coef) = ', round(exp(best.mod.summ$coefficients[2,1]), 4))
  if(p.path.node < 0.0001){p.path.node <- "P < 0.0001"}else{p.path.node <- paste0("P = ", round(p.path.node, 4))}
  
  # Get predicted values
  pred <- best.mod$fitted.values
  meds$pred <- pred
  
  path.nodes.p <- 
    ggplot(data = meds, aes(y = NumPaths, x = NumNodes.InNet)) +
    geom_point(alpha = 0.5, size = 3.5, pch  = 21, aes(fill = Trajectory),
               position = position_jitterdodge(jitter.height = 0.05, 
                                               jitter.width = 0.05, 
                                               dodge.width = 0)) + 
    stat_smooth(aes(y = pred), method = 'glm', color ='black', size =2) + 
    annotate(geom = 'text', label = path.node.form, size = 5,
             y = max(meds$NumPaths), x = min(meds$NumNodes.InNet),
             hjust = 0, vjust = 0, color = 'black') +
    annotate(geom = 'text', label = paste0(p.path.node, '; ', path.node.coef), size = 5,
             y = max(meds$NumPaths), x = min(meds$NumNodes.InNet),
             hjust = 0, vjust = 2, color = 'black') +
    scale_fill_manual(values = c(m.col, p.col)) +
    scale_y_continuous(trans = log10_trans(), breaks = base_breaks()) + 
    scale_x_continuous(trans = log10_trans(), breaks = base_breaks()) + 
    annotation_logticks(sides = 'lb') +
    ylab('# of Accessible Paths') + 
    xlab('# of Nodes in Network') +
    theme_classic(base_size = 18) +
    theme(legend.position = 'none')
  path.nodes.p
  
  # And a little table of the models/their fits. 
  path.node.mods <- 
    data.frame(Model = c(format(path.node.mod1$formula), 
                         format(path.node.mod2$formula), 
                         format(path.node.mod3$formula)),
               Family = rep('Poisson', 3),
               AIC = c(path.node.mod1$aic, path.node.mod2$aic, path.node.mod3$aic),
               Coefficient = c(path.node.mod1[[1]][2][[1]],
                               path.node.mod2[[1]][2][[1]],
                               path.node.mod3[[1]][2][[1]]),
               expCoefficient = exp(c(path.node.mod1[[1]][2][[1]],
                                      path.node.mod2[[1]][2][[1]],
                                      path.node.mod3[[1]][2][[1]])),
               `P-value` = c(summary(path.node.mod1)[[12]][2,4],
                             summary(path.node.mod2)[[12]][2,4],
                             summary(path.node.mod3)[[12]][2,4]))

  ######
  # Lastly, how do the number of accessible paths scale with the number of edges in the network?
  path.edge.mod1 <- 
    glm(NumPaths ~ NumEdges.InNet, 
        data = meds, family = "poisson")
  path.edge.mod2 <- 
    glm(NumPaths ~ NumEdges.InNet + Trajectory, 
        data = meds, family = "poisson")
  path.edge.mod3 <- 
    glm(NumPaths ~ NumEdges.InNet * Trajectory, 
        data = meds, family = "poisson")
  
  aic <- c(AIC(path.edge.mod1), AIC(path.edge.mod2), AIC(path.edge.mod3))
  best.mod <- list(path.edge.mod1, path.edge.mod2, path.edge.mod3)[[which(aic == min(aic))]]
  # Come up with a neat formula to include in the plot
  if(which(aic == min(aic)) == 1){
    path.edge.form <- "glm(# Paths ~ # Edges, family = 'poisson')"
  }else{
    if(which(aic == min(aic)) == 2){
      path.edge.form <- "glm(# Paths ~ # Edges + Trajectory, family = 'poisson')"
    }else{
      path.edge.form <- "glm(# Paths ~ # Edges * Trajectory, family = 'poisson')"
    }
  }
  
  best.mod.summ <- summary(best.mod)
  p.path.edge <- best.mod.summ$coefficients[2,4]
  path.edge.coef <- paste0('exp(Coef) = ', round(exp(best.mod.summ$coefficients[2,1]), 4))
  if(p.path.edge < 0.0001){p.path.edge <- "P < 0.0001"}else{p.path.edge <- paste0("P = ", round(p.path.edge, 4))}
  
  # Get predicted values
  pred <- best.mod$fitted.values
  meds$pred <- pred
  
  path.edges.p <- 
    ggplot(data = meds, aes(y = NumPaths, x = NumEdges.InNet)) +
    geom_point(alpha = 0.5, size = 3.5, pch  = 21, aes(fill = Trajectory),
               position = position_jitterdodge(jitter.height = 0.05, 
                                               jitter.width = 0.05, 
                                               dodge.width = 0)) + 
    stat_smooth(aes(y = pred), method = 'glm', color ='black', size =2) + 
    annotate(geom = 'text', label = path.edge.form, size = 5,
             y = max(meds$NumPaths), x = min(meds$NumEdges.InNet),
             hjust = 0, vjust = 0, color = 'black') +
    annotate(geom = 'text', label = paste0(p.path.edge, '; ', path.edge.coef), size = 5,
             y = max(meds$NumPaths), x = min(meds$NumEdges.InNet),
             hjust = 0, vjust = 2, color = 'black') +
    scale_fill_manual(values = c(m.col, p.col)) +
    scale_y_continuous(trans = log10_trans(), breaks = base_breaks()) + 
    scale_x_continuous(trans = log10_trans(), breaks = base_breaks()) + 
    annotation_logticks(sides = 'lb') +
    ylab('# of Accessible Paths') + 
    xlab('# of Edges in Network') +
    theme_classic(base_size = 18) +
    theme(legend.position = 'none')
  path.edges.p
  
  # And a little table of the models/their fits. 
  path.edge.mods <- 
    data.frame(Model = c(format(path.edge.mod1$formula), 
                         format(path.edge.mod2$formula), 
                         format(path.edge.mod3$formula)),
               Family = rep('Poisson', 3),
               AIC = c(path.edge.mod1$aic, path.edge.mod2$aic, path.edge.mod3$aic),
               Coefficient = c(path.edge.mod1[[1]][2][[1]],
                               path.edge.mod2[[1]][2][[1]],
                               path.edge.mod3[[1]][2][[1]]),
               expCoefficient = exp(c(path.edge.mod1[[1]][2][[1]],
                                      path.edge.mod2[[1]][2][[1]],
                                      path.edge.mod3[[1]][2][[1]])),
               `P-value` = c(summary(path.edge.mod1)[[12]][2,4],
                             summary(path.edge.mod2)[[12]][2,4],
                             summary(path.edge.mod3)[[12]][2,4]))
  
  # Combine them all
  net.path.scaling.mods <-
    rbind(path.edge.mods,
          edge.node.mods,
          path.node.mods)
  
  # And write out
  write.table(net.path.scaling.mods, col.names = T, row.names = F, quote = F, sep = "\t",
              file = paste(fitness, '-NetSize-NumPathScaling-models-AllRedo.txt'))
  
  leg <- 
    get_legend(ggplot(data = meds, aes(y = NumPaths, x = NumNodes.InNet, 
                                       color = Trajectory)) +
                 geom_point(size = 3.5) + 
                 scale_color_manual(values = c(m.col, p.col), 
                                    labels = c("Generalist to Molluscivore",
                                               "Generalist to Scale-eater")) +
                 theme_classic(base_size = 18))
  
  # And plot together 
  path.net.size.p <- 
    plot_grid(edges.nodes.p, path.nodes.p, path.edges.p, leg, nrow = 2, ncol = 2,
              labels = c('a.', 'b.', 'c.', ''), label_size = 18)
  # path.net.size.p <- plot_grid(leg, path.net.size.p, ncol = 1)
  ggsave(path.net.size.p, filename = paste0(fitness, '-NetSize-NumAccessiblePaths-Scaling-AllRedo.pdf'), 
         height = 12, width = 12)
  
  
  #################################################################
  # And now the length of these paths
  mean(meds$MinDistance[which(meds$Trajectory == 'AtoM')])
  mean(meds$MinDistance[which(meds$Trajectory == 'AtoP')])
  
  plts <- prep.plot(meds, 'MinDistance', mode = 'nonparam')
  
  mindist.errs <- Results[,c(1:2,10:11)]
  mindist.errs$X <- c("Generalist to Molluscivore","Generalist to Scale-eater")
  colnames(mindist.errs)[3] <- 'MinDistance'
  mindist.errs$Upper <- mindist.errs$MinDistance + 2*mindist.errs$MinDistance.SE
  mindist.errs$Lower <- mindist.errs$MinDistance - 2*mindist.errs$MinDistance.SE
  
  min.dist.p <- 
    ggplot(meds, aes(x = X, y = MinDistance, fill = Trajectory)) + 
    scale_fill_manual(values = c(m.col, p.col)) +
    geom_sina(alpha = 0.2, aes(color = Trajectory), size = 1.5) +
    scale_color_manual(values = c(m.col, p.col)) +
    # geom_text(x=0.5, y=max(meds$MinDistance), label = plts$p,
    #           hjust = 0, size = 5) +
    #ylim(coords) + 
    ylab('Length of Shortest \n Accessible Path') + 
    theme_classic(base_size = 18) + 
    theme(legend.position = 'none', 
          axis.title.x = element_blank()) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 13)) +
    geom_errorbar(data = mindist.errs, aes(x = X, ymax = Upper, ymin = Lower), 
                  size = 1, width = 0.3) +
    geom_point(data = mindist.errs, aes(x = X, y = MinDistance, fill = Trajectory), 
               size = 5, pch = 21)
  min.dist.p
  
  # Calculate odds ratios for these
  MinDist.or <- pairwiseOR(dat = meds, variable = "MinDistance", by.source = F)

  allsites.access <- 
    plot_grid(RelNumPaths.p, min.dist.p, align = 'H')
  
  ggsave(allsites.access, filename = 'AllSites-Accessibility-A2M-A2P-AllRedo.pdf', width = 8, height = 8, useDingbats = F)
  
  ##################################################################################################################
  # Let's also look at some other statistics
  
  # Like mean fitness effect of steps along paths
  plts <- prep.plot(meds, 'MeanFit', mode = 'param')
  
  MeanFit.p <- 
    ggplot(meds, aes(x = X, y = MeanFit, fill = Trajectory)) + 
    scale_fill_manual(values = c(m.col, p.col)) +
    geom_violin(alpha = 0.5) +
    geom_sina(alpha = 0.2, aes(color = Trajectory), size = 1.5, stroke = 0.0, bins = 50) +
    scale_color_manual(values = c(m.col, p.col)) +
    stat_boxplot(geom ='errorbar', width = 0.3, size = 1) +
    geom_boxplot(width = 0.075, outlier.colour = NA) +
    # geom_point(position = position_jitterdodge(jitter.height = 0.1, jitter.width = 3), alpha = 0.05) + 
    geom_text(x=0.5, y=plts$coords[2], label = plts$p,
              hjust = 0, size = 5) +
    ylim(plts$coords) + 
    ylab('Mean Mutational \n Fitness Effect') + 
    theme_classic(base_size = 18) + 
    theme(legend.position = 'none', 
          axis.title.x = element_blank()) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 13))
  MeanFit.p
  
  ######################################################################
  # How does the mean fitness effect of mutational paths scale along with length of mutational paths?
  fit.scaling <- na.omit(meds[,c('Trajectory', 'MinDistance', 'MeanFit')])
  fit.scaling$MinDistance <- factor(fit.scaling$MinDistance, levels = c(1:20))
  dodge <- position_dodge(width = 0.9)
  MeanFit.LenScale.p <- 
    ggplot(fit.scaling, aes(x = MinDistance, y = MeanFit, fill = Trajectory)) + 
    scale_fill_manual(values = c(m.col, p.col)) +
    scale_color_manual(values = c(m.col, p.col)) +
    geom_violin(alpha = 0.5, position = dodge) +
    geom_sina(alpha = 0.2, aes(color = Trajectory), 
              size = 1, stroke = 0.0, bins = 50, position = dodge) +
    stat_boxplot(geom ='errorbar', width = 0.3, size = 1, position = dodge) +
    geom_boxplot(width = 0.075, outlier.colour = NA, position = dodge) +
    ylim(plts$coords) + 
    ylab('Mean Mutational Fitness Effect') + 
    xlab('Length of \n Accessible Path') + 
    theme_classic(base_size = 18) + 
    theme(legend.position = 'none')
  #scale_x_discrete(labels = function(x) str_wrap(x, width = 7))
  MeanFit.LenScale.p
  
  #Plot these together
  MeanFit.Len <- 
    plot_grid(MeanFit.p, MeanFit.LenScale.p, ncol = 1, rel_heights = c(1, 1), rel_widths = c(0.8,1), 
              align = 'hv', labels = c('a', 'b'), label_size = 24)
  
  ggsave(MeanFit.Len, filename = 'AllSites-MeanFit-LengthScaling-A2M-A2P-AllRedo.pdf', width = 12, height = 12, useDingbats = F)
  
  main.res <- Results[,c(1,31:32,4:5,7:8,10:11)]
  #########################################################################
  # Now, what about some general attributes of the networks, such as the geometry of fitness peaks? 
  # How many peaks are there? How connected are these peaks to the rest of the network?
  # How close are parental species to the peaks?
  
  m.nodes <- read.table(paste0('./sgv-intro/AtoM-MedMeanFit-PerNet-sgv-intro-M-', fit, '.tsv'), header = T)[,10]
  p.nodes <- read.table(paste0('./sgv-denovo-intro/AtoP-MedMeanFit-PerNet-sgv-denovo-intro-P-', fit, '.tsv'), header = T)[,10]
  
  m.all.summs <- read.table(paste0('./sgv-intro/NetworkSummaries-sgv-intro-M-', fit, '.tsv'), header = T)
  p.all.summs <- read.table(paste0('./sgv-denovo-intro/NetworkSummaries-sgv-denovo-intro-P-', fit, '.tsv'), header = T)
  
  m.all.summs$Trajectory <- 'AtoM'
  p.all.summs$Trajectory <- 'AtoP'
  m.all.summs$Source <- 'SGV+Intro'
  p.all.summs$Source <- 'SGV+Intro+DNV'
  m.all.summs$NumNodes.InNet <- m.nodes
  p.all.summs$NumNodes.InNet <- p.nodes
  
  m.all.summs <- na.omit(m.all.summs)
  p.all.summs <- na.omit(p.all.summs)
  
  # summs <- 
  #   rbind(m.all.summs, 
  #         p.all.summs)
  
  # Lets pull out the variables we're most curious about - this includes:
  # 1) Proportion of realized paths
  # 2) Number of peaks
  # 3) Mean peak fitness
  # 4) Mean peak neighboring slope
  # 5) Mean peak degree (# of edges leading to peak)
  # 6) Mean peak degree probability (how 'unusually connected are these peaks'typical' are the connectivity of these peaks)
  # 7) Mean distance of specialist to peak 
  # 8) Min distance of specialist to peak
  
  Results <- 
    data.frame(Trajectory = c('AtoM', 'AtoP'),
               Source = c('SGV + Introgression', 'SGV + Introgression + De Novo'), 
               ProportionRealizedPaths = NA, ProportionRealizedPaths.SE = NA, 
               NumPeaks = NA, NumPeaks.SE = NA, 
               MeanPeakFitness = NA, MeanPeakFitness.SE = NA, 
               MeanPeakNeighborSlope = NA, MeanPeakNeighborSlope.SE = NA, 
               MeanNumPathsToPeaks = NA, MeanNumPathsToPeaks.SE = NA, 
               MeanPeakDegreeProb = NA,  MeanPeakDegreeProb.SE = NA, 
               MaxPeakFitness = NA, MaxPeakFitness.SE = NA, 
               MeanNumPathsToMaxPeak = NA, MeanNumPathsToMaxPeak.SE = NA, 
               MaxPeakDegreeProb = NA,  MaxPeakDegreeProb.SE = NA, 
               MeanDistanceToPeak = NA, MeanDistanceToPeak.SE = NA, 
               MinDistanceToPeak = NA, MinDistanceToPeak.SE = NA,
               MeanDistanceToMaxPeak = NA, MeanDistanceToMaxPeak.SE = NA, 
               MinDistanceToMaxPeak = NA, MinDistanceToMaxPeak.SE = NA,
               RelNumPeaks = NA, RelNumPeaks.SE = NA,
               RelNumPathsToPeaks = NA, RelNumPathsToPeaks.SE = NA)
  
  m.all.summs <- m.all.summs[,c(29:30,12:17,22:24,20,18,27,25,31)]
  p.all.summs <- p.all.summs[,c(29:30,12:17,22:24,21,19,28,26,31)]
  colnames(m.all.summs)[12:15] <- c('MeanDistanceToPeak', 'MinDistanceToPeak', 
                                    'MeanDistanceToMaxPeak', 'MinDistanceToMaxPeak')
  colnames(p.all.summs)[12:15] <- c('MeanDistanceToPeak', 'MinDistanceToPeak',
                                    'MeanDistanceToMaxPeak', 'MinDistanceToMaxPeak')
  m.all.summs$RelNumPeaks <- m.all.summs$NumPeaks / m.all.summs$NumNodes.InNet
  p.all.summs$RelNumPeaks <- p.all.summs$NumPeaks / p.all.summs$NumNodes.InNet
  m.all.summs$RelNumPathsToPeaks <- m.all.summs$MeanNumPathsToPeaks / m.all.summs$NumNodes.InNet
  p.all.summs$RelNumPathsToPeaks <- p.all.summs$MeanNumPathsToPeaks / p.all.summs$NumNodes.InNet
  # Begin by getting means and SE for the permutations, for each measure of fitness:
  m.all <- c()
  p.all <- c()
  
  for(i in c(3:15,17:18)){
    m.all <- c(m.all, summary_func(m.all.summs, i)[c(2,4)])
    p.all <- c(p.all, summary_func(p.all.summs, i)[c(2,4)])
    
    names(m.all) <- NULL
    names(p.all) <- NULL
  }
  
  Results[1:2,-c(1:2)] <- rbind(m.all, p.all)
  write.table(Results, file = paste0('../', fitness, 'AllSites-NetworkSummaries-MeanSE-AllRedo.txt'), sep = '\t', quote = F)
  
  summs <- rbind(m.all.summs, p.all.summs)
  
  summs$X <- 
    gsub(summs$Trajectory, pattern = "AtoM", 
         replacement = "Generalist to Molluscivore") %>%
    gsub(., pattern = "AtoP", replacement = "Generalist to Scale-eater")
  
  # First proportion of realized paths
  plts <- prep.plot(summs, 'ProportionRealizedPaths', mode = 'param')
  
  RealPaths.p <- 
    ggplot(summs, aes(x = X, y = ProportionRealizedPaths, fill = Trajectory)) + 
    scale_fill_manual(values = c(m.col, p.col)) +
    geom_violin(alpha = 0.5) +
    geom_sina(alpha = 0.2, aes(color = Trajectory), size = 1.5, stroke = 0.0, bins = 50) +
    scale_color_manual(values = c(m.col, p.col)) +
    stat_boxplot(geom ='errorbar', width = 0.3, size = 1) +
    geom_boxplot(width = 0.075, outlier.colour = NA, notch = T) +
    # geom_point(position = position_jitterdodge(jitter.height = 0.1, jitter.width = 3), alpha = 0.05) + 
    geom_text(x=0.5, y=plts$coords[2], label = plts$p,
              hjust = 0, size = 5) +
    ylim(plts$coords) + 
    ylab('Proportion of \n Realized Paths') + 
    theme_classic(base_size = 18) + 
    theme(legend.position = 'none', 
          axis.title.x = element_blank()) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 13))
  
  # Then number of peaks
  plts <- prep.plot(summs, 'NumPeaks', mode = 'nonparam')
  
  numpeaks.errs <- Results[,c(1:2,5:6)]
  numpeaks.errs$X <- c("Generalist to Molluscivore", "Generalist to Scale-eater")
  colnames(numpeaks.errs)[3] <- 'NumPeaks'
  numpeaks.errs$Upper <- numpeaks.errs$NumPeaks + 2*numpeaks.errs$NumPeaks.SE
  numpeaks.errs$Lower <- numpeaks.errs$NumPeaks - 2*numpeaks.errs$NumPeaks.SE
  
  
  NumPeaks.p <- 
    ggplot(summs, aes(x = X, y = NumPeaks, fill = Trajectory)) + 
    scale_fill_manual(values = c(m.col, p.col)) +
    geom_sina(alpha = 0.2, aes(color = Trajectory), size = 1.5) +
    scale_color_manual(values = c(m.col, p.col)) +
    # geom_text(x=0.5, y=max(summs$NumPeaks), label = plts$p,
    #           hjust = 0, size = 5) +
    #ylim(coords) + 
    ylab('Number of Peaks') + 
    theme_classic(base_size = 18) + 
    theme(legend.position = 'none',
          axis.title.x = element_blank()) +
    scale_y_continuous(breaks = seq(from = min(na.omit(summs$NumPeaks)), 
                                    to = max(na.omit(summs$NumPeaks)), by = 3)) + 
    scale_x_discrete(labels = function(x) str_wrap(x, width = 13)) +
    geom_errorbar(data = numpeaks.errs, aes(x = X, ymax = Upper, ymin = Lower), 
                  size = 1, width = 0.3) +
    geom_point(data = numpeaks.errs, aes(x = X, y = NumPeaks, fill = Trajectory), 
               size = 5, pch = 21)
  
  # Calculate odds ratios
  NumPeaks.or <- pairwiseOR(dat = summs, variable = "NumPeaks", by.source = F)
  
  # The same for the relative number of peaks
  # Then number of peaks
  plts <- prep.plot(summs, 'RelNumPeaks', mode = 'nonparam')
  
  relnumpeaks.errs <- Results[,c(1:2,29:30)]
  relnumpeaks.errs$X <- c("Generalist to Molluscivore", "Generalist to Scale-eater")
  colnames(relnumpeaks.errs)[3] <- 'RelNumPeaks'
  relnumpeaks.errs$Upper <- relnumpeaks.errs$RelNumPeaks + 2*relnumpeaks.errs$RelNumPeaks.SE
  relnumpeaks.errs$Lower <- relnumpeaks.errs$RelNumPeaks - 2*relnumpeaks.errs$RelNumPeaks.SE
  
  RelNumPeaks.p <- 
    ggplot(summs, aes(x = X, y = RelNumPeaks, fill = Trajectory)) + 
    scale_fill_manual(values = c(m.col, p.col)) +
    geom_violin(alpha = 0.5) +
    geom_sina(alpha = 0.2, aes(color = Trajectory), size = 1.5, stroke = 0.0, bins = 50) +
    scale_color_manual(values = c(m.col, p.col)) +
    stat_boxplot(geom ='errorbar', width = 0.3, size = 1) +
    geom_boxplot(width = 0.075, outlier.colour = NA, notch = T) +
    ylim(plts$coords) + 
    # geom_point(position = position_jitterdodge(jitter.height = 0.1, jitter.width = 3), alpha = 0.05) + 
    # geom_text(x=0.5, y=plts$coords[2], label = plts$p,
    #           hjust = 0, size = 5) +
    #scale_y_continuous(trans = log10_trans(), breaks = base_breaks(), limits = plts$coords) + 
    #annotation_logticks(sides = 'l') +
    ylab(expression(frac("# Peaks", "# Nodes in Network"))) + 
    theme_classic(base_size = 18) + 
    theme(legend.position = 'none', 
          axis.title.x = element_blank()) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 13))
  
  # Calculate odds ratios
  RelNumPeaks.or <- pairwiseOR(dat = summs, variable = "RelNumPeaks", by.source = F)

  # Scaled # of paths to peaks
  plts <- prep.plot(summs, 'RelNumPathsToPeaks', mode = 'nonparam')
  
  numpathstopeaks.errs <- Results[,c(1:2,31:32)]
  numpathstopeaks.errs$X <- c("Generalist to Molluscivore", "Generalist to Scale-eater")
  colnames(numpathstopeaks.errs)[3] <- 'RelNumPathsToPeaks'
  numpathstopeaks.errs$Upper <- numpathstopeaks.errs$RelNumPathsToPeaks + 2*numpathstopeaks.errs$RelNumPathsToPeaks.SE
  numpathstopeaks.errs$Lower <- numpathstopeaks.errs$RelNumPathsToPeaks - 2*numpathstopeaks.errs$RelNumPathsToPeaks.SE
  
  RelNumPathsToPeaks.p <- 
    ggplot(summs, aes(x = X, y = RelNumPathsToPeaks, fill = Trajectory)) + 
    scale_fill_manual(values = c(m.col, p.col)) +
    geom_violin(alpha = 0.5) +
    geom_sina(alpha = 0.2, aes(color = Trajectory), size = 1.5, stroke = 0.0, bins = 50) +
    scale_color_manual(values = c(m.col, p.col)) +
    stat_boxplot(geom ='errorbar', width = 0.3, size = 1) +
    geom_boxplot(width = 0.075, outlier.colour = NA, notch = T) +
    ylim(plts$coords) + 
    # geom_point(position = position_jitterdodge(jitter.height = 0.1, jitter.width = 3), alpha = 0.05) + 
    # geom_text(x=0.5, y=plts$coords[2], label = plts$p,
    #           hjust = 0, size = 5) +
    #scale_y_continuous(trans = log10_trans(), breaks = base_breaks(), limits = plts$coords) + 
    #annotation_logticks(sides = 'l') +
    ylab(expression(frac("# Paths to Peaks", "# Nodes in Network"))) + 
    theme_classic(base_size = 18) + 
    theme(legend.position = 'none', 
          axis.title.x = element_blank()) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 13))
  
  # Calculate odds ratios
  RelNumPathsToPeaks.or <- pairwiseOR(dat = summs, variable = "RelNumPathsToPeaks", by.source = F)

  # Mean peak fitness
  plts <- prep.plot(summs, 'MeanPeakFitness', mode = 'param')
  
  MeanPeakFit.p <- 
    ggplot(summs, aes(x = X, y = MeanPeakFitness, fill = Trajectory)) + 
    scale_fill_manual(values = c(m.col, p.col)) +
    geom_violin(alpha = 0.5) +
    geom_sina(alpha = 0.2, aes(color = Trajectory), size = 1.5, stroke = 0.0, bins = 50) +
    scale_color_manual(values = c(m.col, p.col)) +
    stat_boxplot(geom ='errorbar', width = 0.3, size = 1) +
    geom_boxplot(width = 0.075, outlier.colour = NA, notch = T) +
    # geom_text(x=0.5, y=plts$coords[2], label = plts$p,
    #           hjust = 0, size = 5) +
    ylim(plts$coords) + 
    ylab('Mean Peak Fitness') + 
    theme_classic(base_size = 18) + 
    theme(legend.position = 'none', 
          axis.title.x = element_blank()) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 13))
  
  # Calculate odds ratios
  MeanPeakFit.or <- pairwiseOR(dat = summs, variable = "MeanPeakFitness", by.source = F)

  # Max peak fitness
  plts <- prep.plot(summs, 'MaxPeakFitness', mode = 'param')
  
  MaxPeakFit.p <- 
    ggplot(summs, aes(x = X, y = MaxPeakFitness, fill = Trajectory)) + 
    scale_fill_manual(values = c(m.col, p.col)) +
    geom_violin(alpha = 0.5) +
    geom_sina(alpha = 0.2, aes(color = Trajectory), size = 1.5, stroke = 0.0, bins = 50) +
    scale_color_manual(values = c(m.col, p.col)) +
    stat_boxplot(geom ='errorbar', width = 0.3, size = 1) +
    geom_boxplot(width = 0.075, outlier.colour = NA, notch = T) +
    geom_text(x=0.5, y=plts$coords[2], label = plts$p,
              hjust = 0, size = 5) +
    ylim(plts$coords) + 
    ylab('Max Peak Fitness') + 
    theme_classic(base_size = 18) + 
    theme(legend.position = 'none', 
          axis.title.x = element_blank()) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 13))
  
  # Mean peak neighboring slope
  plts <- prep.plot(summs, 'MeanPeakNeighborSlope', mode = 'param')
  
  MeanPeakSlope.p <- 
    ggplot(summs, aes(x = X, y = MeanPeakNeighborSlope, fill = Trajectory)) + 
    scale_fill_manual(values = c(m.col, p.col)) +
    geom_violin(alpha = 0.5) +
    geom_sina(alpha = 0.2, aes(color = Trajectory), size = 1.5, stroke = 0.0, bins = 50) +
    scale_color_manual(values = c(m.col, p.col)) +
    stat_boxplot(geom ='errorbar', width = 0.3, size = 1) +
    geom_boxplot(width = 0.075, outlier.colour = NA, notch = T) +
    # geom_text(x=0.5, y=plts$coords[2], label = plts$p,
    #           hjust = 0, size = 5) +
    ylim(plts$coords) + 
    ylab('Mean Peak \n Neighboring Slope') + 
    theme_classic(base_size = 18) + 
    theme(legend.position = 'none', 
          axis.title.x = element_blank()) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 13))
  
  # Calculate odds ratios
  MeanPeakSlope.or <- pairwiseOR(dat = summs, variable = "MeanPeakNeighborSlope", by.source = F)

  # Mean peak degree
  plts <- prep.plot(summs, 'MeanNumPathsToPeaks', mode = 'param')
  
  PeakDegree.p <- 
    ggplot(summs, aes(x = X, y = MeanNumPathsToPeaks, fill = Trajectory)) + 
    scale_fill_manual(values = c(m.col, p.col)) +
    geom_violin(alpha = 0.5) +
    geom_sina(alpha = 0.2, aes(color = Trajectory), size = 1.5, stroke = 0.0, bins = 50) +
    scale_color_manual(values = c(m.col, p.col)) +
    stat_boxplot(geom ='errorbar', width = 0.3, size = 1) +
    geom_boxplot(width = 0.075, outlier.colour = NA, notch = T) +
    # geom_text(x=0.5, y=plts$coords[2], label = plts$p,
    #           hjust = 0, size = 5) +
    ylim(plts$coords) + 
    ylab('Mean # Paths \n to Peak') + 
    theme_classic(base_size = 18) + 
    theme(legend.position = 'none', 
          axis.title.x = element_blank()) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 13))
  
  # Calculate odds ratios
  PeakDegree.or <- pairwiseOR(dat = summs, variable = "MeanNumPathsToPeaks", by.source = F)

  # Max peak degree
  plts <- prep.plot(summs, 'MeanNumPathsToMaxPeak', mode = 'param')
  
  maxpeakdeg.errs <- Results[,c(1:2,17:18)]
  maxpeakdeg.errs$X <- c("Generalist to Molluscivore","Generalist to Scale-eater")
  colnames(maxpeakdeg.errs)[3] <- 'MeanNumPathsToMaxPeak'
  maxpeakdeg.errs$Upper <- maxpeakdeg.errs$MeanNumPathsToMaxPeak + 2*maxpeakdeg.errs$MeanNumPathsToMaxPeak.SE
  maxpeakdeg.errs$Lower <- maxpeakdeg.errs$MeanNumPathsToMaxPeak - 2*maxpeakdeg.errs$MeanNumPathsToMaxPeak.SE
  
  MaxPeakDegree.p <- 
    ggplot(summs, aes(x = X, y = MeanNumPathsToMaxPeak, fill = Trajectory)) + 
    scale_fill_manual(values = c(m.col, p.col)) +
    geom_sina(alpha = 0.2, aes(color = Trajectory), size = 1.5) +
    scale_color_manual(values = c(m.col, p.col)) +
    geom_text(x=0.5, y=max(summs$MeanNumPathsToMaxPeak), label = plts$p,
              hjust = 0, size = 5) +
    #ylim(coords) + 
    ylab('# Paths to Highest Peak') + 
    theme_classic(base_size = 18) + 
    theme(legend.position = 'none',
          axis.title.x = element_blank()) +
    scale_y_continuous(breaks = seq(from = min(na.omit(summs$MeanNumPathsToMaxPeak)), 
                                    to = max(na.omit(summs$MeanNumPathsToMaxPeak)), by = 3)) + 
    scale_x_discrete(labels = function(x) str_wrap(x, width = 13)) +
    geom_errorbar(data = maxpeakdeg.errs, aes(x = X, ymax = Upper, ymin = Lower), 
                  size = 1, width = 0.3) +
    geom_point(data = maxpeakdeg.errs, aes(x = X, y = MeanNumPathsToMaxPeak, fill = Trajectory), 
               size = 5, pch = 21)
  
  # Mean peak degree probability
  plts <- prep.plot(summs, 'MeanPeakDegreeProb', mode = 'param')
  
  PeakDegreeProb.p <- 
    ggplot(summs, aes(x = X, y = MeanPeakDegreeProb, fill = Trajectory)) + 
    scale_fill_manual(values = c(m.col, p.col)) +
    geom_violin(alpha = 0.5) +
    geom_sina(alpha = 0.2, aes(color = Trajectory), size = 1.5, stroke = 0.0, bins = 50) +
    scale_color_manual(values = c(m.col, p.col)) +
    stat_boxplot(geom ='errorbar', width = 0.3, size = 1) +
    geom_boxplot(width = 0.075, outlier.colour = NA, notch = T) +
    geom_text(x=0.5, y=plts$coords[2], label = plts$p,
              hjust = 0, size = 5) +
    ylim(plts$coords) + 
    ylab('Mean Peak Degree \n Probability') + 
    theme_classic(base_size = 18) + 
    theme(legend.position = 'none', 
          axis.title.x = element_blank()) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 13))
  
  # Mean peak degree for the highest peak
  plts <- prep.plot(summs, 'MaxPeakDegreeProb', mode = 'nonparam')
  
  MaxPeakDegreeProb.p <- 
    ggplot(summs, aes(x = X, y = MaxPeakDegreeProb, fill = Trajectory)) + 
    scale_fill_manual(values = c(m.col, p.col)) +
    geom_violin(alpha = 0.5) +
    geom_sina(alpha = 0.2, aes(color = Trajectory), size = 1.5, stroke = 0.0, bins = 50) +
    scale_color_manual(values = c(m.col, p.col)) +
    stat_boxplot(geom ='errorbar', width = 0.3, size = 1) +
    geom_boxplot(width = 0.075, outlier.colour = NA, notch = T) +
    geom_text(x=0.5, y=plts$coords[2], label = plts$p,
              hjust = 0, size = 5) +
    ylim(plts$coords) + 
    ylab('Mean Degree Probability \n for Highest Peak') + 
    theme_classic(base_size = 18) + 
    theme(legend.position = 'none', 
          axis.title.x = element_blank()) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 13))
  
  # Mean distance to peak
  plts <- prep.plot(summs, 'MeanDistanceToPeak', mode = 'param')
  
  MeanDistToPeak.p <- 
    ggplot(summs, aes(x = X, y = MeanDistanceToPeak, fill = Trajectory)) + 
    scale_fill_manual(values = c(m.col, p.col)) +
    geom_violin(alpha = 0.5) +
    geom_sina(alpha = 0.2, aes(color = Trajectory), size = 1.5, stroke = 0.0, bins = 50) +
    scale_color_manual(values = c(m.col, p.col)) +
    stat_boxplot(geom ='errorbar', width = 0.3, size = 1) +
    geom_boxplot(width = 0.075, outlier.colour = NA, notch = T) +
    geom_text(x=0.5, y=plts$coords[2], label = plts$p,
              hjust = 0, size = 5) +
    ylim(plts$coords) + 
    ylab('Mean Distance to Peak') + 
    theme_classic(base_size = 18) + 
    theme(legend.position = 'none', 
          axis.title.x = element_blank()) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 13))
  
  # Min distance to peak
  plts <- prep.plot(summs, 'MinDistanceToPeak', mode = 'nonparam')
  
  peakdist.errs <- Results[,c(1:2,23:24)]
  peakdist.errs$X <- c("Generalist to Molluscivore","Generalist to Scale-eater")
  colnames(peakdist.errs)[3] <- 'MinDistanceToPeak'
  peakdist.errs$Upper <- peakdist.errs$MinDistanceToPeak + 2*peakdist.errs$MinDistanceToPeak.SE
  peakdist.errs$Lower <- peakdist.errs$MinDistanceToPeak - 2*peakdist.errs$MinDistanceToPeak.SE
  
  MinDistToPeak.p <- 
    ggplot(summs, aes(x = X, y = MinDistanceToPeak, fill = Trajectory)) + 
    scale_fill_manual(values = c(m.col, p.col)) +
    geom_sina(alpha = 0.2, aes(color = Trajectory), size = 1.5) +
    scale_color_manual(values = c(m.col, p.col)) +
    # geom_text(x=0.5, y=max(summs$MinDistanceToPeak), label = plts$p,
    #           hjust = 0, size = 5) +
    #ylim(coords) + 
    ylab('Minimum Distance \n to Peak') + 
    theme_classic(base_size = 18) + 
    theme(legend.position = 'none',
          axis.title.x = element_blank()) +
    scale_y_continuous(breaks = seq(from = min(na.omit(summs$MinDistanceToPeak)), 
                                    to = max(na.omit(summs$MinDistanceToPeak)), by = 3)) + 
    scale_x_discrete(labels = function(x) str_wrap(x, width = 13)) +
    geom_errorbar(data = peakdist.errs, aes(x = X, ymax = Upper, ymin = Lower), 
                  size = 1, width = 0.3) +
    geom_point(data = peakdist.errs, aes(x = X, y = MinDistanceToPeak, fill = Trajectory), 
               size = 5, pch = 21)
  
  # Calculate odds ratios
  MinDistToPeak.or <- pairwiseOR(dat = summs, variable = "MinDistanceToPeak", by.source = F)

  # Mean distance to peak
  plts <- prep.plot(summs, 'MeanDistanceToMaxPeak', mode = 'param')
  
  MeanDistToMaxPeak.p <- 
    ggplot(summs, aes(x = X, y = MeanDistanceToMaxPeak, fill = Trajectory)) + 
    scale_fill_manual(values = c(m.col, p.col)) +
    geom_violin(alpha = 0.5) +
    geom_sina(alpha = 0.2, aes(color = Trajectory), size = 1.5, stroke = 0.0, bins = 50) +
    scale_color_manual(values = c(m.col, p.col)) +
    stat_boxplot(geom ='errorbar', width = 0.3, size = 1) +
    geom_boxplot(width = 0.075, outlier.colour = NA, notch = T) +
    geom_text(x=0.5, y=plts$coords[2], label = plts$p,
              hjust = 0, size = 5) +
    ylim(plts$coords) + 
    ylab('Mean Distance to \n Highest Peak') + 
    theme_classic(base_size = 18) + 
    theme(legend.position = 'none', 
          axis.title.x = element_blank()) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 13))
  
  # Min distance to peak
  plts <- prep.plot(summs, 'MinDistanceToMaxPeak', mode = 'nonparam')
  
  peakdist.errs <- Results[,c(1:2,27:28)]
  peakdist.errs$X <- c("Generalist to Molluscivore","Generalist to Scale-eater")
  colnames(peakdist.errs)[3] <- 'MinDistanceToMaxPeak'
  peakdist.errs$Upper <- peakdist.errs$MinDistanceToMaxPeak + 2*peakdist.errs$MinDistanceToMaxPeak.SE
  peakdist.errs$Lower <- peakdist.errs$MinDistanceToMaxPeak - 2*peakdist.errs$MinDistanceToMaxPeak.SE
  
  
  MinDistToMaxPeak.p <- 
    ggplot(summs, aes(x = X, y = MinDistanceToMaxPeak, fill = Trajectory)) + 
    scale_fill_manual(values = c(m.col, p.col)) +
    geom_sina(alpha = 0.2, aes(color = Trajectory), size = 1.5) +
    scale_color_manual(values = c(m.col, p.col)) +
    geom_text(x=0.5, y=max(summs$MinDistanceToMaxPeak), label = plts$p,
              hjust = 0, size = 5) +
    #ylim(coords) + 
    ylab('Min Distance to \n Highest Peak') + 
    theme_classic(base_size = 18) + 
    theme(legend.position = 'none',
          axis.title.x = element_blank()) +
    scale_y_continuous(breaks = seq(from = min(na.omit(summs$MinDistanceToMaxPeak)), 
                                    to = max(na.omit(summs$MinDistanceToMaxPeak)), by = 3)) + 
    scale_x_discrete(labels = function(x) str_wrap(x, width = 13)) +
    geom_errorbar(data = peakdist.errs, aes(x = X, ymax = Upper, ymin = Lower), 
                  size = 1, width = 0.3) +
    geom_point(data = peakdist.errs, aes(x = X, y = MinDistanceToMaxPeak, fill = Trajectory), 
               size = 5, pch = 21)
  
  
  # Now plot them all 
  mean.peaks.p <- plot_grid(MeanPeakFit.p, PeakDegree.p, 
                            MeanDistToPeak.p, MinDistToPeak.p, 
                            labels = c('a', 'b', 'c', 'd'),
                            label_size = 24)
  
  max.peaks.p <- plot_grid(MaxPeakFit.p, MaxPeakDegree.p, 
                           MeanDistToMaxPeak.p, MinDistToMaxPeak.p, 
                           labels = c('a', 'b', 'c', 'd'),
                           label_size = 24)
  if(fitness != "Survival"){
    max.peaks.p <- plot_grid(MaxPeakFit.p, MaxPeakDegree.p, 
                             MeanDistToMaxPeak.p, MinDistToMaxPeak.p, 
                             labels = c('a', 'b', 'c', 'd'),
                             label_size = 24)
  }
  
  # Combine the odds ratios that will be plotted together
  main.or <- rbind(RelNumPaths.or, MinDist.or, NumPeaks.or, RelNumPathsToPeaks.or, MinDistToPeak.or)
  supp.or <- rbind(MeanPeakFit.or, MeanPeakSlope.or, PeakDegree.or)
  
  # Rename accessibility measures for plotting
  main.or$Accessibility <- c('d.', 'e.', 'f.', 'g.', 'h.')
  supp.or$Accessibility <- c('Mean peak fitness', 'Mean peak slope',
                             'Mean # of paths to peaks')
  
  # Reorder factor levels for plotting
  main.or$Accessibility <- factor(main.or$Accessibility, 
                                  levels = c(main.or[1,2], main.or[2,2], 
                                             main.or[3,2], main.or[4,2], 
                                             main.or[5,2]))
  supp.or$Accessibility <- factor(supp.or$Accessibility, 
                                  levels = c(supp.or[1,2], supp.or[2,2], 
                                             supp.or[3,2]))
  
  # Plot them
  main.or.p <- 
    ggplot(main.or, aes(x = Accessibility)) + 
    geom_errorbar(aes(ymin = OR.Lower, ymax = OR.Upper), width = 0.25, size = 1) + 
    geom_point(aes(y = OddsRatio), fill = 'grey', pch = 21, size = 4) +
    geom_hline(yintercept = 1, 
               alpha = 0.5, size = 0.3, lty = 2) + 
    scale_y_continuous(trans = log10_trans(), breaks = base_breaks()) + 
    annotation_logticks(sides = 'l') + 
    ylab('Odds Ratio: \n AtoM/AtoP') + 
    theme_classic(base_size = 18) + 
    theme(legend.position = 'none', 
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank(),
          strip.text.x = element_blank()) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 16)) 
  main.or.p
  
  supp.or.p <- 
    ggplot(supp.or, aes(x = Accessibility)) + 
    geom_errorbar(aes(ymin = OR.Lower, ymax = OR.Upper), width = 0.25, size = 1) + 
    geom_point(aes(y = OddsRatio), fill = 'grey', pch = 21, size = 4) +
    geom_hline(yintercept = 1, 
               alpha = 0.5, size = 0.3, lty = 2) + 
    scale_y_continuous(trans = log10_trans(), breaks = base_breaks()) + 
    annotation_logticks(sides = 'l') + 
    ylab('Odds Ratio: \n AtoM/AtoP') + 
    theme_classic(base_size = 18) + 
    theme(legend.position = 'none', 
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank(),
          strip.text.x = element_blank()) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 16))
  supp.or.p
  
  # Plot results together in a nice format to be put into the manuscript
  main.res.upper <- 
    plot_grid(RelNumPaths.p, min.dist.p, 
              NumPeaks.p, align = 'vh', ncol = 3,
              rel_widths = c(1,1,1))
  main.res.lower <- 
    plot_grid(RelNumPathsToPeaks.p, MinDistToPeak.p, main.or.p, 
              align = 'vh', ncol = 3, rel_widths = c(1,1,1))
  main.res.p <- plot_grid(main.res.upper, main.res.lower, nrow = 2, align = 'vh')
  
  ggsave(main.res.p, filename = 'AllSites-A2M-A2P-AccessibilityPeakSummary.pdf', width = 12, height = 7, useDingbats = F)
  ggsave(main.res.p, filename = 'AllSites-A2M-A2P-AccessibilityPeakSummary.png', width = 12, height = 7)
  
  supp.res.p <- 
    plot_grid(MeanPeakFit.p, MeanPeakSlope.p, 
              PeakDegree.p, supp.or.p, nrow = 2, 
              align = 'vh')
  
  ggsave(supp.res.p, filename = 'AllSites-A2M-A2P-AccessibilityPeakSummary-Supps.pdf', width = 12, height = 7, useDingbats = F)
  ggsave(supp.res.p, filename = 'AllSites-A2M-A2P-AccessibilityPeakSummary-Supps.png', width = 12, height = 7)
  
  ggsave(NumPeaks.p, filename = 'AllSites-NumberOfPeaks-A2M-A2P-AllRedo.pdf', width = 6, height = 6, useDingbats = F)
  ggsave(mean.peaks.p, filename = 'AllSites-MeanPeakSummary-A2M-A2P-AllRedo.pdf', width = 12, height = 12, useDingbats = F)
  
  # Max fitness for survival isn't really useful, so only plot for the other measures
  if(fitness != "Survival"){
    ggsave(max.peaks.p, filename = 'AllSites-MaxPeakSummary-A2M-A2P-AllRedo.pdf', width = 12, height = 12, useDingbats = F)
  }
  
  
  or.res <- 
    rbind(NumNodes.or,
          NumPaths.or,
          RelNumPaths.or,
          MinDist.or,
          NumPeaks.or,
          RelNumPathsToPeaks.or,
          MinDistToPeak.or)
  
  or.res[,-c(1:2)] <- 
    apply(or.res[,-c(1:2)], 2, 
          function(x){round(x, digits = 4)})
  or.res$OR.CI <- 
    apply(or.res[,-c(1:2)], 1,
          function(x){paste0(x[1], ": (", x[2], ", ", x[3], ")")})
  or.res$LRT.P <- 
    apply(data.frame(or.res[,6]), 1,
          function(x){if(x == 0){x <- "< 0.0001"}})
  
  or.res <- or.res[,c(1:2, 7, 6)]
  
  # Combine with means/SEs
  main.res <- cbind(main.res, Results[,c(5:6,31:32,23:24)])
  main.res[,-1] <- 
    apply(main.res[,-1], 2, 
          function(x){format(round(x, digits = 4), scientific = F)})
  or.res$AtoP.MeanSE <- NA
  or.res$AtoM.MeanSE <- NA
  
  
  for(n in seq(from = 2, to = 14, by = 2)){
    if(n == 2){
      counter <- 1
    }else{
      counter <- counter + 1
    }
    m.mean <- main.res[1,n]
    m.se <- main.res[1,n+1]
    p.mean <- main.res[2,n]
    p.se <- main.res[2,n+1]
    or.res[counter,5] <- paste0(p.mean, ": ", p.se)
    or.res[counter,6] <- paste0(m.mean, ": ", m.se)
  }
  or.res <- apply(or.res, 2, as.character)
  write.table(or.res, file = paste0(plotdir, fitness, '-BroadAccess-A2M-A2P-MainRes-OddsRatio-Vals.txt'), 
              col.names = T, row.names = F, quote = F, sep = '\t')
}

