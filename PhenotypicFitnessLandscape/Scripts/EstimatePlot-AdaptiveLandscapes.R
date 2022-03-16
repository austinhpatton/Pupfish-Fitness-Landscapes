# Estimate adaptive landscapes using generalized additive models for both survival and growth
# Then, follow up by including genotypes that predict growth to these best-fit models

library(car)
library(fields)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(reshape)
library(ggnewscale)
library(mgcViz)
library(MuMIn)
library(MASS)
library(DAAG)
library(dplyr)
library(purrr)
library(cowplot)
library(tidymv)
library(qpcR)
library(metR)
library(formula.tools)
library(scico)
library(mgcv)
library(ggdist)


# This will estimate and plot adaptive landscapes for the sequenced hybrids. 
setwd('~/Dropbox/Research/Martin-Berkeley-Postdoc/Pupfish-Fitness-Genomics/Results/AdaptiveLandscapes/')

# Read in all data
dat <- read.table('../../Data/Pupfish-final-trait-measurements-growth-admix.txt', 
                  header = T, check.names = F)

a.col <- '#D2AF81FF'
m.col <- '#709AE1FF'
p.col <- '#F05C3BFF'
h.col <- '#8A9197FF' 
cols <- c(a.col, m.col, p.col, h.col)

# Separate the parents and hybrids
parents <- dat[which(dat$Experiment == 'Parental'),]
parents <- droplevels(parents)
hybs <- dat[which(dat$Experiment %in% c('RTM1', 'RTM2')),]
hybs <- droplevels(hybs)
hybs$Survivorship <-
  gsub(pattern = 'NonSurvivor', replacement = 0, hybs$Survivorship) %>% 
  gsub(pattern = 'Survivor', replacement = 1, .) %>%
  as.integer(.)

# Now, perform an LDA to discriminate among species, without accounting for lake
# First RTM1
ld.spp <- 
  lda(Species ~ SL + nose + foresnout + bellylen + snoutlen + jaw2pect + 
        pmx2add + pmxlen + jawlen + foreeyewidth + eyewidth + eyeht +
        headht + suspensorium + adductorht + ad2pect + pectinsertion +
        analtocaudal + caudalpedht + dorsaltocaudal + bodydepth +
        dorsalfacelen + eyetosnout + headlen + innereyetosnout + 
        cranialwidth + hindeyewidth + headwidth + lowereyeangle + 
        nasalangle + topeyeangle,
      parents[,-c(1:4)], CV = T)

confusion(actual = parents[,5], predicted = ld.spp$class)

#Now, try accounting for lake as well
parents$SppLake <- paste0(parents[,4], '.', parents[,5])

ld.spp.lake <- 
  lda(SppLake ~ SL + nose + foresnout + bellylen + snoutlen + jaw2pect + 
        pmx2add + pmxlen + jawlen + foreeyewidth + eyewidth + eyeht +
        headht + suspensorium + adductorht + ad2pect + pectinsertion +
        analtocaudal + caudalpedht + dorsaltocaudal + bodydepth +
        dorsalfacelen + eyetosnout + headlen + innereyetosnout + 
        cranialwidth + hindeyewidth + headwidth + lowereyeangle + 
        nasalangle + topeyeangle,
      parents[,-c(1:5)], CV = T)

confusion(actual = parents$SppLake, predicted = ld.spp.lake$class)

# Not actually any added value here - a bit worse. So we're going to proceed with LDA that doesn't include lake.
ld.spp <- 
  lda(Species ~ SL + nose + foresnout + bellylen + snoutlen + jaw2pect + 
        pmx2add + pmxlen + jawlen + foreeyewidth + eyewidth + eyeht +
        headht + suspensorium + adductorht + ad2pect + pectinsertion +
        analtocaudal + caudalpedht + dorsaltocaudal + bodydepth +
        dorsalfacelen + eyetosnout + headlen + innereyetosnout + 
        cranialwidth + hindeyewidth + headwidth + lowereyeangle + 
        nasalangle + topeyeangle,
      parents[,-c(1:4,37)])

ld.pure <- 
  predict(ld.spp, parents[,-c(1:5, 37)])

# Now, let's predict LD scores for the hybrids
traits <- c('SL', 'nose', 'foresnout', 'bellylen', 'snoutlen', 'jaw2pect', 
            'pmx2add', 'pmxlen', 'jawlen', 'foreeyewidth', 'eyewidth', 'eyeht',
            'headht', 'suspensorium', 'adductorht', 'ad2pect', 'pectinsertion',
            'analtocaudal', 'caudalpedht', 'dorsaltocaudal', 'bodydepth',
            'dorsalfacelen', 'eyetosnout', 'headlen', 'innereyetosnout', 
            'cranialwidth', 'hindeyewidth', 'headwidth', 'lowereyeangle', 
            'nasalangle', 'topeyeangle')

hyb.preds <- predict(ld.spp, newdata=hybs[,traits])

# Now, let's do some simple thin plate splines to see if there's anything here!
source('../../Scripts/GetGamVisFunction.R')


# Now, let's plot
# a.col <- "#91bfdb"
# m.col <- "#7fbf7b"
# p.col <- "#fc8d59"
# parent.cols <- c(a.col, m.col, p.col)

a.col <- '#709AE1FF'
m.col <- '#D2AF81FF'
p.col <- '#FD7446FF'
h.col <- '#8A9197FF' 
parent.cols <- c(m.col, a.col, p.col)


# We're also going to turn this step into a function since we'll use it repeatedly.
prep.parents <- 
  function(parents, parental.preds){
    if(dim(parental.preds$x)[2] > 2){
      preds <- 
        data.frame('Species' = parents$Species, 
                   'LD1' = parental.preds$x[,1],
                   'LD2' = parental.preds$x[,2],
                   'LD3' = parental.preds$x[,3],
                   'LD4' = parental.preds$x[,4])
    } else {
      preds <- 
        data.frame('Species' = parents$Species, 
                   'LD1' = parental.preds$x[,1],
                   'LD2' = parental.preds$x[,2])
    }
    
    
    preds$Species <-
      gsub('variegatus', 'Generalist', preds$Species) %>%
      gsub('brontotheroides', 'Molluscivore', .) %>%
      gsub('desquamator', 'Scale-eater', .)
    
    preds$Lake <-
      parents$Lake
    
    return(preds)
  }

prep.hybs <-
  function(hybs, hyb.preds){
    if(dim(hyb.preds$x)[2] > 2){
      preds <- 
        data.frame('ID' = hybs$ID,
                   'Experiment' = hybs$Experiment, 
                   'LD1' = hyb.preds$x[,1],
                   'LD2' = hyb.preds$x[,2],
                   'LD3' = hyb.preds$x[,3],
                   'LD4' = hyb.preds$x[,4], 
                   'Lake' = hybs$Lake,
                   'Growth' = hybs$Growth,
                   'Molluscivore' = hybs$Molluscivore,
                   'ScaleEater' = hybs$`Scale-eater`)
    } else {
      preds <- 
        data.frame('ID' = hybs$ID,
                   'Experiment' = hybs$Experiment, 
                   'LD1' = hyb.preds$x[,1],
                   'LD2' = hyb.preds$x[,2], 
                   'Lake' = hybs$Lake,
                   'Growth' = hybs$Growth,
                   'Molluscivore' = hybs$Molluscivore,
                   'ScaleEater' = hybs$`Scale-eater`)
    }
    
    preds$survival.ID <- 
      hybs$Survivorship
    
    preds$survival.ID <- 
      gsub(pattern = 1, replacement = 'Survivor', preds$survival.ID) %>%
      gsub(pattern ='0', replacement = 'Non-Survivor', .)
    
    return(preds)
  }

# Now, we're ready to go. Let's plot the simple thin plate splines for LD1/LD2 for each 

# 1) RTM1 traits, all fish
parent.preds <- 
  prep.parents(parents = parents, parental.preds = ld.pure)

hyb.preds <- 
  prep.hybs(hybs = hybs, hyb.preds = hyb.preds)

hyb.preds$survival <- as.integer(
  gsub(pattern ='Non-Survivor', replacement = 0, hyb.preds$survival.ID) %>%
    gsub(pattern = 'Survivor', replacement = 1, .))

# write.table(hyb.preds, file = '../Hybrid-GamDat.tsv', quote = F, sep = '\t', row.names = F)
hyb.preds$Experiment <- factor(hyb.preds$Experiment)
hyb.preds$Lake <- factor(hyb.preds$Lake)

m1 <- gam(survival ~ s(LD1, LD2) + Experiment + Lake, data = hyb.preds, method = 'GCV.Cp', family = binomial)
m2 <- gam(survival ~ s(LD1, LD2) + Experiment * Lake, data = hyb.preds, method = 'GCV.Cp', family = binomial)
m3 <- gam(survival ~ s(LD1, LD2) + s(LD1) + s(LD2) + Experiment + Lake, data = hyb.preds, method = 'GCV.Cp', family = binomial)
m4 <- gam(survival ~ s(LD1, LD2) + s(LD1) + s(LD2) + Experiment * Lake, data = hyb.preds, method = 'GCV.Cp', family = binomial)
m5 <- gam(survival ~ s(LD1, LD2) + s(LD1, Experiment, bs = 'fs') + s(LD2, Experiment, bs = 'fs') + Lake, data = hyb.preds, method = 'GCV.Cp', family = binomial)
m6 <- gam(survival ~ s(LD1, LD2) + s(LD1, Lake, bs = 'fs') + s(LD2, Lake, bs = 'fs') + Experiment, data = hyb.preds, method = 'GCV.Cp', family = binomial)


# And summarize for reference later
aics <- 
  c(AICc(m1), AICc(m2), AICc(m3), 
    AICc(m4), AICc(m5), AICc(m6))

weights  <- akaike.weights(aics)

fits <- 
  list(m1, m2, m3, m4, m5, m6)

mods <- 
  c(as.character(m1$formula), as.character(m2$formula), as.character(m3$formula), 
    as.character(m4$formula), as.character(m5$formula), as.character(m6$formula))

mod.table <- 
  data.frame('Model' = mods,
             'AICc' = aics, 
             'deltaAICc' = weights$deltaAIC,
             'AkWeights' = weights$weights)

mod.table
# write.table(x = mod.table, file = 'Morpho-SurvivalGAM-ModTable.tsv', 
#             col.names = T, quote = F, sep = '\t', row.names = F)

# Now get the best model and prepare to visualize
best.mod <- which(aics == min(aics))
best.mod
best.mod.res <- get.gam.viz(fits[[best.mod]], n.grid = 100, too.far = 0.075)
best.mod <- fits[[best.mod]]
# save(best.mod, file = 'BestMorpho-SurvivalGAM.Rsave')

parent.preds$Species <- 
  factor(parent.preds$Species, 
         levels = c('Generalist', 'Molluscivore', 'Scale-eater'))
colnames(hyb.preds)[c(8,9)] <- c('Scale-eater', 'Survival')
colnames(best.mod.res)[3] <- 'Survival Probability'
best.surv.mod.p <- 
  ggplot(best.mod.res, aes(x = LD1, y = LD2)) +
  scale_fill_scico(palette = 'davos', direction = 1, na.value = 'transparent')+#, limits = c(0.05, 0.4)) +
  geom_raster(aes(fill = `Survival Probability`)) +
  geom_contour(color = 'black', aes(z = `Survival Probability`), size = 2, alpha = 0.5) + 
  geom_contour(color = 'white', aes(z = `Survival Probability`), size = 1, alpha = 0.8) + 
  geom_text_contour(aes(z = `Survival Probability`), stroke = 0.1) + 
  new_scale_fill() +
  geom_point(data = hyb.preds, aes(x = LD1, y = LD2, fill = Survival, size = Survival), 
             pch = 21, alpha = 0.75) +
  scale_fill_manual(values = c('grey40', 'white')) +
  scale_size_manual(values = c(1,2.25)) +
  new_scale_fill() +
  stat_ellipse(data = parent.preds, level = c(0.95), 
               aes(x = LD1, y = LD2, color = Species), 
               geom = 'polygon', size = 1, type = 'norm',
               alpha = 0.2, lty = 3) +
  stat_ellipse(data = parent.preds, level = c(0.5), 
               aes(x = LD1, y = LD2, color = Species), 
               geom = 'polygon', size = 1, type = 'norm',
               alpha = 0.2, lty = 1) +
  geom_point(data = parent.preds, aes(x = LD1, y = LD2),
             size = 3.5) +
  geom_point(data = parent.preds, aes(x = LD1, y = LD2, color = Species),
             size = 2.5) +
  scale_shape_manual(values = 21) +
  scale_fill_manual(values = parent.cols) +
  scale_color_manual(values = parent.cols) +
  # coord_cartesian(xlim = c(-10, 10), 
  #                 ylim = c(-10, 6)) +
  theme_bw(base_size = 18)
best.surv.mod.p
# ggsave(plot = best.surv.mod.p, 'Adapt-Landscape-Survival.pdf', width = 10, height = 8, useDingbats = F)


##################################################################################################
# What if we use growth rate data? Note that these results include non-survivors, coded as zeros #
##################################################################################################
hyb.preds$Experiment <- factor(hyb.preds$Experiment)
hyb.preds$Lake <- factor(hyb.preds$Lake)

m1 <- gam(Growth ~ s(LD1, LD2) + Experiment + Lake, data = hyb.preds, method = 'GCV.Cp')
m2 <- gam(Growth ~ s(LD1, LD2) + Experiment * Lake, data = hyb.preds, method = 'GCV.Cp')
m3 <- gam(Growth ~ s(LD1, LD2) + s(LD1) + s(LD2) + Experiment + Lake, data = hyb.preds, method = 'GCV.Cp')
m4 <- gam(Growth ~ s(LD1, LD2) + s(LD1) + s(LD2) + Experiment * Lake, data = hyb.preds, method = 'GCV.Cp')
m5 <- gam(Growth ~ s(LD1, LD2) + s(LD1, Experiment, bs = 'fs') + s(LD2, Experiment, bs = 'fs') + Lake, data = hyb.preds, method = 'GCV.Cp')
m6 <- gam(Growth ~ s(LD1, LD2) + s(LD1, Lake, bs = 'fs') + s(LD2, Lake, bs = 'fs') + Experiment, data = hyb.preds, method = 'GCV.Cp')

# And summarize for reference later
aics <- 
  c(AICc(m1), AICc(m2), AICc(m3), 
    AICc(m4), AICc(m5), AICc(m6))

weights  <- akaike.weights(aics)

fits <- 
  list(m1, m2, m3, m4, m5, m6)

mods <- 
  c(as.character(m1$formula), as.character(m2$formula), as.character(m3$formula), 
    as.character(m4$formula), as.character(m5$formula), as.character(m6$formula))

mod.table <- 
  data.frame('Model' = mods,
             'AICc' = aics, 
             'deltaAICc' = weights$deltaAIC,
             'AkWeights' = weights$weights)

mod.table
# write.table(x = mod.table, file = 'Morpho-Growth-WithNonSurvs-GAM-ModTable.tsv', 
            # col.names = T, quote = F, sep = '\t', row.names = F)

# Now get the best model and prepare to visualize
best.mod <- which(aics == min(aics))
best.mod
best.mod.res <- get.gam.viz(fits[[best.mod]], n.grid = 100, too.far = 0.075)
colnames(best.mod.res)[3] <- 'Value'
best.mod <- fits[[best.mod]]
# save(best.mod, file = 'BestMorpho-Growth-WithNonSurvs-GAM.Rsave')

parent.preds$Species <- 
  factor(parent.preds$Species, 
         levels = c('Generalist', 'Molluscivore', 'Scale-eater'))

colnames(best.mod.res)[3] <- 'Predicted Growth'
pred <- guide_legend('Predicted Growth')
obs <- guide_legend('Observed Growth')
spp.cols <- guide_legend('Species')

upper <- max(na.omit(best.mod.res$`Predicted Growth`))

best.growth.withdead.mod.p <- 
  ggplot(best.mod.res, aes(x = LD1, y = LD2)) +
  scale_fill_scico(palette = 'davos', direction = 1, na.value = 'transparent', limits = c(0.00, upper)) +
  geom_raster(aes(fill = `Predicted Growth`)) +
  geom_contour(color = 'black', aes(z = `Predicted Growth`), size = 2, alpha = 0.5) + 
  geom_contour(color = 'white', aes(z = `Predicted Growth`), size = 1, alpha = 0.8) + 
  geom_text_contour(aes(z = `Predicted Growth`), stroke = 0.1) + 
  guides(fill = pred) +
  new_scale_fill() +
  geom_point(data = hyb.preds, aes(x = LD1, y = LD2, fill = Growth, size = Growth), 
             pch = 21, alpha = 0.75) +
  scale_fill_continuous(low = 'black', high = 'white', 
                        limits = c(0, 1.5), breaks = c(0,0.5,1,1.5)) +
  scale_size_continuous(range = c(1,3.25), 
                        limits = c(0, 1.5), breaks = c(0,0.5,1,1.5)) +
  guides(fill = obs, size = obs) +
  stat_ellipse(data = parent.preds, level = c(0.95), 
               aes(x = LD1, y = LD2, color = Species), 
               geom = 'polygon', size = 1, type = 'norm',
               alpha = 0.2, lty = 3) +
  stat_ellipse(data = parent.preds, level = c(0.5), 
               aes(x = LD1, y = LD2, color = Species), 
               geom = 'polygon', size = 1, type = 'norm',
               alpha = 0.2, lty = 1) +
  geom_point(data = parent.preds, aes(x = LD1, y = LD2),
             size = 3.5) +
  geom_point(data = parent.preds, aes(x = LD1, y = LD2, color = Species),
             size = 2.5) +
  scale_shape_manual(values = 21) +
  scale_color_manual(values = parent.cols) +
  guides(color = spp.cols) +
  theme_bw(base_size = 18) 
best.growth.withdead.mod.p

# ggsave(plot = best.growth.withdead.mod.p, 'Adapt-Landscape-Growth-wDead.pdf', width = 10, height = 8, useDingbats = F)

####################################################################
# What if we use growth rate data? These data exclude non-survivors #
####################################################################
orig <- hyb.preds
hyb.preds <- hyb.preds[which(hyb.preds$Growth > 0),]

m1 <- gam(Growth ~ s(LD1, LD2) + Experiment + Lake, data = hyb.preds, method = 'GCV.Cp')
m2 <- gam(Growth ~ s(LD1, LD2) + Experiment * Lake, data = hyb.preds, method = 'GCV.Cp')
m3 <- gam(Growth ~ s(LD1, LD2) + s(LD1) + s(LD2) + Experiment + Lake, data = hyb.preds, method = 'GCV.Cp')
m4 <- gam(Growth ~ s(LD1, LD2) + s(LD1) + s(LD2) + Experiment * Lake, data = hyb.preds, method = 'GCV.Cp')
m5 <- gam(Growth ~ s(LD1, LD2) + s(LD1, Experiment, bs = 'fs') + s(LD2, Experiment, bs = 'fs') + Lake, data = hyb.preds, method = 'GCV.Cp')
m6 <- gam(Growth ~ s(LD1, LD2) + s(LD1, Lake, bs = 'fs') + s(LD2, Lake, bs = 'fs') + Experiment, data = hyb.preds, method = 'GCV.Cp')

# And summarize for reference later
aics <- 1
  c(AICc(m1), AICc(m2), AICc(m3), 
    AICc(m4), AICc(m5), AICc(m6))

weights  <- akaike.weights(aics)

fits <- 
  list(m1, m2, m3, m4, m5, m6)

mods <- 
  c(as.character(m1$formula), as.character(m2$formula), as.character(m3$formula), 
    as.character(m4$formula), as.character(m5$formula), as.character(m6$formula))

mod.table <- 
  data.frame('Model' = mods,
             'AICc' = aics, 
             'deltaAICc' = weights$deltaAIC,
             'AkWeights' = weights$weights)

mod.table
# write.table(x = mod.table, file = 'Morpho-Growth-GAM-ModTable.tsv', 
            # col.names = T, quote = F, sep = '\t', row.names = F)

# Now get the best model and prepare to visualize
best.mod <- which(aics == min(aics))
best.mod 
best.mod.res <- get.gam.viz(fits[[best.mod]], n.grid = 100, too.far = 0.075)
colnames(best.mod.res)[3] <- 'Value'
best.mod <- fits[[best.mod]]
# save(best.mod, file = 'BestMorpho-Growth-GAM.Rsave')

parent.preds$Species <- 
  factor(parent.preds$Species, 
         levels = c('Generalist', 'Molluscivore', 'Scale-eater'))


colnames(best.mod.res)[3] <- 'Predicted Growth'
pred <- guide_legend('Predicted Growth')
obs <- guide_legend('Observed Growth')
spp.cols <- guide_legend('Species')
best.growth.mod.p <- 
  ggplot(best.mod.res, aes(x = LD1, y = LD2)) +
  scale_fill_scico(palette = 'davos', direction = 1, na.value = 'transparent')+#, limits = c(0.05, 0.4)) +
  geom_raster(aes(fill = `Predicted Growth`)) +
  geom_contour(color = 'black', aes(z = `Predicted Growth`), size = 2, alpha = 0.5) + 
  geom_contour(color = 'white', aes(z = `Predicted Growth`), size = 1, alpha = 0.8) + 
  geom_text_contour(aes(z = `Predicted Growth`), stroke = 0.1) + 
  guides(fill = pred) +
  new_scale_fill() +
  geom_point(data = hyb.preds, aes(x = LD1, y = LD2, fill = Growth, size = Growth), 
             pch = 21, alpha = 0.75) +
  scale_fill_continuous(low = 'black', high = 'white', 
                        limits = c(0, 1.5), breaks = c(0,0.5,1,1.5)) +
  scale_size_continuous(range = c(1,3.25), 
                        limits = c(0, 1.5), breaks = c(0,0.5,1,1.5)) +
  guides(fill = obs, size = obs) +
  stat_ellipse(data = parent.preds, level = c(0.95), 
               aes(x = LD1, y = LD2, color = Species), 
               geom = 'polygon', size = 1, type = 'norm',
               alpha = 0.2, lty = 3) +
  stat_ellipse(data = parent.preds, level = c(0.5), 
               aes(x = LD1, y = LD2, color = Species), 
               geom = 'polygon', size = 1, type = 'norm',
               alpha = 0.2, lty = 1) +
  geom_point(data = parent.preds, aes(x = LD1, y = LD2),
             size = 3.5) +
  geom_point(data = parent.preds, aes(x = LD1, y = LD2, color = Species),
             size = 2.5) +
  scale_shape_manual(values = 21) +
  scale_color_manual(values = parent.cols) +
  guides(color = spp.cols) +
  theme_bw(base_size = 18) 
best.growth.mod.p

# ggsave(plot = best.growth.mod.p, 'Adapt-Landscape-Growth.pdf', width = 10, height = 8, useDingbats = F)

#############################################################################################
# And if we include genotypes significantly associated with growth from GEMMA LMM?          #
# Note that these SNPs were identified using survivors only                                 #
#############################################################################################
#hyb.preds <- orig
setwd('~/Dropbox/Research/Martin-Berkeley-Postdoc/Pupfish-Fitness-Genomics/Results/AdaptiveLandscapes/')

dat <- read.table('../../Data/Pupfish-final-trait-measurements-growth-admix.txt', 
                  header = T, check.names = F)
hybs <- dat[which(dat$Experiment %in% c('RTM1', 'RTM2')),]

genos <- read.table('../../Results/GEMMA/Traits/SignificantSites/growth/growthOnly-BonfSig-sites-Hybs-012-matrix.txt',
                    sep = '\t', header = T)

sites <- colnames(genos)[-1]

prefix <- "Site"
suffix <- seq(1:4)
colnames(genos) <- c('ID', paste(prefix, suffix, sep =''))

hyb.preds <- merge(hyb.preds, genos, by = 'ID')
hyb.preds[,11:14] <- hyb.preds[,11:14] %>% mutate_if(is.numeric, as.factor)
for(i in 11:14){
  hyb.preds[,i] <-
    gsub('wild', 0, hyb.preds[,i]) %>%
    gsub('single', 1, .) %>%
    gsub('double', 2, .) %>%
    factor(., ordered = T)
}

m1 <- gam(Growth ~ s(LD1, LD2) + Experiment * Lake, data = hyb.preds, method = 'GCV.Cp')
m2 <- gam(Growth ~ s(LD1, LD2) + Experiment * Lake + s(Site1, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
m3 <- gam(Growth ~ s(LD1, LD2) + Experiment * Lake + s(Site2, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
m4 <- gam(Growth ~ s(LD1, LD2) + Experiment * Lake + s(Site3, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
m5 <- gam(Growth ~ s(LD1, LD2) + Experiment * Lake + s(Site4, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
m6 <- gam(Growth ~ s(LD1, LD2) + Experiment * Lake + s(Site1, bs = 'fs') + s(Site2, bs = 'fs') + s(Site3, bs = 'fs') + s(Site4, bs = 'fs'), 
          data = hyb.preds, method = 'GCV.Cp')
# And try alternatives
m7 <- gam(Growth ~ s(LD1, LD2) + Experiment * Lake + s(Site1, bs = 'fs') + s(Site3, bs = 'fs') + s(Site4, bs = 'fs'), 
          data = hyb.preds, method = 'GCV.Cp')
m8 <- gam(Growth ~ s(LD1, LD2) + Experiment * Lake + s(Site1, bs = 'fs') + s(Site2, bs = 'fs') + s(Site4, bs = 'fs'), 
          data = hyb.preds, method = 'GCV.Cp')
m9 <- gam(Growth ~ s(LD1, LD2) + Experiment * Lake + s(Site1, bs = 'fs') + s(Site2, bs = 'fs') + s(Site3, bs = 'fs'), 
          data = hyb.preds, method = 'GCV.Cp')
m10 <- gam(Growth ~ s(LD1, LD2) + Experiment * Lake + s(Site1, bs = 'fs') + s(Site2, bs = 'fs'),
           data = hyb.preds, method = 'GCV.Cp')
m11 <- gam(Growth ~ s(LD1, LD2) + Experiment * Lake + s(Site1, bs = 'fs') + s(Site3, bs = 'fs'),
           data = hyb.preds, method = 'GCV.Cp')
m12 <- gam(Growth ~ s(LD1, LD2) + Experiment * Lake + s(Site1, bs = 'fs') + s(Site4, bs = 'fs'),
           data = hyb.preds, method = 'GCV.Cp')
m13 <- gam(Growth ~ s(LD1, LD2) + Experiment * Lake + s(Site2, bs = 'fs') + s(Site3, bs = 'fs'),
           data = hyb.preds, method = 'GCV.Cp')
m14 <- gam(Growth ~ s(LD1, LD2) + Experiment * Lake + s(Site2, bs = 'fs') + s(Site4, bs = 'fs'),
           data = hyb.preds, method = 'GCV.Cp')
m15 <- gam(Growth ~ s(LD1, LD2) + Experiment * Lake + s(Site3, bs = 'fs') + s(Site4, bs = 'fs'),
           data = hyb.preds, method = 'GCV.Cp')

# And summarize for reference later
aics <- c(AICc(m1), AICc(m2), AICc(m3), 
          AICc(m4), AICc(m5), AICc(m6),
          AICc(m7), AICc(m8), AICc(m9),
          AICc(m10), AICc(m11), AICc(m12),
          AICc(m13), AICc(m14), AICc(m15))


weights  <- akaike.weights(aics)

fits <- 
  list(m1, m2, m3,  
       m4, m5, m6,
       m7, m8, m9,
       m10, m11, m12,
       m13, m14, m15) 

mods <- 
  c(as.character(m1$formula), as.character(m2$formula), as.character(m3$formula), 
    as.character(m4$formula), as.character(m5$formula), as.character(m6$formula),
    as.character(m7$formula), as.character(m8$formula), as.character(m9$formula),
    as.character(m10$formula), as.character(m11$formula), as.character(m12$formula),
    as.character(m13$formula), as.character(m14$formula), as.character(m15$formula))

mod.table <- 
  data.frame('Model' = mods,
             'AICc' = aics, 
             'deltaAICc' = weights$deltaAIC,
             'AkWeights' = weights$weights)

mod.table
# write.table(x = mod.table, file = 'Morpho-Geno-Growth-GAM-ModTable.tsv',
# col.names = T, quote = F, sep = '\t', row.names = F)


# Now get the best model and prepare to visualize
best.mod <- which(aics == min(aics))
best.mod
best.mod.res <- get.gam.viz(fits[[best.mod]], n.grid = 100, too.far = 0.075)
best.mod <- fits[[best.mod]]
# save(best.mod, file = 'BestMorphoGeno-Growth-GAM.Rsave')

viz <- getViz(best.mod)

#pdf('Growth-IndSNP-Splines.pdf')
print(plot(viz, select = 2), pages = 1)
#dev.off()

#pdf('Growth-IndSNP-LD-Contours.pdf', height = 24, width = 12)
par(mfrow = c(1, 2),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
vis.gam(best.mod, view = c('LD1', 'Site3'), theta = 230)
vis.gam(best.mod, view = c('LD2', 'Site3'), theta = 230)
#dev.off()

parent.preds$Species <- 
  factor(parent.preds$Species, 
         levels = c('Generalist', 'Molluscivore', 'Scale-eater'))
colnames(best.mod.res)[3] <- 'Predicted Growth'

pred <- guide_legend('Predicted Growth')
obs <- guide_legend('Observed Growth')
spp.cols <- guide_legend('Species')
best.geno.growth.mod.p <- 
  ggplot(best.mod.res, aes(x = LD1, y = LD2)) +
  scale_fill_scico(palette = 'davos', direction = 1, na.value = 'transparent')+#, limits = c(0.05, 0.4)) +
  geom_raster(aes(fill = `Predicted Growth`)) +
  geom_contour(color = 'black', aes(z = `Predicted Growth`), size = 2, alpha = 0.5) + 
  geom_contour(color = 'white', aes(z = `Predicted Growth`), size = 1, alpha = 0.8) + 
  geom_text_contour(aes(z = `Predicted Growth`), stroke = 0.1) + 
  guides(fill = pred) +
  new_scale_fill() +
  geom_point(data = hyb.preds, aes(x = LD1, y = LD2, fill = Growth, size = Growth), 
             pch = 21, alpha = 0.75) +
  scale_fill_continuous(low = 'black', high = 'white', 
                        limits = c(0, 1.5), breaks = c(0,0.5,1,1.5)) +
  scale_size_continuous(range = c(1,3.25), 
                        limits = c(0, 1.5), breaks = c(0,0.5,1,1.5)) +
  guides(fill = obs, size = obs) +
  stat_ellipse(data = parent.preds, level = c(0.95), 
               aes(x = LD1, y = LD2, color = Species), 
               geom = 'polygon', size = 1, type = 'norm',
               alpha = 0.2, lty = 3) +
  stat_ellipse(data = parent.preds, level = c(0.5), 
               aes(x = LD1, y = LD2, color = Species), 
               geom = 'polygon', size = 1, type = 'norm',
               alpha = 0.2, lty = 1) +
  geom_point(data = parent.preds, aes(x = LD1, y = LD2),
             size = 3.5) +
  geom_point(data = parent.preds, aes(x = LD1, y = LD2, color = Species),
             size = 2.5) +
  scale_shape_manual(values = 21) +
  scale_color_manual(values = parent.cols) +
  guides(color = spp.cols) +
  theme_bw(base_size = 18) 
best.geno.growth.mod.p

#ggsave(plot = best.geno.growth.mod.p, 'Adapt-Landscape-Geno-Growth.pdf', width = 10, height = 8, useDingbats = F)

#############################################################################################
# And if we include genotypes significantly associated with growth from GEMMA LMM?          #
# Note that these SNPs were identified using survivors and non-survivors (coded as 0)       #
#############################################################################################
hyb.preds <- orig
setwd('~/Dropbox/Research/Martin-Berkeley-Postdoc/Pupfish-Fitness-Genomics/Results/AdaptiveLandscapes/')

dat <- read.table('../../Data/Pupfish-final-trait-measurements-growth-admix.txt', 
                  header = T, check.names = F)
hybs <- dat[which(dat$Experiment %in% c('RTM1', 'RTM2')),]

genos <- read.table('../GEMMA/Traits/SignificantSites/growth/growth-BonfSig-sites-mutants.txt',
                    sep = '\t', header = T)
genos <- genos[,-c(4:5, 10:11,13)]

sites <- colnames(genos)[-1]

prefix <- "Site"
suffix <- seq(1:10)
colnames(genos) <- c('ID', paste(prefix, suffix, sep =''))

hyb.preds <- merge(hyb.preds, genos, by = 'ID')
for(i in 11:20){
hyb.preds[,i] <-
  gsub('wild', 0, hyb.preds[,i]) %>%
  gsub('single', 1, .) %>%
  gsub('double', 2, .) %>%
  factor(., ordered = T)
}

m1 <- gam(Growth ~ s(LD1, LD2) + Experiment + Lake, data = hyb.preds, method = 'GCV.Cp')
m2 <- gam(Growth ~ s(LD1, LD2) + Experiment + Lake + s(Site1, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
m3 <- gam(Growth ~ s(LD1, LD2) + Experiment + Lake + s(Site2, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
m4 <- gam(Growth ~ s(LD1, LD2) + Experiment + Lake + s(Site3, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
m5 <- gam(Growth ~ s(LD1, LD2) + Experiment + Lake + s(Site4, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
m6 <- gam(Growth ~ s(LD1, LD2) + Experiment + Lake + s(Site5, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
m7 <- gam(Growth ~ s(LD1, LD2) + Experiment + Lake + s(Site6, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
m8 <- gam(Growth ~ s(LD1, LD2) + Experiment + Lake + s(Site7, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
m9 <- gam(Growth ~ s(LD1, LD2) + Experiment + Lake + s(Site8, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
m10 <- gam(Growth ~ s(LD1, LD2) + Experiment + Lake + s(Site9, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
m11 <- gam(Growth ~ s(LD1, LD2) + Experiment + Lake + s(Site10, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
m12 <- gam(Growth ~ s(LD1, LD2) + Experiment + Lake + s(Site1, bs = 'fs') + s(Site2, bs = 'fs') + 
             s(Site3, bs = 'fs') + s(Site4, bs = 'fs') + s(Site5, bs = 'fs') + s(Site6, bs = 'fs') + s(Site7, bs = 'fs') + 
             s(Site8, bs = 'fs') + s(Site9, bs = 'fs') + s(Site10, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
m13 <- gam(Growth ~ s(LD1, LD2) + Experiment + Lake + s(Site1, bs = 'fs') + s(Site2, bs = 'fs') + s(Site6, bs = 'fs') + 
             s(Site7, bs = 'fs') + s(Site8, bs = 'fs') + s(Site9, bs = 'fs') + s(Site10, bs = 'fs'), 
           data = hyb.preds, method = 'GCV.Cp')

viz <- getViz(m13)
# 
# pdf('Composite-IndSNP-Splines.pdf')
# print(plot(viz, select = 2:8), pages = 1)
# dev.off()

# pdf('Composite-IndSNP-LD-Contours.pdf', height = 24, width = 12)
par(mfrow = c(7, 2),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
vis.gam(m13, view = c('LD2', 'Site1'), theta = 230)

vis.gam(m13, view = c('LD1', 'Site2'), theta = 230)
vis.gam(m13, view = c('LD2', 'Site2'), theta = 230)

vis.gam(m13, view = c('LD1', 'Site6'), theta = 230)
vis.gam(m13, view = c('LD2', 'Site6'), theta = 230)

vis.gam(m13, view = c('LD1', 'Site7'), theta = 230)
vis.gam(m13, view = c('LD2', 'Site7'), theta = 230)

vis.gam(m13, view = c('LD1', 'Site8'), theta = 230)
vis.gam(m13, view = c('LD2', 'Site8'), theta = 230)

vis.gam(m13, view = c('LD1', 'Site9'), theta = 230)
vis.gam(m13, view = c('LD2', 'Site9'), theta = 230)

vis.gam(m13, view = c('LD1', 'Site10'), theta = 230)
vis.gam(m13, view = c('LD2', 'Site10'), theta = 230)
# dev.off()
# 
# 
# m14 <- gam(Growth ~ s(LD1, LD2) + s(Site1, bs = 'fs') + s(Site2, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
# m15 <- gam(Growth ~ s(LD1, LD2) + s(Site1, bs = 'fs') + s(Site6, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
# m16 <- gam(Growth ~ s(LD1, LD2) + s(Site1, bs = 'fs') + s(Site7, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
# m17 <- gam(Growth ~ s(LD1, LD2) + s(Site1, bs = 'fs') + s(Site8, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
# m18 <- gam(Growth ~ s(LD1, LD2) + s(Site1, bs = 'fs') + s(Site9, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
# m19 <- gam(Growth ~ s(LD1, LD2) + s(Site1, bs = 'fs') + s(Site10, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
# m20 <- gam(Growth ~ s(LD1, LD2) + s(Site2, bs = 'fs') + s(Site6, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
# m21 <- gam(Growth ~ s(LD1, LD2) + s(Site2, bs = 'fs') + s(Site7, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
# m22 <- gam(Growth ~ s(LD1, LD2) + s(Site2, bs = 'fs') + s(Site8, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
# m23 <- gam(Growth ~ s(LD1, LD2) + s(Site2, bs = 'fs') + s(Site9, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
# m24 <- gam(Growth ~ s(LD1, LD2) + s(Site2, bs = 'fs') + s(Site10, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
# m25 <- gam(Growth ~ s(LD1, LD2) + s(Site6, bs = 'fs') + s(Site7, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
# m26 <- gam(Growth ~ s(LD1, LD2) + s(Site6, bs = 'fs') + s(Site8, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
# m27 <- gam(Growth ~ s(LD1, LD2) + s(Site6, bs = 'fs') + s(Site9, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
# m28 <- gam(Growth ~ s(LD1, LD2) + s(Site6, bs = 'fs') + s(Site10, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
# m29 <- gam(Growth ~ s(LD1, LD2) + s(Site7, bs = 'fs') + s(Site8, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
# m30 <- gam(Growth ~ s(LD1, LD2) + s(Site7, bs = 'fs') + s(Site9, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
# m31 <- gam(Growth ~ s(LD1, LD2) + s(Site7, bs = 'fs') + s(Site10, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
# m32 <- gam(Growth ~ s(LD1, LD2) + s(Site8, bs = 'fs') + s(Site9, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
# m33 <- gam(Growth ~ s(LD1, LD2) + s(Site8, bs = 'fs') + s(Site10, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
# m34 <- gam(Growth ~ s(LD1, LD2) + s(Site9, bs = 'fs') + s(Site10, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
# 
# m14 <- gam(Growth ~ s(LD1, LD2) + s(Site1, by = Site2, bs = 'fs') + s(Site2, by = Site1, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
# m15 <- gam(Growth ~ s(LD1, LD2) + s(Site1, by = Site6, bs = 'fs') + s(Site6, by = Site1, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
# m16 <- gam(Growth ~ s(LD1, LD2) + s(Site1, by = Site7, bs = 'fs') + s(Site7, by = Site1, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
# m17 <- gam(Growth ~ s(LD1, LD2) + s(Site1, by = Site8, bs = 'fs') + s(Site8, by = Site1, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
# m18 <- gam(Growth ~ s(LD1, LD2) + s(Site1, by = Site9, bs = 'fs') + s(Site9, by = Site1, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
# m19 <- gam(Growth ~ s(LD1, LD2) + s(Site1, by = Site10, bs = 'fs') + s(Site10, by = Site1, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
# m20 <- gam(Growth ~ s(LD1, LD2) + s(Site2, by = Site6, bs = 'fs') + s(Site6, by = Site2, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
# m21 <- gam(Growth ~ s(LD1, LD2) + s(Site2, by = Site7, bs = 'fs') + s(Site7, by = Site2, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
# m22 <- gam(Growth ~ s(LD1, LD2) + s(Site2, by = Site8, bs = 'fs') + s(Site8, by = Site2, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
# m23 <- gam(Growth ~ s(LD1, LD2) + s(Site2, by = Site9, bs = 'fs') + s(Site9, by = Site2, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
# m24 <- gam(Growth ~ s(LD1, LD2) + s(Site2, by = Site10, bs = 'fs') + s(Site10, by = Site2, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
# m25 <- gam(Growth ~ s(LD1, LD2) + s(Site6, by = Site7, bs = 'fs') + s(Site7, by = Site6, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
# m26 <- gam(Growth ~ s(LD1, LD2) + s(Site6, by = Site8, bs = 'fs') + s(Site8, by = Site6, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
# m27 <- gam(Growth ~ s(LD1, LD2) + s(Site6, by = Site9, bs = 'fs') + s(Site9, by = Site6, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
# m28 <- gam(Growth ~ s(LD1, LD2) + s(Site6, by = Site10, bs = 'fs') + s(Site10, by = Site6, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
# m29 <- gam(Growth ~ s(LD1, LD2) + s(Site7, by = Site8, bs = 'fs') + s(Site8, by = Site7, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
# m30 <- gam(Growth ~ s(LD1, LD2) + s(Site7, by = Site9, bs = 'fs') + s(Site9, by = Site7, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
# m31 <- gam(Growth ~ s(LD1, LD2) + s(Site7, by = Site10, bs = 'fs') + s(Site10, by = Site7, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
# m32 <- gam(Growth ~ s(LD1, LD2) + s(Site8, by = Site9, bs = 'fs') + s(Site9, by = Site8, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
# m33 <- gam(Growth ~ s(LD1, LD2) + s(Site8, by = Site10, bs = 'fs') + s(Site10, by = Site8, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
# m34 <- gam(Growth ~ s(LD1, LD2) + s(Site9, by = Site10, bs = 'fs') + s(Site10, by = Site9, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
hyb.preds$haps <- as.factor(hyb.preds$haps[,1])
# One haploytpe 
# And summarize for reference later
aics <- c(AICc(m1), AICc(m2), AICc(m3), AICc(m4), AICc(m5),
          AICc(m6), AICc(m7), AICc(m8), AICc(m9), AICc(m10), 
          AICc(m11), AICc(m12), AICc(m13))


weights  <- akaike.weights(aics)

fits <- 
  list(m1, m2, m3, m4, m5,
       m6, m7, m8, m9, m10, 
       m11, m12, m13)


mods <- 
  c(as.character(m1$formula), as.character(m2$formula), as.character(m3$formula), 
    as.character(m4$formula), as.character(m5$formula), as.character(m6$formula),
    as.character(m7$formula), as.character(m8$formula), as.character(m9$formula),
    as.character(m10$formula), as.character(m11$formula), as.character(m12$formula),
    as.character(m13$formula))

mod.table <- 
  data.frame('Model' = mods,
             'AICc' = aics, 
             'deltaAICc' = weights$deltaAIC,
             'AkWeights' = weights$weights)

mod.table
# write.table(x = mod.table, file = 'Morpho-Geno-Composite-GAM-ModTable.tsv', 
#             col.names = T, quote = F, sep = '\t', row.names = F)

# Now get the best model and prepare to visualize
best.mod <- which(aics == min(aics))
best.mod
best.mod.res <- get.gam.viz(fits[[best.mod]], n.grid = 100)
best.mod <- fits[[best.mod]]
# save(best.mod, file = 'BestMorphoGeno-Composite-GAM.Rsave')

parent.preds$Species <- 
  factor(parent.preds$Species, 
         levels = c('Generalist', 'Molluscivore', 'Scale-eater'))
colnames(best.mod.res)[3] <- 'Predicted Growth'

pred <- guide_legend('Predicted Growth')
obs <- guide_legend('Observed Growth')
spp.cols <- guide_legend('Species')

# Set color scale relative to be the same as the genotype free composite fitness landscape
upper <- max(na.omit(best.growth.withdead.mod.p$data$`Predicted Growth`))
lower <- min(na.omit(best.growth.withdead.mod.p$data$`Predicted Growth`))

best.geno.growth.wDead.mod.p <- 
  ggplot(best.mod.res, aes(x = LD1, y = LD2)) +
  # scale_fill_gradientn(colours = pal, na.value = 'transparent') +
  scale_fill_scico(palette = 'davos', direction = 1, na.value = 'transparent', limits = c(0, upper)) +
  geom_raster(aes(fill = `Predicted Growth`)) +
  geom_contour(color = 'black', aes(z = `Predicted Growth`), size = 2, alpha = 0.5) + 
  geom_contour(color = 'white', aes(z = `Predicted Growth`), size = 1, alpha = 0.8) + 
  geom_text_contour(aes(z = `Predicted Growth`), stroke = 0.1) + 
  guides(fill = pred) +
  new_scale_fill() +
  geom_point(data = hyb.preds, aes(x = LD1, y = LD2, fill = Growth, size = Growth), 
             pch = 21, alpha = 0.75) +
  scale_fill_continuous(low = 'black', high = 'white', 
                        limits = c(0, 1.5), breaks = c(0,0.5,1,1.5)) +
  scale_size_continuous(range = c(1,3.25), 
                        limits = c(0, 1.5), breaks = c(0,0.5,1,1.5)) +
  guides(fill = obs, size = obs) +
  stat_ellipse(data = parent.preds, level = c(0.95), 
               aes(x = LD1, y = LD2, color = Species), 
               geom = 'polygon', size = 1, type = 'norm',
               alpha = 0.2, lty = 3) +
  stat_ellipse(data = parent.preds, level = c(0.5), 
               aes(x = LD1, y = LD2, color = Species), 
               geom = 'polygon', size = 1, type = 'norm',
               alpha = 0.2, lty = 1) +
  geom_point(data = parent.preds, aes(x = LD1, y = LD2),
             size = 3.5) +
  geom_point(data = parent.preds, aes(x = LD1, y = LD2, color = Species),
             size = 2.5) +
  scale_shape_manual(values = 21) +
  scale_color_manual(values = parent.cols) +
  guides(color = spp.cols) +
  theme_bw(base_size = 18) 
best.geno.growth.wDead.mod.p

# "Quick" explore uncertainty of this model
# Fit models to 1000 bootstrap replicates of the data
bs.res <- c()
bs.summ<- c()
ld1.res <- c()
ld2.res <- c()
fit <- c()
ld1.summ <- c()
ld2.summ <- c()
ld1.fit.summ <- c()
ld2.fit.summ <- c()

for(i in 1:5000){
  print(i)
  
  # Resample the data with replacement
  boot <- hyb.preds[sample.int(nrow(hyb.preds), replace = TRUE),]

  # Fit the best-fit model from the observed data
  model <- gam(Growth ~ s(LD1, LD2) + Experiment + Lake + s(Site1, bs = 'fs') + s(Site2, bs = 'fs') + s(Site6, bs = 'fs') + 
                 s(Site7, bs = 'fs') + s(Site8, bs = 'fs') + s(Site9, bs = 'fs') + s(Site10, bs = 'fs'), 
               data = boot, method = 'GCV.Cp')
  
  # Output predictions at each point that we'll want to plot later
  res <- get.gam.viz(model, n.grid = 50, plot = F)
  colnames(res)[3] <- 'Fitness'
  
  # Calculate the mean predicted fitness at each point on the landscape,
  # either for LD1 or LD2
  ld1.mean <- 
    res %>% group_by(LD1) %>% 
    summarise(Mean.Fitness = mean(Fitness, na.rm = TRUE))
  ld2.mean <- 
    res %>% group_by(LD2) %>% 
    summarise(Mean.Fitness = mean(Fitness, na.rm = TRUE))
  
  # Now, pull out the values of each LD used in the resampled dataset,
  # and the predicted values of fitness either at each grid point on the landscape,
  # or the mean fitness calculated for each unique value of each LD
  bs.res <- append(bs.res, rep(i, 2500))
  bs.summ <- append(bs.summ, rep(i, 50))
  
  ld1.res <- append(ld1.res, res$LD1)
  ld2.res <- append(ld2.res, res$LD2)
  fit <- append(fit, res$Fitness)

  ld1.summ <- append(ld1.summ, ld1.mean$LD1)
  ld2.summ <- append(ld2.summ, ld2.mean$LD2)
  ld1.fit.summ <- append(ld1.fit.summ, ld1.mean$Mean.Fitness)
  ld2.fit.summ <- append(ld2.fit.summ, ld2.mean$Mean.Fitness)
}

# Put this all into a dataframe, and save
bs.comp.geno.mod <- 
  data.frame(Bootstrap = bs.res,
             LD1 = ld1.res, 
             LD2 = ld2.res,
             Fitness = fit)
bs.comp.geno.mod.summ <- 
  data.frame(Bootstrap = bs.summ, 
             LD1 = ld1.summ,
             LD2 = ld2.summ,
             LD1.MeanFit = ld1.fit.summ, 
             LD2.MeanFit = ld2.fit.summ)
write_delim(x = bs.comp.geno.mod.summ, 
            file = 'BestFit-Composite-wGeno-GAM-Bootstraps-Summary.txt',
            delim = '\t', col_names = T)
write_delim(x = bs.comp.geno.mod, 
            file = 'BestFit-Composite-wGeno-GAM-Bootstraps.txt',
            delim = '\t', col_names = T)


# First, get the predicted values from the observed dataset
# Because we have values of fitness at each of these points, 
# we want to predict fitness at
obs <- get.gam.viz(m13, n.grid = 50, plot = F)
colnames(obs)[3] <- 'Fitness'
lds.toPred <- 
  data.frame(LD1 = unique(obs$LD1),
             LD2 = unique(obs$LD2))

# Lets take the landscape in 3 slices - plot those three slices for 
# both LD1 and LD2 - should give a better sense of how well the two 
# match one another spatially.
# Here we're using the unsummarized data so we can look at predicted 
# fitness for each LD, in slices across the landscape

# First, using the values of each LD axis for which we have predicted 
# fitness for our observed data, determine break points
ld1.breaks <- 
  seq(from = min(obs$LD1), 
      to = max(obs$LD1), 
      by = dist(range(obs$LD1))/3)
ld1.breaks <- 
  list(c(ld1.breaks[1], ld1.breaks[2]), 
       c(ld1.breaks[2], ld1.breaks[3]),
       c(ld1.breaks[3], ld1.breaks[4]))
ld2.breaks <- 
  seq(from = min(obs$LD2), 
      to = max(obs$LD2), 
      by = dist(range(obs$LD2))/3)
ld2.breaks <- 
  list(c(ld2.breaks[1], ld2.breaks[2]), 
       c(ld2.breaks[2], ld2.breaks[3]),
       c(ld2.breaks[3], ld2.breaks[4]))

# 1) pull out a slice of the landscape. Let's say we want the 
# relationship between LD2 and fitness, but only at low lowest
# values of LD1 - we thus take a "slice" from LD1. 
# 2) for each iteration calculate the average fitness at each 
# unique LD value (in our example, LD2)
# 3) for each BS iteration, fit a GAM between  LD2 and mean Fitness
# 4) predict fitness using this model from a fixed set of positions
# 5) calculate the mean at each position +- 1 SD across bootstrap 
# replicates

bs.slice.preds <- 
  data.frame(Bootstrap = rep(1:5000, each = 150),
             Slice = rep(c(rep(1, 50), rep(2, 50),
                           rep(3, 50)), 5000),
             LD1 = NA,
             LD2 = NA, 
             LD1.PredFit = NA,
             LD2.PredFit = NA)
obs.slice.preds <- 
  data.frame(Slice = c(rep(1, 50), rep(2, 50),
                           rep(3, 50)),
             LD1 = NA,
             LD2 = NA, 
             LD1.PredFit = NA,
             LD2.PredFit = NA)
slice.summs <- list()

ld1.plts <- list()
ld2.plts <- list()
ld1.plts.ci <- list()
ld2.plts.ci <- list()

for(i in 1:3){
  print(paste0("Working on Slice Number ", i))
  #Pull out the 'slice'
  
  # First for the observed data
  ld1.slice.obs <- 
    obs[which(obs$LD2 >= 
                ld2.breaks[[i]][1] &
                obs$LD2 <= 
                ld2.breaks[[i]][2]),]
  ld2.slice.obs <- 
    obs[which(obs$LD1 >= 
                ld1.breaks[[i]][1] &
                obs$LD1 <= 
                ld1.breaks[[i]][2]),]
  
  # Define what you'll be predicting values for
  lds.toPred <-
    data.frame(LD1 = seq(from = min(ld1.slice.obs$LD1),
                         to = max(unique(ld1.slice.obs$LD1)),
                         length.out = 50),
               LD2 = seq(from = min(ld2.slice.obs$LD2),
                         to = max(unique(ld2.slice.obs$LD2)),
                         length.out = 50))
  
  ld1.mean.obs <- 
    ld1.slice.obs %>% group_by(LD1) %>% 
    summarise(Mean.Fitness = mean(Fitness, na.rm = TRUE))
  
  ld2.mean.obs <- 
    ld2.slice.obs %>% group_by(LD2) %>% 
    summarise(Mean.Fitness = mean(Fitness, na.rm = TRUE))
  
  # Now get predicted values at the focal points along each LD that we can thus
  # compare to the bootstrapped replicates
  # We will look at both means and medians - means may be overly sensitive to outliers
  # in the bootstrap replicates
  ld1.fit <- gam(Mean.Fitness ~ s(LD1), data = ld1.mean.obs)
  ld1.pred.obs <- predict(object = ld1.fit, newdata = lds.toPred)
  ld2.fit <- gam(Mean.Fitness ~ s(LD2), data = ld2.mean.obs)
  ld2.pred.obs <- predict(object = ld2.fit, newdata = lds.toPred)
  
  obs.slice.preds[which(obs.slice.preds == 1),2:3] <- lds.toPred
  obs.slice.preds[which(obs.slice.preds == 1),4] <- ld1.pred.obs
  obs.slice.preds[which(obs.slice.preds == 1),5] <- ld2.pred.obs
  
  # Then the bootstrapped replicates
  ld1.slice <- 
    bs.comp.geno.mod[which(bs.comp.geno.mod$LD2 >= 
                             ld2.breaks[[i]][1] &
                             bs.comp.geno.mod$LD2 <= 
                             ld2.breaks[[i]][2]),]
  ld2.slice <- 
    bs.comp.geno.mod[which(bs.comp.geno.mod$LD1 >= 
                             ld1.breaks[[i]][1] &
                             bs.comp.geno.mod$LD1 <= 
                             ld1.breaks[[i]][2]),]

  for(x in 1:100){
    if(x %% 100 == 0){print(x)}

    ld1.tmp.dat <- ld1.slice[which(ld1.slice$Bootstrap == x),]
    if(x %in% ld1.tmp.dat$Bootstrap){
      ld1.mean <-
        ld1.tmp.dat %>% group_by(LD1) %>%
        summarise(Mean.Fitness = mean(Fitness, na.rm = TRUE))
      ld1.fit <- gam(Mean.Fitness ~ s(LD1), data = ld1.mean)
      ld1.pred <- predict(object = ld1.fit, newdata = lds.toPred)
      bs.slice.preds[which(bs.slice.preds$Bootstrap == x &
                             bs.slice.preds$Slice == i),3] <- lds.toPred$LD1
      bs.slice.preds[which(bs.slice.preds$Bootstrap == x &
                             bs.slice.preds$Slice == i),5] <- ld1.pred
    }else{
      bs.slice.preds[which(bs.slice.preds$Bootstrap == x &
                             bs.slice.preds$Slice == i),c(3,5)] <- NA
    }
    ld2.tmp.dat <- ld2.slice[which(ld2.slice$Bootstrap == x),]
    if(x %in% ld2.tmp.dat$Bootstrap){
      ld2.mean <-
        ld2.tmp.dat %>% group_by(LD2) %>%
        summarise(Mean.Fitness = mean(Fitness, na.rm = TRUE))
      ld2.fit <- gam(Mean.Fitness ~ s(LD2), data = ld2.mean)
      ld2.pred <- predict(object = ld2.fit, newdata = lds.toPred)
      bs.slice.preds[which(bs.slice.preds$Bootstrap == x &
                             bs.slice.preds$Slice == i),4] <- lds.toPred$LD2
      bs.slice.preds[which(bs.slice.preds$Bootstrap == x &
                             bs.slice.preds$Slice == i),6] <- ld2.pred
    }else{
      bs.slice.preds[which(bs.slice.preds$Bootstrap == x &
                             bs.slice.preds$Slice == i),c(4,6)] <- NA
    }
  }
  # Calculate the mean value of predicted fitness at each position
  ld1.mean <-
    na.omit(bs.slice.preds[which(bs.slice.preds$Slice == i),]) %>%
    group_by(LD1) %>%
    summarise(PredFitness = mean(LD1.PredFit, na.rm = TRUE),
              SD = sd(LD1.PredFit, na.rm = TRUE),
              Upper = quantile(LD1.PredFit, probs = c(0.75), na.rm = T),
              Lower = quantile(LD1.PredFit, probs = c(0.25), na.rm = T))
  ld2.mean <-
    na.omit(bs.slice.preds[which(bs.slice.preds$Slice == i),]) %>%
    group_by(LD2) %>%
    summarise(PredFitness = mean(LD2.PredFit, na.rm = TRUE),
              SD = sd(LD2.PredFit, na.rm = TRUE),
              Upper = quantile(LD2.PredFit, probs = c(0.75), na.rm = T),
              Lower = quantile(LD2.PredFit, probs = c(0.25), na.rm = T))
  summ <-
    data.frame(LD1 = ld1.mean$LD1,
               LD2 = ld2.mean$LD2,
               LD1.PredFit = ld1.mean$PredFitness,
               LD2.PredFit = ld2.mean$PredFitness,
               LD1.PredFit.SD = ld1.mean$SD,
               LD2.PredFit.SD = ld2.mean$SD,
               LD1.PredFit.Upper = ld1.mean$Upper,
               LD2.PredFit.Upper = ld2.mean$Upper,
               LD1.PredFit.Lower = ld1.mean$Lower,
               LD2.PredFit.Lower = ld2.mean$Lower)
  
  slice.summs[[i]] <- summ

  # Plot using means +- 1 SD
  ld1.plts[[i]] <- 
    ggplot(dat = slice.summs[[i]], 
           aes(x = LD1, 
               y = LD1.PredFit)) +
    geom_line(aes(color = "Bootstrap")) + 
    geom_line(aes(y = LD1.PredFit-LD1.PredFit.SD, color = "Bootstrap"), lty = 2) + 
    geom_line(aes(y = LD1.PredFit+LD1.PredFit.SD, color = "Bootstrap"), lty = 2) + 
    geom_line(data = obs.slice.preds, aes(y = LD1.PredFit, color = "Observed")) +
    scale_color_manual(name = "Set", 
                       values = c("Bootstrap" = "black", 
                                  "Observed" = "red")) +
    new_scale("color") + 
    geom_point(data = parent.lds, 
               aes(x = LD1, y = 0, 
                   color = Species),
               shape = 124, alpha = 0.5, 
               size = 8) +
    scale_color_manual(values = parent.cols[c(2,3,1)]) +
    theme_bw(base_size = 12) + 
    ggtitle(label = paste0("Calculated for LD2 = ",
                           round(ld2.breaks[[i]][1], 3), 
                           " - ", round(ld2.breaks[[i]][2], 3))) +
    theme(plot.title = element_text(size=14)) + 
    ylab('Composite Fitness')
  
  ld2.plts[[i]] <- 
    ggplot(dat = slice.summs[[i]], 
           aes(x = LD2, 
               y = LD2.PredFit)) +
    geom_line(aes(color = "Bootstrap")) + 
    geom_line(aes(y = LD2.PredFit-LD2.PredFit.SD, color = "Bootstrap"), lty = 2) + 
    geom_line(aes(y = LD2.PredFit+LD2.PredFit.SD, color = "Bootstrap"), lty = 2) + 
    geom_line(data = obs.slice.preds, aes(y = LD2.PredFit, color = "Observed")) +
    scale_color_manual(name = "Set", 
                       values = c("Bootstrap" = "black", 
                                  "Observed" = "red")) +
    new_scale("color") + 
    geom_point(data = parent.lds, 
               aes(x = LD2, y = 0, 
                   color = Species),
               shape = 124, alpha = 0.5, 
               size = 8) +
    scale_color_manual(values = parent.cols[c(2,3,1)]) +
    theme_bw(base_size = 14) +
    ggtitle(label = paste0("Calculated for LD1 = ",
                           round(ld1.breaks[[i]][1], 3), 
                           " - ", round(ld1.breaks[[i]][2], 3))) +
    theme(plot.title = element_text(size=12)) + 
    ylab('Composite Fitness')
  
  # Then using means and 95% CIs (upper and lower 2.5th percentiles)
  ld1.plts.ci[[i]] <- 
    ggplot(dat = slice.summs[[i]], 
           aes(x = LD1, 
               y = LD1.PredFit)) +
    geom_line(aes(color = "Bootstrap")) + 
    geom_line(aes(y = LD1.PredFit.Lower, color = "Bootstrap"), lty = 2) + 
    geom_line(aes(y = LD1.PredFit.Upper, color = "Bootstrap"), lty = 2) + 
    geom_line(data = obs.slice.preds, aes(y = LD1.PredFit, color = "Observed")) +
    scale_color_manual(name = "Set", 
                       values = c("Bootstrap" = "black", 
                                  "Observed" = "red")) +
    new_scale("color") + 
    geom_point(data = parent.lds, 
               aes(x = LD1, y = 0, 
                   color = Species),
               shape = 124, alpha = 0.5, 
               size = 8) +
    scale_color_manual(values = parent.cols[c(2,3,1)]) +
    theme_bw(base_size = 12) + 
    ggtitle(label = paste0("Calculated for LD2 = ",
                           round(ld2.breaks[[i]][1], 3), 
                           " - ", round(ld2.breaks[[i]][2], 3))) +
    theme(plot.title = element_text(size=14)) + 
    ylab('Composite Fitness')
  
  ld2.plts.ci[[i]] <- 
    ggplot(dat = slice.summs[[i]], 
           aes(x = LD2, 
               y = LD2.PredFit)) +
    geom_line(aes(color = "Bootstrap")) + 
    geom_line(aes(y = LD2.PredFit.Lower, color = "Bootstrap"), lty = 2) + 
    geom_line(aes(y = LD2.PredFit.Upper, color = "Bootstrap"), lty = 2) + 
    geom_line(data = obs.slice.preds, aes(y = LD2.PredFit, color = "Observed")) +
    scale_color_manual(name = "Set", 
                       values = c("Bootstrap" = "black", 
                                  "Observed" = "red")) +
    new_scale("color") + 
    geom_point(data = parent.lds, 
               aes(x = LD2, y = 0, 
                   color = Species),
               shape = 124, alpha = 0.5, 
               size = 8) +
    scale_color_manual(values = parent.cols[c(2,3,1)]) +
    theme_bw(base_size = 14) +
    ggtitle(label = paste0("Calculated for LD1 = ",
                           round(ld1.breaks[[i]][1], 3), 
                           " - ", round(ld1.breaks[[i]][2], 3))) +
    theme(plot.title = element_text(size=12)) + 
    ylab('Composite Fitness')
  
  
}

leg <- get_legend(ld1.plts[[1]] + theme(legend.position = 'bottom'))
slice.plts <- 
  plot_grid(ld1.plts[[1]] + theme(legend.position = 'none'), 
            ld1.plts[[2]] + theme(legend.position = 'none'), 
            ld1.plts[[3]] + theme(legend.position = 'none'),
            ld2.plts[[1]] + theme(legend.position = 'none'), 
            ld2.plts[[2]] + theme(legend.position = 'none'), 
            ld2.plts[[3]] + theme(legend.position = 'none'),
            labels = c('a.', 'b.', 'c.', 
                            'd.', 'e.', 'f.'), 
            label_size = 18, nrow = 2)
slice.plts <- 
  plot_grid(slice.plts, leg, nrow = 2, rel_heights = c(1, 0.1))

ggsave(slice.plts, filename = 'Bootstrap-Slices-BestMod-CompositeMorphoGeno.pdf',
       height = 8, width = 12)

    



hyb.preds <- hyb.preds[-which(rowSums(is.na(hyb.preds[,11:20])) >= 1),]
hyb.preds$haps <- unite(hyb.preds[,11:20], haps, sep = '')
hyb.preds$haps <- as.factor(hyb.preds$haps[,1])
# One haploytpe is particularly common - give this one a different shape
hyb.preds$shape <- '21'
hyb.preds$shape[which(hyb.preds$haps == '0000000101')] <- '22'
#hyb.preds <- hyb.preds[,-c(11:20)]

best.geno.growth.wDead.mod.haplotypes.p <- 
  ggplot(best.mod.res, aes(x = LD1, y = LD2)) +
  scale_fill_scico(palette = 'davos', direction = 1, na.value = 'transparent', limits = c(0, upper)) +
  geom_raster(aes(fill = `Predicted Growth`)) +
  geom_contour(color = 'black', aes(z = `Predicted Growth`), size = 2, alpha = 0.5) + 
  geom_contour(color = 'white', aes(z = `Predicted Growth`), size = 1, alpha = 0.8) + 
  geom_text_contour(aes(z = `Predicted Growth`), stroke = 0.1) + 
  guides(fill = pred) +
  new_scale_fill() +
  scale_fill_continuous(low = 'black', high = 'white', 
                        limits = c(0, 1.5), breaks = c(0,0.5,1,1.5)) +
  scale_size_continuous(range = c(1,3.25), 
                        limits = c(0, 1.5), breaks = c(0,0.5,1,1.5)) +
  guides(fill = obs, size = obs) +
  # stat_ellipse(data = parent.preds, level = c(0.95), 
  #              aes(x = LD1, y = LD2, color = Species), 
  #              geom = 'polygon', size = 1, type = 'norm',
  #              alpha = 0.2, lty = 3) +
  # stat_ellipse(data = parent.preds, level = c(0.5), 
  #              aes(x = LD1, y = LD2, color = Species), 
  #              geom = 'polygon', size = 1, type = 'norm',
  #              alpha = 0.2, lty = 1) +
  # geom_point(data = parent.preds, aes(x = LD1, y = LD2),
  #            size = 3.5) +
  # geom_point(data = parent.preds, aes(x = LD1, y = LD2, color = Species),
  #            size = 2.5) +
  scale_shape_manual(values = c(21, 22)) +
  scale_color_manual(values = parent.cols) +
  new_scale_fill() +
  geom_point(data = hyb.preds, aes(x = LD1, y = LD2, fill = haps, size = Growth*3, shape = shape)
             , alpha = 0.75) +
  scale_fill_manual(values = rainbow(67)) + 
  guides(color = spp.cols) +
  theme_bw(base_size = 18) +
  theme(legend.position = 'none')

best.geno.growth.wDead.mod.haplotypes.p


# ggsave(plot = best.geno.growth.wDead.mod.p, 'Adapt-Landscape-Geno-Growth-wDead-SigSitesOnly.pdf', width = 10, height = 8, useDingbats = F)
# ggsave(plot = best.geno.growth.wDead.mod.haplotypes.p, 'Adapt-Landscape-Geno-Growth-wDead-SigSitesOnly-Haplotypes.pdf', width = 10, height = 8, useDingbats = F)

# And plot 3d
#pdf('Adapt-Landscape-Geno-Growth-wDead-SigSitesOnly-3D.pdf', height = 10, width = 8)
viz.gam(x = best.mod, view = c('LD1', 'LD2'), theta = 45, n.grid = 30, nCol = 100, color = 'davos')
#dev.off()
pdf('Adapt-Landscape-Geno-Growth-wDead-SigSitesOnly-3D-SameColorScale.pdf', height = 10, width = 8)
viz.gam(x = best.mod, view = c('LD1', 'LD2'), theta = 45, n.grid = 30, nCol = 100, color = 'davos', col.limits = c(0, upper))
dev.off()
########################################################################################
# And if we include genotypes significantly associated with survival from GEMMA LMM?
#############################################################################################
# What about if we include admixture proportion in these models?
setwd('~/Dropbox/Research/Martin-Berkeley-Postdoc/Pupfish-Fitness-Genomics/Results/AdaptiveLandscapes/')
hyb.preds <- orig

dat <- read.table('../../Data/Pupfish-final-trait-measurements-growth-admix.txt', 
                  header = T, check.names = F)
hybs <- dat[which(dat$Experiment %in% c('RTM1', 'RTM2')),]

genos <- read.table('../../Results/GEMMA/Traits/SignificantSites/growth/growthOnly-BonfSig-sites-Hybs-012-matrix.txt',
                    sep = '\t', header = T)

sites <- colnames(genos)[-1]

prefix <- "Site"
suffix <- seq(1:4)
colnames(genos) <- c('ID', paste(prefix, suffix, sep =''))

hyb.preds <- merge(hyb.preds, genos, by = 'ID')
hyb.preds[,11:14] <- hyb.preds[,11:14] %>% mutate_if(is.character, as.factor)
for(i in 11:14){
  hyb.preds[,i] <-
    gsub('wild', 0, hyb.preds[,i]) %>%
    gsub('single', 1, .) %>%
    gsub('double', 2, .) %>%
    factor(., ordered = T)
}

m1 <- gam(survival ~ s(LD1, LD2) + Experiment + Lake, data = hyb.preds, method = 'GCV.Cp', family = 'binomial')
m2 <- gam(survival ~ s(LD1, LD2) + Experiment + Lake + s(Site1, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp', family = 'binomial')
m3 <- gam(survival ~ s(LD1, LD2) + Experiment + Lake + s(Site2, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp', family = 'binomial')
m4 <- gam(survival ~ s(LD1, LD2) + Experiment + Lake + s(Site3, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp', family = 'binomial')
m5 <- gam(survival ~ s(LD1, LD2) + Experiment + Lake + s(Site4, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp', family = 'binomial')
m6 <- gam(survival ~ s(LD1, LD2) + Experiment + Lake + s(Site1, bs = 'fs') + s(Site2, bs = 'fs') + s(Site3, bs = 'fs') + s(Site4, bs = 'fs'), 
          data = hyb.preds, method = 'GCV.Cp', family = 'binomial')
# And try alternatives
m7 <- gam(survival ~ s(LD1, LD2) + Experiment + Lake + s(Site1, bs = 'fs') + s(Site3, bs = 'fs') + s(Site4, bs = 'fs'), 
          data = hyb.preds, method = 'GCV.Cp', family = 'binomial')
m8 <- gam(survival ~ s(LD1, LD2) + Experiment + Lake + s(Site1, bs = 'fs') + s(Site2, bs = 'fs') + s(Site4, bs = 'fs'), 
          data = hyb.preds, method = 'GCV.Cp', family = 'binomial')
m9 <- gam(survival ~ s(LD1, LD2) + Experiment + Lake + s(Site1, bs = 'fs') + s(Site2, bs = 'fs') + s(Site3, bs = 'fs'), 
          data = hyb.preds, method = 'GCV.Cp', family = 'binomial')
m10 <- gam(survival ~ s(LD1, LD2) + Experiment + Lake + s(Site1, bs = 'fs') + s(Site2, bs = 'fs'),
           data = hyb.preds, method = 'GCV.Cp', family = 'binomial')
m11 <- gam(survival ~ s(LD1, LD2) + Experiment + Lake + s(Site1, bs = 'fs') + s(Site3, bs = 'fs'),
           data = hyb.preds, method = 'GCV.Cp', family = 'binomial')
m12 <- gam(survival ~ s(LD1, LD2) + Experiment + Lake + s(Site1, bs = 'fs') + s(Site4, bs = 'fs'),
           data = hyb.preds, method = 'GCV.Cp', family = 'binomial')
m13 <- gam(survival ~ s(LD1, LD2) + Experiment + Lake + s(Site2, bs = 'fs') + s(Site3, bs = 'fs'),
           data = hyb.preds, method = 'GCV.Cp', family = 'binomial')
m14 <- gam(survival ~ s(LD1, LD2) + Experiment + Lake + s(Site2, bs = 'fs') + s(Site4, bs = 'fs'),
           data = hyb.preds, method = 'GCV.Cp', family = 'binomial')
m15 <- gam(survival ~ s(LD1, LD2) + Experiment + Lake + s(Site3, bs = 'fs') + s(Site4, bs = 'fs'),
           data = hyb.preds, method = 'GCV.Cp', family = 'binomial')

# And summarize for reference later
aics <- c(AICc(m1), AICc(m2), AICc(m3), 
          AICc(m4), AICc(m5), AICc(m6),
          AICc(m7), AICc(m8), AICc(m9),
          AICc(m10), AICc(m11), AICc(m12),
          AICc(m13), AICc(m14), AICc(m15))


weights  <- akaike.weights(aics)

fits <- 
  list(m1, m2, m3, 
       m4, m5, m6,
       m7, m8, m9,
       m10, m11, m12,
       m13, m14, m15)

mods <- 
  c(as.character(m1$formula), as.character(m2$formula), as.character(m3$formula), 
    as.character(m4$formula), as.character(m5$formula), as.character(m6$formula),
    as.character(m7$formula), as.character(m8$formula), as.character(m9$formula),
    as.character(m10$formula), as.character(m11$formula), as.character(m12$formula),
    as.character(m13$formula), as.character(m14$formula), as.character(m15$formula))

mod.table <- 
  data.frame('Model' = mods,
             'AICc' = aics, 
             'deltaAICc' = weights$deltaAIC,
             'AkWeights' = weights$weights)

mod.table
# write.table(x = mod.table, file = 'Morpho-Geno-Survival-GAM-ModTable.tsv', 
#             col.names = T, quote = F, sep = '\t', row.names = F)



# Now get the best model and prepare to visualize
best.mod <- which(aics == min(aics))
best.mod
best.mod.res <- get.gam.viz(fits[[best.mod]], n.grid = 100, too.far = 0.075)
best.mod <- fits[[best.mod]]
# save(best.mod, file = 'BestMorphoGeno-Survival-GAM.Rsave')

viz <- getViz(best.mod)
# 
# pdf('Survival-IndSNP-Splines.pdf')
# print(plot(viz, select = 2), pages = 1)
# dev.off()
# 
# pdf('Survival-IndSNP-LD-Contours.pdf', height = 24, width = 12)
# par(mfrow = c(1, 2),
#     oma = c(5,4,0,0) + 0.1,
#     mar = c(0,0,1,1) + 0.1)
# 
# vis.gam(best.mod, view = c('LD1', 'Site2'), theta = 230)
# vis.gam(best.mod, view = c('LD2', 'Site2'), theta = 230)
# dev.off()




colnames(best.mod.res)[3] <- 'Survival Probability'
colnames(hyb.preds)[9] <- 'Survival'

best.geno.surv.mod.p <- 
  ggplot(best.mod.res, aes(x = LD1, y = LD2)) +
  scale_fill_scico(palette = 'davos', direction = 1, na.value = 'transparent')+#, limits = c(0.05, 0.4)) +
  geom_raster(aes(fill = `Survival Probability`)) +
  geom_contour(color = 'black', aes(z = `Survival Probability`), size = 2, alpha = 0.5) + 
  geom_contour(color = 'white', aes(z = `Survival Probability`), size = 1, alpha = 0.8) + 
  geom_text_contour(aes(z = `Survival Probability`), stroke = 0.1) + 
  new_scale_fill() +
  geom_point(data = hyb.preds, aes(x = LD1, y = LD2, fill = Survival, size = Survival), 
             pch = 21, alpha = 0.75) +
  scale_fill_manual(values = c('grey40', 'white')) +
  scale_size_manual(values = c(1,2.25)) +
  new_scale_fill() +
  stat_ellipse(data = parent.preds, level = c(0.95), 
               aes(x = LD1, y = LD2, color = Species), 
               geom = 'polygon', size = 1, type = 'norm',
               alpha = 0.2, lty = 3) +
  stat_ellipse(data = parent.preds, level = c(0.5), 
               aes(x = LD1, y = LD2, color = Species), 
               geom = 'polygon', size = 1, type = 'norm',
               alpha = 0.2, lty = 1) +
  geom_point(data = parent.preds, aes(x = LD1, y = LD2),
             size = 3.5) +
  geom_point(data = parent.preds, aes(x = LD1, y = LD2, color = Species),
             size = 2.5) +
  scale_shape_manual(values = 21) +
  scale_fill_manual(values = parent.cols) +
  scale_color_manual(values = parent.cols) +
  # coord_cartesian(xlim = c(-10, 10), 
  #                 ylim = c(-10, 6)) +
  theme_bw(base_size = 18)
best.geno.surv.mod.p

# ggsave(plot = best.geno.surv.mod.p, 'Adapt-Landscape-Geno-Survival.pdf', width = 10, height = 8, useDingbats = F)


########################
# Now, plot all together in a 2x2 grid
surv.leg <- get_legend(best.geno.surv.mod.p)
grow.leg <- get_legend(best.geno.growth.mod.p)
grow.wdead.leg <- get_legend(best.geno.growth.wDead.mod.p)

surv.lands <- 
  plot_grid(best.surv.mod.p + 
              theme(legend.position = 'none'),
            best.geno.surv.mod.p + 
              theme(legend.position = 'none'),
            labels = c('A', 'B'),
            label_size = 28,
            align = 'H')

grow.lands <- 
  plot_grid(best.growth.mod.p + 
              theme(legend.position = 'none'),
            best.geno.growth.mod.p + 
              theme(legend.position = 'none'),
            labels = c('C', 'D'),
            label_size = 28,
            align = 'H')


grow.wDead.lands <- 
  plot_grid(best.growth.withdead.mod.p + 
              theme(legend.position = 'none'),
            best.geno.growth.wDead.mod.p + 
              theme(legend.position = 'none'),
            labels = c('A', 'B'),
            label_size = 28,
            align = 'H')

pdf('AllLandscapes-PhenoGeno-SurvGrowth.pdf', height = 12, width = 18)
plot_grid(surv.lands, surv.leg,
          grow.lands, grow.leg,
          rel_heights = c(1, 1),
          rel_widths = c(1, 0.15),
          ncol = 2, nrow = 2,
          align = 'HV')
dev.off()

pdf('Growth-wDead-PhenoGeno-Landscapes.pdf', height = 6, width = 12)
plot_grid(grow.wDead.lands, 
          grow.wdead.leg,
          rel_heights = c(1, 1),
          rel_widths = c(1, 0.2),
          ncol = 2,
          align = 'HV')
dev.off()

grow.wDead.lands <- 
  plot_grid(best.growth.withdead.mod.p + 
              theme(legend.position = 'none'),
            best.geno.growth.wDead.mod.p + 
              theme(legend.position = 'none'),
            labels = c('E', 'F'),
            label_size = 28,
            align = 'H')

pdf('AllLandscapes-PhenoGeno-SurvGrowth-PlusCombo.pdf', height = 18, width = 18)
plot_grid(surv.lands, surv.leg,
          grow.lands, grow.leg,
          grow.wDead.lands, grow.wdead.leg,
          rel_heights = c(1, 1, 1),
          rel_widths = c(1, 0.15),
          ncol = 2, nrow = 3,
          align = 'HV')
dev.off()



best.lands <- 
  plot_grid(best.surv.mod.p,
            best.geno.growth.mod.p,
            labels = c('A', 'B'),
            label_size = 28, nrow = 1,
            align = 'H')

ggsave(plot = best.lands, filename = 'BestLandscapes-PhenoGeno-SurvGrowth-PlusCombo.pdf', 
       height = 9, width = 21.33, useDingbats = F)


worst.lands <- 
  plot_grid(best.growth.mod.p,
            best.growth.withdead.mod.p,
            labels = c('A', 'B'),
            label_size = 28, nrow = 1,
            align = 'H')

ggsave(plot = worst.lands, filename = 'WorstLandscapes-Pheno-Growth-Composite.pdf', 
       height = 9, width = 21.33, useDingbats = F)

########################################################################################
# What about using de novo sweeping snps?
#############################################################################################
setwd('~/Dropbox/Research/Martin-Berkeley-Postdoc/Pupfish-Fitness-Genomics/Results/AdaptiveLandscapes/')
hyb.preds <- orig

dat <- read.table('../../Data/Pupfish-final-trait-measurements-growth-admix.txt', 
                  header = T, check.names = F)
hybs <- dat[which(dat$Experiment %in% c('RTM1', 'RTM2')),]

genos <- read.table('../../Results/FixedSweeping/denovo-sites.012-mat.txt',
                    sep = ' ', header = T)
genos <- genos[which(genos$ID %in% hybs$ID),]
sites <- colnames(genos)[-1]

prefix <- "Site"
suffix <- seq(1:8)
colnames(genos) <- c('ID', paste(prefix, suffix, sep =''))

hyb.preds <- merge(hyb.preds, genos, by = 'ID')
hyb.preds[,11:18] <- hyb.preds[,11:18] %>% mutate_if(is.character, as.factor)
for(i in 11:18){
  hyb.preds[,i] <-
    factor(hyb.preds[,i], ordered = T)
}

m1 <- gam(Growth ~ s(LD1, LD2) + Experiment + Lake, data = hyb.preds, method = 'GCV.Cp')
m2 <- gam(Growth ~ s(LD1, LD2) + Experiment + Lake + s(Site1, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
m3 <- gam(Growth ~ s(LD1, LD2) + Experiment + Lake + s(Site2, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
m4 <- gam(Growth ~ s(LD1, LD2) + Experiment + Lake + s(Site3, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
m5 <- gam(Growth ~ s(LD1, LD2) + Experiment + Lake + s(Site4, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
m6 <- gam(Growth ~ s(LD1, LD2) + Experiment + Lake + s(Site5, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
m7 <- gam(Growth ~ s(LD1, LD2) + Experiment + Lake + s(Site6, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
m8 <- gam(Growth ~ s(LD1, LD2) + Experiment + Lake + s(Site7, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
m9 <- gam(Growth ~ s(LD1, LD2) + Experiment + Lake + s(Site8, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
m10 <- gam(Growth ~ s(LD1, LD2) + Experiment + Lake + s(Site1, bs = 'fs') + s(Site2, bs = 'fs') + 
             s(Site3, bs = 'fs') + s(Site4, bs = 'fs') + s(Site5, bs = 'fs') + s(Site6, bs = 'fs') + s(Site7, bs = 'fs') + 
             s(Site8, bs = 'fs'), data = hyb.preds, method = 'GCV.Cp')
m13 <- gam(Growth ~ s(LD1, LD2) + Experiment + Lake + s(Site1, bs = 'fs') + s(Site2, bs = 'fs') + s(Site6, bs = 'fs') + 
             s(Site7, bs = 'fs') + s(Site8, bs = 'fs'), 
           data = hyb.preds, method = 'GCV.Cp')


# And summarize for reference later
aics <- c(AICc(m1), AICc(m2), AICc(m3), 
          AICc(m4), AICc(m5), AICc(m6),
          AICc(m7), AICc(m8), AICc(m9),
          AICc(m10))


weights  <- akaike.weights(aics)

fits <- 
  list(m1, m2, m3, 
       m4, m5, m6,
       m7, m8, m9,
       m10)

mods <- 
  c(as.character(m1$formula), as.character(m2$formula), as.character(m3$formula), 
    as.character(m4$formula), as.character(m5$formula), as.character(m6$formula),
    as.character(m7$formula), as.character(m8$formula), as.character(m9$formula),
    as.character(m10$formula))

mod.table <- 
  data.frame('Model' = mods,
             'AICc' = aics, 
             'deltaAICc' = weights$deltaAIC,
             'AkWeights' = weights$weights)

mod.table
write.table(x = mod.table, file = 'Morpho-Geno-Survival-GAM-ModTable.tsv', 
            col.names = T, quote = F, sep = '\t', row.names = F)



# Now get the best model and prepare to visualize
best.mod <- which(aics == min(aics))
best.mod
best.mod.res <- get.gam.viz(fits[[best.mod]], n.grid = 100, too.far = 0.075)
best.mod <- fits[[best.mod]]
save(best.mod, file = 'BestMorphoGeno-Survival-GAM.Rsave')

colnames(best.mod.res)[3] <- 'Survival Probability'
colnames(hyb.preds)[9] <- 'Survival'

best.geno.surv.mod.p <- 
  ggplot(best.mod.res, aes(x = LD1, y = LD2)) +
  scale_fill_viridis_c(option = "B", na.value = 'transparent')+#, limits = c(0, 0.6)) +
  geom_raster(aes(fill = `Survival Probability`)) +
  geom_contour(color = 'black', aes(z = `Survival Probability`), size = 2, alpha = 0.5) + 
  geom_contour(color = 'white', aes(z = `Survival Probability`), size = 1, alpha = 0.8) + 
  geom_text_contour(aes(z = `Survival Probability`), stroke = 0.1) + 
  new_scale_fill() +
  geom_point(data = hyb.preds, aes(x = LD1, y = LD2, fill = Survival, size = Survival), 
             pch = 21, alpha = 0.75) +
  scale_fill_manual(values = c('grey40', 'white')) +
  scale_size_manual(values = c(1,2.25)) +
  new_scale_fill() +
  stat_ellipse(data = parent.preds, level = c(0.95), 
               aes(x = LD1, y = LD2, color = Species), 
               geom = 'polygon', size = 1, type = 'norm',
               alpha = 0.2, lty = 3) +
  stat_ellipse(data = parent.preds, level = c(0.5), 
               aes(x = LD1, y = LD2, color = Species), 
               geom = 'polygon', size = 1, type = 'norm',
               alpha = 0.2, lty = 1) +
  geom_point(data = parent.preds, aes(x = LD1, y = LD2),
             size = 3.5) +
  geom_point(data = parent.preds, aes(x = LD1, y = LD2, color = Species),
             size = 2.5) +
  scale_shape_manual(values = 21) +
  scale_fill_manual(values = parent.cols) +
  scale_color_manual(values = parent.cols) +
  # coord_cartesian(xlim = c(-10, 10), 
  #                 ylim = c(-10, 6)) +
  theme_bw(base_size = 18)
best.geno.surv.mod.p

ggsave(plot = best.geno.surv.mod.p, 'Adapt-Landscape-Geno-Survival.pdf', width = 10, height = 8)



























