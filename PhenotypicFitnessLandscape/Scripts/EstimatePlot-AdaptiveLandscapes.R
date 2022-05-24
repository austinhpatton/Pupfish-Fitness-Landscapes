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

write.table(hyb.preds, file = '../Hybrid-GamDat.tsv', quote = F, sep = '\t', row.names = F)
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
write.table(x = mod.table, file = 'Morpho-SurvivalGAM-ModTable.tsv',
            col.names = T, quote = F, sep = '\t', row.names = F)

# Now get the best model and prepare to visualize
best.mod <- which(aics == min(aics))
best.mod
best.mod.res <- viz.gam(fits[[best.mod]], n.grid = 100, too.far = 0.075, view = c("LD1", "LD2"), plot.type = "contour")
best.mod <- fits[[best.mod]]
save(best.mod, file = 'BestMorpho-SurvivalGAM.Rsave')

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
  theme_bw(base_size = 18)
best.surv.mod.p
ggsave(plot = best.surv.mod.p, 'Adapt-Landscape-Survival.pdf', width = 10, height = 8, useDingbats = F)


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
write.table(x = mod.table, file = 'Morpho-Growth-WithNonSurvs-GAM-ModTable.tsv',
            col.names = T, quote = F, sep = '\t', row.names = F)

# Now get the best model and prepare to visualize
best.mod <- which(aics == min(aics))
best.mod
best.mod.res <- viz.gam(fits[[best.mod]], n.grid = 100, too.far = 0.075, view = c("LD1", "LD2"), plot.type = "contour")
colnames(best.mod.res)[3] <- 'Value'
best.mod <- fits[[best.mod]]
save(best.mod, file = 'BestMorpho-Growth-WithNonSurvs-GAM.Rsave')

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

ggsave(plot = best.growth.withdead.mod.p, 'Adapt-Landscape-Growth-wDead.pdf', width = 10, height = 8, useDingbats = F)

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
write.table(x = mod.table, file = 'Morpho-Growth-GAM-ModTable.tsv',
            col.names = T, quote = F, sep = '\t', row.names = F)

# Now get the best model and prepare to visualize
best.mod <- which(aics == min(aics))
best.mod 
best.mod.res <- viz.gam(fits[[best.mod]], n.grid = 100, too.far = 0.075, view = c("LD1", "LD2"), plot.type = "contour")
colnames(best.mod.res)[3] <- 'Value'
best.mod <- fits[[best.mod]]
save(best.mod, file = 'BestMorpho-Growth-GAM.Rsave')

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

ggsave(plot = best.growth.mod.p, 'Adapt-Landscape-Growth.pdf', width = 10, height = 8, useDingbats = F)

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
write.table(x = mod.table, file = 'Morpho-Geno-Growth-GAM-ModTable.tsv',
col.names = T, quote = F, sep = '\t', row.names = F)


# Now get the best model and prepare to visualize
best.mod <- which(aics == min(aics))
best.mod
best.mod.res <- viz.gam(fits[[best.mod]], n.grid = 100, too.far = 0.075, view = c("LD1", "LD2"), plot.type = "contour")
best.mod <- fits[[best.mod]]
save(best.mod, file = 'BestMorphoGeno-Growth-GAM.Rsave')

viz <- getViz(best.mod)

pdf('Growth-IndSNP-Splines.pdf')
print(plot(viz, select = 2), pages = 1)
dev.off()

pdf('Growth-IndSNP-LD-Contours.pdf', height = 24, width = 12)
par(mfrow = c(1, 2),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
vis.gam(best.mod, view = c('LD1', 'Site3'), theta = 230)
vis.gam(best.mod, view = c('LD2', 'Site3'), theta = 230)
dev.off()

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

ggsave(plot = best.geno.growth.mod.p, 'Adapt-Landscape-Geno-Growth.pdf', width = 10, height = 8, useDingbats = F)

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

pdf('Composite-IndSNP-Splines.pdf')
print(plot(viz, select = 2:8), pages = 1)
dev.off()

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
write.table(x = mod.table, file = 'Morpho-Geno-Composite-GAM-ModTable.tsv',
            col.names = T, quote = F, sep = '\t', row.names = F)

# Now get the best model and prepare to visualize
best.mod <- which(aics == min(aics))
best.mod
best.mod.res <- viz.gam(fits[[best.mod]], n.grid = 100, view = c("LD1", "LD2"), plot.type = "contour")
best.mod <- fits[[best.mod]]
save(best.mod, file = 'BestMorphoGeno-Composite-GAM.Rsave')

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


ggsave(plot = best.geno.growth.wDead.mod.p, 'Adapt-Landscape-Geno-Growth-wDead-SigSitesOnly.pdf', width = 10, height = 8, useDingbats = F)
ggsave(plot = best.geno.growth.wDead.mod.haplotypes.p, 'Adapt-Landscape-Geno-Growth-wDead-SigSitesOnly-Haplotypes.pdf', width = 10, height = 8, useDingbats = F)

# And plot 3d
pdf('Adapt-Landscape-Geno-Growth-wDead-SigSitesOnly-3D.pdf', height = 10, width = 8)
viz.gam(x = best.mod, view = c('LD1', 'LD2'), theta = 45, n.grid = 30, nCol = 100, color = 'davos', plot.type = 'persp', plot = T)
dev.off()
pdf('Adapt-Landscape-Geno-Growth-wDead-SigSitesOnly-3D-SameColorScale.pdf', height = 10, width = 8)
viz.gam(x = best.mod, view = c('LD1', 'LD2'), theta = 45, n.grid = 30, nCol = 100, color = 'davos', plot.type = 'persp', col.limits = c(0, upper))
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
write.table(x = mod.table, file = 'Morpho-Geno-Survival-GAM-ModTable.tsv',
            col.names = T, quote = F, sep = '\t', row.names = F)



# Now get the best model and prepare to visualize
best.mod <- which(aics == min(aics))
best.mod
best.mod.res <- viz.gam(fits[[best.mod]], n.grid = 100, too.far = 0.075, view = c("LD1", "LD2"), plot.type = "contour")
best.mod <- fits[[best.mod]]
save(best.mod, file = 'BestMorphoGeno-Survival-GAM.Rsave')

viz <- getViz(best.mod)

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

ggsave(plot = best.geno.surv.mod.p, 'Adapt-Landscape-Geno-Survival.pdf', width = 10, height = 8, useDingbats = F)


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
best.mod.res <- viz.gam(fits[[best.mod]], n.grid = 100, too.far = 0.075, view = c("LD1", "LD2"), plot.type = "contour")
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


