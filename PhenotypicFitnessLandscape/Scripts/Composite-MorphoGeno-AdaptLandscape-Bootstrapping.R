# We want to explore uncertainty in our most complex model for an 
# adaptive landscape that includes LD axes constructed from all 
# morphological variables, and fitness associated SNPS.


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

load(file = 'BestMorphoGeno-Composite-GAM.Rsave')

# "Quick" explore uncertainty of this model
# Fit models to 10000 bootstrap replicates of the data
bs.res <- c()
bs.summ<- c()
ld1.res <- c()
ld2.res <- c()
fit <- c()
ld1.summ <- c()
ld2.summ <- c()
ld1.fit.summ <- c()
ld2.fit.summ <- c()

for(i in 1:10000){
  if(i %% 100 == 0){
    done <- (i/10000)*100
    print(paste0(done, "% done: ", 10000-i, " iterations to go"))
  }
  
  # Resample the data with replacement
  boot <- hyb.preds[sample.int(nrow(hyb.preds), replace = TRUE),]
  
  # Fit the best-fit model from the observed data
  model <- gam(Growth ~ s(LD1, LD2) + Experiment + Lake + s(Site1, bs = 'fs') + s(Site2, bs = 'fs') + s(Site6, bs = 'fs') + 
                 s(Site7, bs = 'fs') + s(Site8, bs = 'fs') + s(Site9, bs = 'fs') + s(Site10, bs = 'fs'), 
               data = boot, method = 'GCV.Cp')
  
  # Output predictions at each point that we'll want to plot later
  res <- get.gam.viz(model, n.grid = 30, plot = F)
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
  bs.res <- append(bs.res, rep(i, 900))
  bs.summ <- append(bs.summ, rep(i, 30))
  
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
  data.frame(Bootstrap = rep(1:10000, each = 900),
             LD1 = ld1.res, 
             LD2 = ld2.res,
             Fitness = fit)
bs.comp.geno.mod.summ <- 
  data.frame(Bootstrap = rep(1:10000, each = 30), 
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
obs <- get.gam.viz(best.mod, n.grid = 50, plot = F)
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
  data.frame(Bootstrap = rep(1:10000, each = 90),
             Slice = rep(c(rep(1, 30), rep(2, 30),
                           rep(3, 30)), 10000),
             LD1 = NA,
             LD2 = NA, 
             LD1.PredFit = NA,
             LD2.PredFit = NA)
obs.slice.preds <- 
  data.frame(Slice = c(rep(1, 30), rep(2, 30),
                       rep(3, 30)),
             LD1 = NA,
             LD2 = NA, 
             LD1.PredFit = NA,
             LD2.PredFit = NA)
slice.summs <- list()

ld1.plts <- list()
ld2.plts <- list()
ld1.plts.ci <- list()
ld2.plts.ci <- list()

for(i in 2:3){
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
                         length.out = 30),
               LD2 = seq(from = min(ld2.slice.obs$LD2),
                         to = max(unique(ld2.slice.obs$LD2)),
                         length.out = 30))
  
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
  
  for(x in 1:10000){
    if(x %% 100 == 0){
      done <- (x/10000)*100
      print(paste0(done, "% done: ", 10000-x, " iterations to go"))
    }    
    
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
              SD = sd(LD1.PredFit, na.rm = TRUE))
  ld2.mean <-
    na.omit(bs.slice.preds[which(bs.slice.preds$Slice == i),]) %>%
    group_by(LD2) %>%
    summarise(PredFitness = mean(LD2.PredFit, na.rm = TRUE),
              SD = sd(LD2.PredFit, na.rm = TRUE))
  
  slice <- bs.slice.preds[which(bs.slice.preds$Slice == i),]
    
  summ <-
    data.frame(LD1 = ld1.mean$LD1,
               LD2 = ld2.mean$LD2,
               LD1.PredFit = ld1.mean$PredFitness,
               LD2.PredFit = ld2.mean$PredFitness,
               LD1.PredFit.SD = ld1.mean$SD,
               LD2.PredFit.SD = ld2.mean$SD)
  
  slice.summs[[i]] <- summ
  
  # Plot using means +- 1 SD
  ld1.plts[[i]] <- 
    ggplot(dat = slice.summs[[i]], 
           aes(x = LD1, 
               y = LD1.PredFit)) +
    geom_line(data = slice,
              aes(x = LD1, y = LD1.PredFit, group = Bootstrap),
              alpha = 0.005, color = 'black') +
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
    coord_cartesian(ylim = c(min(slice.summs[[i]]$LD1.PredFit - 
                                slice.summs[[i]]$LD1.PredFit.SD)*1.5,
                             max(slice.summs[[i]]$LD1.PredFit + 
                                slice.summs[[i]]$LD1.PredFit.SD)*1.5)) + 
    ylab('Composite Fitness')
  
  ld2.plts[[i]] <- 
    ggplot(dat = slice.summs[[i]], 
           aes(x = LD2, 
               y = LD2.PredFit)) +
    geom_line(data = slice, 
              aes(x = LD1, y = LD1.PredFit, group = Bootstrap),
              alpha = 0.005, color = 'black') + 
    geom_line(aes(color = "Bootstrap")) + 
    geom_line(aes(y = LD2.PredFit-LD2.PredFit.SD, color = "Bootstrap"), lty = 2) + 
    geom_line(aes(y = LD2.PredFit+LD2.PredFit.SD, color = "Bootstrap"), lty = 2) + 
    geom_line(data = obs.slice.preds, aes(y = LD2.PredFit, color = "Observed")) +
    scale_color_manual(name = "Set", 
                       values = c("Bootstrap" = "grey", 
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
    coord_cartesian(ylim = c(min(slice.summs[[i]]$LD2.PredFit - 
                                   slice.summs[[i]]$LD2.PredFit.SD),
                             max(slice.summs[[i]]$LD2.PredFit + 
                                   slice.summs[[i]]$LD2.PredFit.SD))) + 
    ylab('Composite Fitness')
  
  # # Then using means and 95% CIs (upper and lower 2.5th percentiles)
  # ld1.plts.ci[[i]] <- 
  #   ggplot(dat = slice.summs[[i]], 
  #          aes(x = LD1, 
  #              y = LD1.PredFit)) +
  #   geom_line(aes(color = "Bootstrap")) + 
  #   geom_line(aes(y = LD1.PredFit.Lower, color = "Bootstrap"), lty = 2) + 
  #   geom_line(aes(y = LD1.PredFit.Upper, color = "Bootstrap"), lty = 2) + 
  #   geom_line(data = obs.slice.preds, aes(y = LD1.PredFit, color = "Observed")) +
  #   scale_color_manual(name = "Set", 
  #                      values = c("Bootstrap" = "black", 
  #                                 "Observed" = "red")) +
  #   new_scale("color") + 
  #   geom_point(data = parent.lds, 
  #              aes(x = LD1, y = 0, 
  #                  color = Species),
  #              shape = 124, alpha = 0.5, 
  #              size = 8) +
  #   scale_color_manual(values = parent.cols[c(2,3,1)]) +
  #   theme_bw(base_size = 12) + 
  #   ggtitle(label = paste0("Calculated for LD2 = ",
  #                          round(ld2.breaks[[i]][1], 3), 
  #                          " - ", round(ld2.breaks[[i]][2], 3))) +
  #   theme(plot.title = element_text(size=14)) + 
  #   ylab('Composite Fitness')
  # 
  # ld2.plts.ci[[i]] <- 
  #   ggplot(dat = slice.summs[[i]], 
  #          aes(x = LD2, 
  #              y = LD2.PredFit)) +
  #   geom_line(aes(color = "Bootstrap")) + 
  #   geom_line(aes(y = LD2.PredFit.Lower, color = "Bootstrap"), lty = 2) + 
  #   geom_line(aes(y = LD2.PredFit.Upper, color = "Bootstrap"), lty = 2) + 
  #   geom_line(data = obs.slice.preds, aes(y = LD2.PredFit, color = "Observed")) +
  #   scale_color_manual(name = "Set", 
  #                      values = c("Bootstrap" = "black", 
  #                                 "Observed" = "red")) +
  #   new_scale("color") + 
  #   geom_point(data = parent.lds, 
  #              aes(x = LD2, y = 0, 
  #                  color = Species),
  #              shape = 124, alpha = 0.5, 
  #              size = 8) +
  #   scale_color_manual(values = parent.cols[c(2,3,1)]) +
  #   theme_bw(base_size = 14) +
  #   ggtitle(label = paste0("Calculated for LD1 = ",
  #                          round(ld1.breaks[[i]][1], 3), 
  #                          " - ", round(ld1.breaks[[i]][2], 3))) +
  #   theme(plot.title = element_text(size=12)) + 
  #   ylab('Composite Fitness')
  
}

slice.1 <- bs.slice.preds[which(bs.slice.preds$Slice == 1),]
slice.2 <- bs.slice.preds[which(bs.slice.preds$Slice == 2),]
slice.3 <- bs.slice.preds[which(bs.slice.preds$Slice == 3),]

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

ld2.plts[[3]] + 
  geom_line(data = slice.3, 
            aes(x = LD2, y = LD2.PredFit, group = Bootstrap),
            alpha = 0.005, color = 'black')

























bs.preds <- 
  data.frame(Bootstrap = bs.summ,
             LD1 = rep(unique(obs$LD1), 1000),
             LD2 = rep(unique(obs$LD2), 1000), 
             LD1.PredFit = NA,
             LD2.PredFit = NA)
colnames(bs.preds)[4:5] <- c('LD1.PredFit', 'LD2.PredFit')
for(i in 1:1000){
  print(i)
  tmp.dat <- bs.comp.geno.mod.summ[which(bs.comp.geno.mod.summ$Bootstrap == i),]
  ld1.fit <- gam(LD1.MeanFit ~ s(LD1), data = tmp.dat)
  ld2.fit <- gam(LD2.MeanFit ~ s(LD2), data = tmp.dat)
  ld1.pred <- predict(object = ld1.fit, newdata = lds.toPred)
  ld2.pred <- predict(object = ld2.fit, newdata = lds.toPred)
  
  bs.preds[which(bs.preds$Bootstrap == i),4] <- ld1.pred
  bs.preds[which(bs.preds$Bootstrap == i),5] <- ld2.pred
}

ld1.summ <- bs.preds %>% group_by(LD1) %>% 
  summarise(Pred.Fitness = mean(LD1.PredFit, na.rm = TRUE),
            SD = sd(LD1.PredFit, na.rm = TRUE)) %>%
  mutate(Name = "LD1.PredFit")

ld2.summ <- bs.preds %>% group_by(LD2) %>% 
  summarise(Pred.Fitness = mean(LD2.PredFit, na.rm = TRUE),
            SD = sd(LD2.PredFit, na.rm = TRUE)) %>%
  mutate(Name = "LD2.PredFit")

# pull out the observed
ld1.mean.obs <- 
  obs %>% group_by(LD1) %>% 
  summarise(Mean.Fitness = mean(Fitness, na.rm = TRUE),
            SD = sd(Fitness, na.rm = TRUE))
ld2.mean.obs <- 
  obs %>% group_by(LD2) %>% 
  summarise(Mean.Fitness = mean(Fitness, na.rm = TRUE),
            SD = sd(Fitness, na.rm = TRUE))

obs.summ <- 
  data.frame(LD1 = ld1.mean.obs$LD1,
             LD2 = ld2.mean.obs$LD2,
             LD1.MeanFit = ld1.mean.obs$Mean.Fitness,
             LD2.MeanFit = ld2.mean.obs$Mean.Fitness,
             LD1.SD = ld1.mean.obs$SD,
             LD2.SD = ld2.mean.obs$SD)

# We also want to show the distributions of parentals in these 
# plots
parent.lds <- as.data.frame(ld.pure$x)
parent.lds$Species <- parents$Species

ld2.plot <- 
  ggplot(dat = ld2.summ, 
         aes(x = LD2, 
             y = Pred.Fitness)) +
  geom_line() + 
  geom_line(aes(y = Pred.Fitness-SD), lty = 2) + 
  geom_line(aes(y = Pred.Fitness+SD), lty = 2) + 
  geom_line(data = obs.summ, aes(y = LD2.MeanFit),
            color = 'red') +
  geom_point(data = parent.lds, 
             aes(x = LD2, y = 0, 
                 color = Species),
             shape = 124, alpha = 0.5, 
             size = 8) +
  scale_color_manual(values = parent.cols[c(2,3,1)]) +
  coord_cartesian(xlim = c(min(ld2.summ$LD2), 
                           max(ld2.summ$LD2))) +
  theme_bw() + 
  ylab('Composite Fitness')
ld2.plot

ld1.plot <- 
  ggplot(dat = ld1.summ, 
         aes(x = LD1, 
             y = Pred.Fitness)) +
  geom_line() + 
  geom_line(aes(y = Pred.Fitness-SD), lty = 2) + 
  geom_line(aes(y = Pred.Fitness+SD), lty = 2) + 
  geom_line(data = obs.summ, aes(y = LD1.MeanFit),
            color = 'red') +
  geom_point(data = parent.lds, 
             aes(x = LD1, y = 0, 
                 color = Species),
             shape = 124, alpha = 0.5, 
             size = 8) +
  scale_color_manual(values = parent.cols[c(2,3,1)]) +
  coord_cartesian(xlim = c(min(ld1.summ$LD1), 
                           max(ld1.summ$LD1))) +
  theme_bw() + 
  ylab('Composite Fitness')
ld1.plot
