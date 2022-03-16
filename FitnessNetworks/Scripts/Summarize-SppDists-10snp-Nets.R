# At a broad scale, sampling sites across both sets (either for molluscivore or scale-eater),
# how far are each species from each other on a genotypic network?

library(ggplot2)
library(tidyverse)
library(igraph)
library(maditr)
library(cowplot)
library(ggforce)
library(pracma)
library(gtools)
library(agricolae)
library(MASS)
library(oddsratio)
library(scales)
library(facetscales)
library(plyr)
library(magrittr)
library(data.table)
library(future.apply)
library(progressr)
library(readr)

setwd('~/Dropbox/Research/Martin-Berkeley-Postdoc/Pupfish-Fitness-Genomics/Results/FixedSweeping_AllRedo/')
source.dir <- '~/Dropbox/Research/Martin-Berkeley-Postdoc/Pupfish-Fitness-Genomics/Results/FixedSweeping/'
source('../../Scripts/NetworkSummaryFunctions.R')

# Read in genotypes in 012 notation, where 0 are 'wildtype' (same as reference), 1 is single mutant, 2 is double
# We need to restrict to only those samples that have no missing data so the sequences are of equal length
dnv.p <- read.table(paste0(source.dir, 'denovo-sites-P.012-mat.txt'))
sgv.p <- read.table(paste0(source.dir, 'sgv-sites-P.012-mat.txt'))
intro.p <- read.table(paste0(source.dir, 'intro-sites-P.012-mat.txt'))
sgv.m <- read.table(paste0(source.dir, 'sgv-sites-M.012-mat.txt'))
intro.m <- read.table(paste0(source.dir, 'intro-sites-M.012-mat.txt'))

all.sites <- cbind(dnv.p, sgv.p, intro.p, sgv.m, intro.m)
all.sites <- all.sites[,!duplicated(colnames(all.sites))]

meta <- read.table('../FixedSweeping/Hyb-Parent-Metadata.tsv', header = T)
all.sites <- merge(meta, all.sites, by = 'ID')

# Randomly sample five sites from the the P alleles and 5 from the M alleles and do this 100 times
p.sites <- c(colnames(dnv.p[-1]), colnames(sgv.p[-1]), colnames(intro.p[-1]))
m.sites <- c(colnames(sgv.m[-1]), colnames(intro.m[-1]))
# Will use to plot an example network, and determine mean 
# distance among species
meta <- read.table('../../Data/Pupfish-final-trait-measurements-growth-admix.txt',
                   header = T, check.names = F)
meta <- meta[,c(1,2,37)]
colnames(meta)[2] <- 'Survival'
meta$Survival<- 
  meta$Survival %>% 
  gsub(pattern = 'NonSurvivor', replacement = 0) %>% 
  gsub(pattern = 'Survivor', replacement = 1)
meta$Survival <- as.numeric(meta$Survival)
genos <- list()
sites <- list()
seqs <- list()
fitness <- list()
net <- list()

spp.dists <- 
  data.frame(
    'AtoA' = rep(NA,100),
    'MtoM' = rep(NA,100),
    'PtoP' = rep(NA,100),
    'AtoM' = rep(NA,100),
    'AtoP' = rep(NA,100),
    'MtoP' = rep(NA,100)
  )

# Okay, calculating distance for these large haplotypes is incredibly slow
# As a result, we're going to construct all possible haplotypes that we might 
# observe given out permutation samples (there will be a huge number of them!!) 
# - by itself, this will take a long time, but we'll save that table so we never 
# have to do it again, and then use these pre-calculated distances to speed up 
# these loops/iterations. We have to limit to what we might actually see in our
# permutations since there are over 1.6 billion possible pairs of haplotypes 
# 10 snps long. 
haps <- c()
for(i in 1:100){
  samp <- c()
  while(length(samp) < 10){
    p.samp <- sample(p.sites, size = 5, replace = F)
    m.samp <- sample(m.sites, size = 5, replace = F)
    samp <- unique(c(p.samp, m.samp))
  }
  samp.sites <- which(colnames(all.sites) %in% samp)
  samp.genos <- all.sites[,c(1:4, samp.sites)]
  
  # Remove samps with missing data.
  samp.genos <- samp.genos[-which(rowSums(is.na(samp.genos[,5:14])) >= 1),]
  samp.sites <- colnames(samp.genos)[-c(1:4)]
  
  # Convert the individual genotypes to sequences
  new.haps <- unite(samp.genos[,-c(1:4)], haplotype, sep = '')
  haps <- append(haps, unique(new.haps[,1]))
}
# As is, this is 79 million. Will still take quite some time, 
# but a little better than 1.6 billion.
possible.pairs <- expand.grid(haps, haps)
# and exclude redundant pairs
possible.pairs <- 
  possible.pairs[as.character(possible.pairs$Var1) > 
                   as.character(possible.pairs$Var2),]
# Convert to character in advance to limit the number of 
# functions being carried out in the loops while we're indexing 
# this huge dataframe
possible.pairs[,1] <- as.character(possible.pairs[,1])
possible.pairs[,2] <- as.character(possible.pairs[,2])

# Save this dataframe took a while to make, and you never know when your computer is going to 
# crap the bed. 
# write.table(possible.pairs, './Possible-10SNP-HapPairs.txt', 
#             col.names = T, row.names = F, sep = '\t', quote = F)
possible.pairs <- 
  as.data.frame(
    read_delim('./Possible-10SNP-HapPairs.txt', 
               delim = '\t', 
               col_types = 'cc', 
               progress = T)
  )
#possible.pairs <- setDT(possible.pairs)

all.pairs <-
  data.frame(from = rep(NA, length(possible.pairs)),
             to = rep(NA, length(possible.pairs)),
             dist = rep(NA, length(possible.pairs)))

# We're going to use future_apply to parallelize this to some extent and speed up with C-
# apply is already faster than looping, and parallelizing obviously is even better. 
# As a results, let's make a function that takes in a row from the 
# possible.pairs dataframe, convert the haplotype to se vector of numbers,
# calculates the number of mutational steps between them, and then return both haplotypes 
# and their distance. 
get_dist_fast <- 
  function(focal_pair){
    hap1 <- focal_pair[1]
    hap2 <- focal_pair[2]
    a <- as.numeric(strsplit(hap1, "")[[1]])
    b <- as.numeric(strsplit(hap2, "")[[1]])
    dist <- as.numeric(sum(abs(a-b)))
    # Below is what will spit out the progress bar (see below)
    # Because it's specifief in the "with_progress" function, we don't 
    # specify it as input to the get_dist_fast function
    p()
    # res <- 
    #   data.frame(from = as.factor(hap1), 
    #              to = as.factor(hap2), 
    #              dist = as.numeric(dist))
    return(dist)
  }

# Great, now let's blaze forward, applying this function to every row in the 
# 79338419 row dataframe of possible pairs. Maybe work on something else in the meantime...
# So we have any sort of sense of how quickly this is going and how long it's going to take, let's
# use the with_progress function from progressr to get a fancy lil progress bar, telling us
# how long of an eternity we're going to need to wait. 

# The memory footprint of applying this function to the entire dataframe is just too large, 
# so what we're going to do is to iteratively split the original dataset into ten chunks,
# calculate distances for each, and then take the combined vector and add to the original dataframe.
ints1 <- seq(from = 1, to = nrow(possible.pairs), length.out = 21)
ints2 <- (ints1+1)[-c(1,21)]
intervals <-  
  list(ints1[1:2], c(ints2[1],ints1[3]),
       c(ints2[2], ints1[4]), c(ints2[3],ints1[5]), 
       c(ints2[4], ints1[6]), c(ints2[5],ints1[7]),
       c(ints2[6], ints1[8]), c(ints2[7],ints1[9]), 
       c(ints2[8], ints1[10]), c(ints2[9],ints1[11]),
       c(ints2[10],ints1[12]), c(ints2[11],ints1[13]),
       c(ints2[12],ints1[14]), c(ints2[13],ints1[15]),
       c(ints2[14],ints1[16]), c(ints2[15],ints1[17]),
       c(ints2[16],ints1[18]), c(ints2[17],ints1[19]),
       c(ints2[18],ints1[20]), c(ints2[19],ints1[21]))
hap.dists <- c()
for(i in 1:20){
  print(paste0("Beginning chunk ", i))
  # Get the interval
  int <- intervals[[i]]
  
  # Pull out the chunk. 
  print("Pulling out chunk")
  pair.chunk  <- possible.pairs[int[1]:int[2],]
  
  print("Done! Now the slog - calculating distance between haplotypes.")
  with_progress({
    p <- progressor(nrow(pair.chunk))
    chunk.dist <- 
      future_apply(pair.chunk, 
                   MARGIN = 1, 
                   FUN = get_dist_fast)
    hap.dists <- append(hap.dists, chunk.dist)
  })
  remaining <- 10-i
  print(paste0("Done! Congratualations, you have ", remaining, 
               " chunks left to do!"))
}


for(i in 1:nrow(possible.pairs)){
  if(i %% 100000 == 0){
    print(paste0('Iteration ', i, ' of 79890787'))
    this.much.done <- i/nrow(possible.pairs)
    print(paste0("Congratulations, you're ", this.much.done, '% of the way there!'))
  }
  
}

# Now SAVE!!!! This took for ever, even after attempting to parallelize,
# and you should never ever do this again. 
write.table(possible.pairs, './Possible-10SNP-HapPairs.txt',
            col.names = T, row.names = F, sep = '\t', quote = F)

# For composite fitness
for(i in 1:100){
  print(paste0('Sampling network #', i))
  samp <- c()
  while(length(samp) < 10){
    p.samp <- sample(p.sites, size = 5, replace = F)
    m.samp <- sample(m.sites, size = 5, replace = F)
    samp <- unique(c(p.samp, m.samp))
  }
  samp.sites <- which(colnames(all.sites) %in% samp)
  genos[[i]] <- all.sites[,c(1:4, samp.sites)]
  
  # Remove samps with missing data.
  genos[[i]] <- genos[[i]][-which(rowSums(is.na(genos[[i]][,5:14])) >= 1),]
  sites[[i]] <- colnames(genos[[i]])[-c(1:4)]
  
  # Convert the individual genotypes to sequences
  seqs[[i]] <- unite(genos[[i]][,-c(1:4)], haplotype, sep = '')
  seqs[[i]] <- data.frame('ID' = genos[[i]][,1], 
                          'Species' = genos[[i]][,2],
                          'Lake' = genos[[i]][,3],
                          'haplotype' = seqs[[i]])
  
  # And combine these sequences with growth and survival data
  fitness[[i]] <- merge(meta, seqs[[i]], by = 'ID', all = T)
  fitness[[i]] <- fitness[[i]][-which(is.na(fitness[[i]]$haplotype)),]
  
  # Now we need to collapse duplicate strings, averaging growth and survival for each.
  fitness[[i]]$sites.factor <- as.factor(fitness[[i]]$haplotype)
  
  growth <- with(na.omit(fitness[[i]]), tapply(Growth, sites.factor, mean))
  surv <- with(na.omit(fitness[[i]]), tapply(Survival, sites.factor, mean))
  
  num.a <- length(which(fitness[[i]]$Species == 'Generalist'))
  num.m <- length(which(fitness[[i]]$Species == 'Molluscivore'))
  num.p <- length(which(fitness[[i]]$Species == 'Scale-eater'))
  num.h <- length(which(fitness[[i]]$Species == 'Hybrid'))
  
  fit <- data.frame(
    'seq' = as.character(names(growth)),
    'growth' = growth,
    'survival' = surv
  )
  rownames(fit) <- NULL
  
  # Summarize the proportion of times each haplotype is associated with each species or hybrids
  site.spp <- fitness[[i]] %>% dplyr::count(sites.factor, Species)
  site.spp <- dcast(site.spp, sites.factor ~ Species)
  site.spp[is.na(site.spp)] <- 0
  site.spp$Generalist.prop <- site.spp$Generalist / rowSums(site.spp[,2:5])
  site.spp$Molluscivore.prop <- site.spp$Molluscivore / rowSums(site.spp[,2:5])
  site.spp$ScaleEater.prop <- site.spp$`Scale-eater` / rowSums(site.spp[,2:5])
  site.spp$Hybrid.prop <- site.spp$Hybrid / rowSums(site.spp[,2:5])
  site.spp <- site.spp[,c(1,6:9)]
  colnames(site.spp)[1] <- 'seq'
  
  haps <- merge(fit, site.spp, by = 'seq')
  
  print(paste0('Estimating network #', i))
  net[[i]] <- prep.network(haps, dist = 5, calc.dists = T, fit.measure = 'composite')
  
  net[[i]]$nodes$size <- net[[i]]$nodes$fitness
  net[[i]]$nodes$size[is.na(net[[i]]$nodes$size)] <- 
    max(net[[i]]$nodes$size[!is.na(net[[i]]$nodes$size)])
  
  # Now, determine mean distances between parental nodes
  a.nodes <- net[[i]]$nodes$seq[which(net[[i]]$nodes$Generalist.prop > 0)]
  m.nodes <- net[[i]]$nodes$seq[which(net[[i]]$nodes$Molluscivore.prop > 0)]
  p.nodes <- net[[i]]$nodes$seq[which(net[[i]]$nodes$ScaleEater.prop > 0)]
  
  a.a <- expand.grid(a.nodes, a.nodes)
  m.m <- expand.grid(m.nodes, m.nodes)
  p.p <- expand.grid(p.nodes, p.nodes)
  a.m <- expand.grid(a.nodes, m.nodes)
  a.p <- expand.grid(a.nodes, p.nodes)
  m.p <- expand.grid(m.nodes, p.nodes)
  
  pairs <-
    rbind(a.a, m.m, p.p, 
          a.m, a.p, m.p)
  pairs$pair <- c(
    rep('AtoA', nrow(a.a)),
    rep('MtoM', nrow(m.m)),
    rep('PtoP', nrow(p.p)),
    rep('AtoM', nrow(a.m)),
    rep('AtoP', nrow(a.p)),
    rep('MtoP', nrow(m.p))
  )
  pairs <- pairs[-which(pairs$Var1 == pairs$Var2),]
  pairs$dist <- NA
  for(n in 1:nrow(pairs)){
    a <- as.numeric(strsplit(as.character(pairs[n,1]), "")[[1]])
    b <- as.numeric(strsplit(as.character(pairs[n,2]), "")[[1]])
    pair <- rbind(a,b)
    diffs <- c()
    for(col in 1:ncol(pair)){
      if(max(pair[,col]) == 0){
        diffs[col] <- 0
      }else{
        diffs[col] <- max(pair[,col])- min(pair[,col])
      }
      pairs[n,4] <- sum(diffs)
    }
  }
  
  print('Calculating means')
  spp.dists[i,] <-
    c(mean(pairs$dist[which(pairs$pair == 'AtoA')]), 
      mean(pairs$dist[which(pairs$pair == 'MtoM')]),
      mean(pairs$dist[which(pairs$pair == 'PtoP')]),
      mean(pairs$dist[which(pairs$pair == 'AtoM')]),
      mean(pairs$dist[which(pairs$pair == 'AtoP')]),
      mean(pairs$dist[which(pairs$pair == 'MtoP')]))
  print('Done')
}

spp.dists.melted <- melt(spp.dists)
spp.dists.melted$variable <- 
  spp.dists.melted$variable %>%
  gsub(pattern = 'AtoA', replacement = 'Generalist to Generalist') %>%
  gsub(pattern = 'MtoM', replacement = 'Molluscivore to Molluscivore') %>%
  gsub(pattern = 'PtoP', replacement = 'Scale-eater to Scale-eater') %>%
  gsub(pattern = 'AtoM', replacement = 'Generalist to Molluscivore') %>%
  gsub(pattern = 'AtoP', replacement = 'Generalist to Scale-eater') %>%
  gsub(pattern = 'MtoP', replacement = 'Molluscivore to Scale-eater')

spp.dists.melted$variable <- 
  factor(spp.dists.melted$variable, 
         levels = c('Generalist to Generalist', 
                    'Molluscivore to Molluscivore',
                    'Scale-eater to Scale-eater',
                    'Generalist to Molluscivore', 
                    'Generalist to Scale-eater',
                    'Molluscivore to Scale-eater'))

a.col <- '#D2AF81FF'
m.col <- '#709AE1FF'
p.col <- '#FD7446FF'
a.m.col <- '#A1A5B1'
a.p.col <- '#E89264'
m.p.col <- '#B78794'
colors <- c(a.col, m.col, p.col, 
            a.m.col, a.p.col,
            m.p.col)

# get compact letter display
library(agricolae)
anova <- aov(lm(value ~ variable, data = spp.dists.melted))
tukey <- HSD.test(anova, trt = 'variable', unbalanced = T)
tukey.p <- HSD.test(anova, trt = 'variable', unbalanced = T,
                    group = F)
tukey.p <- 
  data.frame('Comparison' = rownames(tukey.p$comparison),
             'P-value' = tukey.p$comparison$pvalue)
write.table(tukey.p, file = 'SppPair-NumSteps-Tukey-Pvals.txt',
            row.names = F, quote = F)

upper <- 
  ggplot(spp.dists.melted, aes(x = value, y = variable)) + 
  stat_boxplot(geom ='errorbar')

lower <- 
  ggplot_build(upper)$data[[1]][,1]

upper <- 
  ggplot_build(upper)$data[[1]][,5]

lets <-
  data.frame('SppPair' = c(rownames(tukey$groups)),
             'Signif' = c(as.character(tukey$groups$groups)))
lets$SppPair <- 
  factor(lets$SppPair, 
         levels = c('Generalist to Generalist', 
                    'Molluscivore to Molluscivore',
                    'Scale-eater to Scale-eater',
                    'Generalist to Molluscivore', 
                    'Generalist to Scale-eater',
                    'Molluscivore to Scale-eater'))
lets <- rbind(lets[which(lets$SppPair == 'Generalist to Generalist'),],
              lets[which(lets$SppPair == 'Molluscivore to Molluscivore'),],
              lets[which(lets$SppPair == 'Scale-eater to Scale-eater'),],
              lets[which(lets$SppPair == 'Generalist to Molluscivore'),],
              lets[which(lets$SppPair == 'Generalist to Scale-eater'),],
              lets[which(lets$SppPair == 'Molluscivore to Scale-eater'),])
lets$NumSteps <- upper


coords <- c(1, max(upper)*1.05)
# 
# NumPaths.p <- 
#   ggplot(meds, aes(x = Source, y = NumPaths, fill = Source)) + 
#   facet_grid(. ~ Trajectory, scales = 'free', space = 'free') +
#   scale_fill_manual(values = c(sgv.col, intro.col, dnv.col, 
#                                sgv.intro.col, sgv.dnv.col)) +
#   # geom_violin(alpha = 0.5) +
#   geom_dotplot(binaxis='y', stackdir='center', drop = T, width = 0, aes(color = Source),
#                dotsize=0.5, binwidth = binwidth, stackratio = stackratio, 
#                position = position_jitter(width = 0, height = jitter), binpositions = 'all',
#                alpha = 0.25) +
#   scale_color_manual(values = c(sgv.col, intro.col, dnv.col, 
#                                 sgv.intro.col, sgv.dnv.col)) +
#   stat_boxplot(geom ='errorbar', width = 0.3, size = 1) +
#   geom_boxplot(width = 0.075, outlier.colour = NA) +
#   # geom_point(position = position_jitterdodge(jitter.height = 0.1, jitter.width = 3), alpha = 0.05) + 
#   geom_text(data = lets, label = lets$Signif, aes(y = NumPaths, x = Source),
#             vjust = -0.5, size = 6) +
#   ylim(coords) + 
#   ylab('Number of Accessible Paths') + 
#   theme_classic(base_size = 18) + 
#   theme(legend.position = 'none', 
#         axis.title.x = element_blank()) +
#   scale_x_discrete(labels = function(x) str_wrap(x, width = 7))
# NumPaths.p

colnames(spp.dists.melted) <- c('SppPair', 'NumSteps')
NumSteps.p <- 
  ggplot(data = spp.dists.melted, aes(y = SppPair, x = NumSteps, fill = SppPair)) +
  geom_text(data = lets, label = lets$Signif, aes(x = NumSteps, y = SppPair),
            size = 6, hjust = -0.5, vjust = -0.5) +
  geom_violin(alpha = 0.75) + 
  stat_boxplot(geom = 'errorbar', width = 0.25) + 
  geom_boxplot(width = 0.2, outlier.colour = 'NA') + 
  theme_bw(base_size = 16) +
  xlab('Number of Mutational Steps') + 
  scale_fill_manual(values = colors) + 
  coord_cartesian(xlim = c(0, 18)) + 
  theme(legend.position = 'none', 
        axis.title.y = element_blank())

ggsave(NumSteps.p, filename = 'NumSteps-SeparatingSppPairs-10SNP-Net.pdf', height = 7, width = 9)

growth.g <- graph_from_data_frame(net[[1]]$pairs, 
                                  directed=F, 
                                  vertices=net[[1]]$nodes)
edge.width=1/(exp(E(g)$weight)/5)



# fitness <- fitness[-which(is.na(fitness$Lake)),]
#fitness <- fitness[-which(fitness$Species == 'Hybrid' & is.na(fitness$Growth)),]

# Also, add in the LD scores
ld.scores <- read.table('../AdaptiveLandscapes/Sequenced-HybridParents-LD-Scores.txt', header = T)
ld.scores <- ld.scores[which(ld.scores$Species == 'Hybrid'),]
fitness <- merge(fitness, ld.scores[,-4], by = 'ID', all = T)
fitness <- fitness[which(!is.na(fitness$sites)),]


# Now we need to collapse duplicate strings, averaging growth and survival for each.
fitness$sites.factor <- as.factor(fitness$sites)

growth <- with(na.omit(fitness), tapply(Growth, sites.factor, mean))
surv <- with(na.omit(fitness), tapply(Survival, sites.factor, mean))
ld1 <- with(na.omit(fitness), tapply(LD1, sites.factor, mean))
ld2 <- with(na.omit(fitness), tapply(LD2, sites.factor, mean))

# And combine
fit <- data.frame(
  'seq' = as.character(names(growth)),
  'growth' = growth,
  'survival' = surv,
  'LD1' = ld1,
  'LD2' = ld2
)
rownames(fit) <- NULL

# Summarize the proportion of times each haplotype is associated with each species or hybrids
site.spp <- fitness %>% count(sites.factor,Species)
site.spp <- dcast(site.spp, sites.factor ~ Species)
site.spp[is.na(site.spp)] <- 0
site.spp$Generalist.prop <- site.spp$Generalist / rowSums(site.spp[,2:5])
site.spp$Molluscivore.prop <- site.spp$Molluscivore / rowSums(site.spp[,2:5])
site.spp$ScaleEater.prop <- site.spp$ScaleEater / rowSums(site.spp[,2:5])
site.spp$Hybrid.prop <- site.spp$Hybrid / rowSums(site.spp[,2:5])
site.spp <- site.spp[,c(1,6:9)]
colnames(site.spp)[1] <- 'seq'

haps <- merge(fit, site.spp, by = 'seq')

net <- prep.network(haps, dist = 4)
net$growth$nodes$growth.size <- net$growth$nodes$growth
net$growth$nodes$growth.size[is.na(net$growth$nodes$growth.size)] <- 
  max(net$growth$nodes$growth.size[!is.na(net$growth$nodes$growth.size)])
net$survival$nodes$growth.size <- net$survival$nodes$growth
net$survival$nodes$growth.size[is.na(net$survival$nodes$growth.size)] <- 
  max(net$survival$nodes$growth.size[!is.na(net$survival$nodes$growth.size)])



surv.g <- graph_from_data_frame(net$survival$pairs, directed=F, vertices=net$survival$nodes)
growth.g <- graph_from_data_frame(net$growth$pairs, directed=F, vertices=net$growth$nodes)
plot(growth.g)
# Two nodes are off in space: 010122001 and 020122001 - remove
#g <- delete.vertices(g, c('010122001', '020122001'))
#net$pairs <- net$pairs[-which(net$pairs$from %in% c('010122001', '020122001')),] 
#net$nodes <- net$nodes[-which(net$nodes$seq %in% c('010122001', '020122001')),] 


pies <- list()
for(i in 1:nrow(net$growth$nodes)){
  pies[[i]] <- c(net$growth$nodes[i,6], net$growth$nodes[i,7], 
                 net$growth$nodes[i,8], net$growth$nodes[i,9])
}

a.col <- '#D2AF81FF'
m.col <- '#709AE1FF'
p.col <- '#FD7446FF'
h.col <- '#8A9197FF' 

colnames(net$growth$nodes)[6:9] <- c('Generalist', 'Molluscivore', 'Scale-eater', 'Hybrid')

set.seed(26)
#set.seed(23)

layout <- layout_with_dh(growth.g, weight.edge.lengths = net$survival$pairs$dist)

pdf('Growth-Raw-ParentHyb-FixedSNP-Genotypic-Network-Pies-Dist4.pdf', height = 8, width = 8)
if (interactive()) {
  plot(growth.g, vertex.shape="pie", vertex.pie=pies,
       vertex.pie.color=list(c(a.col, m.col, p.col, h.col)), 
       vertex.size=rescale(V(growth.g)$growth.size, to = c(0.5,3))*4,
       vertex.label=NA, layout = layout)
  legend("topleft", legend=as.factor(colnames(net$growth$nodes[,6:9])),
         col = c(a.col, m.col, p.col, h.col), 
         bty = "n", pch=20 , pt.cex = 3, cex = 1.5, 
         text.col='black' , horiz = FALSE)
}
dev.off()

