require(RANN)
require(gdata)
require(igraph)
require(zoo)
require(tidyverse)
require(rlist)
require(scales)
library(robustbase)
library(R.utils)

# These are a series of functions to do the legwork for taking data in 012 format and combining into haplotypes, and preparing to plot as a network
# Also includes functions to sumarize fitness landscapes

summary_func <- function(x, col){
  c(med = median(x[[col]], na.rm=TRUE),
    mean = mean(x[[col]], na.rm=TRUE),
    sd = sd(x[[col]], na.rm=TRUE),
    se = sd(x[[col]], na.rm=TRUE)/sqrt(length(x[[col]])))
}

get.rel <- function(summs){
  summs$RelNumPaths <- summs$NumPaths / summs$NumNodes.InNet
  summs <- summs[,c(1:2,12,3:11)]
}

pairwiseOR <- function(dat, trajectory = NULL, variable, by.source = T){
  if(by.source == F){
    id <- which(colnames(dat) == variable)
    trajectory <- which(colnames(dat) == 'Trajectory')
    dat <- 
      dat[which(is.finite(dat[,id])),c(trajectory, id)]
    
    odds.ratios <- 
      data.frame(Comparison = 'AtoM-AtoP',
                 Accessibility = variable, OddsRatio = NA, 
                 OR.Lower = NA, OR.Upper = NA, LRT.P = NA)
    
    # Replace tajectory with a binary - AtoP = 0.
    # This translates to an odds ratio that describes the effect of the trajectory being 
    # AtoM on accessibility, relative to AtoP.
    dat$Trajectory <- 
      gsub(dat$Trajectory, pattern = "AtoP", replacement = 0)
    dat$Trajectory[-which(dat$Trajectory == "0")] <- 1
    dat$Trajectory <- as.numeric(dat$Trajectory)
    colnames(dat) <- c('Trajectory', 'Accessibility')
    
    
    # log10 transform Accessibility for the measures that are very skewed.
    # Without doing so, you can get some obscenely large(and nonsensical) odds ratios
    # We do this per species, since the analyses are effectively separate, 
    # using different sets of sites and thus different networks to measure accessibility
    if(variable %in% c('RelNumPaths', 'RelNumPathsToPeaks',
                       'MeanFit', 'MeanFit.ByLen', 'MedFit', 'MedFit.ByLen')){
      mod <- glm(Trajectory ~ scale(log10(Accessibility)), 
                 family = binomial(link="logit"), 
                 data = dat)
      odds.ratios[,3] <- exp(mod$coefficients[2])
      odds.ratios[,4:5] <- exp(confint(mod)[c(2,4)])
      odds.ratios[,6] <- anova(mod, test = "Chisq")$`Pr(>Chi)`[2]
    }else{
      mod <- glm(Trajectory ~ Accessibility, 
                 family = binomial(link="logit"), 
                 data = dat)
      odds.ratios[,3] <- exp(mod$coefficients[2])
      odds.ratios[,4:5] <- exp(confint(mod)[c(2,4)])
      odds.ratios[,6] <- anova(mod, test = "Chisq")$`Pr(>Chi)`[2]
      
    }
  }else{
    # Reduce the dataframe to what we are comparing
    id <- which(colnames(dat) == variable)
    source <- which(colnames(dat) == 'Source')
    dat <- 
      dat[which(dat$Trajectory == trajectory &
                  is.finite(dat[,id])),c(source, id)]
    
    if(trajectory == 'AtoM'){
      comparisons <- 
        list(c("SGV", "Introgression"), 
             c("SGV", "SGV + Introgression"))
      sources <- 
        c("SGV", "Introgression", 
          "SGV + Introgression")
    }else{
      comparisons <- 
        list(c("SGV", "Introgression"), 
             c("SGV", "De novo"),
             c("SGV", "SGV + Introgression"),
             c("SGV", "SGV + De novo"),
             c("SGV", "SGV + De novo + Introgression"))
      sources <- 
        c("SGV", "Introgression", "De novo", 
          "SGV + Introgression",
          "SGV + De novo",
          "SGV + De novo + Introgression")
    }
    
    odds.ratios <- 
      data.frame(Trajectory = trajectory, Source = sources, 
                 Accessibility = variable, OddsRatio = NA, 
                 OR.Lower = NA, OR.Upper = NA, LRT.P = NA)
    
    for(x in 1:length(comparisons)){
      tmp.dat <- dat[which(dat$Source %in% comparisons[[x]]),]
      colnames(tmp.dat) <- c('Source', 'Accessibility')
      
      # Replace source with a binary - SGV = 0 since we want to compare the effect of the 
      # other source on accessibility relative to SGV. 
      tmp.dat$Source <- 
        gsub(tmp.dat$Source, pattern = "SGV", replacement = 0)
      tmp.dat$Source[-which(tmp.dat$Source == "0")] <- 1
      tmp.dat$Source <- as.numeric(tmp.dat$Source)
      
      # log10 transform Accessibility for the measures that are very skewed.
      # Without doing so, you can get some obscenely large(and nonsensical) odds ratios
      # We do this per species, since the analyses are effectively separate, 
      # using different sets of sites and thus different networks to measure accessibility
      if(variable %in% c('RelNumPaths', 'RelNumPathsToPeaks',
                         'MeanFit', 'MeanFit.ByLen', 'MedFit', 'MedFit.ByLen')){
        mod <- glm(Source ~ scale(log10(Accessibility)), family = binomial(link="logit"), data = tmp.dat)
        odds.ratios[x+1, 4] <- exp(mod$coefficients[2])
        odds.ratios[x+1, 5:6] <- exp(confint(mod)[c(2,4)])
        odds.ratios[x+1,7] <- anova(mod, test = "Chisq")$`Pr(>Chi)`[2]
      }else{
        mod <- glm(Source ~ Accessibility, family = binomial(link="logit"), data = tmp.dat)
        odds.ratios[x+1, 4] <- exp(mod$coefficients[2])
        odds.ratios[x+1, 5:6] <- exp(confint(mod)[c(2,4)])
        odds.ratios[x+1,7] <- anova(mod, test = "Chisq")$`Pr(>Chi)`[2]
      }
    }
    # And correct p-values for multiple testing. 
    odds.ratios[,7] <- p.adjust(odds.ratios[,7], method = 'fdr')
  }
  return(odds.ratios)
}

# let's quick make a little function to simplify the plotting prep work
prep.plot <- function(summ.stats = NULL, variable = NULL, mode = c('param', 'nonparam'), multiplier = 1.05){
  id <- which(colnames(summ.stats) == variable)
  colnames(summ.stats)[id] <- 'Stat'
  if(mode == 'param'){
    res <- t.test(summ.stats$Stat ~ summ.stats$Trajectory)
  }else{
    res <- kruskal.test(summ.stats$Stat ~ summ.stats$Trajectory)
  }
  
  
  a2m.upper <- 
    ggplot(summ.stats[which(summ.stats$Trajectory == 'AtoM'),], aes(x = Trajectory, y = Stat)) + 
    stat_boxplot(geom ='errorbar')
  a2p.upper <- 
    ggplot(summ.stats[which(summ.stats$Trajectory == 'AtoP'),], aes(x = Trajectory, y = Stat)) + 
    stat_boxplot(geom ='errorbar')
  
  a2m.lower <- 
    ggplot_build(a2m.upper)$data[[1]][,1]
  a2p.lower <- 
    ggplot_build(a2p.upper)$data[[1]][,1]
  
  a2m.upper <- 
    ggplot_build(a2m.upper)$data[[1]][,5]
  a2p.upper <- 
    ggplot_build(a2p.upper)$data[[1]][,5]
  
  
  coords <- c(min(c(a2m.lower, a2p.lower)), max(c(a2m.upper, a2p.upper))*multiplier)
  
  p <- res$p.value
  if(p < 0.0001){
    p <- "P < 0.0001"
  }else{
    p <- format(round(p, 4), scientific = F)
    p <- paste0("P = ", p)
  }
  
  return(list(p = p, coords = coords))
  
}



prep.network <-
  function(dat, dist = 1, calc.dists = F, precalc.dists = NULL, fit.measure = c('survival', 'growth', 'composite')){
    # Construct a dataframe that defines all relationships between nodes
    # On other wise, we want to define the distance between each node, in terms of number of nucleotides

    # and now, calculate the difference (in nucleotides) between each pair of strings
    if(calc.dists == T){
      combs <- 
        combn(x = as.character(dat$seq), m = 2, simplify = F)
      
      pairs <- data.frame(from = rep(NA, length(combs)), 
                          to = rep(NA, length(combs)),
                          dist = rep(NA, length(combs)))
      
      for(i in 1:length(combs)){
        a <- as.numeric(strsplit(combs[[i]][1], "")[[1]])
        b <- as.numeric(strsplit(combs[[i]][2], "")[[1]])
        pair <- rbind(a,b)
        diffs <- c()
        for(col in 1:ncol(pair)){
          if(max(pair[,col]) == 0){
            diffs[col] <- 0
          }else{
            diffs[col] <- max(pair[,col])- min(pair[,col])
          }
          pairs[i,1] <- combs[[i]][1]
          pairs[i,2] <- combs[[i]][2]
          pairs[i,3] <- sum(diffs)
        }
      }
    }else{
      pairs <- precalc.dists
    }
    
    pairs <- pairs[which(pairs$dist <= dist),]
    if(nrow(pairs) > 0){
      pairs$fit1 <- NA
      pairs$fit2 <- NA
      # Fill in with fitness data
      if(fit.measure == 'composite'){
        for(s in 1:nrow(pairs)){
          pairs[s,4] <- dat[which(dat$seq == pairs[s,1]),2]
          pairs[s,5] <- dat[which(dat$seq == pairs[s,2]),2]
        }
      }else{
        if(fit.measure == 'growth'){
          for(s in 1:nrow(pairs)){
            pairs[s,4] <- dat[which(dat$seq == pairs[s,1]),2]
            pairs[s,5] <- dat[which(dat$seq == pairs[s,2]),2]
          }
          pairs <- pairs[-which(pairs$fit1 == 0 | pairs$fit2 == 0),]
        }else{
          for(s in 1:nrow(pairs)){
            pairs[s,4] <- dat[which(dat$seq == pairs[s,1]),3]
            pairs[s,5] <- dat[which(dat$seq == pairs[s,2]),3]
          }
        }
      }
      
      # Now begin getting what we need for a directed graph
      # For neighboring parent pairs, duplicate these with swapped From/To
      parent.pairs <- pairs[which(is.na(pairs$fit1) & is.na(pairs$fit2)),]
      parent.swap <- parent.pairs[,c(2:1,3,5:4)]
      colnames(parent.swap) <- colnames(pairs)
      parent.pairs <- rbind(parent.pairs, parent.swap)
      
      # Now temporarily remove these
      if(nrow(parent.pairs) > 0){
        pairs <- pairs[-which(is.na(pairs$fit1) & is.na(pairs$fit2)),]
      }
      
      # Now, we need to be sure that when there is a parent node, 
      # edges lead to it (assuming parental fitness)
      to.swap <- which(is.na(pairs$fit1) & !is.na(pairs$fit2))
      if(length(to.swap) > 0){
        replacement <- replace(pairs[to.swap,], values = pairs[to.swap,c(2:1,3,5:4)])
        pairs <- pairs[-to.swap,]
        pairs <- rbind(pairs, replacement)
      }
      
      
      # Duplicate the equal fitness pairs, swapping from/to, and temporarily remove
      eq.pairs <- pairs[which(pairs$fit1 == pairs$fit2),]
      eq.pairs <- eq.pairs[,c(2:1,3,5:4)]
      colnames(eq.pairs) <- colnames(pairs)
      if(nrow(eq.pairs) > 0){
        pairs <- pairs[-which(pairs$fit1 == pairs$fit2),]
      }
      
      
      # And fix the from/to so that edges always go to higher fitness
      to.swap <- which(pairs$fit2 < pairs$fit1)
      replacement <- replace(pairs[to.swap,], values = pairs[to.swap,c(2:1,3,5:4)])
      if(length(to.swap) > 0){
        pairs <- pairs[-to.swap,]
        pairs <- rbind(pairs, replacement)
      }
      
      # and combine
      pairs <- rbind(pairs, eq.pairs)
      pairs <- rbind(pairs, parent.pairs)
      
      # remove edges connecting non-survivors
      drop <- which(pairs$fit1 == 0 & pairs$fit2 == 0)
      if(length(drop > 0)){
        pairs <- pairs[-drop,]
      }
      
      nodes <- dat[which(dat$seq %in% c(pairs$from, pairs$to)),]
      if(fit.measure == 'survival'){
        colnames(nodes)[3] <- 'fitness'
      }else{
        colnames(nodes)[2] <- 'fitness'
      }
      
      fit.net <- list('pairs' = pairs, 
                      'nodes' = nodes)
    }else{
      fit.net <- NULL
    }
      
    return(fit.net)
  }

get.dists <-
  function(seq1, seq2){
    a <- as.numeric(strsplit(seq1, "")[[1]])
    b <- as.numeric(strsplit(seq2, "")[[1]])
    pair <- rbind(a,b)
    diffs <- c()
    for(col in 1:ncol(pair)){
      if(max(pair[,col]) == 0){
        diffs[col] <- 0
      }else{
        diffs[col] <- max(pair[,col])- min(pair[,col])
      }
    }
    distance <- sum(diffs)
    return(distance)
  }

get.all.pairs <- 
  function(dat){
    # Construct a dataframe that defines all relationships between nodes
    # On other wise, we want to define the distance between each node, in terms of number of nucleotides
    combs <- 
      combn(x = as.character(dat$seq), m = 2, simplify = F)
    
    pairs <- data.frame(from = rep(NA, length(combs)), 
                        to = rep(NA, length(combs)),
                        dist = rep(NA, length(combs)))
    
    # and now, calculate the difference (in nucleotides) between each pair of strings
    for(i in 1:length(combs)){
      a <- as.numeric(strsplit(combs[[i]][1], "")[[1]])
      b <- as.numeric(strsplit(combs[[i]][2], "")[[1]])
      pair <- rbind(a,b)
      diffs <- c()
      for(col in 1:ncol(pair)){
        if(max(pair[,col]) == 0){
          diffs[col] <- 0
        }else{
          diffs[col] <- max(pair[,col])- min(pair[,col])
        }
        pairs[i,1] <- combs[[i]][1]
        pairs[i,2] <- combs[[i]][2]
        pairs[i,3] <- sum(diffs)
      }
    }
    return(pairs)
  }



get.trajectories <- function(g, spp.pair, parent.haps, net, sampled.generalist, sampled.specialist){
  paths <- list()
  # we are keeping the flanking parental nodes for now
  # keep <- which(parent.haps %in% c(as.character(spp.pair[i,1]),
  #                                  as.character(spp.pair[i,2])))
  # drop.nodes <- as.character(parent.haps[-keep])
  # drop.nodes <- as.character(parent.haps[-which(parent.haps %in% c(sampled.generalist, sampled.specialist))])
  # tmp.g <- delete_vertices(g, drop.nodes)
  # 
  # Only continue on the following steps if the parental nodes are still connected after dropping
  # other parental nodes
  # if(edge.connectivity(tmp.g, as.character(spp.pair[i,1]), as.character(spp.pair[i,2])) > 0){
  if(edge.connectivity(g, sampled.generalist, sampled.specialist) > 0){
    # # also want to drop nodes that are too far away from either focal node. We need to restrict the genotype 
    # # space to something smalled so we can more tractably enumerate all simple paths between nodes. 
    # short.pathLen <- length(shortest_paths(tmp.g, from = sampled.generalist,to = sampled.specialist)$vpath[[1]])
    # dd <- distances(tmp.g, mode = "out")
    # gen.neighbs <- names(which(dd[,sampled.generalist] <= short.pathLen))
    # spec.neighbs <- names(which(dd[,sampled.specialist] <= short.pathLen))
    # neighbs <- intersect(gen.neighbs, spec.neighbs)
    # drop.nodes <- as.character(as_ids(V(g))[-which(as_ids(V(g)) %in% neighbs)])
    # tmp.g <- delete_vertices(g, drop.nodes)
    # 
    # trajectories <-
    #   get.all.shortest.paths(tmp.g,
    #                          from = sampled.generalist,
    #                          to = sampled.specialist)$res
    # we should hypothetically now be able to look at all paths if we consider a directed graph, 
    # only allowing traveral across fitness ridges, or fitness inclines (e.g. only going from one node to another if
    # fitness increases along the way - fitness decreases are not allowed).
    trajectories <-
      all_simple_paths(g,
                       from = sampled.generalist,
                       to = sampled.specialist,
                       mode = 'out')
    #Check that these do no go through haplotypes not found in hybrids (missing fitness data)
    # If any nodes do not, ignore
    for(n in 1:length(trajectories)){
      paths <- 
        list.append(paths, trajectories[[n]])       
    }
  }
  return(paths)
}

get.path.values <- 
  function(short.paths, nodes, fitness, trajectory, return.simple = F){
    
    # Get the actual values along the path
    paths <- short.paths
    cumDist <- list()
    cumDistRescale <- list()
    monoton <- c()
    MeanFit <- c()
    MeanFit.ByLen <- c()
    MedFit <- c()
    MedFit.ByLen <- c()
    abs.dist <- c()
    # too.short <- which(lengths(short.paths) <= 2)
    # if(length(too.short) >=1){
    #   short.paths <- short.paths[-too.short]
    # }
    
    # We need to initialize a vector of paths we'll later need to drop. 
    # These will include paths that are of length two, but with one NA,
    # as well as paths of length three with two NAs, or an internal NA,
    # and lastly paths of any length with internal NAs 
    toDrop <- c()
    # In some cases there are a prohibitive number of paths. If so, draw 10000 at random and calculate for these
    num.paths <- length(short.paths)
    if(num.paths > 5000){
      path.samps <- sample(1:num.paths, size = 5000)
    }else{
      path.samps <- 1:num.paths
    }
    
    for(path in path.samps){
      steps <- c()
      fits <- c()
      origin <- as_ids(short.paths[[path]])[1]
      # Here we're just getting values of fitness at each node
      # Also get the cumulative distance from origin
      for(i in 1:length(short.paths[[path]])){
        seq <- as_ids(short.paths[[path]])[i]
        fits[i] <- nodes[which(nodes$seq == seq),'fitness']
        steps[i] <- i-1
      }
      toDrop[path] <- 0
      # Now, check the conditions of the path to determine if it's suitable for determining accessibility.
      # If not, flag for dropping
      if(length(fits) == 2 & sum(is.na(fits)) >= 1){
        toDrop[path] <- 1
      }
      if(length(fits) == 3 & sum(is.na(fits)) >= 2){
        toDrop[path] <- 1
      }
      if(length(fits) >= 3 & is.na(fits[2])){
        toDrop[path] <- 1
      }
      if(length(fits) >= 4 & sum(is.na(fits[2:c(length(fits)-1)]))){
        toDrop[path] <- 1
      }
      
      # # Now, clean these up. If there are no NAs within (internal to - excluding flanking parentals) the path, continue,
      # # otherwise proceed to the next path
      # len <- length(fits)
      # na <- is.na(fits[-c(1, len)])
      # If there are NAs between the two parental nodes, ignore the path, we can't use it
      if(toDrop[path] == 1){
        paths[[path]] <- NULL
        cumDist[[path]] <- NULL
      }else{
        # otherwise, remove parental nodes if they are NAs for the measure of fitness
        na <- is.na(fits)
        if(sum(na) > 0){
          paths[[path]] <- fits[-which(is.na(fits))]
          cumDist[[path]] <- steps[-which(is.na(fits))]
        }else{
          # take all fitness values if parents aren't NAs
          paths[[path]] <- fits
          cumDist[[path]] <- steps
        }
        cumDistRescale[[path]] <- rescale(cumDist[[path]], to = c(0,1))
      }
    }
    
    # Also calculate accessibility
    if(length(paths) >= 1){
      for(i in 1:length(paths)){
        path <- paths[[i]]
        if(is.null(path)){
          MeanFit[i] <- NA
          MeanFit.ByLen[i] <- NA
          MedFit[i] <- NA
          MedFit.ByLen[i] <- NA
          monoton[i] <- NA
          abs.dist[i] <- NA
        }else{
          
          diffs <- c()
          dist <- max(cumDist[[i]])
          if(length(path) == 1){
            MeanFit[i] <- NA
            MeanFit.ByLen[i] <- NA
            MedFit[i] <- NA
            MedFit.ByLen[i] <- NA
            monoton[i] <- NA
            abs.dist[i] <- NA
          }else{
            for(x in 1:length(path)-1){
              diffs[x] <- path[x+1] - path[x]
            }
            # Here, accessibility is defined by the proportion of positive (or = zero) steps in a path
            # divided by the total length of the path
            if(all(path == cummax(path))){
              monoton[i] <- 1
            }else{
              monoton[i] <- 0
            }
            MeanFit[i] <- mean(diffs)
            MeanFit.ByLen[i] <- mean(diffs) / dist
            MedFit[i] <- median(diffs)
            MedFit.ByLen[i] <- median(diffs) / dist
            abs.dist[i] <- dist
          }
          abs.dist <- as.numeric(gsub(-Inf, 0, abs.dist))
        }
      }
    }
    
    
    access <- 
      list('MeanFit' = MeanFit,
           'MeanFit.ByLen' = MeanFit.ByLen,
           'MedFit' = MedFit,
           'MedFit.ByLen' = MedFit.ByLen,
           'monotonic' = monoton,
           'distance' = abs.dist)

    # If you're going to plot individual trajectories (i.e. in ggplot, then do the following)
    if(return.simple == F){
      if(length(paths) >= 1){
        res <- 
          data.frame(PathID = NA, 
                     Distance = NA,
                     CumulativeDistance = NA,
                     CumulativeDistanceRescale = NA,
                     Fitness = NA,
                     IsMonotonic = NA,
                     Trajectory = NA)
        for(path in which(lengths(paths) > 1)){
          path.id <- rep(path, length(paths[[path]]))
          tmp <- 
            data.frame(PathID = path.id, 
                       Distance = abs.dist[[path]],
                       CumulativeDistance = cumDist[[path]],
                       CumulativeDistanceRescale = cumDistRescale[[path]],
                       Fitness = paths[[path]],
                       IsMonotonic = rep(access$monotonic[[path]], length(paths[[path]])),
                       Trajectory = rep(trajectory, length(paths[[path]])))
          res <- 
            rbind(res, tmp)
        }
        res <- res[-1,]
        res$IsMonotonic <- factor(res$IsMonotonic, levels = c(0,1))
      }
    }else{
      if(length(paths) >= 1){
        res <- 
          data.frame(PathID = NA, 
                     Distance = NA,
                     MeanFit = NA,
                     MeanFit.ByLen = NA,
                     MedFit = NA,
                     MedFit.ByLen = NA,
                     IsMonotonic = NA,
                     Trajectory = NA)
        
        for(path in which(lengths(paths) > 1)){
          tmp <- 
            data.frame(PathID = path, 
                       Distance = abs.dist[[path]],
                       MeanFit = access$MeanFit[[path]],
                       MeanFit.ByLen = access$MeanFit.ByLen[[path]],
                       MedFit = access$MedFit[[path]],
                       MedFit.ByLen = access$MedFit.ByLen[[path]],
                       IsMonotonic = access$monotonic[[path]],
                       Trajectory = trajectory)
          res <- 
            rbind(res, tmp)
        }
        res <- res[-1,]
      }
    }
    if(exists('res') == FALSE){
      res <- NULL
    }
    return(res)
  }

check.connectivity <- function(node.pairs, g, h.haps){
  for(t in 1:nrow(node.pairs)){
    drop.nodes <- unique(c(as.character(node.pairs[,1]), as.character(node.pairs[,2])))
    drop.nodes <- drop.nodes[-which(drop.nodes %in% c(as.character(node.pairs[t,1]), as.character(node.pairs[t,2])))]
    h.keepers <- which(drop.nodes %in% h.haps)
    if(length(h.keepers > 0)){
      drop.nodes <- drop.nodes[-h.keepers]
    }
    if(length(drop.nodes) > 0){
      tmp.g <- delete.vertices(g, v = drop.nodes)
    }else{
      tmp.g <- g
    }
    
    node.pairs$conn[t] <-
      edge.connectivity(tmp.g,
                        source = as.character(node.pairs[t,1]),
                        target = as.character(node.pairs[t,2]))
    res <- list('traj' = node.pairs, 'drop' = drop.nodes)

  }
  return(res)
}

get.peak.dists <- function(peak.nodes = NULL, parent.nodes = NULL, potent.pairs = NULL, mode = c('min', 'mean')){
  dist.to.parent <- c()
  if(mode == 'min'){
    for(n in 1:length(peak.nodes)){
      peak <- peak.nodes[n]
      set.a <- which(potent.pairs$from %in% parent.nodes & 
                       potent.pairs$to == peak)
      set.b <- which(potent.pairs$to %in% parent.nodes & 
                       potent.pairs$from == peak)
      focal <- unique(c(set.a, set.b))
      if(peak %in% parent.nodes){
        dist.to.parent[n] <- min(c(potent.pairs[focal,3], 0))
      }else{
        if(sum(potent.pairs$from %in% parent.nodes) >= 1){
          dist.to.parent[n] <- min(potent.pairs[focal,3])
          
        }else{
          dist.to.parent[n] <- NA
        }
      }
    }
    if(length(na.omit(dist.to.parent)) > 0){
      dist.to.parent <- min(na.omit(dist.to.parent))
    }else{
      dist.to.parent <- NA
    }
  }else{
    if(mode == 'mean'){
      for(n in 1:length(peak.nodes)){
        peak <- peak.nodes[n]
        set.a <- which(potent.pairs$from %in% parent.nodes & 
                         potent.pairs$to %in% peak)
        set.b <- which(potent.pairs$to %in% parent.nodes & 
                         potent.pairs$from %in% peak)
        focal <- unique(c(set.a, set.b))
        if(peak %in% parent.nodes){
          dist.to.parent[n] <- mean(c(potent.pairs[focal,3], 0))
        }else{
          if(sum(potent.pairs$from %in% parent.nodes) >= 1){
            dist.to.parent[n] <- mean(potent.pairs[focal,3])
            
          }else{
            dist.to.parent[n] <- NA
          }
        }
      }
      if(length(na.omit(dist.to.parent)) > 0){
        dist.to.parent <- mean(na.omit(dist.to.parent))
      }else{
        dist.to.parent <- NA
      }
    }
  }
  return(dist.to.parent)
}


summarize.nets <-
  function(dat, genet.source, gam.res, raw = NULL, n, base, dir, spp,
           fitness, returnSumms = T, return.simple = F){
    # First we will begin using all sites, irrespective of their source
    geno.fpath <- list.files(path = base, pattern = paste0(genet.source, '-sites-', spp))
    genos <- read.table(paste0(base, geno.fpath[1]), header = T)
    genos <- genos[which(genos$ID %in% dat$ID),]
    genos <- merge(dat, genos, by = 'ID', all = T)

    #Need to do for both measures of accessibility, for each trajectory
    site.MedFit <- data.frame(colnames(genos)[-c(1:9)])
    colnames(site.MedFit) <- 'Site'
    site.MeanFit <- data.frame(colnames(genos)[-c(1:9)])
    colnames(site.MeanFit) <- 'Site'
    site.MedFit.ByLen <- data.frame(colnames(genos)[-c(1:9)])
    colnames(site.MedFit.ByLen) <- 'Site'
    site.MeanFit.ByLen <- data.frame(colnames(genos)[-c(1:9)])
    colnames(site.MeanFit.ByLen) <- 'Site'

    sites <- dat

    fname <- list.files(path = base, pattern = paste0(genet.source, '-permuations-', spp))[1]
    perms <- as.matrix(read.table(paste0(base, fname)), header = F)

    netSumm <-
      matrix(nrow = n, ncol = 29,
             dimnames = list(1:n, c('Modularity', 'NumModules', 'MeanDegree', 'MeanStrength',
                                    'MeanCloseness', 'EigenCentrality', 'MeanBetweenness',
                                    'MeanDistance', 'Diameter', 'EdgeDensity', 'NumComponents',
                                    'ProportionRealizedPaths', 'NumPeaks', 'MeanPeakFitness',
                                    'MeanPeakNeighborSlope', 'MeanNumPathsToPeaks', 'MeanPeakDegreeProb',
                                    'MinDistToPeakM', 'MinDistToPeakP', 'MeanDistToPeakM', 'MeanDistToPeakP',
                                    'MaxPeakFitness', 'MeanNumPathsToMaxPeak', 'MaxPeakDegreeProb',
                                    'MinDistToMaxPeakM', 'MinDistToMaxPeakP', 'MeanDistToMaxPeakM',
                                    'MeanDistToMaxPeakP', 'NumNodes.InNet')))

    possible.seqs <- expand.grid(0:2, 0:2, 0:2, 0:2, 0:2)
    possible.seqs <- unite(possible.seqs[,1:5], sites, sep = '')
    possible.pairs <-
      combn(x = as.character(possible.seqs$sites), m = 2, simplify = F)
    all.pairs <-
      data.frame(from = rep(NA, length(possible.pairs)),
                 to = rep(NA, length(possible.pairs)),
                 dist = rep(NA, length(possible.pairs)))

    print('Calculating all distances in nucleotides between all possible haplotypes')
    for(i in 1:length(possible.pairs)){
      a <- as.numeric(strsplit(possible.pairs[[i]][1], "")[[1]])
      b <- as.numeric(strsplit(possible.pairs[[i]][2], "")[[1]])
      pair <- rbind(a,b)
      diffs <- c()
      for(col in 1:ncol(pair)){
        if(max(pair[,col]) == 0){
          diffs[col] <- 0
        }else{
          diffs[col] <- max(pair[,col])- min(pair[,col])
        }
        all.pairs[i,1] <- possible.pairs[[i]][1]
        all.pairs[i,2] <- possible.pairs[[i]][2]
        all.pairs[i,3] <- sum(diffs)
      }
    }

    # Construct dataframes of accessibility means for both trajectories.
    AtoB.meds <-
      data.frame('Permutation' = 1:n,
                 'NumPaths' = rep(NA, n),
                 'MinDistance' = rep(NA, n),
                 'MeanDistance' = rep(NA, n),
                 'MedDistance' = rep(NA, n),
                 'MeanFit' = rep(NA, n),
                 'MeanFit.ByLen' = rep(NA, n),
                 'MedFit' = rep(NA, n),
                 'MedFit.ByLen' = rep(NA, n),
                 'NumNodes.InNet' = rep(NA, n),
                 'NumEdges.InNet' = rep(NA, n))
    
    epistasis <- 
      data.frame(
        'Permutation' = 1:n,
        'NetFitnessCorrelation' = NA,
        'NetCorrelationPvalue' = NA,
        'NumSitesSignifCorrs' = NA,
        'MeanSiteFitnessCorr_Signif' = NA,
        'MinSiteFitnessCorr_Signif' = NA,
        'MaxSiteFitnessCorr_Signif' = NA,
        'MeanSiteFitnessCorr_All' = NA,
        'MinSiteFitnessCorr_All' = NA,
        'MaxSiteFitnessCorr_All' = NA,
        'Site1_Corr' = NA,
        'Site2_Corr' = NA,
        'Site3_Corr' = NA,
        'Site4_Corr' = NA,
        'Site5_Corr' = NA,
        'Site1_Corr_P' = NA,
        'Site2_Corr_P' = NA,
        'Site3_Corr_P' = NA,
        'Site4_Corr_P' = NA,
        'Site5_Corr_P' = NA
      )
    
    for(i in 1:n){
      print(paste0("Beginning permutation ", i, " of ", n))
      # pull out the sites for this permutation
      geno <- genos[,c(1:9,perms[i,]+9)]
      
      # remove samples with missing genotypes
      geno <- geno[-which(rowSums(is.na(geno[,10:14])) >= 1),]
      
      # pull out genotype ids
      geno.ids <- data.frame(colnames(geno)[10:14])
      colnames(geno.ids) <- 'Site'
      
      # construct haplotypes
      seqs <- unite(geno[,c(10:14)], sites, sep = '')
      seqs <- cbind(geno[,1:9], seqs)
      
      # these haplotypes need to be treated as factors
      seqs$sites.factor <- as.factor(seqs$sites)
      
      # reduce the full matrix of possible pairs
      potent.pairs <- 
        all.pairs[which(all.pairs[,1] %in% seqs$sites &
                          all.pairs[,2]  %in% seqs$sites),]
      
      growth <- with(na.omit(seqs), tapply(Growth, sites.factor, mean))
      surv <- with(na.omit(seqs), tapply(survival, sites.factor, mean))
      ld1 <- with(na.omit(seqs), tapply(LD1, sites.factor, mean))
      ld2 <- with(na.omit(seqs), tapply(LD2, sites.factor, mean))
      
      # And combine
      fit <- data.frame(
        'seq' = as.character(names(growth)),
        'growth' = growth,
        'survival' = surv,
        'LD1' = ld1,
        'LD2' = ld2
      )
      rownames(fit) <- NULL
      
      if(is.null(raw)){
        fitness.trait <- 'PredictedFitness'
        # Get predicted fitness
        probs<- c()
        for(hap in 1:nrow(fit)){
          if(sum(is.na(fit[hap,4:5])) == 2){
            probs[hap] <- NA
          }else{
            index <- nn2(gam.res[,1:2], query = fit[hap,4:5])[[1]][1]
            probs[hap] <- gam.res[index,3]
          }
        }
        fit$PredictedFitness <- probs
      }else{
        fitness.trait <- raw
      }
      
      site.spp <- seqs %>% dplyr::count(sites.factor,Species)
      site.spp <- dcast(site.spp, sites.factor ~ Species, value.var = 'n')
      site.spp[is.na(site.spp)] <- 0
      site.spp$Generalist.prop <- site.spp$Generalist / rowSums(site.spp[,2:5])
      site.spp$Molluscivore.prop <- site.spp$Molluscivore / rowSums(site.spp[,2:5])
      site.spp$ScaleEater.prop <- site.spp$`Scale-eater` / rowSums(site.spp[,2:5])
      site.spp$Hybrid.prop <- site.spp$Hybrid / rowSums(site.spp[,2:5])
      site.spp <- site.spp[,c(1,6:9)]
      colnames(site.spp)[1] <- 'seq'
      
      fit <- merge(fit, site.spp, by = 'seq')
      
      #print(paste0('Estimating network'))
      # prepare for network construction
      fit.measure <- fitness
      if(fitness == "GrowSurv"){fit.measure <- 'composite'}
      
      net <- prep.network(fit, dist = 1, precalc.dists = potent.pairs, fit.measure = fit.measure)
      
      if(!is.null(net)){
        g <- graph_from_data_frame(net$pairs, directed=T, vertices=net$nodes)
        
        ########################################################################
        # Now we begin to summarize
        # We also want to identify and summarize peaks
        # To do so we need to ignore cases where nodes with data direct into parental 
        # nodes without data, as this will mask observed peaks
        # For this, we need to exclude nodes without fitness data
        dat.nodes <- net$nodes[which(!is.na(net$nodes$fitness)),]
        dat.pairs <- net$pairs[which(net$pairs$from %in% dat.nodes$seq & 
                                       net$pairs$to %in% dat.nodes$seq ),]
        
        a.nodes <- net$nodes$seq[which(net$nodes$Generalist.prop > 0)]
        m.nodes <- net$nodes$seq[which(net$nodes$Molluscivore.prop > 0)]
        p.nodes <- net$nodes$seq[which(net$nodes$ScaleEater.prop > 0)]
        h.nodes <- net$nodes$seq[which(net$nodes$Hybrid.prop > 0)]
        parent.nodes <- unique(c(a.nodes, m.nodes, p.nodes))
        
        dat.g <- graph_from_data_frame(dat.pairs, directed=T, vertices=dat.nodes)
        
        # calculate the effect of epistasis as measured by the pairwise correlation in 
        # fitness effect of mutation between neighbor genotypes. 
        
        for(x in 1:5){
          # Pull out the set of nodes we're focusing on
          obs.nodes <- unique(c(dat.pairs$from, dat.pairs$to))
          
          # Pull out the genotypes with same allele at the focal position
          obs0 <- obs.nodes[which(substr(obs.nodes, x, x) == 0)]
          obs1 <- obs.nodes[which(substr(obs.nodes, x, x) == 1)]
          obs2 <- obs.nodes[which(substr(obs.nodes, x, x) == 2)]
          
          # Get their distance from "00000", "10000", or "20000", and order them by increasing distance
          dists.0 <- all.pairs[which(all.pairs$from == "00000" & all.pairs$to %in% obs0),]
          dists.1 <- all.pairs[which(all.pairs$from == "10000" & all.pairs$to %in% obs1),]
          dists.2 <- all.pairs[which(all.pairs$from == "20000" & all.pairs$to %in% obs2),]
          
          dists.0 <- unique(merge(dists.0, dat.pairs[,c(2,5)], by = 'to'))
          dists.1 <- unique(merge(dists.1, dat.pairs[,c(2,5)], by = 'to'))
          dists.2 <- unique(merge(dists.2, dat.pairs[,c(2,5)], by = 'to'))
          
          # and combine
          tmp <- rbind(dists.0, dists.1, dists.2)[,c(3:4)]
          tmp$position <- x
          
          if(x == 1){
            fit.dists <- tmp
          }else{
            fit.dists <- rbind(fit.dists, tmp)
          }
        }
         
        net.epi <- summary(lm(fit.dists$fit2 ~ fit.dists$dist))
        epistasis$Permutation[i] <- i
        epistasis[i,'NetFitnessCorrelation'] <- net.epi$coefficients[2]
        epistasis[i,'NetCorrelationPvalue'] <- net.epi$coefficients[8]
        
        num.signif <- 0
        site.corrs <- c()
        signif.corrs <- c()
        site.ps <- c()
        for(x in 1:5){
          site.fit.dists <- fit.dists[which(fit.dists$position == x),]
          epi.site <- summary(lm(site.fit.dists$fit2 ~ site.fit.dists$dist))
          if(epi.site$coefficients[8] < 0.05){
            num.signif <- num.signif + 1
            signif.corrs <- append(signif.corrs, epi.site$coefficients[2])
          }
          site.corrs <- append(site.corrs, epi.site$coefficients[2])
          site.ps <- append(site.ps, epi.site$coefficients[8])
        }
        epistasis[i,'NumSitesSignifCorrs'] <- num.signif
        epistasis[i,11:15] <- site.corrs
        epistasis[i,16:20] <- site.ps
        if(num.signif > 0){
          epistasis[i,'MeanSiteFitnessCorr_Signif'] <- mean(signif.corrs)
          epistasis[i,'MinSiteFitnessCorr_Signif'] <- min(signif.corrs)
          epistasis[i,'MaxSiteFitnessCorr_Signif'] <- max(signif.corrs)
        }else{
          epistasis[i,'MeanSiteFitnessCorr_Signif'] <- NA
          epistasis[i,'MinSiteFitnessCorr_Signif'] <- NA
          epistasis[i,'MaxSiteFitnessCorr_Signif'] <- NA
        }
        epistasis[i,'MeanSiteFitnessCorr_All'] <- mean(site.corrs)
        epistasis[i,'MinSiteFitnessCorr_All'] <- min(site.corrs)
        epistasis[i,'MaxSiteFitnessCorr_All'] <- max(site.corrs)

        comm <- cluster_walktrap(g, steps = 5000)
        netSumm[i,1] <- modularity(comm)
        netSumm[i,2] <- length(unique(comm$membership))
        netSumm[i,3] <- mean(degree(g))
        netSumm[i,4] <- mean(strength(g))
        netSumm[i,5] <- suppressWarnings(mean(closeness(g, normalized = T)))
        netSumm[i,6] <- mean(eigen_centrality(g)$vector)
        netSumm[i,7] <- mean(betweenness(g))
        d <- distances(g)
        netSumm[i,8] <- mean(d[is.finite(d)])
        netSumm[i,9] <- diameter(g)
        netSumm[i,10] <- edge_density(g)
        netSumm[i,11] <- components(g)$no
        netSumm[i,12] <- nrow(dat.pairs) / nrow(potent.pairs)
        
        # Find peaks. These could be defined either as nodes that only have edges leading to them 
        # (i.e. all neighboring nodes have reduced fitness), but doing so excludes cases where  
        # you could have a fitness ridge atop the peak (i.e. several equally high-fitness nodes are 
        # connected, but together are higher fitness than all other adjacent nodes). 
        # So, we want to include both the single node peaks, and the fitness ridges that are surrounded 
        # by lower fitness nodes
        
        # Get the single node peaks
        ins <- V(dat.g)[degree(dat.g, mode = 'in') >= 2] # they should at least be connected by two edges
        outs <- V(dat.g)[degree(dat.g, mode = 'out') > 0] # these have edges leading out - exclude
        peaks <- ins[-which(ins %in% outs)]
        peak.nodes <- names(peaks)
        
        # Now get the fitness ridges
        possible.ridges <- names(ins[which(ins %in% outs)])
        ridge.pairs <- 
          dat.pairs[which(dat.pairs$from %in% possible.ridges | 
                            dat.pairs$to %in% possible.ridges),]
        # Fitness for one 'from' node can never equal zero
        ridge.pairs <- ridge.pairs[-which(ridge.pairs$fit1 == 0),]
        # Exclude any node in the 'from' column if it's fitness is less than the node in the 'to' columns
        ridge.pairs <- ridge.pairs[-which(ridge.pairs$fit1 < ridge.pairs$fit2),]
        ridge.nodes <- unique(c(ridge.pairs$from, ridge.pairs$to))

        # And combine with the peaks.
        peak.nodes <- unique(c(peak.nodes, ridge.nodes))
        
        # IF PEAKS EXIST, PUT INTO DATAFRAME, OTHERWISE PUT IN NA
        # Get the number
        netSumm[i,13] <-  length(peak.nodes)
        
        # Identify which nodes are the neighbors to these peak, 
        # calculate the mean fitness step to these peaks,
        # the mean number of paths feeding into them,
        # and the mean distance from peaks to each specialist
        if(length(peak.nodes) > 0){
          # Need to pull out any pair including the peak - in this dataframe, 
          # the "to" node isn't necessarily higher fitness
          peak.pairs <- dat.pairs[which(dat.pairs$to %in% peak.nodes | 
                                          dat.pairs$from %in% peak.nodes),]
          fit.diffs <- c()
          peak.fit <- c()
          for(x in 1:nrow(peak.pairs)){
            from <- dat.nodes$seq[which(dat.nodes$seq == peak.pairs[x,1])]
            to <- dat.nodes$seq[which(dat.nodes$seq == peak.pairs[x,2])]
            fit.diffs[x] <- 
              dat.nodes$fitness[which(dat.nodes$seq == to)] - 
              dat.nodes$fitness[which(dat.nodes$seq == from)]
          }
          peak.dat <- dat.nodes[which(dat.nodes$seq %in% peak.nodes),]
          netSumm[i,14] <- mean(peak.dat$fitness)
          netSumm[i,15] <- mean(fit.diffs)
          
          npaths.to.peak <- c()
          for(x in 1:length(peak.nodes)){
            npaths.to.peak[x] <- nrow(peak.pairs[which(peak.pairs$to == peak.nodes[x]),])
          }
          netSumm[i,16] <- mean(npaths.to.peak)
          
          # Get the mean degree probability for peaks (i.e., are these peaks unusually connected?)
          peak.degree.probs <- c()
          degree.probs <- degree_distribution(dat.g)
          for(x in 1:length(npaths.to.peak)){
            peak.degree.probs[x] <- degree.probs[npaths.to.peak[x]+1]
          }
          netSumm[i,17] <- mean(peak.degree.probs)
          netSumm[i,18] <- get.peak.dists(peak.nodes, m.nodes, potent.pairs, mode = 'min')
          netSumm[i,19] <- get.peak.dists(peak.nodes, p.nodes, potent.pairs, mode = 'min')
          netSumm[i,20] <- get.peak.dists(peak.nodes, m.nodes, potent.pairs, mode = 'mean')
          netSumm[i,21] <- get.peak.dists(peak.nodes, p.nodes, potent.pairs, mode = 'mean')
          
          # Now the greatest peak
          max.peak <- peak.dat$seq[which(peak.dat$fitness == max(peak.dat$fitness))]
          netSumm[i,22] <- max(peak.dat$fitness)
          # Because of ridges/max fitness of 1 for survival, we can't assume a single maximum
          # Need to summarize across
          npaths.to.max <- c()
          max.peak.degree.p <- c()
          min.m.max.peak.dist <- c()
          mean.m.max.peak.dist <- c()
          min.p.max.peak.dist <- c()
          mean.p.max.peak.dist <- c()
          for(p in 1:length(max.peak)){
            npaths.to.max[p] <- nrow(peak.pairs[which(peak.pairs$to == max.peak[p]),])
            max.peak.degree.p[p] <- degree.probs[npaths.to.max+1]
            min.m.max.peak.dist[p] <- get.peak.dists(max.peak[p], m.nodes, potent.pairs, mode = 'min')
            min.p.max.peak.dist[p] <- get.peak.dists(max.peak[p], p.nodes, potent.pairs, mode = 'min')
            mean.m.max.peak.dist[p] <- get.peak.dists(max.peak[p], m.nodes, potent.pairs, mode = 'mean')
            mean.p.max.peak.dist[p] <- get.peak.dists(max.peak[p], p.nodes, potent.pairs, mode = 'mean')
          }
          netSumm[i,23] <- mean(npaths.to.max)
          netSumm[i,24] <- mean(max.peak.degree.p)
          netSumm[i,25] <- min(min.m.max.peak.dist)
          netSumm[i,26] <- min(min.p.max.peak.dist)
          netSumm[i,27] <- mean(mean.m.max.peak.dist)
          netSumm[i,28] <- mean(mean.p.max.peak.dist)
          netSumm[i,29] <- nrow(dat.nodes)
          
        }else{
          netSumm[i,14:28] <- NA
        }
        
        # Now, we need to identify all parental nodes, and come up with all combinations between them (among different species)
        
        
        # Pull out the paths and remove paths that are a single step between parents (these have no data for fitness)
        #print('Extracting and summarizing paths')
        
        # Now, we need to identify all paths between parental species. 
        # So, come up with all pairwise combinations between species
        # We also want to exclude parent nodes that dont have neighboring hybrids with fitness data
        if(spp == 'M'){
          spec.nodes <- 
            unique(net$pairs[which(net$pairs$to %in% m.nodes & 
                                     net$pairs$from %in% h.nodes),2])
          a.nodes  <- 
            unique(net$pairs[which(net$pairs$from %in% a.nodes & 
                                     net$pairs$to %in% h.nodes),1])
          dups <- which(a.nodes %in% as.character(spec.nodes))
          if(length(dups) > 0){
            a.nodes <- a.nodes[-which(a.nodes %in% as.character(spec.nodes))]
          }
          trajectory = 'A2M'
        }else{
          spec.nodes <- 
            unique(net$pairs[which(net$pairs$to %in% p.nodes & 
                                     net$pairs$from %in% h.nodes),2])
          a.nodes  <- 
            unique(net$pairs[which(net$pairs$from %in% a.nodes & 
                                     net$pairs$to %in% h.nodes),1])
          dups <- which(a.nodes %in% as.character(spec.nodes))
          if(length(dups) > 0){
            a.nodes <- a.nodes[-which(a.nodes %in% as.character(spec.nodes))]
          }
          trajectory = 'A2P'
        }
        if(length(a.nodes) > 0 && length(spec.nodes) > 0){
          AtoB.path.summ <- NULL
          AtoB.path.steps <- NULL
          AtoB <- expand.grid(a.nodes, spec.nodes)
          AtoB <- AtoB[which(as.character(AtoB$Var1) %in% as_ids(V(g)) & 
                               as.character(AtoB$Var2) %in% as_ids(V(g))),]
          AtoB$conn <- NA
          # In same cases, generalist and specialist nodes will overlap. Because there are fewer
          # specialist nodes on average than generalist nodes, remove the overlapping node from the 
          # list of generalist nodes, but keep for specialist.
          shared <- which(AtoB$Var1 %in% AtoB$Var2)
          if(length(shared) > 0){
            AtoB <- AtoB[-shared,]
          }
          if(nrow(AtoB) > 0){
            for(node in 1:nrow(AtoB)){
              AtoB$conn[node] <-
                edge.connectivity(g,
                                  source = as.character(AtoB[node,1]),
                                  target = as.character(AtoB[node,2]))
            } 
            AtoB <- AtoB[which(AtoB$conn > 0),]
          }
          
          
          # if(nrow(AtoB) > 0){
          #   # exclude paths that are a single step
          #   AtoB$dist <- NA
          #   for(hap in 1:nrow(AtoB)){
          #     AtoB$dist[hap] <- 
          #       all.pairs[which(all.pairs$from == AtoB$Var1[hap] & 
          #                         all.pairs$to == AtoB$Var2[hap] |
          #                         all.pairs$from == AtoB$Var2[hap] & 
          #                         all.pairs$to == AtoB$Var1[hap]),3]
          #   }
          #   AtoB <- AtoB[which(AtoB$dist > 1),1:3]
          # }
          
          if(nrow(AtoB) > 0){
            if(identical(a.nodes, spec.nodes) == T | length(a.nodes) == 0 | length(spec.nodes) == 0){
              AtoB.paths <- NA
            }else{
              # Determine which parental nodes are connected by a contiguous string of hybrids
              # (nodes may be connected, but with intervening parental nodes without data)
              # We want to make sure that we only choose paths that actually have complete data 
              # so we can determine accessibility
              #print('Sampling parental nodes and obtaining trajectories')
              conn <- check.connectivity(AtoB, g, h.nodes)
              AtoB <- conn$traj
              drop.nodes <- conn$drop
              
              # Now, remove pairs of nodes that are not fully connected
              AtoB <- AtoB[which(AtoB$conn > 0),]
              
              rand.pair <- sample(1:nrow(AtoB), 1)
              rand.A <- as.character(AtoB[rand.pair,1])
              rand.B <- as.character(AtoB[rand.pair,2])
              
              drop.nodes <- drop.nodes[-which(drop.nodes %in% c(rand.A, rand.B))]
              h.keepers <- which(drop.nodes %in% h.nodes)
              
              if(length(h.keepers > 0)){
                drop.nodes <- drop.nodes[-h.keepers]
              }
              if(length(drop.nodes) > 0){
                tmp.g <- delete.vertices(g, v = drop.nodes)
              }else{
                tmp.g <- g
              }
              
              AtoB.paths <-
                get.trajectories(tmp.g, spp.pair = AtoB, parent.haps = parent.nodes, 
                                 sampled.generalist = rand.A, sampled.specialist = rand.B)
              
              # Check which nodes along these paths we have hybrid fitness data for.
              # At a minimum, all interior nodes need to have fitness data
              for(x in 1:length(AtoB.paths)){
                h <- which(names(AtoB.paths[[x]]) %in% h.nodes)
                # We can allow a two node path (single step) only if both parental nodes have hybrid fitness data
                if(length(AtoB.paths[[x]]) == 2){
                  enough.h <- length(h) == 2
                  if(enough.h == FALSE){
                    AtoB.paths[[x]] <- NA
                  }
                }
                # Then a check for three node paths
                # If paths are only three nodes long, we need at least two hybrid nodes (data for one parent) - 
                # other wise, we have no information about change in fitness, and the second (internal) node must be hybrid
                if(length(AtoB.paths[[x]]) == 3){
                  enough.h <- length(h) >= 2
                  if(enough.h == FALSE | 2 %in% h == FALSE){
                    AtoB.paths[[x]] <- NA
                  }
                }
                # If paths are four steps or longer (including two parent nodes), we can allow paths without parent
                # fitness data, but all internal nodes must have fitness data. 
                if(length(AtoB.paths[[x]]) >= 4){
                  internal.h <- names(AtoB.paths[[x]][2:(length(AtoB.paths[[x]])-1)])
                  if(all(internal.h %in% h.nodes) == FALSE){
                    AtoB.paths[[x]] <- NA
                  }
                }
              }
              AtoB.paths[which(is.na(AtoB.paths))] <- NULL
              
              # AtoB.paths <- AtoB.paths[sapply(AtoB.paths, length)>3]
              if(length(AtoB.paths) > 0){
                numPaths <- length(AtoB.paths)
                min.dist <- min(lengths(AtoB.paths))
                med.dist <- median(lengths(AtoB.paths))
                mean.dist <- mean(lengths(AtoB.paths))
              }else{
                AtoB.paths <- NA
              }
            }
          }else{
            AtoB.paths <- NA
          }
          
          # Then, obtain values for growth at each step of these paths
          if(sum(is.na(AtoB.paths)) < length(AtoB.paths)){
            #print('Getting fitness values along paths')
            AtoB.path.summ <-
              get.path.values(short.paths = AtoB.paths, nodes = net$nodes[which(net$nodes$seq %in% names(V(tmp.g))),],
                              fitness = fit.measure, trajectory = trajectory,
                              return.simple = T)
            AtoB.path.summ$Permutation <- i
            
            AtoB.path.steps <-
              get.path.values(short.paths = AtoB.paths, nodes = net$nodes[which(net$nodes$seq %in% names(V(tmp.g))),],
                              fitness = fitness.trait, trajectory = trajectory)
            AtoB.path.steps$Permutation <- i
            
            if(is.null(AtoB.path.summ)){
              numPaths <- NA
            }
          }
        }else{
          AtoB.path.summ <- NULL
          AtoB.path.steps <- NULL
        }
      }else{
        AtoB.path.summ <- NULL
        AtoB.path.steps <- NULL
      }
      
      
      if(is.null(AtoB.path.summ)){
        AtoB.meds[i,c(2:11)] <- NA
        
        # All paths go through parent node, so can't assess fitness - assign as NA
        geno.ids[,2] <- NA
        colnames(geno.ids)[2] <- paste0('MedFit-', i)
        site.MedFit <- merge(site.MedFit, geno.ids, by = 'Site', all = T)
        
        geno.ids[,2] <- NA
        colnames(geno.ids)[2] <- paste0('MeanFit-', i)
        site.MeanFit <- merge(site.MeanFit, geno.ids, by = 'Site', all = T)
        
        # Then with median and mean fitness scaled by path length
        geno.ids[,2] <- NA
        colnames(geno.ids)[2] <- paste0('MedFit.ByLen-', i)
        site.MedFit.ByLen <- merge(site.MedFit.ByLen, geno.ids, by = 'Site', all = T)
        
        geno.ids[,2] <- NA
        colnames(geno.ids)[2] <- paste0('MeanFit.ByLen-', i)
        site.MeanFit.ByLen <- merge(site.MeanFit.ByLen, geno.ids, by = 'Site', all = T)
        
      }else{
        #meds <- colMedians(as.matrix(AtoB.path.summ[,c(3:4)]))
        AtoB.meds[i,c(2:11)] <- c(numPaths, min.dist, mean.dist, med.dist, 
                                 mean(AtoB.path.summ$MeanFit), mean(AtoB.path.summ$MeanFit.ByLen), 
                                 median(AtoB.path.summ$MedFit), median(AtoB.path.summ$MedFit.ByLen),
                                 length(V(g)), length(E(g)))
        
        # First using median and mean fitness
        geno.ids[,2] <- median(AtoB.path.summ$MedFit)
        colnames(geno.ids)[2] <- paste0('MedFit-', i)
        site.MedFit <- merge(site.MedFit, geno.ids, by = 'Site', all = T)
        
        geno.ids[,2] <- mean(AtoB.path.summ$MeanFit)
        colnames(geno.ids)[2] <- paste0('MeanFit-', i)
        site.MeanFit <- merge(site.MeanFit, geno.ids, by = 'Site', all = T)
        
        # Then with median and mean fitness scaled by path length
        geno.ids[,2] <- median(AtoB.path.summ$MedFit.ByLen)
        colnames(geno.ids)[2] <- paste0('MedFit.ByLen-', i)
        site.MedFit.ByLen <- merge(site.MedFit.ByLen, geno.ids, by = 'Site', all = T)
        
        geno.ids[,2] <- mean(AtoB.path.summ$MeanFit.ByLen)
        colnames(geno.ids)[2] <- paste0('MedFit.ByLen-', i)
        site.MedFit.ByLen <- merge(site.MedFit.ByLen, geno.ids, by = 'Site', all = T)
      }

      if(i == 1){
        if(!exists('AtoB.path.summ') | is.null(AtoB.path.summ)){
          AtoB.path.summ <- data.frame(PathID = NA, Distance = NA, MeanFit = NA, MeanFit.ByLen = NA, 
                                        MedFit = NA, MedFit.ByLen = NA, IsMonotonic = NA, Trajectory = trajectory, 
                                        Permutation = i)
          net.summ.res <- AtoB.path.summ
        }else{
          net.summ.res <- AtoB.path.summ
        }
      }
      if(i > 1){
        if(!exists('AtoB.path.summ') | is.null(AtoB.path.summ)){
          AtoB.path.summ <- data.frame(PathID = NA, Distance = NA, MeanFit = NA, MeanFit.ByLen = NA, 
                                        MedFit = NA, MedFit.ByLen = NA, IsMonotonic = NA, Trajectory = trajectory, 
                                        Permutation = i)
          net.summ.res <- rbind(net.summ.res, AtoB.path.summ)
        }else{
          net.summ.res <- rbind(net.summ.res, AtoB.path.summ)
        }
      }

      if(i == 1){
        if(!exists('AtoB.path.steps') | is.null(AtoB.path.steps)){
          AtoB.path.steps <- data.frame(PathID = NA, Distance = NA, CumulativeDistance = NA, CumulativeDistanceRescale = NA, 
                                        Fitness = NA, IsMonotonic = NA, Trajectory = trajectory, 
                                        Permutation = i)
          net.step.res <- AtoB.path.steps
        }else{
          net.step.res <- AtoB.path.steps
        }
      }
      if(i > 1){
        if(!exists('AtoB.path.steps') | is.null(AtoB.path.steps)){
          AtoB.path.steps <- data.frame(PathID = NA, Distance = NA, CumulativeDistance = NA, CumulativeDistanceRescale = NA, 
                                        Fitness = NA, IsMonotonic = NA, Trajectory = trajectory, 
                                        Permutation = i)
          net.step.res <- rbind(net.step.res, AtoB.path.steps)
        }else{
          net.step.res <- rbind(net.step.res, AtoB.path.steps)
        }
      }
        

      if(i %% 50 == 0){
        write.table(net.step.res, paste0(dir, genet.source, '/PathStepSummary-', genet.source, '-', spp, '-', fitness, '.tsv'),
                    col.names = T, row.names = F, quote = F, sep = '\t')
        write.table(net.summ.res, paste0(dir, genet.source, '/PathSummary-', genet.source, '-', spp, '-', fitness, '.tsv'),
                    col.names = T, row.names = F, quote = F, sep = '\t')
        write.table(netSumm, paste0(dir, genet.source, '/NetworkSummaries-', genet.source, '-', spp, '-', fitness,  '.tsv'),
                    col.names = T, row.names = F, quote = F, sep = '\t')

        write.table(AtoB.meds, paste0(dir, genet.source, '/Ato', spp, '-MedMeanFit-PerNet-', genet.source, '-', spp, '-', fitness,  '.tsv'),
                    col.names = T, row.names = F, quote = F, sep = '\t')
        write.table(site.MedFit, paste0(dir, genet.source, '/PerSite-Snp-MedFit-', genet.source, '-', spp, '-', fitness,  '.tsv'),
                    col.names = T, row.names = F, quote = F, sep = '\t')
        write.table(site.MedFit.ByLen, paste0(dir, genet.source, '/PerSite-MedFitByLen-', genet.source, '-', spp, '-', fitness,  '.tsv'),
                    col.names = T, row.names = F, quote = F, sep = '\t')
        write.table(site.MeanFit, paste0(dir, genet.source, '/PerSite-Snp-MeanFit-', genet.source, '-', spp, '-', fitness,  '.tsv'),
                    col.names = T, row.names = F, quote = F, sep = '\t')
        write.table(site.MeanFit.ByLen, paste0(dir, genet.source, '/PerSite-MeanFitByLen-', genet.source, '-', spp, '-', fitness,  '.tsv'),
                    col.names = T, row.names = F, quote = F, sep = '\t')
        write.table(epistasis, paste0(dir, genet.source, '/EpistasisSummaries-', genet.source, '-', spp, '-', fitness,  '.tsv'),
                    col.names = T, row.names = F, quote = F, sep = '\t')
      }
    }
    
    return(netSumm)
  }


















































# get.cumDists <- 
#   function(short.paths, dists, rescale = T){
#     # This function determines how far each step is from the origin
#     cum.dist <- list()
#     for(i in 1:length(short.paths)){
#       path <- short.paths[[i]]
#       dist <- c()
#       for(n in 1:(length(path)-1)){
#         from <- as_ids(short.paths[[i]])[1]
#         to <- as_ids(short.paths[[i]])[n+1]
#         a <- which(dists$from == to | dists$to == to)
#         b <- which(dists$from == from | dists$to == from)
#         dist[n] <- all.dists[a[which(a %in% b)],3]
#       }
#       if(rescale == TRUE){
#         cum.dist[[i]] <- rescale(dist, to = c(0,1))
#       }else{
#         cum.dist[[i]] <- dist
#         }
#       
#     }
#     return(cum.dist)
#   }

# Accessibility = proportion changes increase fitness / length
calc.accessibility <- 
  function(paths, path.ids, dists){
    # Here we're just calculating the relative accessibility of the path (number of postitive steps, penalized for path length)
    access <- c()
    monoton <- c()
    for(i in 1:length(paths)){
      path <- paths[[i]]
      diffs <- c()
      cum.dists <- list()
      from <- as_ids(path.ids[[i]])[1]
      to <- as_ids(path.ids[[i]])[length(path.ids[[i]])]
      a <- which(dists$from == to | dists$to == to)
      b <- which(dists$from == from | dists$to == from)
      dist <- all.dists[a[which(a %in% b)],3]
      if(length(path) == 1){
        access[i] <- NA
        monoton[i] <- NA
      }else{
        for(x in 1:length(path)-1){
          diffs[x] <- path[x+1] - path[x]
        }
        # Here, accessibility is defined by the proportion of positive (or = zero) steps in a path
        # divided by the total length of the path
        if(sum(diffs >= 0)/dist == T){
          monoton[i] <- 1
        } else{
          monoton[i] <- 0
        }
        access[i] <- (sum(diffs >= 0) / dist) / dist
      }
    }
    access <- list('accessibility' = access,
                   'monotonic' = monoton)
    return(access)
  }

get.accessible <- 
  function(paths, min.max.paths, access){
    # Here we're simply identifying paths that are accessible (monotonic)
    # This function returns a list of paths that are monotonic, providing both the actual nodes, and their values along the path
    accessible <- list()
    accessible.path.nodes <- list()
    if(sum(access$monotonic == 1) > 0){
      mono <- which(access$monotonic == 1)
      for(i in 1:length(mono)){
        accessible[[i]] <- paths[[mono[i]]]
        accessible.path.nodes[[i]] <- min.max.paths[[mono[i]]]
      }
      accessible <- list('nodes' = accessible.path.nodes,
                         'vals' = accessible)
      return(accessible)
    }else{print('No accessible paths!')}
  }


get.shared.nodes <- 
  function(accessible){
    set <- length(accessible$nodes)
    
    # Get the number of nodes within each accessible (monotonic) path (excluding the min and max) that are shared,
    # and do this for each pair of accessible paths.
    shared <- c()
    
    comparisons <- combn(x = length(accessible$nodes), m = 2)
    for(i in 1:ncol(comparisons)){
      shared[i] <- sum(accessible$nodes[[comparisons[1,i]]] %in% accessible$nodes[[comparisons[2,i]]])-2
    }
    shared <- list('comparisons' = comparisons,
                   'numShared' = shared)
    return(shared)
  }

get.ruggedness <- 
  function(g, nodes, fitness = 'PredictedGrowth'){
    # Here we want to calculate two metrics of 'ruggedness':
    # number of peaks (one genotype away in any direction is reduced fitness)
    # number of sinks (one genotype away in any direction is increased fitness)
    # Note that these are a property of the graph as a whole, not of individual 
    # paths (# peaks and sinks will be equal otherwise)
    
    num.peaks <- 0
    num.sinks <- 0
    for(i in 1:length(V(g))){
      node <- V(g)[i]
      neighs <- neighbors(g, node)
      diffs <- c()
      for(x in 1:length(neighs)){
        if(length(neighs) >= 2){
          a <- as_ids(node)
          b <- as_ids(neighs[x])
          diffs[x] <- nodes[which(nodes$seq == a),fitness] - nodes[which(nodes$seq == b),fitness]
        }else{diffs <- 0}
      }
      if(all(diffs > 0)){
        num.peaks <- num.peaks + 1
      }
      if(all(diffs < 0)){
        num.sinks <- num.sinks + 1
      }
    }
    
    ratio <- num.peaks / num.sinks
    ruggedness <- list('numPeaks' = num.peaks,
                       'numSinks' = num.sinks,
                       'ratio' = ratio)
    return(ruggedness)
  }

get.chains <- 
  function(network = g, nodes, fitness){
    # Here, we want to find paths in which for a given node, there is only a single 
    # neighboring genotype that can lead to an increase in fitness.
    chains <- 
      data.frame(from = NA,
                 to = NA,
                 from.val = NA,
                 to.val = NA)
    
    for(i in 1:length(V(g))){
      node <- V(g)[i]
      neighs <- neighbors(g, node)
      from <- c()
      to <- c()
      for(x in 1:length(neighs)){
        a <- as_ids(node)
        b <- as_ids(neighs[x])
        
        chain <- nodes[which(nodes$seq == a),fitness] < nodes[which(nodes$seq == b),fitness]
        if(chain == TRUE){
          chains <- rbind(chains, 
                          c(a, b, nodes[which(nodes$seq == a),fitness], 
                            nodes[which(nodes$seq == b),fitness]))
        }
      }
    }
    chains <- chains[-1,]
    chains[,3] <- as.numeric(chains[,3])
    chains[,4] <- as.numeric(chains[,4])
    return(chains)
  }

netCols <- 
  function(net, fitness = 'PredictedGrowth'){
    if(min(net$nodes[,fitness]) == 0){
      cols <- rev(RColorBrewer::brewer.pal(n = 11, name = 'RdYlBu'))[6:11]
    }else{
      cols <- rev(RColorBrewer::brewer.pal(n = 11, name = 'RdYlBu'))
    }
    rgb2hex <- function(vals) rgb(vals[,1], vals[,2], vals[,3], maxColorValue = 255)
    getCols <- colorRamp(cols)
    #cols <- getCols(rescale(fitness, to = c(0,1)))
    cols <- getCols(rescale(net$nodes[,fitness][which(net$nodes[,fitness] != 0)], to = c(0,1)))
    cols <- rgb2hex(cols)
    cols <- c(rep("#313695", length(which(net$nodes[,fitness] == 0))), cols)
    return(cols)
  }

