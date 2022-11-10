# Final script for analysis
#'
#'Analysis process:
#'1. Obtain site trends
#'2. Obtain site diversities
#'3. Obtain species trends

browseVignettes("vegan")

rm(list=ls())

source("R/cwac_functions.r")
library(vegan)

# load relevant data
load('data/robertsDB.RData')
load('data/barberspan.counts.RData')

# do analysis on species with more than 6 counts only
barbernas <- as.vector(apply(bp_counts, 2, function(x){
  
  return()
  
}))

bp_counts <- bp_counts[,barbernas>6]

#### getting site trend ####

# BARBERSPAN

# tallying counts for overall analysis of barberspan
allCounts <- apply(bp_counts[,1:(length(bp_counts)-2)], 1, sum, na.rm = T)
allCounts <- as.data.frame(cbind(allCounts,
                                 'Season' = bp_counts$Season,
                                 'Year' = bp_counts$Year))
allCounts[allCounts==0] <- NA

allCounts$allCounts <- as.numeric(allCounts$allCounts)

# run jags on barberspan area
allJags <- jags_analysis(allCounts)

# plot
ts_jag_plot(allJags, allCounts, "Barberspan Combined Waterbird Population")

# BLESBOK

allBles <- apply(bles_counts[,1:(length(bles_counts)-2)], 1, sum, na.rm = T)

allBles <- as.data.frame(cbind(allBles,
                                 'season' = bles_counts$season,
                                 'startDate' = bles_counts$startDate))

allBles$allBles <- as.numeric(allBles$allBles)

# run jags on barberspan area
blessJags <- jags_analysis(allBles)

# plot
ts_jag_plot(blessJags, allBles, 2)


#### getting site diversity ####

# BARBERSPAN

# generate diversities
hill <- hill_nums(lessc)

# plot
plot_hills(hill, unique(counts$startDate))

# BLESBOK
bles_species <- get_local_species('gauteng', '26162830')

# getting counts of species
bles_counts <- species_counts('gauteng', '26162830', bles_species$id)
bles_counts <- as.data.frame(bles_counts)

blesnas <- as.vector(apply(bles_counts, 2, function(x){
  return(sum(!is.na(x)))
}))

bles_counts <- bles_counts[,blesnas>6]
bleshill <- hill_nums(bles_counts)
plot_hills(bleshill, unique(bles_counts$startDate))


#### species trends ####

# BARBERSPAN

# loop through species and save output in list
barber.jags <- list()

# species with more than 6 counts
aspec <- as.data.frame(cbind(bp_counts[,"56"],
                             'Season' = bp_counts$Season,
                             'Year' = bp_counts$Year))

ajag <- jags_analysis(aspec)
comjag <- jags_com(aspec)

aplot <- ts_jag_plot(ajag, aspec, "a plot")
aplot$s.plot
aplot$w.plot

plot_resident(comjag, aspec, "a title")

ms <- lapply(barber.jags, function(x){
  return(exp(x$mean$mu_t))
})

# apply shannons to each
# gives the best and most likely results but we won't be able to have CI
shanms <- exp(diversity(as.data.frame(ms), index = "shannon"))
simpsms <- diversity(as.data.frame(ms), index = "simpson")


####  ########################################################
# use sims lists to calculate the shannon entropy/diversity of barberspan
summer_sims <- lapply(barber.jags, function(x){
  return(x$samples[,1:26])
})

# test showing mean(exp()s) != exp(means)
test <- c()
real <- c()
testie <- list()
realms <- list()
for(j in 1:66){
  for(i in 1:26){
    test <- c(test, mean(exp(unlist(summer_sims[[j]][,i]))))
    real <- c(real, exp(mean(unlist(summer_sims[[j]][,i]))))
  }
  testie[[j]] <- test
  realms[[j]] <- real
  test <- c()
  real <- c()
}

ms <- c()
lows <- c()
ups <- c()

simps <- c()
lowSimps <- c()
upSimps <- c()

# grouping like columns together
for(i in 1:(nrow(counts)/2)){
  first.year <- lapply(summer_sims, function(x){
    return(unlist(x[,i]))
  })
  
  first.year <- as.data.frame(first.year)
  
  year1.shan <- exp(diversity(exp(first.year), index = "shannon"))
  year1.simps <- diversity(exp(first.year), index = "simpson")
  
  ms <- c(ms, data.frame(mean(year1.shan)))
  lows <- c(lows, quantile(year1.shan, probs = c(0.025, 0.975))[1])
  ups <- c(ups, quantile(year1.shan, probs = c(0.025, 0.975))[2])
  
  simps <- c(simps, mean(year1.simps))
  lowSimps <- c(lowSimps, quantile(year1.simps, props = c(0.025, 0.975))[2])
  upSimps <- c(upSimps, quantile(year1.simps, props = c(0.025, 0.975))[4])
  
}

summer_shans <- data.frame('estimate' = unlist(ms),
                           'simps' = unlist(simps),
                           'upsimps' = unlist(upSimps),
                           'lowsimps' = unlist(lowSimps),
                           'ups' = unlist(ups),
                           'lows' = unlist(lows),
                           'years' = unique(counts$startDate))

shans_p <- ggplot(summer_shans, aes(x = years)) +
  
  #geom_ribbon(aes(ymin = lower, ymax = upper), fill = "gray80") +
  #geom_line(aes(y = estimated, color = "grey1"), lwd = 1, lty = 2) +
  geom_line(aes(y = estimate, color = 'red'), lwd = 1, lty = 1) +
  geom_line(aes(y = lows, colour = 'red'), lwd = 1, lty = 3) +
  geom_line(aes(y = ups, colour = 'red'), lwd = 1, lty = 3) +
  
  # scale_color_identity(guide = "legend",
  #                      name = "",
  #                      labels = c("Equitability", "H")) +
  
  labs(title = "Shannons Index(S)")

simps_p <- ggplot(summer_shans, aes(x = years)) +
  geom_line(aes(y = simps, color = 'red'), lwd = 1, lty = 1) +
  geom_line(aes(y = lowsimps, colour = 'red'), lwd = 1, lty = 3) +
  geom_line(aes(y = upsimps, colour = 'red'), lwd = 1, lty = 3) +
  
  # scale_color_identity(guide = "legend",
  #                      name = "",
  #                      labels = c("Equitability", "H")) +
  
  labs(title = "Simpson Index(S)")

plot_grid(shans_p, simps_p, nrow = 2)  

########################################################################

ts_jag_plot(barber.jags[[63]],
            as.data.frame(cbind(morecounts[,63],
                                'season' = counts$season,
                                'startDate' = counts$startDate)),
            2)


##### Species trends based on nest site #####

robertsDB %>% group_by(`Nest site`) %>% count()

nestSite_trends <- function(traitsdf, countsdf, site, title){

  # traitsdf <- robertsDB
  # countsdf <- bp_counts
  # site <- 'Aquatic'
  # title <- "a title"
  
  # extract species with certain trait
  Ids <- traitsdf[which(traitsdf$`Nest site`==site),]$SppRef
  t <- bp_counts[,Ids]
  
  # combine the counts
  combCounts <- apply(t, 1, sum, na.rm = T)
  combCounts <- as.data.frame(cbind(combCounts,
                                'Season' = countsdf$Season,
                                'Year' = countsdf$Year))
  
  combCounts[combCounts==0] <- NA
  
  # run jags on them
  jagmod <- jags_analysis(combCounts)
  
  # plot them
  jPlot <- ts_jag_plot(jagmod, combCounts, title)
  
  # run jags on each species in trait group
  trait_species <- list()
  
  for (i in 1:length(t)){
    trait_species[[i]] <- jags_analysis(as.data.frame(cbind(t[,i],
                                                      'Season' = countsdf$Season,
                                                      'Year' = countsdf$Year)))
  }
  
  # extract beta values
  betas <- lapply(trait_species, function(x){
    return(x$mean$p)
  })
  lower <- lapply(trait_species, function(x){
    return(x$q2.5$p)
  })
  upper <- lapply(trait_species, function(x){
    return(x$q97.5$p)
  })
  
  betas <- as.data.frame(betas)
  colnames(betas) <- colnames(t)
  
  lower <- as.data.frame(lower)
  colnames(lower) <- colnames(t)
  
  upper <- as.data.frame(upper)
  colnames(upper) <- colnames(t)
  
  # get average trend per year
  betas.avg <- apply(betas, 1, mean, rm.na = T)
  l.avg <- apply(lower, 1, mean, rm.na = T)
  u.avg <- apply(upper, 1, mean, rm.na = T)

  # calc LPI
  # lpi <- c(1)
  # for(i in 2:length(betas.avg)){
  #   lpi <-  c(lpi, (lpi[i-1] * exp(betas.avg[i])))
  # }
  
  lpidf <- as.data.frame(cbind("lpi" = betas.avg,
                               'lower' = l.avg,
                               "upper" = u.avg,
                               "Year" = unique(countsdf$Year)))
  
  lpi.plot <- ggplot(lpidf, aes(x = Year, group = 1)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "gray80") +
    geom_line(aes(y = lpi, color = "grey1"), lwd = 1, lty = 2) +
    scale_x_continuous(breaks = seq(min(lpidf$Year), max(lpidf$Year), by = 1)) +
    theme(axis.title.x=element_blank(),
          axis.text.x = element_text(angle = 90),
          legend.position = 'none')
  
  return(list("jag.overall" = jagmod,
              "jag.mods" = trait_species,
              "jag.plot" = jPlot,
              "ind" = betas.avg,
              "ind2.5" = l.avg,
              "ind97.5" = u.avg,
              "ind.plot" = lpi.plot))
}

# aquatic nest type
aq.nest <- nestSite_trends(robertsDB, bp_counts, 'Aquatic', 'Aquatic Nest Type')

aq.nest$ind.plot
plot_grid(aq.nest$jag.plot$s.plot, aq.nest$jag.plot$w.plot, aq.nest$lpi.plot,
          ncol = 1,
          nrow = 3,
          rel_heights = c(1/2, 1/2, 1/4),
          rel_widths = c(1,1,1),
          left = "Log Population",
          bottom = "Year",
          align = 'v',
          axis = 'l')
#-----------------------------------------------------------------------------
# cavity nest type
cav.nest <- nestSite_trends(robertsDB, bp_counts, 'Cavity', 'Cavity Nest Type')

cav.nest$ind.plot
plot_grid(cav.nest$jag.plot$s.plot, cav.nest$jag.plot$w.plot, cav.nest$ind.plot,
             ncol = 1,
             nrow = 3, rel_heights = c(1/2, 1/2, 1/4),
             left = "Log Population",
             bottom = "Year")
#-----------------------------------------------------------------------------
# diverse nest type
div.nest <- nestSite_trends(robertsDB, bp_counts, 'Diverse', 'Diverse Nest Type')

div.nest$ind.plot
grid.arrange(div.nest$jag.plot$s.plot, div.nest$jag.plot$w.plot,
             nrow = 2, heights = c(1/2, 1/2),
             left = "Log Population",
             bottom = "Year")
#-----------------------------------------------------------------------------
# grass nest type
grass.nest <- nestSite_trends(robertsDB, bp_counts, 'Grass', 'Grass Nest Type')

grass.nest$ind.plot
grid.arrange(grass.nest$jag.plot$s.plot, grass.nest$jag.plot$w.plot,
             nrow = 2, heights = c(1/2, 1/2),
             left = "Log Population",
             bottom = "Year")
#-----------------------------------------------------------------------------
# ground nest type
ground.nest <- nestSite_trends(robertsDB, bp_counts, 'Ground', 'Ground Nest Type')

ground.nest$ind.plot
grid.arrange(ground.nest$jag.plot$s.plot, ground.nest$jag.plot$w.plot,
             nrow = 2, heights = c(1/2, 1/2),
             left = "Log Population",
             bottom = "Year")
#-----------------------------------------------------------------------------
# reed nest type
reed.nest <- nestSite_trends(robertsDB, bp_counts, 'Reed', 'Reed Nest Type')

reed.nest$lpi
grid.arrange(reed.nest$jag.plot$s.plot, reed.nest$jag.plot$w.plot,
             nrow = 2, heights = c(1/2, 1/2),
             left = "Log Population",
             bottom = "Year")
#-----------------------------------------------------------------------------
# shrub nest type
shrub.nest <- nestSite_trends(robertsDB, bp_counts, 'Shrub', 'Shrub Nest Type')

shrub.nest$lpi
grid.arrange(shrub.nest$jag.plot$s.plot, shrub.nest$jag.plot$w.plot,
             nrow = 2, heights = c(1/2, 1/2),
             left = "Log Population",
             bottom = "Year")
#-----------------------------------------------------------------------------
# tree nest type
tree.nest <- nestSite_trends(robertsDB, bp_counts, 'Tree', 'Tree Nest Type')

tree.nest$lpi
grid.arrange(tree.nest$jag.plot$s.plot, tree.nest$jag.plot$w.plot,
             nrow = 2, heights = c(1/2, 1/2),
             left = "Log Population",
             bottom = "Year")

##### Species trends based on foraging substratum #####

robertsDB %>% group_by(numForage) %>% count()

foraging_trends <- function(traitsdf, countsdf, numForagSub, title){
  
  # traitsdf <- robertsDB
  # countsdf <- bp_counts
  # numForagSub <- 1

  # extract species with certain trait
  Ids <- traitsdf[which(traitsdf$numForage==numForagSub),]$SppRef
  t <- bp_counts[,Ids]
  
  # combine the counts
  combCounts <- apply(t, 1, sum, na.rm = T)
  combCounts <- as.data.frame(cbind('counts' = combCounts,
                                    'Season' = countsdf$Season,
                                    'Year' = countsdf$Year))
  
  combCounts[combCounts==0] <- NA
  
  # run jags on them
  jagmod <- jags_analysis(combCounts)
  
  # plot them
  jPlot <- ts_jag_plot(jagmod, combCounts, title)
  
  # run jags on each species in trait group
  trait_species <- list()
  
  for (i in 1:length(t)){
    trait_species[[i]] <- jags_analysis(as.data.frame(cbind(t[,i],
                                                            'Season' = countsdf$Season,
                                                            'Year' = countsdf$Year)))
  }
  
  # extract beta values
  betas <- lapply(trait_species, function(x){
    return(x$mean$beta)
  })
  
  betas <- as.data.frame(betas)
  colnames(betas) <- colnames(t)
  
  # get average trend per year
  betas.avg <- apply(betas, 1, mean, rm.na = T)
  
  # calc LPI
  lpi <- c(1)
  for(i in 2:length(betas.avg)){
    lpi <-  c(lpi, (lpi[i-1] * exp(betas.avg[i])))
    
  }
  
  return(list("jag.plot" = jPlot,
              "lpi"=lpi))
}

one.type <- foraging_trends(robertsDB, bp_counts, 1, "One Type of Foraging Substratum")

one.type$lpi
grid.arrange(one.type$jag.plot$s.plot, one.type$jag.plot$w.plot,
             nrow = 2, heights = c(1/2, 1/2),
             left = "Log Population",
             bottom = "Year")


two.types <- foraging_trends(robertsDB, bp_counts, 2, "Two Type of Foraging Substratum")

two.types$lpi
grid.arrange(two.types$jag.plot$s.plot, two.types$jag.plot$w.plot,
             nrow = 2, heights = c(1/2, 1/2),
             left = "Log Population",
             bottom = "Year")

three.types <- foraging_trends(robertsDB, bp_counts, 3, "Three Type of Foraging Substratum")

three.types$lpi
grid.arrange(three.types$jag.plot$s.plot, three.types$jag.plot$w.plot,
             nrow = 2, heights = c(1/2, 1/2),
             left = "Log Population",
             bottom = "Year")


# plotting individual models
#for(i in length(aq.nest)){
aplot <-  ts_jag_plot(aq.nest$jag.mods[[6]], 
              as.data.frame(cbind(t[,i],
                                  'Season' = bp_counts$Season,
                                  'Year' = bp_counts$Year)),
              title = "a title")
#}





















