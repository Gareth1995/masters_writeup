source("R/cwac_functions.r")
library(knitr)
library(vegan)
library(tidyverse)
library(cowplot)
library(bayesplot)
#
# load relevant data
load('data/robertsDB.RData')
load('data/barberspan.counts.RData')

birdcnts <- bp_counts[,1:103]

# separate into resident and migrant
resid <- robertsDB %>%
  filter(Migrant=="n") %>%
  select(SppRef)

migid <- robertsDB %>%
  filter(Migrant=="y") %>%
  select(SppRef)

residents <- birdcnts[,resid$SppRef]
migrants <- birdcnts[,migid$SppRef]

# Applying the jag model on a random bird
bird_df <- as.data.frame(bp_counts[c('89', 'Season', 'Year')])

bird.jags <- jags_analysis(bird_df)
birdtau <- ts_jag_plot(bird.jags, bird_df, "Egyptian Goose")

sims <- bird.jags$sims.list$mu_t

# using simulation sample to create yearly change plot
for(i in 1:ncol(sims)){
  perc_change <- data.frame(lower = 0, median = 0, upper = 0)
  
  for(x in 2:ncol(sims)){
    # using median and lower and upper quantile
    perc_change[x,] = quantile(((sims[,x] - sims[,x-1]) / (sims[,x-1]))*100)[c(2,3,4)]
  }
}

perc_change <- cbind(perc_change, "Year" = unique(bp_counts$Year))

ggplot(as.data.frame(perc_change), aes(x = Year, group = 1)) +
  
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "gray80") +
  geom_line(aes(y = median, color = "grey1"), lwd = 1, lty = 2)

# finding yearly change per year for migrant group
# extract species with certain trait
traitsdf <- robertsDB
countsdf <- bp_counts

Ids <- traitsdf[which(traitsdf$Migrant=='y'),]$SppRef
t <- countsdf[,Ids]

# run jags on each species in trait group
trait_species <- list()
# 
for (i in 1:length(t)){
  trait_species[[i]] <- jags_analysis(as.data.frame(cbind(t[,i],
                                                          'Season' = countsdf$Season,
                                                          'Year' = countsdf$Year)))
}

# sum together all sim lists
migtot <- data.frame(matrix(0, 9000, 26))
for(i in 1:26){
  
  migtot <- migtot + trait_species[[i]]$sims.list$mu_t
  
}

# loop through sims list to find lower median and upper quantile
perc_change <- data.frame(lower = 0, med = 0, upper = 0)
for(x in 2:26){
  # using median and lower and upper quantile
  perc_change[x,] = quantile(((migtot[,x] - migtot[,x-1]) / (migtot[,x-1]))*100)[c(2,3,4)]
}
perc_change <- cbind(perc_change, "Year" = unique(bp_counts$Year))

# plot the yearly percentage change
yrlychange.plot <- ggplot(perc_change, aes(x = Year)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "gray80") +
  geom_line(aes(y = med, color = "red"), lwd = 1) +
  geom_hline(yintercept=0, linetype = 'dashed') +
  scale_x_continuous(breaks = seq(min(perc_change$Year), max(perc_change$Year), by = 1)) +
  labs(title = title, y = "Average yearly percentage change", x = "") +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = 'none')


#### (Q11) Displaying distributions for variance parameters #####

hist(bird.jags$sims.list$zeta)# showing variance of error term zeta
hist(bird.jags$sims.list$w)   # showing variance of error term w
hist(bird.jags$sims.list$tau.alpha, breaks = 100) # showing variance of error term eps
hist((bird.jags$sims.list$sig.eps)**-1, breaks = 100)

hist(bird.jags$sims.list$mu_t[,2]) # showing distribution of a summer count
hist(bird.jags$sims.list$mu_t[,23])
hist(bird.jags$sims.list$mu_wt[,2]) # showing distribution of a winter count
hist(bird.jags$sims.list$mu_wt[,23])

# (Q12) MCMC output (convergence plots)
traceplot(bird.jags) # to show convergence of estimated counts


#### (Q8) applying shannon and simpson index to sims lists ####

# run jag analysis on resident birds
res.jags <- list()

for (i in 1:length(residents)){
  res.jags[[i]] <- jags_analysis(as.data.frame(cbind(residents[,i],
                                          'Season' = bp_counts$Season,
                                          'Year' = bp_counts$Year)))
}

# run jag analysis on migrant birds
mig.jags <- list()

for (i in 1:length(migrants)){
  mig.jags[[i]] <- jags_analysis(as.data.frame(cbind(migrants[,i],
                                          'Season' = bp_counts$Season,
                                          'startDate' = bp_counts$Year)))
}

# iterate through migrant group and extract year i
shan.mig.sum <- data.frame('lower' = 0, 'med' = 0, 'upper' = 0)

# place all year i's together in one df
for(i in 1:26){
  ayear <- lapply(mig.jags, function(x){
    x$sims.list$mu_t[,i]
  })
  
  # run the diversity() function on the year i df
  # save upper, lower and median to a separate df
  shan.mig.sum <- rbind(shan.mig.sum, quantile(diversity(exp(as.data.frame(ayear)), index = "shannon"),
                                probs = c(0.025,0.5,0.975)))
}

shan.mig.sum <- cbind(exp(shan.mig.sum[2:nrow(shan.mig.sum),]),
                      'Year' = unique(bp_counts$Year))

ggplot(shan.mig.sum, aes(x = Year)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "gray80") +
  geom_line(aes(y = med, color = 'red')) 


#### Q11 plotting distribution of variance parameters ####

# use a common bird with a fair amount of data (Egyptian goose)
egoose <- as.data.frame(bp_counts[c('89', 'Season', 'Year')])
egoose.jags <- jags_analysis(egoose)
egoose.plot <- ts_jag_plot(egoose.jags, egoose, "Egyptian Goose")

hist(egoose.jags$sims.list$sig.alpha, breaks = 100, main = "Sigma_alpha for Egyptian Goose", xlab = "variance") #sigma_alpha
hist(egoose.jags$sims.list$sig.e, breaks = 100, main = "Sigma_e for Egyptian Goose", xlab = "variance") # sigma_e

hist(egoose.jags$sims.list$sig.w2, breaks = 100, main = "Sigma_w for Egyptian Goose", xlab = "variance") # sigma_w
hist(egoose.jags$sims.list$sig.zeta, breaks = 100, main = "Sigma_zeta for Egyptian Goose", xlab = "variance") #sigma_zeta
hist(egoose.jags$sims.list$sig.eps, breaks = 100, main = "Sigma_epsilon for Egyptian Goose", xlab = "variance") # sigma_eps

hist(egoose.jags$sims.list$w, breaks = 100, main = "w for Egyptian Goose", xlab = "variance") # sigma_w
hist(egoose.jags$sims.list$zeta, breaks = 100, main = "zeta for Egyptian Goose", xlab = "variance") #sigma_zeta
hist(egoose.jags$sims.list$eps, breaks = 100, main = "epsilon for Egyptian Goose", xlab = "variance") # sigma_eps

# checking traceplots
traceplot(egoose.jags, "sigma_e")
traceplot(egoose.jags, "sigma_alpha")
traceplot(egoose.jags, "sigma_w")
traceplot(egoose.jags, "sigma_zeta")
traceplot(egoose.jags, "sigma_epsilon")

# first 5 summer counts and first 5 winter counts
traceplot(egoose.jags, "mu_t")
traceplot(egoose.jags, "mu_wt")
#coda::traceplot(coda::as.mcmc(egoose.jags$sims.list$mu_t), col = 1:3)

color_scheme_set("blue")
mcmc_trace(as.array(egoose.jags), pars = c("sig.e", "sig.alpha"))

# and use a rare bird

# convergence for 13 or more counts
eurasian_curlew <- as.data.frame(bp_counts[c('263', 'Season', 'Year')]) # 13 observations
curlew.jags <- jags_analysis(eurasian_curlew)


curlew.plot <- ts_jag_plot(curlew.jags, eurasian_curlew, "Eurasian Curlew")

hist(curlew.jags$sims.list$sig.alpha, breaks = 100, main = "Sigma_alpha for Eurasian Curlew", xlab = "variance") #sigma_alpha
hist(curlew.jags$sims.list$sig.e, breaks = 100, main = "Sigma_e for Eurasian Curlew", xlab = "variance") # sigma_e
hist(curlew.jags$sims.list$sig.w2, breaks = 100, main = "Sigma_w for Eurasian Curlew", xlab = "variance") # sigma_w
hist(curlew.jags$sims.list$sig.zeta, breaks = 100, main = "Sigma_zeta for Eurasian Curlew", xlab = "variance") #sigma_zeta
hist(curlew.jags$sims.list$sig.eps, breaks = 100, main = "Sigma_epsilon for Eurasian Curlew", xlab = "variance") # sigma_eps

hist(curlew.jags$sims.list$w, breaks = 100, main = "w for Eurasian Curlew", xlab = "variance") # sigma_w
hist(curlew.jags$sims.list$zeta, breaks = 100, main = "zeta for Eurasian Curlew", xlab = "variance") #sigma_zeta
hist(curlew.jags$sims.list$eps, breaks = 100, main = "epsilon for Eurasian Curlew", xlab = "variance") # sigma_eps

# checking traceplots
traceplot(curlew.jags, "sig.e")
traceplot(curlew.jags, "sig.alpha")
traceplot(curlew.jags, "sig.w2")
traceplot(curlew.jags, "sig.zeta")
traceplot(curlew.jags, "sig.eps")

# do analysis on species with more than 6 counts only
barbernas <- apply(birdcnts, 2, function(x){
  return(x==0)
})
barbernas

barber0s = apply(barbernas, 2, sum, na.rm = T)
less13 = names(barber0s[barber0s>19])

df <- robertsDB %>% select(SppRef, `Common name`, Migrant)
newdf <- df[df$SppRef %in% less13,]

bp_counts <- bp_counts[,barbernas>6]


