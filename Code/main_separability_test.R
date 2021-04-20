
# load relevant libraries
library(ggplot2)
library(ggpubr)
library(usmap)
library(sp)
library(spacetime)
library(gstat)
library(covatest)
library(IRISSeismic)
library(tidyverse)
library(psd)

# read data of the residuals
df = read.csv("~/Documents/GitHub/COVID_Analysis/Code/covid_withresidual.csv")
df$response = df$resid_log_prev

# calculate cross spectrum for all pairs of locations
locs = as.character(unique(df$state))
pair_cs = data.frame(expand.grid(locs,locs))
colnames(pair_cs) = c("loc1","loc2")
pair_cs = pair_cs %>% filter(loc1 != loc2) %>% mutate_if(is.factor,as.character)
all_cse = list()
for (i in 1:nrow(pair_cs)){
  s1 = as.numeric(subset(df,state == pair_cs$loc1[i])$response)
  x1 = ts((s1 - mean(s1)))
  s2 = as.numeric(subset(df,state == pair_cs$loc2[i])$response)
  x2 = ts((s2 - mean(s2)))
  all_cse[[i]] = crossSpectrum(x = ts.union(x1,x2),kernel("daniell",floor(length(s1)^0.2)),demean = T,detrend = T)
}

# calculate coherence values
coh = do.call(cbind,lapply(all_cse,function(x) atanh(x$coh)))

# conduct an anova like test many times
rpts = 10000
pvals_sep = numeric(rpts)
pvals_stn = numeric(rpts)

for (j in 1:rpts){
  
  # take a random sample of 24 location pairs, which are different from each other
  locidx1 = sample(locs,24)
  locidx2 = setdiff(locs,locidx1)[1:24]
  subsetlocs = cbind(locidx1,locidx2)
  subsetids = unlist(apply(subsetlocs,1,function(x) which(pair_cs$loc1 == x[1] & pair_cs$loc2 == x[2])))
  
  # subset out coherence values based on frequency amd the above random sample
  freq_idx = c(1,6,11,16,21)
  ival = 24
  jval = length(freq_idx)
  subsetdf = data.frame(locs = rep(paste("pair",c(1:ival),sep = ""),each = jval),
                        times = rep(paste("freq",c(1:jval),sep = ""),ival),
                        coherence = c(coh[freq_idx,subsetids]))
  
  # perform the anova tests to find out stationarity and separability
  h1model = lm(coherence ~ locs + times,data = subsetdf)
  h0model_loc = lm(coherence ~ locs,data = subsetdf)
  h0model_time = lm(coherence ~ times,data = subsetdf)
  anova1 = anova(h0model_loc,h1model)
  anova2 = anova(h0model_time,h1model)
  pvals_sep[j] = anova1$`Pr(>F)`[2]
  pvals_stn[j] = anova2$`Pr(>F)`[2]
}

# see how many times the test is not rejecting the separability assumption
mean(pvals_sep > 0.05)

# see how many times the test is not rejecting the spatial stationarity assumption
mean(pvals_stn > 0.05)



#:::::::::::::: PREVIOUS CODES (ROUGH WORK)


ival = ncol(coh)
jval = 6#nrow(coh)
anovadf = data.frame(locs = rep(paste("pair",c(1:ival),sep = ""),each = jval),
                     times = rep(paste("freq",c(1:jval),sep = ""),ival),
                     coherence = c(coh[c(1,5,9,13,17,21),]))

# test using anova
subsetdf = anovadf
h1model = lm(coherence ~ locs + times,data = subsetdf)
h0model_loc = lm(coherence ~ locs,data = subsetdf)
h0model_time = lm(coherence ~ times,data = subsetdf)

anova(h0model_loc,h1model)
anova(h0model_time,h1model)

#::::: WORKING WITH RESIDUALS

df_res = read.csv("~/Documents/GitHub/COVID_Analysis/Code/covid_withresidual.csv")
df_res$response = log((df_res$prevalence+0.01)/(1-df_res$prevalence))

# calculate cross spectrum for all pairs of locations
locs = as.character(unique(df_res$state))
pair_cs = data.frame(expand.grid(locs,locs))
colnames(pair_cs) = c("loc1","loc2")
pair_cs = pair_cs %>% filter(loc1 != loc2) %>% mutate_if(is.factor,as.character)
pair_cs$distance = NA
all_cse = list()
for (i in 1:nrow(pair_cs)){
  loc1data = subset(df_res,state == pair_cs$loc1[i])
  loc2data = subset(df_res,state == pair_cs$loc2[i])
  coord1 = as.numeric(c(loc1data$lat[1],loc1data$long[1]))
  coord2 = as.numeric(c(loc2data$lat[1],loc2data$long[1]))
  pair_cs$distance[i] = pracma::haversine(coord1,coord2)
  s1 = as.numeric(loc1data$response)
  x1 = ts((s1 - mean(s1)))
  s2 = as.numeric(loc2data$response)
  x2 = ts((s2 - mean(s2)))
  # tmp = pspectrum(cbind(x1,x2),detrend = T,verbose = F)
  all_cse[[i]] = crossSpectrum(x = ts.union(x1,x2),kernel("daniell",floor(length(s1)^0.2)),demean = T,detrend = T)
}

# subset out based on distance and time
dist_idx = which(pair_cs$distance > 500)
time_idx = c(1,6,11,16,21)

# calculate coherence values
coh = do.call(cbind,lapply(all_cse[dist_idx],function(x) atanh(x$coh)))

# make an anova like structure
ival = ncol(coh)
jval = length(time_idx)
anovadf = data.frame(locs = rep(paste("pair",c(1:ival),sep = ""),each = jval),
                     times = rep(paste("freq",c(1:jval),sep = ""),ival),
                     coherence = c(coh[time_idx,]))

# test using anova
subsetdf = anovadf
h1model = lm(coherence ~ locs + times,data = subsetdf)
h0model_loc = lm(coherence ~ locs,data = subsetdf)
h0model_time = lm(coherence ~ times,data = subsetdf)

anova1 = anova(h0model_loc,h1model)
anova(h0model_time,h1model)


#::::: SOME ADDITIONAL EXPERIMENTATION

pvals_sep = numeric(3000)
coh = do.call(cbind,lapply(all_cse,function(x) atanh(x$coh)))
for (j in 1:3000){
  locidx1 = sample(locs,24)
  locidx2 = setdiff(locs,locidx1)[1:24]
  subsetlocs = cbind(locidx1,locidx2)
  subsetids = unlist(apply(subsetlocs,1,function(x) which(pair_cs$loc1 == x[1] & pair_cs$loc2 == x[2])))
  
  ival = 24
  subsetdf = data.frame(locs = rep(paste("pair",c(1:ival),sep = ""),each = jval),
                        times = rep(paste("freq",c(1:jval),sep = ""),ival),
                        coherence = c(coh[time_idx,subsetids]))
  
  h1model = lm(coherence ~ locs + times,data = subsetdf)
  h0model_loc = lm(coherence ~ locs,data = subsetdf)
  anova1 = anova(h0model_loc,h1model)
  pvals_sep[j] = anova1$`Pr(>F)`[2]
}
mean(pvals_sep > 0.05)

#::::::::::: ANOTHER METHOD: TRYING IT OUT

# spatial points object
locs = unique(df[,c("state","long","lat")])
sp = SpatialPoints(coords = locs[,2:3])
time = as.POSIXct(sort(unique(df$week)))
mydata = data.frame(values = df[order(df$week,df$state),"log_prevalence"])
stfdf_obj = STFDF(sp,time,mydata)

# find spatial distances
spdist = spDists(sp)
sorted = sort.int(spdist,index.return = T)
minidx = sorted$ix[which(sorted$x > 0)]
rowno = minidx %% 49; rowno[rowno == 0] = 49
colno = ceiling(minidx/49)

nn = length(sorted$x[sorted$x > 0])/2
choosen_locs = data.frame(state1 = rep(NA,nn),state2 = rep(NA,nn),dist = rep(0,nn))
for (i in 1:nn){
  choosen_locs[i,1] = as.character(locs$state[rowno[(2*i-1)]])
  choosen_locs[i,2] = as.character(locs$state[rowno[(2*i)]])
  choosen_locs[i,3] = spdist[rowno[(2*i-1)],rowno[(2*i)]]
}


states_mindist = c(unlist(t(choosen_locs[c(1,2,7,9,15,16),-3])))
locids = as.numeric(sapply(states_mindist,function(x) which(locs$state == x)))
sel.staz.sym = as.character(locids)
sp.couples.in.sym <- matrix(data = sel.staz.sym, ncol = 2, byrow = TRUE)
t.couples.in.sym <- c(1,2)
couples.sym <- couples(sel.staz = sel.staz.sym, 
                       sp.couples.in = sp.couples.in.sym, 
                       t.couples.in = t.couples.in.sym,
                       typetest = "sym", 
                       typecode = character())
rr_13 = stfdf_obj[sel.staz.sym,]
block.sym <- blocks(lb = 40,ls = 10,matdata = rr_13,pardata1 = 1,pardata2 = 1,stpairs = couples.sym)
