library(readr)
library(invgamma)
library(goftest)
library(matlib)
library(geosphere)
library(ape)
library(ggplot2)

state_lat_long = unique(df[,c(2,3,4)])

weight_mat=as.matrix(distm(cbind(state_lat_long$long,state_lat_long$lat), fun = distHaversine))
rho =0.001   #converting to Km
weight_decay_mat=exp(-rho*weight_mat)
diag(weight_decay_mat)=0
moran_i_vec = vector()
for (tmp in seq.Date(min(df$week),max(df$week), 7)) {
  moran_i_vec = c(moran_i_vec, Moran.I(df %>% filter(week == (tmp)) %>% pull(log_prevalence),weight_decay_mat)$observed)
}


