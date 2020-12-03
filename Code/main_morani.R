
# load relevant libraries
library(readr)
library(invgamma)
library(goftest)
library(matlib)
library(geosphere)
library(ape)
library(ggplot2)

# run the main_data_prep script to get df
View(df)

# find the latitude longitude and prepare weight matrix
state_lat_long = unique(df[,c("state","lat","long")])
weight_mat = as.matrix(distm(cbind(state_lat_long$long,state_lat_long$lat), fun = distHaversine))

# find exponential decay matrix
rho = 0.001   # needed converting to Km
weight_decay_mat = exp(-rho*weight_mat)
diag(weight_decay_mat) = 0

# find the Moran's I index for all weeks
all_weeks = unique(df$week)
moran_i_mat = matrix(nrow = length(all_weeks),ncol = 4)
for (tmp in 1:length(all_weeks)){
  moran_test = Moran.I(df %>% filter(week == all_weeks[tmp]) %>% pull(log_prevalence),weight_decay_mat)
  moran_i_mat[tmp,] = unlist(moran_test)
}

# combine all results
moran_result = data.frame(all_weeks,moran_i_mat)
colnames(moran_result) = c("week","observed","expected","sd","pvalue")
moran_result$significant = ifelse(moran_result$pvalue < 0.05,"significant","not significant")

# see the plot
moran_result %>%
  ggplot(aes(x = week,y = observed)) +
  geom_line() +
  geom_point(aes(x = week,y = observed,pch = significant),size = 3) +
  scale_shape_manual(values = c("significant" = 19,"not significant" = 1)) +
  ylab("Moran's I statistic") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_blank()) + 
  ggsave(filename = "~/Documents/GitHub/COVID_Analysis/Code/morans_I.eps",
         width = 300,
         height = 200,
         units = "mm",
         dpi = 800)
