library(spatialreg)
library(stringi)

us_neighboring_states = read.csv(file = "D:/Term5/IS&RA/US_state_neighbors.csv",stringsAsFactors = FALSE)
#Check if all the states match in our df and us_neighboring_states
sum = 0
for (l in (unique(df$state))) {
  sum = sum + (l %in% us_neiboring_states[,1])
}
#Removing space from neighbor column
for (l in (1:length(us_neighboring_states[,1]))) {
  us_neighboring_states[l,2] = gsub(", ",",", us_neighboring_states[l,2])
}
#Getting weight matrix
weight_mat = matrix(0, nrow = n_sp, ncol = n_sp)
for (loc in unique(df$state)) {
  loc_idx = which(unique(df$state) == loc)
  for (neighbor in strsplit(us_neighboring_states[which(us_neighboring_states[,1] == loc) ,2], ",")[[1]]) {
    neighbor_idx = which(unique(df$state) == neighbor)
    weight_mat[loc_idx, neighbor_idx] = 1
  }
}
weight_mat_tmp = matrix(1, nrow = n_tmp_full_training2, ncol = n_tmp_full_training2)

weight_mat_fin = weight_mat %x% weight_mat_tmp

weight_mat_fin = weight_mat_fin/(sum(weight_mat_fin))

slx_model = lmSLX(100*log_prevalence~ log(population) + time + time_sq + prev_log_prevalence,data =  full_training2,
                  listw =spdep::nb2listw(weight_mat_fin,style="W") )
