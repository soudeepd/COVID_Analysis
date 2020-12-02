library(dplyr)
library(utils)
library(readr)
library(invgamma)
library(goftest)
library(matlib)
library(geosphere)
library(ape)
library(zoo)
library(matrixcalc)
library(MASS)
library(mnormt)
library(matrixStats)

#estimating with the first training data
unique_state_lat_long = unique(training1 %>% dplyr::select(state,lat,long))
distance_mat_training1 = as.matrix(distm(cbind(unique_state_lat_long$long,unique_state_lat_long$lat), fun = distHaversine))
#training 1 values
y_training1 = 100*training1$log_prevalence
n_training1 = length(y_training1)
n_sp_training1 = length(unique_state_lat_long$state)
n_tmp_training1 = n_training1/n_sp_training1
#total dataset values
n_total = length(df$state)
n_sp_total = length(unique(df$state))
n_tmp_total = n_total/ n_sp_total
#validation set values
y_validation1 = 100*validation1$log_prevalence
n_validation1 = length(y_validation1)
n_sp_validation1 = length(unique(validation1$state))
n_tmp_validation1 = n_validation1/n_sp_validation1
#testing values
y_testing = 100*testing$log_prevalence
n_testing = length(y_testing)
n_sp_testing = length(unique(testing$state))
n_tmp_testing = n_testing/ n_sp_testing
p=6
#X matrix for the model
X_mat_training1 = cbind(rep(1,n_training1), log(training1$population), training1$time, training1$time_sq, 
                        100*(training1$prev_log_prevalence), training1$prev_log_death,#training1$prev_log_death
                        100*(training1$prev_log_prevalence_country))

#Getting week difference mat for Sigma_T
week_vec = 1:n_tmp_training1
week_diff_mat = as.matrix(dist(week_vec,diag  =  TRUE,upper  =  TRUE))

#Error Vectors
mape_vec_training_training1 = vector()
mape_vec_validation_training1 = vector() 

burn_period_training1 = 1000
sample_size_training1 = 500
diff_in_random_draws_training1 = 20

#Getting sample means for every rho
v_sample_mean_training1 = vector()
theta_sample_mean_training1 = vector()
sigma1_sq_sample_mean_training1 = vector()
sigma2_sq_sample_mean_training1 = vector()
theta_significance_training1 = vector()
#constant parameters values
a_training1=2
lambda_training1=1
#Saving sample of posterior values for each rho values
v_mat_sample_list_training1 = list()
theta_sample_list_training1 = list()
sigma1_sq_sample_vec_training1 = vector()
sigma2_sq_sample_vec_training1 = vector()

rho_vec_tmp_training1 = c(0.10, 0.50, 0.75, 1, 1.5, 3)#seq(0.1, 0.5, 0.1)
rho_vec_sp_training1 = c(0.0001, 0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1)#seq(0.5, 1, 0.1)
rho_idx = 1

for(rho_tmp in rho_vec_tmp_training1){
  for (rho_sp in rho_vec_sp_training1) {
    
    SIGMA_sp_training1 = exp(-rho_sp*distance_mat_training1/(1000)) # dividing by 80000 to get matrix values similar to temporal values
    SIGMA_tmp_training1 = exp(-rho_tmp*week_diff_mat)
    inv_SIGMA_sp_training1 = chol2inv(chol(SIGMA_sp_training1))
    inv_SIGMA_tmp_training1 = chol2inv(chol(SIGMA_tmp_training1))
    inv_SIGMA_sp_kro_inv_SIGMA_tmp_training1 = inv_SIGMA_sp_training1 %x% inv_SIGMA_tmp_training1
    #Giving values to different parameters of our model
    v_mat_training1 = matrix(0,nrow = n_tmp_training1,ncol = n_sp_training1)
    theta_training1 = rep(0,p+1)
    sigma1_sq_training1 = 1
    sigma2_sq_training1 = 1
    
    #time taken analysis for parameters
    total_tt_v = 0
    total_tt_theta = 0
    total_tt_sigma1_sq = 0
    total_tt_sigma2_sq = 0
    
    #Getting inv Sigma_sp_22 to use for every iteration
    inv_SIGMA_sp_trans_22_matlist_training1=list()
    SIGMA_sp_trans_21_mat_training1 = vector()
    for (j in (1:n_sp_training1)) {
      idx = 1:n_sp_training1
      
      if(j != 1){
        idx[1]=j
        idx[j]=1
      }
      SIGMA_sp_trans_training1 = diag(n_sp_training1)[idx,] %*% SIGMA_sp_training1 %*% diag(n_sp_training1)[idx,]
      inv_SIGMA_sp_trans_training1 = diag(n_sp_training1)[idx,] %*% inv_SIGMA_sp_training1 %*% diag(n_sp_training1)[idx,]
      inv_SIGMA_sp_trans_a_training1 = inv_SIGMA_sp_trans_training1[1,1]
      inv_SIGMA_sp_trans_b_training1 = inv_SIGMA_sp_trans_training1[1,-1]
      inv_SIGMA_sp_trans_c_training1 = inv_SIGMA_sp_trans_training1[-1,1]
      inv_SIGMA_sp_trans_d_training1 = inv_SIGMA_sp_trans_training1[-1,-1]
      
      inv_SIGMA_sp_trans_22_training1 = inv_SIGMA_sp_trans_d_training1-(inv_SIGMA_sp_trans_c_training1%*%t(inv_SIGMA_sp_trans_b_training1))/inv_SIGMA_sp_trans_a_training1
      SIGMA_sp_trans_21_training1 = SIGMA_sp_trans_training1[-1,1]
      
      inv_SIGMA_sp_trans_22_matlist_training1[[j]] = inv_SIGMA_sp_trans_22_training1
      SIGMA_sp_trans_21_mat_training1 = cbind(SIGMA_sp_trans_21_mat_training1,SIGMA_sp_trans_21_training1)
      
    }
    
    theta_sample_training1=vector()
    sigma1_sq_sample_training1=vector()
    sigma2_sq_sample_training1=vector()
    v_sample_training1=vector()
    for (i in (1:(burn_period_training1+sample_size_training1*diff_in_random_draws_training1))) {
      #v posterior
      st_v = Sys.time()
      
      for (j in (1:n_sp_training1)) {
        idx=1:n_sp_training1
        
        if(j!=1){
          idx[1] = j
          idx[j] = 1
        }
        
        v_cond_training1 = (v_mat_training1 %*% diag(n_sp_training1)[idx,])[(1+n_tmp_training1):(n_training1)]
        
        #SIGMA' conditional without sigma1_sq as it's included in the posterior 
        SIGMA_cond_training1 = (1 -t(SIGMA_sp_trans_21_mat_training1[,j])%*%inv_SIGMA_sp_trans_22_matlist_training1[[j]] %*%
                                  SIGMA_sp_trans_21_mat_training1[,j])%x% SIGMA_tmp_training1
        mu_cond_training1 = ((t(SIGMA_sp_trans_21_mat_training1[,j])%*%inv_SIGMA_sp_trans_22_matlist_training1[[j]])%x% 
                               diag(n_tmp_training1))%*%v_cond_training1
        inv_SIGMA_cond_training1 = chol2inv(chol(SIGMA_cond_training1))
        
        v_covar_training1 = chol2inv(chol(inv_SIGMA_cond_training1 /sigma1_sq_training1 + diag(n_tmp_training1)/sigma2_sq_training1))
        v_mean_training1  = v_covar_training1 %*% (inv_SIGMA_cond_training1 %*% 
                                                     mu_cond_training1/sigma1_sq_training1 + 
                                                     (y_training1[(1+(j-1)*n_tmp_training1):(j*n_tmp_training1)] -
                                                        X_mat_training1[(1+(j-1)*n_tmp_training1):(j*n_tmp_training1),] %*% 
                                                        theta_training1)/sigma2_sq_training1 )
        v_mat_training1[,j] = v_mean_training1 + t(chol(v_covar_training1))  %*% rnorm(n_tmp_training1)
        j=j+1
      }
      
      total_tt_v = total_tt_v + as.numeric(difftime(Sys.time(), st_v), units="secs")
      
      #theta posterior
      st_theta = Sys.time()
      
      theta_covar_training1 = solve( (t(X_mat_training1) %*% X_mat_training1)/sigma2_sq_training1 + diag(p+1) ) 
      theta_mean_training1 = theta_covar_training1 %*% ( t(X_mat_training1)%*% 
                                                           (y_training1 - v_mat_training1[1:n_training1]))/sigma2_sq_training1
      theta_training1 = theta_mean_training1 + t(chol(theta_covar_training1)) %*%rnorm(p+1)
      
      total_tt_theta = total_tt_theta + as.numeric(difftime( Sys.time(), st_theta), units="secs")
      #Sigma1^2 posterior
      st_sigma1_sq = Sys.time()
      
      sigma1_sq_a_training1 = a_training1 + n_training1/2
      sigma1_sq_lambda_training1 = (t(v_mat_training1[1:n_training1]) %*% inv_SIGMA_sp_kro_inv_SIGMA_tmp_training1 
                                    %*% v_mat_training1[1:n_training1])/2 + lambda_training1
      sigma1_sq_training1 = invgamma::rinvgamma(1,sigma1_sq_a_training1,rate = sigma1_sq_lambda_training1)
      
      total_tt_sigma1_sq = total_tt_sigma1_sq + as.numeric(difftime(Sys.time(), st_sigma1_sq), units="secs")
      #Sigma2^2 posterior
      st_sigma2_sq = Sys.time()
      
      sigma2_sq_a_training1 = a_training1 + n_training1/2
      temp_vec = y_training1 -X_mat_training1%*%theta_training1 - v_mat_training1[1:n_training1]
      sigma2_sq_lambda_training1 = (t(temp_vec) %*% (temp_vec))/2 + lambda_training1
      sigma2_sq_training1 = invgamma::rinvgamma(1,sigma2_sq_a_training1,rate = sigma2_sq_lambda_training1)
      
      total_tt_sigma2_sq = total_tt_sigma2_sq + as.numeric(difftime(Sys.time(), st_sigma2_sq), units="secs")
      
      #Collecting samples for prediction
      if(i> burn_period_training1 & i%%diff_in_random_draws_training1 ==0){
        v_sample_training1 = cbind(v_sample_training1, v_mat_training1[1:n_training1])
        theta_sample_training1 = cbind(theta_sample_training1,theta_training1)
        sigma1_sq_sample_training1 = c(sigma1_sq_sample_training1,sigma1_sq_training1)
        sigma2_sq_sample_training1 = c(sigma2_sq_sample_training1,sigma2_sq_training1)
      }
    }
    #time analysis
    avg_tt_v = total_tt_v / (burn_period_training1+sample_size_training1*diff_in_random_draws_training1)
    avg_tt_theta = total_tt_theta / (burn_period_training1+sample_size_training1*diff_in_random_draws_training1)
    avg_tt_sigma1_sq = total_tt_sigma1_sq / (burn_period_training1+sample_size_training1*diff_in_random_draws_training1)
    avg_tt_sigma2_sq = total_tt_sigma2_sq / (burn_period_training1+sample_size_training1*diff_in_random_draws_training1)
    #estimation of y with the posterior parameters
    estimated_y_training1 = X_mat_training1 %*% rowMeans(theta_sample_training1) + rowMeans(v_sample_training1) + 
      rnorm(n_training1,mean = 0, sd= sqrt(mean(sigma2_sq_sample_training1)))
    mape_training1 = (sum(abs( (y_training1 - estimated_y_training1)/y_training1 )))/n_training1
    mape_vec_training_training1 = c(mape_vec_training_training1, mape_training1)
    
    #Testing for significance
    significance_vec=vector()
    for (ts in (1:(p+1))) {
      significance_vec = c(significance_vec, unname(!(0 >= quantile(theta_sample_training1[ts,], 0.025) & 
                                                        0 <= quantile(theta_sample_training1[ts,], 0.975))))
    }
    theta_significance_training1 = cbind(theta_significance_training1, significance_vec)
    
    #Saving samples for each rho
    v_mat_sample_list_training1[[rho_idx]] = v_sample_training1
    theta_sample_list_training1[[rho_idx]] = theta_sample_training1
    sigma1_sq_sample_vec_training1 = cbind(sigma1_sq_sample_vec_training1, sigma1_sq_sample_training1)
    sigma2_sq_sample_vec_training1 = cbind(sigma2_sq_sample_vec_training1, sigma2_sq_sample_training1)
    
    #Saving value for sample means
    v_sample_mean_training1 = cbind(v_sample_mean_training1, rowMeans(v_sample_training1))
    theta_sample_mean_training1 = cbind(theta_sample_mean_training1, rowMeans(theta_sample_training1))
    sigma1_sq_sample_mean_training1 = c(sigma1_sq_sample_mean_training1, mean(sigma1_sq_sample_training1))
    sigma2_sq_sample_mean_training1 = c(sigma2_sq_sample_mean_training1, mean(sigma2_sq_sample_training1))
    
    #Validation
    v_pred_val_vec = vector()
    y_pred_val_vec = vector()
    y_pred_val_mat = vector()
    pop_vec_validation1 = unique(validation1[,c(2,6)])[,2]
    for (pred_i in (1:n_validation1)) {
      Sigma_12_vec = vector()
      for (j in (1:n_training1)) {
        sp_dist = distm(c(training1$long[j],  training1$lat[j]), c(validation1$long[pred_i], validation1$lat[pred_i]), fun = distHaversine)
        tmp_dist = abs((as.numeric(difftime(validation1$week[pred_i] ,  training1$week[j]), units="days"))/7)
        Sigma_12_vec = c(Sigma_12_vec, exp(-rho_sp*  sp_dist/1000) * exp(-rho_tmp * tmp_dist))
      }
      v_tmp_given_v1_covar = 1 - t(Sigma_12_vec) %*% inv_SIGMA_sp_kro_inv_SIGMA_tmp_training1 %*% Sigma_12_vec
      v_tmp_given_v1_mean = t(Sigma_12_vec) %*% inv_SIGMA_sp_kro_inv_SIGMA_tmp_training1 %*% rowMeans(v_sample_training1)
      v_tmp_given_v1 = v_tmp_given_v1_mean + sqrt(mean(sigma1_sq_sample_training1)) * t(chol(v_tmp_given_v1_covar)) * rnorm(1)
      
      
      if(pred_i %% n_sp_validation1 != 0){#which spatial location s=pred_i %% n_sp_validation1
        tmp = floor(pred_i / n_sp_validation1) +1
        sp = pred_i %% n_sp_validation1
        idx = (sp -1)*n_tmp_validation1 + tmp
        if(tmp == 1){#Implies time t=1
          X_vec_tmp = c(1, log(validation1$population[idx]), validation1$time[idx], validation1$time_sq[idx], 
                        100*validation1$prev_log_prevalence[idx], validation1$prev_log_death[idx], 
                        100*(validation1$prev_log_prevalence_country[idx]))
        }else{#Implies t>1
          X_vec_tmp = c(1, log(validation1$population[idx]), validation1$time[idx], validation1$time_sq[idx],
                        y_pred_val_mat[tmp - 1, sp], validation1$prev_log_death[idx], 
                        100*log(sum((exp(y_pred_val_mat[tmp - 1, ]/100) - 0.1)*pop_vec_validation1)/sum(pop_vec_validation1) +0.1))
        }
        
        y_pred_val = t(X_vec_tmp) %*% rowMeans(theta_sample_training1) + v_tmp_given_v1 + 
          sqrt(mean(sigma2_sq_sample_training1))*rnorm(1)
        y_pred_val_vec = c(y_pred_val_vec, y_pred_val)
        
      }else{#location = n
        sp = n_sp_validation1
        tmp = pred_i / n_sp_validation1
        idx = (sp -1)*n_tmp_validation1 + tmp
        if(tmp == 1){#Implies time t=1
          X_vec_tmp = c(1, log(validation1$population[idx]), validation1$time[idx], validation1$time_sq[idx], 
                        100*validation1$prev_log_prevalence[idx], validation1$prev_log_death[idx], 
                        100*(validation1$prev_log_prevalence_country[idx]))
        }else{#Implies t>1
          X_vec_tmp = c(1, log(validation1$population[idx]), validation1$time[idx], validation1$time_sq[idx],
                        y_pred_val_mat[tmp - 1, sp], validation1$prev_log_death[idx], #log(validation1$prev_new_death[idx] + 0.1)
                        100*log(sum((exp(y_pred_val_mat[tmp - 1, ]/100) - 0.1)*pop_vec_validation1)/sum(pop_vec_validation1) +0.1))
        }
        y_pred_val = t(X_vec_tmp) %*% rowMeans(theta_sample_training1) + v_tmp_given_v1 + 
          sqrt(mean(sigma2_sq_sample_training1))*rnorm(1)
        y_pred_val_vec = c(y_pred_val_vec, y_pred_val)
        y_pred_val_mat = rbind(y_pred_val_mat, unname(y_pred_val_vec))
        y_pred_val_vec = vector()
        
      }
      
    }
    mape_validation = (sum(abs( (y_validation1 - y_pred_val_mat[1:n_validation1])/y_validation1 )))/n_validation1
    mape_vec_validation_training1 = c(mape_vec_validation_training1, mape_validation)
    
    rho_idx = rho_idx + 1
  }
  
}  
#Prediction
#training 1 data values

y_full_training1 = 100*full_training1$log_prevalence
n_full_training1 = length(y_full_training1)
n_sp_full_training1 = length(unique(full_training1$state))
n_tmp_full_training1 = n_full_training1/n_sp_full_training1

#X mat full training
X_mat_full_training1 = cbind(rep(1,n_full_training1), log(full_training1$population), full_training1$time, full_training1$time_sq, 
                             100*(full_training1$prev_log_prevalence), full_training1$prev_log_death,
                             100*(full_training1$prev_log_prevalence_country))

#Getting week difference mat for Sigma_T
week_vec = 1:n_tmp_full_training1
week_diff_mat = as.matrix(dist(week_vec,diag  =  TRUE,upper  =  TRUE))

#Getting distance matrix
unique_state_lat_long_full_training1 = unique(full_training1 %>% dplyr::select(state,lat,long))
distance_mat_full_training1 = as.matrix(distm(cbind(unique_state_lat_long_full_training1$long,unique_state_lat_long_full_training1$lat), 
                                              fun = distHaversine))
#Error vectors
mape_vec_training_full_training1 = vector()
mape_vec_validation_full_training1 = vector() 
#pred_error_vec_full_training1 = vector()

burn_period_full_training1 = 1000
sample_size_full_training1 = 500
diff_in_random_draws_full_training1 = 20
#Getting sample means for every rho
v_sample_mean_full_training1 = vector()
theta_sample_mean_full_training1 = vector()
theta_significance_full_training1 = vector()
sigma1_sq_sample_mean_full_training1 = vector()
sigma2_sq_sample_mean_full_training1 = vector()

#constant parameters values
a_full_training1=10
lambda_full_training1=1000

#Min rho values
rho_sp = 0.01
rho_tmp = 3

SIGMA_sp_full_training1 = exp(-rho_sp*distance_mat_full_training1/1000) #to get data similar to temporal matrix values
SIGMA_tmp_full_training1 = exp(-rho_tmp*week_diff_mat)
inv_SIGMA_sp_full_training1 = chol2inv(chol(SIGMA_sp_full_training1))
inv_SIGMA_tmp_full_training1 = chol2inv(chol(SIGMA_tmp_full_training1))
inv_SIGMA_sp_kro_inv_SIGMA_tmp_full_training1 = inv_SIGMA_sp_full_training1 %x% inv_SIGMA_tmp_full_training1
#Giving values to different parameters of our model
v_mat_full_training1 = matrix(0,nrow = n_tmp_full_training1,ncol = n_sp_full_training1)
theta_full_training1 = rep(0,p+1)
sigma1_sq_full_training1 = 1
sigma2_sq_full_training1 = 1

#time taken analysis for parameters
total_tt_v_full_training1 = 0
total_tt_theta_full_training1 = 0
total_tt_sigma1_sq_full_training1 = 0
total_tt_sigma2_sq_full_training1 = 0

#Getting inv Sigma_sp_22 to use for every iteration
inv_SIGMA_sp_trans_22_matlist_full_training1=list()
SIGMA_sp_trans_21_mat_full_training1 = vector()
for (j in (1:n_sp_full_training1)) {
  idx = 1:n_sp_full_training1
  
  if(j != 1){
    idx[1]=j
    idx[j]=1
  }
  SIGMA_sp_trans_full_training1 = diag(n_sp_full_training1)[idx,] %*% SIGMA_sp_full_training1 %*% diag(n_sp_full_training1)[idx,]
  inv_SIGMA_sp_trans_full_training1 = diag(n_sp_full_training1)[idx,] %*% inv_SIGMA_sp_full_training1 %*% diag(n_sp_full_training1)[idx,]
  inv_SIGMA_sp_trans_a_full_training1 = inv_SIGMA_sp_trans_full_training1[1,1]
  inv_SIGMA_sp_trans_b_full_training1 = inv_SIGMA_sp_trans_full_training1[1,-1]
  inv_SIGMA_sp_trans_c_full_training1 = inv_SIGMA_sp_trans_full_training1[-1,1]
  inv_SIGMA_sp_trans_d_full_training1 = inv_SIGMA_sp_trans_full_training1[-1,-1]
  
  inv_SIGMA_sp_trans_22_full_training1 = inv_SIGMA_sp_trans_d_full_training1-(inv_SIGMA_sp_trans_c_full_training1%*%t(inv_SIGMA_sp_trans_b_full_training1))/inv_SIGMA_sp_trans_a_full_training1
  SIGMA_sp_trans_21_full_training1 = SIGMA_sp_trans_full_training1[-1,1]
  
  inv_SIGMA_sp_trans_22_matlist_full_training1[[j]] = inv_SIGMA_sp_trans_22_full_training1
  SIGMA_sp_trans_21_mat_full_training1 = cbind(SIGMA_sp_trans_21_mat_full_training1,SIGMA_sp_trans_21_full_training1)
  
}

theta_sample_full_training1 = vector()
sigma1_sq_sample_full_training1 = vector()
sigma2_sq_sample_full_training1 = vector()
v_sample_full_training1 = vector()
for (i in (1:(burn_period_full_training1+sample_size_full_training1*diff_in_random_draws_full_training1))) {
  #v posterior
  st_v = Sys.time()
  
  for (j in (1:n_sp_full_training1)) {
    idx=1:n_sp_full_training1
    
    if(j!=1){
      idx[1] = j
      idx[j] = 1
    }
    
    v_cond_full_training1 = (v_mat_full_training1 %*% diag(n_sp_full_training1)[idx,])[(1+n_tmp_full_training1):(n_full_training1)]
    
    #SIGMA' conditional without sigma1_sq as it's included in the posterior 
    SIGMA_cond_full_training1 = (1 -t(SIGMA_sp_trans_21_mat_full_training1[,j])%*%inv_SIGMA_sp_trans_22_matlist_full_training1[[j]] 
                                 %*%SIGMA_sp_trans_21_mat_full_training1[,j])%x% SIGMA_tmp_full_training1
    mu_cond_full_training1 = ((t(SIGMA_sp_trans_21_mat_full_training1[,j])%*%inv_SIGMA_sp_trans_22_matlist_full_training1[[j]])
                              %x% diag(n_tmp_full_training1))%*%v_cond_full_training1
    inv_SIGMA_cond_full_training1 = chol2inv(chol(SIGMA_cond_full_training1))
    
    v_covar_full_training1 = chol2inv(chol(inv_SIGMA_cond_full_training1 /sigma1_sq_full_training1 + diag(n_tmp_full_training1)/sigma2_sq_full_training1))
    v_mean_full_training1  = v_covar_full_training1 %*% (inv_SIGMA_cond_full_training1 %*% 
                                                           mu_cond_full_training1/sigma1_sq_full_training1 + 
                                                           (y_full_training1[(1+(j-1)*n_tmp_full_training1):(j*n_tmp_full_training1)] -
                                                              X_mat_full_training1[(1+(j-1)*n_tmp_full_training1):(j*n_tmp_full_training1),] %*% 
                                                              theta_full_training1)/sigma2_sq_full_training1 )
    v_mat_full_training1[,j] = v_mean_full_training1 + t(chol(v_covar_full_training1))  %*% rnorm(n_tmp_full_training1)
    j=j+1
  }
  
  total_tt_v_full_training1 = total_tt_v_full_training1 + as.numeric(difftime(Sys.time(), st_v), units="secs")
  
  #theta posterior
  st_theta = Sys.time()
  
  theta_covar_full_training1 = solve( (t(X_mat_full_training1) %*% X_mat_full_training1)/sigma2_sq_full_training1 + diag(p+1) ) 
  theta_mean_full_training1 = theta_covar_full_training1 %*% ( t(X_mat_full_training1) %*% (y_full_training1 - v_mat_full_training1[1:n_full_training1] ) )/sigma2_sq_full_training1
  theta_full_training1 = theta_mean_full_training1 + t(chol(theta_covar_full_training1)) %*%rnorm(p+1)
  
  total_tt_theta_full_training1 = total_tt_theta_full_training1 + as.numeric(difftime(Sys.time(), st_theta), units="secs")
  #Sigma1^2 posterior
  st_sigma1_sq = Sys.time()
  
  sigma1_sq_a_full_training1 = a_full_training1 + n_full_training1/2
  sigma1_sq_lambda_full_training1 = (t(v_mat_full_training1[1:n_full_training1]) %*% inv_SIGMA_sp_kro_inv_SIGMA_tmp_full_training1 
                                     %*% v_mat_full_training1[1:n_full_training1])/2 + lambda_full_training1
  sigma1_sq_full_training1 = invgamma::rinvgamma(1,sigma1_sq_a_full_training1,rate = sigma1_sq_lambda_full_training1)
  
  total_tt_sigma1_sq_full_training1 = total_tt_sigma1_sq_full_training1 + as.numeric(difftime(Sys.time(), st_sigma1_sq), units="secs")
  #Sigma2^2 posterior
  st_sigma2_sq = Sys.time()
  
  sigma2_sq_a_full_training1 = a_full_training1 + n_full_training1/2
  temp_vec = y_full_training1 -X_mat_full_training1%*%theta_full_training1 - v_mat_full_training1[1:n_full_training1]
  sigma2_sq_lambda_full_training1 = (t(temp_vec) %*% (temp_vec))/2 + lambda_full_training1
  sigma2_sq_full_training1 = invgamma::rinvgamma(1,sigma2_sq_a_full_training1,rate = sigma2_sq_lambda_full_training1)
  
  total_tt_sigma2_sq_full_training1 = total_tt_sigma2_sq_full_training1 + as.numeric(difftime(Sys.time(), st_sigma2_sq), units="secs")
  
  #Collecting samples for prediction
  if(i> burn_period_full_training1 & i%%diff_in_random_draws_full_training1 ==0){
    v_sample_full_training1 = cbind(v_sample_full_training1,v_mat_full_training1[1:n_full_training1])
    theta_sample_full_training1 = cbind(theta_sample_full_training1,theta_full_training1)
    sigma1_sq_sample_full_training1 = c(sigma1_sq_sample_full_training1,sigma1_sq_full_training1)
    sigma2_sq_sample_full_training1 = c(sigma2_sq_sample_full_training1,sigma2_sq_full_training1)
  }
}
#time analysis
avg_tt_v = total_tt_v_full_training1 / (burn_period_full_training1+sample_size_full_training1*diff_in_random_draws_full_training1)
avg_tt_theta = total_tt_theta_full_training1 / (burn_period_full_training1+sample_size_full_training1*diff_in_random_draws_full_training1)
avg_tt_sigma1_sq = total_tt_sigma1_sq_full_training1 / (burn_period_full_training1+sample_size_full_training1*diff_in_random_draws_full_training1)
avg_tt_sigma2_sq = total_tt_sigma2_sq_full_training1 / (burn_period_full_training1+sample_size_full_training1*diff_in_random_draws_full_training1)
#estimation of y with the posterior parameters
estimated_y_full_training1 = X_mat_full_training1 %*% rowMeans(theta_sample_full_training1) + rowMeans(v_sample_full_training1) + 
  rnorm(n_full_training1,mean = 0, sd= sqrt(mean(sigma2_sq_sample_full_training1)))
mape_full_training1 = (sum(abs( (y_full_training1 - estimated_y_full_training1)/y_full_training1 )))/n_full_training1
mape_vec_training_full_training1 = c(mape_vec_training_full_training1, mape_full_training1)

#Testing for significance
theta_significance_vec_full_training1=vector()
for (ts in (1:(p+1))) {
  theta_significance_vec_full_training1 = c(theta_significance_vec_full_training1, unname(!(0 >= quantile(theta_sample_full_training1[ts,], 0.025) & 
                                                                                              0 <= quantile(theta_sample_full_training1[ts,], 0.975))))
}

#Pred


v_pred_test_vec = vector()
y_pred_test_vec = vector()
y_pred_test_mat = vector()
pop_vec_testing = unique(testing[,c(2,6)])[,2]
for (pred_i in (1:n_testing)) {
  Sigma_12_vec = vector()
  for (j in (1:n_full_training1)) {
    sp_dist = distm(c(full_training1$long[j],  full_training1$lat[j]), c(testing$long[pred_i], testing$lat[pred_i]), fun = distHaversine)
    tmp_dist = abs((as.numeric(difftime(testing$week[pred_i] ,  full_training1$week[j]), units="days"))/7)
    Sigma_12_vec = c(Sigma_12_vec, exp(-rho_sp*  sp_dist/1000) * exp(-rho_tmp * tmp_dist))
  }
  v_tmp_given_v1_covar = 1 - t(Sigma_12_vec) %*% inv_SIGMA_sp_kro_inv_SIGMA_tmp_full_training1 %*% Sigma_12_vec
  v_tmp_given_v1_mean = t(Sigma_12_vec) %*% inv_SIGMA_sp_kro_inv_SIGMA_tmp_full_training1 %*% rowMeans(v_sample_full_training1)
  v_tmp_given_v1 = v_tmp_given_v1_mean + sqrt(mean(sigma1_sq_sample_full_training1)) * t(chol(v_tmp_given_v1_covar)) * rnorm(1)
  v_pred_test_vec = c(v_pred_test_vec, v_tmp_given_v1)
  
  
  if(pred_i %% n_sp_testing != 0){#which spatial location s=pred_i %% n_sp_validation1
    tmp = floor(pred_i / n_sp_testing) +1
    sp = pred_i %% n_sp_testing
    idx = (sp -1)*n_tmp_testing + tmp
    if(tmp == 1){#Implies time t=1
      if(unique(testing[,c(2,6)])[sp,1] %in% test_loc ){#location not trained 
        X_vec_tmp = c(1, log(testing$population[idx]), testing$time[idx], testing$time_sq[idx], 
                      100*log(0.1), testing$prev_log_death[idx], 
                      100*(testing$prev_log_prevalence_country[idx]))
      } else{
        X_vec_tmp = c(1, log(testing$population[idx]), testing$time[idx], testing$time_sq[idx], 
                      100*testing$prev_log_prevalence[idx], testing$prev_log_death[idx], 
                      100*(testing$prev_log_prevalence_country[idx]))
      }
    }else{#Implies t>1
      X_vec_tmp = c(1, log(testing$population[idx]), testing$time[idx], testing$time_sq[idx],
                    y_pred_test_mat[tmp - 1, sp], testing$prev_log_death[idx], 
                    100*log(sum((exp(y_pred_test_mat[tmp - 1, ]/100) - 0.1)*pop_vec_testing)/sum(pop_vec_testing) +0.1))
    }
    
    y_pred_test = t(X_vec_tmp) %*% rowMeans(theta_sample_full_training1) + v_tmp_given_v1 + 
      sqrt(mean(sigma2_sq_sample_full_training1))*rnorm(1)
    y_pred_test_vec = c(y_pred_test_vec, y_pred_test)
    
  }else{#location = n
    sp = n_sp_testing
    tmp = pred_i / n_sp_testing
    idx = (sp -1)*n_tmp_testing + tmp
    if(tmp == 1){#Implies time t=1
      if(unique(testing[,c(2,6)])[sp,1] %in% test_loc ){#location not trained 
        X_vec_tmp = c(1, log(testing$population[idx]), testing$time[idx], testing$time_sq[idx], 
                      100*log(0.1), testing$prev_log_death[idx] , 
                      100*(testing$prev_log_prevalence_country[idx]))
      } else{
        X_vec_tmp = c(1, log(testing$population[idx]), testing$time[idx], testing$time_sq[idx], 
                      100*testing$prev_log_prevalence[idx], testing$prev_log_death[idx], #log(testing$prev_new_death[idx] + 0.1)
                      100*(testing$prev_log_prevalence_country[idx]))
      }
    }else{#Implies t>1
      X_vec_tmp = c(1, log(testing$population[idx]), testing$time[idx], testing$time_sq[idx],
                    y_pred_test_mat[tmp - 1, sp], testing$prev_log_death[idx], 
                    100*log(sum((exp(y_pred_test_mat[tmp - 1, ]/100) - 0.1)*pop_vec_testing)/sum(pop_vec_testing) +0.1))
    }
    y_pred_test = t(X_vec_tmp) %*% rowMeans(theta_sample_full_training1) + v_tmp_given_v1 + 
      sqrt(mean(sigma2_sq_sample_full_training1))*rnorm(1)
    y_pred_test_vec = c(y_pred_test_vec, y_pred_test)
    y_pred_test_mat = rbind(y_pred_test_mat, y_pred_test_vec)
    y_pred_test_vec = vector()
    
  }
}

#prediction errors
mape_prediction_test = (sum(abs( (y_testing - y_pred_test_mat[1:n_testing])/y_testing )))/n_testing  
####################################################################################################################################
#estimating with training 2 data

#training 2 data values
y_training2 = 100*training2$log_prevalence
n_training2 = length(y_training2)
n_sp_training2 = length(unique(training2$state))
n_tmp_training2 = n_training2/n_sp_training2

#validation set values
y_validation2 = 100*validation2$log_prevalence
n_validation2 = length(y_validation2)
n_sp_validation2 = length(unique(validation2$state))
n_tmp_validation2 = n_validation2/n_sp_validation2
#total dataset values
n_total = length(df$state)
n_sp_total = length(unique(df$state))
n_tmp_total = n_total/ n_sp_total
#testing values
y_testing = 100*testing$log_prevalence
n_testing = length(y_testing)
n_sp_testing = length(unique(testing$state))
n_tmp_testing = n_testing/ n_sp_testing
p=4
#X matrix for the model
X_mat_training2 = cbind(rep(1,n_training2), log(training2$population), training2$time, training2$time_sq, 
                        100*(training2$prev_log_prevalence))

#Getting week difference mat for Sigma_T
week_vec = 1:n_tmp_training2
week_diff_mat = as.matrix(dist(week_vec,diag  =  TRUE,upper  =  TRUE))

#Getting distance matrix
unique_state_lat_long_training2 = unique(training2 %>% dplyr::select(state,lat,long))
distance_mat_training2 = as.matrix(distm(cbind(unique_state_lat_long_training2$long,unique_state_lat_long_training2$lat), 
                                         fun = distHaversine))
#Error vectors
mape_vec_training_training2 = vector()
mape_vec_validation_training2 = vector() 
#pred_error_vec_training2 = vector()

burn_period_training2 = 1000
sample_size_training2 = 500
diff_in_random_draws_training2 = 20
#Getting sample means for every rho
v_sample_mean_training2 = vector()
theta_sample_mean_training2 = vector()
theta_significance_training2 = vector()
sigma1_sq_sample_mean_training2 = vector()
sigma2_sq_sample_mean_training2 = vector()

#constant parameters values
a_training2=10
lambda_training2=1000

rho_vec_tmp_training2 = seq(0.1, 0.5, 0.1)#c(0.1, 0.5, 0.75, 1, 1.5, 3)
rho_vec_sp_training2 = seq(0.5, 1, 0.1)#c(0.01, 0.1, 0.5, 0.75, 1)

for(rho_tmp in rho_vec_tmp_training2){
  for (rho_sp in rho_vec_sp_training2) {
    
    SIGMA_sp_training2 = exp(-rho_sp*distance_mat_training2/80000) #to get data similar to temporal matrix values
    SIGMA_tmp_training2 = exp(-rho_tmp*week_diff_mat)
    inv_SIGMA_sp_training2 = chol2inv(chol(SIGMA_sp_training2))
    inv_SIGMA_tmp_training2 = chol2inv(chol(SIGMA_tmp_training2))
    inv_SIGMA_sp_kro_inv_SIGMA_tmp_training2 = inv_SIGMA_sp_training2 %x% inv_SIGMA_tmp_training2
    #Giving values to different parameters of our model
    v_mat_training2 = matrix(0,nrow = n_tmp_training2,ncol = n_sp_training2)
    theta_training2 = rep(0,p+1)
    sigma1_sq_training2 = 1
    sigma2_sq_training2 = 1
    
    #time taken analysis for parameters
    total_tt_v_training2 = 0
    total_tt_theta_training2 = 0
    total_tt_sigma1_sq_training2 = 0
    total_tt_sigma2_sq_training2 = 0
    
    #Getting inv Sigma_sp_22 to use for every iteration
    inv_SIGMA_sp_trans_22_matlist_training2=list()
    SIGMA_sp_trans_21_mat_training2 = vector()
    for (j in (1:n_sp_training2)) {
      idx = 1:n_sp_training2
      
      if(j != 1){
        idx[1]=j
        idx[j]=1
      }
      SIGMA_sp_trans_training2 = diag(n_sp_training2)[idx,] %*% SIGMA_sp_training2 %*% diag(n_sp_training2)[idx,]
      inv_SIGMA_sp_trans_training2 = diag(n_sp_training2)[idx,] %*% inv_SIGMA_sp_training2 %*% diag(n_sp_training2)[idx,]
      inv_SIGMA_sp_trans_a_training2 = inv_SIGMA_sp_trans_training2[1,1]
      inv_SIGMA_sp_trans_b_training2 = inv_SIGMA_sp_trans_training2[1,-1]
      inv_SIGMA_sp_trans_c_training2 = inv_SIGMA_sp_trans_training2[-1,1]
      inv_SIGMA_sp_trans_d_training2 = inv_SIGMA_sp_trans_training2[-1,-1]
      
      inv_SIGMA_sp_trans_22_training2 = inv_SIGMA_sp_trans_d_training2-(inv_SIGMA_sp_trans_c_training2%*%t(inv_SIGMA_sp_trans_b_training2))/inv_SIGMA_sp_trans_a_training2
      SIGMA_sp_trans_21_training2 = SIGMA_sp_trans_training2[-1,1]
      
      inv_SIGMA_sp_trans_22_matlist_training2[[j]] = inv_SIGMA_sp_trans_22_training2
      SIGMA_sp_trans_21_mat_training2 = cbind(SIGMA_sp_trans_21_mat_training2,SIGMA_sp_trans_21_training2)
      
    }
    
    theta_sample_training2 = vector()
    sigma1_sq_sample_training2 = vector()
    sigma2_sq_sample_training2 = vector()
    v_sample_training2 = vector()
    for (i in (1:(burn_period_training2+sample_size_training2*diff_in_random_draws_training2))) {
      #v posterior
      st_v = Sys.time()
      
      for (j in (1:n_sp_training2)) {
        idx=1:n_sp_training2
        
        if(j!=1){
          idx[1] = j
          idx[j] = 1
        }
        
        v_cond_training2 = (v_mat_training2 %*% diag(n_sp_training2)[idx,])[(1+n_tmp_training2):(n_training2)]
        
        #SIGMA' conditional without sigma1_sq as it's included in the posterior 
        SIGMA_cond_training2 = (1 -t(SIGMA_sp_trans_21_mat_training2[,j])%*%inv_SIGMA_sp_trans_22_matlist_training2[[j]] 
                                %*%SIGMA_sp_trans_21_mat_training2[,j]) %x% SIGMA_tmp_training2
        mu_cond_training2 = ((t(SIGMA_sp_trans_21_mat_training2[,j])%*%inv_SIGMA_sp_trans_22_matlist_training2[[j]])
                             %x% diag(n_tmp_training2))%*%v_cond_training2
        inv_SIGMA_cond_training2 = chol2inv(chol(SIGMA_cond_training2))
        
        v_covar_training2 = chol2inv(chol(inv_SIGMA_cond_training2 /sigma1_sq_training2 + diag(n_tmp_training2)/sigma2_sq_training2))
        v_mean_training2  = v_covar_training2 %*% (inv_SIGMA_cond_training2 %*% 
                                                     mu_cond_training2/sigma1_sq_training2 + 
                                                     (y_training2[(1+(j-1)*n_tmp_training2):(j*n_tmp_training2)] -
                                                        X_mat_training2[(1+(j-1)*n_tmp_training2):(j*n_tmp_training2),] %*% 
                                                        theta_training2)/sigma2_sq_training2 )
        v_mat_training2[,j] = v_mean_training2 + t(chol(v_covar_training2))  %*% rnorm(n_tmp_training2)
        j=j+1
      }
      
      total_tt_v_training2 = total_tt_v_training2 + as.numeric(difftime(Sys.time(), st_v), units="secs")
      
      #theta posterior
      st_theta = Sys.time()
      
      theta_covar_training2 = solve( (t(X_mat_training2) %*% X_mat_training2)/sigma2_sq_training2 + diag(p+1) ) 
      theta_mean_training2 = theta_covar_training2 %*% ( t(X_mat_training2) %*% (y_training2 - v_mat_training2[1:n_training2] ) )/sigma2_sq_training2
      theta_training2 = theta_mean_training2 + t(chol(theta_covar_training2)) %*%rnorm(p+1)
      
      total_tt_theta_training2 = total_tt_theta_training2 + as.numeric(difftime(Sys.time(), st_theta), units="secs")
      #Sigma1^2 posterior
      st_sigma1_sq = Sys.time()
      
      sigma1_sq_a_training2 = a_training2 + n_training2/2
      sigma1_sq_lambda_training2 = (t(v_mat_training2[1:n_training2]) %*% inv_SIGMA_sp_kro_inv_SIGMA_tmp_training2 
                                    %*% v_mat_training2[1:n_training2])/2 + lambda_training2
      sigma1_sq_training2 = invgamma::rinvgamma(1,sigma1_sq_a_training2,rate = sigma1_sq_lambda_training2)
      
      total_tt_sigma1_sq_training2 = total_tt_sigma1_sq_training2 + as.numeric(difftime(Sys.time(), st_sigma1_sq), units="secs")
      #Sigma2^2 posterior
      st_sigma2_sq = Sys.time()
      
      sigma2_sq_a_training2 = a_training2 + n_training2/2
      temp_vec = y_training2 -X_mat_training2%*%theta_training2 - v_mat_training2[1:n_training2]
      sigma2_sq_lambda_training2 = (t(temp_vec) %*% (temp_vec))/2 + lambda_training2
      sigma2_sq_training2 = invgamma::rinvgamma(1,sigma2_sq_a_training2,rate = sigma2_sq_lambda_training2)
      
      total_tt_sigma2_sq_training2 = total_tt_sigma2_sq_training2 + as.numeric(difftime(Sys.time(), st_sigma2_sq), units="secs")
      
      #Collecting samples for prediction
      if(i> burn_period_training2 & i%%diff_in_random_draws_training2 ==0){
        v_sample_training2 = cbind(v_sample_training2,v_mat_training2[1:n_training2])
        theta_sample_training2 = cbind(theta_sample_training2,theta_training2)
        sigma1_sq_sample_training2 = c(sigma1_sq_sample_training2,sigma1_sq_training2)
        sigma2_sq_sample_training2 = c(sigma2_sq_sample_training2,sigma2_sq_training2)
      }
    }
    #time analysis
    avg_tt_v = total_tt_v_training2 / (burn_period_training2+sample_size_training2*diff_in_random_draws_training2)
    avg_tt_theta = total_tt_theta_training2 / (burn_period_training2+sample_size_training2*diff_in_random_draws_training2)
    avg_tt_sigma1_sq = total_tt_sigma1_sq_training2 / (burn_period_training2+sample_size_training2*diff_in_random_draws_training2)
    avg_tt_sigma2_sq = total_tt_sigma2_sq_training2 / (burn_period_training2+sample_size_training2*diff_in_random_draws_training2)
    #estimation of y with the posterior parameters
    estimated_y_training2 = X_mat_training2 %*% rowMeans(theta_sample_training2) + rowMeans(v_sample_training2) + 
      rnorm(n_training2,mean = 0, sd= sqrt(mean(sigma2_sq_sample_training2)))
    mape_training2 = (sum(abs( (y_training2 - estimated_y_training2)/y_training2 )))/n_training2
    mape_vec_training_training2 = c(mape_vec_training_training2, mape_training2)
    
    #Testing for significance
    significance_vec_training2 = vector()
    significance_vec_training2 = c(significance_vec_training2, unname(!(0 >= quantile(theta_sample_training2[1,], 0.025) & 
                                                                          0 <= quantile(theta_sample_training2[1,], 0.975))))
    significance_vec_training2 = c(significance_vec_training2, unname(!(0 >= quantile(theta_sample_training2[2,], 0.025) & 
                                                                          0 <= quantile(theta_sample_training2[2,], 0.975))))
    significance_vec_training2 = c(significance_vec_training2, unname(!(0 >= quantile(theta_sample_training2[3,], 0.025) & 
                                                                          0 <= quantile(theta_sample_training2[3,], 0.975))))
    significance_vec_training2 = c(significance_vec_training2, unname(!(0 >= quantile(theta_sample_training2[4,], 0.025) & 
                                                                          0 <= quantile(theta_sample_training2[4,], 0.975))))
    significance_vec_training2 = c(significance_vec_training2, unname(!(0 >= quantile(theta_sample_training2[5,], 0.025) & 
                                                                          0 <= quantile(theta_sample_training2[5,], 0.975))))
    theta_significance_training2 = cbind(theta_significance_training2, significance_vec_training2)
    #Saving value for sample means
    v_sample_mean_training2 = cbind(v_sample_mean_training2, rowMeans(v_sample_training2))
    theta_sample_mean_training2 = cbind(theta_sample_mean_training2, rowMeans(theta_sample_training2))
    sigma1_sq_sample_mean_training2 = c(sigma1_sq_sample_mean_training2, mean(sigma1_sq_sample_training2))
    sigma2_sq_sample_mean_training2 = c(sigma2_sq_sample_mean_training2, mean(sigma2_sq_sample_training2))
    
    #Validation
    v_pred_val_vec_training2 = vector()
    y_pred_val_vec_training2 = vector()
    for (pred_i in (1:n_validation2)) {
      Sigma_12_vec = vector()
      for (j in (1:n_training2)) {
        sp_dist = distm(c(training2$long[j],  training2$lat[j]), c(validation2$long[pred_i], validation2$lat[pred_i]), fun = distHaversine)
        tmp_dist = abs((as.numeric(difftime(validation2$week[pred_i] ,  training2$week[j]), units="days"))/7)
        Sigma_12_vec = c(Sigma_12_vec, exp(-rho_sp*  sp_dist/80000) * exp(-rho_tmp * tmp_dist))
      }
      v_tmp_given_v1_covar = 1 - t(Sigma_12_vec) %*% inv_SIGMA_sp_kro_inv_SIGMA_tmp_training2 %*% Sigma_12_vec
      v_tmp_given_v1_mean = t(Sigma_12_vec) %*% inv_SIGMA_sp_kro_inv_SIGMA_tmp_training2 %*% rowMeans(v_sample_training2)
      v_tmp_given_v1 = v_tmp_given_v1_mean + sqrt(mean(sigma1_sq_sample_training2)) * t(chol(v_tmp_given_v1_covar)) * rnorm(1)
      if(pred_i %% n_tmp_validation2 == 1){
        X_vec_tmp = c(1, log(validation2$population[pred_i]), validation2$time[pred_i], validation2$time_sq[pred_i], 
                      100*validation2$prev_log_prevalence[pred_i])
      }else{
        X_vec_tmp = c(1, log(validation2$population[pred_i]), validation2$time[pred_i], validation2$time_sq[pred_i], 
                      y_pred_val_vec_training2[pred_i - 1])
      }
      y_pred_val = X_vec_tmp %*% rowMeans(theta_sample_training2) + v_tmp_given_v1 + sqrt(mean(sigma2_sq_sample_training2))*rnorm(1)
      v_pred_val_vec_training2 = c(v_pred_val_vec_training2, v_tmp_given_v1)
      y_pred_val_vec_training2 = c(y_pred_val_vec_training2, y_pred_val)
    }
    mape_validation = (sum(abs( (y_validation2 - y_pred_val_vec_training2)/y_validation2 )))/n_validation2
    mape_vec_validation_training2 = c(mape_vec_validation_training2, mape_validation)
  }
}
#Prediction



#training 2 data values

y_full_training2 = 100*full_training2$log_prevalence
n_full_training2 = length(y_full_training2)
n_sp_full_training2 = length(unique(full_training2$state))
n_tmp_full_training2 = n_full_training2/n_sp_full_training2
#X matrix for the model
X_mat_full_training2 = cbind(rep(1,n_full_training2), log(full_training2$population), full_training2$time, full_training2$time_sq, 
                             100*(full_training2$prev_log_prevalence))


#Getting week difference mat for Sigma_T
week_vec = 1:n_tmp_full_training2
week_diff_mat = as.matrix(dist(week_vec,diag  =  TRUE,upper  =  TRUE))

#Getting distance matrix
unique_state_lat_long_full_training2 = unique(full_training2 %>% dplyr::select(state,lat,long))
distance_mat_full_training2 = as.matrix(distm(cbind(unique_state_lat_long_full_training2$long,unique_state_lat_long_full_training2$lat), 
                                              fun = distHaversine))
#Error vectors
mape_vec_training_full_training2 = vector()
mape_vec_validation_full_training2 = vector() 
#pred_error_vec_full_training2 = vector()

burn_period_full_training2 = 1000
sample_size_full_training2 = 500
diff_in_random_draws_full_training2 = 20
#Getting sample means for every rho
v_sample_mean_full_training2 = vector()
theta_sample_mean_full_training2 = vector()
theta_significance_full_training2 = vector()
sigma1_sq_sample_mean_full_training2 = vector()
sigma2_sq_sample_mean_full_training2 = vector()

#constant parameters values
a_full_training2=10
lambda_full_training2=1000

#Min rho values
rho_sp = 1
rho_tmp = 0.1

SIGMA_sp_full_training2 = exp(-rho_sp*distance_mat_full_training2/80000) #to get data similar to temporal matrix values
SIGMA_tmp_full_training2 = exp(-rho_tmp*week_diff_mat)
inv_SIGMA_sp_full_training2 = chol2inv(chol(SIGMA_sp_full_training2))
inv_SIGMA_tmp_full_training2 = chol2inv(chol(SIGMA_tmp_full_training2))
inv_SIGMA_sp_kro_inv_SIGMA_tmp_full_training2 = inv_SIGMA_sp_full_training2 %x% inv_SIGMA_tmp_full_training2
#Giving values to different parameters of our model
v_mat_full_training2 = matrix(0,nrow = n_tmp_full_training2,ncol = n_sp_full_training2)
theta_full_training2 = rep(0,p+1)
sigma1_sq_full_training2 = 1
sigma2_sq_full_training2 = 1

#time taken analysis for parameters
total_tt_v_full_training2 = 0
total_tt_theta_full_training2 = 0
total_tt_sigma1_sq_full_training2 = 0
total_tt_sigma2_sq_full_training2 = 0

#Getting inv Sigma_sp_22 to use for every iteration
inv_SIGMA_sp_trans_22_matlist_full_training2=list()
SIGMA_sp_trans_21_mat_full_training2 = vector()
for (j in (1:n_sp_full_training2)) {
  idx = 1:n_sp_full_training2
  
  if(j != 1){
    idx[1]=j
    idx[j]=1
  }
  SIGMA_sp_trans_full_training2 = diag(n_sp_full_training2)[idx,] %*% SIGMA_sp_full_training2 %*% diag(n_sp_full_training2)[idx,]
  inv_SIGMA_sp_trans_full_training2 = diag(n_sp_full_training2)[idx,] %*% inv_SIGMA_sp_full_training2 %*% diag(n_sp_full_training2)[idx,]
  inv_SIGMA_sp_trans_a_full_training2 = inv_SIGMA_sp_trans_full_training2[1,1]
  inv_SIGMA_sp_trans_b_full_training2 = inv_SIGMA_sp_trans_full_training2[1,-1]
  inv_SIGMA_sp_trans_c_full_training2 = inv_SIGMA_sp_trans_full_training2[-1,1]
  inv_SIGMA_sp_trans_d_full_training2 = inv_SIGMA_sp_trans_full_training2[-1,-1]
  
  inv_SIGMA_sp_trans_22_full_training2 = inv_SIGMA_sp_trans_d_full_training2-(inv_SIGMA_sp_trans_c_full_training2%*%t(inv_SIGMA_sp_trans_b_full_training2))/inv_SIGMA_sp_trans_a_full_training2
  SIGMA_sp_trans_21_full_training2 = SIGMA_sp_trans_full_training2[-1,1]
  
  inv_SIGMA_sp_trans_22_matlist_full_training2[[j]] = inv_SIGMA_sp_trans_22_full_training2
  SIGMA_sp_trans_21_mat_full_training2 = cbind(SIGMA_sp_trans_21_mat_full_training2,SIGMA_sp_trans_21_full_training2)
  
}

theta_sample_full_training2 = vector()
sigma1_sq_sample_full_training2 = vector()
sigma2_sq_sample_full_training2 = vector()
v_sample_full_training2 = vector()
for (i in (1:(burn_period_full_training2+sample_size_full_training2*diff_in_random_draws_full_training2))) {
  #v posterior
  st_v = Sys.time()
  
  for (j in (1:n_sp_full_training2)) {
    idx=1:n_sp_full_training2
    
    if(j!=1){
      idx[1] = j
      idx[j] = 1
    }
    
    v_cond_full_training2 = (v_mat_full_training2 %*% diag(n_sp_full_training2)[idx,])[(1+n_tmp_full_training2):(n_full_training2)]
    
    #SIGMA' conditional without sigma1_sq as it's included in the posterior 
    SIGMA_cond_full_training2 = (1 -t(SIGMA_sp_trans_21_mat_full_training2[,j])%*%inv_SIGMA_sp_trans_22_matlist_full_training2[[j]] 
                                 %*%SIGMA_sp_trans_21_mat_full_training2[,j])%x% SIGMA_tmp_full_training2
    mu_cond_full_training2 = ((t(SIGMA_sp_trans_21_mat_full_training2[,j])%*%inv_SIGMA_sp_trans_22_matlist_full_training2[[j]])
                              %x% diag(n_tmp_full_training2))%*%v_cond_full_training2
    inv_SIGMA_cond_full_training2 = chol2inv(chol(SIGMA_cond_full_training2))
    
    v_covar_full_training2 = chol2inv(chol(inv_SIGMA_cond_full_training2 /sigma1_sq_full_training2 + diag(n_tmp_full_training2)/sigma2_sq_full_training2))
    v_mean_full_training2  = v_covar_full_training2 %*% (inv_SIGMA_cond_full_training2 %*% 
                                                           mu_cond_full_training2/sigma1_sq_full_training2 + 
                                                           (y_full_training2[(1+(j-1)*n_tmp_full_training2):(j*n_tmp_full_training2)] -
                                                              X_mat_full_training2[(1+(j-1)*n_tmp_full_training2):(j*n_tmp_full_training2),] %*% 
                                                              theta_full_training2)/sigma2_sq_full_training2 )
    v_mat_full_training2[,j] = v_mean_full_training2 + t(chol(v_covar_full_training2))  %*% rnorm(n_tmp_full_training2)
    j=j+1
  }
  
  total_tt_v_full_training2 = total_tt_v_full_training2 + as.numeric(difftime(Sys.time(), st_v), units="secs")
  
  #theta posterior
  st_theta = Sys.time()
  
  theta_covar_full_training2 = solve( (t(X_mat_full_training2) %*% X_mat_full_training2)/sigma2_sq_full_training2 + diag(p+1) ) 
  theta_mean_full_training2 = theta_covar_full_training2 %*% ( t(X_mat_full_training2) %*% (y_full_training2 - v_mat_full_training2[1:n_full_training2] ) )/sigma2_sq_full_training2
  theta_full_training2 = theta_mean_full_training2 + t(chol(theta_covar_full_training2)) %*%rnorm(p+1)
  
  total_tt_theta_full_training2 = total_tt_theta_full_training2 + as.numeric(difftime(Sys.time(), st_theta), units="secs")
  #Sigma1^2 posterior
  st_sigma1_sq = Sys.time()
  
  sigma1_sq_a_full_training2 = a_full_training2 + n_full_training2/2
  sigma1_sq_lambda_full_training2 = (t(v_mat_full_training2[1:n_full_training2]) %*% inv_SIGMA_sp_kro_inv_SIGMA_tmp_full_training2 
                                     %*% v_mat_full_training2[1:n_full_training2])/2 + lambda_full_training2
  sigma1_sq_full_training2 = invgamma::rinvgamma(1,sigma1_sq_a_full_training2,rate = sigma1_sq_lambda_full_training2)
  
  total_tt_sigma1_sq_full_training2 = total_tt_sigma1_sq_full_training2 + as.numeric(difftime(Sys.time(), st_sigma1_sq), units="secs")
  #Sigma2^2 posterior
  st_sigma2_sq = Sys.time()
  
  sigma2_sq_a_full_training2 = a_full_training2 + n_full_training2/2
  temp_vec = y_full_training2 -X_mat_full_training2%*%theta_full_training2 - v_mat_full_training2[1:n_full_training2]
  sigma2_sq_lambda_full_training2 = (t(temp_vec) %*% (temp_vec))/2 + lambda_full_training2
  sigma2_sq_full_training2 = invgamma::rinvgamma(1,sigma2_sq_a_full_training2,rate = sigma2_sq_lambda_full_training2)
  
  total_tt_sigma2_sq_full_training2 = total_tt_sigma2_sq_full_training2 + as.numeric(difftime(Sys.time(), st_sigma2_sq), units="secs")
  
  #Collecting samples for prediction
  if(i> burn_period_full_training2 & i%%diff_in_random_draws_full_training2 ==0){
    v_sample_full_training2 = cbind(v_sample_full_training2,v_mat_full_training2[1:n_full_training2])
    theta_sample_full_training2 = cbind(theta_sample_full_training2,theta_full_training2)
    sigma1_sq_sample_full_training2 = c(sigma1_sq_sample_full_training2,sigma1_sq_full_training2)
    sigma2_sq_sample_full_training2 = c(sigma2_sq_sample_full_training2,sigma2_sq_full_training2)
  }
}
#time analysis
avg_tt_v = total_tt_v_full_training2 / (burn_period_full_training2+sample_size_full_training2*diff_in_random_draws_full_training2)
avg_tt_theta = total_tt_theta_full_training2 / (burn_period_full_training2+sample_size_full_training2*diff_in_random_draws_full_training2)
avg_tt_sigma1_sq = total_tt_sigma1_sq_full_training2 / (burn_period_full_training2+sample_size_full_training2*diff_in_random_draws_full_training2)
avg_tt_sigma2_sq = total_tt_sigma2_sq_full_training2 / (burn_period_full_training2+sample_size_full_training2*diff_in_random_draws_full_training2)
#estimation of y with the posterior parameters
estimated_y_full_training2 = X_mat_full_training2 %*% rowMeans(theta_sample_full_training2) + rowMeans(v_sample_full_training2) + 
  rnorm(n_full_training2,mean = 0, sd= sqrt(mean(sigma2_sq_sample_full_training2)))
mape_full_training2 = (sum(abs( (y_full_training2 - estimated_y_full_training2)/y_full_training2 )))/n_full_training2
mape_vec_training_full_training2 = c(mape_vec_training_full_training2, mape_full_training2)

#Testing for significance
significance_vec_full_training2=vector()
significance_vec_full_training2 = c(significance_vec_full_training2, unname(!(0 >= quantile(theta_sample_full_training2[1,], 0.025) & 
                                                  0 <= quantile(theta_sample_full_training2[1,], 0.975))))
significance_vec_full_training2 = c(significance_vec_full_training2, unname(!(0 >= quantile(theta_sample_full_training2[2,], 0.025) & 
                                                  0 <= quantile(theta_sample_full_training2[2,], 0.975))))
significance_vec_full_training2 = c(significance_vec_full_training2, unname(!(0 >= quantile(theta_sample_full_training2[3,], 0.025) & 
                                                  0 <= quantile(theta_sample_full_training2[3,], 0.975))))
significance_vec_full_training2 = c(significance_vec_full_training2, unname(!(0 >= quantile(theta_sample_full_training2[4,], 0.025) & 
                                                  0 <= quantile(theta_sample_full_training2[4,], 0.975))))
significance_vec_full_training2 = c(significance_vec_full_training2, unname(!(0 >= quantile(theta_sample_full_training2[5,], 0.025) & 
                                                  0 <= quantile(theta_sample_full_training2[5,], 0.975))))
#Pred


v_pred_test_vec_training2 = vector()
y_pred_test_vec_training2 = vector()
for (pred_i in (1:n_testing)) {
  Sigma_12_vec = vector()
  for (j in (1:n_full_training2)) {
    sp_dist = distm(c(full_training2$long[j],  full_training2$lat[j]), c(testing$long[pred_i], testing$lat[pred_i]), fun = distHaversine)
    tmp_dist = abs((as.numeric(difftime(testing$week[pred_i] ,  full_training2$week[j]), units="days"))/7)
    Sigma_12_vec = c(Sigma_12_vec, exp(-rho_sp*  sp_dist/80000) * exp(-rho_tmp * tmp_dist))
  }
  v_tmp_given_v1_covar = 1 - t(Sigma_12_vec) %*% inv_SIGMA_sp_kro_inv_SIGMA_tmp_full_training2 %*% Sigma_12_vec
  v_tmp_given_v1_mean = t(Sigma_12_vec) %*% inv_SIGMA_sp_kro_inv_SIGMA_tmp_full_training2 %*% rowMeans(v_sample_full_training2)
  v_tmp_given_v1 = v_tmp_given_v1_mean + sqrt(mean(sigma1_sq_sample_full_training2)) * t(chol(v_tmp_given_v1_covar)) * rnorm(1)
  v_pred_test_vec_training2 = c(v_pred_test_vec_training2, v_tmp_given_v1)
  if(pred_i %% n_tmp_testing == 1){
    X_vec_tmp = c(1, log(testing$population[pred_i]), testing$time[pred_i], testing$time_sq[pred_i], 
                  100*testing$prev_log_prevalence[pred_i])
  }else{
    X_vec_tmp = c(1, log(testing$population[pred_i]), testing$time[pred_i], testing$time_sq[pred_i], 
                  y_pred_test_vec_training2[pred_i - 1])
  }
  y_pred_test = X_vec_tmp %*% rowMeans(theta_sample_full_training2) + v_tmp_given_v1 + sqrt(mean(sigma2_sq_sample_full_training2))*rnorm(1)
  v_pred_val_vec = c(v_pred_val_vec, v_tmp_given_v1)
  y_pred_test_vec_training2 = c(y_pred_test_vec_training2, y_pred_test)
}

#prediction errors
mape_prediction_test_training2 = (sum(abs( (y_testing - y_pred_test_vec_training2)/y_testing )))/n_testing


###################################################################################################################################


#Adding columns for fitted and predicted values for training set 1

log_prev_model_vec = vector()#c(estimated_y_full_training2, y_pred_test_vec_training2)
is_trained_vec = vector()
for (loc in (1:n_sp)) {
  if(unique(df$state)[loc] %in% test_loc){
    log_prev_model_vec = c(log_prev_model_vec, rep(0,n_tmp_full_training2),
                           y_pred_test_vec_training2[(((loc - 1)*n_tmp_testing+1):(loc*n_tmp_testing))])
    is_trained_vec = c(is_trained_vec, rep(0,n_tmp_full_training2),rep(0,n_tmp_testing))
  }else{
    log_prev_model_vec = c(log_prev_model_vec, 
                           estimated_y_full_training2[(((loc - 1)*n_tmp_full_training2+1):(loc*n_tmp_full_training2))],
                           y_pred_test_vec_training2[(((loc - 1)*n_tmp_testing+1):(loc*n_tmp_testing))])
    is_trained_vec = c(is_trained_vec, rep(1,n_tmp_full_training2),rep(0,n_tmp_testing))
  }
  
}
log_prev_model_vec = log_prev_model_vec/100
df$is_trained = is_trained_vec
df$log_prev_model_values = log_prev_model_vec

write.csv(df,"D:/Term5/IS&RA/US_statewise_weekly_withfittedprevnewdeath_train1model.csv",row.names = FALSE)
##Adding columns for fitted and predicted values for training set 2

log_prev_model_vec = vector()#c(estimated_y_full_training2, y_pred_test_vec_training2)
is_trained_vec = vector()
for (loc in (1:n_sp)) {
  log_prev_model_vec = c(log_prev_model_vec, 
                         estimated_y_full_training2[(((loc - 1)*n_tmp_full_training2+1):(loc*n_tmp_full_training2))],
                         y_pred_test_vec_training2[(((loc - 1)*n_tmp_testing+1):(loc*n_tmp_testing))])
  is_trained_vec = c(is_trained_vec, rep(1,n_tmp_full_training2),rep(0,n_tmp_testing))
}
log_prev_model_vec = log_prev_model_vec/100
df$is_trained = is_trained_vec
df$log_prev_model_values = log_prev_model_vec

write.csv(df,"D:/Term5/IS&RA/US_statewise_weekly_withfitted_model.csv",row.names = FALSE)
