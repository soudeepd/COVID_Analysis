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

set.seed(2020)

# read the statewise weekly data
df = read.csv(file = "D:/Term5/IS&RA/gitrepo/COVID_Analysis/Code/US_statewise_weekly.csv")

# check class of the columns in the data and change their class if needed
sapply(df,class) 
df = df %>% 
  mutate_if(is.factor,as.character) %>% 
  mutate(
    week = as.Date(week)
  )

# case1: use all but two random locations as training and last four weeks as test set 
# case2: use all locations for training and last four weeks as test set 
test_loc = sample(unique(df$state),2)
test_wk_min = max(df$week) - 28
training1 = df %>% filter(! state %in% test_loc & week <= test_wk_min)
training2 = df %>% filter(week <= test_wk_min)
testing = df %>% filter(week > test_wk_min)

#estimating with the first training data

unique_state_lat_long = unique(training1 %>% dplyr::select(state,lat,long))
distance_mat_training1 = as.matrix(distm(cbind(unique_state_lat_long$long,unique_state_lat_long$lat), fun = distHaversine))
#training 1 values
y_training1 = 100*training1$log_incidence
n_training1 = length(y_training1)
n_sp_training1 = length(unique_state_lat_long$state)
n_tmp_training1 = n_training1/n_sp_training1
#total dataset values
n_total = length(df$state)
n_sp_total = length(unique(df$state))
n_tmp_total = n_total/ n_sp_total
#testing values
y_testing = 100*testing$log_incidence
n_testing = length(y_testing)
n_sp_testing = length(unique(testing$state))
n_tmp_testing = n_testing/ n_sp_testing
p=2
#X matrix for the model
X_row_vec = vector()
X_mat_training1 = vector()
for (j in (1:n_training1)) {
  if(j %% n_tmp_training1==0){
    X_row_vec = c(log(training1$population[j]),n_tmp_training1/n_tmp_total)
  }else{
    X_row_vec = c(log(training1$population[j]), (j %% n_tmp_training1)/n_tmp_total)
  }
  X_mat_training1 = rbind(X_mat_training1,X_row_vec)
  X_row_vec = vector()
}

#Adding column of 1s for theta_0
X_mat_training1 = cbind(rep(1,n_training1),X_mat_training1)

#Getting week difference mat for Sigma_T
week_vec = 1:n_tmp_training1
week_diff_mat = as.matrix(dist(week_vec,diag  =  TRUE,upper  =  TRUE))

mse_vec_training1=vector()
pred_error_vec_training1=vector()

burn_period_training1 = 1000
sample_size_training1 = 500
diff_in_random_draws_training1 = 20
#Getting sample means for every rho
v_sample_mean_training1 = vector()
theta_sample_mean_training1 = vector()
sigma1_sq_sample_mean_training1 = vector()
sigma2_sq_sample_mean_training1 = vector()

#Getting sample means for every rho
v_sample_mean_training1 = vector()
theta_sample_mean_training1 = vector()
sigma1_sq_sample_mean_training1 = vector()
sigma2_sq_sample_mean_training1 = vector()

#constant parameters values
a_training1=10
lambda_training1=1000

rho_vec_tmp_training1 = seq(0.1, 1, 0.3)
rho_vec_sp_training1 = seq(0.1, 1, 0.3)
                         
for(rho_tmp in rho_vec_tmp_training1){
  for (rho_sp in rho_vec_sp_training1) {
   
    SIGMA_sp_training1 = exp(-rho_sp*distance_mat_training1/(80000)) # dividing by 80000 to get matrix values similar to temporal values
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
        SIGMA_cond_training1 = (1 -t(SIGMA_sp_trans_21_mat_training1[,j])%*%inv_SIGMA_sp_trans_22_matlist_training1[[j]] %*%SIGMA_sp_trans_21_mat_training1[,j])%x% SIGMA_tmp_training1
        mu_cond_training1 = ((t(SIGMA_sp_trans_21_mat_training1[,j])%*%inv_SIGMA_sp_trans_22_matlist_training1[[j]])%x% diag(n_tmp_training1))%*%v_cond_training1
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
    mse_training1 = (sum((y_training1 - estimated_y_training1)^2))/n_training1
    mse_vec_training1 = c(mse_vec_training1, mse_training1)
    
    #Saving value for sample means
    v_sample_mean_training1 = cbind(v_sample_mean_training1, rowMeans(v_sample_training1))
    theta_sample_mean_training1 = cbind(theta_sample_mean_training1, rowMeans(theta_sample_training1))
    sigma1_sq_sample_mean_training1 = c(sigma1_sq_sample_mean_training1, mean(sigma1_sq_sample_training1))
    sigma2_sq_sample_mean_training1 = c(sigma2_sq_sample_mean_training1, mean(sigma2_sq_sample_training1))
    
    #Prediction
    X_mat_pred_training1 = vector()
    for (j in (1:n_testing)) {
      if(j %% n_tmp_testing == 0){
        X_row_vec = c(log(testing$population[j]),(n_tmp_training1+n_tmp_testing)/n_tmp_total)
      }else{
        X_row_vec = c(log(testing$population[j]), (j %% n_tmp_training1 + n_tmp_training1)/n_tmp_total)
      }
      X_mat_pred_training1 = rbind(X_mat_pred_training1,X_row_vec)
      X_row_vec = vector()
    }
    X_mat_pred_training1 = cbind(rep(1, n_testing), X_mat_pred_training1)
    
    
    st_v_pred = Sys.time()
    v_pred_vec = vector()
    for (pred_i in (1:n_testing)) {
      Sigma_12_vec = vector()
      for (j in (1:n_training1)) {
        sp_dist = distm(c(training1$long[j],  training1$lat[j]), c(testing$long[pred_i], testing$lat[pred_i]), fun = distHaversine)
        tmp_dist = abs((as.numeric(difftime(testing$week[pred_i] ,  training1$week[j]), units="days"))/7)
        Sigma_12_vec = c(Sigma_12_vec, exp(-rho_sp*  sp_dist/80000) * exp(-rho_tmp * tmp_dist))#dividng by 8000 to get data similar to temp
      }
      v_tmp_given_v1_covar = 1 - t(Sigma_12_vec) %*% inv_SIGMA_sp_kro_inv_SIGMA_tmp_training1 %*% Sigma_12_vec
      v_tmp_given_v1_mean = t(Sigma_12_vec) %*% inv_SIGMA_sp_kro_inv_SIGMA_tmp_training1 %*% rowMeans(v_sample_training1)
      v_tmp_given_v1 = v_tmp_given_v1_mean + sqrt(mean(sigma1_sq_sample_training1)) * t(chol(v_tmp_given_v1_covar)) * rnorm(1)
      v_pred_vec = c(v_pred_vec, v_tmp_given_v1)
    }
    Sys.time() -st_v_pred
    #prediction errors
    y2_given_y1 = X_mat_pred_training1%*% rowMeans(theta_sample_training1) + v_pred_vec + sqrt(mean(sigma2_sq_sample_training1))*rnorm(n_testing)
    y2_pred_error = (sum((y_testing - y2_given_y1)^2))/n_testing
    pred_error_vec_training1 = c(pred_error_vec_training1, y2_pred_error)
  }
}

####################################################################################################################################
#estimating with training 2 data

#training 2 data values
y_training2 = 100*training2$log_incidence
n_training2 = length(y_training2)
n_sp_training2 = length(unique(training2$state))
n_tmp_training2 = n_training2/n_sp_training2

#X matrix for the model
X_row_vec = vector()
X_mat_training2 = vector()
for (j in (1:n_training2)) {
  if(j %% n_tmp_training2==0){
    X_row_vec = c(log(training2$population[j]),n_tmp_training2/n_tmp_total)
  }else{
    X_row_vec = c(log(training2$population[j]), (j %% n_tmp_training2)/n_tmp_total)
  }
  X_mat_training2 = rbind(X_mat_training2,X_row_vec)
  X_row_vec = vector()
}

#Adding column of 1s for theta_0
X_mat_training2 = cbind(rep(1,n_training2),X_mat_training2)

#Getting week difference mat for Sigma_T
week_vec = 1:n_tmp_training2
week_diff_mat = as.matrix(dist(week_vec,diag  =  TRUE,upper  =  TRUE))

#Getting distance matrix
unique_state_lat_long_training2 = unique(training2 %>% dplyr::select(state,lat,long))
distance_mat_training2 = as.matrix(distm(cbind(unique_state_lat_long_training2$long,unique_state_lat_long_training2$lat), 
                                         fun = distHaversine))

mse_vec_training2=vector()
pred_error_vec_training2=vector()

burn_period_training2 = 1000
sample_size_training2 = 500
diff_in_random_draws_training2 = 20
#Getting sample means for every rho
v_sample_mean_training2 = vector()
theta_sample_mean_training2 = vector()
sigma1_sq_sample_mean_training2 = vector()
sigma2_sq_sample_mean_training2 = vector()

#constant parameters values
a_training2=10
lambda_training2=1000

rho_vec_tmp_training2 = seq(0.1, 1, 0.3)
rho_vec_sp_training2 = seq(0.1, 1, 0.3)

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
        SIGMA_cond_training2 = (1 -t(SIGMA_sp_trans_21_mat_training2[,j])%*%inv_SIGMA_sp_trans_22_matlist_training2[[j]] %*%SIGMA_sp_trans_21_mat_training2[,j])%x% SIGMA_tmp_training2
        mu_cond_training2 = ((t(SIGMA_sp_trans_21_mat_training2[,j])%*%inv_SIGMA_sp_trans_22_matlist_training2[[j]])%x% diag(n_tmp_training2))%*%v_cond_training2
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
    mse_training2 = (sum((y_training2 - estimated_y_training2)^2))/n_training2
    mse_vec_training2 = c(mse_vec_training2, mse_training2)
    
    #Saving value for sample means
    v_sample_mean_training2 = cbind(v_sample_mean_training2, rowMeans(v_sample_training2))
    theta_sample_mean_training2 = cbind(theta_sample_mean_training2, rowMeans(theta_sample_training2))
    sigma1_sq_sample_mean_training2 = c(sigma1_sq_sample_mean_training2, mean(sigma1_sq_sample_training2))
    sigma2_sq_sample_mean_training2 = c(sigma2_sq_sample_mean_training2, mean(sigma2_sq_sample_training2))
    
    
    #Prediction
    
    X_mat_pred_training2 = vector()
    for (j in (1:n_testing)) {
      if(j %% n_tmp_testing == 0){
        X_row_vec = c(log(testing$population[j]),(n_tmp_training2+n_tmp_testing)/n_tmp_total)
      }else{
        X_row_vec = c(log(testing$population[j]), (j %% n_tmp_training2 + n_tmp_training2)/n_tmp_total)
      }
      X_mat_pred_training2 = rbind(X_mat_pred_training2,X_row_vec)
      X_row_vec = vector()
    }
    X_mat_pred_training2 = cbind(rep(1, n_testing), X_mat_pred_training2)
    
    v_pred_vec = vector()
    for (pred_i in (1:n_testing)) {
      Sigma_12_vec = vector()
      for (j in (1:n_training2)) {
        sp_dist = distm(c(training2$long[j],  training2$lat[j]), c(testing$long[pred_i], testing$lat[pred_i]), fun = distHaversine)
        tmp_dist = abs((as.numeric(difftime(testing$week[pred_i] ,  training2$week[j]), units="days"))/7)
        Sigma_12_vec = c(Sigma_12_vec, exp(-rho_sp*  sp_dist/80000) * exp(-rho_tmp * tmp_dist))
      }
      v_tmp_given_v1_covar = 1 - t(Sigma_12_vec) %*% inv_SIGMA_sp_kro_inv_SIGMA_tmp_training2 %*% Sigma_12_vec
      v_tmp_given_v1_mean = t(Sigma_12_vec) %*% inv_SIGMA_sp_kro_inv_SIGMA_tmp_training2 %*% rowMeans(v_sample_training2)
      v_tmp_given_v1 = v_tmp_given_v1_mean + sqrt(mean(sigma1_sq_sample_training2)) * t(chol(v_tmp_given_v1_covar)) * rnorm(1)
      v_pred_vec = c(v_pred_vec, v_tmp_given_v1)
    }
    
    #prediction errors
    y2_given_y1 = X_mat_pred_training2%*% rowMeans(theta_sample_training2) + v_pred_vec + sqrt(mean(sigma2_sq_sample_training2))*rnorm(n_testing)
    y2_pred_error = (sum((y_testing - y2_given_y1)^2))/n_testing
    pred_error_vec_training2 = c(pred_error_vec_training2, y2_pred_error)
  }
}