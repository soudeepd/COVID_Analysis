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

unique_state_lat_long = unique(training1 %>% select("state","lat","long"))
distance_mat_training1 = as.matrix(distm(cbind(unique_state_lat_long$long,unique_state_lat_long$lat), fun = distHaversine))
y_training1 = training1$log_incidence
n_training1 = length(y_training1)
n_sp_training1 = length(unique_state_lat_long$state)
n_tmp_training1 = n_training1/n_sp_training1
p=12
#X matrix for the model
X_row_vec = t(rep(0,12))
X_mat_training1 = vector()
for (j in (1:n_training1)) {
  X_row_vec[strtoi(format(as.Date(training1$week[j]), "%m"))] = 1 #changed to all data
  X_mat_training1 = rbind(X_mat_training1,X_row_vec)
  X_row_vec = t(rep(0,12))
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
#constant parameters values
a_training1=10
lambda_training1=1000

rho_vec_tmp_training1 = c(1)#seq(0.1,1,0.1)
rho_vec_sp_training1 = c(1)

for(rho_tmp in rho_vec_tmp_training1){
  for (rho_sp in rho_vec_sp_training1) {
   
    SIGMA_sp = exp(-rho_sp*distance_mat_training1)
    SIGMA_tmp = exp(-rho_tmp*week_diff_mat)
    inv_SIGMA_sp = chol2inv(chol(SIGMA_sp))
    inv_SIGMA_tmp = chol2inv(chol(SIGMA_tmp))
    inv_SIGMA_sp_kro_inv_SIGMA_tmp = inv_SIGMA_sp %x% inv_SIGMA_tmp
    #Giving values to different parameters of our model
    v_training1 = rep(0,n_training1)
    theta_training1 = rep(0,p+1)
    sigma1_sq_training1 = 1
    sigma2_sq_training1 = 1
    
    theta_sample_training1=vector()
    sigma1_sq_sample_training1=vector()
    sigma2_sq_sample_training1=vector()
    v_sample_training1=vector()
    for (i in (1:(burn_period_training1+sample_size_training1*diff_in_random_draws_training1))) {
      #v posterior
      v_covar_training1 = chol2inv((chol((inv_SIGMA_sp_kro_inv_SIGMA_tmp)/sigma1_sq_training1 + diag(n_training1)/sigma2_sq_training1)))
      v_mean_training1  = v_covar_training1 %*% ((y_training1-X_mat_training1 %*% theta_training1)/sigma2_sq_training1 )
      v_training1 = v_mean_training1 + t(chol(v_covar_training1))  %*% rnorm(n_training1)
      
      #theta posterior
      theta_covar_training1 = solve( (t(X_mat_training1) %*% X_mat_training1)/sigma2_sq_training1 + diag(p+1) ) 
      theta_mean_training1 = theta_covar_training1 %*% ( t(X_mat_training1) %*% (y_training1 - v_training1 ) )/sigma2_sq_training1
      theta_training1 = theta_mean_training1 + t(chol(theta_covar_training1)) %*%rnorm(p+1)
      #Sigma1^2 posterior
      sigma1_sq_a_training1 = a_training1 + n_training1/2
      sigma1_sq_lambda_training1 = (t(v_training1) %*% inv_SIGMA_sp_kro_inv_SIGMA_tmp %*% v_training1)/2 + lambda_training1
      sigma1_sq_training1 = invgamma::rinvgamma(1,sigma1_sq_a_training1,rate = sigma1_sq_lambda_training1) 
      #Sigma2^2 posterior
      sigma2_sq_a_training1 = a_training1 + n_training1/2
      temp_vec = y_training1 -X_mat_training1%*%theta_training1 - v_training1
      sigma2_sq_lambda_training1 = (t(temp_vec) %*% (temp_vec))/2 + lambda_training1
      sigma2_sq_training1 = invgamma::rinvgamma(1,sigma2_sq_a_training1,rate = sigma2_sq_lambda_training1)
      
      if(i> burn_period_training1 & i%%diff_in_random_draws_training1 ==0){
        v_sample_training1 = cbind(v_sample_training1,v_training1)
        theta_sample_training1 = cbind(theta_sample_training1,theta_training1)
        sigma1_sq_sample_training1 = c(sigma1_sq_sample_training1,sigma1_sq_training1)
        sigma2_sq_sample_training1 = c(sigma2_sq_sample_training1,sigma2_sq_training1)
      }
    }
    #estimation of y with the posterior parameters
    estimated_y_training1 = X_mat_training1%*%rowMeans(theta_sample_training1) + rowMeans(v_sample_training1) + 
      rnorm(n_training1,mean = 0, sd= sqrt(mean(sigma2_sq_sample_training1)))
    mse_training1 = (sum((y_training1 - estimated_y_training1)^2))/n_training1
    
    #Prediction
  }
}
