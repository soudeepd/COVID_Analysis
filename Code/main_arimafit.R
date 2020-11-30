# Fitting Auto ARIMA model

library(forecast)
library(tidyverse)

# Use training1, testing and df generated from main_dataprep.R

# get list of states
all_states = unique(df$state)

# create new columns for ARIMA
df = df %>% 
  mutate(
    arima_trained = as.numeric(week <= test_wk_min),
    arima_fit = 999
  ) %>%
  ungroup()

# assign the regressor names for the arima models
regset = c("time","time_sq","prev_log_new_death","prev_log_prevalence_country")

# fit arima to every state separately
all_arima_models = list()
for (i in 1:length(all_states)){
  idx1 = which(df$state == all_states[i] & df$arima_trained == 1)
  idx2 = which(df$state == all_states[i] & df$arima_trained == 0)
  reg_train = cbind(rep(1,length(idx1)),df[idx1,regset])
  reg_test = cbind(rep(1,length(idx2)),df[idx2,regset])
  colnames(reg_train) = colnames(reg_test) = c("beta0",regset)
  fit_auto = auto.arima(df$log_prevalence[idx1],
                        max.p = 3,max.P = 3,max.q = 3,max.Q = 3,max.d = 1,max.D = 0,
                        seasonal = FALSE,
                        trace = FALSE,
                        approximation = FALSE,
                        allowdrift = FALSE,
                        allowmean = FALSE,
                        stepwise = TRUE,
                        biasadj = FALSE,
                        ic = "aicc",
                        lambda = NULL,
                        xreg = as.matrix(reg_train))
  all_arima_models[[i]] = fit_auto
  y_pred = forecast::forecast(fit_auto,xreg = as.matrix(reg_test))
  df$arima_fit[idx1] = fitted(fit_auto)
  df$arima_fit[idx2] = y_pred$mean
}

# find out the errors
df = df %>%
  mutate(
    arima_error = arima_fit - log_prevalence,
    arima_mape = 100*abs(arima_error/log_prevalence),
    arima_smape = 100*abs(arima_error/(abs(log_prevalence)+abs(arima_fit))),
    arima_mse = arima_error^2#,
    # our_error = log_prev_model_values - log_prevalence,
    # our_mape = 100*abs(our_error/log_prevalence),
    # our_smape = 100*abs(our_error/(abs(log_prev_model_values) + abs(log_prevalence))),
    # our_mse = our_error^2
  )

# check state-wise prediction error
error_summary = df %>%
  filter(arima_trained == 0) %>%
  group_by(state) %>%
  summarize(
    arima_MAE = mean(abs(arima_error)),
    arima_MAPE = mean(arima_mape),
    arima_SMAPE = mean(arima_smape),
    arima_MSE = mean(arima_mse)
  ) %>%
  ungroup() %>%
  mutate_if(is.numeric,round,digits = 3) 
View(error_summary)

# check total errors
mean(error_summary$arima_MAPE)
mean(error_summary$arima_MSE)
mean(error_summary$arima_SMAPE)


# see the model performance for a state
df %>%
  filter(state == 'Wisconsin') %>%
  ggplot(aes(x = week)) +
  geom_line(aes(y = log_prevalence,col = "Actual")) +
  geom_line(aes(y = arima_fit,col = "ARIMA")) +
  # geom_line(aes(y = log_prev_model_values,col = "Our model")) +
  geom_vline(aes(xintercept = test_wk_min),lty = 2) +
  scale_color_manual(values = c("Actual" = "black","ARIMA" = "red","Our model" =  "blue")) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(size = 20))




#::::: OLD CODES FROM SIDDHARTH
Xreg  = cbind(rep(1, n_full_training2), full_training2$time, full_training2$time_sq) #rep(1, n_full_training2),
y_pred_arima_mat = vector()
for (loc in (1:n_sp_full_training2)) {
  #fitting auto arima
  fit_auto = auto.arima(y_full_training2[((loc - 1)*n_tmp_full_training2+1):(loc*n_tmp_full_training2)],
                        max.p = 7,max.P = 7,max.q = 7,max.Q = 7,max.d = 1,max.D = 0,
                        seasonal = FALSE,
                        trace = FALSE,
                        approximation = FALSE,
                        allowdrift = FALSE,
                        allowmean = FALSE,
                        stepwise = TRUE,
                        biasadj = FALSE,
                        ic = "aicc",
                        lambda = NULL,
                        xreg = as.matrix(Xreg[(((loc - 1)*n_tmp_full_training2+1):(loc*n_tmp_full_training2)),]))
  
  #Forecasting using the fitted model on forecast at a time
  y_pred_arima_vec = vector()
  for (pred_i in (1:n_tmp_testing)) {
    if(pred_i %% n_tmp_testing == 1){
      X_vec_tmp = t(c(1, testing$time[pred_i], testing$time_sq[pred_i]))
    }else{
      X_vec_tmp = t(c(1, testing$time[pred_i], testing$time_sq[pred_i]))
    }
    y_pred_arima = forecast(fit_auto, xreg = X_vec_tmp)
    y_pred_arima_vec = c(y_pred_arima_vec, unname(y_pred_arima$mean[1]))
  }
  y_pred_arima_mat = cbind(y_pred_arima_mat, y_pred_arima_vec)
}
#Mape for predicted values
mape_pred_arima = (sum(abs( (y_testing - y_pred_arima_mat[1:(n_sp_full_training2*n_tmp_testing)])/y_testing )))/n_testing
