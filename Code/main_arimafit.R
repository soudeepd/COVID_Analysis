#Fitting Auto ARIMA model

library(forecast)

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
