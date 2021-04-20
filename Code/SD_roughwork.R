
library(forecast)
library(tidyverse)
library(spdep)
library(sf)
library(MCS)
library(spatialreg)

setwd("~/OneDrive - Indian Institute of Management/IIMB/Independent-Study/2020-Siddharth/Datasets/")

df = read.csv("Datasets/US_statewise_weekly_training2_2weeks.csv")
df = df %>% 
  mutate_if(is.factor,as.character) %>% 
  mutate(
    week = as.Date(week),
    log_population = log(population)
  ) %>%
  arrange(country,state,week) 

maxdate = max(subset(df,is_trained == 1)$week)

# fit the SLX model
spdf <- df %>% filter(week <= maxdate)
coordinates(spdf) <- ~ long + lat

# create the neighbours list for the SLX code
neib <- knn2nb(knearneigh(coordinates(spdf),longlat = TRUE))
lw <- nb2listw(neib,style = "W") # here see documentation, you can use W/B/C/U/minimax/S

# define the regressors and then run the model
regset = c("time","time_sq","prev_log_new_death","log_population")
model_formula = paste("log_prevalence~",paste(regset,collapse = "+"))
slxmodel <- spatialreg::lagsarlm(formula(model_formula),data = spdf,listw = lw)

# see predictions: first use the whole data to define the spatial object
# and then see predictions only for the training part
preddf = df
coordinates(preddf) <- ~ long + lat
predneib <- knn2nb(knearneigh(coordinates(preddf),longlat = TRUE))
predlw <- spdep::nb2listw(predneib,style = "W")

slxpredict <- spatialreg::predict.sarlm(slxmodel,newdata = preddf,listw = predlw)
df$slxfit = slxpredict


# fitting arima model
# get list of states
all_states = unique(df$state)

# create new columns for ARIMA
df = df %>% 
  mutate(
    arima_trained = as.numeric(week <= maxdate),
    arima_fit = 999
  ) 

# define the regressors and then run the model
regset = c("time","time_sq","prev_log_new_death")

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

# create the error columns
df = df %>%
  mutate(
    arima_error = arima_fit - log_prevalence,
    slx_error = slxfit - log_prevalence,
    our_error = log_prev_model_values - log_prevalence
  )

# only the testing part errors
test_error = df %>% filter(week > maxdate) 
test_error = test_error[complete.cases(test_error),]

# summarize test errors by states
test_error %>%
  group_by(state) %>%
  summarize(
    arima_mape = 100*mean(abs(arima_error/log_prevalence)),
    slx_mape = 100*mean(abs(slx_error/log_prevalence)),
    our_mape = 100*mean(abs(our_error/log_prevalence))
  ) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(
    least_error = min(arima_mape,slx_mape,our_mape),
    best_method = ifelse(least_error == arima_mape,"arima",ifelse(least_error == slx_mape,"slx","ours"))
  ) %>% 
  mutate_if(is.numeric,round,digits = 3) -> state_summary#%>%
  View()


# MCS procedure
testvals = test_error$log_prevalence
predvals = data.frame(#arima = test_error$arima_fit,
                      slx = test_error$slxfit,
                      ours = test_error$log_prev_model_values)
mse_rv = abs(predvals - testvals)^2
mae_rv = abs(predvals - testvals)
MCSprocedure(Loss = mse_rv,alpha = 0.1,statistic = "TR")
MCSprocedure(Loss = mae_rv,alpha = 0.1,statistic = "TR")

test_error %>% 
  dplyr::select(state,week,log_prevalence,log_prev_model_values,arima_error,slx_error,our_error) %>% 
  mutate_if(is.numeric,round,digits=3) %>% 
  View()

# view the plot
df %>%
  filter(state == 'Illinois') %>%
  ggplot(aes(x = week)) +
  geom_line(aes(y = log_prevalence,col = "Actual")) +
  geom_line(aes(y = log_prev_model_values,col = "Modeled")) +
  scale_color_manual(values = c("Actual" = "black","Modeled" = "red")) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(size = 20)) 
  

# error curve for different weeks
allerrors = data.frame(test_wks = c(1:12),
                       training1_error = runif(12,2,4),
                       training2_error = runif(12,2,4))
allerrors %>%
  ggplot(aes(x = test_wks)) +
  geom_line(aes(y = training1_error,lty = "case 1")) + 
  geom_line(aes(y = training2_error,lty = "case 2")) + 
  scale_x_continuous(breaks = c(2,4,6,8,10,12)) +
  xlab("test period (weeks)") +
  ylab("overall MAPE") +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(size = 20)) 


df = read.csv("US_statewise_weekly_training2_withresiduals.csv")
