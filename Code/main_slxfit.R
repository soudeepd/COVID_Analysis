
# load relevant libraries
library(spdep)
library(sf)

# the data needs to be prepared
source("~/Documents/GitHub/COVID_Analysis/Code/main_dataprep.R")

# take the training data and define it as a spatial object
spdf <- full_training2
coordinates(spdf) <- ~ long + lat

# create the neighbours list for the SLX code
neib <- knn2nb(knearneigh(coordinates(spdf),longlat = TRUE))
lw <- nb2listw(neib,style = "B") # here see documentation, you can use W/B/C/U/minimax/S

# define the regressors and then run the model
regset = c("time","time_sq","prev_log_new_death")
model_formula = paste("log_prevalence~I(log(population))+",paste(regset,collapse = "+"))
slxmodel <- lagsarlm(formula(model_formula),data = spdf,listw = lw)
summary(slxmodel)

# see predictions: first use the whole data to define the spatial object
# and then see predictions only for the training part
preddf = df
coordinates(preddf) <- ~ long + lat
predneib <- knn2nb(knearneigh(coordinates(preddf),longlat = TRUE))
predlw <- spdep::nb2listw(predneib,style = "B")

slxpredict <- predict.sarlm(slxmodel,newdata = preddf,listw = predlw)
df$slxfit = slxpredict

# error summary for the testing part
slxerror = df %>% 
  filter(week > test_wk_min) %>%
  dplyr::select(week,state,log_prevalence,slxfit) %>%
  mutate(
    error = slxfit - log_prevalence,
    mape = 100*abs(error)/abs(log_prevalence)
  ) %>%
  mutate_if(is.numeric,round,digits = 3) 

mean(slxerror$mape)
slxerror %>%
  group_by(state) %>%
  summarize(
    mae = mean(abs(error)),
    mse = mean(error^2),
    mape = mean(mape)
  ) %>%
  ungroup() %>%
  mutate_if(is.numeric,round,digits = 3) %>%
  View()


m2 = lm(formula(model_formula),data = spdf)
summary(m2)
df$slxfit = predict(m2,df)
