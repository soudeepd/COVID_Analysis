
# load necessary libraries
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
df = read.csv(file = "~/Documents/GitHub/COVID_Analysis/Code/US_statewise_weekly.csv")
df = read.csv(file = "D:/Term5/IS&RA/gitrepo/COVID_Analysis/Code/US_statewise_weekly.csv")

# check class of the columns in the data and change their class if needed
sapply(df,class) 
df = df %>% 
  mutate_if(is.factor,as.character) %>% 
  mutate(
    week = as.Date(week)
  ) %>%
  arrange(country,state,week) %>%             # arranging first by state then by week
  group_by(state) %>%                         # creating new regressor columns
  mutate(
    time = c(1:n())/n(),
    time_sq = time^2,
    prev_log_prevalence = lag(log_prevalence,1,default = log_prevalence[1]),
    log_death = log(death + 0.1),
    prev_log_death = lag(log_death,1,default = log_death[1]),
    new_death = death - lag(death,1,default = 0),
    new_death = ifelse(new_death >= 0,new_death,0),    # correcting one row where new_death<0
    log_new_death = log(new_death + 0.1),
    prev_log_new_death = lag(log_new_death,1,default = log_new_death[1])
  ) %>%
  ungroup() %>%
  group_by(week) %>%
  mutate(
    prevalence_country = sum(confirmed)/sum(population),
    log_prevalence_country = log(prevalence_country + 0.1),
    prev_log_prevalence_country = lag(log_prevalence_country,1,default = log_prevalence_country[1])
  ) %>%
  ungroup()


# case1: use all but two random locations as training and last four weeks as test set 
# case2: use all locations for training and last four weeks as test set 
test_loc = sample(unique(df$state),2)
test_wk_min = max(df$week) - 56
validation_wk_min = max(df$week) - 84
training1 = df %>% filter(! state %in% test_loc & week <= validation_wk_min)
training2 = df %>% filter(week <= validation_wk_min)
validation1 = df %>% filter(! state %in% test_loc & week > validation_wk_min & week <= test_wk_min)
validation2 = df %>% filter( week > validation_wk_min & week <= test_wk_min)
testing = df %>% filter(week > test_wk_min)
full_training1 = df %>% filter(! state %in% test_loc & week <= test_wk_min)
full_training2 = df %>% filter(week <= test_wk_min)



#:::::: PREVIOUS CODES
# 
# #Adding columns to Data
# n = length(df$country)
# n_sp = length(unique(df$state))
# n_tmp = n/n_sp
# 
# time_vec = vector()
# time_sq_vec = vector()
# y_st_vec = vector()
# for (j in (1:n)) {
#   if(j %% n_tmp==0){
#     temp_var = 1
#     y_st = df$log_prevalence[j-1]
#   }else if(j %% n_tmp == 1){
#     temp_var = (j %% n_tmp)/n_tmp
#     y_st = 0
#   }else{
#     temp_var = (j %% n_tmp)/n_tmp
#     y_st = df$log_prevalence[j-1]
#   }
#   time_vec = c(time_vec, temp_var)
#   time_sq_vec = c(time_sq_vec, temp_var^2)
#   y_st_vec = c(y_st_vec, y_st)
# }
# 
# df$time = time_vec
# df$time_sq = time_sq_vec
# df$prev_log_prevalence = y_st_vec
# 
