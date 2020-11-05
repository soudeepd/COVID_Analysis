
set.seed(2020)

# read the statewise weekly data
df = read.csv(file = "~/Documents/GitHub/COVID_Analysis/Code/US_statewise_weekly.csv")

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
