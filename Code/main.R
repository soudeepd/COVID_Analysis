
# read the statewise weekly data
df = read.csv(file = "~/Documents/GitHub/COVID_Analysis/Code/US_statewise_weekly.csv")

# check class of the columns in the data
sapply(df,class) 

# change week column if needed
df$week = as.Date(df$week)
