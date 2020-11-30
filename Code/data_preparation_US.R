
library(tidyverse)

#:::::::::::::: Read, clean and prepare data

# assign the directory
filepath = "~/Documents/GitHub/COVID-19/csse_covid_19_data/csse_covid_19_time_series/"

# read the data files
df1 = read.csv(paste(filepath,"time_series_covid19_confirmed_US.csv",sep = ""))
df2 = read.csv(paste(filepath,"time_series_covid19_deaths_US.csv",sep = ""))

# find the maximum date
lastcol = colnames(df2)[ncol(df2)]
if (colnames(df1)[ncol(df1)] == lastcol){
  string = paste(substr(lastcol,2,nchar(lastcol)),"20",sep = "")
  maxdate = as.Date(string,format = "%m.%d.%Y")
}

# assign location_ids and dates, extract location info
loc_cols = c("Province_State","Country_Region","Lat","Long_","Combined_Key")
ids = unlist(lapply(seq(1,nrow(df1),by = 1),FUN = function(x) paste("ID",x,sep = "")))
dates = seq.Date(from = as.Date("2020-01-22"),to = maxdate,by = 1)
location_info = cbind(ids,df1[,loc_cols])
names(location_info) = c("ID","state","country","lat","long","key")

# extract data on confirmed cases
date_cols = colnames(df1)[which(substr(colnames(df1),1,1) == 'X')]
covid_conf = data.frame(dates,t(df1[,date_cols]),stringsAsFactors = T)
colnames(covid_conf) = c("dates",as.character(location_info$ID))
rownames(covid_conf) = 1:nrow(covid_conf)

# extract data on death cases
covid_dead = data.frame(dates,t(df2[,date_cols]),stringsAsFactors = T)
colnames(covid_dead) = c("dates",as.character(location_info$ID))
rownames(covid_dead) = 1:nrow(covid_dead)

# merge info on confirmed and death 
N = dim(covid_conf)[1]
conf = c()
de = c()
loc = c()
for (i in 1:nrow(location_info)){
  colno = which(colnames(covid_conf) == location_info$ID[i])
  conf = c(conf,covid_conf[,colno])
  de = c(de,covid_dead[,colno])
  tmp = location_info[i,]
  tmp = tmp %>% slice(rep(1:n(), each = N))
  loc = rbind(loc,tmp)
}
usdata = data.frame(date = as.Date(rep(dates,nrow(location_info))))
usdata = data.frame(usdata,loc,conf,de)
colnames(usdata)[8:9] = c("confirmed","death")

# change type of columns as required
usdata = usdata %>% mutate_if(is.factor,as.character)
usdata$date = as.Date(usdata$date)

# filter out no lat-long info cases
usdata = usdata %>% filter(lat != 0 & long != 0)

# merge population data
usdata = merge(usdata,df2[,c("Lat","Long_","Population")],by.x = c("lat","long"),by.y = c("Lat","Long_")) 
colnames(usdata) = tolower(colnames(usdata))
usdata = usdata %>% filter(population > 0)

# order by date
usdata = usdata[order(usdata$date,usdata$lat,usdata$long),]

# states to be excluded - not in contiguous US
exclusion = c("American Samoa","Alaska","Guam","Hawaii","Northern Mariana Islands","Puerto Rico","Virgin Islands")

# combine statewise data into one data frame
statewise = usdata %>%
  filter(! state %in% exclusion) %>%
  group_by(date,state,country) %>%
  summarize(
    lat = median(lat),
    long = median(long),
    confirmed = sum(confirmed), # adding total numbers for different locations within same state
    death = sum(death),
    population = sum(population)
  ) %>%
  ungroup()

# add week column (last Monday) to the data
lastmon2 <- function(x) x - as.numeric(x-1+4)%%7
statewise = statewise %>%
  mutate(
    week = lastmon2(date)
  )

# make weekly data and create new variables
statewise_weekly = statewise %>%
  group_by(country,state,lat,long,week,population) %>%
  summarize(
    confirmed = max(confirmed), # getting total cases as the maximum in that week
    death = max(death)
  ) %>%
  ungroup() %>%
  group_by(country,state,lat,long) %>%
  mutate(
    prevalence = confirmed/population,
    log_prevalence = log(prevalence + 0.1), # avoiding error for 0 cases
    prev_day_total = dplyr::lag(confirmed,default = min(confirmed)),
    new_cases = confirmed - prev_day_total,
    incidence = ifelse(new_cases > 0,new_cases/population,0),
    log_incidence = log(incidence + 0.1), # avoiding error for 0 cases
    prev_day_actual = dplyr::lag(new_cases,default = min(confirmed))
  ) %>%
  ungroup()

# write the data in a csv file
write.csv(x = statewise_weekly,file = "~/Documents/GitHub/COVID_Analysis/Code/US_statewise_weekly.csv",row.names = F)
