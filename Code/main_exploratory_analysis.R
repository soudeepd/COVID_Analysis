
# load relevant libraries
library(ggplot2)
library(ggpubr)
library(usmap)

# see summary of timeline to find quartiles
summary(df$week)

# divide the data according to quartiles
data_subset1 = df %>%
  filter(week == as.Date("2020-04-06")) %>%
  dplyr::select(state,log_prevalence) 

data_subset2 = df %>%
  filter(week == as.Date("2020-06-22")) %>%
  dplyr::select(state,log_prevalence) 

data_subset3 = df %>%
  filter(week == as.Date("2020-09-07")) %>%
  dplyr::select(state,log_prevalence) 

data_subset4 = df %>%
  filter(week == as.Date("2020-11-23")) %>%
  dplyr::select(state,log_prevalence) 

# set lower and upper limit of every plot
LL = min(df$log_prevalence)
UL = max(df$log_prevalence)

# create four different plots
p1 = plot_usmap(data = data_subset1,values = "log_prevalence",color = "black",exclude = c("AK","HI")) + 
  scale_fill_continuous(name = "log_prevalence",low = "white",high = "black",limits = c(LL,UL),label = scales::comma) + 
  theme(legend.position = "right") +
  ggtitle("6 Apr, 2020") +
  theme(plot.title = element_text(hjust = 0.5,size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20))

p2 = plot_usmap(data = data_subset2,values = "log_prevalence",color = "black",exclude = c("AK","HI")) + 
  scale_fill_continuous(name = "log_prevalence",low = "white",high = "black",limits = c(LL,UL),label = scales::comma) + 
  theme(legend.position = "right") +
  ggtitle("22 Jun, 2020") +
  theme(plot.title = element_text(hjust = 0.5,size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20))

p3 = plot_usmap(data = data_subset3,values = "log_prevalence",color = "black",exclude = c("AK","HI")) + 
  scale_fill_continuous(name = "log_prevalence",low = "white",high = "black",limits = c(LL,UL),label = scales::comma) + 
  theme(legend.position = "right") +
  ggtitle("7 Sep, 2020") +
  theme(plot.title = element_text(hjust = 0.5,size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20))

p4 = plot_usmap(data = data_subset4,values = "log_prevalence",color = "black",exclude = c("AK","HI")) + 
  scale_fill_continuous(name = "log_prevalence",low = "white",high = "black",limits = c(LL,UL),label = scales::comma) + 
  theme(legend.position = "right") +
  ggtitle("23 Nov, 2020") +
  theme(plot.title = element_text(hjust = 0.5,size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20))

# merge them together and save on file
ggpubr::ggarrange(p1,p2,p3,p4) + 
  ggsave(filename = "~/Documents/GitHub/COVID_Analysis/Code/US_map_COVID.eps",
         width = 600,
         height = 400,
         units = "mm",
         dpi = 800)
