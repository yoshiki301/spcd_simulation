library("tidyverse")
library("GGally")
library("config")

config <<- config::get()
datapath <<- config$datapath
destdir <<- config$destdir

simulation_csv_path <- file.path(destdir, datapath)
patient_data <- read.csv(simulation_csv_path)

# summarise
patient_data %>%
  group_by(group) %>%
  summarise(count=n(), across(y01:y2, list(mean=mean, sd=sd)))

# check correlation and visualize
ggpairs(select(patient_data, -X), aes(color = group, alpha = 0.7))