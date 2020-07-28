library("config")

config <<- config::get()
datapath <<- config$datapath
destdir <<- config$destdir

simulation_csv_path <- file.path(destdir, datapath)
patient_df <- read.csv(simulation_csv_path)

source("class/EMAlgorithm.R")

EMModel <- EMAlgorithmModel$new(
  df = patient_df,
  eps = 0.001
)

EMModel$emSimulation(patient_df)

EMModel$printProcess()