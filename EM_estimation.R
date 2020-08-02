library("config")

config <<- config::get()
datapath <<- config$datapath
destdir <<- config$destdir

simulation_csv_path <- file.path(destdir, datapath)
patient_df <- read.csv(simulation_csv_path)

#source("class/EMAlgorithm.R")

#EMModel <- EMAlgorithmModel$new(
#  df = patient_df,
#  eps = 0.001
#)

#EMModel$emSimulation(patient_df)

#EMModel$printProcess()

library("mclust")

placebo_df <- patient_df[patient_df$group != "DD", ]

responder <- c()
for (i in 1:nrow(placebo_df)) {
  g <- placebo_df$group[i]
  if (g == "rPD" || g == "rPP") {
    responder <- append(responder, "responder")
  }
  else {
    responder <- append(responder, "nonresponder")
  }
}

data_columns <- c("y01", "y1", "y2")
placebo_data <- placebo_df[, data_columns]

result <- Mclust(placebo_data, G=2, modelNames="EEE")
classification <- result$classification

estimated_pi <- table(classification)[1] / length(classification)
print(estimated_pi)

