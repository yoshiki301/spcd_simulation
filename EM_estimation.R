library("config")

config <<- config::get()
datapath <<- config$datapath
destdir <<- config$destdir

simulation_csv_path <- file.path(destdir, datapath)
patient_df <- read.csv(simulation_csv_path)

source("class/EMAlgorithm.R")

EMModel <- EMAlgorithmModel$new(
  df = patient_df,
  eps = 0.0001
)

EMModel$emSimulation(patient_df)

EMModel$printProcess()

EMModel$printQlist()

# em algorithm by mclust
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

estimated_label <- c()
for (i in 1:nrow(placebo_df)) {
  truth_g <- placebo_df$group[i]
  estimated_g <- classification[i]
  if (estimated_g == 1) {
    if (truth_g == "rPD" || truth_g == "nPD") {
      estimated_label <- append(estimated_label, "rPD")
    }
    else {
      estimated_label <- append(estimated_label, "rPP")
    }
  }
  else {
    if (truth_g == "rPD" || truth_g == "nPD") {
      estimated_label <- append(estimated_label, "nPD")
    }
    else {
      estimated_label <- append(estimated_label, "nPP")
    }
  }
}

# estimate treatment effect
w <- 0.5
mu11 <- mean(patient_df[patient_df$group == "DD", "y1"])
sigma11 <- var(patient_df[patient_df$group == "DD", "y1"])
mu10 <- estimated_pi*mean(placebo_df[classification == 1, "y1"]) + (1-estimated_pi)*mean(placebo_df[classification == 2, "y1"])
sigma10 <- (estimated_pi^2)*var(placebo_df[classification == 1, "y1"]) + ((1-estimated_pi)^2)*var(placebo_df[classification == 2, "y1"])
delta1 <- mu11 - mu10
sigma_delta1 <- sigma11 + sigma10
delta2_nr <- mean(placebo_df[estimated_label == "nPD", "y2"]) - mean(placebo_df[estimated_label == "nPP", "y2"])
sigma_delta2_nr <- var(placebo_df[estimated_label == "nPD", "y2"]) + var(placebo_df[estimated_label == "nPP", "y2"])

delta_w <- w*delta1 + (1-w)*delta2_nr
t <- delta_w / sqrt((w^2)*sigma_delta1+((1-w)^2)*sigma_delta2_nr)

