library("config")
library("dplyr")
library("mclust")

source("helper.R")

config <<- config::get()
datapath <<- config$datapath
destdir <<- config$destdir
model <<- config$model

# read config file
source(model)
simulation_csv_path <- file.path(destdir, datapath)

seed <- 101
set.seed(seed = seed)

sim_num <- 5000

# fixed parameter
pi <- 0.3
delta1 <- 1.5
delta2_nr <- 1.5
h <- 1 # delta2_r = h * delta2_nr
drug_assign <- 1/3
w <- 0.5

source("class/EMAlgorithm.R")

rho <- c(-0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5)
sample <- c(90, 120, 150, 225, 300, 450, 600)

pi_df <- data.frame(matrix(0, nrow=11, ncol=7))
names(pi_df) <- sample
rownames(pi_df) <- rho

effect_df <- data.frame(matrix(0, nrow=11, ncol=7))
names(effect_df) <- sample
rownames(effect_df) <- rho

#mse_df <- data.frame(matrix(0, nrow=7, ncol=11))
#names(mse_df) <- sample
#rownames(mse_df) <- rho

for (i in 1:sim_num){
  print(i)
  for (rho_12 in rho){
    for (patient_size in sample){
      drug_size <- patient_size*drug_assign
      responder_size <- patient_size*(1 - drug_assign)*pi
      nonresponder_size <- patient_size*(1 - drug_assign)*(1 - pi)
    
      # generate patients data
      dd_group <- rep("DD", round(drug_size))
      rpd_group <- rep("rPD", ceiling(responder_size*drug_assign))
      rpp_group <- rep("rPP", floor(responder_size*(1-drug_assign)))
      npd_group <- rep("nPD", ceiling(nonresponder_size*drug_assign))
      npp_group <- rep("nPP", floor(nonresponder_size*(1-drug_assign)))
    
      group <- c(dd_group, rpd_group, rpp_group, npd_group, npp_group)
    
      GenerativeModel <- NormalGenerativeModel$new(
        pi=pi,
        delta1=delta1,
        delta2_nr=delta2_nr,
        h=h,
        rho_12=rho_12
      )
    
      random_sample <- lapply(group, GenerativeModel$generateSample)
      y01 <- sapply(random_sample, "[[", 1)
      y1 <- sapply(random_sample, "[[", 2)
      y2 <- sapply(random_sample, "[[", 3)
    
      patient_df <- data.frame(
        y01=y01,
        y1=y1,
        y2=y2,
        group=group
      )
      
      #EMModel <- EMAlgorithmModel$new(
      #  df = patient_df,
      #  eps = 0.001,
      #  w = w
      #)
    
      #EMModel$emSimulation(patient_df)
    
      placebo_df <- patient_df[patient_df$group != "DD", ]
      responder <- lapply(placebo_df$group, labeling_responder)
  
    
      data_columns <- c("y01", "y1", "y2")
      placebo_data <- placebo_df[, data_columns]
    
      clustering <- Mclust(placebo_data, G=2, modelNames="EEE", verbose=FALSE)
      classification <- clustering$classification
    
      estimated_pi <- table(classification)[1] / length(classification)
    
      pi_df[as.character(rho_12),as.character(patient_size)] <- pi_df[as.character(rho_12),as.character(patient_size)] + estimated_pi
    
      estimated_label <- mapply(labeling_estimation, placebo_df$group, classification)
    
      # estimate treatment effect
      mu11 <- mean(patient_df[patient_df$group == "DD", "y1"])
      #sigma11 <- var(patient_df[patient_df$group == "DD", "y1"])
      mu10 <- estimated_pi*mean(placebo_df[classification == 1, "y1"]) + (1-estimated_pi)*mean(placebo_df[classification == 2, "y1"])
      #sigma10 <- (estimated_pi^2)*var(placebo_df[classification == 1, "y1"]) + ((1-estimated_pi)^2)*var(placebo_df[classification == 2, "y1"])
      hat_delta1 <- mu11 - mu10
      #sigma_delta1 <- sigma11 + sigma10
      hat_delta2_nr <- mean(placebo_df[estimated_label == "nPD", "y2"]) - mean(placebo_df[estimated_label == "nPP", "y2"])
      #sigma_delta2_nr <- var(placebo_df[estimated_label == "nPD", "y2"]) + var(placebo_df[estimated_label == "nPP", "y2"])
    
      hat_delta_w <- w*hat_delta1 + (1-w)*hat_delta2_nr
      #t <- delta_w / sqrt((w^2)*sigma_delta1+((1-w)^2)*sigma_delta2_nr)

      effect_df[as.character(rho_12),as.character(patient_size)] <- effect_df[as.character(rho_12),as.character(patient_size)] + hat_delta_w
    }
  }
  #print(pi_df)
  #print(effect_df)
}

pi_df <- pi_df / sim_num
effect_df <- effect_df / sim_num

write.csv(pi_df, "./normal/pi.csv")
write.csv(effect_df, "./normal/effect.csv")
