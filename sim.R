library("dplyr")
# library("mclust")

source("helper.R")
source("random_sampling.R")
# source("class/EMAlgorithm.R")
source("EMy2.0.R")

seed <- 101
set.seed(seed = seed)

sim_num <- 100

# fixed paramters
spcd_w <- 0.6

# prepare parameters of EM estimation
em_parameters <- rep(0.5,20)
initial_pi <- 0.4
m01 <- 0.5
m02 <- 0.5

# variable parameter
rho <- c(-0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5)
sample <- c(90, 300, 600, 1200)

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
      
      # generate virtual data
      patient_data <- random_sampling(
        rho_12 = rho_12,
        patient_size = patient_size
      )
      
      patient_data["y02"] <- patient_data["y01"] + patient_data["y1"]
      encode_g <- apply(patient_data["group"], 1, encode_group)
      patient_data["g1"] <- encode_g[1,]
      patient_data["g2"] <- encode_g[2,]
      
      ordered_data <- patient_data[order(patient_data["g1"]),]
      
      np <- sum(ordered_data["g1"] == 0)
      
      # EM estimation
      # the latter 27 parameters are dummies, so set to 0.5
      initial_parameters <- c(em_parameters, initial_pi, m01, m02, rep(0.5,27))
      result <- doEMy(
        d = ordered_data,
        s = initial_parameters
      )
      
      print(result)
    
      # pi_df[as.character(rho_12),as.character(patient_size)] <- pi_df[as.character(rho_12),as.character(patient_size)] + estimated_pi
    
      # estimated_label <- mapply(labeling_estimation, placebo_df$group, classification)
    
      # estimate treatment effect
      # mu11 <- mean(patient_df[patient_df$group == "DD", "y1"])
      # sigma11 <- var(patient_df[patient_df$group == "DD", "y1"])
      # mu10 <- estimated_pi*mean(placebo_df[classification == 1, "y1"]) + (1-estimated_pi)*mean(placebo_df[classification == 2, "y1"])
      # sigma10 <- (estimated_pi^2)*var(placebo_df[classification == 1, "y1"]) + ((1-estimated_pi)^2)*var(placebo_df[classification == 2, "y1"])
      # hat_delta1 <- mu11 - mu10
      # sigma_delta1 <- sigma11 + sigma10
      # hat_delta2_nr <- mean(placebo_df[estimated_label == "nPD", "y2"]) - mean(placebo_df[estimated_label == "nPP", "y2"])
      # sigma_delta2_nr <- var(placebo_df[estimated_label == "nPD", "y2"]) + var(placebo_df[estimated_label == "nPP", "y2"])
    
      # hat_delta_w <- w*hat_delta1 + (1-w)*hat_delta2_nr
      # t <- delta_w / sqrt((w^2)*sigma_delta1+((1-w)^2)*sigma_delta2_nr)

      # effect_df[as.character(rho_12),as.character(patient_size)] <- effect_df[as.character(rho_12),as.character(patient_size)] + hat_delta_w
    }
  }
  #print(pi_df)
  #print(effect_df)
}

#pi_df <- pi_df / sim_num
#effect_df <- effect_df / sim_num

#write.csv(pi_df, "./normal/pi.csv")
#write.csv(effect_df, "./normal/effect.csv")
