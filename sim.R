library("dplyr")

source("helper.R")
source("random_sampling.R")
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
    }
  }
}
