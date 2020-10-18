library("dplyr")

source("helper.R")
source("random_sampling.R")
source("EMy2.0.R")

seed <- 101
set.seed(seed = seed)

sim_num <- 1100

# fixed paramters
spcd_w <- 0.5

# prepare parameters of EM estimation
em_parameters <- c(
  14, 0.2, 49,      # drug stage 1
  8.0, 0.2, 49,     # drug stage 2
  16, 0.1, 4,       # responder stage 1
  9.6, 0.1, 1.5, 4, # responder stage 2
  9.7, 0.1, 4,      # nonresponder stage 1
  5.8, 0.1, 1.5, 4  # nonresponder stage 2
)
initial_pi <- 0.5
m01 <- 1.5
m02 <- 1.5

# variable parameter
rho <- c(-0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5)
sample <- c(90, 300, 600, 1200)

# define result matrix shape
result_shape <- data.frame(
  matrix(
    0,
    nrow = length(sample),
    ncol = length(rho)
  )
)
dimnames(result_shape) <- list(sample, rho)

basedir <- "./sim1/"

for (i in 1:sim_num){
  
  print(i)
  
  pi_matrix <- result_shape
  effect_matrix <- result_shape
  mu10_matrix <- result_shape
  mu11_matrix <- result_shape
  delta1_matrix <- result_shape
  delta2_matrix <- result_shape
  
  save_flag <- TRUE
  
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

      # EM estimation
      # the latter 27 parameters are dummies, so set to 0.5
      initial_parameters <- c(em_parameters, initial_pi, m01, m02, rep(0.5,27))
      result <- doEMy(
        d = patient_data,
        s = initial_parameters
      )
      
      if (result == "solve error"){
        save_flag <- FALSE
      } else {
        stage1_drug <- patient_data[patient_data["g1"]==1, "y01"]
        stage1_placebo <- patient_data[patient_data["g1"]==0, "y01"]
      
        estimated_values <- calculate_estimated_values(
          parameters = result$parameters,
          stage1_drug = stage1_drug,
          stage1_placebo = stage1_placebo,
          spcd_w = spcd_w
        )

        s <- as.character(patient_size)
        r <- as.character(rho_12)
        pi_matrix[s, r] <- estimated_values$pi
        effect_matrix[s, r] <- estimated_values$delta_w
        mu10_matrix[s, r] <- estimated_values$mu10
        mu11_matrix[s, r] <- estimated_values$mu11
        delta1_matrix[s, r] <- estimated_values$delta1
        delta2_matrix[s, r] <- estimated_values$delta2_nr
      }
      
    }
  }
  
  if (save_flag) {
    number <- as.character(i)
  
    pi_filepath <- paste(basedir, "pi_", number, ".csv", sep = "")
    effect_filepath <- paste(basedir, "effect_", number, ".csv", sep = "")
    mu10_filepath <- paste(basedir, "mu10_", number, ".csv", sep = "")
    mu11_filepath <- paste(basedir, "mu11_", number, ".csv", sep = "")
    delta1_filepath <- paste(basedir, "delta1_", number, ".csv", sep = "")
    delta2_filepath <- paste(basedir, "delta2_", number, ".csv", sep = "")
  
    write.csv(pi_matrix, pi_filepath)
    write.csv(effect_matrix, effect_filepath)
    write.csv(mu10_matrix, mu10_filepath)
    write.csv(mu11_matrix, mu11_filepath)
    write.csv(delta1_matrix, delta1_filepath)
    write.csv(delta2_matrix, delta2_filepath)
  } else {
    print("fail solving")
  }
}