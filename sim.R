library("dplyr")

source("helper.R")
source("random_sampling.R")
source("EMy2.0.R")
source("line_notify.R")

seed <- 101
set.seed(seed = seed)

sim_num <- 1020

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
rho <- c(-0.5, 0, 0.5)
sample <- c(90, 300, 600)
pi <- 0.3

# define result matrix shape
result_shape <- data.frame(
  matrix(
    0,
    nrow = length(sample),
    ncol = length(rho)
  )
)
dimnames(result_shape) <- list(sample, rho)

basedir <- "./sim7/"

# start notification
datetime <- Sys.time()
notification_text <- paste(datetime, "start", pi, basedir)
r <- notify(notification_text)

for (i in 1:sim_num){
  
  print(i)
  
  pi_matrix <- result_shape
  effect_matrix <- result_shape
  effect_t_matrix <- result_shape
  
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
        effect_matrix[s, r] <- estimated_values$effect
        effect_t_matrix[s, r] <- (estimated_values$effect / estimated_values$effect_var)
      }
      
    }
  }
  
  if (save_flag) {
    number <- as.character(i)
  
    pi_filepath <- paste(basedir, "pi_", number, ".csv", sep = "")
    effect_filepath <- paste(basedir, "effect_", number, ".csv", sep = "")
    t_filepath <- paste(basedir, "t_", number, ".csv", sep = "")
  
    write.csv(pi_matrix, pi_filepath)
    write.csv(effect_matrix, effect_filepath)
    write.csv(effect_t_matrix, t_filepath)
  } else {
    print("fail solving")
  }
}

# notification
datetime <- Sys.time()
notification_text <- paste(datetime, "end", pi, basedir)
r <- notify(notification_text)