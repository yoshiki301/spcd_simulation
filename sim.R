library("stringr")
library("dplyr")
library("config")

source("helper.R")
source("random_sampling.R")
source("EMy2.0.R")

# read config file
config_params <<- config::get(config = "simulation")
model <<- config_params$generative_model
seed <<- config_params$set_seed
sim_num <<- config_params$sim_num
spcd_w <<- config_params$spcd_w
pi_list <<- config_params$pi_list
sample_list <<- config_params$sample_list
rho_list <<- config_params$rho_list
save_dir <<- config_params$save_dir

set.seed(seed = seed)

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

# define result matrix shape (sample * rho)
result_shape <- data.frame(
  matrix(
    0,
    nrow = length(sample_list),
    ncol = length(rho_list)
  )
)
dimnames(result_shape) <- list(sample_list, rho_list)

for (pi in pi_list) {
  basedir <- paste(save_dir, "/", pi, "/", sep = "")
  dir.create(basedir, showWarnings = F, recursive = T)
  
  make_notification(
    header = "---simulation start---",
    sim_num = sim_num,
    pi = pi,
    spcd_w = spcd_w
  )
  
  i <- 1
  while (i <= sim_num) {
  
    print(i)
  
    data_matrix <- data.frame()
    pi_matrix <- result_shape
    delta1_matrix <- result_shape
    delta2_matrix <- result_shape
    effect_matrix <- result_shape
    effect_t_matrix <- result_shape
  
    save_flag <- TRUE
  
    for (rho_12 in rho_list){
      for (patient_size in sample_list){
      
        # generate virtual data
        patient_data <- random_sampling(
          rho_12 = rho_12,
          patient_size = patient_size,
          pi = pi
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
      
          # bind data with response
          patient_data_with_response <- data.frame(
            patient_data,
            response = result$response,
            rho_12 = rho_12,
            patient_size = patient_size
          )
          data_matrix <- rbind(data_matrix, patient_data_with_response)
          
          # calculate overall pi and effect
          estimated_values <- calculate_estimated_values(
            parameters = result$parameters,
            stage1_drug = stage1_drug,
            stage1_placebo = stage1_placebo,
            spcd_w = spcd_w
          )

          s <- as.character(patient_size)
          r <- as.character(rho_12)
          pi_matrix[s, r] <- estimated_values$pi
          delta1_matrix[s, r] <- estimated_values$delta1
          delta2_matrix[s, r] <- estimated_values$delta2nr
          effect_matrix[s, r] <- estimated_values$effect
        }
      
      }
    }

  
    if (save_flag) {
      number <- str_pad(i, 4, pad="0")
  
      data_filepath <- paste(basedir, "data_", number, ".csv", sep = "")
      pi_filepath <- paste(basedir, "pi_", number, ".csv", sep = "")
      delta1_filepath <- paste(basedir, "delta1_", number, ".csv", sep = "")
      delta2_filepath <- paste(basedir, "delta2_", number, ".csv", sep = "")
      effect_filepath <- paste(basedir, "effect_", number, ".csv", sep = "")
  
      write.csv(data_matrix, data_filepath)
      write.csv(pi_matrix, pi_filepath)
      write.csv(delta1_matrix, delta1_filepath)
      write.csv(delta2_matrix, delta2_filepath)
      write.csv(effect_matrix, effect_filepath)
      
      i <- i + 1
    } else {
      print("fail solving")
    }
  }
  
  make_notification(
    header = "---simulation end---",
    sim_num = sim_num,
    spcd_w = spcd_w,
    pi = pi
  )
}

