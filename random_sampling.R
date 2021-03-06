library("config")

# read config file
config_params <<- config::get(config = "generate")
model <<- config_params$generative_model
seed <<- config_params$set_seed
source(model)
set.seed(seed = seed)

random_sampling = function(rho_12, patient_size, pi) {

  # fixed parameter
  delta1 <- 1.5
  delta2_nr <- 1.5
  h <- 1 # delta2_r = h * delta2_nr
  drug_assign <- 1/3

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

  sample <- lapply(group, GenerativeModel$generateSample)
  y01 <- sapply(sample, "[[", 1)
  y1 <- sapply(sample, "[[", 2)
  y2 <- sapply(sample, "[[", 3)

  patient_data <- data.frame(
    y01 = y01,
    y1 = y1,
    y2 = y2,
    group = group
  )
  
  return(patient_data)
}
