labeling_responder = function(group) {
  if (group == "rPD" || group == "rPP") {
    return("responder")
  }
  else {
    return("nonresponder")
  }
}

labeling_estimation = function(truth_g, estimated_g){
  if (estimated_g == 1) {
    if (truth_g == "rPD" || truth_g == "nPD") {
      return("rPD")
    }
    else {
      return("rPP")
    }
  }
  else {
    if (truth_g == "rPD" || truth_g == "nPD") {
      return("nPD")
    }
    else {
      return("nPP")
    }
  }
}

encode_group = function(group) {
  # encode group to 2 binary labels (0 or 1)
  if (group == "DD") {
    g1 <- 1
    g2 <- 1
  } else {
    g1 <- 0
    if (group == "rPD" || group == "nPD") {
      g2 <- 1
    } else {
      g2 <- 0
    }
  }
  
  res <- c(g1, g2)
  return(res)
}

calculate_estimated_values = function(
  parameters,
  stage1_drug,
  stage1_placebo,
  spcd_w
) {
  # calculate pi and effect
  estimated_pi <- parameters[21]
  effect <- parameters[48]
  effect_var <- parameters[49]
  
  #hat_mu11 <- parameters[1] + parameters[2]*mean(stage1_drug)
  #hat_mu101 <- parameters[7] + parameters[8]*mean(stage1_placebo)
  #hat_mu102 <- parameters[14] + parameters[15]*mean(stage1_placebo)
  #hat_mu10 <- estimated_pi*hat_mu101 + (1-estimated_pi)*hat_mu102
  #hat_delta1 <- hat_mu11 - hat_mu10
  #hat_delta2_nr <- parameters[19]
  #hat_delta_w <- spcd_w*hat_delta1 + (1-spcd_w)*hat_delta2_nr
  
  result_values <- list(
    estimated_pi,
    effect,
    effect_var
  )
  
  names(result_values) <- c(
    "pi",
    "effect",
    "effect_var"
  )
  
  return(result_values)
}