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