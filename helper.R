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