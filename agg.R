library(tidyverse)

target_dirs <- c("./sim3/", "./sim4/", "./sim5/", "./sim6/", "./sim7/")
true_pi <- c(0.1, 0.5, 0.7, 0.9, 0.3)
true_effect <- rep(1.5, 5)

calculate_csv_statistics = function(pattern, target_dir, true_val) {
  csv_list <- list.files(target_dir, pattern = pattern, full.names = T) %>%
    lapply(read.csv)
  mean <- (csv_list %>%
             reduce(`+`)) / length(csv_list)
  bias <- (csv_list %>%
            map(function(df) {return((df - true_val))}) %>%
            reduce(`+`)) / length(csv_list)
  mse <- (csv_list %>%
            map(function(df) {return((df - true_val)^2)}) %>%
            reduce(`+`)) / length(csv_list)
  result <- list(mean, bias, mse)
  names(result) <- c("mean", "bias", "mse")
  return(result)
}

for (i in 1:length(target_dirs)) {
  target_dir <- target_dirs[i]
  pi_val <- true_pi[i]
  effect_val <- true_effect[i]
  pi_stat <- calculate_csv_statistics("pi", target_dir, pi_val)
  effect_stat <- calculate_csv_statistics("effect", target_dir, effect_val)
  print(target_dir)
  print(pi_val)
  print("pi")
  print(pi_stat)
  print("effect")
  print(effect_stat)
}