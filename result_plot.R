library(tidyverse)
library(gridExtra)

source("EMy2.0.R")

plot_responder_scatter = function(
    target_dir,
    file_id,
    rho_12,
    patient_size,
    only_placebo = FALSE
) {
    filepath <- paste(target_dir, "data_", file_id, ".csv", sep = "")
    # パラメータ指定のためfilter_を使用する
    data <- read.csv(filepath) %>%
        filter_(paste("rho_12==", rho_12))  %>%
        filter_(paste("patient_size==", patient_size))
    if (only_placebo) {
        data <- data %>%
            filter(g1 == 0)
    }
    p1 <- ggplot(data = data, mapping = aes(x = y01, y = y1, color = response)) +
        geom_point()
    p2 <- ggplot(data = data, mapping = aes(x = y1, y = y2, color = response)) +
        geom_point()
    gridExtra::grid.arrange(p1, p2, nrow = 1)
}

plot_pi_effect = function(
    merged_data
) {
    p1 <- ggplot(data = merged_data, mapping = aes(x = p_m0.5, y = e_m0.5)) +
        geom_point() +
        geom_hline(yintercept=1.5,color="red") +
        geom_vline(xintercept=0.3,color="blue") 
    p2 <- ggplot(data = merged_data, mapping = aes(x = p_0, y = e_0)) +
        geom_point() +
        geom_hline(yintercept=1.5,color="red") +
        geom_vline(xintercept=0.3,color="blue")
    p3 <- ggplot(data = merged_data, mapping = aes(x = p_0.5, y = e_0.5)) +
        geom_point() +
        geom_hline(yintercept=1.5,color="red") +
        geom_vline(xintercept=0.3,color="blue")
    gridExtra::grid.arrange(p1, p2, p3, nrow = 1)
}

target_dir <- "./sim_result/0.1/"

pi_csv <- list.files(target_dir, pattern = "pi", full.names = T) %>%
    lapply(read.csv)
pi_df <- Reduce(rbind, pi_csv) # TODO: Reduceの順番が適当なので直す
pi_df["id"] <- rep(1:length(pi_csv), each = 3)
names(pi_df) <- c("patient_size", "p_m0.5", "p_0", "p_0.5", "file_id")

effect_csv <- list.files(target_dir, pattern = "effect", full.names = T) %>%
    lapply(read.csv)
effect_df <- Reduce(rbind, effect_csv) # TODO: Reduceの順番が適当なので直す
effect_df["id"] <- rep(1:length(effect_csv), each = 3)
names(effect_df) <- c("patient_size", "e_m0.5", "e_0", "e_0.5", "file_id")

merged_df <- inner_join(pi_df, effect_df, by = c("patient_size", "file_id"))

plot_responder_scatter(target_dir, file_id = "190", rho_12 = 0, patient_size = 60)