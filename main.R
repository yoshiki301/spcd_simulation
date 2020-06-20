source("./class/NormalGenerativeModel.R")

seed <- 101
patient_size <- 10000
pi <- 0.3 # placebo response probability
drug_assign <- 0.4

drug_size <- patient_size*drug_assign
responder_size <- patient_size*(1 - drug_assign)*pi
nonresponder_size <- patient_size*(1 - drug_assign)*(1 - pi)

# baseline parameter
mu0 <- 31
sigma0 <- 25

set.seed(seed = seed)

# generate patients data
y01 <- NULL
group <- NULL
stage2_assign <- NULL

for (i in 1:drug_size) {
  y01 <- append(y01, rnorm(1, mean = mu0, sd = sqrt(sigma0)))
  group <- append(group, "drug")
  stage2_assign <- append(stage2_assign, 0)
}

for (i in 1:responder_size) {
  y01 <- append(y01, rnorm(1, mean = mu0, sd = sqrt(sigma0)))
  group <- append(group, "responder")
  if (i <= responder_size*drug_assign) {
    stage2_assign <- append(stage2_assign, 1)
  } else {
    stage2_assign <- append(stage2_assign, 0)
  }
}

for (i in 1:nonresponder_size) {
  y01 <- append(y01, rnorm(1, mean = mu0, sd = sqrt(sigma0)))
  group <- append(group, "nonresponder")
  if (i <= nonresponder_size*drug_assign) {
    stage2_assign <- append(stage2_assign, 1)
  } else {
    stage2_assign <- append(stage2_assign, 0)
  }
}

# model parameter
b1 <- c(1.5, 1.2)
b2 <- c(1.5, 1.2)
b3 <- c(0.5, 1.0)
b4 <- c(0.5, 1.0, 2.5)
b5 <- c(0.5, 1.2)
b6 <- c(0.5, 1.2, 2.5)

GenerativeModel <- NormalGenerativeModel$new(
  b1=b1,
  b2=b2,
  b3=b3,
  b4=b4,
  b5=b5,
  b6=b6
)

res <- mapply(GenerativeModel$generateSample, y01, group, stage2_assign)

patient_data_baseline <- data.frame(
  y01 = y01,
  y1 = res[1,],
  y2 = res[2,],
  group = group,
  stage2_assign = stage2_assign
)

write.csv(patient_data_baseline, "./normal_sample.csv")
