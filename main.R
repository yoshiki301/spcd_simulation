library("config")

config <<- config::get()
datapath <<- config$datapath
destdir <<- config$destdir
model <<- config$model

# read config file
source(model)
simulation_csv_path <- file.path(destdir, datapath)

seed <- 101
set.seed(seed = seed)

# fixed parameter
pi <- 0.3
delta1 <- 1.5
delta2_nr <- 1.5
h <- 0.5 # delta2_r = h * delta2_nr
drug_assign <- 1/3

# variable parameter
rho_12 <- 0.2 # correlation between y1 and y2
patient_size <- 600

drug_size <- patient_size*drug_assign
responder_size <- patient_size*(1 - drug_assign)*pi
nonresponder_size <- patient_size*(1 - drug_assign)*(1 - pi)

# generate patients data
group <- NULL

for (i in 1:drug_size) {
  group <- append(group, "DD")
}

for (i in 1:responder_size) {
  if (i <= responder_size*drug_assign) {
    group <- append(group, "rPD")
  } else {
    group <- append(group, "rPP")
  }
}

for (i in 1:nonresponder_size) {
  if (i <= responder_size*drug_assign) {
    group <- append(group, "nPD")
  } else {
    group <- append(group, "nPP")
  }
}

GenerativeModel <- NormalGenerativeModel$new(
  pi=pi,
  delta1=delta1,
  delta2_nr=delta2_nr,
  h=h,
  rho_12=rho_12
)

res <- lapply(group, GenerativeModel$generateSample)
y01 <- sapply(res, "[[", 1)
y1 <- sapply(res, "[[", 2)
y2 <- sapply(res, "[[", 3)

patient_data <- data.frame(
  y01 = y01,
  y1 = y1,
  y2 = y2,
  group = group
)

write.csv(patient_data, simulation_csv_path)
