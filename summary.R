library("config")

config <<- config::get()
datapath <<- config$datapath
destdir <<- config$destdir

simulation_csv_path <- file.path(destdir, datapath)

patient_data <- read.csv(simulation_csv_path)

y01 <- patient_data$y01
y1 <- patient_data$y1
y2 <- patient_data$y2

y02 <- y01 + y1
y03 <- y02 + y2

plot(y01, y1)
plot(y1, y2)

plot(y01, y02)
plot(y01, y03)
plot(y02, y03)