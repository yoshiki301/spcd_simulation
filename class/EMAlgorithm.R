library("R6")

EMAlgorithmModel <- R6Class("EMAlgorithmModel",
  public = list(
    initialize = function(df, eps) {
      private$m.pi_i <- rep_len(0.5, nrow(df))
      private$m.eps <- eps
    },
    printPramas = function() {
      print(private$m.pi)
    },
    printProcess = function() {
      print(private$m.pi_process)
    },
    eStep = function(df) {
      pi_list <- c()
      # update pi_i except DD group
      for (i in 1:nrow(df)) {
        row <- df[i,]
        group_label <- as.character(row["group"])
        if (group_label != "DD") {
          y01 <- as.numeric(row["y01"])
          mu101 <- private$m.b03 + private$m.b13 * y01
          mu102 <- private$m.b05 + private$m.b15 * y01
          p101 <- dnorm(y01, mean = mu101, sd = sqrt(private$m.sigma01))
          p102 <- dnorm(y01, mean = mu102, sd = sqrt(private$m.sigma02))
          private$m.pi_i[i] <- (private$m.pi * p101) / ((private$m.pi * p101) + ((1-private$m.pi) * p102))
          pi_list <- append(pi_list, private$m.pi_i[i])
        }
      }
      pi <- mean(pi_list)
      private$m.pi <- pi
      private$m.pi_process <- append(private$m.pi_process, pi)
    },
    mStep = function(df) {
      b03_num <- c()
      b03_den <- c()
      b13_num <- c()
      b13_den <- c()
      sigma01_num <- c()
      sigma01_den <- c()
      b05_num <- c()
      b05_den <- c()
      b15_num <- c()
      b15_den <- c()
      sigma02_num <- c()
      sigma02_den <- c()
      # update only b03, b13, sigma01, b05, b15, sigma02
      for (i in 1:nrow(df)) {
        row <- df[i,]
        group_label <- as.character(row["group"])
        if (group_label != "DD") {
          y01 <- as.numeric(row["y01"])
          y1 <- as.numeric(row["y1"])
          b03_num <- append(b03_num, private$m.pi_i[i] * (y1 - private$m.b13 * y01))
          b03_den <- append(b03_den, private$m.pi_i[i])
          b13_num <- append(b13_num, private$m.pi_i[i] * y01 * (y1 - private$m.b03))
          b13_den <- append(b13_den, private$m.pi_i[i] * (y01^2))
          sigma01_num <- append(sigma01_num, private$m.pi_i[i] * (y1 - private$m.b03 - private$m.b13 * y01)^2)
          sigma01_den <- append(sigma01_den, private$m.pi_i[i])
          b05_num <- append(b05_num, (1 - private$m.pi_i[i]) * (y1 - private$m.b15 * y01))
          b05_den <- append(b05_den, (1 - private$m.pi_i[i]))
          b15_num <- append(b15_num, (1 - private$m.pi_i[i]) * y01 * (y1 - private$m.b05))
          b15_den <- append(b15_den, (1 - private$m.pi_i[i]) * (y01^2))
          sigma02_num <- append(sigma02_num, (1 - private$m.pi_i[i]) * (y1 - private$m.b05 - private$m.b15 * y01)^2)
          sigma02_den <- append(sigma02_den, (1 - private$m.pi_i[i]))
        }
      }
      private$m.b03 <- sum(b03_num) / sum(b03_den)
      private$m.b13 <- sum(b13_num) / sum(b13_den)
      private$m.sigma01 <- sum(sigma01_num) / sum(sigma01_den)
      private$m.b05 <- sum(b05_num) / sum(b05_den)
      private$m.b15 <- sum(b15_num) / sum(b15_den)
      private$m.sigma02 <- sum(sigma02_num) / sum(sigma02_den)
    },
    emSimulation = function(df, limit_step = 10000) {
      if (length(private$m.pi_process) == 0) {
        pi_s <- private$m.pi
        private$m.pi_process <- append(private$m.pi_process, pi_s)
        self$eStep(df)
        pi_ss <- rev(private$m.pi_process)[1]
        s <- 1
        print(s)
        while (abs(pi_s-pi_ss) > private$m.eps && s <= limit_step) {
          self$mStep(df)
          pi_s <- pi_ss
          self$eStep(df)
          pi_ss <- rev(private$m.pi_process)[1]
          s <- s + 1
          print(s)
        }
        print("Simulation has done.")
      }
      else {
        print("Simulation had done. Please check output with printProcess().")
      }
    }
  ),
  private = list(
    # epsilon for convergence
    m.eps = NULL,
    # Latent pi
    m.pi = 0.5,
    # Each pi
    m.pi_i = NULL,
    # Process of iteration
    m.pi_process = c(),
    # Drug Stage I
    m.b01 = 0.1,
    m.b11 = 1.0,
    m.sigma11 = 5.0,
    # Drug Stage II
    m.b02 = 0.1,
    m.b12 = 1.0,
    m.sigma21 = 5.0,
    # Placebo responder Stage I
    m.b03 = 0.1,
    m.b13 = 1.0,
    m.sigma01 = 5.0,
    # Placebo responder Stage II
    m.b04 = 0.1,
    m.b14 = 1.0,
    m.b24 = 1.0,
    m.sigma201 = 5.0,
    # Placebo nonresponder Stage I
    m.b05 = 0.1,
    m.b15 = 1.0,
    m.sigma02= 5.0,
    # Placebo nonresponder Stage II
    m.b06 = 0.1,
    m.b16 = 1.0,
    m.b26 = 1.0,
    m.sigma202 = 5.0
  ),
)