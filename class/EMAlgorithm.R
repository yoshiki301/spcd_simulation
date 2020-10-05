library("R6")

EMAlgorithmModel <- R6Class("EMAlgorithmModel",
  public = list(
    initialize = function(df, eps, w) {
      private$m.pi_i <- rep_len(private$m.pi, nrow(df))
      private$m.eps <- eps
      private$m.w <- w
    },
    printPramas = function() {
      print(private$m.pi)
    },
    printProcess = function() {
      print(private$m.pi_process)
      plot(private$m.pi_process, xlab="step", ylab="π", xlim=c(0,7), type="l")
    },
    printQlist = function() {
      print(private$m.Q_list)
      plot(private$m.Q_list, type="l")
    },
    calculateQfunciton = function(df) {
      Q <- 0
      for (i in 1:nrow(df)) {
        row <- df[i,]
        group_label <- as.character(row["group"])
        y01 <- as.numeric(row["y01"])
        y1 <- as.numeric(row["y1"])
        if (group_label != "DD") {
          mu101 <- private$m.b03 + private$m.b13 * y01
          mu102 <- private$m.b05 + private$m.b15 * y01
          q <- private$m.pi * log(dnorm(y1, mean = mu101, sd = sqrt(private$m.sigma01))) + (1 - private$m.pi) * log(dnorm(y1, mean = mu102, sd = sqrt(private$m.sigma02)))
          Q <- Q + q
        }
      }
      Q <- Q / nrow(df)
      return(Q)
    },
    eStep = function(df) {
      pi_list <- c()
      # update pi_i except DD group
      for (i in 1:nrow(df)) {
        row <- df[i,]
        group_label <- as.character(row["group"])
        if (group_label != "DD") {
          y01 <- as.numeric(row["y01"])
          y1 <- as.numeric(row["y1"])
          mu101 <- private$m.b03 + private$m.b13 * y01
          mu102 <- private$m.b05 + private$m.b15 * y01
          p101 <- dnorm(y1, mean = mu101, sd = sqrt(private$m.sigma01))
          p102 <- dnorm(y1, mean = mu102, sd = sqrt(private$m.sigma02))
          private$m.pi_i[i] <- (private$m.pi * p101) / ((private$m.pi * p101) + ((1 - private$m.pi) * p102))
          pi_list <- append(pi_list, private$m.pi_i[i])
        }
      }
      pi <- mean(pi_list)
      private$m.pi <- pi
      private$m.pi_process <- append(private$m.pi_process, pi)
    },
    mStep = function(df) {
      b01_num <- c()
      b01_den <- c()
      b11_num <- c()
      b11_den <- c()
      sigma11_num <- c()
      sigma11_den <- c()
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
      b06_num <- c()
      b06_den <- c()
      b16_num <- c()
      b16_den <- c()
      b26_num <- c()
      b26_den <- c()
      sigma202_num <- c()
      sigma202_den <- c()
      for (i in 1:nrow(df)) {
        row <- df[i,]
        group_label <- as.character(row["group"])
        y01 <- as.numeric(row["y01"])
        y1 <- as.numeric(row["y1"])
        y2 <- as.numeric(row["y2"])
        if (group_label != "DD") {
          b03_num <- append(b03_num, private$m.pi_i[i] * (y1 - private$m.b13 * y01))
          b03_den <- append(b03_den, private$m.pi_i[i])
          b13_num <- append(b13_num, private$m.pi_i[i] * y01 * (y1 - private$m.b03))
          b13_den <- append(b13_den, private$m.pi_i[i] * (y01**2))
          sigma01_num <- append(sigma01_num, private$m.pi_i[i] * (y1 - private$m.b03 - private$m.b13 * y01)**2)
          sigma01_den <- append(sigma01_den, private$m.pi_i[i])
          b05_num <- append(b05_num, (1 - private$m.pi_i[i]) * (y1 - private$m.b15 * y01))
          b05_den <- append(b05_den, (1 - private$m.pi_i[i]))
          b15_num <- append(b15_num, (1 - private$m.pi_i[i]) * y01 * (y1 - private$m.b05))
          b15_den <- append(b15_den, (1 - private$m.pi_i[i]) * (y01**2))
          sigma02_num <- append(sigma02_num, (1 - private$m.pi_i[i]) * (y1 - private$m.b05 - private$m.b15 * y01)**2)
          sigma02_den <- append(sigma02_den, (1 - private$m.pi_i[i]))
          if (group_label == "rPD" || group_label == "nPD") {
            b06_num <- append(b06_num, (1 - private$m.pi_i[i]) * (y2 - private$m.b16 * y1 - private$m.b26))
            b06_den <- append(b06_den, (1 - private$m.pi_i[i]))
            b16_num <- append(b16_num, (1 - private$m.pi_i[i]) * y1 * (y2 - private$m.b06 - private$m.b26))
            b16_den <- append(b16_den, (1 - private$m.pi_i[i]) * (y1**2))
            b26_num <- append(b26_num, (1 - private$m.pi_i[i]) * (y2 - private$m.b06 - private$m.b16 * y1))
            b26_den <- append(b26_den, (1 - private$m.pi_i[i]))
            sigma202_num <- append(sigma202_num, (1 - private$m.pi_i[i]) * (y2 - private$m.b06 - private$m.b16 * y1 - private$m.b26)**2)
            sigma202_den <- append(sigma202_den, (1 - private$m.pi_i[i]))
          } else {
            b06_num <- append(b06_num, (1 - private$m.pi_i[i]) * (y2 - private$m.b16 * y1))
            b06_den <- append(b06_den, (1 - private$m.pi_i[i]))
            b16_num <- append(b16_num, (1 - private$m.pi_i[i]) * y1 * (y2 - private$m.b06))
            b16_den <- append(b16_den, (1 - private$m.pi_i[i]) * (y1**2))
            sigma202_num <- append(sigma202_num, (1 - private$m.pi_i[i]) * (y2 - private$m.b06 - private$m.b16 * y1)**2)
            sigma202_den <- append(sigma202_den, (1 - private$m.pi_i[i]))
          }
        } else {
          b01_num <- append(b01_num, (y1 - private$m.b11 * y01))
          b01_den <- append(b01_den, 1)
          b11_num <- append(b11_num, y01 * (y1 - private$m.b01))
          b11_den <- append(b11_den, y01**2)
          sigma11_num <- append(sigma11_num, (y1 - private$m.b01 - private$m.b11 * y01)**2)
          sigma11_den <- append(sigma11_den, 1)
        }
      }
      
      private$m.b01 <- sum(b01_num) / sum(b01_den)
      private$m.b11 <- sum(b11_num) / sum(b11_den)
      private$m.sigma11 <- sum(sigma11_num) / sum(sigma11_den)
      private$m.b03 <- sum(b03_num) / sum(b03_den)
      private$m.b13 <- sum(b13_num) / sum(b13_den)
      private$m.sigma01 <- sum(sigma01_num) / sum(sigma01_den)
      private$m.b05 <- sum(b05_num) / sum(b05_den)
      private$m.b15 <- sum(b15_num) / sum(b15_den)
      private$m.sigma02 <- sum(sigma02_num) / sum(sigma02_den)
      private$m.b06 <- sum(b06_num) / sum(b06_den)
      private$m.b16 <- sum(b16_num) / sum(b16_den)
      private$m.b26 <- sum(b26_num) / sum(b26_den)
      private$m.sigma202 <- sum(sigma202_num) / sum(sigma202_den)
    },
    emSimulation = function(df, limit_step = 1000) {
      if (length(private$m.pi_process) == 0) {
        pi_s <- private$m.pi
        private$m.pi_process <- append(private$m.pi_process, pi_s)
        Q <- self$calculateQfunciton(df)
        private$m.Q_list <- append(private$m.Q_list, Q)
        self$eStep(df)
        pi_ss <- rev(private$m.pi_process)[1]
        s <- 1
        # print(s)
        # print(private$m.pi_i[205])
        # print(private$m.pi_i[325])
        show_pi <- c()
        show_sigma <- c()
        while (abs(pi_s-pi_ss) > private$m.eps && s <= limit_step) {
          
          # print("---params---")
          # print(private$m.pi)
          # print(private$m.sigma01)
          # print(private$m.sigma02)
          # print("---")
          show_pi <- append(show_pi, private$m.pi)
          show_sigma <- append(show_sigma, private$m.sigma01)
          
          self$mStep(df)
          Q <- self$calculateQfunciton(df)
          private$m.Q_list <- append(private$m.Q_list, Q)
          pi_s <- pi_ss
          self$eStep(df)
          pi_ss <- rev(private$m.pi_process)[1]
          s <- s + 1
          # print(s)
          # print(private$m.pi_i[205]) # responder
          # print(private$m.pi_i[325]) # nonresponder
        }
        # print("Estimation has been done.")
        # plot(show_pi, type="l")
        # par(new=T)
        # plot(show_sigma)
      }
      else {
        print("Estimation had been done. Please check output with printProcess().")
      }
    },
    calculateEffect = function(patient_df) {
      mu11 <- private$m.b01 + (private$m.b11 * mean(patient_df[patient_df$group == "DD", "y01"]))
      mu101 <- private$m.b03 + (private$m.b13 * mean(patient_df[patient_df$group != "DD", "y01"]))
      mu102 <- private$m.b05 + (private$m.b15 * mean(patient_df[patient_df$group != "DD", "y01"]))
      mu10 <- private$m.pi * mu101 + (1 - private$m.pi) * mu102
      sigma10 <- (private$m.pi**2 * private$m.sigma01) + ((1 - private$m.pi)**2 * private$m.sigma02)
      delta1 <- mu11 - mu10
      var_delta1 <- private$m.sigma11 - sigma10
      delta2_nr <- private$m.b26
      delta_w <- private$m.w * delta1 + (1 - private$m.w) * delta2_nr
      return(delta_w)
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
    # List of Q function changes
    m.Q_list = c(),
    # Drug Stage I
    m.b01 = 18,
    m.b11 = 0.1,
    m.sigma11 = 49,
    # Drug Stage II
    m.b02 = 0.1,
    m.b12 = 1.0,
    m.sigma21 = 5.0,
    # Placebo responder Stage I
    m.b03 = 18,
    m.b13 = 0.1,
    m.sigma01 = 25,
    # Placebo responder Stage II
    m.b04 = 0.1,
    m.b14 = 1.0,
    m.b24 = 1.0,
    m.sigma201 = 5.0,
    # Placebo nonresponder Stage I
    m.b05 = 5,
    m.b15 = 0.05,
    m.sigma02= 25,
    # Placebo nonresponder Stage II
    m.b06 = 6,
    m.b16 = 0.1,
    m.b26 = 1.4,
    m.sigma202 = 4.0,
    # predefined w
    m.w = 0.5
  ),
)