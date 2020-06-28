library("R6")
library("assertthat")

NormalGenerativeModel <- R6Class("NormalGenerativeModel",
  public = list(
    initialize = function(pi, delta1, delta2_nr, h, rho_12) {
      # check assertion and set private member
      assert_that(
        is.numeric(pi),
        is.numeric(delta1),
        is.numeric(delta2_nr),
        is.numeric(h),
        is.numeric(rho_12),
        length(pi) == 1,
        length(delta1) == 1,
        length(delta2_nr) == 1,
        length(h) == 1,
        length(rho_12) == 1
      )
      private$m.pi <- pi
      private$m.delta1 <- delta1
      private$m.delta2_nr <- delta2_nr
      private$m.h <- h
      private$m.rho_12 <- rho_12
      
      # execute cholesky factorization
      corr <- rbind(c(1, 0.1, 0.1),
                    c(0.1, 1, private$m.rho_12),
                    c(0.1,private$m.rho_12,1))
      private$m.L <- chol(corr)
    },
    printPrameters = function() {
      # print private members as parameters
      print(private$m.pi)
      print(private$m.delta1)
      print(private$m.delta2_nr)
      print(private$m.h)
      print(private$m.rho_12)
      print(private$m.L)
    },
    generateSample = function(group) {
      # generate samples following different normal distribution in 2 stages
      # return vector of stage1 outcome y1 and stage2 outcome y2
      assert_that(
        (group == "DD" || group == "rPD" || group == "rPP" || group == "nPD" || group == "nPP")
      )
      
      # generate 3 standard normal random number
      x <- rnorm(3, mean=0, sd=1)
      
      # baseline
      z1 <- private$m.L[1,1] * x[1]
      y01 <- 5*z1 + 31
      
      # change in stage1
      z2 <- private$m.L[1,2] * x[1] + private$m.L[2,2] * x[2]
      switch (group,
        "DD" = y1 <- 7*z2 + (16*private$m.pi + 9.7*(1-private$m.pi) + private$m.delta1),
        "rPD" = y1 <- 2*z2 + 16,
        "rPP" = y1 <- 2*z2 + 16,
        "nPD" = y1 <- 2*z2 + 9.7,
        "nPP" = y1 <- 2*z2 + 9.7,
      )
      
      # change in stage2
      z3 <- private$m.L[1,3] * x[1] + private$m.L[2,3] * x[2] + private$m.L[3,3] * x[3]
      switch (group,
              "DD" = y2 <- 7*z3 + 0.6*y1,
              "rPD" = y2 <- 2*z3 + (0.6*y1 + private$m.delta2_nr*private$m.h),
              "rPP" = y2 <- 2*z3 + 0.6*y1,
              "nPD" = y2 <- 2*z3 + (0.6*y1 + private$m.delta2_nr),
              "nPP" = y2 <- 2*z3 + 0.6*y1,
      )
      
      return (c(y01, y1, y2))
    }
  ),
  private = list(
    m.pi = NULL,
    m.delta1 = NULL,
    m.delta2_nr = NULL,
    m.h = NULL,
    m.rho_12 = NULL,
    m.L = NULL
  )
)
