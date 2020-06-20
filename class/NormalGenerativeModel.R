library("R6")
library("assertthat")

NormalGenerativeModel <- R6Class("NormalGenerativeModel",
  public = list(
    mu1 = NULL,
    mu2 = NULL,
    sigma1 = NULL,
    sigma2 = NULL,
    initialize = function(b1, b2, b3, b4, b5, b6) {
      # check assertion and set private member
      assert_that(
        is.numeric(b1),
        is.numeric(b2),
        is.numeric(b3),
        is.numeric(b4),
        is.numeric(b5),
        is.numeric(b6),
        length(b1) == 2,
        length(b2) == 2,
        length(b3) == 2,
        length(b4) == 3,
        length(b5) == 2,
        length(b6) == 3
      )
      private$m.b1 <<- b1
      names(private$m.b1) <<- c("b01", "b11")
      private$m.b2 <<- b2
      names(private$m.b2) <<- c("b02", "b12")
      private$m.b3 <<- b3
      names(private$m.b3) <<- c("b03", "b13")
      private$m.b4 <<- b4
      names(private$m.b4) <<- c("b04", "b14", "b24")
      private$m.b5 <<- b5
      names(private$m.b5) <<- c("b05", "b15")
      private$m.b6 <<- b6
      names(private$m.b6) <<- c("b06", "b16", "b26")
    },
    printPrameters = function() {
      # print private members as parameters
      print(private$m.b1)
      print(private$m.b2)
      print(private$m.b3)
      print(private$m.b4)
      print(private$m.b5)
      print(private$m.b6)
    },
    generateSample = function(baseline, group, stage2_assign = NULL) {
      # generate samples following different normal distribution in 2 stages
      # return vector of stage1 outcome y1 and stage2 outcome y2
      assert_that(
        is.numeric(baseline),
        length(baseline) == 1,
        (group == "drug" || group == "responder" || group == "nonresponder")
      )
      
      # stage 1
      if (group == "drug") {
        self$mu1 <- private$m.meanFunction(private$m.b1, baseline)
        self$sigma1 <- 49
      } else if (group == "responder") {
        self$mu1 <- private$m.meanFunction(private$m.b3, baseline)
        self$sigma1 <- 4
      } else {
        self$mu1 <- private$m.meanFunction(private$m.b5, baseline)
        self$sigma1 <- 4
      }
      y1 <- rnorm(1, self$mu1, sqrt(self$sigma1))
      
      # stage 2
      if (group == "drug") {
        self$mu2 <- private$m.meanFunction(private$m.b2, y1)
        self$sigma2 <- 49
      } else if (group == "responder") {
        assert_that(stage2_assign == 0 || stage2_assign == 1)
        self$mu2 <- private$m.meanFunction(private$m.b4, y1, assign = stage2_assign)
        self$sigma2 <- 4
      } else {
        assert_that(stage2_assign == 0 || stage2_assign == 1)
        self$mu2 <- private$m.meanFunction(private$m.b6, y1, assign = stage2_assign)
        self$sigma2 <- 4
      }
      y2 <- rnorm(1, self$mu2, sqrt(self$sigma2))
      
      return (c(y1, y2))
    }
  ),
  private = list(
    # Drug parameters
    m.b1 = NULL,
    m.b2 = NULL,
    # Placebo responders parameters
    m.b3 = NULL,
    m.b4 = NULL,
    # Placebo nonresponders parameters
    m.b5 = NULL,
    m.b6 = NULL,
    m.meanFunction = function(b, y, assign = NULL) {
      # calculate mean by regression model
      # the length(b) is 3 if the stage2 assign exists, otherwise 2.
      if (length(b) == 2) {
        return (b[1] + b[2]*y)
      } else {
        return (b[1] + b[2]*y + b[3]*assign)
      }
    }
  )
)
