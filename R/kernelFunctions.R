##--------------------------------------------------------------------##
##                                                                    ##
##                         Kernel Operators                           ##
##                                                                    ##
##--------------------------------------------------------------------##

##
##  Kernel function operators ----
##

kernelExponential <- function(z, r = 2, s = 3) {
  p <- exp(-s * (z ^ r))
  return(p)
}

kernelHyperbolic <- function(z, r = 2, s = 3) {
  p <- 1 / (z ^ s)
  return(p)
}

kernelToDegreeS <- function(z, r = 2, s = 3) {
  p <- (1 - (z ^ r)) ^ s
  return(p)
}

kernelExpHyperbolic <- function(z, r = 2, s = 3) {
  p <- 1 / (s * (z ^ r) + 1)
  return(p)
  }
