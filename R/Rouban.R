#' Rouban: package allows you to use the algorithm for finding global extrema by averaging coordinates
#'
#' @name Rouban
#' @docType package
#'
NULL

nuclearFunction <- function(z, x = 1, r = 2, s = 3) {
  if (all.equal(z, as.double(z), check.attributes = FALSE) != TRUE
      || is.nan(z)
      || length(z) != 1)
    stop("Incorrect value of z. z must be double value")
  if (all.equal(s, as.double(s), check.attributes = FALSE) != TRUE
      || s < 0
      || is.nan(s)
      || length(s) != 1)
    stop("Error, expected s is double and s >= 0")
  if (all.equal(r, as.integer(r), check.attributes = FALSE) != TRUE
      || r < 1
      || is.nan(r)
      || length(r) != 1)
    stop("Error, expected r is integer and r > 0")
  if (all.equal(x, as.integer(x), check.attributes = FALSE) != TRUE
      || x < 1
      || x > 4
      || is.nan(x)
      || length(x) != 1) {
    warning("you choose incorrect number of nuclear function.
            Choose default value: 1")
    x <- 1
  }
  if (x == 1)
    p <- (1 - (z ^ r)) ^ s
  else if (x == 2)
    p <- exp(-s * (z ^ r))
  else if (x == 3)
    p <- 1 / (s * (z ^ r) + 1)
  else if (x == 4)
    p <- 1 / (z ^ s)
  return(p)
}

#' Search for extremum value using coordinate averaging method
#'
#' @param x starting point coordinate
#' @param delta increment x denotes search range
#' @param f search function
#' @param lower lower extremum limits
#' @param upper upper extremum limits
#' @param n amount of points
#' @param e precision constant
#' @param M maximum number of iterations
#' @param y parameter of Rouban's algorithm
#' @param q parameter of Rouban's algorithm
#' @param nuclearType type of nuclear function
#' @param r parameter for nuclear function
#' @param s parameter for nuclear function
#' @return numeric the estimated amount of good honey
#' @export
#' @examples
#' f<-function(x) {
#' z<-7*(abs(x[1])^2) + 7*(abs(x[2])^2)
#' z<-c(z, 5*(abs(x[1]-3))^0.8 + 5*(abs(x[2]-3)^0.6) + 6)
#' z<-c(z, 5*(abs(x[1]-6))^1.3 + 5*(abs(x[2]-6)^1.3) + 2)
#' z<-c(z, 5*(abs(x[1]-6))^1 + 5*(abs(x[2]+6)^1) + 8)
#' z<-c(z, 4*(abs(x[1]+6))^1.5 + 4*(abs(x[2]+6)^1.5) + 7)
#' z<-c(z, 5*(abs(x[1]+3))^1.8 + 5*(abs(x[2])^1.8) + 9)
#' z<-c(z, 6*(abs(x[1]+6))^0.6 + 6*(abs(x[2]-6)^0.9))
#' return(min(z))
#' }
#' x<-Rouban(x=c(10,10),
#' delta=c(20,20),
#' lower=c(-10,-10),
#' upper = c(30,30),
#' f = f,
#' n=500,
#' y=1,
#' q=2,
#' s=100,
#' e=0.0001,
#' r=2,
#' nuclearType = 1)
#' print(x)
Rouban <- function(x,
                  delta,
                  f,
                  lower,
                  upper,
                  n = 500,
                  e = 0.001,
                  M = 1000,
                  y = 1,
                  q = 2,
                  nuclearType = 1,
                  r = 2,
                  s = 100) {
  if (all.equal(x, as.double(x), check.attributes = FALSE) != TRUE
      || testOnNaN(x))
    stop("Incorrect value of x. x must be numeric vector")
  if (all.equal(delta, as.double(delta), check.attributes = FALSE) != TRUE
      || testOnNaN(delta))
    stop("Incorrect value of delta. delta must be numeric vector")
  if (all.equal(upper, as.double(upper), check.attributes = FALSE) != TRUE
      || testOnNaN(upper))
    stop("Incorrect value of upper. upper must be numeric vector")
  if (all.equal(lower, as.double(lower), check.attributes = FALSE) != TRUE
      || testOnNaN(lower))
    stop("Incorrect value of lower. lower must be numeric vector")
  if (length(x) != length(delta)
      || length(x) != length(lower)
      || length(x) != length(upper)) {
    stop("Error. x, delta, upper and lower must have one size")
  }
  if (all.equal(n, as.integer(n), check.attributes = FALSE) != TRUE
      || n < 1
      || is.nan(n)
      || length(n) != 1)
    stop("Error, expected n is integer and n > 0")
  if (all.equal(y, as.double(y), check.attributes = FALSE) != TRUE
      || y < 0
      || is.nan(y)
      || length(y) != 1)
    stop("Error, expected y is double and y > 0")
  if (all.equal(q, as.integer(q), check.attributes = FALSE) != TRUE
      || q < 1
      || is.nan(q)
      || length(q) != 1)
    stop("Error, expected q is integer and q >= 1")
  if (all.equal(e, as.double(e), check.attributes = FALSE) != TRUE
      || e <= 0
      || is.nan(e)
      || length(e) != 1)
    stop("Error, expected e is double and e > 0")
  if (testOnUncorrentUpperAndBound(lower, upper))
    stop("Error, your value of lower vector > than value of upper vector")
  k <- 1
  testX <- matrix(0, n, length(x))
  fValues <- rep(0, n)
  uValues <- matrix(0, n, length(x))
  p <- rep(0, n)
  pNorm <- rep(0, n)
  allX <- matrix(0, M, length(x))
  allX[1, ] <- x
  allResults <- f(x)

  cols <- createColForDF(length(x))
  resultObj  <- sResult(iterations = 0,
                         allX = data.frame(Iteration = 0,  cols, stringsAsFactors=FALSE),
                         allDelta = data.frame(Iteration = 0,  cols, stringsAsFactors=FALSE),
                         x = x,
                         delta = delta,
                         lower = lower,
                         upper = upper,
                         n = n,
                         e = e,
                         M = M,
                         y = y,
                         q = q,
                         nuclearType = nuclearType,
                         r = r,
                         s = s
  )

  while (k < M) {
    for (i in 1:n) {
      for (j in 1:length(x))
        uValues[i, j] <- stats::runif(1, 0, 1) * 2 - 1
      testX[i, ] <- x + delta * uValues[i, ]
      fValues[i] <- f(testX[i, ])
    }
    gmin <- rep(0, n)
    for (i in 1:n) {
      a <- fValues[i] - min(fValues)
      b <- max(fValues) - min(fValues)
      gmin[i] <- a / b
    }
    for (i in 1:n)
      p[i] <- nuclearFunction(x = nuclearType, z = gmin[i], r = r, s = s)
    for (i in 1:n)
      pNorm[i] <- p[i] / sum(p)
    for (i in 1:length(x))
      x[i] <- x[i] + delta[i] * sum(sapply(1:n, function(x) {
        uValues[x, i] * pNorm[x]
        }))
    for (i in 1:length(x)) {
      if (x[i] < lower[i])
        x[i] <- lower[i]
      if (x[i] < lower[i])
        x[i] <- lower[i]
      if (x[i] > upper[i])
        x[i] <- upper[i]
      if (x[i] > upper[i])
        x[i] <- upper[i]
    }
    for (i in 1:length(x))
      delta[i] <- y * delta[i] * ( (sum(sapply(1:n, function(x) {
        (abs(uValues[x, i]) ^ (q)) * pNorm[x]
        }
        ))) ^ (1 / q))
    resultObj @allX = rbind(resultObj @allX, addRowForDF(k, x))
    resultObj @allDelta = rbind(resultObj @allDelta, addRowForDF(k, delta))
    k <- k + 1
    allX[k, ] <- x
    allResults <- c(allResults, f(x))
    flag <- (max(delta) <= e)
    if (flag == TRUE)
      break
  }
  resultObj @iterations = k
  resultObj @results = allX[which.min(allResults), ]
  return(resultObj )
}

createColForDF <- function(ncols) {
  cols <- vector("list", ncols)
  for (i in 1:ncols) {
    value <- paste0("x",i)
    cols[i] <- value
  }
  return(cols)
}

addRowForDF <- function(k, x) {
  new_col <- vector("list", length(x) + 1)
  new_col[1] <- k
  for (i in 2:(length(x) + 1)) {
    new_col[i] <- x[i-1]
  }
  return(new_col)
}

sResult <- methods::setClass("sResult", slots = c(iterations = "numeric",
                                                              x = "numeric",
                                                              delta = "numeric",
                                                              lower = "numeric",
                                                              upper = "numeric",
                                                              n = "numeric",
                                                              e = "numeric",
                                                              M = "numeric",
                                                              y = "numeric",
                                                              q = "numeric",
                                                              nuclearType = "numeric",
                                                              r = "numeric",
                                                              s = "numeric",
                                                              allX = "data.frame",
                                                              allDelta = "data.frame",
                                                              results = "numeric"
                                                              ),
                                   package = "Rouban")

setMethod("summary", "sResult",
          function(object)
          {
            cat(cli::rule(left = crayon::bold("Rouban Algorithm"),
                          width = min(getOption("width"),40)), "\n\n")
            cat("+-----------------------------------+\n")
            cat("|               Rouban              |\n")
            cat("+-----------------------------------+\n\n")
            cat(cli::rule(left = crayon::bold("Algorithm settings"),
                          width = min(getOption("width"),40)), "\n")
            cat("Starting point               = ", object@x, "\n")
            cat("Increment                    = ", object@delta, "\n")
            cat("Lower bound                  = ", object@lower, "\n")
            cat("Upper bound                  = ", object@upper, "\n")
            cat("Accuracy                     = ", object@e, "\n")
            cat("Max.iterations               = ", object@M, "\n")
            cat("y                            = ", object@y, "\n")
            cat("q                            = ", object@q, "\n")
            cat("Selected nuclear function    = ", object@nuclearType, "\n")
            cat("r                            = ", object@r, "\n")
            cat("s                            = ", object@s, "\n")
            cat("+-----------------------------------+\n\n")
            cat(cli::rule(left = crayon::bold("Results"),
                          width = min(getOption("width"),40)), "\n")
            cat("Iterations                   =", object@iterations,  "\n")
            cat("Fitness function value       =", object@results,  "\n")
          }
)

