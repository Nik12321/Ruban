#' Ruban: package allows you to use the algorithm for finding global extrema by averaging coordinates
#'
#' @name Ruban
#' @docType package
#'
NULL

#' Use of nuclear function
#'
#' The choice is made at the discretion of the user
#'
#' @param z value
#' @param x choice of nuclear function
#' @param r nuclear function parameter
#' @param s nuclear function parameter
#' @return result of nuclear function
#' @export
#' @examples
#' z <- 5
#' x <- nuclearFunction(z)
#' print(z)
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
#' @param y parameter of Ruban's algorithm
#' @param q parameter of Ruban's algorithm
#' @param nuclear_function choice of nuclear function
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
#' x<-Ruban(x=c(10,10),
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
#' nuclear_function = 1)
#' print(x)
Ruban <- function(x,
                  delta,
                  f,
                  lower,
                  upper,
                  n = 500,
                  e = 0.001,
                  M = 1000,
                  y = 1,
                  q = 2,
                  nuclear_function = 1,
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
  test_x <- matrix(0, n, length(x))
  f_values <- rep(0, n)
  u_values <- matrix(0, n, length(x))
  p <- rep(0, n)
  p_norm <- rep(0, n)
  all_x <- matrix(0, M, length(x))
  all_x[1, ] <- x
  all_results <- f(x)
  while (k < M) {
    for (i in 1:n) {
      for (j in 1:length(x))
        u_values[i, j] <- stats::runif(1, 0, 1) * 2 - 1
      test_x[i, ] <- x + delta * u_values[i, ]
      f_values[i] <- f(test_x[i, ])
    }
    gmin <- rep(0, n)
    for (i in 1:n) {
      a <- f_values[i] - min(f_values)
      b <- max(f_values) - min(f_values)
      gmin[i] <- a / b
    }
    for (i in 1:n)
      p[i] <- nuclearFunction(x = nuclear_function, z = gmin[i], r = r, s = s)
    for (i in 1:n)
      p_norm[i] <- p[i] / sum(p)
    for (i in 1:length(x))
      x[i] <- x[i] + delta[i] * sum(sapply(1:n, function(x) {
        u_values[x, i] * p_norm[x]
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
        (abs(u_values[x, i]) ^ (q)) * p_norm[x]
        }
        ))) ^ (1 / q))
    k <- k + 1
    all_x[k, ] <- x
    all_results <- c(all_results, f(x))
    flag <- (max(delta) <= e)
    if (flag == TRUE)
      break
  }
  return(all_x[which.min(all_results), ])
}

testOnNaN <- function(value) {
  value <- sapply(value, function(x) {
    !is.nan(x)
    })
  for (i in 1:length(value))
    if (value[i] != TRUE) {
      return(TRUE)
    }
  return(FALSE)
}

testOnUncorrentUpperAndBound <- function(lower, upper) {
  value <- lower > upper
  for (i in 1:length(value))
    if (value[i] == TRUE) {
      return(TRUE)
    }
  return(FALSE)
}
