##--------------------------------------------------------------------##
##                      Global optimization method                    ##
##                              based on                              ##
##     the selective averaging coordinate  with restrictions          ##
##--------------------------------------------------------------------##

#' Search for extremum value using global optimization method based on the selective averaging coordinate  with restrictions.
#'
#' @param x starting point coordinate
#' @param delta increment x denotes search range
#' @param fitness function to search for extremum
#' @param lower lower extremum limits
#' @param upper upper extremum limits
#' @param n amount of test points
#' @param e precision constant
#' @param M maximum number of iterations
#' @param y stretch factor
#' @param q matched fixed parameter
#' @param kernelType type of kernel function
#' @param r matched fixed parameter
#' @param s core selectivity factor
#' @return potential point extremum
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
#' x<-sacExtended(x=c(10,10),
#' delta=c(20,20),
#' lower=c(-10,-10),
#' upper = c(30,30),
#' fitness = f)
#' summary(x)
sacExtended<- function(type = c("sacNormal", "sacExtended", "sacIterative"),
                       x,
                       delta,
                       fitness,
                       lower,
                       upper,
                       n = sacControl(type)$n,
                       e = sacControl(type)$e,
                       M = sacControl(type)$M,
                       y = sacControl(type)$y,
                       q = sacControl(type)$q,
                       kernelType = sacControl(type)$kernelType,
                       r = sacControl(type)$r,
                       s = sacControl(type)$s) {
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

  k <- 1
  testX <- matrix(0, n, length(x))
  fValues <- rep(0, n)
  uValues <- matrix(0, n, length(x))
  p <- rep(0, n)
  pNorm <- rep(0, n)
  allX <- matrix(0, M, length(x))
  allX[1, ] <- x
  allResults <- fitness(x)

  cols <- createColForDF(length(x))
  resultObj  <- extResult(iterations = 0,
                        allX = data.frame(Iteration = 0,  cols, stringsAsFactors=FALSE),
                        allDelta = data.frame(Iteration = 0,  cols, stringsAsFactors=FALSE),
                        uValues = data.frame(Iteration = 0, "number of uValue",  cols, stringsAsFactors=FALSE),

                        testPoint = data.frame(Iteration = 0, "number of testPoint", cols, "fintess value", stringsAsFactors=FALSE),

                        gminValues = data.frame(Iteration = 0,  "number of gmin", "value", stringsAsFactors=FALSE),
                        pValues = data.frame(Iteration = 0,  "number of p", "value", stringsAsFactors=FALSE),
                        pNormValues = data.frame(Iteration = 0, "number of pNorm", "value", stringsAsFactors=FALSE),

                        x = x,
                        delta = delta,
                        lower = lower,
                        upper = upper,
                        n = n,
                        e = e,
                        M = M,
                        y = y,
                        q = q,
                        kernelType = kernelType,
                        r = r,
                        s = s,
                        func = fitness
  )
  if (kernelType == "kernelExponential" ||
      kernelType == "kernelHyperbolic" ||
      kernelType == "kernelToDegreeS" ||
      kernelType == "kernelExpHyperbolic")
    kernelFunction <- get(kernelType)
  else
    kernelFunction <- get("kernelExponential")
  while (k < M) {
    for (i in 1:n) {
      for (j in 1:length(x)) {
        if (x[j] - delta[j] < lower[j])
          a <- (((lower[j] - x[j]) / delta[j]) + 1) / 2
        else
          a <- 0
        if (x[j] + delta[j] > upper[j])
          b <- (((upper[j] - x[j]) / delta[j]) + 1) / 2
        else
          b <- 1
        uValues[i, j] <- stats::runif(1, a, b) * 2 - 1
      }
      resultObj @uValues = rbind(resultObj @uValues, addRowForDF(k, c(i, uValues[i,])))
      testX[i, ] <- x + delta * uValues[i, ]
      fValues[i] <- fitness(testX[i, ])
      resultObj @testPoint = rbind(resultObj @testPoint, addRowForDF(k, c(i, testX[i,], fValues[i])))
    }
    gmin <- rep(0, n)
    for (i in 1:n) {
      a <- fValues[i] - min(fValues)
      b <- max(fValues) - min(fValues) + 0.0000000000000001
      gmin[i] <- a / b
      resultObj @gminValues = rbind(resultObj @gminValues, addRowForDF(k, c(i, gmin[i])))
    }
    for (i in 1:n) {
      p[i] <- kernelFunction(z = gmin[i], r = r, s = s)
      resultObj @pValues = rbind(resultObj @pValues, addRowForDF(k, c(i, p[i])))
    }
    for (i in 1:n) {
      pNorm[i] <- p[i] / sum(p)
      resultObj @pNormValues = rbind(resultObj @pNormValues, addRowForDF(k, c(i, pNorm[i])))
    }
    for (i in 1:length(x))
      x[i] <- x[i] + delta[i] * sum(sapply(1:n, function(x) {
        uValues[x, i] * pNorm[x]
      }))
    for (i in 1:length(x))
      delta[i] <- y * delta[i] * ( (sum(sapply(1:n, function(x) {
        (abs(uValues[x, i]) ^ (q)) * pNorm[x]
      }
      ))) ^ (1 / q))
    resultObj @allX = rbind(resultObj @allX, addRowForDF(k, x))
    resultObj @allDelta = rbind(resultObj @allDelta, addRowForDF(k, delta))
    k <- k + 1
    allX[k, ] <- x
    allResults <- c(allResults, fitness(x))
    flag <- (max(delta) <= e)
    if (flag == TRUE)
      break
  }
  resultObj @iterations <- k
  resultObj @extremePoint <- allX[which.min(allResults), ]
  resultObj @fitnessValue <- allResults[which.min(allResults)]
  return(resultObj )
}

extResult <- methods::setClass("extResult", slots = c(iterations = "numeric",
                                                  x = "numeric",
                                                  delta = "numeric",
                                                  lower = "numeric",
                                                  upper = "numeric",
                                                  n = "numeric",
                                                  e = "numeric",
                                                  M = "numeric",
                                                  y = "numeric",
                                                  q = "numeric",
                                                  kernelType = "character",
                                                  r = "numeric",
                                                  s = "numeric",
                                                  allX = "data.frame",
                                                  allDelta = "data.frame",
                                                  fitnessValue = "numeric",
                                                  extremePoint = "numeric",
                                                  uValues = "data.frame",
                                                  testPoint = "data.frame",
                                                  gminValues = "data.frame",
                                                  pValues = "data.frame",
                                                  pNormValues = "data.frame",
                                                  func = 'function'
),
package = "SAC")

setMethod("summary", "extResult",
          function(object)
          {
            cat(cli::rule(left = crayon::bold("SAC Algorithm"),
                          width = min(getOption("width"),40)), "\n\n")
            cat("+-----------------------------------+\n")
            cat("|               sacExtended         |\n")
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
            cat("Selected kernel function     = ", object@kernelType, "\n")
            cat("r                            = ", object@r, "\n")
            cat("s                            = ", object@s, "\n")
            cat("+-----------------------------------+\n\n")
            cat(cli::rule(left = crayon::bold("Progress"),
                          width = min(getOption("width"),40)), "\n")
            for (i in 2:nrow(object@allX)) {
              cat(cli::rule(left = crayon::bold("Iteration", (i-1)),
                            width = min(getOption("width"),40)), ":\nFound point: ", as.numeric(object@allX[i,2:length(object@allX)]), "\nNew increment:", as.numeric(object@allDelta[i,2:length(object@allX)]),  "\n")
            }
            cat("+-----------------------------------+\n")
            cat(cli::rule(left = crayon::bold("Results"),
                          width = min(getOption("width"),40)), "\n")
            cat("Iterations                   =", object@iterations,  "\n")
            cat("Fitness function value       =", object@fitnessValue,  "\n")
            cat("Extreme Point                =", object@extremePoint,  "\n")
          }
)

setMethod("print", "extResult", function(x, ...) utils::str(x))

setMethod("show", "extResult",
          function(object)
          { cat("An object of class \"extResult\"\n")
            cat("Available slots:\n")
            print(slotNames(object))
          })
