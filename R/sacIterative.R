##--------------------------------------------------------------------##
##                          Function for the                          ##
##                  selection of optimal parameters                   ##
##                        for SAC function                         ##
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
#' x<-sacIterative(x=c(10,10),
#' delta=c(20,20),
#' lower=c(-10,-10),
#' upper = c(30,30),
#' fitness = f)
#' summary(x)
sacIterative <- function(type = c("sacNormal", "sacExtended", "sacIterative"),
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
  resultObj  <- normalResult(iterations = 0,
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
  ITERATION <- readline("Conduct a detailed iteration? (Y/N)")
  while(ITERATION != "N" && ITERATION != "Y")
    ITERATION <- readline()
  if (ITERATION == "Y")
    cat(cli::rule(left = crayon::bold("Detailed process"),
                  width = min(getOption("width"), 40)), "\n")
  while (k < M) {
    pastDelta <- delta
    if (ITERATION == "Y") {
      ANSWER <- readline("To form test points? (Y/N)")
      while(ANSWER != "N" && ANSWER != "Y")
        ANSWER <- readline()
      if (ANSWER == "Y")
        cat(cli::rule(left = crayon::bold("Formation of test points"),
                      width = min(getOption("width"),40)), "\n")
    }
    else
      ANSWER <- "N"
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
          if (ANSWER == "Y")
            cat(paste0("uValues for x", j, " = "), uValues[i, j], "\n")
        }
      testX[i, ] <- x + delta * uValues[i, ]
      fValues[i] <- fitness(testX[i, ])
      if (ANSWER == "Y") {
        cat(paste0("testX (n = ", i, ")  = "), testX[i, ], "\n")
        cat(paste0("fValues        = "),  fValues[i], "\n")
        cat("+-----------------------------------+\n\n")
        ANSWER <- readline("")
        if (ANSWER != "N")
          ANSWER = "Y"
      }
    }
    if (ANSWER == "Y")
      ANSWER <- readline("Test points were formed. Press ENTER to continue.")

    gmin <- rep(0, n)
    if (ITERATION == "Y") {
      ANSWER <- readline("To form gmin? (Y/N)")
      while(ANSWER != "N" && ANSWER != "Y")
        ANSWER <- readline()
    }
    if (ANSWER == "Y")
      cat(cli::rule(left = crayon::bold("Formation of gmin"),
                    width = min(getOption("width"),40)), "\n")
    for (i in 1:n) {
      a <- fValues[i] - min(fValues)
      b <- max(fValues) - min(fValues) + 0.0000000000000001
      gmin[i] <- a / b
      if (ANSWER == "Y") {
        cat(paste0("gmin (n = ", i, ") = "), gmin[i], "\n")
        cat("+-----------------------------------+\n\n")
        ANSWER <- readline("")
        if (ANSWER != "N")
          ANSWER = "Y"
      }
    }
    if (ANSWER == "Y")
      ANSWER <- readline("All gmin were formed. Press ENTER to continue.")

    if (ITERATION == "Y") {
      ANSWER <- readline("To form value of the kernel function? (Y/N)")
      while(ANSWER != "N" && ANSWER != "Y")
        ANSWER <- readline()
    }
    if (ANSWER == "Y")
      cat(cli::rule(left = crayon::bold("Formation of value of the kernel function"),
                    width = min(getOption("width"),40)), "\n")
    for (i in 1:n) {
      p[i] <- kernelFunction(z = gmin[i], r = r, s = s)
      if (ANSWER == "Y") {
        cat(paste0("p (n = ", i, ") = "), p[i], "\n")
        cat("+-----------------------------------+\n\n")
        ANSWER <- readline("")
        if (ANSWER != "N")
          ANSWER = "Y"
      }
    }
    if (ANSWER == "Y")
      ANSWER <- readline("All value of the kernel function were formed. Press ENTER to continue.")

    if (ITERATION == "Y") {
      ANSWER <- readline("To form normalized value of the kernel function? (Y/N)")
      while(ANSWER != "N" && ANSWER != "Y")
        ANSWER <- readline()
    }
    if (ANSWER == "Y")
      cat(cli::rule(left = crayon::bold("Formation of normalized value of the kernel function"),
                    width = min(getOption("width"),40)), "\n")
    for (i in 1:n) {
      pNorm[i] <- p[i] / sum(p)
      if (ANSWER == "Y") {
        cat(paste0("pNorm (n = ", i, ") = "), pNorm[i], "\n")
        cat("+-----------------------------------+\n\n")
        ANSWER <- readline("")
        if (ANSWER != "N")
          ANSWER = "Y"
      }
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
    if (ITERATION == "Y") {
      cat("New center point value = ", x, "\n")
      cat("Modified increment     = ", delta, "\n")
      cat("+-----------------------------------+\n\n")
    }
    if (ITERATION == "Y"){
      ANSWER <- "N"
      ITERATION <- readline("Conduct a new detailed iteration or repeat current iteration? (Y/N/R)")
      while(ITERATION != "N" && ITERATION != "Y" && ITERATION != "R")
        ITERATION <- readline()
    }
    if (ITERATION != "R") {
      resultObj @allX = rbind(resultObj @allX, addRowForDF(k, x))
      resultObj @allDelta = rbind(resultObj @allDelta, addRowForDF(k, delta))
      k <- k + 1
      allX[k, ] <- x
      allResults <- c(allResults, f(x))
      flag <- (max(delta) <= e)
      if (flag == TRUE)
        break
    }
    else {
      ITERATION = "Y"
      x = allX[k,]
      delta <- pastDelta
    }
  }
  resultObj @iterations <- k
  resultObj @extremePoint <- allX[which.min(allResults), ]
  resultObj @fitnessValue <- allResults[which.min(allResults)]
  return(resultObj)
}
