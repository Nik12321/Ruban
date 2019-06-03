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
