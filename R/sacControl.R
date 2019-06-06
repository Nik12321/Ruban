#-------------------------------------------------------------------#
#           Function to search and configure advanced options       #
#-------------------------------------------------------------------#

sacControl <- function(...)
{
  current <- .sac.default

  if(nargs() == 0)
    return(current)

  args <- list(...)

  if(length(args) == 1 && is.null(names(args))) {
    arg <- args[[1]]
    switch(mode(arg),
           list = args <- arg,
           character = return(.sac.default[[arg]]),
           stop("invalid argument: ", dQuote(arg)))
  }

  if (length(args) == 0)
    return(current)

  nargs <- names(args)

  if (is.null(nargs))
    stop("options must be given by name")

  if(is.list(args)) {
    changed <- current[nargs]
  for(i in 1:length(nargs)) {
    if(is.list(args[[i]])) {
      what <- names(args[[i]])
      changed[[i]][what] <- args[[i]][what]
      }
    else {
      changed[i] <- args[[i]]
      }
  }
  current[nargs] <- changed
  }
  else {
    changed <- current[nargs]
    current[nargs] <- args
  }

  if(sys.parent() == 0)
    env <- asNamespace("SAC")
  else
    env <- parent.frame()
  assign(".sac.default", current, envir = env)
  invisible(current)
}

x <- function (...) {
  y =   nargs()
  return(y)
}
.sac.default <- list("sacNormal" = list(kernelType = "kernelExponential",
                                    n = 500,
                                    e = 0.01,
                                    r  = 2,
                                    s  = 100,
                                    y = 1,
                                    q = 2,
                                    M = 150),
                     "sacExtended" = list(kernelType = "kernelExponential",
                                        n = 500,
                                        e = 0.01,
                                        r  = 2,
                                        s  = 100,
                                        y = 1,
                                        q = 2,
                                        M = 150),
                     "sacIterative" = list(kernelType = "kernelExponential",
                                        n = 500,
                                        e = 0.01,
                                        r  = 2,
                                        s  = 100,
                                        y = 1,
                                        q = 2,
                                        M = 150))
