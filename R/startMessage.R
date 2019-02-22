RoubanStartupMessage <- function() {
  # Startup message obtained as
  # > Rouban
  msg <- c(paste0("
\"```````````````````````````````````````````
```___`````____`````_````_````__````````__````_````_``Rouban
``/   \\```/`__`\\```|`|``|`|``|`_`\\`````/``|``| \\``|`|`   Algorithm
`|`````|`|`/``\\`|``|`|` |`|` |___/````/```|``|``\\`|`|`
`|____/``||````||``|`|``|`|``|```\\```/____|``| |\\\\|`|`
`|````\\``|`\\__/`|``|`|__|`|``|`__`|`/-----|``| |`\\``|`
`|`````\\``\\____/````\\____/```|___/`/``````|``|`|``\\`|`
````````````````````````````````````````````\" version ",
                  utils::packageVersion("Rouban")),
    "\nType 'citation(\"Rouban\")' for citing this R package in publications.")
  return(msg)
}

.onAttach <- function(lib, pkg) {
  # startup message
  msg <- RoubanStartupMessage()
  if(!interactive())
    msg[1] <- paste("Package 'Rouban' version", utils::packageVersion("Rouban"))
  packageStartupMessage(msg)
  invisible()
}

