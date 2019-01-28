RubanStartupMessage <- function() {
  # Startup message obtained as
  # > Ruban
  msg <- c(paste0("
\"```````````````````````````````````````````
```___````_````_````__````````__````_````_``Ruban
``/   \\``|`|``|`|``|`_`\\`````/``|``| \\``|`|`   Algorithm
`|`````|`|`|` |`|` |___/````/```|``|``\\`|`|`
`|____/``|`|``|`|``|```\\```/____|``| |\\\\|`|`
`|````\\``|`|__|`|``|`__`|`/-----|``| |`\\``|`
`|`````\\``\\````/```|___/`/``````|``|`|``\\`|`
````````````````````````````````````````````\" version ",
                  utils::packageVersion("Ruban")),
    "\nType 'citation(\"Ruban\")' for citing this R package in publications.")
  return(msg)
}

.onAttach <- function(lib, pkg) {
  # startup message
  msg <- RubanStartupMessage()
  if (!interactive())
    msg[1] <- paste("Package 'Ruban' version", utils::packageVersion("Ruban"))
  packageStartupMessage(msg)
  invisible()
}
