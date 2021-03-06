sacStartupMessage <- function() {
  # Startup message obtained as
  # > SAC
  msg <- c(paste0("
\"```````````````````````````````````````````
`````````````````````````````````````````
`````_____``````````````````````____`````
```/``````\\````````````````````/`___`\\```
``/```__```\\````````__````````/`/```\\`\\``
`|```/``\\```|``````/``\\``````|`|`````|_|`
`|``|````|__|`````/`/\\`\\`````|`|`````````
`|``|````````````/`/``\\`\\````|`|`````````
``\\``\\``````````/`/````\\`\\```|`|`````````
````\\``\\```````/`/``````\\`\\``|`|`````````
``````\\``\\````|`|````````|`|`|`|`````````
````````\\``\\``|`|________|`|`|`|`````````
``__`````|``|`|`|________|`|`|`|`````````
`|``|````|``|`|`|````````|`|`|`|``````_``
`|```\\__/```|`|`|````````|`|`|`|`````|`|`
``\\````````/``|`|````````|`|``\\`\\___/`/``
```\\______/```|_|````````|_|```\\_____/```
`````````````````````````````````````````
`````````````````````````````````````````
````````````````````````````````````````````\" version ",
                  utils::packageVersion("SAC")),
    "\nType 'citation(\"SAC\")' for citing this R package in publications.")
  return(msg)
}

.onAttach <- function(lib, pkg) {
  # startup message
  msg <- sacStartupMessage()
  if(!interactive())
    msg[1] <- paste("Package 'SAC' version", utils::packageVersion("SAC"))
  packageStartupMessage(msg)
  invisible()
}

