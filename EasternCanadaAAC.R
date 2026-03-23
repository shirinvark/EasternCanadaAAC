## Everything in this file and any files in the R directory are sourced during `simInit()`;
## all functions and objects are put into the `simList`.
## To use objects, use `sim$xxx` (they are globally available to all modules).
## Functions can be used inside any function that was sourced in this module;
## they are namespaced to the module, just like functions in R packages.
## If exact location is required, functions will be: `sim$.mods$<moduleName>$FunctionName`.
defineModule(sim, list(
  name = "EasternCanadaAAC",
  description = "",
  keywords = "",
  authors = structure(list(list(given = c("First", "Middle"), family = "Last", role = c("aut", "cre"), email = "email@example.com", comment = NULL)), class = "person"),
  childModules = character(0),
  version = list(EasternCanadaAAC = "0.0.0.9000"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("NEWS.md", "README.md", "EasternCanadaAAC.Rmd"),
  reqdPkgs = list("SpaDES.core (>= 3.0.4)", "ggplot2"),
  parameters = bindrows(
    
    defineParameter(
      name = "replanInterval",
      class = "numeric",
      default = 10,
      desc = "Years between AAC recalculation"
    ),
    
    defineParameter(
      name = ".plots",
      class = "character",
      default = "screen",
      desc = "Plotting option"
    )
  ),
  inputObjects = bindrows(
    #expectsInput("objectName", "objectClass", "input object description", sourceURL, ...),
    expectsInput(objectName = NA, objectClass = NA, desc = NA, sourceURL = NA)
  ),
  outputObjects = bindrows(
    #createsOutput("objectName", "objectClass", "output object description", ...),
    createsOutput(objectName = NA, objectClass = NA, desc = NA)
  )
))

doEvent.EasternCanadaAAC = function(sim, eventTime, eventType, debug = FALSE) {
  switch(
    eventType,
    init = {
      sim <- Init(sim)
      sim <- scheduleEvent(sim, time(sim), "EasternCanadaAAC", "plan")
    },
    plan = {
      sim <- Plan(sim)
      sim <- scheduleEvent(sim, time(sim) + P(sim)$replanInterval, "EasternCanadaAAC", "plan")
    },
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}



calcHanzlik <- function(ytM, sim){

  if (is.matrix(ytM)){
    idx <- dim(ytM)[2]
    yt <- ytM[idx]
  }
  else if (is.numeric(ytM)){
    yt = ytM
  }
  else
    stop("illegal yield table: matrix expected")
  n <- length(yt)
  vt<- yt/1:n
  #plot(vt)
  R <- order(vt,decreasing=TRUE)[1] #trouve le age de la culmination
  inc <- c(0, diff(yt))
  tmp <- yt[R:n]/(P(sim)$rationPeriodMultiplier*R)     #contribution of each age class to Vm/R in m^3/ha
  tmp <- c(inc[1:R-1], tmp)
  
  return(list(R=R,I=inc,hVec=tmp))   #should not alter Vm after initial plan?
  
}

### template initialization
Init <- function(sim) {
  
  sim$hanzlikPars <- lapply(sim$yieldTables, calcHanzlik, sim) # we assume there is a list these structures
  
  return(invisible(sim))
}


Plan <- function(sim) {
  #browser()
  print(sim$analysisUnitMap)
  AUvals <- terra::values(sim$analysisUnitMap)
  print(head(AUvals))
  hVec<-sim$hanzlikPars[[1]]$hVec  
  sim$cellSize <- prod(terra::res(sim$standAgeMap))
  #for now, assume annualCut has at least one object, and take the 1st
  #Comment peut-on le generaliser pour plusieurs courbes de rendemment?
  nh <- length(hVec)
  res<-rep(0,nh)
  ages <- sim$standAgeMap[]
  ages <- ifelse(ages > nh, nh, ages) #silently truncate
  x<-tabulate(ages)  #NAs are silently ignored   
  nx<-length(x)
  if (nx <= nh){
    res[1:nx] <- x
  }
  else {
    res <- x[1:nh]
    x[nh] <- x[nh] + sum(x[nh+1:nx]) #accumulate any "missing ones"
  }
  
  res<-sum(res*hVec) * sim$cellSize  #OOPS
  message(sprintf("Hanzlik:Plan. %d AAC = %5.1f x 1000 m^3\n",as.integer(time(sim)),res/1e3))
  sim$annualCut <- res  
  return(invisible(sim))
}

.inputObjects <- function(sim) {
 
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")

  # ! ----- EDIT BELOW ----- ! #

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

