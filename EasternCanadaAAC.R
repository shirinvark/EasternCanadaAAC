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
      name = "rotationPeriodMultiplier",
      class = "numeric",
      default = 1,
      desc = "Multiplier for Hanzlik rotation"
    ),
    
    defineParameter(
      name = ".plots",
      class = "character",
      default = "screen",
      desc = "Plotting option"
    )
  ),
    inputObjects = bindrows(
      
      expectsInput("analysisUnitMap", "SpatRaster",
                   "Map of analysis units per pixel"),
      
      expectsInput("standAgeMap", "SpatRaster",
                   "Stand age per pixel"),
      
      expectsInput("yieldTables", "list",
                   "Yield tables per AU")
  ),
    outputObjects = bindrows(
  createsOutput("AAC", "numeric",
                "Annual Allowable Cut (m3/year)")

))
)
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
    yt <- ytM[, idx]
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
  tmp <- yt[R:n]/(P(sim)$rotationPeriodMultiplier * R)    #contribution of each age class to Vm/R in m^3/ha
  tmp <- c(inc[1:R-1], tmp)
  
  return(list(R=R,I=inc,hVec=tmp))   #should not alter Vm after initial plan?
  
}

### template initialization
Init <- function(sim) {
  
  sim$hanzlikPars <- lapply(sim$yieldTables, calcHanzlik, sim)
  
  # 🔥 این خط خیلی مهمه
  names(sim$hanzlikPars) <- names(sim$yieldTables)
  
  return(invisible(sim))
}


Plan <- function(sim) {
  
  # 🔹 1. گرفتن AU و سن
  AUvals  <- as.vector(terra::values(sim$analysisUnitMap))
  ageVals <- as.vector(terra::values(sim$standAgeMap))
  
  dt <- data.table(
    AU  = AUvals,
    age = ageVals
  )
  
  dt <- dt[!is.na(AU) & !is.na(age)]
  
  # 🔹 2. گروه‌بندی بر اساس AU
  dt[, age := floor(age)]
  # 🔴 چک حیاتی: آیا همه AUها yield دارن؟
  if (!all(unique(dt$AU) %in% names(sim$hanzlikPars))) {
    stop("Some AU values do not match yieldTables!")
  }
  # 🔹 3. محاسبه AAC برای هر AU
  AAC_by_AU <- dt[, .(
    AAC = {
      pars <- sim$hanzlikPars[[as.character(.BY$AU)]]
      ages <- pmin(pmax(age, 1), length(pars$hVec))
      sum(pars$hVec[ages])
    }
  ), by = AU]
  
  # 🔹 4. مساحت هر پیکسل (ha)
  cellArea_ha <- prod(terra::res(sim$standAgeMap)) / 10000
  
  # 🔹 5. AAC کل
  totalAAC <- sum(AAC_by_AU$AAC) * cellArea_ha  
  # 🔹 6. ذخیره
  sim$AAC <- totalAAC
  
  message(sprintf(
    "Year %d | AAC = %.2f m3/year",
    time(sim), totalAAC
  ))
  
  return(invisible(sim))
}


.inputObjects <- function(sim) {
  
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")
  
  # =========================================================
  # analysisUnitMap
  # =========================================================
  
  if (!is.null(sim$analysisUnitMap) && inherits(sim$analysisUnitMap, "SpatRaster")) {
    
    message("Using supplied analysisUnitMap")
    
  } else {
    
    message("Creating fake analysisUnitMap...")
    
    analysisUnitMap <- terra::rast(
      nrows = 10, ncols = 10,
      xmin = 0, xmax = 1000,
      ymin = 0, ymax = 1000
    )
    
    terra::values(analysisUnitMap) <- sample(1:3, 100, replace = TRUE)
    
    sim$analysisUnitMap <- analysisUnitMap
  }
  # =========================================================
  # standAgeMap
  # =========================================================
  
  if (!is.null(sim$standAgeMap)) {
    
    message("Using supplied standAgeMap")
    
  } else {
    
    message("Creating fake standAgeMap...")
    
    standAgeMap <- terra::rast(
      nrows = 10, ncols = 10,
      xmin = 0, xmax = 1000,
      ymin = 0, ymax = 1000
    )
    
    terra::values(standAgeMap) <- sample(1:80, 100, replace = TRUE)
    
    sim$standAgeMap <- standAgeMap
  }
  
  # =========================================================
  # yieldTables
  # =========================================================
  
  if ("yieldTables" %in% names(sim)) {
    
    message("Using supplied yieldTables")
    
  } else {
    
    message("No yieldTables supplied → creating based on AU")
    
    # 🔹 اول AU ها رو بگیر
    AUvals <- unique(terra::values(sim$analysisUnitMap))
    AUvals <- AUvals[!is.na(AUvals)]
    
    # 🔹 یه yield ساده بساز (یا از .yld اگر خواستی)
    yt <- cumsum(seq(1, 3, length.out = 100))
    
    # 🔹 برای هر AU یک curve بساز
    sim$yieldTables <- setNames(
      replicate(length(AUvals), yt, simplify = FALSE),
      as.character(AUvals)
    )
  }
  }
  
  return(invisible(sim))
}