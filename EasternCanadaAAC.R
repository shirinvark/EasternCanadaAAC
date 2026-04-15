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

  # =========================================================
  # Hanzlik calculation function
  # =========================================================
  calcHanzlik <- function(ytM, sim){
    
    # If yield table is a matrix, take the last column (assumed volume)
    if (is.matrix(ytM)){
      idx <- dim(ytM)[2]
      yt <- ytM[, idx]
      
      # If already numeric vector, use as is
    } else if (is.numeric(ytM)){
      yt = ytM
      
    } else {
      stop("illegal yield table: matrix expected")
    }
    
    # Number of age classes
    n <- length(yt)
    
    # Mean annual increment (MAI)
    vt <- yt / 1:n
    
    # Rotation age (R): age at maximum MAI
    R <- order(vt, decreasing = TRUE)[1]
    
    # Increment (annual volume growth)
    inc <- c(0, diff(yt))
    
    # Contribution of each age class after rotation
    # (based on Hanzlik sustained yield formulation)
    tmp <- yt[R:n] / (P(sim)$rotationPeriodMultiplier * R)
    
    # Combine:
    # - early ages use increment
    # - older ages use normalized volume contribution
    tmp <- c(inc[1:R-1], tmp)
    
    # Return parameters:
    # R = rotation age
    # I = increment vector
    # hVec = contribution vector used in AAC calculation
    return(list(R = R, I = inc, hVec = tmp))
  }
  
  
  # =========================================================
  # Initialization function
  # =========================================================
  Init <- function(sim) {
    
    # Compute Hanzlik parameters for each yield table
    sim$hanzlikPars <- lapply(sim$yieldTables, calcHanzlik, sim)
    
    # 🔥 IMPORTANT:
    # Ensure names match analysisUnit IDs exactly
    names(sim$hanzlikPars) <- names(sim$yieldTables)
    
    return(invisible(sim))
  }
  
  
  # =========================================================
  # Plan function (AAC calculation)
  # =========================================================
  Plan <- function(sim) {
    
    # -------------------------------------------------------
    # 1. Extract Analysis Unit (AU) and stand age values
    # -------------------------------------------------------
    AUvals  <- as.vector(terra::values(sim$analysisUnitMap))
    ageVals <- as.vector(terra::values(sim$standAgeMap))
    
    dt <- data.table(
      AU  = AUvals,
      age = ageVals
    )
    
    # Remove missing values
    dt <- dt[!is.na(AU) & !is.na(age)]
    
    # -------------------------------------------------------
    # 2. Prepare data
    # -------------------------------------------------------
    
    # Convert age to integer (required for indexing)
    dt[, age := floor(age)]
    
    # 🔴 CRITICAL CHECK:
    # Ensure every AU has a corresponding yield table
    if (!all(unique(dt$AU) %in% names(sim$hanzlikPars))) {
      stop("Some AU values do not match yieldTables!")
    }
    
    # -------------------------------------------------------
    # 3. Compute AAC per Analysis Unit
    # -------------------------------------------------------
    AAC_by_AU <- dt[, .(
      AAC = {
        
        # Retrieve Hanzlik parameters for this AU
        pars <- sim$hanzlikPars[[as.character(.BY$AU)]]
        
        # Clamp ages to valid range of hVec
        ages <- pmin(pmax(age, 1), length(pars$hVec))
        
        # Sum contributions across pixels
        sum(pars$hVec[ages])
      }
    ), by = AU]
    
    # -------------------------------------------------------
    # 4. Compute cell area (hectares)
    # -------------------------------------------------------
    cellArea_ha <- prod(terra::res(sim$standAgeMap)) / 10000
    
    # -------------------------------------------------------
    # 5. Total AAC (m3/year)
    # -------------------------------------------------------
    totalAAC <- sum(AAC_by_AU$AAC) * cellArea_ha
    
    # -------------------------------------------------------
    # 6. Store result in sim object
    # -------------------------------------------------------
    sim$AAC <- totalAAC
    
    # Print result
    message(sprintf(
      "Year %d | AAC = %.2f m3/year",
      time(sim), totalAAC
    ))
    
    return(invisible(sim))
  }

# =========================================================
# .inputObjects function
# Responsible for ensuring all required inputs exist in sim
# If not provided, creates fallback (test) data
# =========================================================
.inputObjects <- function(sim) {
  
  # -------------------------------------------------------
  # Resolve data path for module inputs
  # -------------------------------------------------------
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")
  
  # =========================================================
  # analysisUnitMap
  # =========================================================
  
  # Check if analysisUnitMap is already provided and valid
  if (!is.null(sim$analysisUnitMap) && inherits(sim$analysisUnitMap, "SpatRaster")) {
    
    message("Using supplied analysisUnitMap")
    
  } else {
    
    # Create a small dummy raster for testing
    message("Creating fake analysisUnitMap...")
    
    analysisUnitMap <- terra::rast(
      nrows = 10, ncols = 10,
      xmin = 0, xmax = 1000,
      ymin = 0, ymax = 1000
    )
    
    # Assign random AU IDs (1–3)
    terra::values(analysisUnitMap) <- sample(1:3, 100, replace = TRUE)
    
    sim$analysisUnitMap <- analysisUnitMap
  }
  
  # =========================================================
  # standAgeMap
  # =========================================================
  
  # Check if standAgeMap is already provided
  if (!is.null(sim$standAgeMap)) {
    
    message("Using supplied standAgeMap")
    
  } else {
    
    # Create a dummy stand age raster
    message("Creating fake standAgeMap...")
    
    standAgeMap <- terra::rast(
      nrows = 10, ncols = 10,
      xmin = 0, xmax = 1000,
      ymin = 0, ymax = 1000
    )
    
    # Assign random stand ages (1–80 years)
    terra::values(standAgeMap) <- sample(1:80, 100, replace = TRUE)
    
    sim$standAgeMap <- standAgeMap
  }
  
  # =========================================================
  # yieldTables
  # =========================================================
  
  # If yieldTables already exist in sim, use them
  if ("yieldTables" %in% names(sim)) {
    
    message("Using supplied yieldTables")
    
  } else {
    
    message("No yieldTables supplied → building internally")
    
    # -------------------------------------------------------
    # Extract unique Analysis Unit IDs from raster
    # -------------------------------------------------------
    AUvals <- unique(terra::values(sim$analysisUnitMap))
    AUvals <- AUvals[!is.na(AUvals)]
    
    # -------------------------------------------------------
    # Attempt to read real .yld file
    # -------------------------------------------------------
    yldFile <- file.path(dPath, "NL/YTF/BarNS_sub_all.yld")
    
    if (file.exists(yldFile)) {
      
      message("Reading real .yld file...")
      
      # Read raw file lines
      lines <- readLines(yldFile)
      
      # Locate header line containing "Age"
      start <- grep("Age", lines)[1]
      
      # Read table from file (after header)
      tab <- read.table(
        text = lines[(start + 1):length(lines)],
        header = FALSE,
        fill = TRUE
      )
      
      # Assign column names (based on NL yield table structure)
      colnames(tab) <- c("Age", "BSv", "WSv", "BFv", "WPv", "TLv", "YBv")
      
      # Select one yield curve (e.g., black spruce volume)
      yt <- tab$BSv
      
      # Interpolate yield curve to annual resolution (1–100 years)
      yt_interp <- approx(
        x = tab$Age,
        y = yt,
        xout = 1:100
      )$y
      
      # Replace NA values (from interpolation gaps) with zero
      yt_interp[is.na(yt_interp)] <- 0
      
      message("Using .yld-based yield curve")
      
    } else {
      
      # Fallback: create synthetic yield curve
      message("No .yld file found → using fake yield curve")
      
      yt_interp <- cumsum(seq(1, 3, length.out = 100))
    }
    
    # -------------------------------------------------------
    # Build yield table list for each AU
    # -------------------------------------------------------
    # IMPORTANT:
    # Each Analysis Unit must have its own yield curve entry
    sim$yieldTables <- setNames(
      replicate(length(AUvals), yt_interp, simplify = FALSE),
      as.character(AUvals)
    )
  }
  
  return(invisible(sim))
}