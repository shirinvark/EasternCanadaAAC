## Everything in this file and any files in the R directory are sourced during `simInit()`;
## all functions and objects are put into the `simList`.
## To use objects, use `sim$xxx` (they are globally available to all modules).
## Functions can be used inside any function that was sourced in this module;
## they are namespaced to the module, just like functions in R packages.
## If exact location is required, functions will be: `sim$.mods$<moduleName>$FunctionName`.
defineModule(sim, list(
  name = "EasternCanadaAAC",
  description = "Compute Annual Allowable Cut (AAC) using Hanzlik method",
  keywords = c("AAC", "Hanzlik", "forest management", "SpaDES"),
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
      name = "rotationPeriodShift",
      class = "numeric",
      default = 0,
      desc = "Years added to rotation age"
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
      expectsInput("pixelGroupMap", "SpatRaster", "Pixel groups"),
      expectsInput("pixelAreaDT", "data.table", "Effective area"),
      
      expectsInput("standAgeMap", "SpatRaster",
                   "Stand age per pixel"),
      expectsInput("AUtoCurve", "data.table",
                   "Mapping between analysis units and yield curves"),
      expectsInput("yieldTables", "list",
                   "Yield tables by curve"),
  ),
    outputObjects = bindrows(
  createsOutput("AAC", "numeric",
                "Annual Allowable Cut (m3/year)"),
  createsOutput("AAC_by_AU", "data.table", "AAC per AU"),
  createsOutput("hanzlikPars", "list", "Hanzlik parameters")

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
# Hanzlik parameter calculation
# ---------------------------------------------------------
# Converts a yield table into:
# - annual volume (V)
# - increment (I)
# - rotation age (R) based on maximum MAI
#
# This function prepares all parameters required for
# AAC calculation following the Hanzlik method.
# =========================================================  
calcHanzlik <- function(yt, sim){
  
  # -------------------------------------------------------
  # check input
  # -------------------------------------------------------
  
  if (!is.numeric(yt)) {
    stop("Yield table must be numeric")
  }
  
  # -------------------------------------------------------
  # number of age classes
  # -------------------------------------------------------
  
  n <- length(yt)
  
  # -------------------------------------------------------
  # Mean Annual Increment (MAI)
  # -------------------------------------------------------
  
  vt <- yt / (1:n)
  
  # -------------------------------------------------------
  # Rotation age (maximum MAI)
  # -------------------------------------------------------
  
  R <- which.max(vt)
  
  # -------------------------------------------------------
  # annual increment
  # -------------------------------------------------------
  
  inc <- c(
    0,
    diff(yt)
  )
  
  # -------------------------------------------------------
  # mature forest contribution
  # -------------------------------------------------------
  
  mature_part <- yt[R:n] /
    (
      P(sim)$rotationPeriodMultiplier *
        (
          R +
            P(sim)$rotationPeriodShift
        )
    )
  
  # -------------------------------------------------------
  # immature contribution
  # -------------------------------------------------------
  
  young_part <- if (R > 1) {
    inc[1:(R - 1)]
  } else {
    numeric(0)
  }
  
  # -------------------------------------------------------
  # Hanzlik vector
  # -------------------------------------------------------
  
  hVec <- c(
    young_part,
    mature_part
  )
  
  # -------------------------------------------------------
  # return
  # -------------------------------------------------------
  
  return(
    list(
      R = R,
      I = inc,
      V = yt,
      hVec = hVec
    )
  )
}

# =========================================================
# Initialization function
# =========================================================
Init <- function(sim) {
  
  # =====================================================
  # Prepare annual yield curves
  # =====================================================
  
  if (is.null(sim$yieldTables)) {
    
    yt <- prepareYieldTables(sim)
    
    if (!is.null(yt)) {
      sim$yieldTables <- yt
    }
  }  
  # =====================================================
  # Compute Hanzlik parameters
  # =====================================================
  
  sim$hanzlikPars <- lapply(
    sim$yieldTables,
    function(x) {
      
      # ---------------------------------------------------
      # ON tables
      # ---------------------------------------------------
      
      if ("volume" %in% names(x)) {
        
        yt <- x$volume
        
      } else {
        
        # -------------------------------------------------
        # NL species-wise tables
        # -------------------------------------------------
        
        numeric_cols <- names(x)[
          sapply(x, is.numeric)
        ]
        
        numeric_cols <- setdiff(
          numeric_cols,
          "age"
        )
        
        if (is.data.table(x)) {
          
          yt <- rowSums(
            x[, numeric_cols, with = FALSE],
            na.rm = TRUE
          )
          
        } else {
          
          yt <- rowSums(
            as.data.frame(x)[, numeric_cols, drop = FALSE],
            na.rm = TRUE
          )
        }
      }
      
      calcHanzlik(
        yt,
        sim
      )
    }
  )
  # =====================================================
  # Ensure names match analysis units
  # =====================================================
  
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
    pixelGroup = as.vector(terra::values(sim$pixelGroupMap)),
    AU  = AUvals,
    age = ageVals
  )
  # -------------------------------------------------------
  # 🔥 join effective area from classifier
  # -------------------------------------------------------
  dt <- merge(
    dt,
    sim$pixelAreaDT[, .(pixelGroup, effectiveArea)],
    by = "pixelGroup",
    all.x = TRUE
  )
  
  # Check for missing effective area (should not happen)
  if (any(is.na(dt$effectiveArea))) {
    stop("❌ NA in effectiveArea — join with pixelAreaDT failed")
  }
  
  # Remove missing values
  dt <- dt[!is.na(AU) & !is.na(age)]
  
  # -------------------------------------------------------
  # 2. Prepare data
  # -------------------------------------------------------
  
  # Convert age to integer for indexing yield vectors
  dt[, age := floor(age)]
  
  # 🔴 CRITICAL CHECK:
  # Ensure every AU has a corresponding yield table
  curveIDs <- unique(sim$AUtoCurve$curveID)
  print(curveIDs)
  
  print(names(sim$yieldTables))
  
  print(names(sim$hanzlikPars))
  if (!all(curveIDs %in% names(sim$hanzlikPars))) {
    stop("Some curveIDs do not match yieldTables!")
  }
  
  
  # -------------------------------------------------------
  # 3. Compute AAC per Analysis Unit (Hanzlik method)
  # -------------------------------------------------------
  # For each AU:
  # - Mature stands (age >= R): contribute via standing volume
  # - Immature stands (age < R): contribute via growth
  #
  # AAC = (Mature Volume / Rotation Period) + Growth
  # -------------------------------------------------------
  AAC_by_AU <- dt[, .(
    AAC = {
      
      # Retrieve Hanzlik parameters for this AU
      curveID <- sim$AUtoCurve[
        AU == .BY$AU
      ]$curveID
      
      pars <- sim$hanzlikPars[[curveID]]
      # Rotation age derived from maximum MAI (Mean Annual Increment)
      R <- pars$R
      if (isTRUE(getOption("aac.debug", TRUE))) {
        message("AU: ", .BY$AU, " | Rotation age: ", R)
      }        
      # Clamp ages to valid range
      ages <- pmin(pmax(age, 1), length(pars$V))
      ages <- as.integer(ages)
      
      # ----------------------------------------
      # Split stands based on rotation age
      # ----------------------------------------
      # Mature stands: age >= R → available for harvest
      # Immature stands: age < R → still growing
      mature_idx   <- ages >= R
      immature_idx <- ages < R
      
      # Check area exists
      if (!"effectiveArea" %in% names(.SD)) {
        stop("❌ effectiveArea missing — classifier not connected properly")
      }
      
      a <- effectiveArea
      
      # ----------------------------------------
      # Compute mature volume
      # ----------------------------------------
      # Total standing volume in stands older than rotation age
      V_mature <- sum(pars$V[ages[mature_idx]] * a[mature_idx])
      # ----------------------------------------
      # Compute growth from immature stands
      # ----------------------------------------
      # Total annual increment from younger stands
      I_total  <- sum(pars$I[ages[immature_idx]] * a[immature_idx])
      
      # ----------------------------------------
      # Hanzlik AAC formula
      # ----------------------------------------
      # AAC = (Mature Volume / Rotation Period) + Growth
      #
      # - Mature volume is distributed over the rotation period
      # - Growth from immature stands is added directly
      #
      # This ensures sustained yield over time
      (V_mature /
         (P(sim)$rotationPeriodMultiplier *
            (R + P(sim)$rotationPeriodShift))) +
        I_total      
    }
  ), by = AU]
  
  # -------------------------------------------------------
  # 4. Compute cell area (hectares)
  # -------------------------------------------------------
  #cellArea_ha <- prod(terra::res(sim$standAgeMap)) / 10000
  
  # -------------------------------------------------------
  # 5. Total AAC (m3/year)
  # -------------------------------------------------------
  totalAAC <- sum(AAC_by_AU$AAC)    
  # -------------------------------------------------------
  # 6. Store result in sim object
  # -------------------------------------------------------
  sim$AAC <- totalAAC
  sim$AAC_by_AU <- AAC_by_AU
  
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
  if (!is.null(sim$standAgeMap) &&
      inherits(sim$standAgeMap, "SpatRaster")) {
    
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
  # pixelGroupMap (fallback)
  # =========================================================
  
  if (!is.null(sim$pixelGroupMap) && inherits(sim$pixelGroupMap, "SpatRaster")) {
    
    message("Using supplied pixelGroupMap")
    
  } else {
    
    message("Creating fake pixelGroupMap...")
    
    r <- sim$analysisUnitMap
    
    terra::values(r) <- 1:terra::ncell(r)
    
    sim$pixelGroupMap <- r
  }
  # =========================================================
  # pixelAreaDT (fallback)
  # =========================================================
  
  if (!is.null(sim$pixelAreaDT) &&
      inherits(sim$pixelAreaDT, "data.table")) {
    
    message("Using supplied pixelAreaDT")
    
  } else {
    
    message("Creating fake pixelAreaDT...")
    
    cellArea_ha <- prod(terra::res(sim$analysisUnitMap)) / 10000
    
    pg <- as.vector(terra::values(sim$pixelGroupMap))
    
    sim$pixelAreaDT <- data.table(
      pixelGroup = pg,
      #effectiveArea = cellArea_ha
      # TEST ONLY: reduce effective area by 50% to test AAC sensitivity
      effectiveArea = cellArea_ha * 0.5   
    )
  }
  # =========================================================
  # yieldTables
  # =========================================================
  
  # If yieldTables already exist in sim, use them
  if (!is.null(sim$yieldTables)) {
    
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
      
      # Standardize yield curve to annual resolution
      fakeYT <- standardizeYieldCurve(
        ages = tab$Age,
        volumes = yt,
        maxAge = 200
      )
      
      message("Using .yld-based yield curve")
      
    } else {
      
      # Fallback: create synthetic yield curve
      message("No .yld file found → using fake yield curve")
      
      
      ages <- 1:100
      
      yt_interp <- 0.05 * ages^2 * exp(-0.03 * ages)
      yt_interp <- yt_interp / max(yt_interp) * 300 
      fakeYT <- standardizeYieldCurve(
        ages = ages,
        volumes = yt_interp,
        maxAge = 200
      )
      # scale
      # realistic growth curve  
    }
    
    # -------------------------------------------------------
    # Build yield table list for each AU
    # -------------------------------------------------------
    # IMPORTANT:
    # Each Analysis Unit must have its own yield curve entry
    if (!is.null(sim$AUtoCurve) &&
        inherits(sim$AUtoCurve, "data.table"))  {
      
      message("Using supplied AUtoCurve")
      
    } else {
      
      sim$AUtoCurve <- data.table(
        AU = as.character(AUvals),
        curveID = as.character(AUvals)
      )
    }
    
    
    sim$yieldTables <- setNames(
      replicate(length(AUvals), fakeYT, simplify = FALSE),
      sim$AUtoCurve$curveID
    )
  }
  
  return(invisible(sim))
}