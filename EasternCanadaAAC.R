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
# NEWFOUNDLAND PARSER
# =========================================================

parseNL <- function(sim,
                    maxAge = 205) {
  
  library(data.table)
  
  message("Preparing Newfoundland yield tables...")
  
  # =======================================================
  # curves selected by classifier
  # =======================================================
  
  curves_needed <- unique(
    sim$analysisUnitDT$bestCurve
  )
  
  curves_needed <- curves_needed[!is.na(curves_needed)]
  
  # =======================================================
  # output
  # =======================================================
  
  yieldTables <- list()
  
  # =======================================================
  # loop over curves
  # =======================================================
  
  for (curve_name in curves_needed) {
    
    message("Processing: ", curve_name)
    
    # ---------------------------------------------------
    # file path
    # ---------------------------------------------------
    
    yld_path <- file.path(
      "data/NL/YTF",
      paste0(curve_name, ".yld")
    )
    
    if (!file.exists(yld_path)) {
      warning("Missing: ", yld_path)
      next
    }
    
    # ---------------------------------------------------
    # read file
    # ---------------------------------------------------
    
    lines <- readLines(yld_path)
    
    # ---------------------------------------------------
    # parse species rows
    # ---------------------------------------------------
    
    species_vectors <- list()
    
    i <- 1
    
    while (i <= length(lines)) {
      
      line <- lines[i]
      
      if (grepl(
        "^[[:space:]]+[A-Za-z]{2,3}v",
        line
      )) {
        
        combined <- trimws(line)
        
        j <- i + 1
        
        # continuation lines
        while (
          j <= length(lines) &&
          grepl(
            "^[[:space:]]+[0-9]",
            lines[j]
          )
        ) {
          
          combined <- paste(
            combined,
            trimws(lines[j])
          )
          
          j <- j + 1
        }
        
        vals <- strsplit(
          combined,
          "\\s+"
        )[[1]]
        
        vals <- vals[vals != ""]
        
        if (length(vals) > 5) {
          
          sp <- vals[1]
          
          y <- suppressWarnings(
            as.numeric(vals[-c(1,2)])
          )
          
          y <- y[!is.na(y)]
          
          species_vectors[[sp]] <- y
        }
        
        i <- j
        
      } else {
        
        i <- i + 1
      }
    }
    
    # ---------------------------------------------------
    # equalize lengths
    # ---------------------------------------------------
    
    min_len <- min(
      sapply(species_vectors, length)
    )
    
    species_vectors <- lapply(
      species_vectors,
      function(x) x[1:min_len]
    )
    
    # ---------------------------------------------------
    # total volume
    # ---------------------------------------------------
    
    total_volume <- Reduce(
      "+",
      species_vectors
    )
    
    # ---------------------------------------------------
    # ages
    # ---------------------------------------------------
    
    ages <- seq(
      0,
      by = 5,
      length.out = length(total_volume)
    )
    
    # ---------------------------------------------------
    # annual interpolation
    # ---------------------------------------------------
    
    annual <- approx(
      x = ages,
      y = total_volume,
      xout = 1:maxAge,
      method = "linear",
      rule = 2
    )$y
    
    annual[annual < 0] <- 0
    
    # ---------------------------------------------------
    # save
    # ---------------------------------------------------
    
    yieldTables[[curve_name]] <- data.table(
      age = 1:maxAge,
      volume = annual
    )
  }
  
  return(yieldTables)
}


# =========================================================
# ONTARIO PARSER
# =========================================================

parseON <- function(sim,
                    maxAge = 255) {
  
  library(data.table)
  
  message("Preparing Ontario yield tables...")
  
  # =======================================================
  # curves selected by classifier
  # =======================================================
  
  curves_needed <- unique(
    sim$AUtoCurve$curveID
  )
  
  curves_needed <- curves_needed[
    !is.na(curves_needed)
  ]
  
  # =======================================================
  # Ontario YTF files
  # =======================================================
  
  ytf_files <- list.files(
    "data/ON/YTF",
    pattern = "tbl_yield_final",
    full.names = TRUE
  )
  
  # =======================================================
  # species columns
  # =======================================================
  
  species_cols <- c(
    "PW","PR","PJ",
    "SB","SW","BF",
    "CE","OC","HE",
    "PO","PB",
    "BW","MH","QR",
    "YB","OH"
  )
  
  # =======================================================
  # output object
  # =======================================================
  
  yieldTables <- list()
  
  # =======================================================
  # loop over curves
  # =======================================================
  
  for (curve_id in curves_needed) {
    
    message("Processing curve: ", curve_id)
    
    curve_found <- FALSE
    
    # ---------------------------------------------------
    # search all YTF files
    # ---------------------------------------------------
    
    for (f in ytf_files) {
      
      message("Checking file: ", basename(f))
      
      yt <- fread(f)
      
      # -------------------------------------------------
      # skip if curve absent
      # -------------------------------------------------
      
      if (!(curve_id %in% yt$CURVENO)) {
        next
      }
      
      curve_found <- TRUE
      
      # -------------------------------------------------
      # extract curve
      # -------------------------------------------------
      
      curve_dt <- yt[
        CURVENO == curve_id
      ]
      
      # -------------------------------------------------
      # sort by age
      # -------------------------------------------------
      
      setorder(curve_dt, AC10)
      
      # -------------------------------------------------
      # total volume
      # -------------------------------------------------
      
      curve_dt[, total_volume := rowSums(
        .SD,
        na.rm = TRUE
      ), .SDcols = species_cols]
      
      # -------------------------------------------------
      # annual interpolation
      # -------------------------------------------------
      
      annual <- approx(
        x = curve_dt$AC10,
        y = curve_dt$total_volume,
        xout = 1:maxAge,
        method = "linear",
        rule = 2
      )$y
      
      annual[annual < 0] <- 0
      
      # -------------------------------------------------
      # save standardized curve
      # -------------------------------------------------
      
      yieldTables[[as.character(curve_id)]] <- data.table(
        age = 1:maxAge,
        volume = annual
      )
      
      # -------------------------------------------------
      # curve found → stop searching
      # -------------------------------------------------
      
      break
    }
    
    # ---------------------------------------------------
    # warning if missing
    # ---------------------------------------------------
    
    if (!curve_found) {
      
      warning(
        "Curve not found in Ontario YTF files: ",
        curve_id
      )
    }
  }
  
  return(yieldTables)
}


# =========================================================
# WRAPPER
# =========================================================

prepareYieldTables <- function(sim,
                               maxAge = 255) {
  
  jur <- toupper(
    P(sim)$jurisdiction
  )
  
  if (jur == "NL") {
    return(parseNL(sim, maxAge))
  }
  
  if (jur == "ON") {
    return(parseON(sim, maxAge))
  }
  
  stop(
    "Unsupported jurisdiction: ",
    jur
  )
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
  # MAI = V(t) / t (Mean Annual Increment)
  vt <- yt / 1:n
  
  # Rotation age (R): age at maximum MAI
  R <- order(vt, decreasing = TRUE)[1]
  
  # Increment (annual volume growth)
  inc <- c(0, diff(yt))
  
  # Normalize mature volume by rotation period
  # This distributes harvestable volume over time
  tmp <- yt[R:n] /
    (P(sim)$rotationPeriodMultiplier *
       (R + P(sim)$rotationPeriodShift))  
  # Combine:
  # - early ages use increment
  # - older ages use normalized volume contribution
  tmp <- c(inc[1:R-1], tmp)
  
  # Return parameters:
  # R = rotation age
  # I = increment vector
  # hVec = contribution vector used in AAC calculation
  return(list(R = R, I = inc, V = yt, hVec = tmp))
}


# =========================================================
# Initialization function
# =========================================================
Init <- function(sim) {
  
  # =====================================================
  # Prepare annual yield curves
  # =====================================================
  
  sim$yieldTables <- prepareYieldTables(sim)
  
  # =====================================================
  # Compute Hanzlik parameters
  # =====================================================
  
  sim$hanzlikPars <- lapply(
    sim$yieldTables,
    function(x) {
      calcHanzlik(
        x$volume,
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
  
  if (!is.null(sim$pixelAreaDT)) {
    
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
      
      
      ages <- 1:100
      
      yt_interp <- 0.05 * ages^2 * exp(-0.03 * ages)
      yt_interp <- yt_interp / max(yt_interp) * 300  # scale
      # realistic growth curve  
    }
    
    # -------------------------------------------------------
    # Build yield table list for each AU
    # -------------------------------------------------------
    # IMPORTANT:
    # Each Analysis Unit must have its own yield curve entry
    sim$AUtoCurve <- data.table(
      AU = as.character(AUvals),
      curveID = as.character(AUvals)
    )
    
    fakeYT <- data.table(
      age = 1:100,
      volume = yt_interp
    )
    
    sim$yieldTables <- setNames(
      replicate(length(AUvals), fakeYT, simplify = FALSE),
      sim$AUtoCurve$curveID
    )
  }
  
  return(invisible(sim))
}