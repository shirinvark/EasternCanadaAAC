# =========================================================
# ONTARIO PARSER
# =========================================================
# Purpose:
#   Reads raw Ontario yield table (YTF) files,
#   extracts only the curves selected by the
#   classifier, converts species-level yields
#   into total stand volume, and standardizes
#   the curves into annual age increments
#   required by the AAC/Hanzlik workflow.
#
# General workflow:
#
#   AU classifier
#          ↓
#   selected curveIDs
#          ↓
#   search Ontario YTF files
#          ↓
#   extract matching curve
#          ↓
#   sum species volumes
#          ↓
#   standardize annual yield curve
#          ↓
#   AAC-ready yield tables
#
# Output:
#   A named list of standardized annual
#   yield tables indexed by curveID.
#
# Notes:
#   - Ontario YTF files store yields by species.
#   - AAC currently uses total stand volume.
#   - Interpolation is handled before AAC
#     so downstream modules only receive
#     clean annual yield curves.
# =========================================================

parseON <- function(sim,
                    maxAge = 200) {
  
  library(data.table)
  
  # =======================================================
  # project data directory
  # =======================================================
  
  dPath <- dataPath(sim)
  
  message("Preparing Ontario yield tables...")
  
  # =======================================================
  # curves selected by the classifier
  #
  # sim$AUtoCurve links analysis units (AU)
  # to the Ontario yield curve numbers
  # chosen during the classifier step.
  #
  # Only these curves will be extracted
  # from the raw Ontario YTF files.
  # =======================================================
  
  curves_needed <- unique(
    sim$AUtoCurve$curveID
  )
  
  # =======================================================
  # remove missing curve IDs
  # =======================================================
  
  curves_needed <- curves_needed[
    !is.na(curves_needed)
  ]
  
  # =======================================================
  # Ontario yield table files
  #
  # These are the raw provincial yield tables
  # containing species-specific volumes
  # at irregular age intervals.
  # =======================================================
  
  ytf_files <- list.files(
    file.path(dPath, "ON/YTF"),
    pattern = "tbl_yield_final",
    full.names = TRUE
  )
  
  # =======================================================
  # species volume columns expected in
  # Ontario yield table files
  #
  # These columns are used only to build
  # total stand volume curves.
  #
  # AAC/Hanzlik currently uses total volume,
  # not species-specific yield trajectories.
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
  #
  # Stores standardized annual yield tables
  # indexed by curveID.
  # =======================================================
  
  yieldTables <- list()
  
  # =======================================================
  # loop over all required curves
  # =======================================================
  
  for (curve_id in curves_needed) {
    
    message("Processing curve: ", curve_id)
    
    # =====================================================
    # tracking flag
    #
    # used to determine whether the curve
    # was successfully located in any YTF file
    # =====================================================
    
    curve_found <- FALSE
    
    # ---------------------------------------------------
    # search across all Ontario YTF files
    #
    # each Ontario region may store curves
    # in different files, so all files are checked
    # until the requested curve is found
    # ---------------------------------------------------
    
    for (f in ytf_files) {
      
      message(
        "Checking file: ",
        basename(f)
      )
      
      # -------------------------------------------------
      # read Ontario yield table file
      # -------------------------------------------------
      
      yt <- fread(f)
      
      # -------------------------------------------------
      # skip file if the requested curve
      # does not exist in this table
      # -------------------------------------------------
      
      if (!(curve_id %in% yt$CURVENO)) {
        next
      }
      
      # -------------------------------------------------
      # curve successfully found
      # -------------------------------------------------
      
      curve_found <- TRUE
      
      # -------------------------------------------------
      # extract all rows belonging
      # to the selected curve
      # -------------------------------------------------
      
      curve_dt <- yt[
        CURVENO == curve_id
      ]
      
      # -------------------------------------------------
      # sort by age to ensure proper interpolation
      #
      # standardizeYieldCurve assumes ages
      # are ordered correctly
      # -------------------------------------------------
      
      setorder(
        curve_dt,
        AC10
      )
      
      # -------------------------------------------------
      # identify species columns that are
      # actually present in this file
      #
      # not all Ontario files necessarily
      # contain every species column
      # -------------------------------------------------
      
      spp_available <- intersect(
        species_cols,
        names(curve_dt)
      )
      
      # -------------------------------------------------
      # safety check:
      # ensure at least one species
      # volume column exists
      # -------------------------------------------------
      
      if (length(spp_available) == 0) {
        
        warning(
          "No species columns found for curve: ",
          curve_id
        )
        
        next
      }
      
      # -------------------------------------------------
      # calculate total stand volume
      #
      # Ontario YTF files store yields
      # separately by species.
      #
      # Here all species volumes are summed
      # to produce a total stand yield curve
      # required by the AAC workflow.
      # -------------------------------------------------
      
      curve_dt[
        ,
        total_volume := rowSums(
          .SD,
          na.rm = TRUE
        ),
        .SDcols = spp_available
      ]
      
      # -------------------------------------------------
      # remove rows with missing ages
      #
      # interpolation requires valid ages
      # -------------------------------------------------
      
      curve_dt <- curve_dt[
        !is.na(AC10)
      ]
      
      # -------------------------------------------------
      # standardize yield curve
      #
      # converts irregular Ontario age intervals
      # into annual age increments:
      #
      #   age = 1:maxAge
      #
      # using linear interpolation
      # -------------------------------------------------
      
      yt_standard <- standardizeYieldCurve(
        ages = curve_dt$AC10,
        volumes = curve_dt$total_volume,
        maxAge = maxAge
      )
      
      # -------------------------------------------------
      # store standardized yield table
      #
      # output is indexed by curveID
      # -------------------------------------------------
      
      yieldTables[[as.character(curve_id)]] <- list(yt_standard)
      
      # -------------------------------------------------
      # stop searching additional files
      #
      # the requested curve has already
      # been found and processed
      # -------------------------------------------------
      
      break
    }
    
    # ---------------------------------------------------
    # warning if curve could not be found
    # in any Ontario YTF file
    # ---------------------------------------------------
    
    if (!curve_found) {
      
      warning(
        "Curve not found in Ontario YTF files: ",
        curve_id
      )
    }
  }
  
  # =======================================================
  # return standardized Ontario yield tables
  # =======================================================
  
  return(yieldTables)
}