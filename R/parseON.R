# =========================================================
# ONTARIO PARSER
# =========================================================

parseON <- function(sim,
                    maxAge = 200) {
  
  library(data.table)
  
  dPath <- dataPath(sim)
  
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
    file.path(dPath, "ON/YTF"),
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
      
      message(
        "Checking file: ",
        basename(f)
      )
      
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
      
      setorder(
        curve_dt,
        AC10
      )
      
      # -------------------------------------------------
      # available species columns only
      # -------------------------------------------------
      
      spp_available <- intersect(
        species_cols,
        names(curve_dt)
      )
      
      # -------------------------------------------------
      # check species columns
      # -------------------------------------------------
      
      if (length(spp_available) == 0) {
        
        warning(
          "No species columns found for curve: ",
          curve_id
        )
        
        next
      }
      
      # -------------------------------------------------
      # total volume
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
      # remove NA ages
      # -------------------------------------------------
      
      curve_dt <- curve_dt[
        !is.na(AC10)
      ]
      
      # -------------------------------------------------
      # annual interpolation
      # -------------------------------------------------
      
      yt_standard <- standardizeYieldCurve(
        ages = curve_dt$AC10,
        volumes = curve_dt$total_volume,
        maxAge = maxAge
      )
      
      # -------------------------------------------------
      # save standardized curve
      # -------------------------------------------------
      
      yieldTables[
        [as.character(curve_id)]
      ] <- list(yt_standard)
      
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