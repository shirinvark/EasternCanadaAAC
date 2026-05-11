

# =========================================================
# ONTARIO PARSER
# =========================================================

parseON <- function(sim,
                    maxAge = 200) {
  
  library(data.table)
  
  message("Preparing Ontario yield tables...")
  
  # =======================================================
  # curves selected by classifier
  # =======================================================
  
  curves_needed <- c(1, 2, 3, 4)
  
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
      
      yt_standard <- standardizeYieldCurve(
        ages = curve_dt$AC10,
        volumes = curve_dt$total_volume,
        maxAge = maxAge
      )
      # -------------------------------------------------
      # save standardized curve
      # -------------------------------------------------
      
      yieldTables[[as.character(curve_id)]] <- yt_standard
      
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
