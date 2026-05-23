parseON <- function(sim,
                    maxAge = 200) {
  

  message("Preparing Ontario yield tables...")
  
  curves_needed <- unique(
    sim$AUtoCurve$curveID
  )
  
  curves_needed <- curves_needed[
    !is.na(curves_needed)
  ]
  
  ytf_tables <- sim$rawYieldTables$ON
  
  species_cols <- c(
    "PW","PR","PJ",
    "SB","SW","BF",
    "CE","OC","HE",
    "PO","PB",
    "BW","MH","QR",
    "YB","OH"
  )
  
  yieldTables <- list()
  
  for (curve_id in curves_needed) {
    
    message("Processing curve: ", curve_id)
    
    curve_found <- FALSE
    
    for (region_name in names(ytf_tables)) {
      
      yt <- ytf_tables[[region_name]]
      
      if (!(curve_id %in% yt$CURVENO)) {
        next
      }
      
      curve_found <- TRUE
      
      curve_dt <- yt[
        CURVENO == curve_id
      ]
      
      setorder(
        curve_dt,
        AC10
      )
      
      spp_available <- intersect(
        species_cols,
        names(curve_dt)
      )
      
      if (length(spp_available) == 0) {
        
        warning(
          "No species columns found for curve: ",
          curve_id
        )
        
        next
      }
      
      curve_dt[
        ,
        total_volume := rowSums(
          .SD,
          na.rm = TRUE
        ),
        .SDcols = spp_available
      ]
      
      curve_dt <- curve_dt[
        !is.na(AC10)
      ]
      
      yt_standard <- standardizeYieldCurve(
        ages = curve_dt$AC10,
        volumes = curve_dt$total_volume,
        maxAge = maxAge
      )
      
      yieldTables[[as.character(curve_id)]] <- list(
        yt_standard
      )
      
      break
    }
    
    if (!curve_found) {
      
      warning(
        "Curve not found in Ontario yield tables: ",
        curve_id
      )
    }
  }
  
  return(yieldTables)
}