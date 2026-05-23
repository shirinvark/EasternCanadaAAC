parseNL <- function(sim,
                    maxAge = 200) {
  

  message("Preparing Newfoundland yield tables...")
  
  curves_needed <- unique(
    sim$AUtoCurve$curveID
  )
  
  curves_needed <- curves_needed[
    !is.na(curves_needed)
  ]
  
  raw_curves <- sim$rawYieldTables$NL
  
  yieldTables <- list()
  
  for (curve_name in curves_needed) {
    
    message("Processing: ", curve_name)
    
    curve_found <- FALSE
    
    for (region_name in names(raw_curves)) {
      
      region_curves <- raw_curves[[region_name]]
      
      if (!(curve_name %in% names(region_curves))) {
        next
      }
      
      curve_found <- TRUE
      
      species_vectors <- region_curves[[curve_name]]
      
      if (length(species_vectors) == 0) {
        
        warning(
          "No species parsed from: ",
          curve_name
        )
        
        next
      }
      
      min_len <- min(
        sapply(
          species_vectors,
          length
        )
      )
      
      species_vectors <- lapply(
        species_vectors,
        function(x) {
          x[1:min_len]
        }
      )
      
      total_volume <- Reduce(
        "+",
        species_vectors
      )
      
      ages <- seq(
        0,
        by = 5,
        length.out = length(total_volume)
      )
      
      yt_standard <- standardizeYieldCurve(
        ages = ages,
        volumes = total_volume,
        maxAge = maxAge
      )
      
      yieldTables[[curve_name]] <- yt_standard
      
      break
    }
    
    if (!curve_found) {
      
      warning(
        "Curve not found in Newfoundland yield tables: ",
        curve_name
      )
    }
  }
  
  return(yieldTables)
}