
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
    
    yt_standard <- standardizeYieldCurve(
      ages = ages,
      volumes = total_volume,
      maxAge = maxAge
    )
    
    # ---------------------------------------------------
    # save
    # ---------------------------------------------------
    
    yieldTables[[curve_name]] <- yt_standard
  }
  
  return(yieldTables)
}

