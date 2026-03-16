Init <- function(sim) {
  
  message("AAC module started")
  
  if (!is.null(sim$areaByAU)) {
    
    message("Area table exists")
    print(sim$areaByAU)
    
  } else {
    
    message("No area table yet")
    
  }
  
  return(invisible(sim))
}