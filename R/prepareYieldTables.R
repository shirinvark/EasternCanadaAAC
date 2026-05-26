prepareYieldTables <- function(sim,
                               maxAge = 200) {
  
  if (is.null(sim$yieldTables)) {
    
    sim$yieldTables <- standardizeYieldTables(
      sim = sim,
      maxAge = maxAge
    )
  }
  
  return(sim$yieldTables)
}