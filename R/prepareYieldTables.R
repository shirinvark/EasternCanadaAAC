prepareYieldTables <- function(sim,
                               maxAge = 255) {
  
  # -------------------------------------------------------
  # Jurisdiction parameter
  # -------------------------------------------------------
  
  jur <- P(sim)$jurisdiction
  
  # -------------------------------------------------------
  # Standalone mode:
  # if no jurisdiction supplied,
  # keep existing/fake yieldTables
  # -------------------------------------------------------
  
  if (is.null(jur) || length(jur) == 0) {
    
    message(
      "No jurisdiction supplied → using existing/fake yieldTables"
    )
    
    return(NULL)
  }
  
  jur <- toupper(jur)
  
  # -------------------------------------------------------
  # Newfoundland
  # -------------------------------------------------------
  
  if (jur == "NL") {
    return(parseNL(sim, maxAge))
  }
  
  # -------------------------------------------------------
  # Ontario
  # -------------------------------------------------------
  
  if (jur == "ON") {
    return(parseON(sim, maxAge))
  }
  
  # -------------------------------------------------------
  # Unsupported jurisdiction
  # -------------------------------------------------------
  
  stop(
    "Unsupported jurisdiction: ",
    jur
  )
}