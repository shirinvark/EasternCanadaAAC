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