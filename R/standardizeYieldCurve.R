# =========================================================
# STANDARDIZE YIELD CURVE
# =========================================================

standardizeYieldCurve <- function(
    ages,
    volumes,
    maxAge = 200
){
  
  annual <- approx(
    x = ages,
    y = volumes,
    xout = 1:maxAge,
    method = "linear",
    rule = 2
  )$y
  
  annual[annual < 0] <- 0
  
  return(
    data.table(
      age = 1:maxAge,
      volume = annual
    )
  )
}