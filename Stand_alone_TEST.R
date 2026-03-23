rm(list = ls())
gc()

# =========================================================
# LOAD LIBRARIES
# =========================================================

library(SpaDES.core)
library(SpaDES.project)
library(terra)
library(data.table)

# =========================================================
# SET PATHS
# =========================================================

setPaths(
  cachePath   = "E:/EasternCanadaAAC/cache",
  inputPath   = "E:/EasternCanadaAAC/inputs",
  outputPath  = "E:/EasternCanadaAAC/outputs",
  modulePath  = "E:/EasternCanadaAAC/modules",
  scratchPath = "E:/EasternCanadaAAC/scratch"
)

print(getPaths())

# =========================================================
# DOWNLOAD AAC MODULE
# =========================================================

SpaDES.project::getModule(
  modules    = "shirinvark/EasternCanadaAAC",
  modulePath = getPaths()$modulePath,
  overwrite  = TRUE
)

# =========================================================
# CREATE TEST DATA
# =========================================================

# Fake analysisUnitMap

analysisUnitMap <- rast(
  nrows = 10,
  ncols = 10,
  xmin  = 0,
  xmax  = 1000,
  ymin  = 0,
  ymax  = 1000
)

values(analysisUnitMap) <- sample(1:5, 100, replace = TRUE)

# Fake area table

areaByAU <- data.table(
  analysisUnit = 1:5,
  area_ha = runif(5, 1000, 5000)
)

# Fake age summary

ageSummaryByAU <- data.table(
  analysisUnit = 1:5,
  meanAge = sample(20:120, 5),
  nStands = sample(50:200, 5)
)

# Fake yield curves

yieldAges <- seq(0,120,10)

yieldTables <- matrix(
  runif(13*5, 0, 200),
  nrow = 5
)

# =========================================================
# INITIALIZE SIMULATION
# =========================================================

sim <- simInit(
  times   = list(start = 0, end = 1),
  modules = "EasternCanadaAAC",
  
  objects = list(
    analysisUnitMap = analysisUnitMap,
    areaByAU = areaByAU,
    ageSummaryByAU = ageSummaryByAU,
    yieldTables = yieldTables,
    yieldAges = yieldAges
  ),
  
  options = list(
    spades.checkpoint = FALSE,
    spades.save       = FALSE,
    spades.progress   = FALSE
  )
)

# =========================================================
# RUN MODEL
# =========================================================

system.time({
  sim <- spades(sim)
})

# =========================================================
# CHECK OUTPUTS
# =========================================================

print(names(sim))

