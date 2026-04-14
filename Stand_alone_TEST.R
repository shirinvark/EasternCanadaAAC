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
# DOWNLOAD MODULE
# =========================================================

SpaDES.project::getModule(
  modules    = "shirinvark/EasternCanadaAAC",
  modulePath = getPaths()$modulePath,
  overwrite  = TRUE
)

# =========================================================
# CREATE TEST DATA (MINIMAL + CORRECT)
# =========================================================

# 🔹 analysisUnitMap (3 AU فقط برای سادگی)
analysisUnitMap <- rast(
  nrows = 10,
  ncols = 10,
  xmin  = 0,
  xmax  = 1000,
  ymin  = 0,
  ymax  = 1000
)

values(analysisUnitMap) <- sample(1:3, 100, replace = TRUE)

# 🔹 standAgeMap
standAgeMap <- rast(
  nrows = 10,
  ncols = 10,
  xmin  = 0,
  xmax  = 1000,
  ymin  = 0,
  ymax  = 1000
)

values(standAgeMap) <- sample(1:80, 100, replace = TRUE)

# 🔹 yieldTables (مهم‌ترین بخش)
# هر AU یک curve
# باید list باشه + اسم‌ها match AU

yieldTables <- list(
  "1" = cumsum(rep(2, 100)),
  "2" = cumsum(rep(3, 100)),
  "3" = cumsum(rep(4, 100))
)

# =========================================================
# SANITY CHECK (خیلی مهم قبل اجرا)
# =========================================================

print(unique(values(analysisUnitMap)))
print(names(yieldTables))

# =========================================================
# INITIALIZE SIMULATION
# =========================================================

sim <- simInit(
  times   = list(start = 0, end = 1),
  modules = "EasternCanadaAAC",
  
  objects = list(
    analysisUnitMap = analysisUnitMap,
    standAgeMap     = standAgeMap,
    yieldTables     = yieldTables
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
print(sim$AAC)

# optional debug
print(sim$hanzlikPars)