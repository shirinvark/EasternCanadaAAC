rm(list = ls())
gc()

# =========================================================
# LOAD LIBRARIES
# =========================================================

library(SpaDES.core)
library(SpaDES.project)
library(terra)

# =========================================================
# SET PATHS
# =========================================================

setPaths(
  cachePath   = "E:/EasternCanada/cache",
  inputPath   = "E:/EasternCanada/inputs",
  outputPath  = "E:/EasternCanada/outputs",
  modulePath  = "E:/EasternCanada/modules",
  scratchPath = "E:/EasternCanada/scratch"
)

print(getPaths())

# =========================================================
# DOWNLOAD MODULES
# =========================================================

SpaDES.project::getModule(
  modules = c(
    "shirinvark/EasternCanadaClassifier",
    "shirinvark/EasternCanadaAAC"
  ),
  modulePath = getPaths()$modulePath,
  overwrite = TRUE
)

# =========================================================
# INITIALIZE SIMULATION (PIPELINE)
# =========================================================

sim <- simInit(
  times = list(start = 0, end = 1),
  
  modules = c(
    "EasternCanadaClassifier",
    "EasternCanadaAAC"
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

# 🔹 نقشه analysis unit (از classifier)
plot(sim$analysisUnitMap)

# 🔹 مقدار AAC
print(sim$AAC)

# 🔹 debug
print(unique(terra::values(sim$analysisUnitMap)))
print(names(sim$yieldTables))