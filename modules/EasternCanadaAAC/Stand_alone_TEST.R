rm(list = ls())
gc()

library(SpaDES.core)
library(SpaDES.project)
library(terra)

setPaths(
  cachePath   = "E:/EasternCanadaAAC/cache",
  inputPath   = "E:/EasternCanadaAAC/inputs",
  outputPath  = "E:/EasternCanadaAAC/outputs",
  modulePath  = "E:/EasternCanadaAAC/modules",
  scratchPath = "E:/EasternCanadaAAC/scratch"
)

SpaDES.project::getModule(
  modules    = "shirinvark/EasternCanadaAAC",
  modulePath = getPaths()$modulePath,
  overwrite  = TRUE
)

# -----------------------------
# 1️⃣ Build sim
# -----------------------------
sim <- simInit(
  times   = list(start = 0, end = 1),
  modules = "EasternCanadaAAC",
  
  params = list(
    EasternCanadaAAC = list(
      #jurisdiction = "NL",
      rotationPeriodShift = 0
    )
  ),
  
  options = list(
    spades.checkpoint = FALSE,
    spades.save       = FALSE,
    spades.progress   = FALSE
  )
)

# -----------------------------
# 2️⃣ Force young forest
# -----------------------------
terra::values(sim$standAgeMap) <- 5

# -----------------------------
# 3️⃣ Run
# -----------------------------
sim <- spades(sim)

# -----------------------------
# 4️⃣ AAC result
# -----------------------------
sim$AAC
