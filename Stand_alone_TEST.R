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
# 1️⃣ اول sim رو بساز
# -----------------------------
sim <- simInit(
  times   = list(start = 0, end = 1),
  modules = "EasternCanadaAAC",
  options = list(
    spades.checkpoint = FALSE,
    spades.save       = FALSE,
    spades.progress   = FALSE
  )
)

# -----------------------------
# 2️⃣ بعد تغییر بده
# -----------------------------
terra::values(sim$standAgeMap) <- 5
terra::values(sim$standAgeMap) <- 80

# -----------------------------
# 3️⃣ اجرا
# -----------------------------
sim <- spades(sim)

sim$AAC
