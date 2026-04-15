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

#yieldTables <- list(
 # "1" = cumsum(seq(0.5, 2, length.out = 100)),
  #"2" = cumsum(seq(1, 3, length.out = 100)),
  #"3" = cumsum(seq(1, 4, length.out = 100))
#)
# =========================================================
# READ REAL YIELD (.yld)
# =========================================================

lines <- readLines("E:/EasternCanadaClassifier/data/NL/YTF/BarNS_sub_all.yld")

# پیدا کردن شروع جدول
start <- grep("Age", lines)[1]

# حذف header و خواندن بدون header
tab <- read.table(
  text = lines[(start+1):length(lines)],
  header = FALSE,
  fill = TRUE
)

# اسم ستون‌ها رو دستی بده (بر اساس NL structure)
colnames(tab) <- c("Age", "BSv", "WSv", "BFv", "WPv", "TLv", "YBv")

# حالا کار می‌کنه 👇
yt <- tab$BSv

# اگر خواستی کل volume:
# yt <- tab$BSv + tab$WSv + tab$BFv

# =========================================================
# INTERPOLATE TO ANNUAL
# =========================================================

yt_interp <- approx(
  x = tab$Age,
  y = yt,
  xout = 1:100
)$y

yt_interp[is.na(yt_interp)] <- 0
yt_interp <- stats::filter(yt_interp, rep(1/3, 3), sides = 2)
yt_interp <- as.numeric(yt_interp)
yt_interp <- pmax(yt_interp, 0)
# =========================================================
# USE IN MODEL
# =========================================================

yieldTables <- list(
  "1" = yt_interp,
  "2" = yt_interp,
  "3" = yt_interp
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

