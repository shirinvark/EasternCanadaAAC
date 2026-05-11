library(data.table)

yt <- fread(
  "data/ON/YTF/3e_tbl_yield_final.csv"
)
names(yt)
head(yt)
unique(yt$CURVENO)[1:20]
unique(yt$AC10)
curve1 <- yt[CURVENO == 1]
head(curve1)
species_cols <- c(
  "PW","PR","PJ","SB","SW","BF",
  "CE","OC","HE","PO","PB",
  "BW","MH","QR","YB","OH"
)
curve1[, total_volume := rowSums(
  .SD,
  na.rm = TRUE
), .SDcols = species_cols]
plot(
  curve1$AC10,
  curve1$total_volume,
  type = "b"
)
annual <- approx(
  x = curve1$AC10,
  y = curve1$total_volume,
  xout = 1:255,
  rule = 2
)$y
plot(
  1:255,
  annual,
  type = "l"
)
