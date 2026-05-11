---
title: "EasternCanadaAAC Manual"
subtitle: "v.0.0.0.9000"
date: "Last updated: 2026-05-04"
output:
  bookdown::html_document2:
    toc: true
    toc_float: true
    theme: sandstone
    number_sections: true
    df_print: paged
    keep_md: yes
editor_options:
  chunk_output_type: console
link-citations: true
always_allow_html: true
---

# EasternCanadaAAC Module

(ref:EasternCanadaAAC) *EasternCanadaAAC*



---

## Authors


``` r
paste(as.character(SpaDES.core::moduleMetadata(
  module = "EasternCanadaAAC",
  path = ".."
)$authors), collapse = ", ")
```

```
## [1] "First Middle Last <email@example.com> [aut, cre]"
```

---

# 1. Overview

## Summary

The **EasternCanadaAAC** module calculates Annual Allowable Cut (AAC) using a simplified implementation of the Hanzlik method.

The module translates forest structure (age, yield, and area) into a sustainable harvest level.

---

## What the module does

* Computes rotation age from yield tables (based on MAI)
* Separates forest into:

  * Mature stands (harvestable)
  * Immature stands (growing)
* Integrates effective harvestable area
* Computes AAC per Analysis Unit (AU)
* Aggregates total AAC

---

## Key idea

AAC is calculated as:

**AAC = (Mature Volume / Rotation Period) + Growth**

Where:

* Mature volume comes from stands older than rotation age
* Growth comes from younger stands

---

# 2. Inputs


``` r
df_inputs <- SpaDES.core::moduleInputs("EasternCanadaAAC", path = "..")
knitr::kable(df_inputs, caption = "Module input objects.")
```



Table: (\#tab:unnamed-chunk-2)Module input objects.

|objectName      |objectClass |desc                            |sourceURL |
|:---------------|:-----------|:-------------------------------|:---------|
|analysisUnitMap |SpatRaster  |Map of analysis units per pixel |NA        |
|pixelGroupMap   |SpatRaster  |Pixel groups                    |NA        |
|pixelAreaDT     |data.table  |Effective area                  |NA        |
|standAgeMap     |SpatRaster  |Stand age per pixel             |NA        |
|yieldTables     |list        |Yield tables per AU             |NA        |

### Description of inputs

* **analysisUnitMap** → assigns each pixel to an AU
* **standAgeMap** → stand age per pixel
* **pixelGroupMap** → unique pixel identifier
* **pixelAreaDT** → effective (harvestable) area per pixel
* **yieldTables** → yield curves per AU

---

# 3. Parameters


``` r
df_params <- SpaDES.core::moduleParams("EasternCanadaAAC", path = "..")
knitr::kable(df_params, caption = "Module parameters.")
```



Table: (\#tab:unnamed-chunk-3)Module parameters.

|paramName                |paramClass |default |min |max |paramDesc                       |
|:------------------------|:----------|:-------|:---|:---|:-------------------------------|
|replanInterval           |numeric    |10      |NA  |NA  |Years between AAC recalculation |
|rotationPeriodMultiplier |numeric    |1       |NA  |NA  |Multiplier for Hanzlik rotation |
|.plots                   |character  |screen  |NA  |NA  |Plotting option                 |

---

# 4. Outputs


``` r
df_outputs <- SpaDES.core::moduleOutputs("EasternCanadaAAC", path = "..")
knitr::kable(df_outputs, caption = "Module outputs.")
```



Table: (\#tab:unnamed-chunk-4)Module outputs.

|objectName  |objectClass |desc                           |
|:-----------|:-----------|:------------------------------|
|AAC         |numeric     |Annual Allowable Cut (m3/year) |
|AAC_by_AU   |data.table  |AAC per AU                     |
|hanzlikPars |list        |Hanzlik parameters             |

---

## Key outputs

* **AAC** → total allowable cut (m³/year)
* **AAC_by_AU** → breakdown per analysis unit
* **hanzlikPars** → rotation age and growth parameters

---

# 5. Methodology

## Step 1: Yield processing

Each yield table is converted into:

* Volume (V)
* Increment (I)
* Rotation age (R)

Rotation age is defined as:

> Age at which Mean Annual Increment (MAI = V(t)/t) is maximized

---

## Step 2: Pixel-level data

Each pixel contributes:

* Analysis Unit (AU)
* Age
* Effective area

---

## Step 3: Separation of stands

* Mature stands: age ≥ R
* Immature stands: age < R

---

## Step 4: AAC calculation

For each AU:

* Mature volume is summed across pixels
* Growth is summed from immature stands

Then:

**AAC = (V_mature / Rotation Period) + I_total**

---

# 6. Example Run


``` r
library(SpaDES.core)

sim <- simInit(
  times   = list(start = 0, end = 1),
  modules = "EasternCanadaAAC",
  paths   = list(modulePath = ".."),
  options = list(
    spades.checkpoint = FALSE,
    spades.save       = FALSE,
    spades.progress   = FALSE
  )
)

sim <- spades(sim)

sim$AAC
```

```
## [1] 365.3845
```

``` r
sim$AAC_by_AU
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["AU"],"name":[1],"type":["int"],"align":["right"]},{"label":["AAC"],"name":[2],"type":["dbl"],"align":["right"]}],"data":[{"1":"1","2":"132.7765"},{"1":"2","2":"123.0728"},{"1":"3","2":"109.5352"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

---

# 7. Notes

* If inputs are missing, synthetic data are generated
* Yield curves can be:

  * Real (.yld files)
  * Synthetic (for testing)
* Effective area integrates:

  * protected areas
  * riparian buffers
  * landbase constraints

---

# 8. Integration

This module is designed to work with:

* Classifier module (providing AU and area)
* Landbase module (constraints)
* Yield preparation module

---

# 9. Interpretation

* Higher AAC → more mature forest or higher growth
* Lower AAC → younger forest or constrained area

---

# 10. References
