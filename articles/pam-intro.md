# Getting Started with pam

## Overview

The **pam** package builds geographic presence-absence matrices (PAMs)
from species range polygons. A PAM is a sites-by-species binary table
where rows are grid cells and columns are species; a cell contains `1`
when the species range overlaps that cell and `0` otherwise.

PAMs are a fundamental data structure in macroecology and biogeography,
used for computing species richness, range-diversity plots,
co-occurrence analyses, and input to conservation-prioritisation
algorithms.

## Installation

``` r

# Install from GitHub
# install.packages("remotes")
remotes::install_github("alrobles/PAM")
```

## Quick Start

The core workflow has two steps:

1.  Load species range polygons as a
    [`terra::SpatVector`](https://rspatial.github.io/terra/reference/SpatVector-class.html).
2.  Call
    [`get_geo_pam()`](https://alrobles.github.io/pam/reference/get_geo_pam.md)
    with a desired grid resolution.

``` r

library(pam)

# Step 1: load example data (Rodentia of Mexico)
shp <- get_rodentia_mexico()
shp

# Step 2: build the PAM at 1-degree resolution
pam_result <- get_geo_pam(shp, res = 1)
head(pam_result)
```

The result is a tibble. The first two columns (`x`, `y`) are cell
centroids; the remaining columns are species, coded as `1`/`0`.

## Adjusting Resolution

The `res` argument controls the grid cell size in the same units as the
CRS (decimal degrees for EPSG:4326): Higher resolution → more grid cells
→ finer spatial detail.

``` r

# Fine resolution (0.5 degree ~ 55 km at the equator)
pam_fine <- get_geo_pam(shp, res = 0.5)
dim(pam_fine)

# Coarse resolution (2 degrees ~ 220 km)
pam_coarse <- get_geo_pam(shp, res = 2)
dim(pam_coarse)
```

## Computing Species Richness

Once you have a PAM, species richness per cell is just a row sum:

``` r

pam_result$richness <- rowSums(pam_result[, -(1:2)])

# Quick map with base R
plot(pam_result$x, pam_result$y,
     cex  = sqrt(pam_result$richness) * 0.3,
     pch  = 15,
     col  = hcl.colors(max(pam_result$richness))[pam_result$richness],
     xlab = "Longitude",
     ylab = "Latitude",
     main = "Rodentia Species Richness — Mexico")
```

## Environmental Space

The package also provides a WorldClim PCA raster for environmental
analyses:

``` r

world_pca <- get_world_pca()
terra::plot(world_pca)
```

## Using Your Own Data

You can pass any
[`terra::SpatVector`](https://rspatial.github.io/terra/reference/SpatVector-class.html)
that has a `sciname` column: each unique value of `sciname` becomes a
column in the output PAM.

``` r

library(terra)

my_shapes <- vect("path/to/my_ranges.shp")
names(my_shapes)
# [1] "sciname" "family" "order"  ...

my_pam <- get_geo_pam(my_shapes, res = 0.25)
```
