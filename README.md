# pam

<!-- badges: start -->
[![R-CMD-check](https://github.com/alrobles/PAM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/alrobles/PAM/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/alrobles/PAM/actions/workflows/pkgdown.yaml/badge.svg)](https://alrobles.github.io/pam/)
<!-- badges: end -->

**pam** builds geographic presence-absence matrices (PAMs) from species
range polygons stored as `terra::SpatVector` objects.

A PAM is a sites-by-species binary table: rows are grid cells, columns
are species, and each entry is `1` (present) or `0` (absent). PAMs are a
foundational data structure in macroecology and biogeography for
computing species richness, range-diversity plots, co-occurrence
analysis, and conservation prioritisation.

## Installation

```r
# Install from GitHub
# install.packages("remotes")
remotes::install_github("alrobles/PAM")
```

## Quick Start

```r
library(pam)

# Load example data: Rodentia of Mexico
shp <- get_rodentia_mexico()

# Build PAM at 1-degree resolution
pam_result <- get_geo_pam(shp, res = 1)
head(pam_result)
#> # A tibble: 6 × 42
#>       x     y Chaetodipus arenarius ...
#>   <dbl> <dbl>                 <dbl>
#> 1 -118.  23.5                     0
#> ...
```

## Main Functions

| Function | Description |
|---|---|
| `get_geo_pam()` | Build a PAM from a `SpatVector` at a given resolution |
| `get_rodentia_mexico()` | Example data: rodent ranges for Mexico |
| `get_world_pca()` | WorldClim bioclimatic PCA raster (2 layers) |

## Citation

If you use this package, please cite:

```
Robles-Fernandez, A. L. (2023). pam: Presence-absence matrix from shapefiles.
R package. doi:10.5281/zenodo.8189393
```

## License

GPL (>= 3)
