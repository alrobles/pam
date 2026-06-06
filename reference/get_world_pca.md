# Retrieve a WorldClim PCA Environmental Raster

Downloads a
[`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
containing the first two principal components derived from WorldClim 2.0
bioclimatic variables at global extent.

## Usage

``` r
get_world_pca()
```

## Source

WorldClim 2.0 project: <https://www.worldclim.org/>

## Value

A
[`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
object, or `"Can't access data."` if the remote server is unreachable.

## Examples

``` r
# \donttest{
ras <- get_world_pca()
# }
```
