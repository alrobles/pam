# Build a Presence-Absence Matrix from Species Range Polygons

Rasterizes a set of species range polygons onto a regular grid and
returns a sites-by-species presence-absence matrix (PAM) as a tibble.

## Usage

``` r
get_geo_pam(shp, res = 0.5, ...)
```

## Arguments

- shp:

  A
  [`SpatVector`](https://rspatial.github.io/terra/reference/SpatVector-class.html)
  with a `sciname` attribute column identifying each species.

- res:

  Numeric. Grid cell size in the same units as the coordinate reference
  system of `shp` (typically decimal degrees). Default `0.5`.

- ...:

  Additional arguments passed to
  [`rasterize`](https://rspatial.github.io/terra/reference/rasterize.html).

## Value

A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html) with
columns `x`, `y` (cell centroids) followed by one column per species
containing `1` (present) or `0` (absent).

## Examples

``` r
# \donttest{
shp <- get_rodentia_mexico()
pam <- get_geo_pam(shp, res = 1)
head(pam)
#> # A tibble: 6 × 27
#>       x     y `Orthogeomys grandis` `Peromyscus gambelii` Reithrodontomys fulv…¹
#>   <dbl> <dbl>                 <dbl>                 <dbl>                  <dbl>
#> 1 -117.  32.0                     0                     0                      0
#> 2 -116.  30.0                     0                     0                      0
#> 3 -116.  31.0                     0                     0                      0
#> 4 -116.  32.0                     0                     0                      0
#> 5 -115.  30.0                     0                     0                      0
#> 6 -115.  32.0                     0                     0                      0
#> # ℹ abbreviated name: ¹​`Reithrodontomys fulvescens`
#> # ℹ 22 more variables: `Reithrodontomys megalotis` <dbl>,
#> #   `Sciurus alleni` <dbl>, `Sciurus niger` <dbl>, `Sigmodon alleni` <dbl>,
#> #   `Sigmodon arizonae` <dbl>, `Sigmodon hirsutus` <dbl>,
#> #   `Sigmodon hispidus` <dbl>, `Sigmodon leucotis` <dbl>,
#> #   `Sigmodon mascotensis` <dbl>, `Sigmodon ochrognathus` <dbl>,
#> #   `Sigmodon planifrons` <dbl>, `Sigmodon toltecus` <dbl>, …
# }
```
