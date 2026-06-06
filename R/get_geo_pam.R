#' Build a Presence-Absence Matrix from Species Range Polygons
#'
#' Rasterizes a set of species range polygons onto a regular grid and
#' returns a sites-by-species presence-absence matrix (PAM) as a tibble.
#'
#' @param shp A \code{\link[terra]{SpatVector}} with a \code{sciname}
#'   attribute column identifying each species.
#' @param res Numeric.
#'   Grid cell size in the same units as the coordinate reference system
#'   of \code{shp} (typically decimal degrees). Default \code{0.5}.
#' @param ... Additional arguments passed to \code{\link[terra]{rasterize}}.
#'
#' @return A \code{\link[tibble]{tibble}} with columns \code{x}, \code{y}
#'   (cell centroids) followed by one column per species containing
#'   \code{1} (present) or \code{0} (absent).
#'
#' @export
#'
#' @examples
#' \donttest{
#' shp <- get_rodentia_mexico()
#' pam <- get_geo_pam(shp, res = 1)
#' head(pam)
#' }
get_geo_pam <- function(shp, res = 0.5, ...) {
  checkmate::assert_class(shp, "SpatVector")
  checkmate::assert_number(res, lower = 0, finite = TRUE)

  r <- terra::rast(shp, res = res)
  z <- terra::rasterize(shp, r, "sciname")
  z_df <- terra::as.data.frame(z, xy = TRUE)

  z_df |>
    dplyr::mutate(value = 1) |>
    tidyr::spread(key = 3, value = 4, fill = 0) |>
    tibble::as_tibble()
}
