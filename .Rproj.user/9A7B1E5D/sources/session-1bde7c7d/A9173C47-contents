
#' Title
#'
#' @param shp A terra vect object
#' @param res The resolution of the pam
#' @param ... ellipsis
#'
#' @return A data frame with sites by species
#' @export
#'
#' @examples
#' get_geo_pam(get_rodentia_mexico())
get_geo_pam <- function(shp, res = 0.5, ... ){
  #columns <- rlang::enquos(...)
  r <- terra::rast(shp, res = res)
  z <- terra::rasterize(shp, r, "sciname")
  z_df <- terra::as.data.frame(z, xy = TRUE)

  z_df %>%
    dplyr::mutate(value = 1) %>%
    tidyr::spread(key = 3, value = 4, fill = 0) %>%
    tibble::as_tibble()
}
