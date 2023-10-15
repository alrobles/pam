#' Function to obtain the environmental space of the world projected
#' to two dimensions through Principal Component Analysis
#' @importFrom utils download.file
#' @importFrom readr read_rds
#' @importFrom terra unwrap
#' @return A SpatRaster object with World PCA environmental space
#' @export
#' @source WorldClim 2.0 project: \url{https://www.worldclim.org/}
#'
#' @examples
#' get_world_pca()
get_world_pca <- function(){
  url <- "https://phymaps.nyc3.digitaloceanspaces.com/worldPCA.rds"

  # check if the url exists
  if (!url_exists(url)) {
    return("Can't access data.")
  }
  ras <- readr::read_rds(url)
  return(terra::unwrap(ras))
}
