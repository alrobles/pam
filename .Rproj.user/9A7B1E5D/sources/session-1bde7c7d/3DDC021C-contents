#' Function to return mammal shapefiles of afrosoricida Order
#' @importFrom utils download.file
#' @return A SpatVector object with mammal shapefiles
#' @export
#' @source A phylogeny-informed characterization of global tetrapod traits
#'  addresses data gaps and biases
#'  \url{https://doi.org/10.1101/2023.03.04.531098}
#'
#'
#' @examples
#' get_rodentia_mexico()
get_rodentia_mexico <- function(){
 url <- "https://phymaps.nyc3.digitaloceanspaces.com/rodentia_mexico.rds"

    # check if the url exists
    if (!url_exists(url)) {
      return("Can't access data.")
    }
  shp <- readr::read_rds(url)
  return(terra::unwrap(shp))
}
