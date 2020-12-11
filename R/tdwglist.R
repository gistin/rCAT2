#' Extract wgsrpd regions from a set of points in lat long
#' @param thepoints dataframe of points ie c(x,y) in lat long
#' @return dataframe or sf object with wgsrpd regions
#' @examples
#' #built set of points
#' thepoints <- ptsNormal(29,0.1)
#'#shift to Peru, Colombia and Brazil tripoint
#' thepoints <- data.frame(long = thepoints$X - 69.9465972, lat = thepoints$Y - 4.2260167)
#' thepoints <- sf::st_as_sf(thepoints, coords = c("long", "lat"), crs = 4326)
#' tdwglist(thepoints)
#' @import sf
#' @import dplyr
#' @export

tdwglist <- function(thepointsll, returnV = 'S') {
  
  # read in the wgsrpd level 3 polygon data
  wgsrpd_l3 <- readRDS("data/wgsrpd_l3.rds")
  
  # intersect the points with the wgsrpd level 3 polygons
  l3_points <- sf::st_intersection(thepointsll, wgsrpd_l3)
  
  # df <- st_drop_geometry(df)
  l3_points_df <- sf::st_drop_geometry(l3_points)
  
  # remove duplicates from df
  l3_unique <- unique(l3_points_df)
  
  if (returnV == "SF") {
    
    # merge the table of unique with tdwg regions with sf
    l3_unique_sf <- merge(wgsrpd_l3, l3_unique)
    
    return(l3_unique_sf)
  }
  
  return(l3_unique)
}
