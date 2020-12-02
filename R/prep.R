#data preparations scripts to help upload to IUCN etc

##########################################################
#Country list from a set of points in geographic coordinates#
##########################################################
#' calculates the MER of a set of numbers'
#' @title PRE-ALPHA Country list from a set of points in geographic coordinates
#' @description 
#' Calculates the country occurances from a set of points in lat long
#' @author Justin Moat. J.Moat@kew.org
#' @param thepoints dataframe of points ie c(x,y) in lat long
#' @return list of ISO3 codes
#' @examples
#' #built set of points
#' thepoints <- ptsNormal(29,0.1)
#'#shift to Peru, Colombia and Brazil tripoint
#' thepoints <- data.frame(long = thepoints$X - 69.9465972, lat = thepoints$Y - 4.2260167)
#' countrylist(thepoints)
#' @import sf
#' @import dplyr
#' @import rnaturalearth
#' @export



countrylist <- function(thepointsll){
  #check to see if it's been downloaded
  if (!exists('countries10')){
    countries10 <- new.env()
    #note below does not always work, ie NE site can be down!!!
    countries10 <- ne_download(scale = 10, type='countries', category = 'cultural', returnclass='sf')
  }
  thepoints <- sf::st_as_sf(thepointsll, coords = c("long", "lat"), crs = 4326)
  clist <- sf::st_intersects(thepoints,countries10)
  iso2l <- unique(countries10$ISO_A2[unlist(clist)])
  # convert to data frame
  iso2_df <- data.frame(ISO2 = iso2l)
  # join on ISO2 to get IUCN codes
  iucn_df <- dplyr::left_join(iso2_df,tdwg_to_iucn, 
                       by = c("ISO2" = "countryoccurrence.countryoccurrencesubfield.countryoccurrencelookup"), 
                       keep = TRUE)
  # reformat
  iucn_df$countryoccurrence.countryoccurrencesubfield.presence = "Extant"
  iucn_df$countryoccurrence.countryoccurrencesubfield.origin = "Native"
  iucn_df$countryoccurrence.countryoccurrencesubfield.seasonality = "Resident"
  #iucn_df$internal_taxon_id = POWO_ID
  
  countries = dplyr::select(
    iucn_df,
    #internal_taxon_id,
    countryoccurrence.countryoccurrencesubfield.countryoccurrencename,
    countryoccurrence.countryoccurrencesubfield.presence,
    countryoccurrence.countryoccurrencesubfield.origin,
    countryoccurrence.countryoccurrencesubfield.seasonality,
    countryoccurrence.countryoccurrencesubfield.countryoccurrencelookup
  )
  
  countries = dplyr::distinct(countries)
  return(countries)
}
