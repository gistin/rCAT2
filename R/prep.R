#data preparations scripts to help upload to IUCN etc
#note these will need work to get them more useful

##########################################################
#Country list from a set of points in geographic coordinates#
##########################################################
#' calculates the MER of a set of numbers'
#' @title PRE-ALPHA Country list from a set of points in geographic coordinates
#' @description 
#' Calculates the country Occurances from a set of points in lat long
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
#' @export



countrylist <- function(thepointsll){
  #check to see if it's been downloaded
  if (!exists('countries10')){
    countries10 <- new.env()
    countries10 <- ne_download_rcat()
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


#downloads rnatural earth dataset for countries
#borrowed but much reduced from https://github.com/ropensci/rnaturalearth/blob/master/R/ne_download.R

ne_download_rcat <- function() {
  
  file_name <- "ne_10m_admin_0_countries"
  address <- "https://naturalearth.s3.amazonaws.com/10m_cultural/ne_10m_admin_0_countries.zip"
  # download zip to temporary location, unzipped files are saved later
  # tryCatch catches error, returns NUll if no error
  download_failed <- tryCatch(
    utils::download.file(file.path(address), zip_file <- tempfile()),
    error = function(e) {
      message(paste("Natural earth download failed"))
      return(TRUE)
    }
  )
  # return from this function if download error was caught by tryCatch
  if (download_failed) {
    return()
  }
  
  #unzip it
  utils::unzip(zip_file, exdir = tempdir())
  #read as sf
  sf_object <- read_sf(tempdir(), file_name)
  library (spdep)
  sf_object <- poly2nb(st_make_valid(sf_object))
  return(sf_object)
}




