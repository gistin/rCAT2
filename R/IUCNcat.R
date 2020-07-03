###########################################################
#calculates the IUCN rating based on AOO and EOO area
###########################################################
#' Calculates IUCN rating from EOO
#' @title IUCN rating based from EOO Area
#' @description 
#' Calculates IUCN rating based on Extent of Occurance (EOO) Area in km2
#' @author Justin Moat. J.Moat@kew.org
#' @param EOOArea Area in km2
#' @param abb TRUE or FALSE , TRUE = 2 letter code, FALSE = full text (see value), default = TRUE
#' @return Text
#' one of CR, EN, VU, NT, LC or Critically Endangered, Endangered, Vulnerable, Near Threatened, Least Concern
#' @note Any negative values are assumpted to be positive. Near Threatened is set at 30,000 km2, follow example in IUCN petition 2014
#' @examples
#' EOOArea <- 25 
#' EOORtext <- EOORating(EOOArea,TRUE)
#' @seealso \code{\link{EOOarea}} for EOO calculations
#' @export
#' @references
#' Bachman, S., Moat, J., Hill, A.W., de Torre, J., Scott, B., 2011. Supporting Red List threat assessments with GeoCAT: geospatial conservation assessment tool. Zookeys 126, 117–26. doi:10.3897/zookeys.150.2109 
#' 
#' Moat, J., Bachman, S., n.d. GeoCAT Geospatial Conservation Assessment Tool [WWW Document]. URL http://geocat.kew.org/
#' 
#' IUCN, 2012. IUCN RED LIST CATEGORIES AND CRITERIA, 2nd ed. IUCN, Gland, Switzerland. 
#' 
#' IUCN Standards and Petitions Subcommittee, 2014. Guidelines for Using the IUCN Red List Categories and Criteria. Version 11.
#' 
#' IUCN Standards and Petitions Subcommittee, 2016. Guidelines for Using the IUCN Red List Categories and Criteria. Version 12.
#' 
#' Joppa, L.N., Butchart, S.H.M., Hoffmann, M., Bachman, S.P., Akçakaya, H.R., Moat, J.F., Böhm, M., Holland, R.A., Newton, A., Polidoro, B., Hughes, A., 2016. Impact of alternative metrics on estimates of extent of occurrence for extinction risk assessment. Conserv. Biol. 30, 362–370. doi:10.1111/cobi.12591

EOORating <- function(EOOArea,abb){
  if(missing(abb)){
    abb = TRUE
  }
#  EOOArea <- 250
#  abb <- FALSE
  #make positive
  EOOArea <- sqrt(EOOArea * EOOArea)
  cat <- NA
  if (identical(abb,FALSE)){
    if (EOOArea < 100){
      cat <- "Critically Endangered"
    } else if (EOOArea < 5000){
      cat <- "Endangered"
    } else if (EOOArea < 20000){
      cat <- "Vulnerable"
    } else if (EOOArea < 30000){
      cat <- "Near Threatened"
    } else
      cat <- "Least Concern"
    
  } else {
  if (EOOArea < 100){
    cat <- "CR"
  } else if (EOOArea < 5000){
    cat <- "EN"
  } else if (EOOArea < 20000){
    cat <- "VU"
  } else if (EOOArea < 30000){
    cat <- "NT"
  } else
    cat <- "LC"
  }
  return (cat)
}

#' @title IUCN rating based from AOO Area
#' @description 
#' Calculates IUCN rating based on Area of occupancy (AOO) Area in km2
#' @author Justin Moat. J.Moat@kew.org
#' @param AOOArea Area in km2
#' @param abb TRUE or FALSE , TRUE = 2 letter code, FALSE = full text (see value), default = TRUE
#' @return Text one of CR, EN, VU, NT, LC or Critically Endangered, Endangered, Vulnerable, Near Threatened, Least Concern
#' @note Any negative values are assumpted to be positive. 
#' Near Threatened is set at 3,000 km2, follow example in IUCN Guidelines version 11. 2014
#' @examples
#' AOOArea <- 25 
#' AOORtext <- EOORating(AOOArea,FALSE)
#' @seealso  \code{\link{AOOsimp}} for AOO calculations
#' @export
#' @references
#' Bachman, S., Moat, J., Hill, A.W., de Torre, J., Scott, B., 2011. Supporting Red List threat assessments with GeoCAT: geospatial conservation assessment tool. Zookeys 126, 117–26. doi:10.3897/zookeys.150.2109 
#' 
#' Moat, J., Bachman, S., n.d. GeoCAT Geospatial Conservation Assessment Tool [WWW Document]. URL http://geocat.kew.org/
#' 
#' IUCN, 2012. IUCN Red List Categories and Criteria, 2nd ed. IUCN, Gland, Switzerland. 
#' 
#' IUCN Standards and Petitions Subcommittee, 2014. Guidelines for Using the IUCN Red List Categories and Criteria. Version 11.
#' 
#' IUCN Standards and Petitions Subcommittee, 2016. Guidelines for Using the IUCN Red List Categories and Criteria. Version 12.
#' 
#' Joppa, L.N., Butchart, S.H.M., Hoffmann, M., Bachman, S.P., Akçakaya, H.R., Moat, J.F., Böhm, M., Holland, R.A., Newton, A., Polidoro, B., Hughes, A., 2016. Impact of alternative metrics on estimates of extent of occurrence for extinction risk assessment. Conserv. Biol. 30, 362–370. doi:10.1111/cobi.12591

AOORating <- function(AOOArea,abb){
  if(missing(abb)){
    abb = TRUE
  }
  cat <- NA
  cat <- EOORating(AOOArea*10,abb)
}
