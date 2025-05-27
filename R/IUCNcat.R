###########################################################
#calculates the IUCN rating based on AOO and EOO area
###########################################################
#' Calculates IUCN rating on EOO
#' @title IUCN rating based on EOO Area
#' @description 
#' Calculates IUCN rating based on Extent of Occurrence (EOO) Area in km2
#' @author Justin Moat. J.Moat@kew.org
#' @param EOOArea Area in km2
#' @param abb abbreviation TRUE or FALSE , TRUE = 2 letter code, FALSE = full text (see value), default = TRUE
#' @return Text
#' one of CR, EN, VU, NT, LC or Critically Endangered, Endangered, Vulnerable, Near Threatened, Least Concern
#' @note Any negative values are assumed to be positive. Near Threatened is set at 30,000 km2, follow example in IUCN petition 2014
#' @examples
#' ratingEoo(250,TRUE)
#' ratingEoo(250,FALSE)
#' @seealso \code{\link{eoo}} for EOO calculations
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

ratingEoo <- function(EOOArea,abb=TRUE){
  #  EOOArea <- 250
  #  abb <- FALSE
  #make positive
  EOOArea <- sqrt(EOOArea * EOOArea)
  cat <- NA
  if (identical(abb,FALSE)){
    if (is.na(EOOArea)) {
      cat <- "Data Deficient"
    } else if (EOOArea < 100){
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
    if (is.na(EOOArea)){
      cat <- "DD"
    } else if (EOOArea < 100){
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

#kept for consistency with old version 0.1.6
ratingEoo <- ratingEoo

#' @title IUCN rating based on AOO Area
#' @description 
#' Calculates IUCN rating based on Area of occupancy (AOO) in km2
#' @author Justin Moat. J.Moat@kew.org
#' @param AOOArea Area in km2
#' @param abb abbreviation TRUE or FALSE , TRUE = 2 letter code, FALSE = full text (default = TRUE)
#' @return Text one of CR, EN, VU, NT, LC or Critically Endangered, Endangered, Vulnerable, Near Threatened, Least Concern
#' @note  
#' Near Threatened is set at 3,000 km2, follow example in IUCN Guidelines version 11. 2014
#' @examples
#' #will return "Endangered"
#' ratingAoo(25,FALSE)
#' #will return "EN"
#' ratingAoo(25,TRUE)
#' @seealso  \code{\link{aoo}} for AOO calculations
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

ratingAoo <- function(AOOArea,abb=TRUE){
  if(missing(abb)){
    abb = TRUE
  }
  cat <- NA
  cat <- ratingEoo(AOOArea*10,abb)
  return(cat)
}

#kept for consistency with old version 0.1.6
ratingAoo <- ratingAoo


###########################################################
#calculates the IUCN rating based on population reduction
###########################################################
#' Calculates IUCN rating based on population reduction
#' @title IUCN rating based on population reduction
#' @description 
#' Calculates IUCN rating based on based on population reduction as a percentage
#' @author Justin Moat. J.Moat@kew.org
#' @param pReduction reduction as a percentage
#' @param subCr sub Criteria category: 1 or 2 or 3 or 4 (2 default)
#' @note subCr 1 use the limits 90%,70%,50%,25% for CR,EN,VU,NT
#' @note subCr 2,3,4, use the limits 80%,50%,30%,10% for CR,EN,VU,NT
#' @note NT category is a guideline 
#' @param abb abbreviation TRUE or FALSE , TRUE = 2 letter code, FALSE = full text (see value), default = TRUE
#' @return Text
#' one of CR, EN, VU, NT, LC or Critically Endangered, Endangered, Vulnerable, Near Threatened, Least Concern
#' @examples
#' ratingPop(65,2,TRUE)
#' ratingPop(65,2,FALSE)
#' @export
#' @references
#' IUCN, 2012. IUCN RED LIST CATEGORIES AND CRITERIA, 2nd ed. IUCN, Gland, Switzerland. 
#' 
#' IUCN Standards and Petitions Subcommittee, 2014. Guidelines for Using the IUCN Red List Categories and Criteria. Version 11.
#' 
#' IUCN Standards and Petitions Subcommittee, 2016. Guidelines for Using the IUCN Red List Categories and Criteria. Version 12.


ratingPop <- function(pReduction,subCr=2,abb=TRUE){
  #  EOOArea <- 250
  #  abb <- FALSE
  cat <- NA
  if ( subCr > 1){ #2 or 3 or 4
  if (identical(abb,FALSE)){
    if (pReduction >= 80){
      cat <- "Critically Endangered"
    } else if (pReduction >= 50){
      cat <- "Endangered"
    } else if (pReduction >= 30){
      cat <- "Vulnerable"
    } else if (pReduction >= 10){
      cat <- "Near Threatened"
    } else
      cat <- "Least Concern"
    
  } else {
    if (pReduction >= 80){
      cat <- "CR"
    } else if (pReduction >= 50){
      cat <- "EN"
    } else if (pReduction >= 30){
      cat <- "VU"
    } else if (pReduction >= 10){
      cat <- "NT"
    } else
      cat <- "LC"
  }}
  if ( subCr == 1){ #2 or 3 or 4
    if (identical(abb,FALSE)){
      if (pReduction >= 90){
        cat <- "Critically Endangered"
      } else if (pReduction >= 70){
        cat <- "Endangered"
      } else if (pReduction >= 50){
        cat <- "Vulnerable"
      } else if (pReduction >= 25){
        cat <- "Near Threatened"
      } else
        cat <- "Least Concern"
      
    } else {
      if (pReduction >= 90){
        cat <- "CR"
      } else if (pReduction >= 70){
        cat <- "EN"
      } else if (pReduction >= 50){
        cat <- "VU"
      } else if (pReduction >= 25){
        cat <- "NT"
      } else
        cat <- "LC"
    }}
  return (cat)
}


