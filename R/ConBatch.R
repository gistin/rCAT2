###########################################################
#Simple fuction to process multiple species #
###########################################################
#' @title Batch process, preliminary conservation assessments
#' @description 
#' Combines all of routines in rCAT to process multiple species for AOO, EOO etc.
#' @author Justin Moat. J.Moat@kew.org
#' @note Has a switch to either project all data as a whole or each taxa separately. I would suggest you use this switch if data is all from a simlar area (i.e. all from one country/region)
#' @details This function expects a list of taxa and latitudes and longitudes.ie\cr
#' \tabular{lll}{
#' species_w \tab 85.388000  \tab   84.33100\cr 
#' species_w \tab -45.467000   \tab  88.41500\cr
#' species_w \tab -34.339000 \tab -149.52600\cr
#' species_x \tab -29.620000   \tab  79.11900\cr
#' species_x \tab 33.409000  \tab  -33.94700\cr
#' species_x \tab 64.692000 \tab -149.18900\cr
#' species_y \tab 2.308000 \tab -140.21900\cr
#' species_y \tab 41.452000  \tab -3.65600\cr
#' species_y \tab -30.475000 \tab -129.99600\cr}
#' etc
#' @param taxa field which defines a list of species or taxa
#' @param lat field which defines the latitude set of points
#' @param long field which defines the longtitude set of points
#' @param cellsize cell length in metres used to for AOO projection N.B. IUCN recommend 2000 m (ie 2 km)
#' @param project2gether TRUE or FALSE, TRUE all data is projected together using the centre of all latitudes and longtitudes. FALSE each species is projected separately. Default = FALSE
#' @seealso \code{\link{MER}} for Minimum Enclosing Rectangle calculations, \code{\link{EOOarea}} for EOO calculations, \code{\link{EOORating}} for EOO Ratings, \code{\link{AOOsimp}} for AOO calculations, \code{\link{AOORating}} for AOO Ratings,
#' 
#' @return dataframe with; taxa ,Number of points, Area of the enclosing recetangle, EOO Area in km2, AOO area in km2, EOO IUCN category, AOO IUCN category\cr
#' \tabular{llllllll}{
#'  \tab      taxa\tab  NOP\tab              MER\tab           EOOkm2\tab AOO2km\tab EOOcat\tab AOOcat\cr
#' 1 \tab species x \tab  14 \tab 918562.259811711 \tab 585915.607417865 \tab     14 \tab     LC\tab     EN\cr
#' 2 \tab species y \tab 124 \tab 1717224.64389286 \tab 634149.482670821  \tab  124   \tab  LC  \tab   EN\cr
#' 3 \tab species z \tab  61 \tab 22622717.2314339 \tab 17839113.1220552  \tab   61   \tab  LC  \tab   EN\cr
#' 4 \tab species w \tab 1130 \tab 509390660.388499 \tab 506445176.073246 \tab  1130  \tab   LC \tab    VU\cr}
#' 
#' 
#' @examples
#' lat <- runif (200,-24,-12)
#' long <- runif (200,43,51)
#' spa <- rep('aa',50)
#' spb <- rep('bb',150)
#' mydata <- data.frame(species=c(spa,spb),lat,long)
#' resultsdf <- ConBatch(mydata$species,mydata$lat,mydata$long,2000,FALSE)
#'
#' @references 
#' Bachman, S., Moat, J., Hill, A.W., de Torre, J., Scott, B., 2011. Supporting Red List threat assessments with GeoCAT: geospatial conservation assessment tool. Zookeys 126, 117–26. doi:10.3897/zookeys.150.2109 
#' 
#' Moat, J. F. 2007. Conservation assessment tools extension for ArcView 3.x, version 1.2. Retrieved from http://www.kew.org/gis/projects/cats/catsdoc.pdf
#' 
#' Moat, J., Bachman, S.,(2017) GeoCAT Geospatial Conservation Assessment Tool [WWW Document]. URL http://geocat.kew.org/
#' 
#' Joppa, L.N., Butchart, S.H.M., Hoffmann, M., Bachman, S.P., Akçakaya, H.R., Moat, J.F., Böhm, M., Holland, R.A., Newton, A., Polidoro, B., Hughes, A., 2016. Impact of alternative metrics on estimates of extent of occurrence for extinction risk assessment. Conserv. Biol. 30, 362–370. doi:10.1111/cobi.12591
#' @export

#' 
#test values

#mydata <- read.csv(system.file("extdata","multiple_species.csv",package="rCAT"))
#ConBatch(mydata$scientificname,mydata$latitude,mydata$longitude,2000,FALSE)
#taxa <- mydata$scientificname
#lat <- mydata$latitude
#long <- mydata$longitude
#cellsize <- 2000
#project2gether <- TRUE

################################
ConBatch <- function(taxa,lat,long,cellsize,project2gether){
  #catch the default
  if(missing(project2gether)){
    project2gether = FALSE
  }
  
  #sort data out
  mypointsll <- data.frame(taxa,lat,long)
  specieslist <- unique(mypointsll$taxa)
  
  #project all together
  if (project2gether){
    mypointsint <- simProjWiz(mypointsll[,c("lat","long")],trueCOGll(mypointsll[,c("lat","long")]))
    mypointsxy <- data.frame(taxa=mypointsll$taxa, x=mypointsint$x, y=mypointsint$y)
  }
  
  #make a dataframe to store results
  resultsdf<- data.frame(taxa=character(),NOP=integer(), MER=double(),EOOkm2=double(),AOOkm=double(), EOOcat=character(), AOOcat=character(), stringsAsFactors=FALSE)
  #rename the AOOfield to use the cellsize in km
  names(resultsdf)[5] <- paste("AOO",cellsize/1000,"km",sep="")
  
  #loop thought all taxa
  for (taxa in specieslist){
    print(paste("Processing",taxa))
    #get just one taxa to work on, if already projected just taxa get if not projected, then project and then get taxa
    if (project2gether){
      taxapointsxy <- (mypointsxy[mypointsxy$taxa==taxa,c("x","y")])
    } else {
      taxapointsll <- (mypointsll[mypointsll$taxa==taxa,c("lat","long")])
      #reproject points
      taxapointsxy <- simProjWiz(taxapointsll,trueCOGll(taxapointsll))
    }
    
    #CALCULATE METRICS
    #number of points
    nop <- nrow(taxapointsxy)
    #Minimium enclosing rectangle
    MERps <- MER(taxapointsxy)/1000
    MERarea <- (MERps[2] - MERps[1]) * (MERps[4]- MERps[3])
    #EOO
    EOOa <- EOOarea(taxapointsxy)/-1000000
    #AOO with cellsize
    AOOa <- AOOsimp(taxapointsxy,cellsize) * (cellsize/1000)^2
    #EOO IUCN category
    EOOcat <- EOORating(EOOa)
    #AOO IUCn category
    AOOcat <- AOORating(AOOa)
    #population the results data frame
    resultsdf[nrow(resultsdf)+1,] <- c(taxa,nop,MERarea,EOOa,AOOa,EOOcat,AOOcat)
    } #end of loop
  return (resultsdf)
} #end of function