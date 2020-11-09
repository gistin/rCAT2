#NB just reset to new calls only

###########################################################
#Simple function to process multiple species #
###########################################################
#' @title Batch process, preliminary conservation assessments
#' @description 
#' Combines the main of routines in rCAT to process multiple species for AOO, EOO etc.
#' @author Justin Moat. J.Moat@kew.org
#' @note Has a switch to either project all data as a whole or each taxa separately. We would suggest you use this switch if data is all from a similar area (i.e. all from one continent/country/region)
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
#' @param long field which defines the longitude set of points
#' @param project2gether TRUE or FALSE, TRUE all data is projected together using the centre of all latitudes and longitudes. FALSE each species is projected separately. Default = TRUE
#' @param cellsize cell length in metres used to for AOO projection N.B. IUCN recommend 2000 m (default 2000)
#' @param aooMin calls the aooMin routines as well as simple aoo, be warned with lots of species and points this can take some time to run (default=FALSE)
#' @param it if aooMin=TRUE this determines the number of iterations it will run to find aooMin (default=1296)
#' @param returnV switches to return different sets of results: \cr
#' S = simple returns a dataframe of results = (taxa ,Number of points,EOO in km2, Simple AOO in km2,Minimum AOO, EOO category, AOOcategory, Cellwidth, projection parameters) \cr
#' SF = simple features dataframe will all results, taxa in taxa field, geometryclass=(EOO,AOO,points). NB all points will be projected together and aooMin is ignored
#' 
#' 
#' 
#' @examples
#' 
#' #getting batch working
#'#set up simple data within Cyprus
#'lat <- runif (200,32.8294,32.9394)
#'long <- runif (200,34.8720,34.9720)
#'mydata <- data.frame(taxa=c('aa','bb','cc','dd',"xx"),lat,long)
#'
#'#default get a dataframe of results and project all with the same project
#'resultsdf <- conBatch(mydata$taxa,mydata$long,mydata$lat)
#'#project each species individually
#'resultsdf <- conBatch(mydata$taxa,mydata$long,mydata$lat,project2gether = FALSE)
#'#switch on aooMin
#'resultsdf <- conBatch(mydata$taxa,mydata$long,mydata$lat,aooMin=TRUE)
#'
#'#default to return Simple feature objects
#'resultsf <- conBatch(mydata$taxa,mydata$long,mydata$lat,returnV = "SF")
#'#plot all the EOO results
#'library(ggplot2)
#'ggplot(data=resultsf[resultsf$geom_cat=="eoo",]) + geom_sf(fill=NA)
#'#pullone species and plot
#'oneSp <- resultsf[resultsf$taxa=="aa",]
#'ggplot(data=oneSp) + geom_sf(fill=NA)

#' @references 
#' Bachman, S., Moat, J., Hill, A.W., de Torre, J., Scott, B., 2011. Supporting Red List threat assessments with GeoCAT: geospatial conservation assessment tool. Zookeys 126, 117–26. doi:10.3897/zookeys.150.2109 
#' 
#' Moat, J. F. 2007. Conservation assessment tools extension for ArcView 3.x, version 1.2. Retrieved from http://www.kew.org/gis/projects/cats/catsdoc.pdf
#' 
#' Moat, J., Bachman, S.,(2017) GeoCAT Geospatial Conservation Assessment Tool [WWW Document]. URL http://geocat.kew.org/
#' 
#' Joppa, L.N., Butchart, S.H.M., Hoffmann, M., Bachman, S.P., Akçakaya, H.R., Moat, J.F., Böhm, M., Holland, R.A., Newton, A., Polidoro, B., Hughes, A., 2016. Impact of alternative metrics on estimates of extent of occurrence for extinction risk assessment. Conserv. Biol. 30, 362–370. doi:10.1111/cobi.12591
#' @export

batchCon <- function(taxa,long,lat,project2gether=TRUE,cellsize=2000,aooMin=FALSE,it=1296, returnV='S'){
  if (returnV == "SF") {
    project2gether = TRUE
  }
  
  points <- data.frame(lat, long)
  
  if (project2gether) {
    points <- simProjWiz(points)
  }
  
  split_points <- split(points, f=taxa)
  
  if(returnV=='SF'){
    crs <- attr(points, "crs")
    ntaxa <- length(unique(taxa))
    
    eoo_geoms <- lapply(split_points, function(p) eoo(p, "SF"))
    aoo_geoms <- lapply(split_points, function(p) aoo(p, cellsize, "SF"))
    aoo_geoms <- lapply(aoo_geoms, function(g) st_sfc(st_multipolygon(g), crs=crs))
    point_geoms <- lapply(split_points, function(p) st_sfc(st_multipoint(data.matrix(p)), crs=crs))
    
    geoms <- c(
      do.call(c, eoo_geoms),
      do.call(c, aoo_geoms),
      do.call(c, point_geoms)
    )
    
    results <- st_sf(
      taxon=rep(unique(taxa), 3),
      type=c(rep("eoo", ntaxa), rep("aoo", ntaxa), rep("points", ntaxa)),
      geometry=geoms
    )
    
  } else {
    if (! project2gether) {
      split_points <- lapply(split_points, simProjWiz)
      proj_strings <- lapply(split_points, function(p) attr(p, "crs"))
      proj_strings <- do.call(c, proj_strings)
    } else {
      proj_strings <- attr(points, "crs")
    }
    
    n_points <- lapply(split_points, nrow)
    eoo_areas <- lapply(split_points, eoo)
    eoo_ratings <- lapply(eoo_areas, ratingEoo)
    aoo_areas <- lapply(split_points, function(p) aoo(p, cellsize))
    
    if (aooMin) {
      min_aoo_areas <- lapply(split_points, function(p) aooFixedRotation(p, cellsize, it))
      aoo_ratings <- lapply(min_aoo_areas, ratingAoo)
    } else {
      aoo_ratings <- lapply(aoo_areas, ratingAoo)
    }
    
    results <- data.frame(
      taxon=unique(taxa),
      NOP=do.call(c, n_points),
      EOOkm2=do.call(c, eoo_areas),
      AOOkm=do.call(c, aoo_areas),
      EOOcat=do.call(c, eoo_ratings),
      AOOcat=do.call(c, aoo_ratings),
      cellwidth=cellsize,
      proj_metadata=proj_strings
    )
    
    if (aooMin) {
      min_aoo_areas <- lapply(split_points, function(p) aooFixedRotation(p, cellsize, it))
      results$MinAOO <- do.call(c, min_aoo_areas)
    }
    
  }
  
  results
}
##############



#old one left for compatablity with 1.6
#update for new function names

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
    MERps <- mer(taxapointsxy)/1000
    MERarea <- (MERps[2] - MERps[1]) * (MERps[4]- MERps[3])
    #EOO
    EOOa <- EOOarea(taxapointsxy)/-1000000
    #AOO with cellsize
    AOOa <- aoo(taxapointsxy,cellsize) * (cellsize/1000)^2
    #EOO IUCN category
    EOOcat <- ratingEoo(EOOa)
    #AOO IUCn category
    AOOcat <- ratingAoo(AOOa)
    #population the results data frame
    resultsdf[nrow(resultsdf)+1,] <- c(taxa,nop,MERarea,EOOa,AOOa,EOOcat,AOOcat)
    } #end of loop
  return (resultsdf)
} #end of function