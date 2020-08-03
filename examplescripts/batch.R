#getting batch working
#set up simple data within Cyprus
lat <- runif (200,32.8294,32.9394)
long <- runif (200,34.8720,34.9720)
mydata <- data.frame(taxa=c('aa','bb','cc','dd',"xx"),lat,long)

#default get a dataframe of results and project all with the same project
resultsdf <- conBatch(mydata$taxa,mydata$long,mydata$lat)
#project each species individually
resultsdf <- conBatch(mydata$taxa,mydata$long,mydata$lat,project2gether = FALSE)
#switch on aooMin
resultsdf <- conBatch(mydata$taxa,mydata$long,mydata$lat,aooMin=TRUE)
#default to return Simple feature objects
resultsf <- conBatch(mydata$taxa,mydata$long,mydata$lat,returnV = "SF")
#plot all the EOO results
library(ggplot2)
ggplot(data=resultsf[resultsf$geom_cat=="eoo",]) + geom_sf(color='green', fill=NA)
#pullone species and plot
oneSp <- resultsf[resultsf$taxa=="aa",]
ggplot(data=oneSp) + geom_sf(fill=NA)




#various switches
conBatch(mydata$taxa,mydata$long,mydata$lat,project2gether = TRUE,eooR=FALSE,aooR=FALSE,aooFixedRotationR = TRUE)


#############################functions

batchCon <- function(taxa,long,lat,project2gether=TRUE,cellsize=2000,aooMin=FALSE,it=1296, returnV='S'){
  #taxa <- mydata$taxa
  #long <- mydata$long
  #lat <- mydata$lat
  #project2gether=FALSE
  
  #project all either when specified of when returning sf's
  if(project2gether || returnV=="SF"){
    mypointsxy <- simProjWiz(data.frame(lat,long))
    mypointsxy$taxa <- taxa
  } else {
    mypointsll <- data.frame(taxa,lat,long)
    }
  #taxa list for loop
  taxalist <- unique(taxa)

  #return SF if asked for
  if(returnV=='SF'){
  resultsf <- st_sf(st_sfc(crs = attr(mypointsxy,'crs')))  
    for (thetaxa in taxalist){
      #project one at a time
        ppts <- mypointsxy[mypointsxy$taxa==thetaxa,]
        ppts <- ppts[1:2]
        attr(ppts,'crs') <- attr(mypointsxy,'crs')
      #do the work
      eoopoly <- eoo(ppts,"SF")
      aoopolys <- aoo(ppts,cellsize,"SF")
      #build the points
      mpts = st_sf(geometry = st_sfc(st_multipoint(data.matrix(ppts))), crs = attr(mypointsxy,'crs'))
      eoopoly$geom_cat <- "eoo"
      aoopolys$geom_cat <- "aoo"
      #setup the points data so it matches above
      mpts$geom_cat <- "points"
      mpts$Group.1 <- 1
      #bind
      temp <- rbind(eoopoly,aoopolys,mpts)
      #add species name
      temp$taxa <- thetaxa
      resultsf <- rbind(resultsf,temp)
    
    }
  return(resultsf)
  }
  
  
    #dataframe to store results
    resultsdf<- data.frame(taxa=character(),NOP=integer(),EOOkm2=double(),
                           AOOkm=double(),MinAOO=double(),EOOcat=character(),
                           AOOcat=character(),cellwidth=integer(),proj_metadata=character(),
                           stringsAsFactors=FALSE)
  
  
  for (thetaxa in taxalist){
    #project one at a time
    if(!project2gether){
      pptll <- mypointsll[mypointsll$taxa==thetaxa,]
      ppts <- simProjWiz(data.frame(lat=pptll$lat,long=pptll$long))
    } else {
      ppts <- mypointsxy[mypointsxy$taxa==thetaxa,]
      ppts <- ppts[1:2]
      attr(ppts,'crs') <- attr(mypointsxy,'crs')
       }
    #do the work
      eooarea <- eoo(ppts)
      eooRatingtxt <- eooRating(eooarea)
      aooarea <- aoo(ppts,cellsize)
      aooRatingtxt <- aooRating(aooarea)
    if(aooMin){
      minaooarea <- aooFixedRotation(ppts,cellsize,it)
      aooRatingtxt <- aooRating(minaooarea)
    } else {
      minaooarea <- NA
      }
    resultsdf[nrow(resultsdf)+1,] <- c(thetaxa, nrow(ppts), eooarea,
                                       aooarea, minaooarea, eooRatingtxt,
                                       aooRatingtxt, cellsize, attr(ppts,'crs'))
  }
  return(resultsdf)
}
##############