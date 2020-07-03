###########################################################
#Simple example code for multiple species. ie batch processing #
###########################################################
#' @title Example code for multiple species
#' @description 
#' Simple example code for multiple species. ie batch processing 
#' @author Justin Moat. J.Moat@kew.org


library ("rCAT")

#read in simple csv data
mydata <- read.csv(system.file("extdata","multiple_species.csv",package="rCAT"))
mypointsll <- data.frame(taxa=mydata$scientificname,lat=mydata$latitude,long=mydata$longitude)
specieslist <- unique(mypointsll$taxa)

#make a dataframe to store results
resultsdf<- data.frame(taxa=character(),NOP=integer(), MER=double(),EOOkm2=double(),AOO2km=double(),A00100km=double(),stringsAsFactors=FALSE)

for (taxa in specieslist){
  print(paste("Processing",taxa))
  #get just one taxa to work on
  taxapointsll <- (mypointsll[mypointsll$taxa==taxa,c("lat","long")])
  #reproject points
  taxapointsxy <- simProjWiz(taxapointsll,trueCOGll(taxapointsll))
  
  #CALCULATE METRICS
  #number of points
  nop <- nrow(taxapointsxy)
  #Minimium enclosing rectangle
  MERps <- MER(taxapointsxy)/1000
  MERarea <- (MERps[2] - MERps[1]) * (MERps[4]- MERps[3])
  #EOO
  EOOa <- EOOarea(taxapointsxy)/-1000000
  #AOO with 2 x 2 km cells
  AOO2km <- AOOsimp(taxapointsxy,2000)
  #AOO with 100 x 100 km cell
  AOO10km <- AOOsimp(taxapointsxy,100000)
  
  #population the results data frame
  resultsdf[nrow(resultsdf)+1,] <- c(taxa,nop,MERarea,EOOa,AOO2km,AOO10km)
}

#finally to export the results
write.csv(resultsdf,file="results.csv")

