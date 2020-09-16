#trying to Increasing speed
library (rCAT)
library(ggplot2)
library(profvis)

#Build some normally distributed point data around the Troodos mountains ~ 10 km diametre
thepoints <- normalofPs(29,0.1)
#TODO maybe think about make all upper case as SF objects???
#shift to Troodos mountaions
thepoints <- data.frame(long = thepoints$X + 32.8794, lat = thepoints$Y + 34.9220)
#check to see where they are on a map
ggplot(data=Cyprus) + geom_sf() + geom_point(data=thepoints,aes(long,lat))
#project the points
ppts <- simProjWiz(thepoints)
##check crs
attr(ppts,'crs')

#speed test old version
start_time <- Sys.time()
for (i in 1:100){
  aooFixedGrid(ppts)  
}
Sys.time() - start_time
#takes ~ 28-32 secs

#new version
start_time <- Sys.time()
for (i in 1:100){
  aooFixedGridm(ppts)  
}
Sys.time() - start_time
#25 sec

thepoints <- ppts
#http://adv-r.had.co.nz/Profling.html
#matrix  version

profvis({
#set up to do the main work
cellsize=2000
  #Warning for large dataset, user can at least be warned and kill
  if (nrow(thepoints) > 70){ warning(paste("This will run",nrow(thepoints)^2,"times - it may take some time!"),immediate. = TRUE)}
  #starting grid 0,0 just to get  marker (not used later)
  ######################
  #work with matrix
  #####################
  thepoints <- as.matrix(thepoints)
  #####################
  cpoints <- thepoints/cellsize
  bestpoints <- unique(floor(cpoints)) 
  minN <- nrow(bestpoints)
  #setup varibles to populate with the 2x loops
  resultsdf <- data.frame(nofcells=numeric(),rotation=numeric(),xshift=numeric(),yshift=numeric())
  t <- 0
  shiftx <- 0
  shifty <- 0
  #the sequences from left of point
  xpon <- cpoints[,"X"] - floor(cpoints[,"X"])
  ypon <- cpoints[,"Y"] - floor(cpoints[,"Y"])
  for (i in xpon){
    for (j in ypon){
      testps <- cbind(cpoints[,"X"] - i,cpoints[,"Y"] - j)
      testcps <- unique(floor(testps))
      t <- t+1
      resultsdf[t,] <- c(nrow(testcps),0,i*cellsize,j*cellsize)
      if (nrow(testcps)< minN){
        bestpoints<- cbind(testcps[,1] + i, testcps[,2] +j) 
        minN <- nrow(bestpoints)
      }
    }
  }
  
})

#################
profvis({
  #set up to do the main work
  cellsize <- 2000
  #Warning for large dataset, user can at least be warned and kill
  if (nrow(thepoints) > 70){ warning(paste("This will run",nrow(thepoints)^2,"times - it may take some time!"),immediate. = TRUE)}
  #starting grid 0,0 just to get  marker (not used later)
  ######################
  #work with matrix
  #####################
  thepoints <- as.matrix(thepoints)
  #####################
  cpoints <- thepoints/cellsize
  bestpoints <- unique(floor(cpoints)) 
  minN <- nrow(bestpoints)
  #setup varibles to populate with the 2x loops
  resultsdf <- data.frame(nofcells=numeric(),rotation=numeric(),xshift=numeric(),yshift=numeric())
  t <- 0
  shiftx <- 0
  shifty <- 0
  #the sequences from left of point
  xpon <- cpoints[,"X"] - floor(cpoints[,"X"])
  ypon <- cpoints[,"Y"] - floor(cpoints[,"Y"])
  ###################################################
  #this bit need to be moved to a map or apply
  ###################################################
  work <-function(i,j){
    testps <- cbind(cpoints[,"X"] - i,cpoints[,"Y"] - j)
    testcps <- unique(floor(testps))
    c(nrow(testcps),0,i*cellsize,j*cellsize)
  }
  
  #only gives 29 solutions I think I need to build the actual list
  xponl <- rep(xpon,each=29)
  yponl <- rep(ypon,29)
  
  ###############
  t <- mapply (work,xponl,yponl,SIMPLIFY = T)
  #ok seems to work makes a matrix which is 
  
  #ok seems to work
  
  for (i in xpon){
    for (j in ypon){
      testps <- cbind(cpoints[,"X"] - i,cpoints[,"Y"] - j)
      testcps <- unique(floor(testps))
      t <- t+1
      resultsdf[t,] <- c(nrow(testcps),0,i*cellsize,j*cellsize)
      if (nrow(testcps)< minN){
        bestpoints<- cbind(testcps[,1] + i, testcps[,2] +j) 
        minN <- nrow(bestpoints)
      }
    }
  }
  
})



