#' @title Extent of occurrence (EOO) from a set of points - accounting for uncertainty
#'
#' @description Calculates the extent of occurrence (EOO) for a set of points using a minimum convex polygon. 
#' Takes into account uncertainty in the point data due to georeference or identification uncertainty
#' Error buffers are drawn with a set buffer distance for all points or a user defined distance for each point. 
#' Minimum, maximum, average and standard deviation are calculated by repeated uniform sampling within each error buffer
#' Uncertainty can also calculated for points that are deemed uncertain 
#' @keywords extent of occurrence
#' @import sf
#' @return dataframe with EOO statistics
#' @export


# parameters: 
# x y points
# number of samples
# add option if user wants to define single distance error?
# add option if user wants to see impact of uncertain occurrences being removed

#build set of points
#thepoints <- ptsNormal(50,10)
#shift to Peru, Colombia and Brazil
#thepoints <- data.frame(long = thepoints$X - 69.9465972, lat = thepoints$Y - 4.2260167)
#thepoints$error <- "50"
#thepointsEooError <- eooError(thepoints$lat, thepoints$long, thepoints$error, runs = 100)

#plot(thepoints$long, thepoints$lat)

#lat = dypsis_raw$latitude
#long = dypsis_raw$longitude
#error = dypsis_raw$coordinateuncertaintyinmeters
#runs = 1000

#thepointsEooError <- eooError(lat, long, error, runs)

# eooError - generate the min, max and average EOO based on sampling within error parameters 
eooError <- function(long, lat, error, runs){
  
  # project the points 
  mydf <- simProjWiz(data.frame(long,lat))
  
  # add the error back in
  mydf <- cbind(mydf, error = as.integer(error))
  
  # generate replicates of points data shifted to random error
  reps <- replicate(n = runs, pointShift(mydf, mydf$err), simplify = F)
  
  # run replicates against eoo function
  #testapply = lapply(test, eoo)
  runEoo = purrr::map_dbl(reps, eoo)
  
  print(summary(runEoo))
  
  # histogram of results
  h <- hist(runEoo, col="blue", xlab="EOO area", 
       main="Histogram of EOO")
  xfit<-seq(min(runEoo),max(runEoo),length=40) 
  yfit<-dnorm(xfit,mean=mean(runEoo),sd=sd(runEoo)) 
  yfit <- yfit*diff(h$mids[1:2])*length(runEoo) 
  lines(xfit, yfit, col="red", lwd=2)
  
  
  return(summary(runEoo))
  
}
