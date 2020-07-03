###########################################################
#Simple example code for one species #
###########################################################
#' @title Example code for one species
#' @description 
#' Calculates EOO and AOO for one species
#' @author Justin Moat. J.Moat@kew.org
#' @note
#' @details 
#' @examples
#' library ("rCAT")
#' 
#' #read in csv data
#' mydata <- read.csv("species_x.csv")
#' mypointsll <- data.frame(lat=mydata$latitude,long=mydata$longitude)
#' 
#' #find the true centre of the points
#' centreofpoints <- trueCOGll(mypointsll)
#' #project to a equal area projection
#' mypointsxy <- simProjWiz(mypointsll,centreofpoints)
#' 
#' #Calculate the Extent of Occurance EOO
#' EOOm2 <- EOOarea(mypointsxy)
#' EOOkm2 <- -(EOOm2/1000000)
#' 
#' #Calculate the Area of Occupancy AOO
#' cellsizem <- 2000
#' AOOnocells <- AOOsimp (mypointsxy,cellsizem)
#' AOOkm2 <- AOOnocells * (cellsizem/1000)^2
#' x <- runif (20,0,10)
#' y <- runif (20,0,10)
#' df <- data.frame(x,y) 
#' AOOsimp (df)


library ("rCAT")

#read in csv data
mydata <- read.csv("species_x.csv")
mypointsll <- data.frame(lat=mydata$latitude,long=mydata$longitude)

#find the true centre of the points
centreofpoints <- trueCOGll(mypointsll)
#project to a equal area projection
mypointsxy <- simProjWiz(mypointsll,centreofpoints)

#Calculate the Extent of Occurance EOO
EOOm2 <- EOOarea(mypointsxy)
EOOkm2 <- -(EOOm2/1000000)

#Calculate the Area of Occupancy AOO
cellsizem <- 2000
AOOnocells <- AOOsimp (mypointsxy,cellsizem)
AOOkm2 <- AOOnocells * (cellsizem/1000)^2
