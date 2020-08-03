#TODO
#maybe bring into one function ie aoo(...........,method="simple,fixed,iterative")
#look at removing dataframes from rotation and only using matrix will give x 4 speed increase
#look a lapply for the both loops https://gist.github.com/dsparks/3706541
#move to functions
#add meta data to the SF so they have the cellsize,rotation,x and y shift embeded just incase it's needed somewhere



################Functions#####################

###########################################################
#calculates the optimum AOO (smallest) by shifting the grid
###########################################################
#' @title Optimal shifting grid, Area of Occupancy (AOO) 
#' @description
#' Calculates the optimal (smallest) Area of Occupancy AOO  by shifting the grid in x and y direction only.
#' The minimum solution will be achieved but large point datasets (i.e. over 70 points) will take some time to process.
#' Processing time is proportional to nop^2 (number of points squared).
#' Please cite below if using this algorithm:  
#' Moat, J., Bachman, S. P., Field, R., & Boyd, D. S. (2018). Refining area of occupancy to address the modifiable areal unit problem in ecology and conservation. Conservation biology, 32(6), 1278-1289.
#' @author Justin Moat. J.Moat@kew.org
#' @param thepoints dataframe of points in metres i.e. c(x,y)
#' @param cellsize width of cell in metres (default 2000 m) 
#' @param returnV, switches to return different sets of results:  
#' S = Simple, returns just the minimum are in km2, (DEFAULT)  
#' E = Expended simple, returns the solution for the smallest AOO  as list of (area,number of cells, rotation (0 in this case), shift in x direction, shift in y direction) 
#' ALL = returns a dataframe of all of the results with (number of cells, rotation (0 in this case), shift in x direction (metres), shift in y direction (metres))  
#' SF = returns a polygon simple feature for mapping, plotting in ggplot or export to GIS systems  
#' @return dependent on switch, default is area in km2
#' @note will give a warning if large dataset are used, if this take too long think about using: aooFixedRotation
#' @examples
#'#Build and project some points
#'thepoints <- squareOfPs(19,0.1)
#'names(thepoints) <- c("lat","long")
#'thepoints <- simProjWiz(thepoints)
#'attr(thepoints,'crs')
#'cellsize = 2000
#'
#'#just get area in km2
#'aooFixedGrid(thepoints,cellsize)
#'#extended return
#'aooFixedGrid(thepoints,cellsize,"E")
#'#Build polygons for plotting
#'library(ggplot2)
#'polys <- aooFixedGrid(thepoints,cellsize,"SF")
#'ggplot(data=polys) + geom_sf(color="black",fill="red") +geom_point(data=thepoints,aes(x,y))
#'#get it all results for all iterations
#'resultsdf <- aooFixedGrid(thepoints,cellsize,"ALL")
#'#NB as this is the optimal algorithm, do not expect the histogram to be normal
#'ggplot(data=resultsdf,aes(nofcells)) + geom_histogram()
#'#range of cellsize returned
#'range(resultsdf$nofcells)
#'#Example of building polygons from one set of results
#'#plotting the result with the highest AOO
#'worstgrid <- resultsdf[which.max(resultsdf$nofcells),]
#'worstpoly <- buildCellPolys_rxy(thepoints,cellsize,worstgrid$rotation,worstgrid$xshift,worstgrid$yshift)
#'ggplot(data=worstpoly) + geom_sf(color="black",fill="red") +geom_point(data=thepoints,aes(x,y))
#' @seealso \code{\link{aooRating}} for AOO Ratings
#' @seealso \code{\link{aoo}} for simple AOO method
#' @seealso \code{\link{aooFixedRotation}} for systematic methods
#' @seealso \code{\link{buildCellPolys_rxy}} for building grid polygons from points, rotation and shift
#' @export
#' @import sf
#' @references
#' Moat, J., Bachman, S. P., Field, R., & Boyd, D. S. (2018). Refining area of occupancy to address the modifiable areal unit problem in ecology and conservation. Conservation biology, 32(6), 1278-1289.
#' 
#' Bachman, S., Moat, J., Hill, A.W., de Torre, J., Scott, B., 2011. Supporting Red List threat assessments with GeoCAT: geospatial conservation assessment tool. Zookeys 126, 117–26. doi:10.3897/zookeys.150.2109 


aooFixedGrid <- function(thepoints,cellsize=2000,returnV="S"){
  #Warning for large dataset, user can at least be warned and kill
  if (nrow(thepoints) > 70){ warning(paste("This will run",nrow(thepoints)^2,"times - it may take some time!"),immediate. = TRUE)}
  #starting grid 0,0 just to get  marker (not used later)
  cpoints <- thepoints/cellsize
  bestpoints <- unique(floor(cpoints)) 
  minN <- nrow(bestpoints)
  #setup varibles to populate with the 2x loops
  resultsdf <- data.frame(nofcells=numeric(),rotation=numeric(),xshift=numeric(),yshift=numeric())
  t <- 0
  shiftx <- 0
  shifty <- 0
  #the sequences from left of point
  xpon <- cpoints$X - floor(cpoints$X)
  ypon <- cpoints$Y - floor(cpoints$Y)
  for (i in xpon){
    for (j in ypon){
      testps <- cbind(cpoints$X - i,cpoints$Y - j)
      testcps <- unique(floor(testps))
      t <- t+1
      resultsdf[t,] <- c(nrow(testcps),0,i*cellsize,j*cellsize)
      if (nrow(testcps)< minN){
        bestpoints<- cbind(testcps[,1] + i, testcps[,2] +j) 
        minN <- nrow(bestpoints)
      }
    }
  }
  
  #get the first minimum grid for results and returns
  bestgrid <- resultsdf[which.min(resultsdf$nofcells),]
  if(returnV == "E"){
    return(list(area=bestgrid$nofcells * (cellsize^2)/1000000,nocells=bestgrid$nofcells, 
                rotation = rad2deg(bestgrid$rotation),
                xshift= bestgrid$xshift, yshift = bestgrid$yshift))
    }
  if(returnV == "SF"){
    #build df of best points for SF production
    bestpoints <- data.frame (bestpoints * cellsize)
    colnames(bestpoints)<-c("x","y")
    return(buildCells(bestpoints,cellsize,0,shiftx,shifty,attr(thepoints,'crs')))
    }
  if(returnV == "ALL"){return(resultsdf)}
  else {return(minN * (cellsize^2)/1000000)}
}

###########################################################
#Calculates  rotation and shift, using it iterations      #
#Split 50:50  rotation : shift                            #
###########################################################
#' @title Systematic shifting and rotating grid for Area of Occupancy (AOO) 
#' @description
#' Calculates the Area of Occupancy AOO (smallest) by shifting and rotating the grid in x and y direction only.
#' In a very few occasions the minimum solution will not always be achieved but it is quick and consistent (not driven by the number of points).
#' If your species is near a threshold you may want to increase the number of iterations.
#' Please cite below if using this algorithm:  
#' Moat, J., Bachman, S. P., Field, R., & Boyd, D. S. (2018). Refining area of occupancy to address the modifiable areal unit problem in ecology and conservation. Conservation biology, 32(6), 1278-1289.
#' @author Justin Moat. J.Moat@kew.org
#' @param thepoints dataframe of points in metres i.e. c(x,y)
#' @param it the number of iterations you wish it to run, (default 1296)
#' @param cellsize width of cell in metres (default 2000 m) 
#' @param returnV, switches to return different sets of results:  
#' S = Simple, returns just the minimum are in km2, (DEFAULT)  
#' E = Expended simple, returns the solution for the smallest AOO  as list of (area,number of cells, rotation (degrees), shift in x direction, shift in y direction)  
#' ALL = returns a dataframe of all of the results with (number of cells, rotation (radians), shift in x direction (metres), shift in y direction (metres))  
#' SF = returns a polygon simple feature for mapping, plotting in ggplot or export to GIS systems.
#' @param rotation allow rotation of grids? (default = TRUE). If rotations are selected iterations are shared 50:50 rotation:shift
#' @return dependent on switch, default is area in km2
#' @note
#' @examples
#'#Build and project some points
#'thepoints <- squareOfPs(19,0.1)
#'names(thepoints) <- c("lat","long")
#'thepoints <- simProjWiz(thepoints)
#'attr(thepoints,'crs')
#'cellsize = 2000
#'
#'#just get minimum AOO area in km2
#'aooFixedRotation(thepoints,cellsize,1296,"S")
#'#get extended results for minimum AOO
#'aooFixedRotation(thepoints,cellsize,1296,"E")
#'#get results for all iterations
#'a <- aooFixedRotation(thepoints,cellsize,1296,"ALL")
#'hist(a$nofcells)
#'#Build polygon for plotting
#'polysrot <- aooFixedRotation(thepoints,cellsize,1296,"SF")
#'ggplot(data=polysrot) + geom_sf(color="black",fill="red") +
#'   geom_point(data=thepoints,aes(X,Y))
#'#Build polygons and plot the worst case
#'worstgrid <- all[which.max(all$nofcells),]
#'worstpoly <- buildCellPolys_rxy(thepoints,cellsize,worstgrid$rotation,worstgrid$xshift,worstgrid$yshift)
#'ggplot(data=worstpoly) + geom_sf(color="black",fill="red") +
#'   geom_point(data=thepoints,aes(X,Y))
#'
#' @seealso \code{\link{aooRating}} for AOO Ratings
#' @seealso \code{\link{aooFixedGrid}} for fixed grid optimal method
#' @seealso \code{\link{aoo}} for simple AOO method
#' @seealso \code{\link{buildCellPolys_rxy}} for building grid polygons from points, rotation and shift
#' @export
#' @import sf
#' @references
#' Moat, J., Bachman, S. P., Field, R., & Boyd, D. S. (2018). Refining area of occupancy to address the modifiable areal unit problem in ecology and conservation. Conservation biology, 32(6), 1278-1289.
#' 
#' Bachman, S., Moat, J., Hill, A.W., de Torre, J., Scott, B., 2011. Supporting Red List threat assessments with GeoCAT: geospatial conservation assessment tool. Zookeys 126, 117–26. doi:10.3897/zookeys.150.2109 

aooFixedRotation <- function(thepoints,cellsize=2000,it=1296,returnV="S",rotation=TRUE){
  #Build iterative lists
  if (rotation){
    rotlist <- seq(from = 0, to = pi/2, length=it^(1/2))
    xpon <- seq(from=0, to=cellsize, length=it^(1/4))
    ypon <- seq(from=0, to=cellsize, length=it^(1/4))
    } else {
    rotlist <- 0
    xpon <- seq(from=0, to=cellsize, length=it^(1/2))
    ypon <- seq(from=0, to=cellsize, length=it^(1/2))
    }
  #setup variables for storage and testing within loops
  minn <- nrow(thepoints)
  #bestpoints <- unique(floor(thepoints/cellsize))
  shiftx <- 0
  shifty <- 0
  minr <- 0
  t <- 0
  resultsdf <- data.frame(nofcells=numeric(),rotation=numeric(),xshift=numeric(),yshift=numeric())
  #the big ^3 loop
  for (i in xpon){
    for (j in ypon){
      for (r in rotlist){
        #shift points
        testps <- cbind(thepoints$X-i,thepoints$Y-j)
        #rotate points
        if(rotation){
          rps <- rotateP(testps,r)
        } else {
          rps<-data.frame(x=testps[,1],y=testps[,2])
          }
        rcells <- unique(floor(rps/cellsize)) * cellsize 
        t <- t + 1
        resultsdf[t,] <- c(nrow(rcells),r,i,j)
        if (nrow(rcells) < minn){
          minn <- nrow(rcells)
          bestpoints<- cbind(rcells[,1], rcells[,2])
          shiftx <- i 
          shifty <- j 
          minr <- r
        }
        
      }
    }
  }
  bestpoints <- data.frame (bestpoints)
  colnames(bestpoints)<-c("x","y")
  bestgrid <- resultsdf[which.min(resultsdf$nofcells),]
  if(returnV == "E"){return(list(area=bestgrid$nofcells * (cellsize^2)/1000000,nocells=bestgrid$nofcells, rotation =rad2deg(bestgrid$rotation),
                                 xshift= bestgrid$xshift, yshift = bestgrid$yshift))}
  if(returnV == "ALL"){return(resultsdf)}
  if(returnV == "SF"){buildCells(bestpoints,cellsize,-minr,shiftx,shifty,attr(thepoints,'crs'))}
  else{return(bestgrid$nofcells * (cellsize^2)/1000000)}
}



###########################################################
#calculates the intial AOO, with simple grid 0,0 #
###########################################################
#' calculates a very simple AOO area from a set of points
#' @title Area of Occupancy (AOO)
#' @description Calculates the number of occupied cells for Area of Occupancy from a set of points (x,y), in metres, with origin 0,0.
#'  cite below if using this algorithm: Moat, J., Bachman, S. P., Field, R., & Boyd, D. S. (2018). Refining area of occupancy to address the modifiable areal unit problem in ecology and conservation. Conservation biology, 32(6), 1278-1289.
#' @author Justin Moat. J.Moat@kew.org
#' @param thepoints set of points in metres i.e. c(x,y)
#' @param cellsize size of cell (length) in metres
#' @param returnV, switches to return different sets of results:
#' S = Simple, returns just the minimum are in km2, (DEFAULT)
#' E = Expended simple, returns the solution for the smallest AOO as list of (area,number of cells, rotation (degrees), shift in x direction, shift in y direction)
#' SF = returns a polygon simple feature for mapping, plotting in ggplot or export to GIS systems
#' @return as returnV, default is area in km2
#' @examples
#'library(ggplot2)
#'#Build and project some points
#'thepoints <- squareOfPs(19,0.1)
#'names(thepoints) <- c("lat","long")
#'thepoints <- simProjWiz(thepoints)
#'attr(thepoints,'crs')
#'cellsize = 2000
#'
#'#return area in km2
#'aoo (thepoints,cellsize)
#'#return list of parameters
#'aoo (thepoints,cellsize,returnV="E")
#'#return polygon for plotting
#'gridpoly <- aoo(thepoints,cellsize,returnV="SF")
#'ggplot(data=gridpoly) 
#'   + geom_sf(color="black",fill="red") 
#'   + geom_point(data=thepoints,aes(X,Y))
#' @seealso \code{\link{aooRating}} for AOO Ratings from IUCN categories
#' @export
#' @references
#' Moat, J., Bachman, S. P., Field, R., & Boyd, D. S. (2018). Refining area of occupancy to address the modifiable areal unit problem in ecology and conservation. Conservation biology, 32(6), 1278-1289.
#' 
#' Bachman, S., Moat, J., Hill, A.W., de Torre, J., Scott, B., (2011). Supporting Red List threat assessments with GeoCAT: geospatial conservation assessment tool. Zookeys 126, 117–26. doi:10.3897/zookeys.150.2109 

aoo <- function(thepoints,cellsize=2000,returnV="S"){
  bottomleftpoints <- unique(floor(thepoints/cellsize))
  cellp <- data.frame(x=(bottomleftpoints$X * cellsize), y=(bottomleftpoints$Y * cellsize))
  if(returnV == "E"){return(list(area=nrow(cellp) * (cellsize^2)/1000000,nrow(cellp), rotation =0,
                                 xshift= 0, yshift = 0))}
  if(returnV == "SF"){buildCells(cellp,cellsize,0,0,0,attr(thepoints,'crs'))}
  else{return(nrow(cellp) * (cellsize^2)/1000000)}
}
###for consistency with old rCAT 0.1.6
AOOsimp <- aoo


###############building blocks for polygon production######



###########################################################
#builds polygons from points and rotation, shift in X and y
#returns polygons for ggplot2 and mapping
###########################################################
#' @title Build cell polygons from points, rotation and shift in x and y
#' @description 
#' Builds square cell polygons from points and rotation, shift in X and y returns polygons for ggplot2 and mapping.
#' 
#' @author Justin Moat. J.Moat@kew.org
#' @param thepoints set of points in metres i.e. c(x,y)
#' @param cellsize size of cell (length) in metres
#' @param rot rotation of the grid in radian
#' @param shiftx shift in the x direction in metres
#' @param shifty shift in the y direction in metres
#' @return Simple Feature of polygons
#' @examples
#'#Build and project some points
#'thepoints <- squareOfPs(19,0.1)
#'names(thepoints) <- c("lat","long")
#'thepoints <- simProjWiz(thepoints)
#'#Check projection information is attributed
#'attr(thepoints,'crs')
#'cellsize = 2000
#'all <- aooFixedRotation(thepoints,cellsize,1296,"ALL")
#'worstgrid <- all[which.max(all$nofcells),]
#'worstpoly <- buildCellPolys_rxy(thepoints,cellsize,worstgrid$rotation,worstgrid$xshift,worstgrid$yshift)
#'ggplot(data=worstpoly) 
#'   + geom_sf(color="black",fill="red") 
#'   + geom_point(data=thepoints,aes(X,Y))
#'
#' @export
#' @references
#' Moat, J., Bachman, S. P., Field, R., & Boyd, D. S. (2018). Refining area of occupancy to address the modifiable areal unit problem in ecology and conservation. Conservation biology, 32(6), 1278-1289.

buildCellPolys_rxy<- function(thepoints,cellsize,rot,shiftx,shifty){
  #shift first
  testps <- cbind(thepoints$X - shiftx,thepoints$Y - shifty)
  #then rotate
  rpoints <- rotateP(testps,rot)
  testcps <- unique(floor(rpoints/cellsize))*cellsize
  colnames(testcps)<-c("x","y")
  buildCells(testcps,cellsize,-rot,shiftx,shifty,attr(thepoints,'crs'))
  
}





###########################################################
#Rotates a set of points                                  #
###########################################################
#Note angle in radians and only needed between 0 and 2pi for 360's
#but if using with shift you really only need 0 and pi/2
#' @title Rotates a set of points
#' @description 
#' Rotates a set of point by an angle in radians. Used as part of the AOO rotation calculations.
#' 
#' @author Justin Moat. J.Moat@kew.org
#' @param thepoints set of points in metres i.e. c(x,y)
#' @param angle in radians
#' @return dataframe of points
#' @examples
#'#Build and project some points
#'thepoints <- squareOfPs(200,0.1)
#'#rotate by 45 degrees
#'rpoints <- rotateP(thepoints,deg2rad(45))
#'plot(rpoints,asp=1)
#' @export
#' @references
#' Moat, J., Bachman, S. P., Field, R., & Boyd, D. S. (2018). Refining area of occupancy to address the modifiable areal unit problem in ecology and conservation. Conservation biology, 32(6), 1278-1289.

rotateP <- function(thepoints, angle){
  pointmatrix <- as.matrix(thepoints)
  rotationmatrix <- matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)),byrow = TRUE, 2, 2)
  pr <- pointmatrix %*% rotationmatrix
  pdf <- as.data.frame(pr)
  colnames(pdf) <- c("x","y")
  return(pdf)
}

###########################################################
#builds all corners for a square from the lower left corner, returns df id,x,y for use in ggplot
###########################################################
buildCellPolys <- function (llcorners,cellsize){
  mydf <- data.frame(id=integer(),x=double(),y=double())
  for (i in 1:nrow(llcorners)){
    mydf[nrow(mydf)+1,] <- c(i,llcorners[i,]$x,llcorners[i,]$y)
    mydf[nrow(mydf)+1,] <- c(i,llcorners[i,]$x+cellsize,llcorners[i,]$y)
    mydf[nrow(mydf)+1,] <- c(i,llcorners[i,]$x+cellsize,llcorners[i,]$y+cellsize)
    mydf[nrow(mydf)+1,] <- c(i,llcorners[i,]$x,llcorners[i,]$y+cellsize)
    mydf[nrow(mydf)+1,] <- c(i,llcorners[i,]$x,llcorners[i,]$y)
  }
  return(mydf)
}

###########################################################
#builds all corners for a square from the lower left corner,
#rotation, shift in X and y
#returns polygons
###########################################################
#internal called from main scripts

buildCells <- function (llcorners,cellsize,rot=0,shiftx=0,shifty=0,crstxt){
  #build cells 
  mincells <- buildCellPolys(as.data.frame(llcorners),cellsize)
  #rotate these back to original point orientation
  cells <- rotateP(mincells[,2:3],rot)
  #shift
  cells$x <- cells$x + shiftx
  cells$y <- cells$y + shifty
  #Add ID
  cells$id <- mincells$id
  ###cells is a df which ggplot can read, below make them sf spatial objects
  xys <- st_as_sf(cells,coords=c("x","y"))
  polys <- st_sf(
    aggregate(
      xys$geometry,
      list(xys$id),
      function(g){
        st_cast(st_combine(g),"POLYGON")
      }
    ))
  st_crs(polys) <- crstxt
  return(polys)
}

