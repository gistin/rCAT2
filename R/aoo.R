#TODO
# Once tested remove old (aooFixedGrid & aooFixedRotation) and update with aooFixedRotationo aooFixedGrido
#possibles
#maybe bring into one function ie aoo(...........,method="simple,fixed,iterative")
#Maybe add meta data to the SF so they have the cellsize,rotation,x and y shift embeded just incase it's needed somewhere



################Functions#####################

###########################################################
#calculates the optimum AOO (smallest) by shifting the grid
###########################################################
#' @title Area of Occupancy (AOO), optimal shifting grid,  
#' @description
#' Calculates the optimal (smallest) Area of Occupancy AOO  by shifting the grid in the x and y direction only.
#' The minimum solution will be achieved but large point datasets (i.e. over 100 points) will take some time to process.
#' Processing time is proportional to (number of points squared (n^2)).
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
#'thepoints <- ptsSquare(19,0.1)
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
#' @seealso \code{\link{ratingAoo}} for AOO Ratings
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
#' @title Area of Occupancy (AOO) calculated by Systematic shifting and rotating of the grid 
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
#'thepoints <- ptsSquare(19,0.1)
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
#' @seealso \code{\link{ratingAoo}} for AOO Ratings
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
#calculates the initial AOO, with simple grid 0,0 #
###########################################################
#' calculates a very simple AOO area from a set of points
#' @title Area of Occupancy (AOO), grid orgin 0,0
#' @description 
#' Calculates the number area the of occupied cells for (Area of Occupancy AOO) from a set of points (x,y), projected into metres, with origin 0,0. 
#' Please cite: Moat, J., Bachman, S. P., Field, R., & Boyd, D. S. (2018). Refining area of occupancy to address the modifiable areal unit problem in ecology and conservation. Conservation biology, 32(6), 1278-1289. if using this algorithm:
#' @author Justin Moat. J.Moat@kew.org
#' @param thepoints set of points in metres i.e. c(x,y)
#' @param cellsize size of cell (length) in metres
#' @param returnV, switches to return different sets of results: \cr
#' S = Simple, returns just the AOO area in km2, (DEFAULT) \cr
#' E = Expended simple, returns the solution for the AOO as list of (area,number of cells, rotation (0 degrees), shift in x direction(0), shift in y direction(0)). This is used so as be compatiable with other AOO calculations. \cr
#' SF = returns a polygon simple feature for mapping and plotting in ggplot/plot or export to GIS format.
#' @return as returnV, default is area in km2
#' @examples
#'library(ggplot2)
#'#Build and project some points
#'thepoints <- ptsSquare(19,0.1)
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
#' @seealso \code{\link{ratingAoo}} for AOO Ratings from IUCN categories
#' @export
#' @references
#' Moat, J., Bachman, S. P., Field, R., & Boyd, D. S. (2018). Refining area of occupancy to address the modifiable areal unit problem in ecology and conservation. Conservation biology, 32(6), 1278-1289.
#' 
#' Bachman, S., Moat, J., Hill, A.W., de Torre, J., Scott, B., (2011). Supporting Red List threat assessments with GeoCAT: geospatial conservation assessment tool. Zookeys 126, 117–26. doi:10.3897/zookeys.150.2109 

aoo <- function(thepoints, cellsize=2000, returnV="S"){
  bottomleftpoints <- unique(floor(thepoints/cellsize))
  
  cellp <- data.frame(
    x=(bottomleftpoints$X * cellsize), 
    y=(bottomleftpoints$Y * cellsize)
  )
  
  if (returnV == "E") {
    return(list(
      area=nrow(cellp) * (cellsize^2)/1000000,
      nocells=nrow(cellp), 
      rotation=0,
      xshift=0, 
      yshift=0
    ))
  }
  
  
  if (returnV == "SF") {
    buildCells(cellp, cellsize, 0, 0, 0, attr(thepoints,'crs'))
  } else {
    return(nrow(cellp) * (cellsize^2)/1000000)
  }
}
###for consistency with old rCAT 0.1.6
AOOsimp <- aoo

########################################################################################################
#optimised AOO algorithms, I am sure they can be faster, but they are much better than above and will do

##########################################################
#BETA Area of Occupancy (AOO), from systematic grid rotation and shifting
###########################################################
#' @title BETA Area of Occupancy (AOO), from systematic grid rotation and shifting
#' @description 
#' Calculates the Area of Occupancy AOO (smallest) by rotating the grid and shifting in x and y direction. 
#' On a very few occasions the minimum solution will not always be achieved, but this solution is quick and consistent (not driven by the number of points). 
#' If your species is near a threshold you may want to increase the number of iterations. 
#' Please cite if using this algorithm: Moat, J., Bachman, S. P., Field, R., & Boyd, D. S. (2018). Refining area of occupancy to address the modifiable areal unit problem in ecology and conservation. Conservation biology, 32(6), 12781289. \cr
#' Works the same as aooFixedRotation, but much faster. In BETA until fully tested
#' On a very few occasions the minimum solution will not always be achieved but it is quick and consistent (not driven by the number of points).
#' If your species is near a threshold you may want to increase the number of iterations.
#' @author Justin Moat. J.Moat@kew.org
#' @param thepoints dataframe of points in metres i.e. c(x,y)
#' @param it the number of iterations you wish it to run, (default 1296)
#' @param cellsize width of cell in metres (default 2000 m) 
#' @param returnV, switches to return different sets of results: \cr  \cr
#' S = Simple, returns just the minimum AOO area in km2, (DEFAULT)  \cr
#' E = Expended simple, returns the solution for the smallest AOO  as list of (area,number of cells, rotation (degrees), shift in x direction, shift in y direction)  \cr
#' ALL = returns a dataframe of all of the results from all the trials with (number of cells, rotation (0 in this case), shift in x direction (metres), shift in y direction (metres))  \cr
#' SF = returns a polygon simple feature for mapping, plotting in ggplot/plot or to export to a GIS format.
#' @param rotation allow rotation of grids? (default = TRUE). If rotations are trigger selected iterations are shared 50:50 rotation:shift(both in x and y direction)
#' @return dependent on switch, default is area in km2
#' @note Number of Iteration as calculated to the nearest whole number i.e. round(it^0.25), so for 1296 iterations you get 36 rotations, 6 x shifts and 6 y shifts \cr
#' Up to 10,000 iterations (324r, 18x, 18y) is fine (slight pause), probably not worth going further unless you really really want to push for a smaller AOO (you still may not achieve it) \cr
#' Number of point will additional increase processing time \cr
#' If you are not interest in rotation the split is it^0.5 in x and y ie 1296:36x 36y or 104,976: 324x,324y \cr
#' If you have below 100 points and only interested in shifts then running \code{\link{aooFixedGrid}} will give you the optimal (true minimum) result in a reasonable processing time, but this method is systemic, which maybe more useful \cr
#' The break point for aooFixedRotationo vs aooFixedGrid default iterations is 36 points (i.e. both run with 1296 iterations) \cr
#' At 10,000 iterations the breakpoint is ~ 100 points \cr
#' Works the same as \code{\link\{aooFixedRotation}}, but much faster. In BETA until fully tested
#' On a very few occasions the minimum solution will not always be achieved but it is quick and consistent (not driven by the number of points).
#' If your species is near a threshold you may want to increase the number of iterations.
#' 
#' @examples
#'#Build and project some points
#'thepoints <- ptsSquare(19,0.1)
#'names(thepoints) <- c("lat","long")
#'thepoints <- simProjWiz(thepoints)
#'attr(thepoints,'crs')
#'cellsize = 2000
#'
#'#just get minimum AOO area in km2
#'aooFixedRotationo(thepoints,cellsize,6^4,"S")
#'#get extended results for minimum AOO and run more iterations
#'aooFixedRotationo(thepoints,cellsize,8^4,"E")
#'#get results for all iterations
#'a <- aooFixedRotationo(thepoints,cellsize,6^4,"ALL")
#'hist(a$nofcells)
#'#Build polygon for plotting
#'polysrot <- aooFixedRotationp(thepoints,cellsize,6^4,"SF")
#'ggplot(data=polysrot) + geom_sf(color="black",fill="red") +
#'   geom_point(data=thepoints,aes(X,Y))
#'#Build polygons and plot the worst case
#'worstgrid <- all[which.max(all$nofcells),]
#'worstpoly <- buildCellPolys_rxy(thepoints,cellsize,worstgrid$rotation,worstgrid$xshift,worstgrid$yshift)
#'ggplot(data=worstpoly) + geom_sf(color="black",fill="red") +
#'   geom_point(data=thepoints,aes(X,Y))
#'
#' @seealso \code{\link{ratingAoo}} for AOO Ratings
#' @seealso \code{\link{aooFixedGrid}} for fixed grid optimal method
#' @seealso \code{\link{aooFixedRotation}} for original method method
#' @seealso \code{\link{aoo}} for simple AOO method
#' @seealso \code{\link{buildCellPolys_rxy}} for building grid polygons from points, rotation and shift
#' @export
#' @import sf
#' @references
#' Moat, J., Bachman, S. P., Field, R., & Boyd, D. S. (2018). Refining area of occupancy to address the modifiable areal unit problem in ecology and conservation. Conservation biology, 32(6), 1278-1289.
#' 
#' Bachman, S., Moat, J., Hill, A.W., de Torre, J., Scott, B., 2011. Supporting Red List threat assessments with GeoCAT: geospatial conservation assessment tool. Zookeys 126, 117–26. doi:10.3897/zookeys.150.2109 

#' 

aooFixedRotationo <- function(thepoints,cellsize=2000,it=1296,returnV="S",rotation=TRUE){
  #get whole number for the iterations
  itf <- round(it^0.25)
  if (rotation){
    rotlist <- seq(from = 0, to = pi/2, length=itf^2)
    xpon <- seq(from=0, to=cellsize, length=itf)
    ypon <- seq(from=0, to=cellsize, length=itf)
  } else {
    rotlist <- 0
    xpon <- seq(from=0, to=cellsize, length=itf^2)
    ypon <- seq(from=0, to=cellsize, length=itf^2)
  }
  #setup mapply lists
  rl <- rep(rotlist,each=length(rotlist))
  xl <- rep(rep(xpon,length(rotlist)),each=length(xpon))
  yl <- rep(rep(ypon,length(rotlist)),length(ypon))
  #unique (cbind(rl,xl,yl))
  shiftrotgrid <- function(i,j,r){
    #shift
    testps <- cbind(thepoints$X-i,thepoints$Y-j)
    #rotate
    if (rotation){rps <- rotatePm(testps,r)} else {rps<-cbind(testps[,1],testps[,2])}
    cells <- floor(rps/cellsize)
    celltxt <- paste(cells[,1],cells[,2])
    rcells <- unique(celltxt)
    c(length(rcells),r,i,j)
  }
  mresults <- mapply (shiftrotgrid,xl,yl,rl)
  #convert to df
  resultsdf <- as.data.frame(t(mresults))
  names(resultsdf) <- c('nofcells','rotation','xshift','yshift')
  #get the first minimum grid for results and returns
  bestgrid <- resultsdf[which.min(resultsdf$nofcells),]
  if(returnV == "E"){return(list(area=bestgrid$nofcells * (cellsize^2)/1000000,nocells=bestgrid$nofcells, rotation =rad2deg(bestgrid$rotation),
                                 xshift= bestgrid$xshift, yshift = bestgrid$yshift))}
  if(returnV == "ALL"){return(resultsdf)}
  if(returnV == "SF"){buildCellPolys_rxy(thepoints,cellsize,bestgrid$rotation,bestgrid$xshift,bestgrid$yshift)}
  else{return(bestgrid$nofcells * (cellsize^2)/1000000)}
}




##############
#rotates a matrix of points
#much quicker that rotateP, but uses and return a matrix
#used for aooFixedRotationo
rotatePm <- function(thepoints, angle){
  #build rotation matrix
  rotationmatrix <- matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)),byrow = TRUE, 2, 2)
  pr <- thepoints %*% rotationmatrix
  return(pr)
}


###########################################################
#calculates the optimum AOO (smallest) by shifting the grid
###########################################################
#' @title BETA Area of Occupancy (AOO), optimal shifting grid,  
#' @description
#' Calculates the optimal (smallest) Area of Occupancy AOO  by shifting the grid in x and y direction only.
#' The minimum solution will be achieved but large point datasets (i.e. over 70 points) will take some time to process.
#' Processing time is proportional to nop^2 (number of points squared).
#' Please cite below if using this algorithm:  
#' Moat, J., Bachman, S. P., Field, R., & Boyd, D. S. (2018). Refining area of occupancy to address the modifiable areal unit problem in ecology and conservation. Conservation biology, 32(6), 1278-1289.
#' @author Justin Moat. J.Moat@kew.org
#' @param thepoints dataframe of points in metres i.e. c(x,y)
#' @param cellsize width of cell in metres (default 2000 m) 
#' @param returnV, switches to return different sets of results:  \cr \cr
#' S = Simple, returns just the minimum are in km2, (DEFAULT)  \cr
#' E = Expended simple, returns the solution for the smallest AOO  as list of (area,number of cells, rotation (0 in this case), shift in x direction, shift in y direction) \cr
#' ALL = returns a dataframe of all of the results with (number of cells, rotation (0 in this case), shift in x direction (metres), shift in y direction (metres))  \cr
#' SF = returns a polygon simple feature for mapping, plotting in ggplot or export to GIS systems  \cr
#' @return dependent on switch, default is area in km2
#' @note will give a warning if large numbers of points are used, if this take too long think about using: \code{\link{aooFixedRotation}}
#' @examples
#'#Build and project some points
#'thepoints <- ptsSquare(19,0.1)
#'names(thepoints) <- c("lat","long")
#'thepoints <- simProjWiz(thepoints)
#'attr(thepoints,'crs')
#'cellsize = 2000
#'
#'#just get area in km2
#'aooFixedGrido(thepoints,cellsize)
#'#extended return
#'aooFixedGrido(thepoints,cellsize,"E")
#'#Build polygons for plotting
#'library(ggplot2)
#'polys <- aooFixedGrido(thepoints,cellsize,"SF")
#'ggplot(data=polys) + geom_sf(color="black",fill="red") +geom_point(data=thepoints,aes(x,y))
#'#get it all results for all iterations
#'resultsdf <- aooFixedGrido(thepoints,cellsize,"ALL")
#'#NB as this is the optimal algorithm, do not expect the histogram to be normal
#'ggplot(data=resultsdf,aes(nofcells)) + geom_histogram()
#'#range of cellsize returned
#'range(resultsdf$nofcells)
#'#Example of building polygons from one set of results
#'#plotting the result with the highest AOO
#'worstgrid <- resultsdf[which.max(resultsdf$nofcells),]
#'worstpoly <- buildCellPolys_rxy(thepoints,cellsize,worstgrid$rotation,worstgrid$xshift,worstgrid$yshift)
#'ggplot(data=worstpoly) + geom_sf(color="black",fill="red") +geom_point(data=thepoints,aes(x,y))
#' @seealso \code{\link{ratingAoo}} for AOO Ratings
#' @seealso \code{\link{aoo}} for simple AOO method
#' @seealso \code{\link{aooFixedRotation}} for systematic methods with rotation
#' @seealso \code{\link{aooFixedGrid}} for orginal method
#' @seealso \code{\link{buildCellPolys_rxy}} for building grid polygons from points, rotation and shift
#' @export
#' @import sf
#' @references
#' Moat, J., Bachman, S. P., Field, R., & Boyd, D. S. (2018). Refining area of occupancy to address the modifiable areal unit problem in ecology and conservation. Conservation biology, 32(6), 1278-1289.


aooFixedGrido <- function(thepoints,cellsize=2000,returnV="S"){
  #Warning for large dataset, user can at least be warned and kill
  if (nrow(thepoints) > 100){ warning(paste("This will run",nrow(thepoints)^2,"times - it may take some time!"),immediate. = TRUE)}
  cpoints <- thepoints/cellsize
  ####################################
  #################
  #main function for mapply
  shiftgrid <- function(i,j){
    testps <- cbind(X=cpoints$X - i,Y=cpoints$Y - j)
    testcpi <- floor(testps)
    testcpi <- paste(testcpi[,1],testcpi[,2])
    testcps <- unique(testcpi)
    c(length(testcps),0,i*cellsize,j*cellsize)
  }
  #builds the actual lists for mapply
  xpon <- cpoints$X - floor(cpoints$X)
  ypon <- cpoints$Y - floor(cpoints$Y)
  xponl <- rep(xpon,each=nrow(thepoints))
  yponl <- rep(ypon,nrow(thepoints))
  mresults <- mapply (shiftgrid,xponl,yponl) 
  #reformat to df returning results etc
  resultsdf <- as.data.frame(t(mresults))
  ###################################
  names(resultsdf) <- c('nofcells','rotation','xshift','yshift')
  #get the first minimum grid for results and returns
  bestgrid <- resultsdf[which.min(resultsdf$nofcells),]
  if(returnV == "E"){
    return(list(area=bestgrid$nofcells * (cellsize^2)/1000000,nocells=bestgrid$nofcells, 
                rotation = rad2deg(bestgrid$rotation),
                xshift= bestgrid$xshift, yshift = bestgrid$yshift))
  }
  if(returnV == "SF"){
    bestpoly <- buildCellPolys_rxy(thepoints,cellsize,bestgrid$rotation,bestgrid$xshift,bestgrid$yshift)
    return(bestpoly)
  }
  ############################
  if(returnV == "ALL"){return(resultsdf)}
  else {return(bestgrid$nofcells * (cellsize^2)/1000000)}
}




###############building blocks for polygon production######
###########################################################
#builds polygons from points and rotation, shift in X and y
#returns polygons for ggplot2 and mapping
###########################################################
#' @title Build simple feature polygons from point data, rotation and shift in x and y direction
#' @description 
#' Builds cell polygons (as simple features) from points and rotation, shift in X and y returns polygons for ggplot2 and mapping.
#' Generally used to plot data from AOO calculations.
#' @author Justin Moat. J.Moat@kew.org
#' @param thepoints set of points in metres i.e. c(x,y)
#' @param cellsize size of cell (length) in metres
#' @param rot rotation of the grid in radian
#' @param shiftx shift in the x direction in metres
#' @param shifty shift in the y direction in metres
#' @return Simple Feature of polygons
#' @examples
#'#Build and project some points
#'thepoints <- ptsSquare(19,0.1)
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
#'thepoints <- ptsSquare(200,0.1)
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

buildCells <- function (llcorners, cellsize, rot=0, shiftx=0, shifty=0, crs=""){
  #build cells 
  mincells <- buildCellPolys(as.data.frame(llcorners),cellsize)
  #rotate these back to original point orientation
  cells <- rotateP(mincells[, 2:3], rot)
  #shift
  cells$x <- cells$x + shiftx
  cells$y <- cells$y + shifty
  
  cell_list <- split(cells, f=mincells$id)
  poly_list <- lapply(cell_list, function(x) constructPolygon(x$x, x$y, crs))
  
  do.call(c, poly_list)
}



##################################################################################
#calculates the number of Sub-population or number of locations from buffer method
##################################################################################
#' @title Area of Occupancy (AOO) point buffer method
#' @description 
#' Calculates the AOO area using the point buffer method. After Breiner, F. T., & Bergamini, A. (2018), they suggest this avoids issues with grid origin.
#' @author Justin Moat. J.Moat@kew.org
#' @param thepoints dataframe of points of x,y
#' @param bufferradius in metres, (default = 1128 m, therefore circular area = 2 km2)
#' @param returnV three switches for return values 
#' S = returns the area of the buffered output
#' SF = return the simple feature for plotting and mapping
#' 
#' @return AOO area in km2
#' @examples 
#'#Build some normally distributed point data around the Troodos mountains ~ 10 km diameter
#'thepoints <- ptsNormal(50,0.1)
#'#shift to Troodos mountaions
#'thepoints <- data.frame(long = thepoints$X + 32.8794, lat = thepoints$Y + 34.9220)
#'#project the points
#'ppts <- simProjWiz(thepoints)
#'#get the number of sub-pop/locations using default method
#'aooBuf(ppts)
#'#get the sf object and plot
#'aoobufs <- aooBuf(ppts,returnV="SF")
#'plot(aoobufs)
#'#user defined buffer distance in this case 2 km
#'aooBuf(ppts,bufferradius=2000,returnV="S")
#' @export
#' @import sf
#' @references 
#' Breiner, F. T., & Bergamini, A. (2018). Improving the estimation of area of occupancy for IUCN Red List assessments by using a circular buffer approach. Biodiversity and Conservation, 27(9), 2443-2448. 
#' @seealso \code{\link{subLocGrid}} which uses the same method, but different default radius and also calculates the number of sub-populations or locations

aooBuf <- function (thepoints,bufferradius=1128,returnV="S"){
  #just call subLocBuf 
  if (returnV == "S") {returnV <- "AREA"}
  subLocBuf(thepoints,bufferradius,returnV)
}
