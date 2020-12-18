#TD EOOMin not quite returning what I expect

############################
#introductory examples
###########################

#Prepare your data
#In this case we'll make some data for Cyprus
library(rCAT)
library(sf)
library(ggplot2)
Cyprus <-   tdwg3[tdwg3$LEVEL3_COD=="CYP",]
pts <- as.data.frame(st_coordinates(st_sample(Cyprus,29)))
#set names to lat long
names(pts) <- c('long','lat')
#simple plot of this
ggplot(data=Cyprus) + geom_sf() + geom_point(data=pts,aes(long,lat))
#project your points into something useful i.e. not lat long
#project and center point are automatically set by the wizard
#you can set the center your self ie propts <- simProjWiz(pts,c(33.25,34),returnV='S')
##TD seems to set lat_ts = 0 it should be the centre or maybe not check original reference, I don't think it matters for equatorial projections, as it's just a shift
thepoints <- simProjWiz(pts,returnV='S')
thepointsSF <- simProjWiz(pts,returnV='SF')
#check project details are stored
attr(thepoints,'crs')
#lets check the EOO, it should be approximately the area of Cyprus's ~ 9,000 km2
eoo(thepoints)
#plot this
eeoploy <- eoo(thepoints,"SF")
cplot <- ggplot(data=Cyprus) + geom_sf() + geom_sf(data=thepointsSF) + geom_sf(data=eeoploy,fill=NA,col='green')
cplot
#EOO rating
ratingEoo(eoo(thepoints))
#get AOO should be ~ 116 km2
aoo(thepoints)
#get AOO as cell polygons
aoopolys <- aoo(thepoints,2000,"SF")
#plot it
eooaooplot <- cplot + geom_sf(data=aoopolys,fill=NA,col='red')
eooaooplot
#get rating from aoo
ratingAoo(aoo(thepoints),FALSE)

#interested in change in AOO or EOO?
1-aooFixedGrid(thepoints[1:22,])/aooFixedGrid(thepoints)
1-eoo(thepoints[1:22,])/eoo(thepoints)
#rating from pop change note needs to be percentage
ratingPop((1-aooFixedGrid(thepoints[1:10,])/aooFixedGrid(thepoints))*100)


#how about the alpha hull? the default is to drop triangles with 2x the mean
aHullMean(thepoints)
#plotting it with 1.5 x mean
apoly = aHullMean(thepoints,1.5,'SF')
cplot + geom_sf(data=apoly,fill=NA,col='blue')
#want to see the detail of the triangles in the Alpha hull?
cplot + geom_sf(data=aHullMean(thepoints,returnV = 'ALL'),fill=NA,col='blue')

#######################################################
#getting down to some detailed analysis
#######################################################
#build some data

#Build some normally distributed point data around the Troodos mountains ~ 10 km diameter
thepoints <- ptsNormal(29,0.1)
#TODO maybe think about make all upper case as SF objects???
#shift to Troodos mountaions
thepoints <- data.frame(long = thepoints$X + 32.8794, lat = thepoints$Y + 34.9220)
#check to see where they are on a map
ggplot(data=Cyprus) + geom_sf() + geom_point(data=thepoints,aes(long,lat))
#project the points
ppts <- simProjWiz(thepoints)
##check crs
attr(ppts,'crs')

####################
#AOO default
####################
#standard AOO
aoo(ppts)
simpleaoopoly <- aoo(ppts,returnV = "SF")
myplot <- ggplot(data=simpleaoopoly) + geom_sf() + geom_point(data=ppts,aes(X,Y))
myplot

####################
#min AOO with shift only
####################
#getting minimum AOO, but only by shifting
minfixaoo <- aooFixedGrid(ppts,returnV = "SF")
myplot + geom_sf(data=minfixaoo, fill="#00FF0044",colour='green')
#interested in all the AOO results to see the range
aoofixall <- aooFixedGrid(ppts,returnV = "ALL")
range(aoofixall$nofcells)
ggplot(data=aoofixall) + geom_histogram(aes(nofcells))

#want to review the maximum grid for fun? - sad I know
maxgrid <- aoofixall[which.max(aoofixall$nofcells),]
#we need to build the cells from the points, rotation and shifts
maxpoly <- buildCellPolys_rxy(ppts,2000,maxgrid$rotation,maxgrid$xshift,maxgrid$yshift)
#plot it
myplot + geom_sf(data=maxpoly, fill="#FF000044",colour='red')
#OK it was not that much fun, but you get the idea
#report
paste("Default grid give", length(simpleaoopoly), "cells. Minimum grid gave", length(minfixaoo), "cells Whilst max grid gave", length(maxpoly),"cells")

############################
#AOO with shift and rotation
############################
#OK lets look at AOO but with rotation as well
#Also this method is more systematic so AOO range and histograms are statistically more meaningful
#just get the area
aooFixedRotation(ppts)
#return the polygon
minpoly <- aooFixedRotation(ppts,returnV="SF")
myplot + geom_sf(data=minpoly, fill="#00FF0044",colour='green') 
#review the range of results
allresults <- aooFixedRotation(ppts,returnV="ALL")
ggplot(data=allresults) + geom_histogram(aes(nofcells))
#you can pull out any results from the dataframe and review
#here we will look at one of the maximum AOO
maxgrid <- allresults[which.max(allresults$nofcells),]
maxpolyrot <- buildCellPolys_rxy(ppts,2000,maxgrid$rotation,maxgrid$xshift,maxgrid$yshift)
myplot + geom_sf(data=maxpolyrot, fill="#FF000044",colour='red')

#let push it harder to find a smaller combination if possible. Increasing iteration by 10x will take a good 10+ seconds to run
minpolyall <- aooFixedRotation(ppts,2000,10^4,returnV = "ALL")
#you may have been able to squeeze out one extra cell, not sure it was worth it!
# plot the min and max
maxgrids <- minpolyall[which.max(minpolyall$nofcells),]
mingrids <- minpolyall[which.min(minpolyall$nofcells),]
maxpolyrots <- buildCellPolys_rxy(ppts,2000,maxgrids$rotation,maxgrids$xshift,maxgrids$yshift)
minpolyrots <- buildCellPolys_rxy(ppts,2000,mingrids$rotation,mingrids$xshift,mingrids$yshift)
ggplot (data=maxpolyrot) + geom_sf(fill="#FF000044",colour='red') + geom_sf(data=minpolyrots, fill="#00ff0044",colour='green') + geom_point(data=ppts,aes(X,Y))

#report
paste("first iterations gave", length(minpoly),". Second with x10 iterations gave", length(minpolyrots),". Max was", length(maxpolyrots))

############################
#AOO using point buffer method
############################
#aoo area using default 2 km2 circles
aooBuf(ppts)
#aoo with 8 km2 area
aooBuf(ppts,bufferradius=1596)
#plot
polysbuf <- aooBuf(ppts,bufferradius=1596,returnV="SF")
myplot <- ggplot(data=polysbuf) + geom_sf() + geom_point(data=ppts,aes(X,Y))
myplot


###############################
#EOO
###############################
eoo(ppts)

###############################
#EOO with point error
###############################

#Looking at EOO with error incorporated
#same data but with a content error of 1km
#report area
eooMin(ppts,defaultRadius = 1000)
eoopolys <- eooMin(ppts,defaultRadius = 1000, returnV="SFA")
ggplot (data=eoopolys$geometry[2]) + geom_sf(fill="#FF000044",colour='red') + geom_sf(data=eoopolys$geometry[1], fill="#00ff0044",colour='green') + geom_point(data=ppts,aes(X,Y))
#report, note getting area from st_area which is in m2
paste("Standard EOO:", eoo(ppts), "Min EOO:",eoopolys$eoo[1],"Max EOO:",  eoopolys$eoo[2])

#and with some variable error ratings between 0 and 3 km
errorv <- runif (nrow(ppts),0,3000)
ppts$error <- errorv
eoopolys <- eooMin(ppts,errorfield = "error",returnV="SFA")

#build the buffers so we can see the error circles on our plots
ps <- st_as_sf(ppts,coords=c('X','Y'))
psbuff <- st_buffer(ps,ps$error)
#need to set the projection
st_crs(psbuff) <- st_crs(attr(ppts,'crs'))
#plot max
ggplot (data=eoopolys$geometry[2]) + geom_sf(fill="#FF000044",colour='red') + geom_sf(data=eoopolys$min, fill="#00ff0044",colour='green') + geom_point(data=ppts,aes(X,Y)) + geom_sf(data=psbuff, fill=NA)
#plot min
ggplot (data=eoopolys$geometry[1]) + geom_sf(fill="#FF000044",colour='red') + geom_sf(data=eoopolys$min, fill="#00ff0044",colour='green') + geom_point(data=ppts,aes(X,Y)) + geom_sf(data=psbuff, fill=NA)

#################################
#looking at subpops / locations
#################################
#get the number of sub-pop/locations using buffer method
subLocBuf(ppts)
#and area
subLocBuf(ppts,returnV="AREA")
#get the sf object and plot
subpop <- subLocBuf(ppts,returnV="SF")
ggplot(data=subpop) + geom_sf() + geom_point(data=ppts,aes(X,Y))
#user defined buffer distance in this case 2 km
subLocBuf(ppts,bufferradius=2000,returnV="S")

#with cell adjacency, use 1/10 width for default
subLocGrid(ppts)
#get area
subLocGrid(ppts,returnV="AREA")
#get the sf object and plot
subpop <- subLocGrid(ppts,returnV="SF")
ggplot(data=subpop) + geom_sf() + geom_point(data=ppts,aes(X,Y))
#Using rook neighborhood
subLocGrid(ppts,neighborhood = "rook")
#with user defined cellwidth
subLocGrid(ppts,neighborhood = "rook",cellwidth = 2000)

#with Alpha hulls - not recogmented
subLocAlpha(ppts)
#get area
subLocAlpha(ppts,returnV = "AREA")
#user defined distance
subLocAlpha(ppts,barrierDis=2000)
#plotting the results
subpop <- subLocAlpha(ppts,returnV = "SF")
ggplot(data=subpop) + geom_sf() + geom_point(data=ppts,aes(X,Y))

#using Rapoports mean propinquity
#number of sub-populations/locations using default distances
#have a look at the EMST first (just)
subLocRapoport(ppts)
#defining your own distances
subLocRapoport(ppts,barrierDis = 2000,bufferDis = 1000)
#area of sub populations/locations using default distances
subLocRapoport(ppts,returnV = "AREA")
#results as simple features for plotting
sfs <- subLocRapoport(ppts,returnV = "SF")
ggplot (data=sfs$tree) + geom_sf(aes(color=barrier)) + geom_sf(data=sfs$buffers,colour='green', fill=NA) + geom_point(data=ppts,aes(X,Y)) 



#############################
#point building scripts
#used for testing and simulations
#############################
#TD need to update points makings scripts to delivery X Y not x y
#random square of points
plot(ptsSquare(500,1000),asp=1)
#random oval of points
plot(ptsOval(500,1000,0.78,0.5),asp=1)
#Doughnut of points
plot(ptsDoughnut(500,1000,0.6,0.5,0.4),asp=1)
# Normally disturbed (but random within) set of points
plot(ptsNormal(500,1000),asp=1)
#rotating points, mainly used for AOO, but maybe useful to others
thepoints <- ptsSquare(200,0.1)
rpoints <- rotateP(thepoints,deg2rad(45))
plot(rpoints,asp=1)

#############################
#Scripts not used yet
#############################
#Distance of the longest Axis Willis et al. 2003 suggest 1/10 as possible grid/scale size
longestAxis(ppts)
dp <- longestAxis(ppts,'P')
plot(x=ppts$X,y=ppts$Y, asp=1)
points(dp,pch=16)

#Minimum enclosing rectangle = simple metric for area
mer(ppts)

#EMST MStree
euMST <- eMST(ppts)
ggplot (data=euMST) + geom_sf(colour="blue")+ geom_point(data=ppts,aes(X,Y)) 
#ID Barriers using rapoport recommendations (mean x 2)
meandis <- mean(euMST$distance)
euMST$barrier <- euMST$distance > meandis * 2
ggplot (data=euMST) + geom_sf(aes(colour=barrier)) + geom_point(data=ppts,aes(X,Y))


#############################
#Batch processing examples
############################
#getting batch working
#set up simple data within Cyprus
lat <- runif (200,32.8294,32.9394)
long <- runif (200,34.8720,34.9720)
mydata <- data.frame(taxa=c('aa','bb','cc','dd',"xx"),lat,long)

#default get a dataframe of results and project all with the same project
resultsdf <- batchCon(mydata$taxa,mydata$long,mydata$lat)
#project each species individually
resultsdf <- batchCon(mydata$taxa,mydata$long,mydata$lat,project2gether = FALSE)
#switch on aooMin
resultsdf <- batchCon(mydata$taxa,mydata$long,mydata$lat,aooMin=TRUE)
#default to return Simple feature objects
resultsf <- batchCon(mydata$taxa,mydata$long,mydata$lat,returnV = "SF")

#plot all the EOO results
ggplot(data=resultsf[resultsf$type=="eoo",]) + geom_sf(fill=NA)
#pull one species and plot
oneSp <- resultsf[resultsf$taxon=="aa",]
ggplot(data=oneSp) + geom_sf(fill=NA)


################################
#Alphas and betas etc
################################
aooFixedGrido(ppts)
aooFixedRotationo(ppts)
#presently below crashes as NE is not working!!!!
countrylist(ppts)


