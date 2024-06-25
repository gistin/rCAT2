#notes:

Need to update with Steve's country work from https://github.com/gistin/rCAT2/tree/update_countrylist_function https://github.com/gistin/rCAT2/blob/update_countrylist_function/R/countrylist.R https://github.com/gistin/rCAT2/blob/update_countrylist_function/R/get_ne.R https://github.com/gistin/rCAT2/blob/update_countrylist_function/R/get_tdwgl3.R https://github.com/gistin/rCAT2/blob/update_countrylist_function/R/prep.R https://github.com/gistin/rCAT2/blob/update_countrylist_function/R/prep_countries_data.R https://github.com/gistin/rCAT2/blob/update_countrylist_function/R/tdwglist.R


# rCAT: Species Conservation assessment tools in R 
# Version 2.0

## Authors and contributors

Author: Justin Moat
Contributors: Steve Bachman

## 1.0 Introduction

Red Lists are widely used to indicate species at risk of extinction. The tools presented here are scripts that integrate with primary occurrence data, to produce preliminary conservation assessments and to help assessors generate full species conservation assessments.

These tools are to be used with specimen data or sighting data with latitude and longitude co-ordinates and have been primarily designed to work with herbaria data which provides an important source of data relevant for Red List assessments (Willis et al 2003).

The tools follow the IUCN Red Listing guidelines (IUCN 2012, 2016), and primarily produces metrics for Extent of Occurrence (EOO), Area of Occupancy (AOO) and Alpha Hulls. Additionally it adds utilities for the geographic projection of lat, long point data, as well as routinees to build random point dataset for testing.

Any  conservation ratings given here, although calculated based on IUCN categories and criteria, are not equivalent to full IUCN ratings submissible for publication on the IUCN Red List. The majority of measures presented here relate to geometry and species range measures which focus on just one or two aspects of threat considered in the IUCN categories and criteria. Use of these tools assumes the user is familiar with the Red List Criteria (IUCN 2012, IUCN 2016). 

However, the measures calculated do provide a good initial estimate of threat and can be used as a baseline for more detailed assessments that can incorporate population data, further literature, additional GIS analysis and consultation with experts and specialists. 
These tools are primarily produced for quick assessment with the emphasis on speed of processing with map outputs.

### 1.0.1 Brief history

* 2003, first tool produced for ArcView 3.0 as used in Willis et al. 2003 
* 2007, public release of cats.avx = Extension for ESRI ArcView 3.1 (Moat. 2007) 
* 2011, release of GeoCAT online tools geocat.kew.org (Bachman et al 2011) 
* 2017, release of [rCAT 1.5](https://cran.r-project.org/src/contrib/Archive/rCAT/)
* 2020, minor update to [rCAT 1.6](https://cran.r-project.org/web/packages/rCAT/index.html)
* 2020, major update to rCAT 2.0

### 1.0.2	Script development  

Further improvements to the tool will be included in future releases, if you wish to add to this tool, please do fork it or send suggestion to the Authors

### 1.0.3 Installing

Release versions are on cran. Note below is the old version and just for reference
```
install.packages('rCAT')

```
This major update is not on CRAN yes but you can install the development version from GitHub with (you will need the devtools package):

```
install.packages('devtools')
devtools::install_github("gistin/rCAT")
#this will fail until it's opened up
```

### 1.0.4 Quick example

As a quick example, we’ll import and produce metrics for a single species. The data we’ll use is dummy data, but you can easily import a csv.

```
library(rCAT)
#build a square of points of width 0.1 degree
thepoints <- squareOfPs(19,0.1)
#rename field to lat long
names(thepoints) <- c("long","lat")
```
To project the data automatically to an area projection we are going to use the projection wizard. 
```
#this gives you a dataframe of X and Y
thepoints <- simProjWiz(thepoints)
```
To see what projection has been selected
```
attr(thepoints,'crs')
```
This is a base data format for all the analysis
```
#get the EOO area in km2 for these points
eoo(thepoints)
#get the AOO area
aoo(thepoints)
#alpha hull with default distance
aHullMean(thepoints)
```
## Detailed walkthrough
For this detailed walkthrough we will produce an imagery species on the island of Cyprus. To map and visualise the results we will use ggplot2.

```
prepare your data
#In this case we'll make some data for Cyprus
library(sf)
library(ggplot2)
Cyprus <-   tdwg3[tdwg3$LEVEL3_COD=="CYP",]
pts <- as.data.frame(st_coordinates(st_sample(Cyprus,29)))
#set names to lat long
names(pts) <- c('long','lat')
#simple plot of this
ggplot(data=Cyprus) + geom_sf() + geom_point(data=pts,aes(long,lat))
#project your points into something useful ie not lat long
#project and center point are automatically set by the wizard
#you can set the center your self ie propts <- simProjWiz(pts,c(33.25,34),returnV='S')
##TD seems to set lat_ts = 0 it should be the centre or maybe not check orginal reference, I don't think it matters for equatorial projections, as it's just a shift
thepoints <- simProjWiz(pts,,returnV='S')
thepointsSF <- simProjWiz(pts,,returnV='SF')
#check project details are stored
attr(thepoints,'crs')
#lets check the eoo, it should be approximately the area of Cyprus's ~ 9,000 km2
eoo(thepoints)
#plot this
eeoploy <- eoo(thepoints,"SF")
cplot <- ggplot(data=Cyprus) + geom_sf() + geom_sf(data=thepointsSF) + geom_sf(data=eeoploy,fill=NA,col='green')
#EOO rating
eooRating(eoo(thepoints))
#get AOO should be ~ 116 km2
aoo(thepoints)
#get AOO as cell polygons
aoopolys <- aoo(thepoints,2000,"SF")
#plot it
eooaooplot <- cplot + geom_sf(data=aoopolys,fill=NA,col='red')
eooaooplot
#get rating from aoo
aooRating(aoo(thepoints),FALSE)

#interested in change in AOO or EOO?
1-aooFixedGrid(thepoints[1:22,])/aooFixedGrid(thepoints)
1-eoo(thepoints[1:22,])/eoo(thepoints)

#how about the alpha hull? the default is to drop triangles with 2x the mean
aHullMean(thepoints)
#plotting it
apoly = aHullMean(thepoints,1.5,'SF')
cplot + geom_sf(data=apoly,fill=NA,col='blue')
#want to see the detail in the Alpha hull?
cplot + geom_sf(data=aHullMean(thepoints,,'ALL'),fill=NA,col='blue')
```

