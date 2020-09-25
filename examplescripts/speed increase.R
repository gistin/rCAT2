#trying to Increasing speed
library (rCAT)
library(ggplot2)
library(profvis)
library(sf)
library(microbenchmark)

#Build some normally distributed point data around the Troodos mountains ~ 10 km diametre
thepoints <- ptsNormal(29,0.15)
#TODO maybe think about make all upper case as SF objects???
#shift to Troodos mountaions
thepoints <- data.frame(long = thepoints$X + 32.8794, lat = thepoints$Y + 34.9220)
#project the points
ppts <- simProjWiz(thepoints)
##check crs
attr(ppts,'crs')
##speed
#speed test old version
microbenchmark(aooFixedRotation(ppts),times=10L)
microbenchmark(aooFixedGrid (ppts),times=10L)
#1.3-1.4 sec
#NEW
microbenchmark(aooFixedRotationo(ppts,it=1296),times=10L)
microbenchmark(aooFixedGrido (ppts),times=10L)
#0.13
