## investigating covariance between two environmental drivers in two locations
## using PRISM seasonal average temperature

rm(list=ls())

library(raster)
library(spdep)
library(rgdal)
library(ncf)
library(stringr)
library(wsyn)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

summer_tavg <- as.matrix(read.csv("../data/summer_tavg.csv"))
winter_tavg <- as.matrix(read.csv("../data/winter_tavg.csv"))
summer_ppt <- as.matrix(read.csv("../data/summer_ppt.csv"))
winter_ppt <- as.matrix(read.csv("../data/winter_ppt.csv"))

years <- 1990:2009

summer_tavg.cln <- summer_tavg[complete.cases(summer_tavg),] #remove locations with missing data
winter_tavg.cln <- winter_tavg[complete.cases(winter_tavg),]

summer_tavg.cln <- cleandat(summer_tavg.cln, years, clev=3)$cdat #clean data by detrending and standardizing variance
winter_tavg.cln <- cleandat(winter_tavg.cln, years, clev=3)$cdat


summer_ppt.cln <- summer_ppt[complete.cases(summer_ppt),]
winter_ppt.cln <- winter_ppt[complete.cases(winter_ppt),]

summer_ppt.cln <- cleandat(summer_ppt.cln, years, clev=3)$cdat
winter_ppt.cln <- cleandat(winter_ppt.cln, years, clev=3)$cdat


nn <- nrow(summer_tavg.cln)


## analyses


## summer and winter temperature
cor_b1w2<-matrix(NA, nn, nn)

for(ii in 1:nrow(summer_tavg.cln)){
  for(jj in 1:ii){
    #correlate breeding season with next winter
    cor_b1w2[ii,jj]<-cor(summer_tavg.cln[ii,1:19], winter_tavg.cln[jj,2:20])
  }
}

hist(cor_b1w2)
summary(c(cor_b1w2))

quantile(cor_b1w2[lower.tri(cor_b1w2)], c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975), na.rm=T)

quantile(diag(cor_b1w2),c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975), na.rm=T)




## summer temperature and winter precipitation
cor_b1w2<-matrix(NA, nn, nn)

for(ii in 1:nrow(summer_tavg.cln)){
  for(jj in 1:ii){
    #correlate breeding season with next winter
    cor_b1w2[ii,jj]<-cor(summer_tavg.cln[ii,1:19], winter_ppt.cln[jj,2:20])
  }
}

hist(cor_b1w2)
summary(c(cor_b1w2))

quantile(cor_b1w2[lower.tri(cor_b1w2)], c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975), na.rm=T)

quantile(diag(cor_b1w2),c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975), na.rm=T)

