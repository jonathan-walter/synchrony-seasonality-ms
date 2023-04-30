## Theoretical analysis of impacts of seasonality in population processes 
## and environmental drivers on spatial synchrony and metapopulation variability

## setup ------------------------------------------------------------------------------------------

#rm(list=ls())

#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(rootSolve)
library(parallel)
library(lhs)
library(dplyr)
library(viridis)

## source and define helper functions -------------------------------------------------------------
source("simmod_main.R")
source("simmod_alt_nowinter.R")
source("simmod_alt_sameenv.R")
source("simmod_main_deterministic.R")
source("simmod_main_arEnv.R")

getNtcor<-function(inlist){ #extract synchrony in Nt from model output
  out<-rep(NA, length(inlist))
  for(ii in 1:length(inlist)){
    if(length(inlist[[ii]]$Nt)==1){
      out[ii] <- NA
    }
    else{
      tmp<-inlist[[ii]]$Nt[-c(1:burn),]
      out[ii]<-cor(tmp)[2,1]
    }
  }
  return(out)
}

getBtcor<-function(inlist){ #extract synchrony in Bt from model output
  out<-rep(NA, length(inlist))
  for(ii in 1:length(inlist)){
    tmp<-inlist[[ii]]$Bt[-c(1:burn),]
    out[ii]<-cor(tmp)[2,1]
  }
  return(out)
}


#helper function for checking regime
regime.check <- function(inlist){
  out <- rep(NA, length(inlist))
  for(ii in 1:length(out)){
    tmp <- round(inlist[[ii]]$Nt, digits=3)
    if(tmp[length(tmp)] <= 0.1){
      out[ii] <- "check"
    }
    else if(any(tmp > tmp[length(tmp)])){
      out[ii] <- "overcompensatory"
    }
    else{
      out[ii] <- "undercompensatory"
    }
  }
  return(out)
}

# compute analytical solution to linearized model
analytical.solution<-function(f0, kB, s0, kW, cor.ebij, cor.ewij, cor.ebew, sd.e){
  
  #find equilibrium value of N; this is when rate = 0
  rate <- function(N){
    exp(f0)*exp(-N/kB) * exp(s0)*exp(-(N*exp(f0)*exp(-N/kB))/kW) - 1
  }
  
  Eq <- uniroot.all(rate, c(1, max(c(kB,kW))*1.5))
  
  #test stability of equilibrium
  eig <- vector()
  for (i in 1:length(Eq)){
    eig[i] <- sign (gradient(rate, Eq[i]))
  }
  
  if(!any(eig==-1)){
    stop("No stable equilibria between N=1 and N=1.5*max(kB, kW)")
  }
  
  Neq <- Eq[eig==-1] #take the stable equilibrium if multiple
  
  
  # #compute derivatives at N
  g <- expression(Neq*exp(f0)*exp(-Neq/kB)*exp(eB)*
                    exp(s0)*exp(-(Neq*exp(f0)*exp(-Neq/kB)*exp(eB))/kW)*exp(eW))
  dgdeb <- D(g, 'eB')
  dgdew <- D(g, 'eW')
  eB <- eW <- 0
  PB <- as.numeric(eval(dgdeb))
  PW <- as.numeric(eval(dgdew))
  
  #convert between correlation/sd and covariance/variance
  cov.ebij <- cor.ebij*(sd.e^2)
  cov.ewij <- cor.ewij*(sd.e^2)
  cov.ebew <- cor.ebew*(sd.e^2)
  var.e <- sd.e^2
  
  # this is equation 26 from Dan's math document
  return((PB^2*cov.ebij + PW^2*cov.ewij + 2*PB*PW*cov.ebew)/
           (PB^2*var.e + PW^2*var.e + 2*PB*PW*cov.ebew))
  
}


## Simulation -- sensitivity of synchrony to model parameters --------------------------------------------------------

# quick test of model
# test <- simmod_main(tmax=500, f0=2.45, kB=200, s0=-1, kW=200, cor.ebij=0.25, cor.ewij=0.25, cor.ebew=0.1,
#                     sd.e=0.1, dfrac=0)
# plot(test$Nt[-c(1:100),1], type="l")

nn<-300
set.seed(666)
hypercube<-optimumLHS(n=nn, k=9)
#saveRDS(hypercube, "hypercube.rds") #this takes a long time, so optionally save it and load it using the next line
#hypercube <- readRDS("/Users/jonathanwalter/Desktop/synchrony-seasonality code/hypercube.rds")

#maps hypercube onto distributions for parameter choices
f0=qunif(hypercube[,1],0.3,2.45) #density-independent growth rate
kB=qunif(hypercube[,2],10,200) #breeding season density dependence parameter
s0=qunif(hypercube[,3],-0.9, 0)*f0 #density-independent overwintering survival rate
kW=kB*qunif(hypercube[,4],0.1,1) #overwintering density dependence paramter
cor.ebij=qunif(hypercube[,5],0,1) #correlation in breeding season environments
cor.ewij=qunif(hypercube[,6],0,1) #correlation in overwintering season environments
sd.e=qunif(hypercube[,7],0,0.1) #environmental variability
dfrac=qunif(hypercube[,8],0,0.2) #dispersal praction
cor.ebew=qunif(hypercube[,9],-0.5,0.5) #cross-season environmental correlation
cor.eij <- cor.ebij #cross-season environmental correlation for same-environment formulation


tmax=10000 #set max simulation timesteps
burn=1000 #set length of burn-in period

ncores=detectCores()-4

cl<-makeCluster(ncores) #note: this works on a Mac and may not work on a Windows system.
main.sample<-mcmapply(simmod_main, tmax, f0, kB, s0, kW, cor.ebij, cor.ewij, cor.ebew,
                      sd.e, dfrac, getBt=TRUE, SIMPLIFY=FALSE)

nowint.sample<-mcmapply(simmod_nowinter, tmax, f0, kB, cor.ebij, cor.ewij, cor.ebew,
                        sd.e, dfrac, getBt=TRUE, SIMPLIFY=FALSE)

sameenv.sample<-mcmapply(simmod_sameenv, tmax, f0, kB, s0, kW, cor.eij, sd.e, dfrac, 
                         getBt=TRUE, SIMPLIFY=FALSE)
main.determ<-mcmapply(simmod_det, tmax=1000, f0, kB, s0, kW, getBt=FALSE, SIMPLIFY=FALSE)

stopCluster(cl)



table(regime.check(main.determ))



#make data frames for analysis  
sampres.main <- data.frame(rhoN=getNtcor(main.sample), 
                           #scale predictors to standardize effect sizes
                           f0=scale(f0), kB=scale(kB), s0=scale(s0), kW=scale(qunif(hypercube[,4],0.1,1)), cor.ebij=scale(cor.ebij),
                           cor.ewij=scale(cor.ewij), cor.ebew=scale(cor.ebew), sd.e=scale(sd.e), dfrac=scale(dfrac), 
                           regime=regime.check(main.determ))

sampres.nowint <- data.frame(rhoN=getNtcor(nowint.sample), 
                             #scale predictors to standardize effect sizes
                             f0=scale(f0), kB=scale(kB), cor.ebij=scale(cor.ebij),
                             cor.ewij=scale(cor.ewij), cor.ebew=scale(cor.ebew), sd.e=scale(sd.e), dfrac=scale(dfrac),
                             regime=regime.check(main.determ))

sampres.sameenv <- data.frame(rhoN=getNtcor(sameenv.sample), 
                              #scale predictors to standardize effect sizes
                              f0=scale(f0), kB=scale(kB), s0=scale(s0), kW=scale(qunif(hypercube[,4],0.1,1)), cor.ebij=scale(cor.eij),
                              sd.e=scale(sd.e), dfrac=scale(dfrac),
                              regime=regime.check(main.determ))


#regression analysis -- undercompensatory
lm.main.rhoN.uc<-lm(rhoN ~ f0 + kB + s0 + kW + cor.ebij + cor.ewij + cor.ebew + sd.e + dfrac,
                 data=sampres.main[sampres.main$regime=="undercompensatory" & sampres.main$cor.ebij > sampres.main$cor.ewij,], na.action=na.exclude)
summary(lm.main.rhoN.uc)
#hist(resid(lm.main.rhoN))
#plot(lm.main.rhoN.uc)
modsumm.main.uc<-summary(lm.main.rhoN.uc)

lm.nowint.rhoN.uc<-lm(rhoN ~ f0 + kB + cor.ebij + cor.ewij + sd.e + dfrac + cor.ebew,
                   data=sampres.nowint[sampres.nowint$regime=="undercompensatory" & sampres.nowint$cor.ebij > sampres.nowint$cor.ewij,], na.action=na.exclude)
summary(lm.nowint.rhoN.uc)
#hist(resid(lm.nowint.rhoN))
#plot(lm.nowint.rhoN.uc)
modsumm.nowint.uc<-summary(lm.nowint.rhoN.uc)

lm.sameenv.rhoN.uc<-lm(rhoN ~ f0 + kB + s0 + kW + cor.ebij + sd.e + dfrac,
                    data=sampres.sameenv[sampres.sameenv$regime=='undercompensatory',], na.action=na.exclude)
summary(lm.sameenv.rhoN.uc)
#hist(resid(lm.sameenv.rhoN.uc))
#plot(lm.sameenv.rhoN.uc)
modsumm.sameenv.uc<-summary(lm.sameenv.rhoN.uc)

#compile results into a figure
effects.main.uc.eb <- modsumm.main.uc$coefficients[-1,1]
effects.main.uc.eb <- as.data.frame(effects.main.uc.eb)
effects.main.uc.eb$param <- row.names(effects.main.uc.eb)

effects.sameenv.uc.eb <- modsumm.sameenv.uc$coefficients[-1,1]
effects.sameenv.uc.eb <- as.data.frame(effects.sameenv.uc.eb)
effects.sameenv.uc.eb$param <- row.names(effects.sameenv.uc.eb)

effects.nowint.uc.eb <- modsumm.nowint.uc$coefficients[-1,1]
effects.nowint.uc.eb <- as.data.frame(effects.nowint.uc.eb)
effects.nowint.uc.eb$param <- row.names(effects.nowint.uc.eb)

effects.comb.uc.eb <- left_join(effects.main.uc.eb, effects.nowint.uc.eb)
effects.comb.uc.eb <- left_join(effects.comb.uc.eb, effects.sameenv.uc.eb)


#regression analysis -- overcompensatory
lm.main.rhoN.oc<-lm(rhoN ~ f0 + kB + s0 + kW + cor.ebij + cor.ewij + cor.ebew + sd.e + dfrac, 
                    data=sampres.main[sampres.main$regime=="overcompensatory" & sampres.main$cor.ebij > sampres.main$cor.ewij,], na.action=na.exclude)
summary(lm.main.rhoN.oc)
#plot(lm.main.rhoN.oc)
modsumm.main.oc<-summary(lm.main.rhoN.oc)

lm.nowint.rhoN.oc<-lm(rhoN ~ f0 + kB + cor.ebij + cor.ewij + sd.e + dfrac + cor.ebew,
                      data=sampres.nowint[sampres.nowint$regime=="overcompensatory" & sampres.nowint$cor.ebij > sampres.nowint$cor.ewij,], na.action=na.exclude)
summary(lm.nowint.rhoN.oc)
#plot(lm.nowint.rhoN.oc)
modsumm.nowint.oc<-summary(lm.nowint.rhoN.oc)

lm.sameenv.rhoN.oc<-lm(rhoN ~ f0 + kB + s0 + kW + cor.ebij + sd.e + dfrac,
                       data=sampres.sameenv[sampres.sameenv$regime=='overcompensatory',], na.action=na.exclude)
summary(lm.sameenv.rhoN.oc)
#plot(lm.sameenv.rhoN.oc)
modsumm.sameenv.oc<-summary(lm.sameenv.rhoN.oc)

#compile results
effects.main.oc.eb <- modsumm.main.oc$coefficients[-1,1]
effects.main.oc.eb <- as.data.frame(effects.main.oc.eb)
effects.main.oc.eb$param <- row.names(effects.main.oc.eb)

effects.sameenv.oc.eb <- modsumm.sameenv.oc$coefficients[-1,1]
effects.sameenv.oc.eb <- as.data.frame(effects.sameenv.oc.eb)
effects.sameenv.oc.eb$param <- row.names(effects.sameenv.oc.eb)

effects.nowint.oc.eb <- modsumm.nowint.oc$coefficients[-1,1]
effects.nowint.oc.eb <- as.data.frame(effects.nowint.oc.eb)
effects.nowint.oc.eb$param <- row.names(effects.nowint.oc.eb)

effects.comb.oc.eb <- left_join(effects.main.oc.eb, effects.nowint.oc.eb)
effects.comb.oc.eb <- left_join(effects.comb.oc.eb, effects.sameenv.oc.eb)



## cor.ew > cor.eb

#regression analysis -- undercompensatory
lm.main.rhoN.uc<-lm(rhoN ~ f0 + kB + s0 + kW + cor.ebij + cor.ewij + cor.ebew + sd.e + dfrac,
                    data=sampres.main[sampres.main$regime=="undercompensatory" & sampres.main$cor.ebij < sampres.main$cor.ewij,], na.action=na.exclude)
modsumm.main.uc<-summary(lm.main.rhoN.uc)

lm.nowint.rhoN.uc<-lm(rhoN ~ f0 + kB + cor.ebij + cor.ewij + sd.e + dfrac + cor.ebew,
                      data=sampres.nowint[sampres.nowint$regime=="undercompensatory" & sampres.nowint$cor.ebij < sampres.nowint$cor.ewij,], na.action=na.exclude)
modsumm.nowint.uc<-summary(lm.nowint.rhoN.uc)

lm.sameenv.rhoN.uc<-lm(rhoN ~ f0 + kB + s0 + kW + cor.ebij + sd.e + dfrac,
                       data=sampres.sameenv[sampres.sameenv$regime=='undercompensatory',], na.action=na.exclude)
modsumm.sameenv.uc<-summary(lm.sameenv.rhoN.uc)

#compile results into a figure
effects.main.uc.ew <- modsumm.main.uc$coefficients[-1,1]
effects.main.uc.ew <- as.data.frame(effects.main.uc.ew)
effects.main.uc.ew$param <- row.names(effects.main.uc.ew)

effects.sameenv.uc.ew <- modsumm.sameenv.uc$coefficients[-1,1]
effects.sameenv.uc.ew <- as.data.frame(effects.sameenv.uc.ew)
effects.sameenv.uc.ew$param <- row.names(effects.sameenv.uc.ew)

effects.nowint.uc.ew <- modsumm.nowint.uc$coefficients[-1,1]
effects.nowint.uc.ew <- as.data.frame(effects.nowint.uc.ew)
effects.nowint.uc.ew$param <- row.names(effects.nowint.uc.ew)

effects.comb.uc.ew <- left_join(effects.main.uc.ew, effects.nowint.uc.ew)
effects.comb.uc.ew <- left_join(effects.comb.uc.ew, effects.sameenv.uc.ew)


#regression analysis -- overcompensatory
lm.main.rhoN.oc<-lm(rhoN ~ f0 + kB + s0 + kW + cor.ebij + cor.ewij + cor.ebew + sd.e + dfrac, 
                    data=sampres.main[sampres.main$regime=="overcompensatory" & sampres.main$cor.ebij < sampres.main$cor.ewij,], na.action=na.exclude)
modsumm.main.oc<-summary(lm.main.rhoN.oc)

lm.nowint.rhoN.oc<-lm(rhoN ~ f0 + kB + cor.ebij + cor.ewij + sd.e + dfrac + cor.ebew,
                      data=sampres.nowint[sampres.nowint$regime=="overcompensatory" & sampres.nowint$cor.ebij < sampres.nowint$cor.ewij,], na.action=na.exclude)
modsumm.nowint.oc<-summary(lm.nowint.rhoN.oc)

lm.sameenv.rhoN.oc<-lm(rhoN ~ f0 + kB + s0 + kW + cor.ebij + sd.e + dfrac,
                       data=sampres.sameenv[sampres.sameenv$regime=='overcompensatory',], na.action=na.exclude)
modsumm.sameenv.oc<-summary(lm.sameenv.rhoN.oc)

#compile results into a figure
effects.main.oc.ew <- modsumm.main.oc$coefficients[-1,1]
effects.main.oc.ew <- as.data.frame(effects.main.oc.ew)
effects.main.oc.ew$param <- row.names(effects.main.oc.ew)

effects.sameenv.oc.ew <- modsumm.sameenv.oc$coefficients[-1,1]
effects.sameenv.oc.ew <- as.data.frame(effects.sameenv.oc.ew)
effects.sameenv.oc.ew$param <- row.names(effects.sameenv.oc.ew)

effects.nowint.oc.ew <- modsumm.nowint.oc$coefficients[-1,1]
effects.nowint.oc.ew <- as.data.frame(effects.nowint.oc.ew)
effects.nowint.oc.ew$param <- row.names(effects.nowint.oc.ew)

effects.comb.oc.ew <- left_join(effects.main.oc.ew, effects.nowint.oc.ew)
effects.comb.oc.ew <- left_join(effects.comb.oc.ew, effects.sameenv.oc.ew)


## Make a unified figure showing parameter sensitivity--this is Figure 3 in the main text

param.names <- c(expression(italic('f')[0])
                 ,expression(italic('k'['B']))
                 ,expression(italic(hat('s'))[0])
                 ,expression(italic('k'['W']))
                 ,expression(paste('cor(',italic(epsilon['B,i']),',', italic(epsilon['B,j']),')'))
                 ,expression(paste('cor(',italic(epsilon['W,i']),',', italic(epsilon['W,j']),')'))
                 ,expression(paste('cor(',italic(epsilon['B']),',', italic(epsilon['W']),')'))
                 ,expression(paste('sd(', italic(epsilon), ')'))
                 ,expression(italic('d'))
                 )


b <- barplot(t(as.matrix(effects.comb.uc.eb[,-2])), beside=T, plot=FALSE, space=c(0,0.9))
pal <- c("grey20","grey80","grey90")

png("../fig3_simsensitivity_barplot.png", units="in", res=300, width=6.5, height=6)

par(mfcol=c(2,2), mar=c(5.6,2.3,1.6,0.5), mgp=c(2.8,0.8,0), oma=c(0,1.1,0,0))

barplot(t(as.matrix(effects.comb.uc.eb[,-2])), beside=T, names.arg=param.names,
        legend.text=c("main", "no overwintering", "same environ."), ylim=c(-0.15,0.3),
        args.legend=list(x="topleft",bty="n",cex=0.9), las=2, ylab="", col=pal, space=c(0,0.9))
mtext(expression(paste('Undercompensatory, ','cor(',italic(epsilon['B,i']),',', italic(epsilon['B,j']),')','>cor(',italic(epsilon['W,i']),',', italic(epsilon['W,j']),')'))
      , line=0.5, cex=2/3)
mtext("a)", at=0, cex=2/3, line=0.5)
text(x=c(b[2,3:4],b[3,6:7]), y=0, "na", cex=0.4, pos=3, offset=0)

barplot(t(as.matrix(effects.comb.oc.eb[,-2])), beside=T, names.arg=param.names, 
        las=2, ylim=c(-0.15,0.3), ylab="", col=pal, space=c(0,0.9))
mtext(expression(paste('Overcompensatory, ','cor(',italic(epsilon['B,i']),',', italic(epsilon['B,j']),')','>cor(',italic(epsilon['W,i']),',', italic(epsilon['W,j']),')'))
      , line=0.5, cex=2/3)
mtext("b)", at=0, cex=2/3, line=0.5)
text(x=c(b[2,3:4],b[3,6:7]), y=0, "na", cex=0.4, pos=3, offset=0)


barplot(t(as.matrix(effects.comb.uc.ew[,-2])), beside=T, names.arg=param.names,
        ylim=c(-0.15,0.3), space=c(0,0.9),
        args.legend=list(x="topleft",bty="n",cex=0.9), las=2, ylab="", col=pal)
mtext(expression(paste('Undercompensatory, ','cor(',italic(epsilon['W,i']),',', italic(epsilon['W,j']),')','>cor(',italic(epsilon['B,i']),',', italic(epsilon['B,j']),')'))
      , line=0.5, cex=2/3)
mtext("c)", at=0, cex=2/3, line=0.5)
text(x=c(b[2,3:4],b[3,6:7]), y=0, "na", cex=0.4, pos=3, offset=0)

barplot(t(as.matrix(effects.comb.oc.ew[,-2])), beside=T, names.arg=param.names, 
        las=2, ylim=c(-0.15,0.3), ylab="", col=pal, space=c(0,0.9))
mtext(expression(paste('Overcompensatory, ','cor(',italic(epsilon['W,i']),',', italic(epsilon['W,j']),')','>cor(',italic(epsilon['B,i']),',', italic(epsilon['B,j']),')'))
      , line=0.5, cex=2/3)
text(x=c(b[2,3:4],b[3,6:7]), y=0, "na", cex=0.4, pos=3, offset=0)
mtext("d)", at=0, cex=2/3, line=0.5)

mtext("Effect size (standardized regression coefficient)", side=2, outer=TRUE, cex=0.9)

dev.off()




## Simulations -- vary environmental synchrony, parameter set 1 -----------------------------------

scentxt<-"set1" #undercompensatory

#parameter set 1
rho<-seq(0,1,by=0.05) #range of synchrony values to test

tmax = 2000 #simulation length
burn = 1000 #burn-in period
f0 = 1 #density-independent growth rate
kB = 100 #breeding season density dependence parameter
s0 = -0.1 #density-independent overwintering survival rate
kW = 80 #overwintering season density dependence parameter
cor.ebij = rep(expand.grid(rho,rho)[,1],each=75) #breeding season environmental correlations (replicated)
cor.ewij = rep(expand.grid(rho,rho)[,2],each=75) #overwintering season environmental correlations (replicated)
cor.ebew = 0 #cross-season environmental correlation
sd.e = 0.1 #environmental variability
dfrac = 0 #dispersal fraction
cor.eij = cor.ebij #environmental correlations for same-environment model



ncores=detectCores()-4

cl<-makeCluster(ncores)

main.varrhos.uc<-mcmapply(simmod_main, tmax, f0, kB, s0, kW, cor.ebij, cor.ewij, cor.ebew,
                       sd.e, dfrac, getBt=TRUE, SIMPLIFY=FALSE)
nowint.varrhos.uc<-mcmapply(simmod_nowinter, tmax, f0, kB, cor.ebij, cor.ewij, cor.ebew,
                         sd.e, dfrac, getBt=FALSE, SIMPLIFY=FALSE)
sameenv.varrhos.uc<-mcmapply(simmod_sameenv, tmax, f0, kB, s0, kW, cor.eij, sd.e, dfrac, 
                          getBt=TRUE, SIMPLIFY=FALSE)

stopCluster(cl)



#overcompensatory

#parameter set 1
rho<-seq(0,1,by=0.05)

tmax = 2000
burn = 1000
f0 = 2
kB = 100
s0 = -0.1
kW = 80
cor.ebij = rep(expand.grid(rho,rho)[,1],each=75)
cor.ewij = rep(expand.grid(rho,rho)[,2],each=75)
cor.ebew = 0
sd.e = 0.1
dfrac = 0
cor.eij = cor.ebij



ncores=detectCores()-4

cl<-makeCluster(ncores)

main.varrhos.oc<-mcmapply(simmod_main, tmax, f0, kB, s0, kW, cor.ebij, cor.ewij, cor.ebew,
                          sd.e, dfrac, getBt=TRUE, SIMPLIFY=FALSE)
nowint.varrhos.oc<-mcmapply(simmod_nowinter, tmax, f0, kB, cor.ebij, cor.ewij, cor.ebew,
                            sd.e, dfrac, getBt=FALSE, SIMPLIFY=FALSE)
sameenv.varrhos.oc<-mcmapply(simmod_sameenv, tmax, f0, kB, s0, kW, cor.eij, sd.e, dfrac, 
                             getBt=TRUE, SIMPLIFY=FALSE)

stopCluster(cl)



Ntcor.main.uc<-data.frame(cor.ebij=cor.ebij, cor.ewij=cor.ewij, Ntcor=getNtcor(main.varrhos.uc))
Ntcor.main.uc<-aggregate(Ntcor.main.uc, by=list(cor.ebij, cor.ewij), FUN=mean, na.rm=T)
Ntcor.main.uc<-matrix(Ntcor.main.uc$Ntcor, nrow=length(rho), ncol=length(rho))

Ntcor.nowint.uc<-data.frame(cor.ebij=cor.ebij, cor.ewij=cor.ewij, Ntcor=getNtcor(nowint.varrhos.uc))
Ntcor.nowint.uc<-aggregate(Ntcor.nowint.uc, by=list(cor.ebij, cor.ewij), FUN=mean, na.rm=T)
Ntcor.nowint.uc<-matrix(Ntcor.nowint.uc$Ntcor, nrow=length(rho), ncol=length(rho))

Ntcor.sameenv.uc<-data.frame(cor.ebij=cor.ebij, Ntcor=getNtcor(sameenv.varrhos.uc))
Ntcor.sameenv.uc<-aggregate(Ntcor.sameenv.uc, by=list(cor.ebij), FUN=mean, na.rm=T)$Ntcor


Ntcor.main.oc<-data.frame(cor.ebij=cor.ebij, cor.ewij=cor.ewij, Ntcor=getNtcor(main.varrhos.oc))
Ntcor.main.oc<-aggregate(Ntcor.main.oc, by=list(cor.ebij, cor.ewij), FUN=mean, na.rm=T)
Ntcor.main.oc<-matrix(Ntcor.main.oc$Ntcor, nrow=length(rho), ncol=length(rho))

Ntcor.nowint.oc<-data.frame(cor.ebij=cor.ebij, cor.ewij=cor.ewij, Ntcor=getNtcor(nowint.varrhos.oc))
Ntcor.nowint.oc<-aggregate(Ntcor.nowint.oc, by=list(cor.ebij, cor.ewij), FUN=mean, na.rm=T)
Ntcor.nowint.oc<-matrix(Ntcor.nowint.oc$Ntcor, nrow=length(rho), ncol=length(rho))

Ntcor.sameenv.oc<-data.frame(cor.ebij=cor.ebij, Ntcor=getNtcor(sameenv.varrhos.oc))
Ntcor.sameenv.oc<-aggregate(Ntcor.sameenv.oc, by=list(cor.ebij), FUN=mean, na.rm=T)$Ntcor


Btcor.main.uc<-data.frame(cor.ebij=cor.ebij, cor.ewij=cor.ewij, Btcor=getBtcor(main.varrhos.uc))
Btcor.main.uc<-aggregate(Btcor.main.uc, by=list(cor.ebij, cor.ewij), FUN=mean, na.rm=T)
Btcor.main.uc<-matrix(Btcor.main.uc$Btcor, nrow=length(rho), ncol=length(rho))
diff.NtBtcor.main.uc <- Ntcor.main.uc - Btcor.main.uc

Btcor.main.oc<-data.frame(cor.ebij=cor.ebij, cor.ewij=cor.ewij, Btcor=getBtcor(main.varrhos.oc))
Btcor.main.oc<-aggregate(Btcor.main.oc, by=list(cor.ebij, cor.ewij), FUN=mean, na.rm=T)
Btcor.main.oc<-matrix(Btcor.main.oc$Btcor, nrow=length(rho), ncol=length(rho))
diff.NtBtcor.main.oc <- Ntcor.main.oc - Btcor.main.oc

pal<-colorRampPalette(colors=c("red","white","blue"))


#make figure comparing main model to alternate

png(paste0("../figS2_compare_alternates_",scentxt,".png"), 
    units="in", res=300, width=6.5, height=4.71)

par(mfrow=c(2,3), mar=c(1.8,3.1,2.5,0), mgp=c(1.7,0.5,0), tcl=-0.3, oma=c(1.25,0,0,1))

#main model
image(rho, rho, Ntcor.main.uc, xlab="", ylab="Overwintering synchrony", asp=1,
      main="", col=pal(50), zlim=c(-1,1))
contour(rho, rho, Ntcor.main.uc, add=T)
mtext("Main model",3,line=0.1,cex=0.67)
mtext("a)", at=0.01, cex=0.67, line=0.1)

#alternate - no overwintering
image(rho, rho, Ntcor.nowint.uc, xlab="", ylab="Overwintering synchrony", asp=1,
      main="", col=pal(50), zlim=c(-1,1))
contour(rho, rho, Ntcor.nowint.uc, add=T)
mtext("Alt: no overwintering",3,line=0.1,cex=0.67)
mtext("Undercompensatory",3,line=1.3, cex=0.67)
mtext("b)", at=0.01, cex=0.67, line=0.1)


#alternate - same environment
plot(rho, Ntcor.sameenv.uc, pch=16, xlab="", ylab="Population synchrony")
mtext("Alt: same environment",3,line=0.1,cex=0.67)
mtext("c)", at=0.01, cex=0.67, line=0.1)


#main model
image(rho, rho, Ntcor.main.oc, xlab="", ylab="Overwintering synchrony", asp=1,
      main="", col=pal(50), zlim=c(-1,1))
contour(rho, rho, Ntcor.main.oc, add=T)
mtext("Main model",3,line=0.1,cex=0.67)
mtext("d)", at=0.01, cex=0.67, line=0.1)

#alternate - no overwintering
image(rho, rho, Ntcor.nowint.oc, xlab="", ylab="Overwintering synchrony", asp=1,
      main="", col=pal(50), zlim=c(-1,1))
contour(rho, rho, Ntcor.nowint.oc, add=T)
mtext("Alt: no overwintering",3,line=0.1,cex=0.67)
mtext("Overcompensatory",3,line=1.3, cex=0.67)
mtext("e)", at=0.01, cex=0.67, line=0.1)

#alternate - same environment
plot(rho, Ntcor.sameenv.oc, pch=16, xlab="", ylab="Population synchrony")
mtext("Alt: same environment",3,line=0.1,cex=0.67)
mtext("f)", at=0.01, cex=0.67, line=0.1)


mtext("Spatial synchrony of breeding season environment",1,outer=T,cex=0.67,line=0.1)
#mtext("Spatial synchrony of overwintering season environment",2,outer=T,cex=0.8,line=0)

dev.off()



## Make figure comparing syncrhony in different seasons
png(paste0("../figS3_compare_seasons_",scentxt,".png"), units="in", res=300, width=5.98, height=4.6)

par(mfrow=c(2,3), mar=c(1.8,1.8,2.5,0), mgp=c(2.7,0.5,0), tcl=-0.3, oma=c(1.1,1,0,1))
image(rho, rho, Ntcor.main.uc, xlab="", ylab="", asp=1,
      col=pal(50), zlim=c(-1,1))
contour(rho, rho, Ntcor.main.uc, add=T)
mtext("cor(Nt)", cex=0.67, line=0.1)
mtext("a)", at=0.01, cex=0.67, line=0.1)

image(rho, rho, Btcor.main.uc, xlab="", ylab="", asp=1,
      col=pal(50), zlim=c(-1,1))
contour(rho, rho, Btcor.main.uc, add=T)
mtext("b)", at=0.01, cex=0.67, line=0.1)
mtext("cor(Bt)", cex=0.67, line=0.1)
mtext("Undercompensatory", line=1.3, cex=0.67)

image(rho, rho, diff.NtBtcor.main.uc, xlab="", ylab="", asp=1,
      col=pal(50), zlim=c(-.6,.6))
contour(rho, rho, diff.NtBtcor.main.uc, add=T)
mtext("cor(Nt)-cor(Bt)", cex=0.67, line=0.1)
mtext("c)", at=0.01, cex=0.67, line=0.1)


image(rho, rho, Ntcor.main.oc, xlab="", ylab="", asp=1,
      col=pal(50), zlim=c(-1,1))
contour(rho, rho, Ntcor.main.oc, add=T)
mtext("cor(Nt)", cex=0.67, line=0.1)
mtext("d)", at=0.01, cex=0.67, line=0.1)

image(rho, rho, Btcor.main.oc, xlab="", ylab="", asp=1,
      col=pal(50), zlim=c(-1,1))
contour(rho, rho, Btcor.main.oc, add=T)
mtext("e)", at=0.01, cex=0.67, line=0.1)
mtext("cor(Bt)", cex=0.67, line=0.1)
mtext("Overcompensatory", line=1.3, cex=0.67)

image(rho, rho, diff.NtBtcor.main.oc, xlab="", ylab="", asp=1,
      col=pal(50), zlim=c(-.6,.6))
contour(rho, rho, diff.NtBtcor.main.oc, add=T)
mtext("cor(Nt)-cor(Bt)", cex=0.67, line=0.1)
mtext("f)", at=0.01, cex=0.67, line=0.1)

mtext("Spatial synchrony of breeding season environment",1,outer=T,cex=0.67)
mtext("Spatial synchrony of overwintering environment",2,outer=T,cex=0.67)

dev.off()




## Simulations -- vary environmental synchrony, parameter set 2 -----------------------------------

scentxt<-"set2"

rho<-seq(0,1,by=0.05)

tmax = 2000
burn = 1000
f0 = 1
kB = 100
s0 = -0.1
kW = 80
cor.ebij = rep(expand.grid(rho,rho)[,1],each=75)
cor.ewij = rep(expand.grid(rho,rho)[,2],each=75)
cor.ebew = 0.2
sd.e = 0.1
dfrac = 0
cor.eij = cor.ebij



ncores=detectCores()-4

cl<-makeCluster(ncores)

main.varrhos.uc<-mcmapply(simmod_main, tmax, f0, kB, s0, kW, cor.ebij, cor.ewij, cor.ebew,
                          sd.e, dfrac, getBt=TRUE, SIMPLIFY=FALSE)
nowint.varrhos.uc<-mcmapply(simmod_nowinter, tmax, f0, kB, cor.ebij, cor.ewij, cor.ebew,
                            sd.e, dfrac, getBt=FALSE, SIMPLIFY=FALSE)
sameenv.varrhos.uc<-mcmapply(simmod_sameenv, tmax, f0, kB, s0, kW, cor.eij, sd.e, dfrac, 
                             getBt=TRUE, SIMPLIFY=FALSE)

stopCluster(cl)



#overcompensatory
rho<-seq(0,1,by=0.05)

tmax = 2000
burn = 1000
f0 = 2
kB = 100
s0 = -0.1
kW = 80
cor.ebij = rep(expand.grid(rho,rho)[,1],each=75)
cor.ewij = rep(expand.grid(rho,rho)[,2],each=75)
cor.ebew = 0.2
sd.e = 0.1
dfrac = 0
cor.eij = cor.ebij



ncores=detectCores()-4

cl<-makeCluster(ncores)

main.varrhos.oc<-mcmapply(simmod_main, tmax, f0, kB, s0, kW, cor.ebij, cor.ewij, cor.ebew,
                          sd.e, dfrac, getBt=TRUE, SIMPLIFY=FALSE)
nowint.varrhos.oc<-mcmapply(simmod_nowinter, tmax, f0, kB, cor.ebij, cor.ewij, cor.ebew,
                            sd.e, dfrac, getBt=FALSE, SIMPLIFY=FALSE)
sameenv.varrhos.oc<-mcmapply(simmod_sameenv, tmax, f0, kB, s0, kW, cor.eij, sd.e, dfrac, 
                             getBt=TRUE, SIMPLIFY=FALSE)

stopCluster(cl)



Ntcor.main.uc<-data.frame(cor.ebij=cor.ebij, cor.ewij=cor.ewij, Ntcor=getNtcor(main.varrhos.uc))
Ntcor.main.uc<-aggregate(Ntcor.main.uc, by=list(cor.ebij, cor.ewij), FUN=mean, na.rm=T)
Ntcor.main.uc<-matrix(Ntcor.main.uc$Ntcor, nrow=length(rho), ncol=length(rho))

Ntcor.nowint.uc<-data.frame(cor.ebij=cor.ebij, cor.ewij=cor.ewij, Ntcor=getNtcor(nowint.varrhos.uc))
Ntcor.nowint.uc<-aggregate(Ntcor.nowint.uc, by=list(cor.ebij, cor.ewij), FUN=mean, na.rm=T)
Ntcor.nowint.uc<-matrix(Ntcor.nowint.uc$Ntcor, nrow=length(rho), ncol=length(rho))

Ntcor.sameenv.uc<-data.frame(cor.ebij=cor.ebij, Ntcor=getNtcor(sameenv.varrhos.uc))
Ntcor.sameenv.uc<-aggregate(Ntcor.sameenv.uc, by=list(cor.ebij), FUN=mean, na.rm=T)$Ntcor


Ntcor.main.oc<-data.frame(cor.ebij=cor.ebij, cor.ewij=cor.ewij, Ntcor=getNtcor(main.varrhos.oc))
Ntcor.main.oc<-aggregate(Ntcor.main.oc, by=list(cor.ebij, cor.ewij), FUN=mean, na.rm=T)
Ntcor.main.oc<-matrix(Ntcor.main.oc$Ntcor, nrow=length(rho), ncol=length(rho))

Ntcor.nowint.oc<-data.frame(cor.ebij=cor.ebij, cor.ewij=cor.ewij, Ntcor=getNtcor(nowint.varrhos.oc))
Ntcor.nowint.oc<-aggregate(Ntcor.nowint.oc, by=list(cor.ebij, cor.ewij), FUN=mean, na.rm=T)
Ntcor.nowint.oc<-matrix(Ntcor.nowint.oc$Ntcor, nrow=length(rho), ncol=length(rho))

Ntcor.sameenv.oc<-data.frame(cor.ebij=cor.ebij, Ntcor=getNtcor(sameenv.varrhos.oc))
Ntcor.sameenv.oc<-aggregate(Ntcor.sameenv.oc, by=list(cor.ebij), FUN=mean, na.rm=T)$Ntcor


Btcor.main.uc<-data.frame(cor.ebij=cor.ebij, cor.ewij=cor.ewij, Btcor=getBtcor(main.varrhos.uc))
Btcor.main.uc<-aggregate(Btcor.main.uc, by=list(cor.ebij, cor.ewij), FUN=mean, na.rm=T)
Btcor.main.uc<-matrix(Btcor.main.uc$Btcor, nrow=length(rho), ncol=length(rho))
diff.NtBtcor.main.uc <- Ntcor.main.uc - Btcor.main.uc

Btcor.main.oc<-data.frame(cor.ebij=cor.ebij, cor.ewij=cor.ewij, Btcor=getBtcor(main.varrhos.oc))
Btcor.main.oc<-aggregate(Btcor.main.oc, by=list(cor.ebij, cor.ewij), FUN=mean, na.rm=T)
Btcor.main.oc<-matrix(Btcor.main.oc$Btcor, nrow=length(rho), ncol=length(rho))
diff.NtBtcor.main.oc <- Ntcor.main.oc - Btcor.main.oc

pal<-colorRampPalette(colors=c("red","white","blue"))

#make figure comparing main model to alternate

png(paste0("../fig4_compare_alternates_",scentxt,".png"), 
    units="in", res=300, width=6.5, height=4.71)

par(mfrow=c(2,3), mar=c(1.8,3.1,2.5,0), mgp=c(1.7,0.5,0), tcl=-0.3, oma=c(1.25,0,0,1))

#main model
image(rho, rho, Ntcor.main.uc, xlab="", ylab="Overwintering synchrony", asp=1,
      main="", col=pal(50), zlim=c(-1,1))
contour(rho, rho, Ntcor.main.uc, add=T)
mtext("Main model",3,line=0.1,cex=0.67)
mtext("a)", at=0.01, cex=0.67, line=0.1)

#alternate - no overwintering
image(rho, rho, Ntcor.nowint.uc, xlab="", ylab="Overwintering synchrony", asp=1,
      main="", col=pal(50), zlim=c(-1,1))
contour(rho, rho, Ntcor.nowint.uc, add=T)
mtext("Alt: no overwintering",3,line=0.1,cex=0.67)
mtext("Undercompensatory",3,line=1.3, cex=0.67)
mtext("b)", at=0.01, cex=0.67, line=0.1)


#alternate - same environment
plot(rho, Ntcor.sameenv.uc, pch=16, xlab="", ylab="Population synchrony")
mtext("Alt: same environment",3,line=0.1,cex=0.67)
mtext("c)", at=0.01, cex=0.67, line=0.1)


#main model
image(rho, rho, Ntcor.main.oc, xlab="", ylab="Overwintering synchrony", asp=1,
      main="", col=pal(50), zlim=c(-1,1))
contour(rho, rho, Ntcor.main.oc, add=T)
mtext("Main model",3,line=0.1,cex=0.67)
mtext("d)", at=0.01, cex=0.67, line=0.1)

#alternate - no overwintering
image(rho, rho, Ntcor.nowint.oc, xlab="", ylab="Overwintering synchrony", asp=1,
      main="", col=pal(50), zlim=c(-1,1))
contour(rho, rho, Ntcor.nowint.oc, add=T)
mtext("Alt: no overwintering",3,line=0.1,cex=0.67)
mtext("Overcompensatory",3,line=1.3, cex=0.67)
mtext("e)", at=0.01, cex=0.67, line=0.1)

#alternate - same environment
plot(rho, Ntcor.sameenv.oc, pch=16, xlab="", ylab="Population synchrony")
mtext("Alt: same environment",3,line=0.1,cex=0.67)
mtext("f)", at=0.01, cex=0.67, line=0.1)


mtext("Spatial synchrony of breeding season environment",1,outer=T,cex=0.67,line=0.1)
#mtext("Spatial synchrony of overwintering season environment",2,outer=T,cex=0.8,line=0)

dev.off()



## Make figure comparing syncrhony in different seasons
png(paste0("../fig5_compare_seasons_",scentxt,".png"), units="in", res=300, 
    width=5.98, height=4.6)

par(mfrow=c(2,3), mar=c(1.8,1.8,2.5,0), mgp=c(2.7,0.5,0), tcl=-0.3, oma=c(1.1,1,0,1))
image(rho, rho, Ntcor.main.uc, xlab="", ylab="", asp=1,
      col=pal(50), zlim=c(-1,1))
contour(rho, rho, Ntcor.main.uc, add=T)
mtext("cor(Nt)", cex=0.67, line=0.1)
mtext("a)", at=0.01, cex=0.67, line=0.1)

image(rho, rho, Btcor.main.uc, xlab="", ylab="", asp=1,
      col=pal(50), zlim=c(-1,1))
contour(rho, rho, Btcor.main.uc, add=T)
mtext("b)", at=0.01, cex=0.67, line=0.1)
mtext("cor(Bt)", cex=0.67, line=0.1)
mtext("Undercompensatory", line=1.3, cex=0.67)

image(rho, rho, diff.NtBtcor.main.uc, xlab="", ylab="", asp=1,
      col=pal(50), zlim=c(-.6,.6))
contour(rho, rho, diff.NtBtcor.main.uc, add=T)
mtext("cor(Nt)-cor(Bt)", cex=0.67, line=0.1)
mtext("c)", at=0.01, cex=0.67, line=0.1)


image(rho, rho, Ntcor.main.oc, xlab="", ylab="", asp=1,
      col=pal(50), zlim=c(-1,1))
contour(rho, rho, Ntcor.main.oc, add=T)
mtext("cor(Nt)", cex=0.67, line=0.1)
mtext("d)", at=0.01, cex=0.67, line=0.1)

image(rho, rho, Btcor.main.oc, xlab="", ylab="", asp=1,
      col=pal(50), zlim=c(-1,1))
contour(rho, rho, Btcor.main.oc, add=T)
mtext("e)", at=0.01, cex=0.67, line=0.1)
mtext("cor(Bt)", cex=0.67, line=0.1)
mtext("Overcompensatory", line=1.3, cex=0.67)

image(rho, rho, diff.NtBtcor.main.oc, xlab="", ylab="", asp=1,
      col=pal(50), zlim=c(-.6,.6))
contour(rho, rho, diff.NtBtcor.main.oc, add=T)
mtext("cor(Nt)-cor(Bt)", cex=0.67, line=0.1)
mtext("f)", at=0.01, cex=0.67, line=0.1)

mtext("Spatial synchrony of breeding season environment",1,outer=T,cex=0.67)
mtext("Spatial synchrony of overwintering environment",2,outer=T,cex=0.67)

dev.off()


## Simulations -- vary environmental synchrony, parameter set 3 -----------------------------------

scentxt<-"set3"


#parameter set 1
rho<-seq(0,1,by=0.05)

tmax = 2000
burn = 1000
f0 = 1
kB = 100
s0 = -0.1
kW = 80
cor.ebij = rep(expand.grid(rho,rho)[,1],each=75)
cor.ewij = rep(expand.grid(rho,rho)[,2],each=75)
cor.ebew = -0.2
sd.e = 0.1
dfrac = 0
cor.eij = cor.ebij



ncores=detectCores()-4

cl<-makeCluster(ncores)

main.varrhos.uc<-mcmapply(simmod_main, tmax, f0, kB, s0, kW, cor.ebij, cor.ewij, cor.ebew,
                          sd.e, dfrac, getBt=TRUE, SIMPLIFY=FALSE)
nowint.varrhos.uc<-mcmapply(simmod_nowinter, tmax, f0, kB, cor.ebij, cor.ewij, cor.ebew,
                            sd.e, dfrac, getBt=FALSE, SIMPLIFY=FALSE)
sameenv.varrhos.uc<-mcmapply(simmod_sameenv, tmax, f0, kB, s0, kW, cor.eij, sd.e, dfrac, 
                             getBt=TRUE, SIMPLIFY=FALSE)

stopCluster(cl)



#overcompensatory

#parameter set 1
rho<-seq(0,1,by=0.05)

tmax = 2000
burn = 1000
f0 = 2
kB = 100
s0 = -0.1
kW = 80
cor.ebij = rep(expand.grid(rho,rho)[,1],each=75)
cor.ewij = rep(expand.grid(rho,rho)[,2],each=75)
cor.ebew = -0.2
sd.e = 0.1
dfrac = 0
cor.eij = cor.ebij



ncores=detectCores()-4

cl<-makeCluster(ncores)

main.varrhos.oc<-mcmapply(simmod_main, tmax, f0, kB, s0, kW, cor.ebij, cor.ewij, cor.ebew,
                          sd.e, dfrac, getBt=TRUE, SIMPLIFY=FALSE)
nowint.varrhos.oc<-mcmapply(simmod_nowinter, tmax, f0, kB, cor.ebij, cor.ewij, cor.ebew,
                            sd.e, dfrac, getBt=FALSE, SIMPLIFY=FALSE)
sameenv.varrhos.oc<-mcmapply(simmod_sameenv, tmax, f0, kB, s0, kW, cor.eij, sd.e, dfrac, 
                             getBt=TRUE, SIMPLIFY=FALSE)

stopCluster(cl)



Ntcor.main.uc<-data.frame(cor.ebij=cor.ebij, cor.ewij=cor.ewij, Ntcor=getNtcor(main.varrhos.uc))
Ntcor.main.uc<-aggregate(Ntcor.main.uc, by=list(cor.ebij, cor.ewij), FUN=mean, na.rm=T)
Ntcor.main.uc<-matrix(Ntcor.main.uc$Ntcor, nrow=length(rho), ncol=length(rho))

Ntcor.nowint.uc<-data.frame(cor.ebij=cor.ebij, cor.ewij=cor.ewij, Ntcor=getNtcor(nowint.varrhos.uc))
Ntcor.nowint.uc<-aggregate(Ntcor.nowint.uc, by=list(cor.ebij, cor.ewij), FUN=mean, na.rm=T)
Ntcor.nowint.uc<-matrix(Ntcor.nowint.uc$Ntcor, nrow=length(rho), ncol=length(rho))

Ntcor.sameenv.uc<-data.frame(cor.ebij=cor.ebij, Ntcor=getNtcor(sameenv.varrhos.uc))
Ntcor.sameenv.uc<-aggregate(Ntcor.sameenv.uc, by=list(cor.ebij), FUN=mean, na.rm=T)$Ntcor


Ntcor.main.oc<-data.frame(cor.ebij=cor.ebij, cor.ewij=cor.ewij, Ntcor=getNtcor(main.varrhos.oc))
Ntcor.main.oc<-aggregate(Ntcor.main.oc, by=list(cor.ebij, cor.ewij), FUN=mean, na.rm=T)
Ntcor.main.oc<-matrix(Ntcor.main.oc$Ntcor, nrow=length(rho), ncol=length(rho))

Ntcor.nowint.oc<-data.frame(cor.ebij=cor.ebij, cor.ewij=cor.ewij, Ntcor=getNtcor(nowint.varrhos.oc))
Ntcor.nowint.oc<-aggregate(Ntcor.nowint.oc, by=list(cor.ebij, cor.ewij), FUN=mean, na.rm=T)
Ntcor.nowint.oc<-matrix(Ntcor.nowint.oc$Ntcor, nrow=length(rho), ncol=length(rho))

Ntcor.sameenv.oc<-data.frame(cor.ebij=cor.ebij, Ntcor=getNtcor(sameenv.varrhos.oc))
Ntcor.sameenv.oc<-aggregate(Ntcor.sameenv.oc, by=list(cor.ebij), FUN=mean, na.rm=T)$Ntcor



Btcor.main.uc<-data.frame(cor.ebij=cor.ebij, cor.ewij=cor.ewij, Btcor=getBtcor(main.varrhos.uc))
Btcor.main.uc<-aggregate(Btcor.main.uc, by=list(cor.ebij, cor.ewij), FUN=mean, na.rm=T)
Btcor.main.uc<-matrix(Btcor.main.uc$Btcor, nrow=length(rho), ncol=length(rho))
diff.NtBtcor.main.uc <- Ntcor.main.uc - Btcor.main.uc

Btcor.main.oc<-data.frame(cor.ebij=cor.ebij, cor.ewij=cor.ewij, Btcor=getBtcor(main.varrhos.oc))
Btcor.main.oc<-aggregate(Btcor.main.oc, by=list(cor.ebij, cor.ewij), FUN=mean, na.rm=T)
Btcor.main.oc<-matrix(Btcor.main.oc$Btcor, nrow=length(rho), ncol=length(rho))
diff.NtBtcor.main.oc <- Ntcor.main.oc - Btcor.main.oc

pal<-colorRampPalette(colors=c("red","white","blue"))


#make figure comparing main model to alternate

png(paste0("../figS4_compare_alternates_",scentxt,".png"), 
    units="in", res=300, width=6.5, height=4.71)

par(mfrow=c(2,3), mar=c(1.8,3.1,2.5,0), mgp=c(1.7,0.5,0), tcl=-0.3, oma=c(1.25,0,0,1))

#main model
image(rho, rho, Ntcor.main.uc, xlab="", ylab="Overwintering synchrony", asp=1,
      main="", col=pal(50), zlim=c(-1,1))
contour(rho, rho, Ntcor.main.uc, add=T)
mtext("Main model",3,line=0.1,cex=0.67)
mtext("a)", at=0.01, cex=0.67, line=0.1)

#alternate - no overwintering
image(rho, rho, Ntcor.nowint.uc, xlab="", ylab="Overwintering synchrony", asp=1,
      main="", col=pal(50), zlim=c(-1,1))
contour(rho, rho, Ntcor.nowint.uc, add=T)
mtext("Alt: no overwintering",3,line=0.1,cex=0.67)
mtext("Undercompensatory",3,line=1.3, cex=0.67)
mtext("b)", at=0.01, cex=0.67, line=0.1)


#alternate - same environment
plot(rho, Ntcor.sameenv.uc, pch=16, xlab="", ylab="Population synchrony")
mtext("Alt: same environment",3,line=0.1,cex=0.67)
mtext("c)", at=0.01, cex=0.67, line=0.1)


#main model
image(rho, rho, Ntcor.main.oc, xlab="", ylab="Overwintering synchrony", asp=1,
      main="", col=pal(50), zlim=c(-1,1))
contour(rho, rho, Ntcor.main.oc, add=T)
mtext("Main model",3,line=0.1,cex=0.67)
mtext("d)", at=0.01, cex=0.67, line=0.1)

#alternate - no overwintering
image(rho, rho, Ntcor.nowint.oc, xlab="", ylab="Overwintering synchrony", asp=1,
      main="", col=pal(50), zlim=c(-1,1))
contour(rho, rho, Ntcor.nowint.oc, add=T)
mtext("Alt: no overwintering",3,line=0.1,cex=0.67)
mtext("Overcompensatory",3,line=1.3, cex=0.67)
mtext("e)", at=0.01, cex=0.67, line=0.1)

#alternate - same environment
plot(rho, Ntcor.sameenv.oc, pch=16, xlab="", ylab="Population synchrony")
mtext("Alt: same environment",3,line=0.1,cex=0.67)
mtext("f)", at=0.01, cex=0.67, line=0.1)


mtext("Spatial synchrony of breeding season environment",1,outer=T,cex=0.67,line=0.1)
#mtext("Spatial synchrony of overwintering season environment",2,outer=T,cex=0.8,line=0)

dev.off()



## Make figure comparing syncrhony in different seasons
png(paste0("../figS5_compare_seasons_",scentxt,".png"), units="in", res=300, width=5.99, height=4.6)

par(mfrow=c(2,3), mar=c(1.8,1.8,2.5,0), mgp=c(2.7,0.5,0), tcl=-0.3, oma=c(1.1,1,0,1))
image(rho, rho, Ntcor.main.uc, xlab="", ylab="", asp=1,
      col=pal(50), zlim=c(-1,1))
contour(rho, rho, Ntcor.main.uc, add=T)
mtext("cor(Nt)", cex=0.67, line=0.1)
mtext("a)", at=0.01, cex=0.67, line=0.1)

image(rho, rho, Btcor.main.uc, xlab="", ylab="", asp=1,
      col=pal(50), zlim=c(-1,1))
contour(rho, rho, Btcor.main.uc, add=T)
mtext("b)", at=0.01, cex=0.67, line=0.1)
mtext("cor(Bt)", cex=0.67, line=0.1)
mtext("Undercompensatory", line=1.3, cex=0.67)

image(rho, rho, diff.NtBtcor.main.uc, xlab="", ylab="", asp=1,
      col=pal(50), zlim=c(-.7,.7))
contour(rho, rho, diff.NtBtcor.main.uc, add=T)
mtext("cor(Nt)-cor(Bt)", cex=0.67, line=0.1)
mtext("c)", at=0.01, cex=0.67, line=0.1)


image(rho, rho, Ntcor.main.oc, xlab="", ylab="", asp=1,
      col=pal(50), zlim=c(-1,1))
contour(rho, rho, Ntcor.main.oc, add=T)
mtext("cor(Nt)", cex=0.67, line=0.1)
mtext("d)", at=0.01, cex=0.67, line=0.1)

image(rho, rho, Btcor.main.oc, xlab="", ylab="", asp=1,
      col=pal(50), zlim=c(-1,1))
contour(rho, rho, Btcor.main.oc, add=T)
mtext("e)", at=0.01, cex=0.67, line=0.1)
mtext("cor(Bt)", cex=0.67, line=0.1)
mtext("Overcompensatory", line=1.3, cex=0.67)

image(rho, rho, diff.NtBtcor.main.oc, xlab="", ylab="", asp=1,
      col=pal(50), zlim=c(-.4,.4))
contour(rho, rho, diff.NtBtcor.main.oc, add=T)
mtext("cor(Nt)-cor(Bt)", cex=0.67, line=0.1)
mtext("f)", at=0.01, cex=0.67, line=0.1)

mtext("Spatial synchrony of breeding season environment",1,outer=T,cex=0.67)
mtext("Spatial synchrony of overwintering environment",2,outer=T,cex=0.67)

dev.off()




## Simulations -- vary environmental synchrony, parameter set 4 -----------------------------------

scentxt<-"set4"

rho<-seq(0,1,by=0.05)

tmax = 2000
burn = 1000
f0 = 1
kB = 100
s0 = -0.1
kW = 80
cor.ebij = rep(expand.grid(rho,rho)[,1],each=75)
cor.ewij = rep(expand.grid(rho,rho)[,2],each=75)
cor.ebew = 0.2
sd.e = 0.1
dfrac = 0.1
cor.eij = cor.ebij



ncores=detectCores()-4

cl<-makeCluster(ncores)

main.varrhos.uc<-mcmapply(simmod_main, tmax, f0, kB, s0, kW, cor.ebij, cor.ewij, cor.ebew,
                          sd.e, dfrac, getBt=TRUE, SIMPLIFY=FALSE)
nowint.varrhos.uc<-mcmapply(simmod_nowinter, tmax, f0, kB, cor.ebij, cor.ewij, cor.ebew,
                            sd.e, dfrac, getBt=FALSE, SIMPLIFY=FALSE)
sameenv.varrhos.uc<-mcmapply(simmod_sameenv, tmax, f0, kB, s0, kW, cor.eij, sd.e, dfrac, 
                             getBt=TRUE, SIMPLIFY=FALSE)

stopCluster(cl)



#overcompensatory
rho<-seq(0,1,by=0.05)

tmax = 2000
burn = 1000
f0 = 2
kB = 100
s0 = -0.1
kW = 80
cor.ebij = rep(expand.grid(rho,rho)[,1],each=75)
cor.ewij = rep(expand.grid(rho,rho)[,2],each=75)
cor.ebew = 0.2
sd.e = 0.1
dfrac = 0.1
cor.eij = cor.ebij



ncores=detectCores()-4

cl<-makeCluster(ncores)

main.varrhos.oc<-mcmapply(simmod_main, tmax, f0, kB, s0, kW, cor.ebij, cor.ewij, cor.ebew,
                          sd.e, dfrac, getBt=TRUE, SIMPLIFY=FALSE)
nowint.varrhos.oc<-mcmapply(simmod_nowinter, tmax, f0, kB, cor.ebij, cor.ewij, cor.ebew,
                            sd.e, dfrac, getBt=FALSE, SIMPLIFY=FALSE)
sameenv.varrhos.oc<-mcmapply(simmod_sameenv, tmax, f0, kB, s0, kW, cor.eij, sd.e, dfrac, 
                             getBt=TRUE, SIMPLIFY=FALSE)

stopCluster(cl)



Ntcor.main.uc<-data.frame(cor.ebij=cor.ebij, cor.ewij=cor.ewij, Ntcor=getNtcor(main.varrhos.uc))
Ntcor.main.uc<-aggregate(Ntcor.main.uc, by=list(cor.ebij, cor.ewij), FUN=mean, na.rm=T)
Ntcor.main.uc<-matrix(Ntcor.main.uc$Ntcor, nrow=length(rho), ncol=length(rho))

Ntcor.nowint.uc<-data.frame(cor.ebij=cor.ebij, cor.ewij=cor.ewij, Ntcor=getNtcor(nowint.varrhos.uc))
Ntcor.nowint.uc<-aggregate(Ntcor.nowint.uc, by=list(cor.ebij, cor.ewij), FUN=mean, na.rm=T)
Ntcor.nowint.uc<-matrix(Ntcor.nowint.uc$Ntcor, nrow=length(rho), ncol=length(rho))

Ntcor.sameenv.uc<-data.frame(cor.ebij=cor.ebij, Ntcor=getNtcor(sameenv.varrhos.uc))
Ntcor.sameenv.uc<-aggregate(Ntcor.sameenv.uc, by=list(cor.ebij), FUN=mean, na.rm=T)$Ntcor


Ntcor.main.oc<-data.frame(cor.ebij=cor.ebij, cor.ewij=cor.ewij, Ntcor=getNtcor(main.varrhos.oc))
Ntcor.main.oc<-aggregate(Ntcor.main.oc, by=list(cor.ebij, cor.ewij), FUN=mean, na.rm=T)
Ntcor.main.oc<-matrix(Ntcor.main.oc$Ntcor, nrow=length(rho), ncol=length(rho))

Ntcor.nowint.oc<-data.frame(cor.ebij=cor.ebij, cor.ewij=cor.ewij, Ntcor=getNtcor(nowint.varrhos.oc))
Ntcor.nowint.oc<-aggregate(Ntcor.nowint.oc, by=list(cor.ebij, cor.ewij), FUN=mean, na.rm=T)
Ntcor.nowint.oc<-matrix(Ntcor.nowint.oc$Ntcor, nrow=length(rho), ncol=length(rho))

Ntcor.sameenv.oc<-data.frame(cor.ebij=cor.ebij, Ntcor=getNtcor(sameenv.varrhos.oc))
Ntcor.sameenv.oc<-aggregate(Ntcor.sameenv.oc, by=list(cor.ebij), FUN=mean, na.rm=T)$Ntcor

Btcor.main.uc<-data.frame(cor.ebij=cor.ebij, cor.ewij=cor.ewij, Btcor=getBtcor(main.varrhos.uc))
Btcor.main.uc<-aggregate(Btcor.main.uc, by=list(cor.ebij, cor.ewij), FUN=mean, na.rm=T)
Btcor.main.uc<-matrix(Btcor.main.uc$Btcor, nrow=length(rho), ncol=length(rho))
diff.NtBtcor.main.uc <- Ntcor.main.uc - Btcor.main.uc

Btcor.main.oc<-data.frame(cor.ebij=cor.ebij, cor.ewij=cor.ewij, Btcor=getBtcor(main.varrhos.oc))
Btcor.main.oc<-aggregate(Btcor.main.oc, by=list(cor.ebij, cor.ewij), FUN=mean, na.rm=T)
Btcor.main.oc<-matrix(Btcor.main.oc$Btcor, nrow=length(rho), ncol=length(rho))
diff.NtBtcor.main.oc <- Ntcor.main.oc - Btcor.main.oc

pal<-colorRampPalette(colors=c("red","white","blue"))

#make figure comparing main model to alternate

png(paste0("../figS6_compare_alternates_",scentxt,".png"), 
    units="in", res=300, width=6.5, height=4.71)

par(mfrow=c(2,3), mar=c(1.8,3.1,2.5,0), mgp=c(1.7,0.5,0), tcl=-0.3, oma=c(1.25,0,0,1))

#main model
image(rho, rho, Ntcor.main.uc, xlab="", ylab="Overwintering synchrony", asp=1,
      main="", col=pal(50), zlim=c(-1,1))
contour(rho, rho, Ntcor.main.uc, add=T)
mtext("Main model",3,line=0.1,cex=0.67)
mtext("a)", at=0.01, cex=0.67, line=0.1)

#alternate - no overwintering
image(rho, rho, Ntcor.nowint.uc, xlab="", ylab="Overwintering synchrony", asp=1,
      main="", col=pal(50), zlim=c(-1,1))
contour(rho, rho, Ntcor.nowint.uc, add=T)
mtext("Alt: no overwintering",3,line=0.1,cex=0.67)
mtext("Undercompensatory",3,line=1.3, cex=0.67)
mtext("b)", at=0.01, cex=0.67, line=0.1)


#alternate - same environment
plot(rho, Ntcor.sameenv.uc, pch=16, xlab="", ylab="Population synchrony")
mtext("Alt: same environment",3,line=0.1,cex=0.67)
mtext("c)", at=0.01, cex=0.67, line=0.1)


#main model
image(rho, rho, Ntcor.main.oc, xlab="", ylab="Overwintering synchrony", asp=1,
      main="", col=pal(50), zlim=c(-1,1))
contour(rho, rho, Ntcor.main.oc, add=T)
mtext("Main model",3,line=0.1,cex=0.67)
mtext("d)", at=0.01, cex=0.67, line=0.1)

#alternate - no overwintering
image(rho, rho, Ntcor.nowint.oc, xlab="", ylab="Overwintering synchrony", asp=1,
      main="", col=pal(50), zlim=c(-1,1))
contour(rho, rho, Ntcor.nowint.oc, add=T)
mtext("Alt: no overwintering",3,line=0.1,cex=0.67)
mtext("Overcompensatory",3,line=1.3, cex=0.67)
mtext("e)", at=0.01, cex=0.67, line=0.1)

#alternate - same environment
plot(rho, Ntcor.sameenv.oc, pch=16, xlab="", ylab="Population synchrony")
mtext("Alt: same environment",3,line=0.1,cex=0.67)
mtext("f)", at=0.01, cex=0.67, line=0.1)


mtext("Spatial synchrony of breeding season environment",1,outer=T,cex=0.67,line=0.1)
#mtext("Spatial synchrony of overwintering season environment",2,outer=T,cex=0.8,line=0)

dev.off()



## Make figure comparing syncrhony in different seasons
png(paste0("../figS7_compare_seasons_",scentxt,".png"), units="in", res=300, 
    width=5.98, height=4.6)

par(mfrow=c(2,3), mar=c(1.8,1.8,2.5,0), mgp=c(2.7,0.5,0), tcl=-0.3, oma=c(1.1,1,0,1))
image(rho, rho, Ntcor.main.uc, xlab="", ylab="", asp=1,
      col=pal(50), zlim=c(-1,1))
contour(rho, rho, Ntcor.main.uc, add=T)
mtext("cor(Nt)", cex=0.67, line=0.1)
mtext("a)", at=0.01, cex=0.67, line=0.1)

image(rho, rho, Btcor.main.uc, xlab="", ylab="", asp=1,
      col=pal(50), zlim=c(-1,1))
contour(rho, rho, Btcor.main.uc, add=T)
mtext("b)", at=0.01, cex=0.67, line=0.1)
mtext("cor(Bt)", cex=0.67, line=0.1)
mtext("Undercompensatory", line=1.3, cex=0.67)

image(rho, rho, diff.NtBtcor.main.uc, xlab="", ylab="", asp=1,
      col=pal(50), zlim=c(-.6,.6))
contour(rho, rho, diff.NtBtcor.main.uc, add=T)
mtext("cor(Nt)-cor(Bt)", cex=0.67, line=0.1)
mtext("c)", at=0.01, cex=0.67, line=0.1)


image(rho, rho, Ntcor.main.oc, xlab="", ylab="", asp=1,
      col=pal(50), zlim=c(-1,1))
contour(rho, rho, Ntcor.main.oc, add=T)
mtext("cor(Nt)", cex=0.67, line=0.1)
mtext("d)", at=0.01, cex=0.67, line=0.1)

image(rho, rho, Btcor.main.oc, xlab="", ylab="", asp=1,
      col=pal(50), zlim=c(-1,1))
contour(rho, rho, Btcor.main.oc, add=T)
mtext("e)", at=0.01, cex=0.67, line=0.1)
mtext("cor(Bt)", cex=0.67, line=0.1)
mtext("Overcompensatory", line=1.3, cex=0.67)

image(rho, rho, diff.NtBtcor.main.oc, xlab="", ylab="", asp=1,
      col=pal(50), zlim=c(-.6,.6))
contour(rho, rho, diff.NtBtcor.main.oc, add=T)
mtext("cor(Nt)-cor(Bt)", cex=0.67, line=0.1)
mtext("f)", at=0.01, cex=0.67, line=0.1)

mtext("Spatial synchrony of breeding season environment",1,outer=T,cex=0.67)
mtext("Spatial synchrony of overwintering environment",2,outer=T,cex=0.67)

dev.off()




## Simulations -- vary environmental synchrony, parameter set 5 -----------------------------------

scentxt<-"set5"

rho<-seq(0,1,by=0.05)

tmax = 2000
burn = 1000
f0 = 1
kB = 100
s0 = -0.1
kW = 80
cor.ebij = rep(expand.grid(rho,rho)[,1],each=75)
cor.ewij = rep(expand.grid(rho,rho)[,2],each=75)
cor.ebew = 0.2
sd.e = 0.1
dfrac = 0.2
cor.eij = cor.ebij



ncores=detectCores()-4

cl<-makeCluster(ncores)

main.varrhos.uc<-mcmapply(simmod_main, tmax, f0, kB, s0, kW, cor.ebij, cor.ewij, cor.ebew,
                          sd.e, dfrac, getBt=TRUE, SIMPLIFY=FALSE)
nowint.varrhos.uc<-mcmapply(simmod_nowinter, tmax, f0, kB, cor.ebij, cor.ewij, cor.ebew,
                            sd.e, dfrac, getBt=FALSE, SIMPLIFY=FALSE)
sameenv.varrhos.uc<-mcmapply(simmod_sameenv, tmax, f0, kB, s0, kW, cor.eij, sd.e, dfrac, 
                             getBt=TRUE, SIMPLIFY=FALSE)

stopCluster(cl)



#overcompensatory
rho<-seq(0,1,by=0.05)

tmax = 2000
burn = 1000
f0 = 2
kB = 100
s0 = -0.1
kW = 80
cor.ebij = rep(expand.grid(rho,rho)[,1],each=75)
cor.ewij = rep(expand.grid(rho,rho)[,2],each=75)
cor.ebew = 0.2
sd.e = 0.1
dfrac = 0.2
cor.eij = cor.ebij



ncores=detectCores()-4

cl<-makeCluster(ncores)

main.varrhos.oc<-mcmapply(simmod_main, tmax, f0, kB, s0, kW, cor.ebij, cor.ewij, cor.ebew,
                          sd.e, dfrac, getBt=TRUE, SIMPLIFY=FALSE)
nowint.varrhos.oc<-mcmapply(simmod_nowinter, tmax, f0, kB, cor.ebij, cor.ewij, cor.ebew,
                            sd.e, dfrac, getBt=FALSE, SIMPLIFY=FALSE)
sameenv.varrhos.oc<-mcmapply(simmod_sameenv, tmax, f0, kB, s0, kW, cor.eij, sd.e, dfrac, 
                             getBt=TRUE, SIMPLIFY=FALSE)

stopCluster(cl)



Ntcor.main.uc<-data.frame(cor.ebij=cor.ebij, cor.ewij=cor.ewij, Ntcor=getNtcor(main.varrhos.uc))
Ntcor.main.uc<-aggregate(Ntcor.main.uc, by=list(cor.ebij, cor.ewij), FUN=mean, na.rm=T)
Ntcor.main.uc<-matrix(Ntcor.main.uc$Ntcor, nrow=length(rho), ncol=length(rho))

Ntcor.nowint.uc<-data.frame(cor.ebij=cor.ebij, cor.ewij=cor.ewij, Ntcor=getNtcor(nowint.varrhos.uc))
Ntcor.nowint.uc<-aggregate(Ntcor.nowint.uc, by=list(cor.ebij, cor.ewij), FUN=mean, na.rm=T)
Ntcor.nowint.uc<-matrix(Ntcor.nowint.uc$Ntcor, nrow=length(rho), ncol=length(rho))

Ntcor.sameenv.uc<-data.frame(cor.ebij=cor.ebij, Ntcor=getNtcor(sameenv.varrhos.uc))
Ntcor.sameenv.uc<-aggregate(Ntcor.sameenv.uc, by=list(cor.ebij), FUN=mean, na.rm=T)$Ntcor


Ntcor.main.oc<-data.frame(cor.ebij=cor.ebij, cor.ewij=cor.ewij, Ntcor=getNtcor(main.varrhos.oc))
Ntcor.main.oc<-aggregate(Ntcor.main.oc, by=list(cor.ebij, cor.ewij), FUN=mean, na.rm=T)
Ntcor.main.oc<-matrix(Ntcor.main.oc$Ntcor, nrow=length(rho), ncol=length(rho))

Ntcor.nowint.oc<-data.frame(cor.ebij=cor.ebij, cor.ewij=cor.ewij, Ntcor=getNtcor(nowint.varrhos.oc))
Ntcor.nowint.oc<-aggregate(Ntcor.nowint.oc, by=list(cor.ebij, cor.ewij), FUN=mean, na.rm=T)
Ntcor.nowint.oc<-matrix(Ntcor.nowint.oc$Ntcor, nrow=length(rho), ncol=length(rho))

Ntcor.sameenv.oc<-data.frame(cor.ebij=cor.ebij, Ntcor=getNtcor(sameenv.varrhos.oc))
Ntcor.sameenv.oc<-aggregate(Ntcor.sameenv.oc, by=list(cor.ebij), FUN=mean, na.rm=T)$Ntcor

Btcor.main.uc<-data.frame(cor.ebij=cor.ebij, cor.ewij=cor.ewij, Btcor=getBtcor(main.varrhos.uc))
Btcor.main.uc<-aggregate(Btcor.main.uc, by=list(cor.ebij, cor.ewij), FUN=mean, na.rm=T)
Btcor.main.uc<-matrix(Btcor.main.uc$Btcor, nrow=length(rho), ncol=length(rho))
diff.NtBtcor.main.uc <- Ntcor.main.uc - Btcor.main.uc

Btcor.main.oc<-data.frame(cor.ebij=cor.ebij, cor.ewij=cor.ewij, Btcor=getBtcor(main.varrhos.oc))
Btcor.main.oc<-aggregate(Btcor.main.oc, by=list(cor.ebij, cor.ewij), FUN=mean, na.rm=T)
Btcor.main.oc<-matrix(Btcor.main.oc$Btcor, nrow=length(rho), ncol=length(rho))
diff.NtBtcor.main.oc <- Ntcor.main.oc - Btcor.main.oc

pal<-colorRampPalette(colors=c("red","white","blue"))

#make figure comparing main model to alternate

png(paste0("../figS8_compare_alternates_",scentxt,".png"), 
    units="in", res=300, width=6.5, height=4.71)

par(mfrow=c(2,3), mar=c(1.8,3.1,2.5,0), mgp=c(1.7,0.5,0), tcl=-0.3, oma=c(1.25,0,0,1))

#main model
image(rho, rho, Ntcor.main.uc, xlab="", ylab="Overwintering synchrony", asp=1,
      main="", col=pal(50), zlim=c(-1,1))
contour(rho, rho, Ntcor.main.uc, add=T)
mtext("Main model",3,line=0.1,cex=0.67)
mtext("a)", at=0.01, cex=0.67, line=0.1)

#alternate - no overwintering
image(rho, rho, Ntcor.nowint.uc, xlab="", ylab="Overwintering synchrony", asp=1,
      main="", col=pal(50), zlim=c(-1,1))
contour(rho, rho, Ntcor.nowint.uc, add=T)
mtext("Alt: no overwintering",3,line=0.1,cex=0.67)
mtext("Undercompensatory",3,line=1.3, cex=0.67)
mtext("b)", at=0.01, cex=0.67, line=0.1)


#alternate - same environment
plot(rho, Ntcor.sameenv.uc, pch=16, xlab="", ylab="Population synchrony")
mtext("Alt: same environment",3,line=0.1,cex=0.67)
mtext("c)", at=0.01, cex=0.67, line=0.1)


#main model
image(rho, rho, Ntcor.main.oc, xlab="", ylab="Overwintering synchrony", asp=1,
      main="", col=pal(50), zlim=c(-1,1))
contour(rho, rho, Ntcor.main.oc, add=T)
mtext("Main model",3,line=0.1,cex=0.67)
mtext("d)", at=0.01, cex=0.67, line=0.1)

#alternate - no overwintering
image(rho, rho, Ntcor.nowint.oc, xlab="", ylab="Overwintering synchrony", asp=1,
      main="", col=pal(50), zlim=c(-1,1))
contour(rho, rho, Ntcor.nowint.oc, add=T)
mtext("Alt: no overwintering",3,line=0.1,cex=0.67)
mtext("Overcompensatory",3,line=1.3, cex=0.67)
mtext("e)", at=0.01, cex=0.67, line=0.1)

#alternate - same environment
plot(rho, Ntcor.sameenv.oc, pch=16, xlab="", ylab="Population synchrony")
mtext("Alt: same environment",3,line=0.1,cex=0.67)
mtext("f)", at=0.01, cex=0.67, line=0.1)


mtext("Spatial synchrony of breeding season environment",1,outer=T,cex=0.67,line=0.1)
#mtext("Spatial synchrony of overwintering season environment",2,outer=T,cex=0.8,line=0)

dev.off()



## Make figure comparing syncrhony in different seasons
png(paste0("../figS9_compare_seasons_",scentxt,".png"), units="in", res=300, 
    width=5.98, height=4.6)

par(mfrow=c(2,3), mar=c(1.8,1.8,2.5,0), mgp=c(2.7,0.5,0), tcl=-0.3, oma=c(1.1,1,0,1))
image(rho, rho, Ntcor.main.uc, xlab="", ylab="", asp=1,
      col=pal(50), zlim=c(-1,1))
contour(rho, rho, Ntcor.main.uc, add=T)
mtext("cor(Nt)", cex=0.67, line=0.1)
mtext("a)", at=0.01, cex=0.67, line=0.1)

image(rho, rho, Btcor.main.uc, xlab="", ylab="", asp=1,
      col=pal(50), zlim=c(-1,1))
contour(rho, rho, Btcor.main.uc, add=T)
mtext("b)", at=0.01, cex=0.67, line=0.1)
mtext("cor(Bt)", cex=0.67, line=0.1)
mtext("Undercompensatory", line=1.3, cex=0.67)

image(rho, rho, diff.NtBtcor.main.uc, xlab="", ylab="", asp=1,
      col=pal(50), zlim=c(-.6,.6))
contour(rho, rho, diff.NtBtcor.main.uc, add=T)
mtext("cor(Nt)-cor(Bt)", cex=0.67, line=0.1)
mtext("c)", at=0.01, cex=0.67, line=0.1)


image(rho, rho, Ntcor.main.oc, xlab="", ylab="", asp=1,
      col=pal(50), zlim=c(-1,1))
contour(rho, rho, Ntcor.main.oc, add=T)
mtext("cor(Nt)", cex=0.67, line=0.1)
mtext("d)", at=0.01, cex=0.67, line=0.1)

image(rho, rho, Btcor.main.oc, xlab="", ylab="", asp=1,
      col=pal(50), zlim=c(-1,1))
contour(rho, rho, Btcor.main.oc, add=T)
mtext("e)", at=0.01, cex=0.67, line=0.1)
mtext("cor(Bt)", cex=0.67, line=0.1)
mtext("Overcompensatory", line=1.3, cex=0.67)

image(rho, rho, diff.NtBtcor.main.oc, xlab="", ylab="", asp=1,
      col=pal(50), zlim=c(-.6,.6))
contour(rho, rho, diff.NtBtcor.main.oc, add=T)
mtext("cor(Nt)-cor(Bt)", cex=0.67, line=0.1)
mtext("f)", at=0.01, cex=0.67, line=0.1)

mtext("Spatial synchrony of breeding season environment",1,outer=T,cex=0.67)
mtext("Spatial synchrony of overwintering environment",2,outer=T,cex=0.67)

dev.off()





## Explore when the environment has interannual temporal autocorrelation ----------------------------------



rho<-seq(0,1,by=0.05)

tmax = 2000
burn = 1000
f0 = 1
kB = 100
s0 = -0.1
kW = 80
cor.ebij = rep(expand.grid(rho,rho)[,1],each=75)
cor.ewij = rep(expand.grid(rho,rho)[,2],each=75)
cor.ebew = 0.2
sd.e = 0.1
dfrac = 0.2
cor.eij = cor.ebij




ncores=detectCores()-4

cl<-makeCluster(ncores)

main.varrhos.uc.red<-mcmapply(simmod_main_arEnv, tmax, f0, kB, s0, kW, cor.ebij, cor.ewij, cor.ebew,
                          sd.e, arEnv=0.2, dfrac, getBt=TRUE, SIMPLIFY=FALSE)
main.varrhos.uc.blu<-mcmapply(simmod_main_arEnv, tmax, f0, kB, s0, kW, cor.ebij, cor.ewij, cor.ebew,
                              sd.e, arEnv=-0.2, dfrac, getBt=TRUE, SIMPLIFY=FALSE)

stopCluster(cl)



#overcompensatory

f0 = 2


ncores=detectCores()-4

cl<-makeCluster(ncores)

main.varrhos.oc.red<-mcmapply(simmod_main_arEnv, tmax, f0, kB, s0, kW, cor.ebij, cor.ewij, cor.ebew,
                              sd.e, arEnv=0.2, dfrac, getBt=TRUE, SIMPLIFY=FALSE)
main.varrhos.oc.blu<-mcmapply(simmod_main_arEnv, tmax, f0, kB, s0, kW, cor.ebij, cor.ewij, cor.ebew,
                              sd.e, arEnv=-0.2, dfrac, getBt=TRUE, SIMPLIFY=FALSE)

stopCluster(cl)



Ntcor.red.uc<-data.frame(cor.ebij=cor.ebij, cor.ewij=cor.ewij, Ntcor=getNtcor(main.varrhos.uc.red))
Ntcor.red.uc<-aggregate(Ntcor.red.uc, by=list(cor.ebij, cor.ewij), FUN=mean, na.rm=T)
Ntcor.red.uc<-matrix(Ntcor.red.uc$Ntcor, nrow=length(rho), ncol=length(rho))

Ntcor.blu.uc<-data.frame(cor.ebij=cor.ebij, cor.ewij=cor.ewij, Ntcor=getNtcor(main.varrhos.uc.blu))
Ntcor.blu.uc<-aggregate(Ntcor.blu.uc, by=list(cor.ebij, cor.ewij), FUN=mean, na.rm=T)
Ntcor.blu.uc<-matrix(Ntcor.blu.uc$Ntcor, nrow=length(rho), ncol=length(rho))

Ntcor.red.oc<-data.frame(cor.ebij=cor.ebij, cor.ewij=cor.ewij, Ntcor=getNtcor(main.varrhos.oc.red))
Ntcor.red.oc<-aggregate(Ntcor.red.oc, by=list(cor.ebij, cor.ewij), FUN=mean, na.rm=T)
Ntcor.red.oc<-matrix(Ntcor.red.oc$Ntcor, nrow=length(rho), ncol=length(rho))

Ntcor.blu.oc<-data.frame(cor.ebij=cor.ebij, cor.ewij=cor.ewij, Ntcor=getNtcor(main.varrhos.oc.blu))
Ntcor.blu.oc<-aggregate(Ntcor.blu.oc, by=list(cor.ebij, cor.ewij), FUN=mean, na.rm=T)
Ntcor.blu.oc<-matrix(Ntcor.blu.oc$Ntcor, nrow=length(rho), ncol=length(rho))



pal<-colorRampPalette(colors=c("red","white","blue"))



png("./figS10_EnvAr1a.png", units="in", res=300, 
    width=6.5, height=6.41)

par(mfrow=c(2,2), mar=c(3,3.2,1.5,1), mgp=c(1.9,0.5,0), tcl=-0.3, oma=c(0,1.25,0,0))

image(rho, rho, Ntcor.red.uc, asp=1,
      col=pal(50), zlim=c(-1,1), xlab="Breeding synchrony", ylab="Overwintering synchrony")
contour(rho, rho, Ntcor.red.uc, add=T)
mtext("Undercompensatory", cex=0.75)
mtext("a)", at=0.01, cex=0.75, line=0.1)


image(rho, rho, Ntcor.red.oc, asp=1,
      col=pal(50), zlim=c(-1,1), xlab="Breeding synchrony", ylab="Overwintering synchrony")
contour(rho, rho, Ntcor.red.oc, add=T)
mtext("Overcompensatory", cex=0.75)
mtext("b)", at=0.01, cex=0.75, line=0.1)


image(rho, rho, Ntcor.blu.uc, asp=1,
      col=pal(50), zlim=c(-1,1), xlab="Breeding synchrony", ylab="Overwintering synchrony")
contour(rho, rho, Ntcor.blu.uc, add=T)
mtext("Undercompensatory", cex=0.75)
mtext("c)", at=0.01, cex=0.75, line=0.1)

image(rho, rho, Ntcor.blu.oc, asp=1,
      col=pal(50), zlim=c(-1,1), xlab="Breeding synchrony", ylab="Overwintering synchrony")
contour(rho, rho, Ntcor.blu.oc, add=T)
mtext("Overcompensatory", cex=0.75)
mtext("d)", at=0.01, cex=0.75, line=0.1)

mtext("Red noise environment", side=2, at=0.77, line=0.1, outer=TRUE)
mtext("Blue noise environment", side=2, at=0.27, line=0.1, outer=TRUE)

dev.off()




## Explore when the environment has interannual temporal autocorrelation




rho<-seq(0,1,by=0.05)

tmax = 2000
burn = 1000
f0 = 1
kB = 100
s0 = -0.1
kW = 80
cor.ebij = rep(expand.grid(rho,rho)[,1],each=75)
cor.ewij = rep(expand.grid(rho,rho)[,2],each=75)
cor.ebew = 0.2
sd.e = 0.1
dfrac = 0.2
cor.eij = cor.ebij




ncores=detectCores()-4

cl<-makeCluster(ncores)

main.varrhos.uc.red<-mcmapply(simmod_main_arEnv, tmax, f0, kB, s0, kW, cor.ebij, cor.ewij, cor.ebew,
                              sd.e, arEnv=0.5, dfrac, getBt=TRUE, SIMPLIFY=FALSE)
main.varrhos.uc.blu<-mcmapply(simmod_main_arEnv, tmax, f0, kB, s0, kW, cor.ebij, cor.ewij, cor.ebew,
                              sd.e, arEnv=-0.5, dfrac, getBt=TRUE, SIMPLIFY=FALSE)

stopCluster(cl)



#overcompensatory

f0 = 2


ncores=detectCores()-4

cl<-makeCluster(ncores)

main.varrhos.oc.red<-mcmapply(simmod_main_arEnv, tmax, f0, kB, s0, kW, cor.ebij, cor.ewij, cor.ebew,
                              sd.e, arEnv=0.5, dfrac, getBt=TRUE, SIMPLIFY=FALSE)
main.varrhos.oc.blu<-mcmapply(simmod_main_arEnv, tmax, f0, kB, s0, kW, cor.ebij, cor.ewij, cor.ebew,
                              sd.e, arEnv=-0.5, dfrac, getBt=TRUE, SIMPLIFY=FALSE)

stopCluster(cl)



Ntcor.red.uc<-data.frame(cor.ebij=cor.ebij, cor.ewij=cor.ewij, Ntcor=getNtcor(main.varrhos.uc.red))
Ntcor.red.uc<-aggregate(Ntcor.red.uc, by=list(cor.ebij, cor.ewij), FUN=mean, na.rm=T)
Ntcor.red.uc<-matrix(Ntcor.red.uc$Ntcor, nrow=length(rho), ncol=length(rho))

Ntcor.blu.uc<-data.frame(cor.ebij=cor.ebij, cor.ewij=cor.ewij, Ntcor=getNtcor(main.varrhos.uc.blu))
Ntcor.blu.uc<-aggregate(Ntcor.blu.uc, by=list(cor.ebij, cor.ewij), FUN=mean, na.rm=T)
Ntcor.blu.uc<-matrix(Ntcor.blu.uc$Ntcor, nrow=length(rho), ncol=length(rho))

Ntcor.red.oc<-data.frame(cor.ebij=cor.ebij, cor.ewij=cor.ewij, Ntcor=getNtcor(main.varrhos.oc.red))
Ntcor.red.oc<-aggregate(Ntcor.red.oc, by=list(cor.ebij, cor.ewij), FUN=mean, na.rm=T)
Ntcor.red.oc<-matrix(Ntcor.red.oc$Ntcor, nrow=length(rho), ncol=length(rho))

Ntcor.blu.oc<-data.frame(cor.ebij=cor.ebij, cor.ewij=cor.ewij, Ntcor=getNtcor(main.varrhos.oc.blu))
Ntcor.blu.oc<-aggregate(Ntcor.blu.oc, by=list(cor.ebij, cor.ewij), FUN=mean, na.rm=T)
Ntcor.blu.oc<-matrix(Ntcor.blu.oc$Ntcor, nrow=length(rho), ncol=length(rho))



pal<-colorRampPalette(colors=c("red","white","blue"))



png("./figS11_EnvAr1b.png", units="in", res=300, 
    width=6.5, height=6.41)

par(mfrow=c(2,2), mar=c(3,3.2,1.5,1), mgp=c(1.9,0.5,0), tcl=-0.3, oma=c(0,1.25,0,0))

image(rho, rho, Ntcor.red.uc, asp=1,
      col=pal(50), zlim=c(-1,1), xlab="Breeding synchrony", ylab="Overwintering synchrony")
contour(rho, rho, Ntcor.red.uc, add=T)
mtext("Undercompensatory", cex=0.75)
mtext("a)", at=0.01, cex=0.75, line=0.1)


image(rho, rho, Ntcor.red.oc, asp=1,
      col=pal(50), zlim=c(-1,1), xlab="Breeding synchrony", ylab="Overwintering synchrony")
contour(rho, rho, Ntcor.red.oc, add=T)
mtext("Overcompensatory", cex=0.75)
mtext("b)", at=0.01, cex=0.75, line=0.1)


image(rho, rho, Ntcor.blu.uc, asp=1,
      col=pal(50), zlim=c(-1,1), xlab="Breeding synchrony", ylab="Overwintering synchrony")
contour(rho, rho, Ntcor.blu.uc, add=T)
mtext("Undercompensatory", cex=0.75)
mtext("c)", at=0.01, cex=0.75, line=0.1)

image(rho, rho, Ntcor.blu.oc, asp=1,
      col=pal(50), zlim=c(-1,1), xlab="Breeding synchrony", ylab="Overwintering synchrony")
contour(rho, rho, Ntcor.blu.oc, add=T)
mtext("Overcompensatory", cex=0.75)
mtext("d)", at=0.01, cex=0.75, line=0.1)

mtext("Red noise environment", side=2, at=0.77, line=0.1, outer=TRUE)
mtext("Blue noise environment", side=2, at=0.27, line=0.1, outer=TRUE)

dev.off()
