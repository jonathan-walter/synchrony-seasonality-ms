#rm(list=ls())

library(rootSolve)
library(mvtnorm)


## Set up simulation model ------------------------------------------------------------------------

## tmax is the length of the simulation; f0 is the fecundity rate at low density; kB is the
## saturation parameter during the breeding season; s0 is the overwintering survival rate at low density;
## kW is the saturation parameter during the overwintering season; cor.ebij is the spatial synchrony of
## breeding season environment; cor.ewij is the spatial synchrony of the overwintering environment;
## cor.ebew is the correlation between breeding and overwintering environments;
## sd.e is the standard deviation of the environmental noises; N0 are initial population abundances.

simmod_main<-function(tmax, f0, kB, s0, kW, cor.ebij, cor.ewij, cor.ebew,
                      sd.e, dfrac=0, N0=NULL, getBt=FALSE){

  
  if(exp(s0) > 1){
    stop("exp(s0) must be <= 1")
  }
  
  #create environmental noise time series
  sigma<-matrix(0, 4, 4)
  
  sigma[2,1] <- cor.ebij*(sd.e^2)
  sigma[3:4,1] <- cor.ebew*(sd.e^2)
  sigma[3:4,2] <- cor.ebew*(sd.e^2)
  sigma[4,3] <- cor.ewij*(sd.e^2)  
  sigma <- sigma+t(sigma)
  diag(sigma) <- rep(sd.e^2, 4)
  
  env <- rmvnorm(tmax,sigma=sigma)
  colnames(env) <- c("ebi","ebj","ewi","ewj")
  eb <- env[,1:2]
  ew <- env[,3:4]
  
  #if needed, create dispersal matrix
  if(dfrac > 0){
    dmat <- matrix(c(1-dfrac, dfrac, dfrac, 1-dfrac), 2, 2)
    
  }
  
 # if(is.null(N0)){
    N0 <- rnorm(2, min(c(kB, kW)), sd.e)
 # }
  
  Nt <- matrix(NA, tmax, 2)
  Nt[1,] <- N0
  
  if(getBt){
    Bt <- matrix(NA, tmax, 2)
  }
  
  mymin <- function(x){
    out <- rep(NA, length(x))
    for(ii in 1:length(x)){
      out[ii] <- min(c(x[ii]),1)
    }
    return(out)
  }
  
  for(tt in 2:tmax){
    fN <- exp(f0)*exp(-Nt[tt-1,]/kB)*exp(eb[tt,])
    Bt[tt,] <- fN*Nt[tt-1,]
    #sN <- mymin(exp(s0)*exp(-(Nt[tt-1,]*fN)/kW)*exp(ew[tt,]))
    sN <- mymin(exp(s0)*exp(-Bt[tt,]/kW)*exp(ew[tt,]))
    Nt[tt,] <- Nt[tt-1,]*fN*sN
    
    if(dfrac > 0){
      Nt[tt,] <- colSums(dmat*Nt[tt,])
    }
    
 #   if(getBt){
 #     Bt[tt,] <- Nt[tt-1,]*fN
 #  }
  }
  if(getBt){
    out <- list(Nt=Nt, Bt=Bt, env=env)
  }
  else{
    out <- list(Nt=Nt, env=env)
  }
  return(out)
}



##-------------------------------------------------------------------------------------------------
## Set up function for analytical solution

## parameter definitions as above.


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
  g <- expression(Neq*exp(f0)*exp(-Neq/kB)*exp(eB)*exp(s0)*exp(-(Neq*exp(f0)*exp(-Neq/kB)*exp(eB))/kW)*exp(eW))
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
  return((PB^2*cov.ebij + PW^2*cov.ewij + 2*PB*PW*cov.ebew)/(PB^2*var.e + PW^2*var.e + 2*PB*PW*cov.ebew))
  
}

## single scenario run example --------------------------------------------------------------------
# define parameters and run 

###
# parameters across runs
tmax = 2000 #length of sims
burn = 1000 #burn-in period
cor.ebij = seq(0, 1, .05) #breeding season environmental correlation
cor.ewij = seq(0, 1, .05) #overwintering environmental correlation
replicates <- 75  # change to 75 for final run

## multiple scenario run 1 --------------------------------------------------------------------
# define parameters and run 
f0 = 1.05 #density-independent growth rate
kB = 50 #breeding season density dependence
s0 = -0.2 #survival rate
kW = 40 #overwintering season density dependence
cor.ebew = -.1 #cross-season environmental correlation
sd.e = 0.02 #environmental variability
dfrac = 0 #dispersal

results.a <- matrix(NA, nrow=length(cor.ebij), ncol=length(cor.ewij))
simresults.a <- matrix(NA, nrow=length(cor.ebij), ncol=length(cor.ewij))

for(xx in 1:length(cor.ebij)) {
  for(yy in 1:length(cor.ewij)) {
    temp <- rep(NA, replicates)
    for(zz in 1:replicates) {
      linear.test<-simmod_main(tmax, f0, kB, s0, kW, cor.ebij[xx], cor.ewij[yy], cor.ebew, sd.e, dfrac, getBt=TRUE)
      #plot(linear.test$Nt[,1])
      cor_results <- cor(linear.test$Nt[-c(1:burn),])
      temp[zz] <- cor_results[2,1]
    }
    simresults.a[xx,yy] <- mean(temp)
  }
  print(xx)
}


for(xx in 1:length(cor.ebij)) {
  for(yy in 1:length(cor.ewij)) {
    results.a[xx,yy] <- analytical.solution(f0, kB, s0, kW, cor.ebij[xx], cor.ewij[yy], cor.ebew, sd.e)
  }
}

#max(results-simresults)
#min(results-simresults)

## multiple scenario run 2 --------------------------------------------------------------------
# define parameters and run 

f0 = 1.6
kB = 100
s0 = 0
kW = 50
cor.ebew = 0.2
sd.e = 0.01
dfrac = 0

results.b <- matrix(NA, nrow=length(cor.ebij), ncol=length(cor.ewij))
simresults.b <- matrix(NA, nrow=length(cor.ebij), ncol=length(cor.ewij))

for(xx in 1:length(cor.ebij)) {
  for(yy in 1:length(cor.ewij)) {
    temp <- rep(NA, replicates)
    for(zz in 1:replicates) {
      linear.test<-simmod_main(tmax, f0, kB, s0, kW, cor.ebij[xx], cor.ewij[yy], cor.ebew, sd.e, dfrac, getBt=TRUE)
      #plot(linear.test$Nt[,1])
      cor_results <- cor(linear.test$Nt[-c(1:burn),])
      temp[zz] <- cor_results[2,1]
    }
    simresults.b[xx,yy] <- mean(temp)
  }
  print(xx)
}


for(xx in 1:length(cor.ebij)) {
  for(yy in 1:length(cor.ewij)) {
    results.b[xx,yy] <- analytical.solution(f0, kB, s0, kW, cor.ebij[xx], cor.ewij[yy], cor.ebew, sd.e)
  }
}

## multiple scenario run 3 --------------------------------------------------------------------
# define parameters and run 

f0 = 2.2
kB = 100
s0 = -0.1
kW = 80
cor.ebew = 0
sd.e = 0.1
dfrac = 0

results.c <- matrix(NA, nrow=length(cor.ebij), ncol=length(cor.ewij))
simresults.c <- matrix(NA, nrow=length(cor.ebij), ncol=length(cor.ewij))

for(xx in 1:length(cor.ebij)) {
  for(yy in 1:length(cor.ewij)) {
    temp <- rep(NA, replicates)
    for(zz in 1:replicates) {
      linear.test<-simmod_main(tmax, f0, kB, s0, kW, cor.ebij[xx], cor.ewij[yy], cor.ebew, sd.e, dfrac, getBt=TRUE)
      #plot(linear.test$Nt[,1])
      cor_results <- cor(linear.test$Nt[-c(1:burn),])
      temp[zz] <- cor_results[2,1]
    }
    simresults.c[xx,yy] <- mean(temp)
  }
  print(xx)
}


for(xx in 1:length(cor.ebij)) {
  for(yy in 1:length(cor.ewij)) {
    results.c[xx,yy] <- analytical.solution(f0, kB, s0, kW, cor.ebij[xx], cor.ewij[yy], cor.ebew, sd.e)
  }
}


### Plotting

#quartz(height=6, width=4)
pdf("sim_analysis_comparison.pdf", height=6, width=4)

pal<-colorRampPalette(colors=c("red","white","blue"))
par(mfrow=c(3,2), mar=c(1,.5,.5,.5), mgp=c(2.7,0.5,0), tcl=-0.3, oma=c(3,3,3,1))
# sim 1
# analytical
image(cor.ebij, cor.ewij, results.a, zlim=c(-1,1), col=pal(50),
      xlab="", ylab="", xaxt="n", cex=1.25)
contour(cor.ebij, cor.ewij, results.a, add=T)
text(0.0,1.08,"a)", xpd=NA)

#simulation
image(cor.ebij, cor.ewij, simresults.a, zlim=c(-1,1), col=pal(50),
      xlab="", ylab="", cex=1.25, xaxt="n", yaxt="n")
contour(cor.ebij, cor.ewij, simresults.a, add=T)
text(0.0,1.08,"b)", xpd=NA)

# difference
#image(cor.ebij, cor.ewij, results.a-simresults.a, zlim=c(-1,1), col=pal(50),
#      xlab=expression(epsilon[b]), ylab=expression(epsilon[w]), cex=1.25)
#contour(cor.ebij, cor.ewij, results.a-simresults.a, add=T)

# sim 2
# analytical
image(cor.ebij, cor.ewij, results.b, zlim=c(-1,1), col=pal(50),
      xlab="", ylab="", xaxt="n", cex=1.25)
contour(cor.ebij, cor.ewij, results.b, add=T)
text(0.0,1.08,"c)", xpd=NA)

#simulation
image(cor.ebij, cor.ewij, simresults.b, zlim=c(-1,1), col=pal(50),
      xlab="", ylab="", xaxt="n", yaxt="n", cex=1.25)
contour(cor.ebij, cor.ewij, simresults.b, add=T)
text(0.0,1.08,"d)", xpd=NA)

# difference
#image(cor.ebij, cor.ewij, results.b-simresults.b, zlim=c(-1,1), col=pal(50),
#      xlab=expression(epsilon[b]), ylab=expression(epsilon[w]), cex=1.25)
#contour(cor.ebij, cor.ewij, results.b-simresults.b, add=T)

# sim 3
# analytical
image(cor.ebij, cor.ewij, results.c, zlim=c(-1,1), col=pal(50),
      xlab="", ylab="", cex=1.25)
contour(cor.ebij, cor.ewij, results.c, add=T)
text(0.0,1.08,"e)", xpd=NA)

#simulation
image(cor.ebij, cor.ewij, simresults.c, zlim=c(-1,1), col=pal(50),
      xlab="", ylab="", yaxt="n", cex=1.25)
contour(cor.ebij, cor.ewij, simresults.c, add=T)
text(0.0,1.08,"f)", xpd=NA)

# difference
#image(cor.ebij, cor.ewij, results.c-simresults.c, zlim=c(-1,1), col=pal(50),
#      xlab=expression(epsilon[b]), ylab=expression(epsilon[w]), cex=1.25)
#contour(cor.ebij, cor.ewij, results.c-simresults.c, add=T)

mtext(expression(paste("Spatial synchrony of breeding season environment (", epsilon[B], ")")), 
      1, outer=T,cex=0.8, line=1.2)
mtext(expression(paste("Spatial synchrony of overwintering season environment (", epsilon[W], ")")),
      2,outer=T,cex=0.8, line=1.2)
mtext(expression(paste("Analytical")), 3, outer=T, cex=1, line=.3, adj=.2)
mtext(expression(paste("Simuation")), 3, outer=T,cex=1, line=.5, adj=.85)

dev.off()
