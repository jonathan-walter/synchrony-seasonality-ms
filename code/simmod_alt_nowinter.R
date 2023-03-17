## Main simulation model

simmod_nowinter<-function(tmax, f0, kB, cor.ebij, cor.ewij, cor.ebew,
                      sd.e, dfrac=0, getBt=FALSE){
  
  library(mvtnorm)
  
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
  
  N0 <- rnorm(2, mean(c(kB, kW)), sd.e)
  
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
    fN <- exp(f0)*exp(-Nt[tt-1,]/kB)*exp(eb[tt,])*exp(ew[tt,])
    #sN <- mymin(exp(s0)*exp(-(Nt[tt-1,]*fN)/kW)*exp(ew[tt,]))
    Nt[tt,] <- Nt[tt-1,]*fN#*sN
    
    if(dfrac > 0){
      Nt[tt,] <- colSums(dmat*Nt[tt,])
    }
    
    if(getBt){
      Bt[tt,] <- Nt[tt-1,]*fN
    }
  }
  if(getBt){
    out <- list(Nt=Nt, Bt=Bt, env=env)
  }
  else{
    out <- list(Nt=Nt, env=env)
  }
  return(out)
}


