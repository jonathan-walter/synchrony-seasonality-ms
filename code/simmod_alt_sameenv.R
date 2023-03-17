## Main simulation model

simmod_sameenv<-function(tmax, f0, kB, s0, kW, cor.eij,
                      sd.e, dfrac=0, getBt=FALSE){
  
  library(mvtnorm)
  
  if(exp(s0) > 1){
    stop("exp(s0) must be <= 1")
  }
  
  #create environmental noise time series
  sigma<-matrix(0, 2, 2)
  
  sigma[2,1] <- cor.eij*(sd.e^2)
  sigma <- sigma+t(sigma)
  diag(sigma) <- rep(sd.e^2, 2)
  
  env <- rmvnorm(tmax,sigma=sigma)
  colnames(env) <- c("ei","ej")
  
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
    fN <- exp(f0)*exp(-Nt[tt-1,]/kB)*exp(env[tt,])
    sN <- mymin(exp(s0)*exp(-(Nt[tt-1,]*fN)/kW)*exp(env[tt,]))
    Nt[tt,] <- Nt[tt-1,]*fN*sN
    
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


