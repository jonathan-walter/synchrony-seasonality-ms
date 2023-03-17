## Main simulation model

simmod_main_arEnv<-function(tmax, f0, kB, s0, kW, cor.ebij, cor.ewij, cor.ebew,
                      sd.e, arEnv, dfrac=0, getBt=FALSE){
  
  
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
  
  if("try-error" %in% class(try(chol(sigma)))){
    
    if(getBt){
      out <- list(Nt=NA, Bt=NA, env=NA)
    }
    else{
      out <- list(Nt=NA, env=NA)
    }
    
  }
  
  else{
    env <- matrix(NA, 4, tmax)
    for(nn in 1:4){
      env[nn,] <- arima.sim(list(ar=arEnv), tmax)
    }
    for(tt in 1:tmax){
      env[,tt] <- env[,tt] %*% chol(sigma)
    }
    
    env <- t(env)
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
      fN <- exp(f0)*exp(-Nt[tt-1,]/kB)*exp(eb[tt,])
      sN <- mymin(exp(s0)*exp(-(Nt[tt-1,]*fN)/kW)*exp(ew[tt,]))
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
  
  }
  return(out)
}

  
