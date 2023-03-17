simmod_det<-function(tmax, f0, kB, s0, kW, getBt=FALSE){
  
  if(exp(s0) > 1){
    stop("exp(s0) must be <= 1")
  }
  
  N0 <- 0.1
  
  Nt <- rep(NA, tmax)
  Nt[1] <- N0
  
  if(getBt){
    Bt <- rep(NA, tmax)
  }
  
  for(tt in 2:tmax){
    fN <- exp(f0)*exp(-Nt[tt-1]/kB)
    sN <- min(exp(s0)*exp(-(Nt[tt-1]*fN)/kW))
    Nt[tt] <- Nt[tt-1]*fN*sN
    
    if(getBt){
      Bt[tt] <- Nt[tt-1,]*fN
    }
  }
  if(getBt){
    out <- list(Nt=Nt, Bt=Bt)
  }
  else{
    out <- list(Nt=Nt)
  }
  return(out)
}