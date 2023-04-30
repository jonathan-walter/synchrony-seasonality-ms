#rm(list=ls())

library(rootSolve)
library(mvtnorm)

##-------------------------------------------------------------------------------------------------
## Set up function for analytical solution


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

#---------------------------------------------------------------
# set up deterministic version of model

simmod_det<-function(tmax, f0, kB, s0, kW, getBt=FALSE){
  
  if(exp(s0) > 1){
    stop("exp(s0) must be <= 1")
  }
  
  N0 <- 1
  
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


#helper function for checking regime
regime.check <- function(inlist){
  out <- rep(NA, length(inlist))
  for(ii in 1:length(out)){
    tmp <- inlist[[ii]]
    tmp <- round(tmp, digits=3)
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



## generate figure 1 ---------------------------------------------------------------------------------------------


#### define parameters across runs ####
tmax = 2000 #simulation length
burn = 1000 #burn-in period

#### Plotting set up ####

#quartz(height=6, width=4)
pdf(file="fig1_analytical.pdf", width=3, height=6)
pal<-colorRampPalette(colors=c("red","white","blue"))
par(mfrow=c(3,1), mar=c(3,1.5,1.5,.5), mgp=c(2.7,0.5,0), tcl=-0.3, oma=c(2,3,1,1))

# text(0.0,1.08,"a)", xpd=NA)

#### eB eW heat map ####
f0 = 1.6 #density-independent growth rate
kB = 100 #breeding season density dependence
s0 = 0 #survival rate
kW = 50 #overwintering density dependence
cor.ebew = 0.2 #cross-season environmental correlation
sd.e = 0.01 #environmental variability
cor.ebij = seq(0, 1, .05) #breeding season environmental correlation
cor.ewij = seq(0, 1, .05) #overwintering season environmental correlation
dfrac = 0 #dispersal fraction

results.a <- matrix(NA, nrow=length(cor.ebij), ncol=length(cor.ewij))

for(xx in 1:length(cor.ebij)) {
  for(yy in 1:length(cor.ewij)) {
    results.a[xx,yy] <- analytical.solution(f0, kB, s0, kW, cor.ebij[xx], cor.ewij[yy], cor.ebew, sd.e)
  }
}

image(cor.ebij, cor.ewij, results.a, zlim=c(-1,1), col=pal(50),
      xlab="", ylab="", cex=1.25)
contour(cor.ebij, cor.ewij, results.a, add=T)
text(0.02,1.08,"a)", xpd=NA)
mtext(expression(paste("Synchrony of breeding season environ. (", epsilon[B], ")")), 
      1, outer=F,cex=0.75, line=1.7)
mtext(expression(paste("Synchrony of overwintering")),
      2,outer=F,cex=0.75, line=2.5)
mtext(expression(paste("season environment (", epsilon[W], ")")),
      2,outer=F,cex=0.75, line=1.4)

#### cor(eB,eW) figure ####

f0_1 <- 1.2
prop_1 <- -.4
s0_1 <- f0_1*prop_1

f0_2 <- 1.8
prop_2 <- 0
s0_2 <- f0_2*prop_2

f0_3 <- 2.4
prop_3 <- -.1
s0_3 <- f0_3*prop_3

f0_4 <- .6
prop_4 <- -.7
s0_4 <- f0_4*prop_4


kB = 100
kW = 85
cor.ebij = .3
cor.ewij = .1
cor.ebew = seq(-.5,.5,.05)
sd.e = 0.05
dfrac = 0

results.c <- rep(NA, length=length(cor.ebew))
results.d <- rep(NA, length=length(cor.ebew))
results.e <- rep(NA, length=length(cor.ebew))
results.f <- rep(NA, length=length(cor.ebew))

for(xx in 1:length(cor.ebew)) {
  results.c[xx] <- analytical.solution(f0_1, kB, s0_1, kW, cor.ebij, cor.ewij, cor.ebew[xx], sd.e)
  results.d[xx] <- analytical.solution(f0_2, kB, s0_2, kW, cor.ebij, cor.ewij, cor.ebew[xx], sd.e)
  results.e[xx] <- analytical.solution(f0_3, kB, s0_3, kW, cor.ebij, cor.ewij, cor.ebew[xx], sd.e)
  results.f[xx] <- analytical.solution(f0_4, kB, s0_4, kW, cor.ebij, cor.ewij, cor.ebew[xx], sd.e)
  
}

my_min <- min(c(results.c, results.d, results.e, results.f))
my_max <- max(c(results.c, results.d, results.e, results.f))

plot(cor.ebew, results.c, ylim=c(my_min-0.05, my_max+0.1), col="black", lwd=2, type="l",
     xlab="", ylab="")
lines(cor.ebew, results.d,  col="blue", lwd=2)
lines(cor.ebew, results.e,  col="red", lwd=2)
lines(cor.ebew, results.f,  col="darkgreen", lwd=2)
text(-.5,.7,"b)", xpd=NA)
#text(0.35,.08,"b)", xpd=NA)

mtext(expression(paste("Cross season synchrony, cor(", epsilon[B], ",", epsilon[W], ")")), 1, outer=F,cex=0.75, line=1.7)
mtext(expression(paste("Population synchrony")), 2, outer=F,cex=0.75, line=1.7)

#### f0 and s0 heat map ####

f0 <- seq(0.3, 2.45, length.out=21)
prop <- seq(-0.9, 0, length.out=21)

kB = 100
kW = 85
cor.ebij = .8
cor.ewij = .2
cor.ebew = 0
sd.e = 0.05
dfrac = 0

results.b <- matrix(NA, nrow=length(f0), ncol=length(prop))

for(xx in 1:length(f0)) {
  for(yy in 1:length(prop)) {
    s0 <- f0[xx]*prop[yy]
    results.b[xx,yy] <- analytical.solution(f0[xx], kB, s0, kW, cor.ebij, cor.ewij, cor.ebew, sd.e)
  }
}

image(f0, prop, results.b, zlim=c(-1,1), col=pal(50),
      xlab="", ylab="", cex=1.25)
contour(f0, prop, results.b, add=T)
text(0.35,.08,"c)", xpd=NA)
#text(1.5, .08, paste(expression("cor.ebij="), cor.ebij, 
#                     expression("& cor.ewij="), cor.ewij), xpd=NA)

mtext(expression(paste("Growth rate (", italic('f')[0], ")")), 1, outer=F,cex=0.75, line=1.7)
mtext(expression(paste("Relative survival (", italic(hat('s'))[0], ")")), 2, outer=F,cex=0.75, line=1.7)



dev.off()



### changing f0 and s0 --  Figure 2 --------------------------------------------------------------------------

#### define parameters across runs ####
tmax = 2000
burn = 1000

#### 1. f0 and s0 heat map ####

f0 <- seq(0.3, 2.45, length.out=21)
#f0 <- c(0.5150, 0.6225, 0.7300, 0.8375)
prop <- seq(-0.9, 0, length.out=21)

kB = 100
kW = 85
cor.ebij = .8
cor.ewij = .2
cor.ebew = 0
sd.e = 0.05
dfrac = 0

results.a <- matrix(NA, nrow=length(f0), ncol=length(prop))
regimes <- matrix(NA, nrow=length(f0), ncol=length(prop))

for(xx in 1:length(f0)) {
  for(yy in 1:length(prop)) {
    s0 <- f0[xx]*prop[yy]
    results.a[xx,yy] <- analytical.solution(f0[xx], kB, s0, kW, cor.ebij, cor.ewij, cor.ebew, sd.e)
    temp <- simmod_det(tmax, f0[xx], kB, s0, kW)
    regimes[xx,yy] <- regime.check(temp)
    #if(regimes[xx,yy] == "overcompensatory") {
    #  paste(expression("xx="), xx, expression("and yy="), yy)}
  }
}
regimes_clean <- ifelse(regimes == "overcompensatory", 1, 0)

#quartz(height=6, width=4)
pdf(file="fig2_f0_s0.pdf", width=4, height=6)
layout(matrix(c(1, 2, 2, 3,
                4, 4, 5, 5,
                6, 6, 7, 7,
                8, 8, 9, 9), nrow=4, byrow=TRUE))
#layout.show(n=9)

pal<-colorRampPalette(colors=c("red","white","blue"))
#par(mfrow=c(3,2), mar=c(1,.5,.5,.5), mgp=c(2.7,0.5,0), tcl=-0.3, oma=c(3,3,3,1))
par(mar=c(1,1,1,.5), mgp=c(2.7,0.5,0), tcl=-0.3, oma=c(3,3,.5,1))
# sim 1
# analytical
plot.new()
image(f0, prop, regimes_clean, zlim=c(-1,1), col=c("white", "darkolivegreen3"),
      xlab="", ylab="", xaxt="n", cex=1.25)
axis(side=1, tick=TRUE, labels=FALSE)
text(x=1, y=-.6, "Undercompensatory", cex=.75)
text(x=1.96, y=-.1, "Overcompensatory", cex=.75)
text(0.35,.08,"a)", xpd=NA)
#mtext(expression(paste("Survival Rate (", italic('s')[0], ")")), 2, outer=F,cex=0.8, line=1.4)
plot.new()

image(f0, prop, results.a, zlim=c(-1,1), col=pal(50),
      xlab="", ylab="", xaxt="n", cex=1.25)
contour(f0, prop, results.a, add=T)
text(0.35,.08,"b)", xpd=NA)
text(1.5, .08, paste(expression("cor.ebij="), cor.ebij, 
                     expression("& cor.ewij="), cor.ewij), xpd=NA)

#### 2. f0 and s0 heat map ####

f0 <- seq(0.3, 2.45, length.out=21)
prop <- seq(-0.9, 0, length.out=21)

kB = 100
kW = 85
cor.ebij = .2
cor.ewij = .8
cor.ebew = 0
sd.e = 0.05
dfrac = 0

results.b <- matrix(NA, nrow=length(f0), ncol=length(prop))
regimes.b <- matrix(NA, nrow=length(f0), ncol=length(prop))


for(xx in 1:length(f0)) {
  for(yy in 1:length(prop)) {
    s0 <- f0[xx]*prop[yy]
    results.b[xx,yy] <- analytical.solution(f0[xx], kB, s0, kW, cor.ebij, cor.ewij, cor.ebew, sd.e)
    temp <- simmod_det(tmax, f0[xx], kB, s0, kW)
    regimes.b[xx,yy] <- regime.check(temp)
  }
}

regimes_clean.b <- ifelse(regimes.b == "overcompensatory", 1, 0)

image(f0, prop, results.b, zlim=c(-1,1), col=pal(50),
      xlab="", ylab="", xaxt="n", yaxt="n", cex=1.25)
contour(f0, prop, results.b, add=T)
text(0.35,.08,"c)", xpd=NA)
text(1.5, .08, paste(expression("cor.ebij="), cor.ebij, 
                     expression("& cor.ewij="), cor.ewij), xpd=NA)

#### 3. f0 and s0 heat map ####

f0 <- seq(0.3, 2.45, length.out=21)
prop <- seq(-0.9, 0, length.out=21)

kB = 100
kW = 85
cor.ebij = .2
cor.ewij = .1
cor.ebew = 0
sd.e = 0.05
dfrac = 0

results.c <- matrix(NA, nrow=length(f0), ncol=length(prop))
regimes.c <- matrix(NA, nrow=length(f0), ncol=length(prop))


for(xx in 1:length(f0)) {
  for(yy in 1:length(prop)) {
    s0 <- f0[xx]*prop[yy]
    results.c[xx,yy] <- analytical.solution(f0[xx], kB, s0, kW, cor.ebij, cor.ewij, cor.ebew, sd.e)
    temp <- simmod_det(tmax, f0[xx], kB, s0, kW)
    regimes.c[xx,yy] <- regime.check(temp)
  }
}

regimes_clean.c <- ifelse(regimes.c == "overcompensatory", 1, 0)


image(f0, prop, results.c, zlim=c(-1,1), col=pal(50),
      xlab="", ylab="", xaxt="n", cex=1.25)
contour(f0, prop, results.c, add=T)
text(0.35,.08,"d)", xpd=NA)
text(1.5, .08, paste(expression("cor.ebij="), cor.ebij, 
                     expression("& cor.ewij="), cor.ewij), xpd=NA)

#image(f0, prop, regimes_clean.c, zlim=c(-1,1), col=pal(50),
#      xlab="", ylab="", cex=1.25)


#### 4. f0 and s0 heat map ####

f0 <- seq(0.3, 2.45, length.out=21)
prop <- seq(-0.9, 0, length.out=21)

kB = 100
kW = 85
cor.ebij = .1
cor.ewij = .2
cor.ebew = 0
sd.e = 0.05
dfrac = 0

results.d <- matrix(NA, nrow=length(f0), ncol=length(prop))

for(xx in 1:length(f0)) {
  for(yy in 1:length(prop)) {
    s0 <- f0[xx]*prop[yy]
    results.d[xx,yy] <- analytical.solution(f0[xx], kB, s0, kW, cor.ebij, cor.ewij, cor.ebew, sd.e)
  }
}

image(f0, prop, results.d, zlim=c(-1,1), col=pal(50),
      xlab="", ylab="", xaxt="n", yaxt="n", cex=1.25)
contour(f0, prop, results.d, add=T)
text(0.35,.08,"e)", xpd=NA)
text(1.5, .08, paste(expression("cor.ebij="), cor.ebij, 
                     expression("& cor.ewij="), cor.ewij), xpd=NA)


#### 5. f0 and s0 heat map ####

f0 <- seq(0.3, 2.45, length.out=21)
prop <- seq(-0.9, 0, length.out=21)

kB = 100
kW = 85
cor.ebij = .6
cor.ewij = .1
cor.ebew = 0
sd.e = 0.05
dfrac = 0

results.e <- matrix(NA, nrow=length(f0), ncol=length(prop))

for(xx in 1:length(f0)) {
  for(yy in 1:length(prop)) {
    s0 <- f0[xx]*prop[yy]
    results.e[xx,yy] <- analytical.solution(f0[xx], kB, s0, kW, cor.ebij, cor.ewij, cor.ebew, sd.e)
  }
}

image(f0, prop, results.e, zlim=c(-1,1), col=pal(50),
      xlab="", ylab="", cex=1.25)
contour(f0, prop, results.e, add=T)
text(0.35,.08,"f)", xpd=NA)
text(1.5, .08, paste(expression("cor.ebij="), cor.ebij, 
                     expression("& cor.ewij="), cor.ewij), xpd=NA)


#### 6. f0 and s0 heat map ####

f0 <- seq(0.3, 2.45, length.out=21)
prop <- seq(-0.9, 0, length.out=21)

kB = 100
kW = 85
cor.ebij = .1
cor.ewij = .6
cor.ebew = 0
sd.e = 0.05
dfrac = 0

results.f <- matrix(NA, nrow=length(f0), ncol=length(prop))

for(xx in 1:length(f0)) {
  for(yy in 1:length(prop)) {
    s0 <- f0[xx]*prop[yy]
    results.f[xx,yy] <- analytical.solution(f0[xx], kB, s0, kW, cor.ebij, cor.ewij, cor.ebew, sd.e)
  }
}

image(f0, prop, results.f, zlim=c(-1,1), col=pal(50),
      xlab="", ylab="", yaxt="n", cex=1.25)
contour(f0, prop, results.f, add=T)
text(0.35,.08,"g)", xpd=NA)
text(1.5, .08, paste(expression("cor.ebij="), cor.ebij, 
                     expression("& cor.ewij="), cor.ewij), xpd=NA)

mtext(expression(paste("Growth Rate (", italic('f')[0], ")")), 1, outer=T,cex=0.8, line=1.2)
mtext(expression(paste("Relative Survival (", italic(hat('s'))[0], ")")), 2, outer=T,cex=0.8, line=1.2)
#mtext(expression(paste("Survival Rate (", italic('s')[0], ")")), 2, outer=T,cex=0.8, line=1.2)
dev.off()
