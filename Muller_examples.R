#dir <- "/home/pmueller/pap/12/ba/eig121/"
#setwd(dir)
setwd("~/Documents/GitHub/gjamed")
############################################
## Example 4 - Gene expression
## DPM
## implementing MCMC with finite DP
## Sec 2.4.6.

## MODEL:
## y_i | th_i ~ N(th_i, sig^2)
## th_i       ~ G and
## G          ~ DP(M, G0) with G0=N(0,4)
## hyperprior:
## 1/sig ~ Ga(a,b), a=1, b=1
##
## We use the usual notation for unique values ths[1..k] and
##        s_i = j if th_i = ths_j
## Denote with F(y) = \int N(y; th,sig) dG(th)
##        and  f(y) = pdf
##

## Posterior MCMC for DPM models
## (b) using the blocked Gibbs sampler
##         H
##     G = sum w_h delta(m_h)
##         h=1
##     recall: w_h = v_h (1-\sum_{g<h} w_g)
##             r_i = h iff th_i = m_h
##
##     We iterate over the following transition probabilities, sampling
##     from complete conditional posterior distributions:
##     1. r_i   ~ p(r_i   | ... ), i=1..n
##     2. m_h   ~ p(m_h   | ... ), h=1..H
##     3. v_h   ~ p(v_h   | ... ), h=1..H
##     4. sig   ~ p(sig   | ... )


## EIG 121 data
read.dta <- function()
{
  X <- read.table("EIG.txt",header=1)
  y <- X$recorded
  y <- log( y[!is.na(y)] )  # log EIG121 expression
  n <- length(y)
  return(dta=list(y=y, n=n))
}


## hyperparameters
a <- 1;   b <- 1     # 1/sig ~ Ga(a,b)
m0 <- -3;  B0 <- 4    # G0 = N(m0,B0)
M <- 1
H <- 10

# ##################################################################
# initialize clustering..
# ##################################################################

init.DPk <- function()
{ ## inital EDA estimate of G = sum_{h=1..10} w_h delta(m_h)
  ## returns:
  ##   list(mh,wh)
  ## use (mh,wh) to initialize the blocked Gibbs
  
  ## cluster data, and cut at height H=10, to get 10 clusters
  hc <- hclust(dist(y)^2, "cen")
  r  <- cutree(hc, k = 10)
  ## record cluster specific means, order them 
  mh1 <- sapply(split(y,r),mean)    # cluster specific means == m_h
  wh1 <- table(r)/n
  idx <- order(wh1,decreasing=T)    # re-arrange in deceasing order
  mh <- mh1[idx]
  wh <- wh1[idx]
  return(list(mh=mh,wh=wh,r=r))
}   




# ##################################################################
# 2. Blocked GS
# ##################################################################


gibbs.H <- function(n.iter=500)
{
  
  G <- init.DPk()
  sig <- 0.11
  
  ## data structures to save imputed F ~ p(F | ...)
  xgrid <- seq(from= -10, to=2,length=50)
  fgrid <- NULL
  plot(density(y),xlab="X",ylab="Y",bty="l",type="l",
       xlim=c(-10, 2),ylim=c(0,0.4), main="")
  ## Gibbs
  for(iter in 1:n.iter){
    G$r <-  sample.r(G$wh,G$mh,sig)   # 1. r_i ~ p(r_i | ...), i=1..n
    G$mh <- sample.mh(G$wh,G$r,sig)   # 2. m_h ~ p(m_h | ...), h=1..H
    G$vh <- sample.vh(G$r)            # 3. v_h ~ p(v_h | ...), h=1..H
    th <- G$mh[G$r]                   # record implied th[i] = mh[r[i]]
    sig <- sample.sig(th)       # 4. sig ~ p(sig | ...)
    
    ## record draw F ~ p(F | th,sig,y) (approx)
    f   <- fbar.H(xgrid,G$wh,G$mh,sig)
    lines(xgrid,f,col=iter,lty=3)
    fgrid <- rbind(fgrid,f)
  }
  ## add overall average (= posterior mean) to the plot
  fbar <- apply(fgrid,2,mean)
  lines(xgrid,fbar,lwd=3,col=2)
  return(fgrid)
}

sample.r <- function(wh,mh,sig)
{ ## samle allocation indicators
  
  r <- rep(0,n)
  for(i in 1:n){
    ph <-   dnorm(y[i],m=mh,sd=sig)*wh # likelihood   * prior
    ## p(yi | ri=h) * w_h
    r[i] <- sample(1:H,1,prob=ph)
  }
  return(r)
}


sample.mh <- function(wh,r,sig)
{ ## sample mh ~ p(mh | ...)
  ##
  
  mh <- rep(0,H)     # initialize
  for(h in 1:H){
    if(any(r==h)){      # some data assigned to h-th pointmass
      Sh <- which(r==h) # Sh = {i: r[i]=h
      nh <- length(Sh)
      ybarh <- mean(y[Sh])
      varh   <- 1.0/(1/B0 + nh/sig^2)
      meanh  <- varh*(1/B0*m0 + nh/sig^2*ybarh)
    } else {            # no data assinged to h-th pointmass
      varh  <- B0       # sample from base measure
      meanh <- m0
    }
    mh[h] <- rnorm(1,m=meanh,sd=sqrt(varh))
  }
  return(mh)
}

sample.vh <- function(r)
{## sample vh ~ p(vh | ...)
  ## returns: wh
  
  vh <- rep(0,H)  # initialize
  wh <- rep(0,H)
  V <-  1         # record prod_{g<h} (1-vh_h)
  for(h in 1:(H-1)){
    Ah <- which(r==h)
    Bh <- which(r>h)
    vh[h] <-  rbeta(1, 1+length(Ah), M+length(Bh))
    wh[h] <- vh[h]*V
    V <- V*(1-vh[h])
  }
  vh[H] <- 1.0
  wh[H] <- V
  return(wh)
}

fbar.H <- function(xgrid,wh,mh,sig)
{ ## return a draw F ~ p(F | ...) (approx)
  
  fx <- rep(0,length(xgrid))
  for(h in 1:H)
    fx <- fx + wh[h]*dnorm(xgrid,m=mh[h],sd=sig)
  return(fx)
}


sample.sig <- function(th)
{ ## sample
  ##   sig ~ p(sig | ...)
  ## returns: sig
  
  s2 <- sum( (y-th)^2 )    # sum of squared residuals
  a1 <- a+0.5*n
  b1 <- b+0.5*s2
  s2.inv <- rgamma(1,shape=a1,rate=b1)
  return(1/sqrt(s2.inv))
}

plt.all <- function(fgrid,sim=T,dens=T)
{
  xgrid <- seq(from= -10, to=2,length=50)
  M <- nrow(fgrid)
  idx0 <- 21:M
  fgrid0 <- fgrid[idx0,]            # drop initial transient
  idx1 <- which(idx0 %% 5 == 0 )    # thin out for plotting
  fgrid1 <- fgrid[idx1,]
  fbar <- apply(fgrid0,2,mean)
  plot(xgrid,fbar,xlab="log(EIG121)", ylab="G",
       type="l",lwd=3,bty="l",ylim=c(0,0.9))
  if(sim){
    matlines(xgrid,t(fgrid1),col=1)
    ## matlines(xgrid,t(fgrid0),col=2)
    lines(xgrid,fbar,type="l",lwd=4,col="grey")
  }
  if (dens){
    lines(density(y),col="yellow",lty=2,lwd=4)
  }
}

plt.dta <- function()
{
  hist(y,main="",xlab="log(EIG121)",ylab="FREQ",prob=T)
}
##################################################################
## RUN:
## execute the commands below -- best line by line

ex <- function()
{
  dta <- read.dta()
  attach(dta)
  
  ## run MCMC
  fgrid <- gibbs.H()
  plt.dta()
  plt.all(fgrid)
}





z_scores <- seq(-3, 3, by = .1)
pvalues <- pnorm(z_scores)

# Now we'll plot these values
plot(pvalues, # Plot where y = values and x = index of the value in the vector
     xaxt = "n", # Don't label the x-axis
     type = "l", # Make it a line plot
     main = "cdf of the Standard Normal",
     xlab= "Quantiles",
     ylab="Probability Density") 

# These commands label the x-axis
axis(1, at=which(pvalues == pnorm(-2)), labels=round(pnorm(-2), 2))
axis(1, at=which(pvalues == pnorm(-1)), labels=round(pnorm(-1), 2))
axis(1, at=which(pvalues == pnorm(0)), labels=c(.5))
axis(1, at=which(pvalues == pnorm(1)), labels=round(pnorm(1), 2))
axis(1, at=which(pvalues == pnorm(2)), labels=round(pnorm(2), 2))

