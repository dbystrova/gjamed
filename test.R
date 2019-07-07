############################################
## Example 3
############################################

## SAMPLING MODEL
## p(yi | F) = G(yi), yi>0

## PRIOR
## F = DP(M G0)

#############################################
## without censoring -- DP for y[i], y>0 only
#############################################

## data
## healthy Tconv mouse 2
yf <- c(37,11,5,2) # frequencies
xf <- c(1,2,3,4)   # counts
y <- rep(xf,yf)    # raw data
n <- length(y)
k <- n             # initialize

## hyperparameters
## a0 <- 1;   b0 <- 1   # hyperprior b ~ Ga(a0,b0) 
a <- 1; b <- .05     # G0(mu) = Ga(a,b)
lambda <- 300        # k ~ Poi(lambda)

M <- 1
H <- 10
N <- 25; p=8

rdiric<- function(n,a) {
  ## generates x ~ Dir(a,...a) (n-dim)
  p <- length(a)
  m <- matrix(nrow=n,ncol=p)
  for (i in 1:p) {
    m[,i] <- rgamma(n,a[i])
  }
  sumvec <- m %*% rep(1,p)
  m / as.vector(sumvec)
}


sample.dponly <- function(N=10,M=1,p=8,
                          plt.G0=F,plt.spl=F,plt.Ghat=F,
                          cdf=T)
{
  ## generates posterior p(G | y)
  ## for yi ~ G
  ## G0:   prior mean
  ## Ghat: empirical
  ## Gbar: posterior mean
  ## G:    posterior sample
  xgrid <- 1:p
  r <-  1/(1-dpois(0,lambda=2)) # trunction to x>=1
  G0 <- dpois(xgrid,lambda=2)*r # prior base measure
  G1 <- M*G0                    # post base measure  
  G1[xf] <- G1[xf]+yf     # +1 because xgrid starts at 0
  G <- rdiric(N,G1)
  Gcdf <- apply(G,1,cumsum)
  n <- sum(yf)
  Gbar <- G1/(n+M)
  Gbarcdf <- cumsum(Gbar)
  if (cdf)
    matplot(xgrid,Gcdf,type="n",bty="l",
            xlim=c(1,10),xlab="X",ylab="G",ylim=c(0,1))
  else
    matplot(xgrid,t(G),xlim=c(1,8),type="n",bty="l",
            xlab="COUNT",ylab="G") 
  if (plt.spl){
    for(i in 1:N){
      if (cdf)
        cdfplt(xgrid,Gcdf[,i],lw=1,lt=i,hi=10)
      else
        lines(xgrid,G[i,],lw=1,lt=i)
    }
  }
  G0cdf <- cumsum(G0)
  if (plt.G0){
    if (cdf)
      cdfplt(xgrid,G0cdf,  lt=3,lw=3, cl=1)
    else
      lines(xgrid,G0,  lt=3,lw=3, col=1)
  }
  Ghat <- rep(0,p) # initialize
  n <- sum(yf)
  Ghat[xf] <- yf/n
  Ghatcdf <- cumsum(Ghat)
  if (plt.Ghat){
    if (cdf){
      cdfplt(xgrid,Ghatcdf, lt=1, lw=3, cl=1)
      cdfplt(xgrid,Gbarcdf, lt=2, lw=3, cl=1)
    } else{
      xg <- as.numeric(names(table(y)))
      lines(table(y)/n,lwd=1)
      points(xg,table(y)/n,pch=19)
      lines(xgrid,Gbar,lty=1,lwd=3)
    }
  }# plt.Ghat
}

#############################################
## plotting a cdf -- aux function
#############################################

cdfplt <- function(x,y,lo=NULL,hi=NULL,
                   lw=1,lt=1,cl=1)
{
  p <- length(x)
  if (!is.null(hi)) # final line segment
    if(hi > x[p])
      lines(c(x[p],hi), c(1,1),col=cl,lwd=lw,lty=lt)
  if (!is.null(lo)){ # initial line segment
    if (lo<x[1]){
      lines(c(lo,x[1]), c(0,0),col=cl,lwd=lw,lty=lt)
      if (y[1]>0)   # prob mass at x[1]
        points(x[1],0,pch=1)
    }
  }# initial seg
  ylag <- c(0,y)    # prev prob
  for(i in 1:p){
    if (i<p)
      lines(x[c(i,i+1)],y[c(i,i)],col=cl,lwd=lw,lty=lt)
    if (y[i]>ylag[i]){ # prob mass @ x[i]
      points(x[i],y[i],pch=19)
      points(x[i],ylag[i],pch=1)
    }
  }# for i
}



## RUN IT: 
## run the commands below - best line by line

ex <- function()
{
  sample.dponly(plt.spl=T,cdf=F,plt.Ghat=T)
  legend(3,0.65, legend=c("E(G | y)", "G ~ p(G | y)", "Ghat (pins)"),
         lty=c(1,2,1), lwd=c(3,1,1),bty="n")
}


############################################
## Example 3
## Tcell
## DP mixture model
############################################

## SAMPLING MODEL
## p(yi | F) = F(yi)/(1-F(0))                   yi
## p(n | k)  = Bin(n; k, 1-F(0))                n
## p(y,n | k,F) = (k n) F(0)^(k-n) prod F(yi)   Splg model

## PRIOR
## F = \int Poi(mu) dG(mu),
## F0 = F(0); F+ = 1-F0
## G ~ DP(G0,M), G0(mu) = Ga(a,b)

## complete conditionals
## p(k | F,n) \propto (k n) F0^(k-n) p(k)
## p(ri | G), i=1..n
##   r0h = sum(ri==h, i=n+1..k), h=1..H
## p(mh | r)
## p(wh | r)

## data
## healthy Tconv mouse 2
yf <- c(37,11,5,2) # frequencies
xf <- c(1,2,3,4)   # counts
y <- rep(xf,yf)    # raw data
n <- length(y)
k <- n             # initialize

## hyperparameters
a <- 1; b <- .05     # G0(mu) = Ga(a,b)
lambda <- 300        # k ~ Poi(lambda)

M <- 1
H <- 10

# ##################################################################
# initialize clustering..
# ##################################################################

init.DPk <- function(H=10)
{ ## inital EDA estimate of G = sum_{h=1..10} w_h delta(m_h)
  ## returns:
  ##   list(mh,wh)
  ## use (mh,wh) to initialize the blocked Gibbs
  ## r[i]: latent cluster membership
  ## For i=n+1..k, summarize cluster memberships by h=1..H:
  ##     Sh0[h]: number of i=h, for i=n+1...k
  
  ## cluster data, and cut at height H=10, to get 10 clusters
  H1 <- min(round(length(y)/2),H)
  hc <- hclust(dist(y)^2, "cen")
  r  <- cutree(hc, k = H1)
  ## record cluster specific means, order them 
  mh1 <- sapply(split(y,r),mean)    # cluster specific means == m_h
  wh1 <- table(r)/n
  idx <- order(wh1,decreasing=T)    # re-arrange in deceasing order
  mh <- mh1[idx]
  wh <- wh1[idx]
  if (H1<H){
    H0 <- H-H1
    mh0 <- rpois(H0,lambda=mean(mh))+1 # arbitrary way to generate more..
    wh0 <- rdiric(1,rep(1,H0))
    mh <- c(mh,mh0)
    W0 <- 0.1
    wh <- c((1-W0)*wh, W0*wh0)
  }
  return(list(mh=mh,wh=wh,r=r,Sh0=rep(0,H)))
}   




# ##################################################################
#  Blocked GS
# ##################################################################


gibbs <- function(n.iter=100,diter=5,n.iter0=10)
{ ## returns fgrid: (n.iter x 9) matrix of imputed F
  ## klist:         posterior on k = latent sample size
  ## blist:         hyper-parameter of DP of Poi mixture
  
  G <- init.DPk(H)
  ## b <- a0/b0
  k <- n
  ## prepare for sample.k(.) below
  kgrid <- n:(lambda)
  lp0 <- lgamma(kgrid+1)-lgamma(kgrid-n+1)
  ## + dpois(kgrid,lambda=lambda,log=T) [using k ~ U[0,lambda]
  ## missing: ... +(ktgrid-n)*lF0+
  ## will add in MCMC (F0 changes each iter..)
  
  ## data structures to save imputed F ~ p(F | ...)
  xgrid <- 0:8
  fgrid <- NULL
  klist <- NULL
  F0list <- NULL
  pkgrid <- 0*kgrid
  blist <- NULL
  plot(table(y)/n,xlim=c(0,8),bty="l",xlab="COUNT",ylab="FREQ")
  # plot(density(y),xlab="X",ylab="Y",bty="l",type="l",
  #   xlim=c(0,6),ylim=c(0,0.7), main="")
  ## Gibbs
  for(iter in 1:n.iter){
    kp <-  sample.k(G,k,lp0,kgrid)
    k <- kp$k
    G$r  <- sample.r(G)   # 1. r_i ~ p(r_i | ...), i=1..n
    G$Sh0 <- sample.r0(G,k)
    G$mh <- sample.mh(G,b)     # 2. m_h ~ p(m_h | ...), h=1..H
    G$vh <- sample.vh(G)       # 3. v_h ~ p(v_h | ...), h=1..H
    
    ## record draw F ~ p(F | th,sig,y) (approx)
    if ( (iter %% diter == 0) & (iter>n.iter0) ){
      f   <- fbar.H(xgrid,G)
      lines(xgrid,f,col=iter,lty=3)
      fgrid <- rbind(fgrid,f)
      klist <- c(klist,k)
      F0list <- c(F0list, kp$F0)
      pkgrid <- pkgrid+kp$p
      blist <- c(blist,b)
    }
  }
  ## add overall average (= posterior mean) to the plot
  fbar <- apply(fgrid,2,mean)
  lines(xgrid,fbar,lwd=3,col=2)
  pk <- pkgrid/n.iter
  return(list(fgrid=fgrid,klist=klist,blist=blist,pk=pk,
              kgrid=kgrid, F0list=F0list))
}

plt.f <- function(fgrid,pltfgrid=T,pltsim=T,initial=T,cl=NULL,lt=3,
                  norm=F)
{
  ## norm=T: normalize G to (1-G(0))=1 (to make comparable with data..)
  xgrid <- 0:8
  plot(0:8,seq(0,0.7,length=9),type="n",bty="l",xlab="COUNT",ylab="FREQ")
  lines(table(y)/n)  ##,xlim=c(0,8),bty="l",xlab="COUNT",ylab="FREQ")
  ## plot(table(y)/n,xlim=c(0,8),bty="l",xlab="COUNT",ylab="FREQ")
  if (pltsim){
    if (pltfgrid){
      nsim <- nrow(fgrid)
      if (initial)
        nsim0 <- floor(nsim/2)
      else
        nsim0 <- 1
      for(i in nsim0:nsim){
        fi <- fgrid[i,]
        if (norm)
          fi <- fi/(sum(fi[-1])) # normalize to (1-G(0)) = 1
        if (is.null(cl))
          lines(xgrid,fi,col=i,lty=lt)
        else
          lines(xgrid,fi,col=cl,lty=lt)
      }
      ## add overall average (= posterior mean) to the plot
      fbar <- apply(fgrid[nsim0:nsim,],2,mean)
      lines(xgrid,fbar,lwd=3,col=1)
    }
  } #plotsim
  lines(table(y)/n,lwd=1)  ##,xlim=c(0,8),bty="l",xlab="COUNT",ylab="FREQ")
  xg <- as.numeric(names(table(y)))
  points(xg,table(y)/n,pch=19)  ##,xlim=c(0,8),bty="l",xlab="COUNT",ylab="FREQ")
}

sample.k <- function(G,k,lp0,kgrid)
{ ## lp0: part of log(p(k | y) that does not involve F0
  F0 <- sum(exp(-G$mh)*G$wh)
  lF0 <- log(F0)
  lF <- log(1-F0)
  n <-  length(y)
  lp <- lp0+(kgrid-n)*lF0
  p <- exp(lp-max(lp))
  p <- p / sum(p) ## normalize
  k <- sample(kgrid,1,prob=p)
  return(list(k=k,F0=F0,p=p))
}

sample.r <- function(G)
{ ## samle allocation indicators
  
  r <- rep(0,n)
  for(i in 1:n){
    ph <-   dpois(y[i],lambda=G$mh)*G$wh # likelihood   * prior
    ## p(yi | ri=h) * w_h
    r[i] <- sample(1:H,1,prob=ph)
  }
  return(r)
}

sample.r0 <- function(G,k)
{ ## samle allocation indicators
  
  p0 <- exp(-G$mh)*G$wh
  r0 <- sample(1:H, k-n, replace=T, prob=p0)
  Sh0 <- rep(0,H)
  for(h in 1:H)
    Sh0[h] <- sum(r0==h)
  return(Sh0)
}


sample.mh <- function(G,b)
{ ## sample mh ~ p(mh | ...)
  ##
  
  mh <- rep(0,H)     # initialize
  for(h in 1:H){
    Sh <- which(G$r==h) # Sh = {i: r[i]=h
    nh <- length(Sh)+G$Sh0[h]
    sumy <- sum(y[Sh])
    mh[h] <- rgamma(1,shape=a+sumy, rate=b+nh)
  }
  return(mh)
}


sample.vh <- function(G)
{## sample vh ~ p(vh | ...)
  ## returns: wh
  
  vh <- rep(0,H)  # initialize
  wh <- rep(0,H)
  V <-  1         # record prod_{g<h} (1-vh_h)
  for(h in 1:(H-1)){
    Ah <- which(G$r==h)
    nAh <- length(Ah)+G$Sh0[h]
    Bh <- which(G$r>h)
    nBh <- length(G$Bh) + sum(G$Sh0[(h+1):H])
    vh[h] <-  rbeta(1, 1+nAh, M+nBh)
    wh[h] <- vh[h]*V
    V <- V*(1-vh[h])
  }
  vh[H] <- 1.0
  wh[H] <- V
  return(wh)
}

fbar.H <- function(xgrid,G)
{ ## return a draw F ~ p(F | ...) (approx)
  
  fx <- rep(0,length(xgrid))
  for(h in 1:H)
    fx <- fx + G$wh[h]*dpois(xgrid,lambda=G$mh[h])
  return(fx)
}



plt.k <- function(fk,pltk=F,traj=F,hist=F,trans=.5,
                  pltkF0=F)
{ ##  makes plots for p(k | y) (pltk=T)
  ##  discard first trans% of iterations as initial transient
  
  klist <- fk$klist
  F0list <- fk$F0list
  nsim <- length(klist)
  its <- 1:nsim                 # iterations (for x-axis)
  ## not currently used
  if (trans>0){
    nsim0 <- floor( trans*nsim )
    klist <- klist[nsim0:nsim]
    F0list <- F0list[nsim0:nsim]
    its <- nsim0:nsim
  }
  
  if (pltk){
    if (hist)
      plot(table(klist)/M,xlab="k",ylab="p(k|y)",bty="l")
    else if (traj)
      plot(klist,type="l",bty="l",xlab="ITERATION",ylab="k")
    else
      plot(fk$kgrid, fk$pk, xlab="k", type="l",
           ylab="p(k | y)", bty="n")
  }
  if (pltkF0){
    plot(klist,F0list,pch="x",xlab="k",ylab="F0",bty="l")
  }
}

#############################################
## step function
#############################################

cdfplt <- function(x,y,lo=NULL,hi=NULL,
                   lw=1,lt=1,cl=1)
{
  p <- length(x)
  if (!is.null(hi)) # final line segment
    if(hi > x[p])
      lines(c(x[p],hi), c(1,1),col=cl,lwd=lw,lty=lt)
  if (!is.null(lo)){ # initial line segment
    if (lo<x[1]){
      lines(c(lo,x[1]), c(0,0),col=cl,lwd=lw,lty=lt)
      if (y[1]>0)   # prob mass at x[1]
        points(x[1],0,pch=1)
    }
  }# initial seg
  ylag <- c(0,y)    # prev prob
  for(i in 1:p){
    if (i<p)
      lines(x[c(i,i+1)],y[c(i,i)],col=cl,lwd=lw,lty=lt)
    if (y[i]>ylag[i]){ # prob mass @ x[i]
      points(x[i],y[i],pch=19)
      points(x[i],ylag[i],pch=1)
    }
  }# for i
}








##################################################################
## execute the macros
## RUN it: 
##    run the commands below - best line by line

ex <- function()
{
  ## -------------------------------------------------------
  ## density estimation with observed data only (Sec 2.1.2)
  
  sample.dponly(plt.spl=T,cdf=F,plt.Ghat=T)
  legend(5,0.65, legend=c("E(G | y)", "G ~ p(G | y)", "Ghat (pins)"),
         lty=c(1,2,1), lwd=c(3,1,1),bty="n")
  
  ## same, showing also G0 and plotting it as pdf
  sample.dponly(plt.G0=T,plt.Ghat=T)
  legend(7,0.2,lty=1:3,legend=c("data","E(G|y)","G0"),
         bty="n",lwd=3, seg.len=3)
  
  ## same, showing G0 and Ghat as p.m.f.
  sample.dponly(plt.G0=T,plt.Ghat=T,cdf=F,plt.spl=F)
  legend(5.5,0.7,lty=c(1:2,1),
         legend=c("E(G|y)","G0","data"),bty="n",
         seg.len=3.5,lwd=c(3,3,1))
  
  ## -------------------------------------------------------
  ## density estimation with censoring (Sec 2.2.1)
  
  
  lambda <<-  234
  fk234 <- gibbs(5000,diter=20,n.iter0=1000)
  plt.f(fk$fgrid,F,F)
  plt.f(fk234$fgrid,cl="darkgrey",lt=2)
  legend(3,0.6,lty=c(1,2,1),lwd=c(3,1,1),bty="n",
         legend=c("E(G | data)","G ~ p(G | data)","data (pins)"))
  
  ## trajectory of k
  plt.k(fk234,pltk=T,traj=T)
  ## marg posterior in k
  plt.k(fk234,pltk=T)
}
