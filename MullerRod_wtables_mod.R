rm(list=ls())
###Data generation
gen_data<- function(mu_vec=c(1,5,3,5), sigma_vec=rep(0.1, 4), ns=50, prob_vec=c(0.25,0.25,0.25,0.25) ){
  l<- length(mu_vec)
  components <- sample(1:l,prob=prob_vec,size=ns,replace=TRUE)
  mus <- mu_vec
  sds <- sqrt(sigma_vec)
  samples <- rnorm(n=ns,mean=mus[components],sd=sds[components])
  return(samples)
}

set.seed(123)
sample<- gen_data(mu_vec=c(1,10),sigma_vec=rep(0.2, 2),ns=50,prob_vec=c(0.5,0.5))
hist(sample,probability=TRUE,breaks=20,col=grey(.9))
plot(density(sample,bw = "sj"))



## Dirichlet process BM ############################################################


#pdf("Output.pdf")
#y <- sample
#n <- length(y)
## hyperparameters
a  <- 1;   b <- 1   # 1/sig ~ Ga(a,b)
m0 <- 0;  B0 <- 1  # G0 = N(m0,B0)
#M  <- 1

##########INIT##########
init.DPk <- function(ydata)
{ ## inital EDA estimate of th[1..n]
  ## cluster data, and cut at height H=10, to get 10 clusters
  y<- ydata
  hc <- hclust(dist(y)^2, "cen")
  s  <- cutree(hc, k = 10)
  ths <- sapply(split(y,s),mean)
  th <-  ths[s]
  return(th)
}
# cluster membership indicators
# cluster specific means
# return th_i = ths[ s_i ]

##########INIT##########
sample.th <- function(th,sig)
{ ## sample
  ##    th[i]  ~ p(th_i | th[-i],sig,y)
  ## returns updated th vector
  for(i in 1:n){
    ## unique values and counts
    nj <- table(th[-i])
    ths <- as.numeric(names(nj))
    k <- length(nj)
    ## likelihood
    fj <- dnorm(y[i], m=ths, s=sig)
    f0 <- dnorm(y[i], m=m0,   s=sqrt(B0+sig^2)) # q0
    pj <- c(fj*nj,f0*M)                       # p(s[i]=j | ...), j=1..k, k+1
    # pj_norm<- pj/sum(pj)
    s <- sample(1:(k+1), 1, prob=pj)          # sample s[i]
    if (s==k+1){ ## generate new th[i] value
      v.new <- 1.0/(1/B0 + 1/sig^2)
      m.new <- v.new*(1/B0*m0 + 1/sig^2*y[i])
      thi.new   <- rnorm(1,m=m.new,sd=sqrt(v.new))
      ths <- c(ths,thi.new)
    }   
    th[i] <- ths[s]
  }
  return(th) 
} 





sample.ths <- function(th,sig)
{ ## sample ths[j] ~ p(ths[j] | ....)
  ##               = N(ths[j]; m0,B0) * N(ybar[j]; ths[j], sig2/n[j])
  ##               = N(ths[j]; mj,vj)
  ## unique values and counts
  nj <- table(th)                 # counts
  ths <- sort(unique(th))         # unique values
  ##use sort(.) to match table counts
  k <- length(nj)
  for(j in 1:k){
    ## find Sj={i: s[i]=j} and compute sample average over Sj
    idx <- which(th==ths[j])
    ybarj <- mean(y[idx])
    ## posterior moments for p(ths[j] | ...)
    vj <- 1.0/(1/B0 + nj[j]/sig^2)
    mj <- vj*(1/B0*m0 + nj[j]/sig^2*ybarj)
    thsj <- rnorm(1,m=mj,sd=sqrt(vj))
    ## record the new ths[j] by replacing all th[i], i in Sj.
    th[idx] <- thsj
  }
  return(th) 
}  



sample.sig <- function(th)
{ ## sample
  ##   sig ~ p(sig | ...)
  ## returns: sig
  s2 <- sum( (y-th)^2 )
  a1 <- a+0.5*n
  b1 <- b+0.5*s2
  s2.inv <- rgamma(1,shape=a1,rate=b1)
  return(1/sqrt(s2.inv))
}


fbar <- function(x,th,sig)
{ ## conditional draw F ~ p(F | th,sig,y) (approx -- will talk about this..)
  nj <- table(th)              # counts
  ths <- as.numeric(names(nj)) # unique values
  k <- length(nj)
  fx <- M/(n+M)*dnorm(xgrid,m=m0,sd=sqrt(B0+sig^2))
  for(j in 1:k)
    fx <- fx + nj[j]/(n+M)*dnorm(xgrid,m=ths[j],sd=sig)
  return(fx)
}



ebar <- function(x,th,sig)
{ ## conditional draw F ~ p(F | th,sig,y) (approx -- will talk about this..)
  nj <- table(th)              # counts
  ths <- as.numeric(names(nj)) # unique values
  k <- length(nj)
  fx <- M/(n+M)*pnorm(xgrid,m=m0,sd=sqrt(B0+sig^2))
  for(j in 1:k)
    fx <- fx + nj[j]/(n+M)*pnorm(xgrid,m=ths[j],sd=sig)
  return(fx)
}


gibbs <- function(n.iter=1000)
{
  th <- init.DPk()              ## initialize th[1..n]
  sig <- sqrt( mean((y-th)^2))  ## and sig
  ## set up data structures to record imputed posterior draws..1
  xgrid <- seq(from=-2,to=15,length=100)
  fgrid <- NULL     ## we will record imputed draws of f
  ecfgrid <- NULL     ## we will record imputed draws of f
  njlist <- NULL    ## record sizes of 8 largest clusters
  klist <- NULL
  sig_trace<- NULL
  ## start with a plot of a kernel density estimate of the data
  d<-density(y,bw = "sj")
  plot(d,xlab="X",ylab="Y",bty="l",type="l",xlim=c(-2,15),ylim=c(0,1), main="")
  ## now the Gibbs sampler
  for(iter in 1:n.iter){
    th  <- sample.th(th,sig)
    sig <- sample.sig(th)
    th  <- sample.ths(th,sig)
    ## update running summaries #################
    sig_trace<- rbind(sig_trace,sig)
    f   <- fbar(xgrid,th,sig)
    ecf<- ebar(xgrid,th,sig)
    lines(xgrid,f,col=iter,lty=3)
    fgrid <- rbind(fgrid,f)
    ecfgrid<- rbind(ecfgrid,ecf)
    nj <- table(th)                 # counts
    njlist <- rbind(njlist,sort(nj,decr=T)[1:8])
    klist <- c(klist,length(nj))
  }
  ## report summaries ###########################
  fbar <- apply(fgrid,2,mean)
  lines(xgrid,fbar,lwd=3,col=2)
  ecfbar <- apply(ecfgrid,2,mean)
  plot(ecdf(y),xlab="X",ylab="Y",bty="l",xlim=c(-2,15),ylim=c(0,1), main="")
  lines(xgrid,ecfbar,lwd=3,col=2)
  for(iter in 1:n.iter){ lines(xgrid,ecfgrid[iter,],col=iter,lty=3) }
  #ks_test<- ks.test(ecfbar,ecdf(y))
  njbar <- apply(njlist,2,mean,na.rm=T)
  cat("Average cluster sizes:\n",format(njbar),"\n")
  pk <- table(klist)/length(klist)
  cat("Posterior probs p(k): (row1 = k, row2 = p(k) \n ")
  print(pk/sum(pk))
  return(list(fgrid=fgrid,klist=klist,njlist=njlist, sig=sig_trace ))
}

#mc<- gibbs()
#plot(mc$sig)
#ks.test(f_approx$ecf, ecdf(X))
# 
# k<- ecdf(y)
# mcmc <- gibbs()
# names(mcmc)
# ## report summaries
# njbar <- apply(mcmc$njlist,2,mean,na.rm=T)
# cat("Average cluster sizes:\n",format(njbar,digits=1),"\n")
# pk <- table(mcmc$klist)/length(mcmc$klist)
# cat("Posterior probs p(k): (row1 = k, row2 = p(k) \n ")
# print(pk/sum(pk))
# 





#####################PY  process sampling ##########################################################


#initial parameters
#M<-1
#sigma<-0.25

sample.th_py <- function(th,sig,M,sigma,ydata)
{ ## sample
  ##    th[i]  ~ p(th_i | th[-i],sig,y)
  ## returns updated th vector
  y<- ydata
  n<- length(y)
  for(i in 1:n){
    ## unique values and counts
    nj <- table(th[-i])
    ths <- as.numeric(names(nj))
    k <- length(nj)
    ## likelihood
    fj <- dnorm(y[i], m=ths, s=sig)
    f0 <- dnorm(y[i], m=m0,   s=sqrt(B0+sig^2)) # q0
    pj <- c(fj*(nj- sigma) ,f0*(M+ k*sigma))                       # p(s[i]=j | ...), j=1..k, k+1
    # pj_norm<- pj/sum(pj)
    s <- sample(1:(k+1), 1, prob=pj)          # sample s[i]
    if (s==k+1){ ## generate new th[i] value
      v.new <- 1.0/(1/B0 + 1/sig^2)
      m.new <- v.new*(1/B0*m0 + 1/sig^2*y[i])
      thi.new   <- rnorm(1,m=m.new,sd=sqrt(v.new))
      ths <- c(ths,thi.new)
    }   
    th[i] <- ths[s]
  }
  return(th) 
} 


sample.ths_py <- function(th,sig,ydata)
{ ## sample ths[j] ~ p(ths[j] | ....)
  ##               = N(ths[j]; m0,B0) * N(ybar[j]; ths[j], sig2/n[j])
  ##               = N(ths[j]; mj,vj)
  ## unique values and counts
  y<- ydata
  nj <- table(th)                 # counts
  ths <- sort(unique(th))         # unique values
  ##use sort(.) to match table counts
  k <- length(nj)
  for(j in 1:k){
    ## find Sj={i: s[i]=j} and compute sample average over Sj
    idx <- which(th==ths[j])
    ybarj <- mean(y[idx])
    ## posterior moments for p(ths[j] | ...)
    vj <- 1.0/(1/B0 + nj[j]/sig^2)
    mj <- vj*(1/B0*m0 + nj[j]/sig^2*ybarj)
    thsj <- rnorm(1,m=mj,sd=sqrt(vj))
    ## record the new ths[j] by replacing all th[i], i in Sj.
    th[idx] <- thsj
  }
  return(th) 
}  



sample.sig_py <- function(th,ydata)
{ ## sample
  ##   sig ~ p(sig | ...)
  ## returns: sig
  y<- ydata
  n<- length(y)
  s2 <- sum( (y-th)^2 )
  a1 <- a+0.5*n
  b1 <- b+0.5*s2
  s2.inv <- rgamma(1,shape=a1,rate=b1)
  return(1/sqrt(s2.inv))
}


fbar_py <- function(x,th,sig,M,sigma)
{ ## conditional draw F ~ p(F | th,sig,y) (approx -- will talk about this..)
  xgrid<-x
  n<- length(xgrid)
  nj <- table(th)              # counts
  ths <- as.numeric(names(nj)) # unique values
  k <- length(nj)
  fx <-(M+sigma*k)/(n+M)*dnorm(xgrid,m=m0,sd=sqrt(B0+sig^2))
  for(j in 1:k)
    fx <- fx + (nj[j]-sigma)/(n+M)*dnorm(xgrid,m=ths[j],sd=sig)
  return(fx)
}


ebar_py <- function(x,th,sig,M,sigma,ns)
{ ## conditional draw F ~ p(F | th,sig,y) (approx -- will talk about this..)
  xgrid<-x
  n<- ns
  nj <- table(th)              # counts
  ths <- as.numeric(names(nj)) # unique values
  k <- length(nj)
  fx <- (M+sigma*k)/(n+M)*pnorm(xgrid,m=m0,sd=sqrt(B0+sig^2))
  for(j in 1:k)
    fx <- fx + (nj[j]-sigma)/(n+M)*pnorm(xgrid,m=ths[j],sd=sig)
  return(fx)
}


gibbs_py <- function(n.iter=1000,M=1,sigma=0.5, n_s=50,ydata=y, burn=500){
  y<- ydata
  th <- init.DPk(y)              ## initialize th[1..n]
  sig <- sqrt( mean((y-th)^2))  ## and sig
  ## set up data structures to record imputed posterior draws..1
  xgrid <- seq(from=-2,to=15,length=n_s)
  #cgrid<- sort(y)
  cgrid<-  seq(from=-2,to=15,length=500)
  fgrid <- NULL     ## we will record imputed draws of f
  ecfgrid <- NULL     ## we will record imputed draws of f
  njlist <- NULL    ## record sizes of 8 largest clusters
  klist <- NULL
  sig_trace<- NULL
  ## start with a plot of a kernel density estimate of the data
  #d<-density(y,bw = "sj")
  true_dens<- 0.5*dnorm(xgrid, mean=1, sd= sqrt(0.2)) +0.5*dnorm(xgrid, mean=10,sd=sqrt(0.2)) 
 # plot(xgrid,true_dens,xlab="X",ylab="Y",bty="l",type="l",xlim=c(-2,15),ylim=c(0,1), main="")
  ## now the Gibbs sampler
  for(iter in 1:n.iter){
    th  <- sample.th_py(th,sig,M,sigma,y)
    sig <- sample.sig_py(th,y)
    th  <- sample.ths_py(th,sig,y)
    ## update running summaries #################
    sig_trace<- rbind(sig_trace,sig)
    f   <- fbar_py(xgrid,th,sig,M,sigma)
    ecf<- ebar_py(cgrid,th,sig,M,sigma, ns=length(y))
 #   lines(xgrid,f,col=iter,lty=3)
    fgrid <- rbind(fgrid,f)
    ecfgrid<- rbind(ecfgrid,ecf)
    nj <- table(th)                 # counts
    njlist <- rbind(njlist,sort(nj,decr=T)[1:8])
    klist <- c(klist,length(nj))
  }
  ## report summaries ###########################
  fbar <- apply(fgrid[(burn+1):n.iter,],2,mean)
 # lines(xgrid,fbar,lwd=3,col=2)
  ecfbar <- apply(ecfgrid[(burn+1):n.iter,],2,mean)
  true_ecdf<- 0.5*pnorm(cgrid, 1,sqrt(0.2)) +0.5*pnorm(cgrid, 10,sqrt(0.2)) 
  #plot(ecdf(y),xlab="X",ylab="Y",bty="l",xlim=c(min(y),max(y)),ylim=c(0,1), main="")
 # plot(cgrid, true_ecdf,xlab="X",ylab="Y",bty="l",type="l",ylim=c(0,1), main="")
#  lines(cgrid,ecfbar,lwd=3,col=2)
 # for(iter in 1:n.iter){ lines(cgrid,ecfgrid[iter,],col=iter,lty=3) }
  #ks_test<- ks.test(ecfbar,ecdf(y))
  njbar <- apply(njlist,2,mean,na.rm=T)
  cat("Average cluster sizes:\n",format(njbar),"\n")
  pk <- table(klist)/length(klist)
  cat("Posterior probs p(k): (row1 = k, row2 = p(k) \n ")
  #print(pk/sum(pk))
  #cgrid<- sort(y)
  #K_d_py <- max(abs(ecfbar - ecdf(y)(cgrid)))
  K_d_py <- max(abs(ecfbar -true_ecdf))
  return(list(fgrid=fgrid,klist=klist,njlist=njlist, sig=sig_trace, ecf=ecfbar, Kd=K_d_py))
}


sample<- gen_data(mu_vec=c(1,10),sigma_vec=rep(0.2, 2),ns=50,prob_vec=c(0.5,0.5))

#mc_py<- gibbs_py(n.iter=1000,M=1,sigma=0.5, n_s=50,ydata=sample, burn=500)


#xgrid <- seq(from=-2,to=15,length=100)
#cgrid<- sort(y)
#plot(ecdf(y),xlab="X",ylab="Y",bty="l",xlim=c(-2,15),ylim=c(0,1), main="")
#lines(cgrid,mc_py$ecf,lwd=3,col=2)
#K_d_py <- max(abs(mc_py$ecf - ecdf(y)(cgrid)))
#K_d_py


########################################################################PY approx

#M<-1
#sigma<-0.5

sample.th_py_ap1 <- function(th,sig,M,sigma,ydata)
{ ## sample
  ##    th[i]  ~ p(th_i | th[-i],sig,y)
  ## returns updated th vector
  y<- ydata
  n<- length(y)
  for(i in 1:n){
    ## unique values and counts
    nj <- table(th[-i])
    ths <- as.numeric(names(nj))
    k <- length(nj)
    ## likelihood
    fj <- dnorm(y[i], m=ths, s=sig)
    f0 <- dnorm(y[i], m=m0,   s=sqrt(B0+sig^2)) # q0
    #kn sigma/n  (n_i -sigma)/n
    pj <- c(fj*(nj- sigma) ,f0*(k*sigma))                       # p(s[i]=j | ...), j=1..k, k+1
    # pj_norm<- pj/sum(pj)
    s <- sample(1:(k+1), 1, prob=pj)          # sample s[i]
    if (s==k+1){ ## generate new th[i] value
      v.new <- 1.0/(1/B0 + 1/sig^2)
      m.new <- v.new*(1/B0*m0 + 1/sig^2*y[i])
      thi.new   <- rnorm(1,m=m.new,sd=sqrt(v.new))
      ths <- c(ths,thi.new)
    }   
    th[i] <- ths[s]
  }
  return(th) 
} 




fbar_py_ap1 <- function(x,th,sig,M,sigma)
{ ## conditional draw F ~ p(F | th,sig,y) (approx -- will talk about this..)
  xgrid<-x
  n<- length(xgrid)
  nj <- table(th)              # counts
  ths <- as.numeric(names(nj)) # unique values
  k <- length(nj)
  #kn sigma/n  (n_i -sigma)/n
  fx <-(sigma*k)/(n)*dnorm(xgrid,m=m0,sd=sqrt(B0+sig^2))
  for(j in 1:k)
    fx <- fx + (nj[j]-sigma)/(n)*dnorm(xgrid,m=ths[j],sd=sig)
  return(fx)
}



ebar_py_ap1 <- function(x,th,sig,M,sigma,n_ss)
{ ## conditional draw F ~ p(F | th,sig,y) (approx -- will talk about this..)
  xgrid<-x
  n<- n_ss
  nj <- table(th)              # counts
  ths <- as.numeric(names(nj)) # unique values
  k <- length(nj)
  fx <- (sigma*k)/(n)*pnorm(xgrid,m=m0,sd=sqrt(B0+sig^2))
  for(j in 1:k)
    fx <- fx + (nj[j]-sigma)/(n)*pnorm(xgrid,m=ths[j],sd=sig)
  return(fx)
}


gibbs_py_ap1 <- function(n.iter=1000,M=1,sigma=0.5, n_s=50,ydata=y, burn=500)
{
  y<- ydata
  th <- init.DPk(y)              ## initialize th[1..n]
  sig <- sqrt( mean((y-th)^2))  ## and sig
  ## set up data structures to record imputed posterior draws..1
  xgrid <- seq(from=-2,to=15,length=n_s)
  #cgrid<- sort(y)
  cgrid<-  seq(from=-2,to=15,length=500)
  fgrid <- NULL     ## we will record imputed draws of f
  ecfgrid <- NULL     ## we will record imputed draws of f
  njlist <- NULL    ## record sizes of 8 largest clusters
  klist <- NULL
  sig_trace<- NULL
  ## start with a plot of a kernel density estimate of the data
  #d<-density(y,bw = "sj")
  #plot(d,xlab="X",ylab="Y",bty="l",type="l",xlim=c(-2,15),ylim=c(0,1), main="")
  true_dens<- 0.5*dnorm(xgrid, mean=1, sd= sqrt(0.2)) +0.5*dnorm(xgrid, mean=10,sd=sqrt(0.2)) 
  #plot(xgrid,true_dens,xlab="X",ylab="Y",bty="l",type="l",xlim=c(-2,15),ylim=c(0,1), main="")
  ## now the Gibbs sampler
  for(iter in 1:n.iter){
    th  <- sample.th_py_ap1(th,sig,M,sigma,y)
    sig <- sample.sig_py(th,y)
    th  <- sample.ths_py(th,sig,y)
    ## update running summaries #################
    sig_trace<- rbind(sig_trace,sig)
    f   <- fbar_py_ap1(xgrid,th,sig,M,sigma)
    ecf<- ebar_py_ap1(cgrid,th,sig,M,sigma, n_ss=length(y))
  #  lines(xgrid,f,col=iter,lty=3)
    fgrid <- rbind(fgrid,f)
    ecfgrid<- rbind(ecfgrid,ecf)
    nj <- table(th)                 # counts
    njlist <- rbind(njlist,sort(nj,decr=T)[1:8])
    klist <- c(klist,length(nj))
  }
  ## report summaries ###########################
  fbar <- apply(fgrid[(burn+1):n.iter,],2,mean)
 # lines(xgrid,fbar,lwd=3,col=2)
  ecfbar <- apply(ecfgrid[(burn+1):n.iter,],2,mean)
  true_ecdf<- 0.5*pnorm(cgrid, 1,sqrt(0.2)) +0.5*pnorm(cgrid, 10,sqrt(0.2)) 
  #plot(ecdf(y),xlab="X",ylab="Y",bty="l",xlim=c(min(y),max(y)),ylim=c(0,1), main="")
  #plot(cgrid, true_ecdf,xlab="X",ylab="Y",bty="l",type="l",ylim=c(0,1), main="")
  #plot(ecdf(y),xlab="X",ylab="Y",bty="l",xlim=c(min(y),max(y)),ylim=c(0,1), main="")
  #lines(cgrid,ecfbar,lwd=3,col=2)
#  for(iter in 1:n.iter){ lines(cgrid,ecfgrid[iter,],col=iter,lty=3) }
  #ks_test<- ks.test(ecfbar,ecdf(y))
  njbar <- apply(njlist,2,mean,na.rm=T)
  cat("Average cluster sizes:\n",format(njbar),"\n")
  pk <- table(klist)/length(klist)
  cat("Posterior probs p(k): (row1 = k, row2 = p(k) \n ")
  print(pk/sum(pk))
#  cgrid<- sort(y)
 # K_d_py_ap <- max(abs(ecfbar - ecdf(y)(cgrid)))
  K_d_py_ap <- max(abs(ecfbar -true_ecdf))
  return(list(fgrid=fgrid,klist=klist,njlist=njlist, sig=sig_trace, ecf=ecfbar,Kd=K_d_py_ap))
}


#mc_py_ap1<- gibbs_py_ap1(n.iter=1000,M=1,sigma=0.5, n_s=50,ydata=sample)


######################################################NG
######################################################NG(2)

########################################################################PY approx

#initial parameters
#sigma<- 0.5
#tau<- 1

sample.th_ng_ap2 <- function(th,sig,sigma,tau,ydata)
{ ## sample
  ##    th[i]  ~ p(th_i | th[-i],sig,y)
  ## returns updated th vector
  y<- ydata
  n<- length(y)
  for(i in 1:n){
    ## unique values and counts
    nj <- table(th[-i])
    ths <- as.numeric(names(nj))
    k <- length(nj)
    ## likelihood
    
    fj <- dnorm(y[i], m=ths, s=sig)
    f0 <- dnorm(y[i], m=m0,   s=sqrt(B0+sig^2)) # q0
    #kn sigma/n  (n_i -sigma)/n
    #beta= tau n/kn^{1/alpha}
    beta<- (tau*n)/( k^(1/sigma))
    pj <- c(fj*(nj- sigma) ,f0*(k*sigma+ beta))                       # p(s[i]=j | ...), j=1..k, k+1
    # pj_norm<- pj/sum(pj)
    s <- sample(1:(k+1), 1, prob=pj)          # sample s[i]
    if (s==k+1){ ## generate new th[i] value
      v.new <- 1.0/(1/B0 + 1/sig^2)
      m.new <- v.new*(1/B0*m0 + 1/sig^2*y[i])
      thi.new   <- rnorm(1,m=m.new,sd=sqrt(v.new))
      ths <- c(ths,thi.new)
    }   
    th[i] <- ths[s]
  }
  return(th) 
} 




fbar_ng_ap2 <- function(x,th,sig,sigma, tau)
{ ## conditional draw F ~ p(F | th,sig,y) 
  xgrid<-x
  n<- length(xgrid) 
  nj <- table(th)              # counts
  ths <- as.numeric(names(nj)) # unique values
  k <- length(nj)
  #kn sigma/n  (n_i -sigma)/n
  beta<- (tau*n)/( k^(1/sigma))
  fx <-(k*sigma+ beta)/(n+beta)*dnorm(xgrid,m=m0,sd=sqrt(B0+sig^2))
  for(j in 1:k)
    fx <- fx + (nj[j]-sigma)/(n+beta)*dnorm(xgrid,m=ths[j],sd=sig)
  return(fx)
}



ebar_ng_ap2 <- function(x,th,sig,sigma, tau,n_ss)
{ ## conditional draw F ~ p(F | th,sig,y) 
  xgrid<-x
  n<- n_ss
  nj <- table(th)              # counts
  ths <- as.numeric(names(nj)) # unique values
  k <- length(nj)
  beta<- (tau*n)/( k^(1/sigma))
  fx <- (k*sigma+ beta)/(n+beta)*pnorm(xgrid,m=m0,sd=sqrt(B0+sig^2))
  for(j in 1:k)
    fx <- fx + (nj[j]-sigma)/(n+beta)*pnorm(xgrid,m=ths[j],sd=sig)
  return(fx)
}


gibbs_ng_ap2 <- function(n.iter=1000,sigma=0.5, tau=1, n_s=50, ydata=y, burn=500)
{
  y<- ydata
  th <- init.DPk(y)              ## initialize th[1..n]
  sig <- sqrt( mean((y-th)^2))  ## and sig
  ## set up data structures to record imputed posterior draws..1
  xgrid <- seq(from=-2,to=15,length=n_s)
  #cgrid<- sort(y)
  cgrid<-  seq(from=-2,to=15,length=500)
  fgrid <- NULL     ## we will record imputed draws of f
  ecfgrid <- NULL     ## we will record imputed draws of f
  njlist <- NULL    ## record sizes of 8 largest clusters
  klist <- NULL
  sig_trace<- NULL
  ## start with a plot of a kernel density estimate of the data
  #d<-density(y,bw = "sj")
  #plot(d,xlab="X",ylab="Y",bty="l",type="l",xlim=c(-2,15),ylim=c(0,1), main="")
  true_dens<- 0.5*dnorm(xgrid, mean=1, sd= sqrt(0.2)) +0.5*dnorm(xgrid, mean=10,sd=sqrt(0.2)) 
  #plot(xgrid,true_dens,xlab="X",ylab="Y",bty="l",type="l",xlim=c(-2,15),ylim=c(0,1), main="")
  
  ## now the Gibbs sampler
  for(iter in 1:n.iter){
    th  <- sample.th_ng_ap2(th,sig,sigma,tau,y)
    sig <- sample.sig_py(th,y)
    th  <- sample.ths_py(th,sig,y)
    ## update running summaries #################
    sig_trace<- rbind(sig_trace,sig)
    f   <- fbar_ng_ap2(xgrid,th,sig,sigma,tau)
    ecf<- ebar_ng_ap2(cgrid,th,sig,sigma,tau,n_ss=length(y))
   # lines(xgrid,f,col=iter,lty=3)
    fgrid <- rbind(fgrid,f)
    ecfgrid<- rbind(ecfgrid,ecf)
    nj <- table(th)                 # counts
    njlist <- rbind(njlist,sort(nj,decr=T)[1:8])
    klist <- c(klist,length(nj))
  }
  ## report summaries ###########################
  fbar <- apply(fgrid[(burn+1):n.iter,],2,mean)
 # lines(xgrid,fbar,lwd=3,col=2)
  ecfbar <- apply(ecfgrid[(burn+1):n.iter,],2,mean)
 # plot(ecdf(y),xlab="X",ylab="Y",bty="l",xlim=c(min(y),max(y)),ylim=c(0,1), main="")
  true_ecdf<- 0.5*pnorm(cgrid, 1,sqrt(0.2)) +0.5*pnorm(cgrid, 10,sqrt(0.2)) 
  #plot(ecdf(y),xlab="X",ylab="Y",bty="l",xlim=c(min(y),max(y)),ylim=c(0,1), main="")
 # plot(cgrid, true_ecdf,xlab="X",ylab="Y",bty="l",type="l",ylim=c(0,1), main="")
  #lines(cgrid,ecfbar,lwd=3,col=2)
#  for(iter in 1:n.iter){ lines(cgrid,ecfgrid[iter,],col=iter,lty=3) }
  #ks_test<- ks.test(ecfbar,ecdf(y))
  njbar <- apply(njlist,2,mean,na.rm=T)
  cat("Average cluster sizes:\n",format(njbar),"\n")
  pk <- table(klist)/length(klist)
  cat("Posterior probs p(k): (row1 = k, row2 = p(k) \n ")
  #print(pk/sum(pk))
  #cgrid<- sort(y)
  
  K_d_ng_ap2 <- max(abs(ecfbar -true_ecdf))
  return(list(fgrid=fgrid,klist=klist,njlist=njlist, sig=sig_trace, ecf=ecfbar, Kd=K_d_ng_ap2))
}

#mc_ng_ap2<- gibbs_ng_ap2(n.iter=1000,sigma=0.5, tau=1, n_s=50, ydata=sample,burn=500)
# cgrid<- sort(y)
# plot(ecdf(y),xlab="X",ylab="Y",bty="l",xlim=c(-2,15),ylim=c(0,1), main="")
# lines(cgrid,mc_ng_ap2$ecf,lwd=3,col=2)
# 
# cgrid<- sort(y)
# plot(ecdf(y),xlab="X",ylab="Y",bty="l",xlim=c(-2,15),ylim=c(0,1), main="")
# lines(cgrid,mc_py_ap1$ecf,lwd=3,col=2)
# lines(cgrid,mc_py$ecf,lwd=3,col=3)
# lines(cgrid, mc_ng_ap2$ecf,lwd=3,col=4)

# 
# K_d_ng <- max(abs( mc_ng_ap2$ecf - ecdf(y)(cgrid)))
# K_d_ng
# 




##############################################n=50##########################################################################

set.seed(123)
sample<- gen_data(mu_vec=c(1,10),sigma_vec=rep(0.2, 2),ns=50,prob_vec=c(0.5,0.5))
save(sample, file = "sample_data_n50.Rda")


#hist(sample,probability=TRUE,breaks=20,col=grey(.9))
#plot(density(sample,bw = "sj"))
y_data <- sample
n <- length(y_data)
## hyperparameters
a  <- 1;   b <- 1   # 1/sig ~ Ga(a,b)
m0 <- 0;  B0 <- 1  # G0 = N(m0,B0)
F_dmat_PY<- matrix(NA, nrow =3*3*1000, ncol =53)
F_dmat_PY1<- matrix(NA, nrow =3*3*1000, ncol =53)
F_dmat_NG2<- matrix(NA, nrow =3*3*1000, ncol =53)
E_dmat_PY<- matrix(NA, nrow =3*3, ncol =53)
E_dmat_PY1<- matrix(NA, nrow =3*3, ncol =53)
E_dmat_NG2<- matrix(NA, nrow =3*3, ncol =53)

K_dmat<- matrix(NA, nrow =3*3, ncol =5)
sigma_vec<- c(0.25, 0.5, 0.75)
M_vec<- c(1,3,10)
#ns_vec<- c(50,100,200,500) 
i<-1

for(m in 1:length(M_vec)){
  for(j in 1:length(sigma_vec)){
    K_dmat[i,1]<- M_vec[m]
    K_dmat[i,2]<- sigma_vec[j]
    mc_py<- gibbs_py(n.iter=1000,M=M_vec[m],sigma=sigma_vec[j], n_s=length(y_data),ydata=y_data)
    mc_py_ap1<- gibbs_py_ap1(n.iter=1000,M=M_vec[m],sigma=sigma_vec[j], n_s=length(y_data),ydata=y_data)
    mc_ng_ap2<- gibbs_ng_ap2(n.iter=1000,sigma=sigma_vec[j], tau=M_vec[m], n_s=length(y_data), ydata=y_data)
    K_dmat[i,3]<- mc_py$Kd
    K_dmat[i,4]<- mc_py_ap1$Kd
    K_dmat[i,5]<- mc_ng_ap2$Kd
    ind_1<- (i-1)*1000 +1
    ind_2<- i*1000
    F_dmat_PY[ind_1: ind_2,1:50]<- mc_py$fgrid
    F_dmat_PY1[ind_1: ind_2,1:50]<- mc_py_ap1$fgrid
    F_dmat_NG2[ind_1: ind_2,1:50]<- mc_ng_ap2$fgrid
    F_dmat_PY[ind_1: ind_2,51]<- F_dmat_PY1[ind_1: ind_2,51]<-  F_dmat_NG2[ind_1: ind_2,51]<- M_vec[m]
    F_dmat_PY[ind_1: ind_2,52]<- F_dmat_PY1[ind_1: ind_2,52]<-  F_dmat_NG2[ind_1: ind_2,52]<- sigma_vec[j]
    F_dmat_PY[ind_1: ind_2,53]<- F_dmat_PY1[ind_1: ind_2,53]<-  F_dmat_NG2[ind_1: ind_2,53]<- length(y_data)
    E_dmat_PY[i,1:50]<- mc_py$ecf
    E_dmat_PY1[i,1:50]<- mc_py_ap1$ecf
    E_dmat_NG2[i,1:50]<- mc_ng_ap2$ecf
    E_dmat_PY[i,51]<- E_dmat_PY1[i,51]<-  E_dmat_NG2[i,51]<- M_vec[m]
    E_dmat_PY[i,52]<- E_dmat_PY1[i,52]<-  E_dmat_NG2[i,52]<- sigma_vec[j]
    E_dmat_PY[i,53]<- E_dmat_PY1[i,53]<-  E_dmat_NG2[i,53]<- length(y_data)
    i<- i+1
  }
}
k_mat_50<- as.data.frame(K_dmat)
k_mat_50$n<- 50
names(k_mat_50)<- c("alpha","sigma","PY","PY_A1","NG_A2","N")
write.csv(k_mat_50, file = "K_50.csv")
formattable(k_mat_50)

F_mat_50_PY<- as.data.frame(F_dmat_PY)
F_mat_50_PY1<- as.data.frame(F_dmat_PY1)
F_mat_50_NG2<- as.data.frame(F_dmat_NG2)
save(F_mat_50_PY, file = "FmatPY_n50.Rda")
save(F_mat_50_PY1, file = "FmatPY1_n50.Rda")
save(F_mat_50_NG2, file = "FmatNG2_n50.Rda")

E_mat_50_PY<- as.data.frame(E_dmat_PY)
E_mat_50_PY1<- as.data.frame(E_dmat_PY1)
E_mat_50_NG2<- as.data.frame(E_dmat_NG2)
save(E_mat_50_PY, file = "EmatPY_n50.Rda")
save(E_mat_50_PY1, file = "EmatPY1_n50.Rda")
save(E_mat_50_NG2, file = "EmatNG2_n50.Rda")


##############################################n=200##########################################################################
set.seed(123)
sample200<- gen_data(mu_vec=c(1,10),sigma_vec=rep(0.2, 2),ns=200,prob_vec=c(0.5,0.5))
save(sample200, file = "sample_data_n200.Rds")

#hist(sample,probability=TRUE,breaks=20,col=grey(.9))
#plot(density(sample,bw = "sj"))
y_data <- sample200
n <- length(y_data)
## hyperparameters
a  <- 1;   b <- 1   # 1/sig ~ Ga(a,b)
m0 <- 0;  B0 <- 1  # G0 = N(m0,B0)
F_dmat_PY<- matrix(NA, nrow =10000, ncol =203)
F_dmat_PY1<- matrix(NA, nrow =10000, ncol =203)
F_dmat_NG2<- matrix(NA, nrow =10000, ncol =203)
E_dmat_PY<- matrix(NA, nrow =1, ncol =503)
E_dmat_PY1<- matrix(NA, nrow =1, ncol =503)
E_dmat_NG2<- matrix(NA, nrow =1, ncol =503)

K_dmat<- matrix(NA, nrow =3*3, ncol =5)
sigma_vec<- c(0.5)
M_vec<- c(1)
i<-1
niter<- 10000
burnin<- 2000
for(m in 1:length(M_vec)){
  for(j in 1:length(sigma_vec)){
    n_grid<-500
    n_sm<- length(y_data)
    K_dmat[i,1]<- M_vec[m]
    K_dmat[i,2]<- sigma_vec[j]
    mc_py<- gibbs_py(n.iter=niter,M=M_vec[m],sigma=sigma_vec[j], n_s=length(y_data),ydata=y_data, burn=burnin)
    mc_py_ap1<- gibbs_py_ap1(n.iter=niter,M=M_vec[m],sigma=sigma_vec[j], n_s=length(y_data),ydata=y_data, burn=burnin)
    mc_ng_ap2<- gibbs_ng_ap2(n.iter=niter,sigma=sigma_vec[j], tau=M_vec[m], n_s=length(y_data), ydata=y_data, burn=burnin)
    K_dmat[i,3]<- mc_py$Kd
    K_dmat[i,4]<- mc_py_ap1$Kd
    K_dmat[i,5]<- mc_ng_ap2$Kd
    ind_1<- (i-1)*niter +1
    ind_2<- i*niter
    F_dmat_PY[ind_1: ind_2,1:n_sm]<- mc_py$fgrid
    F_dmat_PY1[ind_1: ind_2,1:n_sm]<- mc_py_ap1$fgrid
    F_dmat_NG2[ind_1: ind_2,1:n_sm]<- mc_ng_ap2$fgrid
    F_dmat_PY[ind_1: ind_2,n_sm+1]<- F_dmat_PY1[ind_1: ind_2,n_sm+1]<-  F_dmat_NG2[ind_1: ind_2,n_sm+1]<- M_vec[m]
    F_dmat_PY[ind_1: ind_2,n_sm+2]<- F_dmat_PY1[ind_1: ind_2,n_sm+2]<-  F_dmat_NG2[ind_1: ind_2,n_sm+2]<- sigma_vec[j]
    F_dmat_PY[ind_1: ind_2,n_sm+3]<- F_dmat_PY1[ind_1: ind_2,n_sm+3]<-  F_dmat_NG2[ind_1: ind_2,n_sm+3]<- length(y_data)
    E_dmat_PY[i,1:n_grid]<- mc_py$ecf
    E_dmat_PY1[i,1:n_grid]<- mc_py_ap1$ecf
    E_dmat_NG2[i,1:n_grid]<- mc_ng_ap2$ecf
    E_dmat_PY[i,n_grid+1]<- E_dmat_PY1[i,n_grid+1]<-  E_dmat_NG2[i,n_grid+1]<- M_vec[m]
    E_dmat_PY[i,n_grid+2]<- E_dmat_PY1[i,n_grid+2]<-  E_dmat_NG2[i,n_grid+2]<- sigma_vec[j]
    E_dmat_PY[i,n_grid+3]<- E_dmat_PY1[i,n_grid+3]<-  E_dmat_NG2[i,n_grid+3]<- length(y_data)
    i<- i+1
  }
}



k_mat_200<- as.data.frame(K_dmat)
k_mat_200$n<- length(y_data)
names(k_mat_200)<- c("alpha","sigma","PY","PY_A1","NG_A2","N")
write.csv(k_mat_200, file = "K_200.csv")
formattable(k_mat_200)

F_mat_200_PY<- as.data.frame(F_dmat_PY)
F_mat_200_PY1<- as.data.frame(F_dmat_PY1)
F_mat_200_NG2<- as.data.frame(F_dmat_NG2)
save(F_mat_200_PY, file = "FmatPY_n200.Rds")
save(F_mat_200_PY1, file = "FmatPY1_n200.Rds")
save(F_mat_200_NG2, file = "FmatNG2_n200.Rds")

E_mat_200_PY<- as.data.frame(E_dmat_PY)
E_mat_200_PY1<- as.data.frame(E_dmat_PY1)
E_mat_200_NG2<- as.data.frame(E_dmat_NG2)
save(E_mat_200_PY, file = "EmatPY_n200.Rds")
save(E_mat_200_PY1, file = "EmatPY1_n200.Rds")
save(E_mat_200_NG2, file = "EmatNG2_n200.Rds")

xgrid<-  seq(from=-2,to=15,length=length(y_data)) 
cgrid<-  seq(from=-2,to=15,length=500)
true_ecdf<- 0.5*pnorm(cgrid, 1,sqrt(0.2)) +0.5*pnorm(cgrid, 10,sqrt(0.2)) 

plot(cgrid, true_ecdf,xlab="X",ylab="Y",bty="l",type="l",ylim=c(0,1), main="")
lines(cgrid,E_mat_200_PY[1,1:500],lwd=3,col="red")
lines(cgrid,E_mat_200_PY1[1,1:500],lwd=3,col="blue")
lines(cgrid,E_mat_200_NG2[1,1:500],lwd=3,col="green")


true_dens<- 0.5*dnorm(xgrid, mean=1, sd= sqrt(0.2)) +0.5*dnorm(xgrid, mean=10,sd=sqrt(0.2)) 
plot(xgrid,true_dens,xlab="X",ylab="Y",bty="l",type="l",xlim=c(-2,15),ylim=c(0,1), main="")
for(iter in burnin:niter){ lines(xgrid, F_mat_200_PY[iter,1:200],col=iter,lty=3) }

true_dens<- 0.5*dnorm(xgrid, mean=1, sd= sqrt(0.2)) +0.5*dnorm(xgrid, mean=10,sd=sqrt(0.2)) 
plot(xgrid,true_dens,xlab="X",ylab="Y",bty="l",type="l",xlim=c(-2,15),ylim=c(0,1), main="")
for(iter in burnin:niter){ lines(xgrid, F_mat_200_PY2[iter,1:200],col=iter,lty=3) }

true_dens<- 0.5*dnorm(xgrid, mean=1, sd= sqrt(0.2)) +0.5*dnorm(xgrid, mean=10,sd=sqrt(0.2)) 
plot(xgrid,true_dens,xlab="X",ylab="Y",bty="l",type="l",xlim=c(-2,15),ylim=c(0,1), main="")
for(iter in burnin:niter){ lines(xgrid, F_mat_200_NG2[iter,1:200],col=iter,lty=3) }








##############################################n=500##########################################################################
set.seed(1234)
sample500<- gen_data(mu_vec=c(1,10),sigma_vec=rep(0.2, 2),ns=500,prob_vec=c(0.5,0.5))
save(sample500, file = "sample_data_n500.Rds")

#hist(sample,probability=TRUE,breaks=20,col=grey(.9))
#plot(density(sample,bw = "sj"))
y_data <- sample500
n <- length(y_data)
## hyperparameters
a  <- 1;   b <- 1   # 1/sig ~ Ga(a,b)
m0 <- 0;  B0 <- 1  # G0 = N(m0,B0)
F_dmat_PY<- matrix(NA, nrow =10000, ncol =503)
F_dmat_PY1<- matrix(NA, nrow =10000, ncol =503)
F_dmat_NG2<- matrix(NA, nrow =10000, ncol =503)
E_dmat_PY<- matrix(NA, nrow =1, ncol =503)
E_dmat_PY1<- matrix(NA, nrow =1, ncol =503)
E_dmat_NG2<- matrix(NA, nrow =1, ncol =503)

K_dmat<- matrix(NA, nrow =3*3, ncol =5)
sigma_vec<- c(0.5)
M_vec<- c(1)
#ns_vec<- c(50,100,200,500) 
i<-1

for(m in 1:length(M_vec)){
  for(j in 1:length(sigma_vec)){
    n_grid<-500
    n_sm<-500
    n_sm<- length(y_data)
    K_dmat[i,1]<- M_vec[m]
    K_dmat[i,2]<- sigma_vec[j]
    mc_py<- gibbs_py(n.iter=10000,M=M_vec[m],sigma=sigma_vec[j], n_s=length(y_data),ydata=y_data, burn=2000)
    mc_py_ap1<- gibbs_py_ap1(n.iter=10000,M=M_vec[m],sigma=sigma_vec[j], n_s=length(y_data),ydata=y_data, burn=2000)
    mc_ng_ap2<- gibbs_ng_ap2(n.iter=10000,sigma=sigma_vec[j], tau=M_vec[m], n_s=length(y_data), ydata=y_data, burn=2000)
    K_dmat[i,3]<- mc_py$Kd
    K_dmat[i,4]<- mc_py_ap1$Kd
    K_dmat[i,5]<- mc_ng_ap2$Kd
    ind_1<- (i-1)*10000 +1
    ind_2<- i*10000
    F_dmat_PY[ind_1: ind_2,1:n_sm]<- mc_py$fgrid
    F_dmat_PY1[ind_1: ind_2,1:n_sm]<- mc_py_ap1$fgrid
    F_dmat_NG2[ind_1: ind_2,1:n_sm]<- mc_ng_ap2$fgrid
    F_dmat_PY[ind_1: ind_2,n_sm+1]<- F_dmat_PY1[ind_1: ind_2,n_sm+1]<-  F_dmat_NG2[ind_1: ind_2,n_sm+1]<- M_vec[m]
    F_dmat_PY[ind_1: ind_2,n_sm+2]<- F_dmat_PY1[ind_1: ind_2,n_sm+2]<-  F_dmat_NG2[ind_1: ind_2,n_sm+2]<- sigma_vec[j]
    F_dmat_PY[ind_1: ind_2,n_sm+3]<- F_dmat_PY1[ind_1: ind_2,n_sm+3]<-  F_dmat_NG2[ind_1: ind_2,n_sm+3]<- length(y_data)
    E_dmat_PY[i,1:n_grid]<- mc_py$ecf
    E_dmat_PY1[i,1:n_grid]<- mc_py_ap1$ecf
    E_dmat_NG2[i,1:n_grid]<- mc_ng_ap2$ecf
    E_dmat_PY[i,n_grid+1]<- E_dmat_PY1[i,n_grid+1]<-  E_dmat_NG2[i,n_grid+1]<- M_vec[m]
    E_dmat_PY[i,n_grid+2]<- E_dmat_PY1[i,n_grid+2]<-  E_dmat_NG2[i,n_grid+2]<- sigma_vec[j]
    E_dmat_PY[i,n_grid+3]<- E_dmat_PY1[i,n_grid+3]<-  E_dmat_NG2[i,n_grid+3]<- length(y_data)
    i<- i+1
  }
}



k_mat_500<- as.data.frame(K_dmat)
k_mat_500$n<- 500
names(k_mat_500)<- c("alpha","sigma","PY","PY_A1","NG_A2","N")
write.csv(k_mat_500, file = "K_500.csv")
formattable(k_mat_500)

F_mat_500_PY<- as.data.frame(F_dmat_PY)
F_mat_500_PY1<- as.data.frame(F_dmat_PY1)
F_mat_500_NG2<- as.data.frame(F_dmat_NG2)
save(F_mat_500_PY, file = "FmatPY_n500.Rds")
save(F_mat_500_PY1, file = "FmatPY1_n500.Rds")
save(F_mat_500_NG2, file = "FmatNG2_n500.Rds")
F_mat_500_PY_short<- F_mat_500_PY[burnin:niter,]
save(F_mat_500_PY_short, file = "FmatPY_n500_burn.Rds")
F_mat_500_PY1_short<- F_mat_500_PY1[burnin:niter,]
save(F_mat_500_PY1_short, file = "FmatPY1_n500_burn.Rds")
F_mat_500_NG2_short<- F_mat_500_NG2[burnin:niter,]
save(F_mat_500_NG2_short, file = "FmatNG2_n500_burn.Rds")



E_mat_500_PY<- as.data.frame(E_dmat_PY)
E_mat_500_PY1<- as.data.frame(E_dmat_PY1)
E_mat_500_NG2<- as.data.frame(E_dmat_NG2)
save(E_mat_500_PY, file = "EmatPY_n500.Rds")
save(E_mat_500_PY1, file = "EmatPY1_n500.Rds")
save(E_mat_500_NG2, file = "EmatNG2_n500.Rds")

xgrid<-  seq(from=-2,to=15,length=500) 
cgrid<-  seq(from=-2,to=15,length=500)
true_ecdf<- 0.5*pnorm(cgrid, 1,sqrt(0.2)) +0.5*pnorm(cgrid, 10,sqrt(0.2)) 

plot(cgrid, true_ecdf,xlab="X",ylab="Y",bty="l",type="l",ylim=c(0,1), main="")
lines(cgrid,E_mat_500_PY[1,1:500],lwd=3,col="red")
lines(cgrid,E_mat_500_PY1[1,1:500],lwd=3,col="blue")
lines(cgrid,E_mat_500_NG2[1,1:500],lwd=3,col="green")


true_dens<- 0.5*dnorm(xgrid, mean=1, sd= sqrt(0.2)) +0.5*dnorm(xgrid, mean=10,sd=sqrt(0.2)) 
plot(xgrid,true_dens,xlab="X",ylab="Y",bty="l",type="l",xlim=c(-2,15),ylim=c(0,1), main="")
for(iter in 2000:10000){ lines(xgrid, F_mat_500_PY[iter,1:500],col=iter,lty=3) }

true_dens<- 0.5*dnorm(xgrid, mean=1, sd= sqrt(0.2)) +0.5*dnorm(xgrid, mean=10,sd=sqrt(0.2)) 
plot(xgrid,true_dens,xlab="X",ylab="Y",bty="l",type="l",xlim=c(-2,15),ylim=c(0,1), main="")
for(iter in 2000:10000){ lines(xgrid, F_mat_500_PY2[iter,1:500],col=iter,lty=3) }

true_dens<- 0.5*dnorm(xgrid, mean=1, sd= sqrt(0.2)) +0.5*dnorm(xgrid, mean=10,sd=sqrt(0.2)) 
plot(xgrid,true_dens,xlab="X",ylab="Y",bty="l",type="l",xlim=c(-2,15),ylim=c(0,1), main="")
for(iter in 2000:10000){ lines(xgrid, F_mat_500_NG2[iter,1:500],col=iter,lty=3) }

#############################################Load object#########################################################################


load_object <- function(file) {
  tmp <- new.env()
  load(file = file, envir = tmp)
  tmp[[ls(tmp)[1]]]
}

F1<- load_object("FmatPY1_n200.Rds")
true_dens<- 0.5*dnorm(xgrid, mean=1, sd= sqrt(0.2)) +0.5*dnorm(xgrid, mean=10,sd=sqrt(0.2)) 
plot(xgrid,true_dens,xlab="X",ylab="Y",bty="l",type="l",xlim=c(-2,15),ylim=c(0,1), main="")
for(iter in 1:1000){ lines(xgrid, F1[iter,1:200],col=iter,lty=3) }
