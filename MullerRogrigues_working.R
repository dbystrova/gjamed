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
sample<- gen_data(mu_vec=c(1,10),sigma_vec=rep(0.2, 2),ns=100,prob_vec=c(0.5,0.5))
hist(sample,probability=TRUE,breaks=20,col=grey(.9))
plot(density(sample,bw = "sj"))



## Dirichlet process BM ############################################################


#pdf("Output.pdf")
y <- sample
n <- length(y)
## hyperparameters
a  <- 1;   b <- 1   # 1/sig ~ Ga(a,b)
m0 <- 0;  B0 <- 1  # G0 = N(m0,B0)
M  <- 1

##########INIT##########
init.DPk <- function()
{ ## inital EDA estimate of th[1..n]
  ## cluster data, and cut at height H=10, to get 10 clusters
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

mc<- gibbs()




plot(mc$sig)
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
M<-1
sigma<-0.25

sample.th_py <- function(th,sig)
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


sample.ths_py <- function(th,sig)
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



sample.sig_py <- function(th)
{ ## sample
  ##   sig ~ p(sig | ...)
  ## returns: sig
  s2 <- sum( (y-th)^2 )
  a1 <- a+0.5*n
  b1 <- b+0.5*s2
  s2.inv <- rgamma(1,shape=a1,rate=b1)
  return(1/sqrt(s2.inv))
}


fbar_py <- function(x,th,sig)
{ ## conditional draw F ~ p(F | th,sig,y) (approx -- will talk about this..)
  xgrid<-x
  nj <- table(th)              # counts
  ths <- as.numeric(names(nj)) # unique values
  k <- length(nj)
  fx <-(M+sigma*k)/(n+M)*dnorm(xgrid,m=m0,sd=sqrt(B0+sig^2))
  for(j in 1:k)
    fx <- fx + (nj[j]-sigma)/(n+M)*dnorm(xgrid,m=ths[j],sd=sig)
  return(fx)
}


ebar_py <- function(x,th,sig)
{ ## conditional draw F ~ p(F | th,sig,y) (approx -- will talk about this..)
  xgrid<-x
  nj <- table(th)              # counts
  ths <- as.numeric(names(nj)) # unique values
  k <- length(nj)
  fx <- (M+sigma*k)/(n+M)*pnorm(xgrid,m=m0,sd=sqrt(B0+sig^2))
  for(j in 1:k)
    fx <- fx + (nj[j]-sigma)/(n+M)*pnorm(xgrid,m=ths[j],sd=sig)
  return(fx)
}


gibbs_py <- function(n.iter=1000)
{
  th <- init.DPk()              ## initialize th[1..n]
  sig <- sqrt( mean((y-th)^2))  ## and sig
  ## set up data structures to record imputed posterior draws..1
  xgrid <- seq(from=-2,to=15,length=100)
  cgrid<- sort(y)
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
    th  <- sample.th_py(th,sig)
    sig <- sample.sig_py(th)
    th  <- sample.ths_py(th,sig)
    ## update running summaries #################
    sig_trace<- rbind(sig_trace,sig)
    f   <- fbar_py(xgrid,th,sig)
    ecf<- ebar_py(cgrid,th,sig)
    lines(xgrid,f,col=iter,lty=3)
    fgrid <- rbind(fgrid,f)
    ecfgrid<- rbind(ecfgrid,ecf)
    nj <- table(th)                 # counts
    njlist <- rbind(njlist,sort(nj,decr=T)[1:8])
    klist <- c(klist,length(nj))
  }
  ## report summaries ###########################
  fbar <- apply(fgrid[floor(n.iter/2):n.iter,],2,mean)
  lines(xgrid,fbar,lwd=3,col=2)
  ecfbar <- apply(ecfgrid[floor(n.iter/2):n.iter,],2,mean)
  plot(ecdf(y),xlab="X",ylab="Y",bty="l",xlim=c(min(y),max(y)),ylim=c(0,1), main="")
  lines(cgrid,ecfbar,lwd=3,col=2)
  for(iter in 1:n.iter){ lines(cgrid,ecfgrid[iter,],col=iter,lty=3) }
  #ks_test<- ks.test(ecfbar,ecdf(y))
  njbar <- apply(njlist,2,mean,na.rm=T)
  cat("Average cluster sizes:\n",format(njbar),"\n")
  pk <- table(klist)/length(klist)
  cat("Posterior probs p(k): (row1 = k, row2 = p(k) \n ")
  print(pk/sum(pk))
  return(list(fgrid=fgrid,klist=klist,njlist=njlist, sig=sig_trace, ecf=ecfbar))
}

mc_py<- gibbs_py()
#xgrid <- seq(from=-2,to=15,length=100)
cgrid<- sort(y)
plot(ecdf(y),xlab="X",ylab="Y",bty="l",xlim=c(-2,15),ylim=c(0,1), main="")
lines(cgrid,mc_py$ecf,lwd=3,col=2)
K_d_py <- max(abs(mc_py$ecf - ecdf(y)(cgrid)))
K_d_py


########################################################################PY approx

M<-1
sigma<-0.5




sample.th_py_ap1 <- function(th,sig)
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




fbar_py_ap1 <- function(x,th,sig)
{ ## conditional draw F ~ p(F | th,sig,y) (approx -- will talk about this..)
  xgrid<-x
  nj <- table(th)              # counts
  ths <- as.numeric(names(nj)) # unique values
  k <- length(nj)
  #kn sigma/n  (n_i -sigma)/n
  fx <-(sigma*k)/(n)*dnorm(xgrid,m=m0,sd=sqrt(B0+sig^2))
  for(j in 1:k)
    fx <- fx + (nj[j]-sigma)/(n)*dnorm(xgrid,m=ths[j],sd=sig)
  return(fx)
}



ebar_py_ap1 <- function(x,th,sig)
{ ## conditional draw F ~ p(F | th,sig,y) (approx -- will talk about this..)
  xgrid<-x
  nj <- table(th)              # counts
  ths <- as.numeric(names(nj)) # unique values
  k <- length(nj)
  fx <- (sigma*k)/(n)*pnorm(xgrid,m=m0,sd=sqrt(B0+sig^2))
  for(j in 1:k)
    fx <- fx + (nj[j]-sigma)/(n)*pnorm(xgrid,m=ths[j],sd=sig)
  return(fx)
}


gibbs_py_ap1 <- function(n.iter=1000)
{
  th <- init.DPk()              ## initialize th[1..n]
  sig <- sqrt( mean((y-th)^2))  ## and sig
  ## set up data structures to record imputed posterior draws..1
  xgrid <- seq(from=-2,to=15,length=100)
  cgrid<- sort(y)
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
    th  <- sample.th_py_ap1(th,sig)
    sig <- sample.sig_py(th)
    th  <- sample.ths_py(th,sig)
    ## update running summaries #################
    sig_trace<- rbind(sig_trace,sig)
    f   <- fbar_py_ap1(xgrid,th,sig)
    ecf<- ebar_py_ap1(cgrid,th,sig)
    lines(xgrid,f,col=iter,lty=3)
    fgrid <- rbind(fgrid,f)
    ecfgrid<- rbind(ecfgrid,ecf)
    nj <- table(th)                 # counts
    njlist <- rbind(njlist,sort(nj,decr=T)[1:8])
    klist <- c(klist,length(nj))
  }
  ## report summaries ###########################
  fbar <- apply(fgrid[floor(n.iter/2):n.iter,],2,mean)
  lines(xgrid,fbar,lwd=3,col=2)
  ecfbar <- apply(ecfgrid[floor(n.iter/2):n.iter,],2,mean)
  plot(ecdf(y),xlab="X",ylab="Y",bty="l",xlim=c(min(y),max(y)),ylim=c(0,1), main="")
  lines(cgrid,ecfbar,lwd=3,col=2)
  for(iter in 1:n.iter){ lines(cgrid,ecfgrid[iter,],col=iter,lty=3) }
  #ks_test<- ks.test(ecfbar,ecdf(y))
  njbar <- apply(njlist,2,mean,na.rm=T)
  cat("Average cluster sizes:\n",format(njbar),"\n")
  pk <- table(klist)/length(klist)
  cat("Posterior probs p(k): (row1 = k, row2 = p(k) \n ")
  print(pk/sum(pk))
  return(list(fgrid=fgrid,klist=klist,njlist=njlist, sig=sig_trace, ecf=ecfbar))
}

mc_py_ap1<- gibbs_py_ap1()
xgrid <- seq(from=-2,to=15,length=100)
plot(ecdf(y),xlab="X",ylab="Y",bty="l",xlim=c(-2,15),ylim=c(0,1), main="")
lines(xgrid,mc_py_ap1$ecf,lwd=3,col=2)

cgrid<- sort(y)
plot(ecdf(y),xlab="X",ylab="Y",bty="l",xlim=c(min(y),max(y)),ylim=c(0,1), main="")
lines(cgrid,mc_py_ap1$ecf,lwd=3,col=2)
lines(cgrid,mc_py$ecf,lwd=3,col=4)


K_d_py_1 <- max(abs(mc_py_ap1$ecf - ecdf(y)(cgrid)))
K_d_py_1


######################################################NG




######################################################NG(2)

########################################################################PY approx

#initial parameters
sigma<- 0.25
tau<- 1

sample.th_ng_ap2 <- function(th,sig)
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




fbar_ng_ap2 <- function(x,th,sig)
{ ## conditional draw F ~ p(F | th,sig,y) 
  xgrid<-x
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



ebar_ng_ap2 <- function(x,th,sig)
{ ## conditional draw F ~ p(F | th,sig,y) 
  xgrid<-x
  nj <- table(th)              # counts
  ths <- as.numeric(names(nj)) # unique values
  k <- length(nj)
  beta<- (tau*n)/( k^(1/sigma))
  fx <- (k*sigma+ beta)/(n+beta)*pnorm(xgrid,m=m0,sd=sqrt(B0+sig^2))
  for(j in 1:k)
    fx <- fx + (nj[j]-sigma)/(n+beta)*pnorm(xgrid,m=ths[j],sd=sig)
  return(fx)
}


gibbs_ng_ap2 <- function(n.iter=1000)
{
  th <- init.DPk()              ## initialize th[1..n]
  sig <- sqrt( mean((y-th)^2))  ## and sig
  ## set up data structures to record imputed posterior draws..1
  xgrid <- seq(from=-2,to=15,length=100)
  cgrid<- sort(y)
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
    th  <- sample.th_ng_ap2(th,sig)
    sig <- sample.sig_py(th)
    th  <- sample.ths_py(th,sig)
    ## update running summaries #################
    sig_trace<- rbind(sig_trace,sig)
    f   <- fbar_ng_ap2(xgrid,th,sig)
    ecf<- ebar_ng_ap2(cgrid,th,sig)
    lines(xgrid,f,col=iter,lty=3)
    fgrid <- rbind(fgrid,f)
    ecfgrid<- rbind(ecfgrid,ecf)
    nj <- table(th)                 # counts
    njlist <- rbind(njlist,sort(nj,decr=T)[1:8])
    klist <- c(klist,length(nj))
  }
  ## report summaries ###########################
  fbar <- apply(fgrid[floor(n.iter/2):n.iter,],2,mean)
  lines(xgrid,fbar,lwd=3,col=2)
  ecfbar <- apply(ecfgrid[floor(n.iter/2):n.iter,],2,mean)
  plot(ecdf(y),xlab="X",ylab="Y",bty="l",xlim=c(min(y),max(y)),ylim=c(0,1), main="")
  lines(cgrid,ecfbar,lwd=3,col=2)
  for(iter in 1:n.iter){ lines(cgrid,ecfgrid[iter,],col=iter,lty=3) }
  #ks_test<- ks.test(ecfbar,ecdf(y))
  njbar <- apply(njlist,2,mean,na.rm=T)
  cat("Average cluster sizes:\n",format(njbar),"\n")
  pk <- table(klist)/length(klist)
  cat("Posterior probs p(k): (row1 = k, row2 = p(k) \n ")
  print(pk/sum(pk))
  return(list(fgrid=fgrid,klist=klist,njlist=njlist, sig=sig_trace, ecf=ecfbar))
}

mc_ng_ap2<- gibbs_ng_ap2()
cgrid<- sort(y)
plot(ecdf(y),xlab="X",ylab="Y",bty="l",xlim=c(-2,15),ylim=c(0,1), main="")
lines(cgrid,mc_ng_ap2$ecf,lwd=3,col=2)

cgrid<- sort(y)
plot(ecdf(y),xlab="X",ylab="Y",bty="l",xlim=c(-2,15),ylim=c(0,1), main="")
lines(cgrid,mc_py_ap1$ecf,lwd=3,col=2)
lines(cgrid,mc_py$ecf,lwd=3,col=3)
lines(cgrid, mc_ng_ap2$ecf,lwd=3,col=4)


K_d_ng <- max(abs( mc_ng_ap2$ecf - ecdf(y)(cgrid)))
K_d_ng







dev.off()


############################################################### original NG


