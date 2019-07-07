library(distrEx)

gen_data<- function(mu_vec=c(1,5,3,5), sigma_vec=rep(0.1, 4), ns=50, prob_vec=c(0.25,0.25,0.25,0.25) ){
  l<- length(mu_vec)
  components <- sample(1:l,prob=prob_vec,size=ns,replace=TRUE)
  mus <- mu_vec
  sds <- sqrt(sigma_vec)
  samples <- rnorm(n=ns,mean=mus[components],sd=sds[components])
  return(samples)
}

set.seed(123)
sample<- gen_data(mu_vec=c(1,5),sigma_vec=rep(0.2, 2),ns=200,prob_vec=c(0.5,0.5))
hist(sample,probability=TRUE,breaks=20,col=grey(.9))




y <- sample
n <- length(y)
## hyperparameters
a  <- 1;   b <- 1    # 1/sig ~ Ga(a,b)
m0 <- 0;  B0 <- 4   # G0 = N(m0,B0)
M  <- 1


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
    f0 <- dnorm(y[i], m=0,   s=sqrt(4+sig^2)) # q0
    pj <- c(fj*nj,f0*M)                       # p(s[i]=j | ...), j=1..k, k+1
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
  ##
  nj <- table(th)              # counts
  ths <- as.numeric(names(nj)) # unique values
  k <- length(nj)
  fx <- M/(n+M)*dnorm(xgrid,m=m0,sd=sqrt(B0+sig))
  for(j in 1:k)
    fx <- fx + nj[j]/(n+M)*dnorm(xgrid,m=ths[j],sd=sig)
  return(fx)
}
  
  
gibbs <- function(n.iter=5000){
    th <- init.DPk()              ## initialize th[1..n]
    sig <- sqrt( mean((y-th)^2))  ## and sig
    ## set up data structures to record imputed posterior draws..1
    xgrid <- seq(from=-2,to=10,length=50)
    fgrid <- NULL     ## we will record imputed draws of f
    njlist <- NULL    ## record sizes of 8 largest clusters
    klist <- NULL
    ## start with a plot of a kernel density estimate of the data
    plot(density(y),xlab="X",ylab="Y",bty="l",type="l",
         xlim=c(-2,10),ylim=c(0,0.7), main="")
    ## now the Gibbs sampler
    for(iter in 1:n.iter){
      th  <- sample.th(th,sig)
      sig <- sample.sig(th)
      th  <- sample.ths(th,sig)
      ## update running summaries #################
      if (iter>n.iter/2){
      f   <- fbar(xgrid,th,sig)
      lines(xgrid,f,col=iter,lty=3)
      fgrid <- rbind(fgrid,f)
      nj <- table(th)                 # counts
      njlist <- rbind(njlist,sort(nj,decr=T)[1:8])
      klist <- c(klist,length(nj))
      }
    }
    ## report summaries ###########################
    fbar <- apply(fgrid,2,mean)
    lines(xgrid,fbar,lwd=3,col=2)
    njbar <- apply(njlist,2,mean,na.rm=T)
    cat("Average cluster sizes:\n",format(njbar),"\n")
    pk <- table(klist)/length(klist)
    cat("Posterior probs p(k): (row1 = k, row2 = p(k) \n ")
    print(pk/sum(pk))
    return(list(fgrid=fgrid,klist=klist,njlist=njlist))
}
  


mcmc <- gibbs()


gibbs <- function(n.iter=500)
{
  #th <- init.DPk()              ## initialize th[1..n]
  th<- runif(length(y),-3,3)
  sig <- sqrt( mean((y-th)^2))  ## and sig
  ## set up data structures to record imputed posterior draws..1
  xgrid <- seq(from=0,to=6,length=50)
  fgrid <- NULL     ## we will record imputed draws of f
  njlist <- NULL    ## record sizes of 8 largest clusters
  klist <- NULL
  ## start with a plot of a kernel density estimate of the data
  plot(density(y),xlab="X",ylab="Y",bty="l",type="l",
       xlim=c(0,6),ylim=c(0,0.7), main="")
  ## now the Gibbs sampler
  for(iter in 1:n.iter){
    th  <- sample.th(th,sig)
    sig <- sample.sig(th)
    th  <- sample.ths(th,sig)
    ## update running summaries #################
    f   <- fbar(xgrid,th,sig)
    lines(xgrid,f,col=iter,lty=3)
    fgrid <- rbind(fgrid,f)
    nj <- table(th)                 # counts
    njlist <- rbind(njlist,sort(nj,decr=T)[1:8])
    klist <- c(klist,length(nj))
  }
  ## report summaries ###########################
  fbar <- apply(fgrid,2,mean)
  lines(xgrid,fbar,lwd=3,col=2)
  njbar <- apply(njlist,2,mean,na.rm=T)
  cat("Average cluster sizes:\n",format(njbar),"\n")
  pk <- table(klist)/length(klist)
  cat("Posterior probs p(k): (row1 = k, row2 = p(k) \n ")
  print(pk/sum(pk))
  return(list(fgrid=fgrid,klist=klist,njlist=njlist))
}


gibbs()


names(mcmc)
## report summaries
njbar <- apply(mcmc$njlist,2,mean,na.rm=T)
cat("Average cluster sizes:\n",format(njbar,digits=1),"\n")
pk <- table(mcmc$klist)/length(mcmc$klist)
cat("Posterior probs p(k): (row1 = k, row2 = p(k) \n ")
print(pk/sum(pk))

