

metrop_PY_discount <- function(theta, #previous iteration alpha.DP
                               pvec,
                               lik.fun,
                               ro.disc,
                               V=diag(theta), 
                               N, alpha.PY) { 
  accept<-FALSE
  ndim <- length(theta)
  last.lik<-lik.fun(theta,pvec=pvec,N=N,ro.disc,alpha.PY) 
  last=theta
  #sample from proposal p(a)=1/2dirac(0)+1/2U[0,0.5]
  c<-rbinom(n=1, size=1,prob=0.5)
  if(c==1){proposal<-0}else{proposal<-runif(1,min=0,max=0.5)}
  dproposal<-function(x){if(x==0) {return(0.5)}else{return(1)}}
  proposal.prior <-  log(dproposal(last)) #q(x)
  last.prior <-  log(dproposal(proposal)) #q(y)
  
  proposal.lik <- lik.fun(proposal,pvec,N,ro.disc,alpha.PY)
  alpha <- exp(proposal.lik+proposal.prior-last.lik-last.prior)
  if (alpha > runif(1) & !is.nan(alpha)) accept <- TRUE
  if (accept) {
    last <- proposal
  }
  return(last)
}



metrop_PY_alpha <- function(theta, #previous iteration alpha.DP
                            pvec,
                            lik.fun,
                            V=diag(theta), 
                            N, shape,rate,discount) {
  accept<-FALSE
  ndim <- length(theta)
  last.lik<-lik.fun(alpha=theta,pvec=pvec,N=N,shape=shape,rate=rate,discount=discount) 
  last=theta
  proposal <-.tnorm(1,lo=0,hi=Inf,mu=theta,sig=V)
  proposal.prior <-  log(dtruncnorm(last,a=0,b=Inf,mean=proposal,sd=V)) #q(x,y)
  last.prior <-  log(dtruncnorm(proposal,a=0,b=Inf,mean=last,sd=V)) #q(y,x)
  proposal.lik <- lik.fun(proposal,pvec,N,shape,rate,discount)
  alpha <- exp(proposal.lik+proposal.prior-last.lik-last.prior)
  if (alpha > runif(1) & !is.nan(alpha)) accept <- TRUE
  if (accept) {
    last <- proposal
  }
  return(last)
}




lik.alpha.fun<-function(alpha,pvec,N,shape,rate,discount){
  tmp<-sum(log(gamma(alpha+1+discount*(c(1:(N-1))-1)))-log(gamma((alpha+discount*c(1:(N-1))))))+alpha*log(pvec[N])+(shape-1)*log(alpha)-rate*alpha
  #tmp<-g_func(alpha,discount,N)*pvec[length(pvec)]^(alpha)*alpha^(shape-1)*exp(-rate*alpha)
  return(tmp)
}


lik.disc.fun<-function(discount,pvec,N,ro.disc,alpha){
  tmp<- (-N*log(gamma(1-discount)))+sum(log(gamma(alpha+1+discount*(c(1:(N-1))-1)))-log(gamma((alpha+discount*c(1:(N-1))))))-discount*sum(log(pvec[1:(N-1)]))+discount*(N-1)*pvec[N]+log(ro.disc*ifelse(discount==0,1,0)+2*(1-ro.disc)*ifelse((discount<=0.5 & discount>0),1,0))
  #tmp<-(1/(gamma(1-discount)^N))*g_func(alpha,discount,N)*prod(pvec[1:(N-1)]^(-discount))*(pvec[length(pvec)]^(discount*(N-1)))*(ro.disc*ifelse(discount==0,1,0)+2*(1-ro.disc)*ifelse((discount<=0.5 & discount>0),1,0))
  return(tmp)
}


