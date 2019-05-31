
metrop_DP <- function(theta, #previous iteration alpha.DP
                      pvec,
                      lik.fun,
                      #prior.fun,
                      V=diag(theta), 
                      N, shape,rate #custom
) {
  accept<-FALSE
  ndim <- length(theta)
  last.lik<-lik.fun(alpha=theta,pvec=pvec,N=N,shape=shape,rate=rate)
  last=theta
  proposal <-.tnorm(1,lo=0,hi=Inf,mu=theta,sig=V)
  proposal.prior <-  log(dtruncnorm(last,a=0,b=Inf,mean=proposal,sd=V)) #q(x,y)
  last.prior <-  log(dtruncnorm(proposal,a=0,b=Inf,mean=last,sd=V)) #q(y,x)
  
  proposal.lik <- lik.fun(proposal,pvec,N,shape,rate)
  alpha <- exp(proposal.lik+proposal.prior-last.lik-last.prior)
  if (alpha > runif(1) & !is.nan(alpha) ) accept <- TRUE
  if (accept) {
    last <- proposal
  }
  return(last)
}

