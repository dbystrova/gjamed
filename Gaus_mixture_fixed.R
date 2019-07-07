mu<-c(-3,0,3)
s2<-c(.33,.33,.33)
w<-c(.45,.1,.45)

ths<-seq(-5,5,length=100)
plot(ths, w[1]*dnorm(ths,mu[1],sqrt(s2[1])) +
       w[2]*dnorm(ths,mu[2],sqrt(s2[2])) +
       w[3]*dnorm(ths,mu[3],sqrt(s2[3])) ,type="l" )



#### MC Sampling
set.seed(1)
S<-2000
d<-sample(1:3,S, prob=w,replace=TRUE)
th<-rnorm(S,mu[d],sqrt(s2[d]))
THD.MC<-cbind(th,d)
####

#### MCMC sampling
th<-0
THD.MCMC<-NULL
S<-10000
set.seed(1)
for(s in 1:S) {
  d<-sample(1:3 ,1,prob= w*dnorm(th,mu,sqrt(s2))   )
  th<-rnorm(1,mu[d],sqrt(s2[d]) )
  THD.MCMC<-rbind(THD.MCMC,c(th,d) )
}

plot(THD.MCMC[,1])
lines( mu[THD.MCMC[,2]])

###
#pdf("fig6_4.pdf",family="Times",height=3.5,width=7)
par(mfrow=c(1,1),mar=c(3,3,1,1),mgp=c(1.75,.75,0))
ths<-seq(-6,6,length=1000)
plot(ths, w[1]*dnorm(ths,mu[1],sqrt(s2[1])) +
       w[2]*dnorm(ths,mu[2],sqrt(s2[2])) +
       w[3]*dnorm(ths,mu[3],sqrt(s2[3])) ,type="l" , xlab=expression(theta),ylab=
       expression( paste( italic("p("),theta,")",sep="") ),lwd=2 ,ylim=c(0,.40))
hist(THD.MC[,1],add=TRUE,prob=TRUE,nclass=20,col="gray")
lines( ths, w[1]*dnorm(ths,mu[1],sqrt(s2[1])) +
         w[2]*dnorm(ths,mu[2],sqrt(s2[2])) +
         w[3]*dnorm(ths,mu[3],sqrt(s2[3])),lwd=2 )
#dev.off()