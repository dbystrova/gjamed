
library(plotly)
#Sys.setenv("plotly_username"="Daria_Bystrova")
#Sys.setenv("plotly_api_key"="IIYIoDs5TIMnLJQhlMrR")
#######################################################
## Algorithm 3

##G0 normal N(0,1)
###Fix parameters
theta<-1
alpha<-0.4
b<- theta/alpha
eps<-0.1


################################
B_function<- function(x,al){
  if (x==0) {
    res<- 1/((al^al)*((1-al)^(1-al)))
  } else {
    numerator<- sin(x)
    denumerator<- (sin(al*x)^al)*(sin((1-al)*x)^(1-al))
    res<- numerator/denumerator
  }
  return(res)
}

################################

Calculate_T<- function(b, alpha){
  sigma<- 1/sqrt(b*(1-alpha)*alpha)
  k<-0
  if (sigma >= sqrt(2*pi)) {
    success<-FALSE
    while (!success){
      U<-runif(1,0,pi)
      V<- runif(1,0,1)
      X<-U
      W<-B_function(X,alpha)
      success<-(V<=(W/B_function(0,alpha))^b)
      #print(0)
    }
  } else {
    success<-FALSE
    while(!success){
      N<- rnorm(1,0,1)
      V<-runif(1,0,1)
      X<-sigma*abs(N)
      W<-B_function(X,alpha)
      success<-((X<=pi)&(V*exp(-N^2/2) <= (W/B_function(0,alpha))^b)) 
      k<-k+1
      #print(k)
    }
  }
  return(W)
}
################################

Calculate_tau <- function(theta=1,alpha=0.5,eps=0.1){
  b<- theta/alpha
  W<- Calculate_T(b,alpha)
  G<- rgamma(1,1+b*(1-alpha))
  T<- 1/((W*(G^(1-alpha)))^(1/alpha))
  tau<- 1+ floor( (eps*T/alpha)^(-alpha/(1-alpha)))
  return (tau)
}




################################
###Samples from e-PY, with base measure G0 
PY_eps<-function(theta=1, alpha=0.4, e=0.1){
  tau<- Calculate_tau(theta,alpha,e)
  p<- vector(mode="numeric", length=tau)
  v<- vector(mode="numeric", length=tau)
  for(i in 1:tau){
    v[i] <- rbeta(1,1-alpha,theta+ i*alpha)
    if (i==1) {p[i]=v[1]
    }else {
    p[i] <- prod(1 - v[1:(i - 1)])*v[i]
    }
  }  
  R_tau<- 1-sum(p)
  #equal to prod(1-v)
  #sample tau+1 random variables from G0
  rvec<- rnorm(tau +1, 0,1)
  prob_vec<- c(R_tau,p)
  rList <- list("prob" = prob_vec, "atom_loc" = rvec)
  return(rList) 
}

xres<-PY_eps(1,0.5,0.1)
plot(xres$atom_loc,xres$prob, ylim=c(0,0.3),col=gray(.5),xlab=expression(italic(phi)),ylab=expression(italic(c)) ,lwd=3)



########Sample from tau(eps) for fixed(different) theta, alpha
n <- 10000
theta<-1
alpha<-0.4


a_list <- as.list(rep(alpha,n))
t_list<-  as.list(rep(theta,n))
eps_list<- as.list(rep(eps,n))


ltau<- mapply(Calculate_tau,t_list, a_list,eps_list)
a_list2 <- as.list(rep(0.5,n))
t_list2<-  as.list(rep(1,n))
eps_list2<- as.list(rep(0.1,n))
ltau2<- mapply(Calculate_tau,t_list2, a_list2,eps_list2)
a_list3 <- as.list(rep(0.6,n))
t_list3<-  as.list(rep(1,n))
eps_list3<- as.list(rep(0.1,n))
ltau3<- mapply(Calculate_tau,t_list3, a_list3,eps_list3)

t<- unlist(ltau, recursive = TRUE, use.names = TRUE)
t2<- unlist(ltau2, recursive = TRUE, use.names = TRUE)
t3<- unlist(ltau3, recursive = TRUE, use.names = TRUE)



############### Density for tau(epsilon) *manual
plot(density(t,bw = 2),type="l",col="blue",xlab=expression(tau),
     ylab="Density",xlim=c(0,150), ylim=c(0,0.05))
lines(density(t2,bw = 3),col="red",lty=2,lwd=2)
lines(density(t3),type="l",lty=3,lwd=2)
legend("topright",
       c("alpha=0.4", " alpha=0.5",
         "alpha=0.6"),
       col = c("blue","red","black"), lty =c(1,2,3))
abline(h=0,col="black")

#############Compute estimate on a grid

#####Grid for alpha[limited to 0.6], theta 
xx_m <- seq(0.1,0.6 , by = 0.05)
yy_m<- seq(1,11 , by = 1)


########compute 95% quantile for tau by rbase function *manual
compute_tau_est_q<- function(alpha,theta, eps=0.01){
  n<-1000
  a_list <- as.list(rep(alpha,n))
  t_list<-  as.list(rep(theta,n))
  eps_list<- as.list(rep(eps,n))
  ltau<- mapply(Calculate_tau,t_list, a_list,eps_list)
  q<-floor(quantile(ltau,0.95))
  return(q)
}

########compute 95%  for tau with sort function *manual
compute_tau_est_s<- function(alpha,theta, eps=0.01){
  n<-1000
  a_list <- as.list(rep(alpha,n))
  t_list<-  as.list(rep(theta,n))
  eps_list<- as.list(rep(eps,n))
  ltau<- mapply(Calculate_tau,t_list, a_list,eps_list)
  ltaus<-sort(ltau, decreasing=FALSE)
  q<-ltaus[round(length(ltaus)*0.95)]
  return(q)
}

z_q<- outer(xx_m,yy_m,Vectorize(compute_tau_est_q),eps=0.1)
z_s<-outer(xx_m,yy_m,Vectorize(compute_tau_est_s),eps=0.1)

z_s<-outer(xx_m,yy_m,Vectorize(compute_tau_est_s),eps=0.05)
z_s2<-outer(xx_m,yy_m,Vectorize(compute_tau_est_s),eps=0.03)
z_k<-outer(xx_m,yy_m,Vectorize(compute_tau_est_s),eps=0.01)


####compare qunatile(base) and  quantile manual
p_comp <- plot_ly() %>%
  add_surface(x=xx, y=yy,z = ~z_q, cmin = min(z_q), cmax = max(z_s),colorbar=list(title='Q R'), colorscale = list(c(0,1),c("rgb(255,112,184)","rgb(128,0,64)")),opacity = 0.98) %>%
  add_surface(x=xx, y=yy,z = ~z_s, cmin = min(z_q), cmax = max(z_s),colorbar=list(title='Q M'), colorscale = list(c(0,1),c("rgb(107,184,214)","rgb(0,90,124)"))) %>%
  layout(title="95% quantile for N and mean for N", scene = list(xaxis= list(title="alpha"),
                                                                 yaxis= list(title="theta"),
                                                                 zaxis= list(title="N",range = c(min(z_q),max(z_s)))))
  
p_comp
#COMPARE -OK

################# Plot for different errors
# 
# p <- plot_ly() %>%
#   add_surface(z = ~z_q, colorbar=list(title='Example 1')) %>%
#   add_surface(z = ~z_s, opacity = 0.80, colorbar=list(title='Example 2'))%>%
#   add_surface(z = ~z_s2, opacity = 0.80, colorbar=list(title='Example 2'))%>%
#   add_surface(z = ~z_k, opacity = 0.70, colorbar=list(title='Example 3'))
# 
# p
################# 3D plot on a grid

library("stabledist")


# Zolotaref function
zol_A = function(u,alpha){
  ((sin(alpha*u))^alpha*(sin((1-alpha)*u))^(1-alpha)/sin(u))^(1/(1-alpha))
}
zol_A = Vectorize(zol_A, vectorize.args = "u")

zol_B = function(u,alpha){
  sin(u)/((sin(alpha*u))^alpha*(sin((1-alpha)*u))^(1-alpha))
}
zol_B0 = function(alpha){
  1/(alpha^alpha*(1-alpha)^(1-alpha))
}

# random generation of random variable Z_{alpha,b}
# cf Section 3 of Devroye (2009)
rzol_inter = function(alpha,b){
  # sigma = 1/sqrt(b*alpha*(1-alpha))
  # if(sigma>sqrt(2*pi)){ # 2.506628
  V=1000
  X=.1
  while(V*zol_B0(alpha)^b>zol_B(X,alpha)^b){
    V = runif(1)
    X = runif(1,0,pi)
  }
  # }
  # something needs be corrected in that option: does not lead to correct rv
  # else{
  # V=1000
  # X=2*pi
  # while(X>pi | V*zol_B0(alpha)^b*exp(-X^2/2)>zol_B(X,alpha)^b){
  # V = runif(1)
  # X = abs(rnorm(1,0,sd = sigma))
  # }
  # }
  X
}

rzol = function(n,alpha,b){
  if(b==0){
    runif(n,0,pi)
  }
  else{
    Zvect = rep(0,n)
    for(i in 1:n)
      Zvect[i] = rzol_inter(alpha,b)
    Zvect
  }
}

# density function
dzol = function(x, alpha, b){
  if(b==0){
    dunif(x,0,pi)
  }
  else{
    gamma_constant = exp(lgamma(1+b*alpha)+lgamma(1+b*(1-alpha))-lgamma(1+b))
    gamma_constant/pi/zol_A(x,alpha)^(b*(1-alpha))
  }
}
dzol = Vectorize(dzol, vectorize.args = "x")


# random generation of polynomially tilted unilateral stable random variable T_{alpha,beta}
# cf Section 2 and Theorem 1 of Devroye (2009)
rtstable = function(n, alpha, beta){
  if(beta==0){
    Z = runif(n,0,pi)
    G = rexp(n,1)
  }
  else{
    Z = rzol(n, alpha, beta/alpha)
    G = rgamma(n,1+beta*(1-alpha)/alpha)
  }
  (zol_A(Z, alpha)/G)^((1-alpha)/alpha)
}

n<-10000
compute_tau_vec<- function(alpha,theta){
  #alpha<-0.6
  #theta<-1
  z<-rtstable(n,alpha,theta)
  #plot(density(z))
  T<-rtstable(n,alpha,theta)
  eps<- 0.1
  tau<-1+ floor(((eps*T)/alpha)^(-alpha/(1-alpha)))
  return(tau)
}

tau1<-compute_tau_vec(0.4,1)
tau2<-compute_tau_vec(0.5,1)
tau3<-compute_tau_vec(0.6,1)
plot(density(tau3),xlim=c(0,150), ylim=c(0,0.05))
lines(density(tau2),col="blue")
lines(density(tau1),col="red")
#abline(v=q[2])



#par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
#par(mfrow=c(1,2))
#pdf("density.pdf")
plot(density(tau1,bw = 2),type="l",col="blue",xlab=expression(tau),
     ylab="Probability",xlim=c(0,150), ylim=c(0,0.05), main=expression(paste("Density plot for ",italic(tau),"(",italic(epsilon),")")))
lines(density(tau2,bw = 3),col="red",lty=2,lwd=2)
lines(density(tau3),type="l",lty=3,lwd=2)
legend("topright",
       c("alpha=0.4", "alpha=0.5",
         "alpha=0.6"),
       col = c("blue","red","black"), lty =c(1,2,3))
#abline(h=0,col="black")
#dev.off()
###################Generating 3D plots.

#####Grid for alpha[limited to 0.6], theta 
xx <- seq(0.1,0.6 , by = 0.05)
yy<- seq(1,11 , by = 1)


########compute 95% quantile for tau by rbase function
compute_tau_est<- function(alpha,theta, eps=0.01){
  T<-rtstable(n,alpha,theta)
 # eps<- 0.1
  tau<-1+ floor(((eps*T)/alpha)^(-alpha/(1-alpha)))
  q<-floor(quantile(tau,0.95))
  return(q)
}

########compute 95%  for tau with sort function 
compute_tau_est_b<- function(alpha,theta, eps=0.01){
  T<-rtstable(n,alpha,theta)
  tau<-1+ floor(((eps*T)/alpha)^(-alpha/(1-alpha)))
  taus<-taus<-sort(tau, decreasing=FALSE)
  q<-taus[round(length(taus)*0.95)]
  return(floor(q))
}

########compute mean  for tau
compute_tau_mean<- function(alpha,theta, eps=0.01){
  N_eps<- (eps/alpha)^(-alpha/(1-alpha))
  gamma<- (gamma(1+theta)*gamma(1+ theta/alpha + 1/(1-alpha)))/(gamma(1+theta/alpha)*gamma(1+ alpha/(1-alpha) + theta))
  N<- N_eps*gamma
  return(N)
}


########compute k-th moment  for tau
compute_tau_km<- function(alpha,theta, eps=0.1,k=1){
  N_eps<- (eps/alpha)^(-alpha/(1-alpha))
  gamma<- (gamma(1+theta)*gamma(1+ theta/alpha + k/(1-alpha)))/(gamma(1+theta/alpha)*gamma(1+ (k*alpha)/(1-alpha) + theta))
  N<- N_eps*gamma
  return(N)
}

############## constant gamma function in density
C_const<-function(alpha,theta){
  gamma_c<-(gamma(1+theta))/(gamma(1+theta/alpha))
  return(gamma_c)
}
##############  gamma function depending on kth moment
M<-function(alpha,theta,m=1){
  gamma_m<-(gamma(1+ theta/alpha + m/(1-alpha)))/(gamma(1+ (m*alpha)/(1-alpha) + theta))
  return(gamma_m)
}

##############  compute variance
compute_tau_var<- function(alpha,theta, eps=0.1){
  N_eps<- (eps/alpha)^(-alpha/(1-alpha))
  gamma2<- (gamma(1+theta)*gamma(1+ theta/alpha + 2/(1-alpha)))/(gamma(1+theta/alpha)*gamma(1+ (2*alpha)/(1-alpha) + theta))
  gamma1<- (gamma(1+theta)*gamma(1+ theta/alpha + 1/(1-alpha)))/(gamma(1+theta/alpha)*gamma(1+ (1*alpha)/(1-alpha) + theta))
  gamma<-(gamma2-gamma1*gamma1)
  N<- (N_eps^2)*gamma
  return(sqrt(N))
}

##############  compute skewness
compute_tau_sk<- function(alpha,theta, eps=0.01){
  gamma2<-M(alpha,theta,2)
  gamma1<-M(alpha,theta,1)
  gamma<-(C_const(alpha,theta))*M(alpha,theta,3) - 3*(C_const(alpha,theta)^2)*(gamma1*(gamma2-gamma1))- ((C_const(alpha,theta))^3)*(gamma1)^3
  gamma_fin<-gamma/((C_const(alpha,theta)*(gamma2-gamma1))^3/2)
  return(gamma_fin)
}


############################################################


compute_assymp<- function(alpha,theta, eps=0.1,k=1){
  val<- exp(1*k*log(10)*(alpha/(1-alpha)))*((theta/(1-alpha))^k)
  return(val)
}


##Vectorization
####compute quantile by R
z<- outer(xx,yy,Vectorize(compute_tau_est),eps=0.1)
######compute quantile by sort 
z01<-outer(xx,yy,Vectorize(compute_tau_est_b),eps=0.1)
z005<-outer(xx,yy,Vectorize(compute_tau_est_b),eps=0.05)
z001<-outer(xx,yy,Vectorize(compute_tau_est_b),eps=0.01)


p_comp_qq <- plot_ly(showscale = TRUE) %>%
  add_surface(x=xx, y=yy,z = ~z, cmin = min(z), cmax = max(z01),colorbar=list(title='Q R'), colorscale = list(c(0,1),c("rgb(255,112,184)","rgb(128,0,64)"))) %>%
  add_surface(x=xx, y=yy,z = ~z01, cmin = min(z), cmax = max(z01),colorbar=list(title='Q S'),  colorscale = list(c(0,1),c("rgb(107,184,214)","rgb(0,90,124)"))) %>%
  add_surface(x=xx, y=yy,z = ~z_q, cmin = min(z), cmax = max(z01), colorbar=list(title='Q MM'), colorscale = list(c(0,1),c("rgb(100,154,214)","rgb(0,90,124)"))) %>%
  layout(title="95% quantile for number of iterations by dif methods", scene = list(xaxis= list(title="alpha"),
                                                                                            yaxis= list(title="theta"),
                                                                                            zaxis= list(title="N",range = c(min(z),max(z01)))))
p_comp_qq

chart_link = api_create(p_comp_qq, filename="CHECK")
chart_link
##########Verification that they coincide -DONE
################################################################################################################################################
#######Surface for different errors
p_comp_qe <- plot_ly(showscale = TRUE) %>%
  add_surface(x=xx, y=yy,z = ~z01, cmin = min(z01), cmax = max(z01),colorbar=list(title='e=0.1'), colorscale = list(c(0,1),c("rgb(255,112,184)","rgb(128,0,64)")),opacity = 0.98) %>%
  add_surface(x=xx, y=yy,z = ~z005, cmin = min(z005), cmax = max(z005),colorbar=list(title='e=0.1'),  colorscale = list(c(0,1),c("rgb(107,184,214)","rgb(0,90,124)")),opacity = 0.96) %>%
  add_surface(x=xx, y=yy,z = ~z001, cmin = min(z001), cmax = max(z001),colorbar=list(title='e=0.01'),  colorscale = list(c(0,1),c("rgb(100,154,60)","rgb(150,50,0)")),opacity = 0.78) %>%
  layout(title="95% quantile for number of iterations for different epsilon ", scene = list(xaxis= list(title="alpha"),
                                                                                            yaxis= list(title="theta"),
                                                                                            zaxis= list(title="N",range = c(min(z01),max(z005)))))
p_comp_qe
chart_link = api_create(p_comp_qe, filename="Different error")
chart_link
plotly_IMAGE(p_comp_qe, format = "png", out_file = "Errors.png")

################################################################################################################################################

z_mean<- outer(xx,yy,Vectorize(compute_tau_mean),eps=0.1)


#######Surface for mean and quantile
p_comp_mq <- plot_ly(showscale = TRUE) %>%
  add_surface(x=xx, y=yy,z = ~z01, cmin = min(z01), cmax = max(z01),colorbar=list(title='Quantile'), colorscale = list(c(0,1),c("rgb(255,112,184)","rgb(128,0,64)"))) %>%
  add_surface(x=xx, y=yy,z = ~z_mean, cmin = min(z_mean), cmax = max(z_mean),colorbar=list(title='Mean'),  colorscale = list(c(0,1),c("rgb(107,184,214)","rgb(0,90,124)"))) %>%
   layout(title="Plot for 95% quantile and mean for e=0.1", scene = list(xaxis= list(title="alpha"),
                                                                                            yaxis= list(title="theta"),
                                                                                            zaxis= list(title="N",range = c(min(z_mean),max(z01)))))
p_comp_mq
chart_link = api_create(p_comp_mq, filename="Mean and Quantile")
chart_link
plotly_IMAGE(p_comp_mq, format = "png", out_file = "MeanQ.png")


################################################################################################################################################

z4<-outer(xx,yy,Vectorize(compute_tau_var),eps=0.1)

p_var <- plot_ly(showscale = TRUE) %>%
  add_surface(x=xx, y=yy,z = ~z4, cmin = min(z4), cmax = max(z4),colorbar=list(title='Quantile'), colorscale = list(c(0,1),c("rgb(100,112,184)","rgb(128,50,64)"))) %>%
  layout(title="Plot for Variance, e=0.1", scene = list(xaxis= list(title="alpha"),
                                                                        yaxis= list(title="theta"),
                                                                        zaxis= list(title="Variance",range = c(min(z4),max(z4)))))


p_var
chart_link = api_create(p_var, filename="Variance")
chart_link
plotly_IMAGE(p_var, format = "png", out_file = "Variance.png")



################################################################################################################################################



#######Surface for mean +2 sd and quantile for error 0.05

z_var005<-outer(xx,yy,Vectorize(compute_tau_var),eps=0.05)


z_mean005<- outer(xx,yy,Vectorize(compute_tau_mean),eps=0.05)

z_2sd_05<-z_mean005+2*z_var005


p_comp_2sd05 <- plot_ly(showscale = TRUE) %>%
  add_surface(x=xx, y=yy,z = ~z005, cmin = min(z005), cmax = max(z005),colorbar=list(title='Quantile'), colorscale = list(c(0,1),c("rgb(255,112,184)","rgb(128,0,64)"))) %>%
  add_surface(x=xx, y=yy,z = ~z_2sd_05, cmin = min(z_2sd_05), cmax = max(z_2sd_05),colorbar=list(title='Mean+2sd'),  colorscale = list(c(0,1),c("rgb(107,184,214)","rgb(0,90,124)"))) %>%
  layout(title="Plot for 95% quantile and mean for e=0.05", scene = list(xaxis= list(title="alpha"),
                                                                        yaxis= list(title="theta"),
                                                                        zaxis= list(title="N",range = c(min(z005),max(z005)))))

p_comp_2sd05
chart_link = api_create(p_comp_2sd05 , filename="Mean+2sd and Quantile 05")
chart_link
plotly_IMAGE(p_comp_2sd05 , format = "png", out_file = "Mean2sdQ05.png")


################################################################################################################################################

z_2sd<-z_mean + 2*z4

#######Surface for mean+2 sd   and quantile
p_comp_2sd_q<- plot_ly(showscale = TRUE) %>%
  add_surface(x=xx, y=yy,z = ~z01, cmin = min(z01), cmax = max(z01),colorbar=list(title='Quantile'), colorscale = list(c(0,1),c("rgb(255,112,184)","rgb(128,0,64)"))) %>%
  add_surface(x=xx, y=yy,z = ~z_2sd, cmin = min(z_2sd), cmax = max(z_2sd),colorbar=list(title='Mean'),  colorscale = list(c(0,1),c("rgb(107,184,214)","rgb(0,90,124)"))) %>%
  layout(title="Plot for 95% quantile and mean + 2 sd for e=0.1", scene = list(xaxis= list(title="alpha"),
                                                                        yaxis= list(title="theta"),
                                                                        zaxis= list(title="N",range = c(min(z_2sd),max(z_2sd)))))



p_comp_2sd_q



chart_link = api_create(p_comp_2sd_q, filename="Mean+2sd and Quantile")
chart_link
plotly_IMAGE(p_comp_2sd_q, format = "png", out_file = "Mean2sdQ.png")






################################################################################################################################################


z_km<- outer(xx,yy,Vectorize(compute_tau_km),eps=0.1,k=2)


################################################################################################################################################
z_ass<- outer(xx,yy,Vectorize(compute_assymp),eps=0.1,k=2)





################################################################################################################################################


calc_f<-function(alpha,theta,k=1,eps=0.05){
  val<- alpha^(-k)*(theta + alpha/(1-alpha))^k
  c<- (eps)^(-alpha/(1-alpha))
  val<-val*c
  return (val)
}

zk1<-outer(xx,yy,Vectorize(calc_f),k=1)


p <- plot_ly(x=xx, y=yy, z=zk1) %>% add_surface() %>% 
  layout(scene= list(xaxis= list(title="alpha"),
                     yaxis= list(title="theta"),
                     zaxis= list(title="N")))
p



pk <- plot_ly(showscale = TRUE) %>%
  add_surface(x=xx, y=yy,z = ~z1, cmin = min(zk1), cmax = max(z1), colorscale = list(c(0,1),c("rgb(255,112,184)","rgb(128,0,64)"))) %>%
  add_surface(x=xx, y=yy,z = ~zk1, cmin = min(zk1), cmax = max(z1), colorscale = list(c(0,1),c("rgb(107,184,214)","rgb(0,90,124)"))) %>%
   layout(title="95% quantile for number of iterations for different epsilon ", scene = list(xaxis= list(title="alpha"),
                                                                                            yaxis= list(title="theta"),
                                                                                            zaxis= list(title="N")))
pk


############################################################################################################################################
f
gamma_fun <- function(alpha) {
  (cos(pi*alpha/2))^(1/alpha)
}



lissage  =function(m){
  nr = dim(m)[1]
  nc = dim(m)[2]
  for (i in 1:(nr)){
    for (j in 2:(nc-1)){
      v = m[i,(j-1):(j+1)]
      m[i,j] = mean(v)
    }
  }  
  for (j in 1:(nc)){
    for (i in 2:(nr-1)){
      v = m[(i-1):(i+1),j]
      m[i,j] = mean(v)
    }
  }
  m
}

# applying smoothing to z matrix of integers

z1 = lissage(z1) # do it about ten times :-)

# 3D-plot
build3ds1<-function(x,y,z,par1="",par2=""){
  z<-pmin(z,1520)
  nrz <- nrow(z)
  ncz <- ncol(z)
  jet.colors <- colorRampPalette( c("blue", "red") )
  nbcol <- 100
  color <- jet.colors(nbcol)
  zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
  facetcol <- cut(zfacet, nbcol)
  
  
  sigma_param<-x
  theta_param<-y
  #pdf()
  par(mar = c(2,1,2,.1), mgp = c(1,1,0), xaxs = "i", yaxs = "i", las= 1)
  persp(
    x = sigma_param, 
    y = theta_param, 
    z = z, 
    zlim = c(0,1520), 
    col = color[facetcol], 
    theta = -25, 
    phi = 25,
    ltheta = 120, 
    ticktype = "detailed", 
    shade = 0.3,
    xlab = "", ylab = "", zlab = "", 
    d = 5, r = 10,
    cex.axis = 1, cex.lab = 1.3, nticks = 3, main = expression(paste(italic(N^1),"(",italic(epsilon),",",italic(alpha),",",italic(theta),") ,",
                                                                     paste(epsilon),"=",0.1))
  )
  text(.13,-.37,expression(alpha), cex = 1.5)
  text(-.3,-.35,expression(theta), cex = 1.5)
  #dev.off()
}

#nrz <- nrow(z1)
# ncz <- ncol(z1)
# jet.colors <- colorRampPalette( c("blue", "red") )
# nbcol <- 100
# color <- jet.colors(nbcol)
# zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
# facetcol <- cut(zfacet, nbcol)
# 
# 
# sigma_param<-xx
# theta_param<-yy
# #pdf()
# par(mar = c(2,1,2,.1), mgp = c(1,1,0), xaxs = "i", yaxs = "i", las= 1)
# persp(
#   x = sigma_param, 
#   y = theta_param, 
#   z = z, 
#   zlim = c(0,1520), 
#   col = color[facetcol], 
#   theta = -25, 
#   phi = 25,
#   ltheta = 120, 
#   ticktype = "detailed", 
#   shade = 0.3,
#   xlab = "", ylab = "", zlab = "", 
#   d = 5, r = 10,
#   cex.axis = 1, cex.lab = 1.3, nticks = 3, main = expression(italic(tau)*(italic(epsilon)))
# )
# text(.15,-.37,expression(alpha), cex = 1.5)
# text(-.3,-.35,expression(theta), cex = 1.5)
#dev.off()
####Plot errors##################################################################
#z01<-lissage(z01)
#z005<-lissage(z005)
#z001<-lissage(z001)

pdf("plot_e1.pdf")
build3ds1(xx,yy,z01,par1=0.1)
build3ds1(xx,yy,z005,par1=0.05)
build3ds1(xx,yy,z001,par1=0.01)
dev.off()


build3ds<-function(x,y,z,par1="",par2=""){
  z<-pmin(z,1520)
  nrz <- nrow(z)
  ncz <- ncol(z)
  jet.colors <- colorRampPalette( c("blue", "red") )
  nbcol <- 100
  color <- jet.colors(nbcol)
  zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
  facetcol <- cut(zfacet, nbcol)
  
  
  sigma_param<-x
  theta_param<-y
  #pdf()
  par(mar = c(2,1,2,.1), mgp = c(1,1,0), xaxs = "i", yaxs = "i", las= 1)
  persp(
    x = sigma_param, 
    y = theta_param, 
    z = z, 
    zlim = c(0,320), 
    col = color[facetcol], 
    theta = -25, 
    phi = 25,
    ltheta = 120, 
    ticktype = "detailed", 
    shade = 0.3,
    xlab = "", ylab = "", zlab = "", 
    d = 5, r = 10,
    cex.axis = 1, cex.lab = 1.3, nticks = 3, main = expression(paste("Variance for ",
                                                                     paste(epsilon),"=",0.1))
  )
  text(.13,-.37,expression(alpha), cex = 1.5)
  text(-.3,-.35,expression(theta), cex = 1.5)
  #dev.off()
}






pdf("plot_variance.pdf")
build3ds(xx,yy,z4,par1=0.1)
dev.off()




build3d_diff<-function(x,y,z,par1="",lim=max(z)){
  z<-pmin(z,1520)
  nrz <- nrow(z)
  ncz <- ncol(z)
  jet.colors <- colorRampPalette( c("blue", "red") )
  nbcol <- 100
  color <- jet.colors(nbcol)
  zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
  facetcol <- cut(zfacet, nbcol)
  
  
  sigma_param<-x
  theta_param<-y
  par(mar = c(2,1,2,.1), mgp = c(1,1,0), xaxs = "i", yaxs = "i", las= 1)
  persp(
    x = sigma_param, 
    y = theta_param, 
    z = z, 
    zlim = c(min(z),lim), 
    col = color[facetcol], 
    theta = -25, 
    phi = 25,
    ltheta = 120, 
    ticktype = "detailed", 
    shade = 0.3,
    xlab = "", ylab = "", zlab = "", 
    d = 5, r = 10,
    cex.axis = 1, cex.lab = 1.3, nticks = 3, main = par1)
  
  text(.13,-.37,expression(alpha), cex = 1.5)
  text(-.3,-.35,expression(theta), cex = 1.5)
}



z_2sd<-z_mean + 2*z4
z_diff<-z_2sd-z01 
pdf("difference_msd_q.pdf")
build3d_diff(xx,yy,z_diff,par1=expression(paste("Difference for ",paste(mu),"+2*",paste(sigma)," and 95% quantile for ",
                                                paste(epsilon),"=",0.1)))
dev.off()

#################################################################################"


##########################################DIFFERENCE FOR e=0.05########################################################################################


z005<-outer(xx,yy,Vectorize(compute_tau_est_b),eps=0.05)

z_mean005<- outer(xx,yy,Vectorize(compute_tau_mean),eps=0.05)

z_var005<-outer(xx,yy,Vectorize(compute_tau_var),eps=0.05)


z_2sd_05<-z_mean005+2*z_var005

z_diff05<-z_2sd_05-z005 


pdf("difference_msd_q05.pdf")
build3d_diff(xx,yy,z_diff05,par1=expression(paste("Difference for ",paste(mu),"+2*",paste(sigma)," and 95% quantile for ",
                                                  paste(epsilon),"=",0.05)))
dev.off()


#################################################################################"


z_trunc01<-pmin(z_2sd,1000)
z_trunc05<-pmin(z_2sd_05,1000)





pdf("plot_trunc01.pdf")
build3d_diff(xx,yy,z_trunc01,par1=expression(paste("Truncation number for ",
                                                    paste(epsilon),"=",0.1)),lim=1500)
dev.off()

pdf("plot_trunc05.pdf")
build3d_diff(xx,yy,z_trunc05,par1=expression(paste("Truncation number for ",
                                                paste(epsilon),"=",0.05)),lim=1500)
dev.off()











paste("Jaccard vs. Tsim for depths",  min.depth, "to",max.depth,"m", sep=" "))













####SOME DRAFT TESTS##################################################################
# alpha<-0.5
# theta<-1
# x_s<-seq(1/10000,0.1,by=0.001)
# y1<-1/(x_s^(alpha/(1-alpha)))
# y2<-theta*log(1/x_s)
# plot(x_s,y1,col=gray(0),xlab=expression(italic(epsilon)),ylab=expression(italic(tau(epsilon))),xlim=c(0.1,0.0001),lwd=3, type="l")
# text(1.1,5.1,expression(sigma), cex = 1.5)
# lines(x_s,y2,col="red")
# 

# 
# eps<- 0.1
# alpha<-0.5
# theta<-1
# 
# N_eps<- (eps/alpha)^(-alpha/(1-alpha))
# gamma<- (gamma(1+theta)*gamma(1+ theta/alpha + 1/(1-alpha)))/(gamma(1+theta/alpha)*gamma(1+ alpha/(1-alpha) + theta))
# N<- N_eps*gamma
# 






