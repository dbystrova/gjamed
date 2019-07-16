
rm(list=ls())
#Packages for gamma function estimation
library(pracma)
library(Brobdingnag)
library(copula)
library(plotly)
library(distrEx)
library(expint)

############### V for PY / Dir / NGG########################
v_py<- function(kval, sigma,theta,npoints){
  c_v<-1:(kval-1)
  v_nk<- (theta +sigma*c_v)
  Vup<- prod(v_nk)
  n_vec<- 1:(npoints-1)
  Vlow<- prod(theta +n_vec)
  V_t_nk<-Vup/Vlow
  return(V_t_nk)
}

############### V for NGG########################
v_ng<- function(beta, sigma, kval, npoints){
  sum<-0
  coef_low<-as.brob(gamma(npoints))
  coef_high<-as.brob(exp(beta)* sigma^(kval-1))
  coef<- coef_high/coef_low
  incv<- as.brob(vector(length=npoints))
  for(i in (0:(npoints-1))){
    gn<- as.brob(gammainc(kval - i/sigma,beta))
    #gn<- gamma_inc(kval - i/sigma,beta)
    ckn<- as.brob(choose(npoints-1,i))
    sum<- sum + ((-1)^i)*(beta^(i/sigma))*ckn*gn
    incv[i+1]<- gn
  }
  sumf<- sum*coef
  sumn<- as.numeric(sumf)
  return(list(sum=sumn, incg= incv)) 
}


v_ng2<- function(beta, sigma, kval, npoints){
  sum<-0
  incv<- as.brob(vector(length=npoints))
  for(i in (0:(npoints-1))){
    coef_low<-as.brob(gamma(npoints))
    coef_high<-as.brob(exp(beta)* sigma^(kval-1))
    coef<- coef_high/coef_low
    gn<- as.brob(gammainc(kval - i/sigma,beta))
    ckn<- as.brob(choose(npoints-1,i))
    sum<- sum + ((-1)^i)*(beta^(i/sigma))*ckn*gn*coef
    incv[i+1]<- gn
  }
  sumf<- sum
  sumn<- as.numeric(sumf)
  return(list(sum=sumn, incg= incv)) 
}

###############Generalized coefficient########################

gen_fac_coef<-function(kval,sigma,npoints){
  sum<-0
  kfac<-factorial(kval)
  for(i in (0:kval)){
    n_vec<- 0:(npoints-1)
    sn<- prod(-i*sigma +n_vec)
    ckn<- choose(kval,i)
    #print((-1)^i*ckn*sn)
    sum<- sum + ((-1)^i)*ckn*sn 
  }
  sumf<- sum/kfac
  return(sumf)
}


########### density for  NGG #############################
prob_ng<- function(kg, sigma, npoints, beta){
  pb_v_all<- v_ng(beta, sigma, kg, npoints)
  pb_v<-as.brob(pb_v_all$sum)
  pb_gen<- as.brob(gen_fac_coef(kg,sigma, npoints))
  prob<- (pb_v*pb_gen)/(as.brob(sigma^kg))
  prob_num<- as.numeric(prob)
  return(prob_num)
}

########### density for PY#############################

prob_py<- function(kg, npoints, sigma, theta){
  pb_v<- v_py(kg,sigma, theta,npoints)
  pb_gen<- gen_fac_coef(kg,sigma, npoints)
  prob<- (pb_v*pb_gen)/(sigma^kg)
  return(prob)
}

########### density for Dirichlet#############################
prob_dir<- function(k, npoints, theta){
  n_vec<- 0:(npoints-1)
  theta_n<- prod(theta +n_vec)
  prob<- ((theta^k) *(abs(Stirling1(npoints,k))))/theta_n
  return(prob)
}

prob_dir_large_dim<- function(k, npoints, theta){
  n_vec<-as.brob( 0:(npoints-1))
  theta_n<- prod(theta +n_vec)
  stir<- as.brob(abs(Stirling1(npoints,k)))
  powerk<- as.brob((theta^k))
  prob_brob<- powerk*(stir/theta_n)
  prob<- as.numeric(prob_brob)
  return(prob)
}

#####################################################################################################################################
#prib_ng(kg, npoints, sigma, beta)
k_vec<-seq(1,50,by=5)
sigma_vec<-seq(0.2,0.8, by=0.1)
z<- outer(k_vec,sigma_vec,Vectorize(prob_ng),npoints=50, beta=1)

#p<- plot_ly(showscale = TRUE) %>%
#  add_surface(x=k_vec, y=sigma_vec,z =z, cmin = min(z), cmax = max(z),colorbar=list(title='PY'), colorscale = list(c(0,1),c("rgb(255,112,184)","rgb(128,0,64)")),opacity = 0.98) %>%
#  layout(title="Prior distribution", scene = list(xaxis= list(title="K"),yaxis= list(title="sigma"),zaxis= list(title="N",range = c(min(z),max(z)))))
#p


sigma_df <- rep(seq(0.2,0.8, by=0.6/9), each=50)
k_df<- rep(1:50, 10)
df<- as.data.frame(matrix(NA, nrow=500, ncol=1))

df$k<- k_df
df$sigma<- sigma_df
for(l in (1: nrow(df))){
  df$val[l]<- prob_ng(df$k[l],sigma=df$sigma[l],npoints=50,beta=5)
}

df$sigma<- as.factor(df$sigma)
p<- plot_ly(df, x =df$k , y = df$sigma, z = df$val, split = df$sigma, type = "scatter3d", mode = "lines") %>% 
  layout(title="Prior density", scene = list(xaxis= list(title="K"),yaxis= list(title="sigma"), zaxis= list(title="Probability")))                                                                                          
p



##################################################Approximation for the PY##############
##PY and PY 2 coincide 

#For PY Vnk/Vn = sigma*kn/n


