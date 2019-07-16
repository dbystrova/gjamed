set.seed(123)
sample<- gen_data(mu_vec=c(1,10),sigma_vec=rep(0.2, 2),ns=500,prob_vec=c(0.5,0.5))
save(sample, file = "sample_data_n500_1.Rds")

#hist(sample,probability=TRUE,breaks=20,col=grey(.9))
#plot(density(sample,bw = "sj"))
y_data <- sample
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
n_grid<-500
n_sm<-500
niter<-1000 
burnin<- 200
for(m in 1:length(M_vec)){
  for(j in 1:length(sigma_vec)){
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




k_mat_5001k<- as.data.frame(K_dmat)
k_mat_5001k$n<- length(y_data)
names(k_mat_5001k)<- c("alpha","sigma","PY","PY_A1","NG_A2","N")
write.csv(k_mat_5001k, file = "K_500_1K.csv")
formattable(k_mat_5001k)

F_mat_5001k_PY<- as.data.frame(F_dmat_PY)
F_mat_5001k_PY1<- as.data.frame(F_dmat_PY1)
F_mat_5001k_NG2<- as.data.frame(F_dmat_NG2)
save(F_mat_5001k_PY, file = "FmatPY_n5001k.Rds")
save(F_mat_5001k_PY1, file = "FmatPY1_n5001k.Rds")
save(F_mat_5001k_NG2, file = "FmatNG2_n5001k.Rds")

E_mat_5001k_PY<- as.data.frame(E_dmat_PY)
E_mat_5001k_PY1<- as.data.frame(E_dmat_PY1)
E_mat_5001k_NG2<- as.data.frame(E_dmat_NG2)
save(E_mat_5001k_PY, file = "EmatPY_n5001k.Rds")
save(E_mat_5001k_PY1, file = "EmatPY1_n5001k.Rds")
save(E_mat_5001k_NG2, file = "EmatNG2_n5001k.Rds")

x
xgrid<-  seq(from=-2,to=15,length=length(y_data)) 
cgrid<-  seq(from=-2,to=15,length=500)
true_ecdf<- 0.5*pnorm(cgrid, 1,sqrt(0.2)) +0.5*pnorm(cgrid, 10,sqrt(0.2)) 

plot(cgrid, true_ecdf,xlab="X",ylab="Y",bty="l",type="l",ylim=c(0,1), main="")
lines(cgrid,E_mat_500_PY[1,1:500],lwd=3,col="red")
lines(cgrid,E_mat_500_PY1[1,1:500],lwd=3,col="blue")
lines(cgrid,E_mat_500_NG2[1,1:500],lwd=3,col="green")


true_dens<- 0.5*dnorm(xgrid, mean=1, sd= sqrt(0.2)) +0.5*dnorm(xgrid, mean=10,sd=sqrt(0.2)) 
plot(xgrid,true_dens,xlab="X",ylab="Y",bty="l",type="l",xlim=c(-2,15),ylim=c(0,1), main="")
for(iter in burnin:niter){ lines(xgrid, F_mat_500_PY[iter,1:200],col=iter,lty=3) }

true_dens<- 0.5*dnorm(xgrid, mean=1, sd= sqrt(0.2)) +0.5*dnorm(xgrid, mean=10,sd=sqrt(0.2)) 
plot(xgrid,true_dens,xlab="X",ylab="Y",bty="l",type="l",xlim=c(-2,15),ylim=c(0,1), main="")
for(iter in burnin:niter){ lines(xgrid, F_mat_500_PY2[iter,1:200],col=iter,lty=3) }

true_dens<- 0.5*dnorm(xgrid, mean=1, sd= sqrt(0.2)) +0.5*dnorm(xgrid, mean=10,sd=sqrt(0.2)) 
plot(xgrid,true_dens,xlab="X",ylab="Y",bty="l",type="l",xlim=c(-2,15),ylim=c(0,1), main="")
for(iter in burnin:niter){ lines(xgrid, F_mat_500_NG2[iter,1:200],col=iter,lty=3) }

