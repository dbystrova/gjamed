


load_object <- function(file) {
  tmp <- new.env()
  load(file = file, envir = tmp)
  tmp[[ls(tmp)[1]]]
}


setwd("~/Documents/GitHub/gjamed")
library(ggplot2)
LtGJ <- load_object("Sim_smallSK4_gjam.Rda")
Ltgj0<- load_object( "Sim_smallSK4_gjam0.Rda")
LtT1<-load_object("Sim_smallSK4_type1.Rda")
LtT2<-load_object("Sim_smallSK4_type2.Rda")
LtT3<-load_object("Sim_smallSK4_type3.Rda")
LtT4<-load_object("Sim_smallSK4_type4.Rda")
LtGJ <- load_object("ODSim_smallSK8_gjam.Rda")
Ltgj0<- load_object( "ODSim_smallSK8_gjam0.Rda")
LtT1<-load_object("ODSim_smallSK8_type1.Rda")
LtT2<-load_object("ODSim_smallSK8_type2.Rda")
LtT3<-load_object("ODSim_smallSK8_type3.Rda")
LtT4<-load_object("ODSim_smallSK8_type4.Rda")
 
S_vec<-c(20,50,80,100)
r_vec<-c(5,10,15,20)



table_comp<-as.data.frame(matrix(NA, nrow=length(r_vec)*length(S_vec)*5, ncol=1))
table_comp$S<- rep(S_vec, each=4*5)
table_comp$mod<- rep(c("gjam","gjam1","gjam2","gjam3","gjam4"), 4*4)
table_comp$num<- rep(1:16, each=5)
table_comp$rv<- rep(c(5,10,15,20), each=5, 4)
burn<-2000
it<- 5000
nsamples<-500
for(i in (1: nrow(table_comp))){
  j<-table_comp$num[i]
  table_comp$Schar[i]<- paste0("S=",table_comp$S[i])
  if((table_comp$mod[i]=="gjam")&(length(grep(paste0("r_",table_comp$rv[i]),names(LtGJ)[j]))>0)){ table_comp$res[i]<- mean(LtGJ[[j]]$trace[-c(1:burn)])}
  if((table_comp$mod[i]=="gjam1")&(length(grep(paste0("r_",table_comp$rv[i]),names(LtT1)[j]))>0)){ table_comp$res[i]<-  mean(LtT1[[j]]$trace[-c(1:burn)])}
  if((table_comp$mod[i]=="gjam2")&(length(grep(paste0("r_",table_comp$rv[i]),names(LtT2)[j]))>0)){ table_comp$res[i]<- mean(LtT2[[j]]$trace[-c(1:burn)])}
  if((table_comp$mod[i]=="gjam3")&(length(grep(paste0("r_",table_comp$rv[i]),names(LtT3)[j]))>0)){ table_comp$res[i]<-  mean(LtT3[[j]]$trace[-c(1:burn)])}
  if((table_comp$mod[i]=="gjam4")&(length(grep(paste0("r_",table_comp$rv[i]),names(LtT4)[j]))>0)){ table_comp$res[i]<- mean(LtT4[[j]]$trace[-c(1:burn)])}

  if((table_comp$mod[i]=="gjam")&(length(grep(paste0("r_",table_comp$rv[i]),names(Ltgj0)[j]))>0)){ 
    table_comp$lweight[i]<- apply(Ltgj0[[j]]$pk_chain[-c(1:burn),],2,mean)[length(apply(Ltgj0[[j]]$pk_chain,2,mean))]
    table_comp$fit_err[i]<- Ltgj0[[j]]$fit
    table_comp$err[i]<- Ltgj0[[j]]$err
    }
  if((table_comp$mod[i]=="gjam1")&(length(grep(paste0("r_",table_comp$rv[i]),names(LtT1)[j]))>0)){
    table_comp$lweight[i]<-  apply(LtT1[[j]]$pk_chain[-c(1:burn),],2,mean)[length(apply(LtT1[[j]]$pk_chain,2,mean))]
    table_comp$fit_err[i]<- LtT1[[j]]$fit
    table_comp$err[i]<- LtT1[[j]]$err
    }
  if((table_comp$mod[i]=="gjam2")&(length(grep(paste0("r_",table_comp$rv[i]),names(LtT2)[j]))>0)){
    table_comp$lweight[i]<- apply(LtT2[[j]]$pk_chain[-c(1:burn),],2,mean)[length(apply(LtT2[[j]]$pk_chain,2,mean))]
    table_comp$fit_err[i]<- LtT2[[j]]$fit
    table_comp$err[i]<- LtT2[[j]]$err
    }
  if((table_comp$mod[i]=="gjam3")&(length(grep(paste0("r_",table_comp$rv[i]),names(LtT3)[j]))>0)){ 
    table_comp$lweight[i]<-  apply(LtT3[[j]]$pk_chain[-c(1:burn),],2,mean)[length(apply(LtT3[[j]]$pk_chain,2,mean))]
    table_comp$fit_err[i]<- LtT3[[j]]$fit
    table_comp$err[i]<- LtT3[[j]]$err
    }
  if((table_comp$mod[i]=="gjam4")&(length(grep(paste0("r_",table_comp$rv[i]),names(LtT4)[j]))>0)){ 
    table_comp$lweight[i]<- apply(LtT4[[j]]$pk_chain[-c(1:burn),],2,mean)[length(apply(LtT4[[j]]$pk_chain,2,mean))]
    table_comp$fit_err[i]<- LtT4[[j]]$fit
    table_comp$err[i]<- LtT4[[j]]$err
    }

}

table_comp$Schar <- factor(table_comp$Schar, levels=c('S=20','S=50','S=80','S=100'))
table_comp<-table_comp[,-1]
p <- ggplot(data= table_comp)
q<- p +aes(x=rv,y= res, colour = mod) + geom_point(size = 2) +facet_wrap( ~ Schar, strip.position = "bottom", scales = "free_x")+ ylab("Posterior mean") + xlab("S and r values")+ylim(c(0,10))+
  labs(title="Ability to recover true number of groups for all models", caption=paste0("Number of iterations: ",it," burnin: ",burn," number of samples: ",nsamples))+
  geom_hline(yintercept = 8,color = "red")+theme_bw()+theme(panel.spacing = unit(0, "lines"),
                                                            strip.background = element_blank(),
                                                            strip.placement = "outside",legend.position = "top", plot.title = element_text(hjust = 0.5))
q


p <- ggplot(data= table_comp)
q<- p +aes(x=rv,y= lweight, colour = mod) + geom_point( size = 2) +facet_wrap( ~ Schar, strip.position = "bottom", scales = "free_x")+ ylab("Posterior mean") + xlab("S and r values")+ylim(c(0,0.3))+
  labs(title="Last p_k weight value for all models", caption=paste0("Number of iterations: ",it," burnin: ",burn," number of samples: ",nsamples))+theme_bw()+theme(panel.spacing = unit(0, "lines"),
                                                            strip.background = element_blank(),
                                                            strip.placement = "outside",legend.position = "top", plot.title = element_text(hjust = 0.5))
q






p <- ggplot(data= table_comp)
q<- p +aes(x=rv,y= fit_err, colour = mod) + geom_point( size = 2) +facet_wrap( ~ Schar, strip.position = "bottom", scales = "free_x")+ ylab("Posterior mean") + xlab("S and r values")+ylim(c(0,max(table_comp$fit_err)))+
  labs(title="Fit error value for all models", caption=paste0("Number of iterations: ",it," burnin: ",burn," number of samples: ",nsamples))+theme_bw()+theme(panel.spacing = unit(0, "lines"),
                                                                                                                                                          strip.background = element_blank(),                                                                                                                                                                    strip.placement = "outside",legend.position = "top", plot.title = element_text(hjust = 0.5))
q


p <- ggplot(data= table_comp)
q<- p +aes(x=rv,y= err, colour = mod) + geom_point( size = 2) +facet_wrap( ~ Schar, strip.position = "bottom", scales = "free_x")+ ylab("Posterior mean") + xlab("S and r values")+ylim(c(0,max(table_comp$err)))+
  labs(title="RMSE error between estimated and true covariance matrix for all models", caption=paste0("Number of iterations: ",it," burnin: ",burn," number of samples: ",nsamples))+theme_bw()+theme(panel.spacing = unit(0, "lines"),
                                                                                                                                                                    strip.background = element_blank(),                                                                                                                                                                    strip.placement = "outside",legend.position = "top", plot.title = element_text(hjust = 0.5))
q
########################################################################################################################






