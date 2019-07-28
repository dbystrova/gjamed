data(cls.draw2)
tru.class <- rep(1:8,each=50)
psm2 <- comp.psm(cls.draw2)
# optimize criteria based on PSM
mbind2 <- minbinder(psm2)
mpear2 <- maxpear(psm2)
# Relabelling
k <- apply(cls.draw2,1, function(cl) length(table(cl)))
max.k <- as.numeric(names(table(k))[which.max(table(k))])
relab2 <- relabel(cls.draw2[k==max.k,])
# compare clusterings found by different methods with true grouping
arandi(mpear2$cl, tru.class)
arandi(mbind2$cl, tru.class)
arandi(relab2$cl, tru.class)

vi.dist(mpear2$cl, tru.class)
vi.dist(mbind2$cl, tru.class)
vi.dist(relab2$cl, tru.class)




f<- load_object("")