


f   <- gjamSimData(S=120, Q=10, typeNames = 'PA')
rl <- list(r = 8, N = 50)
ml  <- list(ng = 200, burnin = 50, typeNames = 'PA', reductList = rl)
out <- gjam(f$formula, f$xdata, f$ydata, modelList = ml)

xdata_test<- f$xdata[1:20,]

new <- list(xdata =xdata_test,  nsim = 100) # effort unchanged 
p1  <- gjamPredict(output = out, newdata = new)








plot(f$ydata[1:20,], p1$sdList$yMu )
abline(0,1)



S   <- 200
f   <- gjamSimData(n = 100, S = S, typeNames='PA')
rl  <- list(r = 5, N = 20)
ml  <- list(ng = 2000, burnin = 500, typeNames = f$typeNames, 
            reductList = rl, holdoutIndex=1:10)  
#PREDICTX = F
out <- gjam(f$formula, f$xdata, f$ydata, modelList = ml)
pl  <- list(trueValues = f$trueValues, SMALLPLOTS = F)
gjamPlot(output = out, plotPars = pl)



