.gjam00 <- function(formula, xdata, ydata, modelList){
  
  holdoutN      <-  0
  holdoutIndex  <- numeric(0)
  modelSummary  <- betaPrior  <- traitList <- effort <- NULL
  specByTrait   <- traitTypes <- breakList <- notStandard <- NULL
  censor <- censorCA <- censorDA <- CCgroups <- FCgroups <- intMat <- NULL
  reductList <- y0 <- N  <- r <- otherpar <- pg <- NULL
  ng     <- 2000
  burnin <- 500
  REDUCT <- TRAITS <- FULL <- F
  PREDICTX <- T
  lambdaPrior <- betaPrior <- NULL
  
  RANDOM <- F              # random group intercepts
  
  TIME <- F
  timeList <- timeZero <- timeLast <- timeIndex <- groupIndex <- 
    rowInserts <- Lmat <- Amat <- beta <- NULL
  
  ematAlpha <- .5
  
  # PY alpha.DP <- ncol(ydata)          # large values give more variation
  #alpha.DP <- 1
  
  #PY if(alpha.DP == 1)
  #PY   stop('multivariate model: at least 2 columns needed in ydata')
  
  for(k in 1:length(modelList))assign( names(modelList)[k], modelList[[k]] )
  
  if('CCgroups' %in% names(modelList))attr(typeNames,'CCgroups')  <- CCgroups
  if('FCgroups' %in% names(modelList))attr(typeNames,'FCgroups')  <- FCgroups
  if('CATgroups' %in% names(modelList))attr(typeNames,'CATgroups') <- CATgroups
  
  if(!is.null(timeList)){
    if("betaPrior" %in% names(timeList)){
      colnames(timeList$betaPrior$lo) <- 
        colnames(timeList$betaPrior$hi) <- 
        .cleanNames(colnames(timeList$betaPrior$lo))
    }
    if("lambdaPrior" %in% names(timeList)){
      colnames(timeList$lambdaPrior$lo) <- colnames(timeList$lambdaPrior$hi) <- 
        .cleanNames(colnames(timeList$lambdaPrior$lo))
    }
    
    for(k in 1:length(timeList))assign( names(timeList)[k], timeList[[k]] )
    TIME <- T
    REDUCT <- T
    BPRIOR <- T
    holdoutN      <-  0
    holdoutIndex  <- numeric(0)
  }
  
  if(!is.null(traitList)){
    TRAITS <- T
    for(k in 1:length(traitList))assign( names(traitList)[k], traitList[[k]] )
    
    stt <- .replaceString(colnames(specByTrait),'_','')
    colnames(specByTrait) <- stt
    colnames(plotByTrait) <- stt
    colnames(traitList$specByTrait) <- stt
    colnames(traitList$plotByTrait) <- stt
    modelList$traitList <- traitList
  }
  
  if(burnin >= ng) stop( 'burnin must be < no. MCMC steps, ng' )
  if('censor' %in% names(modelList)){
    for(k in 1:length(censor)){
      if( nrow(censor[[k]]$partition) != 3 )
        stop('censor matrix: 3 rows for value, lo, hi')
      rownames(censor[[k]]$partition) <- c('value','lo','hi')
    }
  }
  
  if(missing(xdata)) xdata <- environment(formula)
  
  S <- ncol(ydata)
  if(length(typeNames) == 1)typeNames <- rep(typeNames,S)
  if(length(typeNames) != S) 
    stop('typeNames must be one value or no. columns in y')
  
  ############### factors in y
  
  tmp <- .checkYfactor(ydata, typeNames)
  ydata <- tmp$ydata; yordNames <- tmp$yordNames
  
  if(TRAITS){
    if(!all( typeNames %in% c('CC','FC') ) )
      stop('trait prediction requires composition data (CC or FC)')
    if(nrow(plotByTrait) != nrow(ydata))
      stop('nrow(plotByTrait) must equal nrow(ydata)')
    if(ncol(plotByTrait) != length(traitTypes))
      stop('ncol(plotByTrait) must equal length(traitTypes)')
    if(ncol(plotByTrait) != length(traitTypes))
      stop('ncol(plotByTrait) must equal length(traitTypes)')
    ii <- identical(rownames(specByTrait),colnames(ydata))
    if(!ii){
      ww <- match(colnames(ydata),rownames(specByTrait) )
      if( is.finite(min(ww)) ){
        specByTrait <- specByTrait[ww,]
      } else {
        stop( 'rownames(specByTrait) must match colnames(ydata)' )
      }
    }
    if(typeNames[1] == 'CC'){
      ytmp <- round(ydata,0)
      ytmp[ytmp == 0 & ydata > 0] <- 1
      ydata <- ytmp
      rm(ytmp)
    }
  }
  
  tmp <- .buildYdata(ydata, typeNames)
  y   <- tmp$y
  ydataNames <- tmp$ydataNames
  typeNames  <- tmp$typeNames
  CCgroups   <- tmp$CCgroups
  FCgroups   <- tmp$FCgroups
  CATgroups  <- tmp$CATgroups
  if(TRAITS) rownames(specByTrait) <- colnames(y)
  
  S <- ncol(y)
  n <- nrow(y)
  
  cat("\nObservations and responses:\n")
  print(c(n, S))
  
  tmp    <- .buildEffort(y, effort, typeNames)
  effort <- tmp
  effMat <- effort$values
  modelList$effort <- effort
  re <- floor( diff( range(log10(effMat),na.rm=T) ) )
  if(re > 2)
    message(paste('sample effort > ', re, ' orders of magnitude--consider units near 1',sep='') )
  
  
  tmp      <- .gjamGetTypes(typeNames)
  typeCols <- tmp$typeCols
  typeFull <- tmp$typeFull
  typeCode <- tmp$TYPES[typeCols]
  allTypes <- sort(unique(typeCols))
  
  tmp <- .gjamXY(formula, xdata, y, typeNames, notStandard)
  x      <- tmp$x; y <- tmp$y; snames <- tmp$snames
  xdata  <- tmp$xdata; xnames <- tmp$xnames
  interBeta   <- tmp$interaction 
  factorBeta  <- tmp$factorAll
  designTable <- tmp$designTable;    xscale <- tmp$xscale
  predXcols   <- tmp$predXcols
  standMat    <- tmp$standMat;      standMu <- tmp$standMu  
  standRows   <- tmp$standRows;    
  xdataNames  <- tmp$xdataNames
  notStandard <- tmp$notStandard[tmp$notStandard %in% xnames]
  
  factorLambda <- interLambda <- NULL
  
  if(!is.null(lambdaPrior)){
    
    lformula <- attr(lambdaPrior$lo,'formula')
    
    tmp <- .gjamXY(lformula, xdata, y, typeNames, notStandard)
    xl   <- tmp$x
    mm <- match(colnames(xl),colnames(xdata))
    wm <- which(is.finite(mm))
    if(length(wm) > 0){
      xdata[,mm[wm]] <- xl[,wm]
    }
    
    xlnames <- tmp$xnames
    interLambda   <- tmp$interaction
    factorLambda <- tmp$factorAll
    
    designTable <- list(beta = designTable, lambda = tmp$designTable)
    
    standMatL    <- tmp$standMat;      standMuL <- tmp$standMu
    standRowsL    <- tmp$standRows;    
    notStandardL <- tmp$notStandard[tmp$notStandard %in% xlnames]
  }
  
  modelList     <- append(modelList, list('formula' = formula,
                                          'notStandard' = notStandard))
  
  Q <- ncol(x)
  
  tmp <- .gjamMissingValues(x, y, factorBeta$factorList, typeNames)
  xmiss  <- tmp$xmiss;   xbound <- tmp$xbound; 
  ymiss  <- tmp$ymiss;   missY <- tmp$missY
  xprior <- tmp$xprior;  yprior <- tmp$yprior
  nmiss  <- nrow(xmiss); mmiss  <- nrow(ymiss)
  x  <- tmp$x; y <- tmp$y
  
  if(TIME){
    tmp <- .gjamMissingValues(xl, y, factorLambda$factorList, typeNames)
    xlmiss  <- tmp$xmiss;   xlbound <- tmp$xbound; 
    xlprior <- tmp$xprior
    nlmiss  <- nrow(xmiss)
    xl <- tmp$x
  }
  
  reductList <- .setupReduct(modelList, S, Q, n) ##########
  N <- reductList$N; r <- reductList$r
  alpha.DP<-reductList$alpha.DP
  if(!is.null(reductList$N))REDUCT <- T
  
  
  
  
  tmp <- .gjamHoldoutSetup(holdoutIndex, holdoutN, n)
  holdoutIndex <- tmp$holdoutIndex; holdoutN <- tmp$holdoutN
  inSamples    <- tmp$inSamples;         nIn <- tmp$nIn
  
  tmp <- .gjamSetup(typeNames, x, y, breakList, holdoutN, holdoutIndex,
                    censor=censor, effort=effort) 
  w <- tmp$w; z <- tmp$z; y <- tmp$y; other <- tmp$other; cuts <- tmp$cuts
  cutLo       <- tmp$cutLo; cutHi <- tmp$cutHi; plo <- tmp$plo; phi <- tmp$phi
  ordCols     <- tmp$ordCols; disCols <- tmp$disCols; compCols <- tmp$compCols 
  conCols     <- which(typeNames == 'CON')
  classBySpec <- tmp$classBySpec; breakMat <- tmp$breakMat
  minOrd      <- tmp$minOrd; maxOrd <- tmp$maxOrd; censorCA <- tmp$censorCA
  censorDA    <- tmp$censorDA; censorCON <- tmp$censorCON; 
  ncut <- ncol(cuts);  corCols <- tmp$corCols
  catCols     <- which(attr(typeNames,'CATgroups') > 0)
  sampleW     <- tmp$sampleW
  ordShift    <- tmp$ordShift
  
  sampleW[censorCA] <- 1
  sampleW[censorDA] <- 1
  sampleW[censorCON] <- 1
  sampleWhold <- tgHold <- NULL
  wHold <- NULL
  wmax  <- apply(y/effMat,2,max,na.rm=T)
  pmin  <- -2*abs(wmax)
  
  if(mmiss > 0){
    phi[ ymiss ] <- wmax[ ymiss[,2] ]
    plo[ ymiss ] <- pmin[ ymiss[,2] ]
    sampleW[ ymiss ] <- 1
  }
  
  ploHold <- phiHold <- NULL
  
  if(holdoutN > 0){
    sampleWhold <- sampleW[holdoutIndex,]  #to predict X
    sampleW[holdoutIndex,] <- 1
    tgHold  <- cuts
    wHold   <- w[drop=F,holdoutIndex,]
    
    ploHold <- plo[drop=F,holdoutIndex,]   # if LOHI: updated to current yp
    phiHold <- phi[drop=F,holdoutIndex,]
  }
  
  byCol <- byRow <- F
  if(attr(sampleW,'type') == 'cols')byCol <- T
  if(attr(sampleW,'type') == 'rows')byRow <- T
  indexW <- attr(sampleW,'index')
  
  notCorCols <- c(1:S)
  if(length(corCols) > 0)notCorCols <- notCorCols[-corCols]
  
  ############ 'other' columns
  sigmaDf  <- nIn - Q + S - 1
  sg <- diag(.1,S)
  SO <- S
  
  notOther <- c(1:S)
  sgOther  <- NULL
  if(length(other) > 0){                     
    notOther   <- notOther[!notOther %in% other]
    SO         <- length(notOther)
    sg[other,] <- sg[,other] <- 0
    sgOther    <- matrix( cbind(other,other),ncol=2 )
    sg[sgOther] <- .1
  }
  
  ############## prior on beta
  loB <- hiB <- NULL
  beta <- bg <- matrix(0,Q,S)
  rownames(beta) <- colnames(x)
  BPRIOR <- F
  
  if( !is.null(betaPrior) ){
    colnames(betaPrior$lo) <- .cleanNames(colnames(betaPrior$lo))
    colnames(betaPrior$hi) <- .cleanNames(colnames(betaPrior$hi))
    loB <- betaPrior$lo
    hiB <- betaPrior$hi
    
    bg <- (loB + hiB)/2
    bg[is.nan(bg)] <- 0
    
    wB <- which(!is.na(t(loB[,notOther])), arr.ind=T)[,c(2,1)]
    wB <- rbind(wB, which(!is.na(t(hiB[,notOther])), arr.ind=T)[,c(2,1)])
    colnames(wB) <- c('row','col')
    
    tmp <- .betaPrior(bg, notOther, loB, hiB)
    bg <- tmp$beta; loB <- tmp$loB; hiB <- tmp$hiB
    wB <- tmp$wB; BPRIOR <- tmp$BPRIOR
    bg[is.nan(bg)] <- 0
    
    tmp <- .getPattern(bg[,notOther], wB)
    Brows <- tmp$rows
    Bpattern <- tmp$pattern
    BPRIOR <- T
    bg[!is.finite(bg)] <- 0
  }
  
  zeroBeta <- .factorCoeffs2Zero(factorBeta, snames, betaPrior)  # max zero is missing factor level
  zeroLambda <- NULL
  
  ############### time 
  if( TIME ){
    
    BPRIOR <- T
    
    tmp <- .getTimeIndex(timeList, other, notOther, xdata, x, xl, y, w)
    Lmat   <- tmp$Lmat; Lpattern <- tmp$Lpattern;  wL <- tmp$wL
    Vmat   <- tmp$Vmat;    Lrows <- tmp$Lrows; gindex <- tmp$gindex
    loLmat <- tmp$loLmat; hiLmat <- tmp$hiLmat; Arows <- tmp$Arows
    Amat   <- tmp$Amat; Apattern <- tmp$Apattern; wA <- tmp$wA
    Umat   <- tmp$Umat;   uindex <- tmp$uindex
    loAmat <- tmp$loAmat; hiAmat <- tmp$hiAmat; aindex <- tmp$aindex
    Brows  <- tmp$Brows;      bg <- tmp$bg; Bpattern <- tmp$Bpattern
    wB     <- tmp$wB;        loB <- tmp$loB; hiB <- tmp$hiB
    timeZero <- tmp$timeZero; timeLast <- tmp$timeLast
    maxTime  <- tmp$maxTime; inSamples <- tmp$inSamples 
    tindex   <- tmp$tindex; sindex <- tmp$sindex; i1 <- tmp$i1; i2 <- tmp$i2
    
    if(is.null(loB))BPRIOR <- F
    
    Unew <- Umat
    Vnew <- Vmat
    mua  <- mub <- mug <- muw <- w*0
    
    zeroLambda <- .factorCoeffs2Zero(factorLambda, snames, lambdaPrior)
    timeList$lambdaPrior$hi[zeroLambda] <- lambdaPrior$hi[zeroLambda] <- 0
    timeList$betaPrior$hi[zeroBeta]     <- betaPrior$hi[zeroBeta] <- 0
    
    standMatLmat <- Lmat*0
    notStandardLmat <- numeric(0)
    
    if(length(standRowsL) > 0){
      csl <- paste('_',names(standRowsL),sep='')
      for(j in 1:length(csl)){
        wj <- grep(csl[j],rownames(Lmat))
        standMatLmat[wj,] <- standMatL[standRowsL[j],]
        notStandardLmat <- c(notStandardLmat,wj)
      }
    }
  } 
  
  if(byCol){
    inw <- intersect( colnames(y)[indexW], colnames(y)[notOther] )
    indexW <- match(inw,colnames(y)[notOther])
  }
  
  IXX <- NULL
  if(nmiss == 0){
    XX    <- crossprod(x)
    IXX <- chol2inv(chol( XX ) )
  }
  
  
  updateBeta <- .betaWrapper(REDUCT, TIME, BPRIOR, notOther, IXX, 
                             betaLim=max(wmax)/2)
  
  ############ dimension reduction
  
  inSamp <- inSamples
  if(TIME)inSamp <- tindex[,1]     # index for x
  
  CLUST <- T   # dirichlet 
  
  .param.fn <- .paramWrapper00(REDUCT, inSamp, SS=length(notOther))
  sigmaerror <- .1
  otherpar   <- list(S = S, Q = Q, sigmaerror = sigmaerror, 
                     Z = NA, K =rep(1,S), sigmaDf = sigmaDf)
  sigErrGibbs <- rndEff <- NULL
  
  yp <- y
  wmax <- ymax <- apply(y,2,max)
  wmax <- wmax/effMat
  
  if(REDUCT){
    cat( paste('\nDimension reduced from',S,'X',S,'->',N,'X',r,'responses\n') )
    otherpar$N <- N; otherpar$r <- r; otherpar$sigmaerror <- 0.1
    otherpar$Z <- rmvnormRcpp(N,rep(0,r),1/S*diag(r))
    otherpar$D <- .riwish(df = (2 + r + N), 
                          S = (crossprod(otherpar$Z) +
                                 2*2*diag(rgamma(r,shape=1,rate=0.001))))
    otherpar$K <- sample(1:N,length(notOther),replace=T)
    
    otherpar$alpha.DP <- alpha.DP
    
    otherpar$pvec     <- .sampleP(N=N, avec=rep(1,(N-1)),
                                  bvec=rep(alpha.DP,(N-1)), K=otherpar$K)
    kgibbs <- matrix(1,ng,S)
    sgibbs <- matrix(0,ng, N*r)
    nnames <- paste('N',1:N,sep='-')
    rnames <- paste('r',1:r,sep='-')
    colnames(sgibbs) <- .multivarChainNames(nnames,rnames)
    sigErrGibbs <- rep(0,ng)   
    
    rndEff <- w*0
    
  } else {
    Kindex <- which(as.vector(lower.tri(diag(S),diag=T)))
    nK     <- length(Kindex)
    sgibbs <- matrix(0,ng,nK)
    colnames(sgibbs) <- .multivarChainNames(snames,snames)[Kindex] # half matrix
  }
  
  out <- .param.fn(CLUST=T, x, beta = bg[,notOther], Y = w[,notOther], otherpar)  
  sg[notOther,notOther]    <- out$sg
  otherpar      <- out$otherpar
  
  muw <- w
  
  if(!TIME){
    
    Y <- w[inSamp,notOther]
    sig <- sg[notOther,notOther]
    
    if(REDUCT){
      Y <- Y - rndEff[inSamp,notOther]
      sig <- sigmaerror
    }
    
    bg[,notOther] <- updateBeta(X = x[inSamp,], Y, sig, beta = bg[,notOther],
                                loB, hiB)
    muw <- x%*%bg
    
  }else{
    
    mua <- Umat%*%Amat
    mug <- Vmat%*%Lmat
    Y <- w - mua - mug - rndEff
    
    if(REDUCT){
      sig <- sigmaerror
    }else{ sig <- sg[notOther,notOther] }
    
    bg[,notOther] <- updateBeta(X = x[tindex[,2],], Y = Y[tindex[,2],notOther], 
                                sig = sig, beta = bg[,notOther], 
                                lo = loB[,notOther], hi = hiB[,notOther], 
                                rows=Brows, pattern=Bpattern)
    mub <- x%*%bg
    muw <- mub + mug + mua 
    wpropTime <- .001 + .1*abs(w)
  }
  
  sg[other,] <- sg[,other] <- 0
  diag(sg)[other]          <- 1
  rownames(bg)  <- xnames
  rownames(sg)  <- colnames(sg) <- colnames(bg) <- snames
  colnames(x)   <- xnames
  
  
  
  ############ ordinal data
  
  cutg <- tg <- numeric(0)
  
  if('OC' %in% typeCode){
    tg       <- cutg <- cuts
    cnames   <- paste('C',1:ncut,sep='-')
    nor      <- length(ordCols)
    cgibbs   <- matrix(0,ng,(ncut-3)*nor)
    colnames(cgibbs) <- as.vector( outer(snames[ordCols],
                                         cnames[-c(1,2,ncut)],paste,sep='_') )
    tmp   <- .gjamGetCuts(y+1,ordCols)
    cutLo <- tmp$cutLo
    cutHi <- tmp$cutHi
    plo[,ordCols] <- tg[cutLo]                                        
    phi[,ordCols] <- tg[cutHi]
    lastOrd <- ncol(tg)
  }
  
  ############ setup w
  tmp <- .gjamGetTypes(typeNames)
  typeFull <- tmp$typeFull
  typeCols <- tmp$typeCols
  allTypes <- unique(typeCols)
  Y <- w
  
  LOHI <- F
  if(!LOHI & holdoutN > 0){
    minlo <- apply(plo,2,min)
    minlo[minlo > 0] <- 0
    maxhi <- apply(phi,2,max)
  }
  
  if(!TIME){
    .updateW <- .wWrapper(REDUCT, RANDOM, S, effMat, corCols, notCorCols, typeNames, 
                          typeFull, typeCols, 
                          allTypes, holdoutN, holdoutIndex, censor, 
                          censorCA, censorDA, censorCON, notOther, sampleW, 
                          byRow, byCol,
                          indexW, ploHold, phiHold, sampleWhold, inSamp)
  }else{
    
    .updateW <- .wWrapperTime(sampleW, y, timeZero, i1, i2, tindex, gindex,
                              uindex, notOther, n, S, REDUCT, RANDOM)
    Y <- w - mua -  mug - rndEff
  }
  
  ycount <- rowSums(y)
  if('CC' %in% typeCode)ycount <- rowSums(y[,compCols])
  
  ############ X prediction
  
  tmp <- .xpredSetup(Y, x, bg, interBeta$isNonLinX, factorBeta, 
                     factorBeta$intMat, 
                     standMat, standMu, notOther, notStandard) 
  factorBeta$linFactor <- tmp$linFactor; xpred <- tmp$xpred; px <- tmp$px
  lox <- tmp$lox; hix <- tmp$hix
  
  priorXIV  <- diag(1e-5,ncol(x))
  priorX    <- colMeans(x)
  priorX[abs(priorX) < 1e-10] <- 0
  
  linFactor <- NULL
  
  ################## random groups
  
  if('random' %in% names(modelList)){
    
    RANDOM <- T
    rname  <- modelList$random
    randGroupTab <- table( as.character(xdata[,rname]) )
    
    wss <- names(randGroupTab[randGroupTab <= 2])
    if(length(wss) > 0){
      xdata[,rname] <- .combineFacLevels(xdata[,rname], fname=wss, 
                                         aname = 'rareGroups', vminF=1)
      randGroupTab <- table( as.character(xdata[,rname]) )
    }
    
    randGroups <- names( randGroupTab )
    G <- length(randGroups)
    
    groupIndex  <- match(as.character(xdata[,rname]),randGroups)
    rmm <- matrix(groupIndex,length(groupIndex), S)
    smm <- matrix(1:S, length(groupIndex), S, byrow=T)
    
    randGroupIndex <- cbind( as.vector(smm), as.vector(rmm) )
    colnames(randGroupIndex) <- c('species','group')
    xdata[,rname] <- as.factor(xdata[,rname])
    alphaRandGroup <- matrix(0, S, G)
    rownames(alphaRandGroup) <- snames
    colnames(alphaRandGroup) <- randGroups
    Cmat <- var(w[,notOther]/2)
    Cmat <- Cmat + diag(.1*diag(Cmat))
    Cprior <- Cmat
    CImat <- solve(Cprior)
    Ckeep <- diag(S)
    
    alphaRanSums <- alphaRandGroup*0
    groupRandEff <- w*0
    
    Kindex <- which(as.vector(lower.tri(diag(S),diag=T)))
    nK     <- length(Kindex)
    alphaVarGibbs <- matrix(0,ng,nK)
    colnames(alphaVarGibbs) <- .multivarChainNames(snames,snames)[Kindex] # half matrix
  }
  
  
  ################################## XL prediction: variables in both
  
  if(TIME){
    
    tmp <- .xpredSetup(Y, xl, lambdaPrior$lo, 
                       interLambda$isNonLinX, factorLambda, interLambda$intMat, standMatL, 
                       standMuL, notOther, notStandardL) 
    factorLambda$linFactor <- tmp$linFactor
    lox <- c(lox,tmp$lox[!names(tmp$lox) %in% names(lox)])
    hix <- c(hix,tmp$lox[!names(tmp$hix) %in% names(hix)])
    
    ################ or
    xpred <- cbind(xpred,xl[,!colnames(xl) %in% colnames(x)])
    Qall <- ncol(xpred) - 1
    
    intMat <- numeric(0)
    if( length(interBeta$intMat) > 0 ){
      intMat <- match(xnames[interBeta$intMat],colnames(xpred))
      intMat <- matrix(intMat,nrow(interBeta$intMat),3)
    }
    if( length(interLambda$intMat) > 0){
      ib <- match(xlnames[interLambda$intMat],colnames(xpred))
      ib <- matrix(ib,nrow(interLambda$intMat),3)
      intMat <- rbind(intMat,ib)
    }
    
    linFactor <- numeric(0)
    lf <- factorBeta$linFactor
    if( length(lf) > 0 ){
      for(k in 1:length(lf)){
        kf <- match(xnames[lf[[k]]],colnames(xpred))
        linFactor <- append(linFactor,list(kf))
      }
    }
    lf <- factorLambda$linFactor
    if( length(lf) > 0 ){
      for(k in 1:length(lf)){
        kf <- match(xlnames[lf[[k]]],colnames(xpred))
        linFactor <- append(linFactor,list(kf))
      }
    }
  }
  
  ############  contrasts, predict F matrix
  
  tmp <- .setupFactors(xdata, xnames, factorBeta)
  ff  <- factorBeta[names(factorBeta) != 'factorList']
  factorBeta <- append(ff,tmp)
  
  ############ E matrix
  emat <- matrix(0,S,S)
  colnames(emat) <- rownames(emat) <- snames
  lo <- hi <- lm <- hm <- ess <- emat
  
  fmat <- factorBeta$fmat
  fnames <- rownames( factorBeta$lCont )
  q2 <- nrow(fmat)
  
  if(TIME){
    tmp <- .setupFactors(xdata, xlnames, factorLambda)
    ff <- factorLambda[names(factorLambda) != 'factorList']
    if(length(tmp) > 0)factorLambda <- append(ff,tmp)
    factorLambda$LCONT <- rep(TRUE, factorLambda$nfact)
    flnames <- rownames( factorLambda$lCont )
    
    ############ E matrix TIME
    ematL <- matrix(0,S,S)
    colnames(ematL) <- rownames(ematL) <- snames
    essL <- ematL
  }
  
  ############ sp richness
  richness <- richFull <- NULL
  RICHNESS <- F
  
  inRichness <- which(!typeNames %in% c('CON','CAT','OC'))
  inRichness <- inRichness[!inRichness %in% other]
  if(length(inRichness) > 2)RICHNESS  <- T
  
  wrich <- y*0 
  wrich[,inRichness] <- 1
  wrich[ymiss] <- 0
  
  presence <- w*0
  
  covx <- cov(x)
  
  ############ sums
  predx  <- predx2 <- xpred*0
  yerror <- ypred  <- ypred2 <- wpred  <- wpred2 <- ymissPred <- ymissPred2 <- y*0
  sumDev <- 0   #for DIC
  sMean  <- sg*0
  ntot   <- 0
  
  if(nmiss > 0){
    xmissSum <- xmissSum2 <- rep(0,nmiss)
  }
  
  if(TIME)predxl <- predxl2 <- xl*0
  
  ############ gibbs chains
  
  q2 <- length(fnames)
  fSensGibbs <- matrix(0,ng,q2)
  colnames(fSensGibbs) <- fnames
  
  bFacGibbs <- matrix(0,ng,q2*SO)
  colnames(bFacGibbs) <- .multivarChainNames(fnames,snames[notOther])
  
  bgibbs <- matrix(0,ng,S*Q)
  colnames(bgibbs) <- .multivarChainNames(xnames,snames)
  bgibbsUn <- bgibbs                   # unstandardized
  
  covE <- cov( x%*%factorBeta$dCont )  # note that x is standardized
  
  if(TRAITS){
    
    specTrait <- specByTrait[colnames(y),]
    tnames    <- colnames(specTrait)
    M         <- ncol(specTrait)
    specTrait <- t(specTrait)
    
    tpred  <- tpred2 <- matrix(0,n,M)
    
    missTrait <- which(is.na(specTrait),arr.ind=T)
    if(length(missTrait) > 0){
      traitMeans <- rowMeans(specTrait,na.rm=T)
      specTrait[missTrait] <- traitMeans[missTrait[,2]]
      warning( paste('no. missing trait values:',nrow(missTrait)) )
    }
    
    bTraitGibbs <- matrix(0,ng,M*Q)
    colnames(bTraitGibbs) <- .multivarChainNames(xnames,tnames)
    
    bTraitFacGibbs <- matrix(0,ng,q2*M)
    colnames(bTraitFacGibbs) <- .multivarChainNames(fnames,tnames)
    
    mgibbs <- matrix(0,ng,M*M)
    colnames(mgibbs) <- .multivarChainNames(tnames,tnames)
  }
  
  if(TIME){
    
    yy <- y*0
    yy[rowInserts,] <- 1
    ymiss <- which(yy == 1, arr.ind=T)
    rm(yy)
    mmiss <- length(ymiss)
    
    covL <- cov( xl%*%factorLambda$dCont )  # note x is standardized
    
    ggibbs <- matrix(0,ng,nrow(wL))
    colnames(ggibbs) <- rownames(wL)
    
    wnames <- apply(wA,1,paste0,collapse='-')  #locations in Amat, not alpha
    alphaGibbs <- matrix(0,ng,nrow(wA))
    colnames(alphaGibbs) <- wnames
    
    nl <- nrow(lambda)
    lgibbs <- matrix(0,ng,length(lambda[,notOther]))
    colnames(lgibbs) <- .multivarChainNames(xlnames,snames[notOther])
    
    gsensGibbs <- matrix(0,ng,nl)
    colnames(gsensGibbs) <- rownames(lambda)
    
    asensGibbs <- matrix(0,ng,nrow(Amat))
    colnames(asensGibbs) <- rownames(Amat)
    
    ni  <- length(i1)
    
    nA <- nrow(wA)
    nL <- nrow(wL)
    
    spA <- rep(.001, nA)
    spL <- rep(.01, nL)
    g1 <- 1
    gcheck <- c(50, 100, 200, 400, 800)
    tinyg <- 1e-6
  }
  
  pbar <- txtProgressBar(min=1,max=ng,style=1)
  
  # unstandardize
  
  tmp <- .getUnstandX(x, standRows, standMu[,1],standMat[,1],
                      interBeta$intMat)
  S2U      <- tmp$S2U
  xUnstand <- tmp$xu
  
  if(TIME){
    tmp <- .getUnstandX(xl, standRowsL, standMuL[,1],standMatL[,1],
                        interLambda$intMat)
    S2UL      <- tmp$S2U
    xlUnstand <- tmp$xu
  }
  
  if(REDUCT){
    rndTot <- w*0 
  }
  notPA <- which(!typeNames == 'PA' & !typeNames == 'CON')
  
  
  if(length(y) < 10000 | FULL) FULL <- T
  
  if(FULL){
    ygibbs <- matrix(0,ng,length(y))
  }
  if(RICHNESS){
    ypredPres <- ypredPres2 <- ypredPresN <- y*0
    shannon   <- rep(0,n)
  }
  
  for(g in 1:ng){ ########################################################
    
    if(REDUCT){
      
      #   if(g > burnin)CLUST <- F
      
      Y <- w[,notOther]
      if(RANDOM)Y <- Y - groupRandEff[,notOther] 
      if(TIME)  Y <- Y - mua[,notOther] - mug[,notOther] 
      
      tmp <- .param.fn(CLUST=T, x, beta = bg[,notOther], Y = Y, otherpar)
      sg[notOther,notOther] <- tmp$sg
      otherpar            <- tmp$otherpar
      rndEff[,notOther]   <- tmp$rndEff
      sigmaerror          <- otherpar$sigmaerror
      kgibbs[g,notOther]  <- otherpar$K
      sgibbs[g,]          <- as.vector(otherpar$Z)
      sigErrGibbs[g]      <- sigmaerror
      
      if(length(corCols) > 0){
        if(max(diag(sg)[corCols]) > 5){  #overfitting covariance
          stop(
            paste('\noverfitted covariance, reductList$N = ',N, 
                  'reductList$r = ',r, '\nreduce N, r\n')
          )
        }
      }
      
      sg[sgOther]         <- .1*sigmaerror
      
      sinv <- .invertSigma(sg[notOther,notOther],sigmaerror,otherpar,REDUCT)
      sdg  <- sqrt(sigmaerror)
      
      if(!TIME){
        Y <- w[inSamp,notOther] - rndEff[inSamp,notOther]
        if(RANDOM)Y <- Y - groupRandEff[inSamp,notOther]
        bg[,notOther] <- updateBeta(X = x[inSamp,], Y, 
                                    sig = sigmaerror, beta = bg[,notOther],
                                    lo=loB[,notOther], hi=hiB[,notOther])
        muw[inSamp,] <- x[inSamp,]%*%bg
        
      } else {
        
        mua  <- Umat%*%Amat 
        mug  <- Vmat%*%Lmat
        
        Y <- w[,notOther] - mua[,notOther] - mug[,notOther] - rndEff[,notOther]
        if(RANDOM)Y <- Y - groupRandEff[,notOther]
        bg[,notOther] <- updateBeta(X = x[tindex[,2],], Y = Y[tindex[,2],], 
                                    sig = sigmaerror, beta = bg[,notOther],
                                    rows = Brows, pattern = Bpattern,
                                    lo=loB[,notOther], hi=hiB[,notOther])
        mub <- x%*%bg
        Y   <- w - mub - mua - rndEff
        if(RANDOM)Y <- Y - groupRandEff
        
        Lmat[,notOther] <- updateBeta(X = Vmat[tindex[,2],], 
                                      Y = Y[tindex[,2],notOther], sig=sigmaerror, 
                                      beta = Lmat[,notOther],
                                      rows = Lrows, pattern = Lpattern, 
                                      lo=loLmat, hi=hiLmat, ixx=F)
        
        #     Lmat[,notOther] <- .updateBetaMet(X = Vmat[tindex[,2],], 
        #                                       Y[tindex[,2],notOther], 
        #                                       B = Lmat[,notOther],
        #                           lo=loLmat, hi=hiLmat, loc = wL, REDUCT, 
        #                           sig=sigmaerror,sp=spL)
        mug  <- Vmat%*%Lmat
        Y    <- w - mub - mug - rndEff
        if(RANDOM)Y <- Y - groupRandEff
        Amat <- updateBeta(X = Umat[tindex[,2],], Y[tindex[,2],], sig=sigmaerror, 
                           rows = Arows, pattern = Apattern, 
                           beta = Amat,
                           lo=loAmat, hi=hiAmat, ixx=F)
        #     Amat <- .updateBetaMet(X = Umat[tindex[,2],], Y[tindex[,2],notOther], 
        #                                       B = Amat,
        #                                       lo=loAmat, hi=hiAmat, loc = wA, REDUCT, 
        #                                      sig=sigmaerror,sp=rexp(nA,1/spA))
        mua <- Umat%*%Amat
        
        #     if(g %in% gcheck){
        #       g2   <- g - 1
        #       spA <- apply(alphaGibbs[g1:g2,],2,sd)/2 + tinyg
        #       spL <- apply(ggibbs[g1:g2,],2,sd)/2 + tinyg
        #       if(g < 200)g1 <- g
        #     }
        
        muw <- mub + mug + mua + rndEff
      }
      
    } else {
      Y <- w[inSamp,notOther]
      if(RANDOM)Y <- Y - groupRandEff[inSamp,notOther]
      bg[,notOther] <- updateBeta(X = x[inSamp,], Y, 
                                  sig = sg[notOther,notOther], 
                                  beta = bg[,notOther], 
                                  lo=loB, hi=hiB)
      
      muw[inSamp,] <- x[inSamp,]%*%bg
      
      SS   <- crossprod(w[inSamp,] - muw[inSamp,])
      SI   <- solveRcpp(SS[notOther,notOther])
      sinv <- .rwish(sigmaDf,SI)
      
      sg[notOther,notOther] <- solveRcpp(sinv)
      sgibbs[g,] <- sg[Kindex]
    }
    
    # muw does not include rndEff or groupRandEff
    
    alphaB <- .sqrtRootMatrix(bg,sg,DIVIDE=T)
    
    if( 'OC' %in% typeCode ){
      tg <- .updateTheta(w,tg,cutLo,cutHi,ordCols,
                         holdoutN,holdoutIndex,minOrd,maxOrd) # var scale
      cutg <- .gjamCuts2theta(tg,ss = sg[ordCols,ordCols]) # corr scale
      breakMat[ordCols,1:lastOrd] <- cutg
      cgibbs[g,] <- as.vector( cutg[,-c(1,2,ncut)] )
      
      plo[,ordCols] <- cutg[cutLo]
      phi[,ordCols] <- cutg[cutHi]
    }
    
    if(RANDOM){
      
      cw <- w - muw
      if(REDUCT){
        cw <- cw - rndEff
        v  <- 1/sigmaerror*.byGJAM(as.vector(cw), randGroupIndex[,1], 
                                   randGroupIndex[,2], alphaRandGroup*0, 
                                   fun='sum')[notOther,]
        sinv <- diag(1/sigmaerror, SO)
      }else{
        v <- .byGJAM(as.vector(cw), randGroupIndex[,1], 
                     randGroupIndex[,2], alphaRandGroup*0, fun='sum')[notOther,]
        v <- sinv%*%v
      }
      
      alphaRandGroup[notOther,] <- randEffRcpp(v, randGroupTab, 
                                               sinv, CImat)
      if(length(other) > 0)alphaRandGroup[other,] <- 0
      if(g < 100){
        alphaRandGroup[notOther,] <- 
          sweep( alphaRandGroup[notOther,], 2, 
                 colMeans(alphaRandGroup[notOther,]), '-')
      }
      SS  <- crossprod(t(alphaRandGroup[notOther,]))
      SS  <- S*SS + Cmat
      
      testv <- try( chol(SS) ,T)
      if( inherits(testv,'try-error') ){
        tiny  <- .1*diag(SS)
        SS  <- SS + diag(diag(SS + tiny))
      }
      
      Ckeep[notOther,notOther] <- .riwish( df = S*G + 1, SS )
      CImat <- solveRcpp(Ckeep[notOther,notOther])
      
      alphaVarGibbs[g,] <- Ckeep[Kindex]
      groupRandEff <- t(alphaRandGroup)[groupIndex,]
    }
    
    if(TIME){
      
      #    muw does not include groupRandEff
      tmp <- .updateW(w,plo,phi,wpropTime,xl,yp,Lmat,Amat,mub,rndEff, groupRandEff,
                      sdg,muw,Umat,Vmat,sinv)
      w <- tmp$w; muw <- tmp$muw; yp <- tmp$yp; Umat <- tmp$Umat; Vmat <- tmp$Vmat
      
      groups <- NULL
      
      for(k in allTypes){
        
        wk <- which(typeCols == k)
        nk <- length(wk)
        wo <- which(wk %in% notOther)
        wu <- which(typeCols[notOther] == k)
        wp <- w[, wk, drop=F]
        yp <- yp[, wk, drop=F]
        
        if(typeFull[wk[1]] == 'countComp')groups <- CCgroups
        if(typeFull[wk[1]] == 'fracComp')groups  <- FCgroups
        if(typeFull[wk[1]] == 'categorical')groups <- CATgroups
        
        glist <- list(wo = wo, type = typeFull[wk[1]], yy = y[,wk,drop=F], 
                      wq = wp, yq = yp, cutg = cutg, 
                      censor = censor, censorCA = censorCA, 
                      censorDA = censorDA, censorCON  = censorCON, 
                      eff = effMat[,wk,drop=F], groups = groups, 
                      k = k, typeCols = typeCols, notOther = notOther, 
                      wk = wk, sampW = sampleW[,wk])
        tmp <- .gjamWLoopTypes( glist )
        w[,wk]  <- tmp[[1]]
        yp[,wk] <- tmp[[2]]
      }
      
      #predict X
      
      ww <- w
      ww[ww < 0] <- 0
      mua  <- Umat%*%Amat
      mug  <- Vmat%*%Lmat
      muw <- mua + mub + mug + rndEff
      
      xtmp <- xpred
      xtmp[,-1] <- .tnorm(n*Qall,-3,3,xpred[,-1],.1)
      
      # factors
      if( length(linFactor) > 0 ){
        
        for(k in 1:length(linFactor)){
          
          mm  <- linFactor[[k]]
          wcol <- sample(mm,n,replace=T)
          xtmp[,mm[-1]] <- 0
          xtmp[ cbind(1:n, wcol) ] <- 1
          
        }
      }
      
      if(length(intMat) > 0){     #  interactions
        xtmp[,intMat[,1]] <- xtmp[,intMat[,2]]*xtmp[,intMat[,3]]
      }
      
      ae     <- mua + rndEff
      Vnow   <- Vmat
      mubNow <- xpred[,xnames]%*%bg
      mubNew <- xtmp[,xnames]%*%bg
      
      Vnow[tindex[,2],] <- ww[tindex[,1],gindex[,'colW']]*
        xpred[tindex[,2],xlnames][,gindex[,'rowG']]
      Vnow[timeZero+1,] <- ww[timeZero,gindex[,'colW']]*
        xpred[timeZero+1,xlnames][,gindex[,'rowG']]
      mugNow <- Vnow%*%Lmat
      muNow  <- mubNow + mugNow + ae
      
      Vnew[tindex[,2],] <- ww[tindex[,1],gindex[,'colW']]*
        xtmp[tindex[,2],xlnames][,gindex[,'rowG']]
      Vnew[timeZero+1,] <- ww[timeZero,gindex[,'colW']]*
        xtmp[timeZero+1,xlnames][,gindex[,'rowG']]
      mugNew <- Vnew%*%Lmat
      muNew  <- mubNew + mugNew + ae
      
      if(REDUCT){
        pnow <- dnorm(w[,notOther],muNow[,notOther],sdg,log=T)
        pnew <- dnorm(w[,notOther],muNew[,notOther],sdg,log=T)
        a1   <- exp( rowSums(pnew - pnow) )
      }else{
        pnow <- .dMVN(w[tindex[,2],notOther],muNow,sinv=sinv,log=T) 
        pnew <- .dMVN(w[tindex[,2],notOther],muNew,sinv=sinv,log=T) 
        a1   <- exp(pnew - pnow)
      }
      z    <- runif(length(a1),0,1)
      za   <- which(z < a1)
      if(length(za) > 0){
        xpred[za,] <- xtmp[za,]
        Vmat[za,] <- Vnew[za,]
        muw[za,]  <- muNew[za,]
        mub[za,]  <- mubNew[za,]
        mug[za,]  <- mugNew[za,]
      }
      
      if(nlmiss > 0)xl[xlmiss] <- xpred[xmiss]
      
      if(nmiss > 0){
        
        x[xmiss] <- xpred[xmiss]
        
        tmp    <- .getUnstandX(x, standRows, standMu[,1],
                               standMat[,1], intMat)            
        S2U    <- tmp$S2U
        XX     <- crossprod(x)
        IXX    <- solveRcpp(XX)
      }
      
      ggibbs[g,]     <- Lmat[wL]
      alphaGibbs[g,] <- Amat[wA]
      
    } else{ #############not TIME
      
      tmp   <- .updateW( rows=1:n, x, w, y, bg, sg, alpha=alphaB, 
                         cutg, plo, phi, rndEff, groupRandEff, 
                         sigmaerror, wHold )
      w     <- tmp$w
      yp    <- tmp$yp
      plo   <- tmp$plo
      phi   <- tmp$phi
      wHold <- tmp$wHold    #values for w if not held out
      
      
      Y <- w[,notOther]
      if(holdoutN > 0) Y[holdoutIndex,] <- wHold[,notOther]  # if w not held out
      if(RANDOM)Y <- Y - groupRandEff[,notOther]
      
      if(nmiss > 0){
        
        x[xmiss] <- .imputX_MVN(x,Y,bg[,notOther],xmiss,sinv,xprior=xprior,
                                xbound=xbound)[xmiss]
        tmp      <- .getUnstandX(x, standRows, standMu[,1],
                                 standMat[,1], intMat)            
        S2U    <- tmp$S2U
        XX     <- crossprod(x)
        IXX    <- solveRcpp(XX)
      }
      
      if( PREDICTX & length(predXcols) > 0){
        
        if( length(interBeta$isNonLinX) > 0 ){
          
          xpred <- .predictY2X_nonLinear(xpred, yy=Y,bb=bg[,notOther],
                                         ss=sg[notOther,notOther],
                                         priorIV = priorXIV,priorX=priorX,
                                         factorObject = factorBeta, interObject = interBeta,
                                         lox, hix)$x
        }
        
        if( length(px) > 0 ){
          wn <- which(!is.finite(xpred),arr.ind=T)
          if(length(wn) > 0){
            tmp <- matrix(priorX,Q,nrow(wn))
            xpred[wn[,1],] <- t(tmp)
          }
          xpred[,px] <- .predictY2X_linear(xpred, yy=Y, bb=bg[,notOther],
                                           ss=sg[notOther,notOther], sinv = sinv,
                                           priorIV = priorXIV, 
                                           priorX=priorX,predCols=px, 
                                           REDUCT=REDUCT, lox, hix)[,px]
          wn <- which(!is.finite(xpred),arr.ind=T)
          if(length(wn) > 0){
            tmp <- matrix(priorX,Q,nrow(wn))
            xpred[wn[,1],] <- t(tmp)
          }
        }
        
        if( length(factorBeta$linFactor) > 0 ){
          
          # predict all factors
          xtmp <- xpred
          xtmp[,factorBeta$findex] <- 
            .predictY2X_linear(xpred, yy=Y, 
                               bb=bg[,notOther],
                               ss=sg[notOther,notOther], sinv = sinv,
                               priorIV = priorXIV, 
                               priorX=priorX,predCols=factorBeta$findex, 
                               REDUCT=REDUCT, lox, hix)[,factorBeta$findex]
          for(k in 1:length(factorBeta$linFactor)){
            
            mm  <- factorBeta$linFactor[[k]]
            tmp <- xtmp[,mm]
            
            tmp[,1] <- 0
            ix  <- apply(tmp,1,which.max)   
            
            tmp <- tmp*0
            tmp[ cbind(1:n,ix) ] <- 1
            tmp <- tmp[,-1,drop=F]
            xpred[,mm[-1]] <- tmp
          }
        }
        xpred[,1] <- 1
      }
    }
    
    setTxtProgressBar(pbar,g)
    
    bgu <- bg                    # unstandardize beta
    if(length(standRows) > 0){
      if(TIME){
        bgu <- S2U%*%mub
        lambda[ gindex[,c('rowG','colW')]] <- Lmat[wL]
        lambdas <- S2UL%*%mug      # unstandardized lambda
        lgibbs[g,] <- lambdas[,notOther]
      }else{
        bgu <- S2U%*%x%*%bg
      }
    }
    
    bgibbsUn[g,] <- bgu          # unstandardized
    bgibbs[g,]   <- bg           # standardized
    
    # Fmatrix centered for factors, 
    # bg is standardized by x, bgu is unstandardized
    
    tmp <- .contrastCoeff(beta=bg[,notOther], 
                          notStand = notStandard[notStandard %in% xnames], 
                          sigma = sg[notOther,notOther], sinv = sinv,
                          stand = standMat, factorObject=factorBeta )
    agg   <- tmp$ag
    beg   <- tmp$eg
    fsens <- tmp$sens
    
    fSensGibbs[g,]  <- sqrt(diag(fsens))
    bFacGibbs[g,] <- agg       # stand for X and W, centered for factors
    
    if(TRAITS){
      Atrait <- bg%*%t(specTrait[,colnames(yp)])  # standardized
      Strait <- specTrait[,colnames(yp)]%*%sg%*%t(specTrait[,colnames(yp)])
      bTraitGibbs[g,] <- Atrait
      mgibbs[g,] <- Strait
      
      minv <- ginv(Strait)
      
      tmp <- .contrastCoeff(beta=Atrait, 
                            notStand = notStandard[notStandard %in% xnames], 
                            sigma = Strait, sinv = minv,
                            stand = standMat, factorObject=factorBeta )
      tagg   <- tmp$ag
      bTraitFacGibbs[g,] <- tagg # stand for X and W, centered for factors
    }
    
    if(TIME){
      
      tmp <- .contrastCoeff(beta=lambda[,notOther], 
                            notStand = notStandardL[notStandardL %in% xlnames], 
                            sigma = sg[notOther,notOther],sinv = sinv,
                            stand=standMatL, factorObject=factorLambda)
      lgg   <- tmp$ag
      leg   <- tmp$eg
      lsens <- tmp$sens
      
      lss <- sqrt(diag(lsens))
      
      if(g == 1){
        if( !all(names(lss) %in% colnames(gsensGibbs)) )
          colnames(gsensGibbs) <- names(lss)
      }
      
      gsensGibbs[g,names(lss)] <- lss
      
      alpha[ aindex[,c('toW','fromW')] ] <- Amat[wA]
      asens <- Amat[,notOther]%*%sinv%*%t(Amat[,notOther])
      asens <- sqrt(diag(asens))
      asensGibbs[g,] <- asens
    }
    
    if(FULL)ygibbs[g,] <- as.vector(yp)
    
    if(g > burnin){
      
      ntot   <- ntot + 1
      ypred  <- ypred + yp
      ypred2 <- ypred2 + yp^2
      
      tmp <- .dMVN(w[,notOther], muw[,notOther], sg[notOther,notOther], log=T)
      
      sumDev <- sumDev - 2*sum(tmp) 
      yerror <- yerror + (yp - y)^2
      
      fmat <- fmat + fsens
      
      sMean  <- sMean + sg
      
      wpred  <- wpred + w
      wpred2 <- wpred2 + w^2
      
      if(RICHNESS){
        
        yy <- yp
        
        if('PA' %in% typeNames){
          wpa <- which(typeNames[inRichness] == 'PA')
          yy[,inRichness[wpa]] <- round(yp[,inRichness[wpa]]) #######
        }
        
        if(length(notPA) > 0){
          w0 <- which(yy[,notPA] <= 0)
          w1 <- which(yy[,notPA] > 0)
          yy[,notPA][w0] <- 0
          yy[,notPA][w1] <- 1
        }
        
        shan <- sweep(yy[,inRichness], 1, rowSums(yy[,inRichness]), '/')
        shan[shan == 0] <- NA
        shan <- -rowSums(shan*log(shan),na.rm=T)
        shannon <- shannon + shan
        
        wpp <- which(yy > 0)
        ypredPres[wpp]  <- ypredPres[wpp] + yp[wpp]
        ypredPres2[wpp] <- ypredPres2[wpp] + yp[wpp]^2
        ypredPresN[wpp] <- ypredPresN[wpp] + 1
        
        presence[,inRichness] <- presence[,inRichness] + yy[,inRichness]
        ones <- round(rowSums(yy[,inRichness]))
        more <- round(rowSums(yy[,inRichness]*wrich[,inRichness,drop=F]))
        richFull <- .add2matrix(ones,richFull)
        richness <- .add2matrix(more,richness)  # only for non-missing
      }
      
      if(RANDOM){
        alphaRanSums <- alphaRanSums + alphaRandGroup
      }
      
      if(mmiss > 0){
        ymissPred[ymiss]  <- ymissPred[ymiss] + y[ymiss]
        ymissPred2[ymiss] <- ymissPred2[ymiss] + y[ymiss]^2
      }
      if(nmiss > 0){
        xmissSum  <- xmissSum + x[xmiss]
        xmissSum2 <- xmissSum2 + x[xmiss]^2
      }
      
      if(PREDICTX & length(predXcols) > 0){
        predx  <- predx + xpred
        predx2 <- predx2 + xpred^2
      }
      
      wa0 <- which(colSums(agg) != 0)
      ess[notOther[wa0],notOther[wa0]]  <- 
        t(agg[,wa0,drop=F])%*%covE%*%agg[,wa0,drop=F] 
      if(TIME){
        wa0 <- which(colSums(lgg) != 0)
        ess[notOther[wa0],notOther[wa0]]  <- 
          ess[notOther[wa0],notOther[wa0]] +
          t(lgg[,wa0,drop=F])%*%covL%*%lgg[,wa0,drop=F] 
      }
      
      emat[notOther[wa0],notOther[wa0]] <- 
        emat[notOther[wa0],notOther[wa0]] + 
        .cov2Cor( ess[notOther[wa0],notOther[wa0]] )
      
      lo[ ess < 0 ] <- lo[ ess < 0 ] + 1
      hi[ ess > 0 ] <- hi[ ess > 0 ] + 1
      
      ess[notOther,notOther] <- ginv(ess[notOther,notOther])
      
      lm[ ess < 0 ] <- lm[ ess < 0 ] + 1  # neg values
      hm[ ess > 0 ] <- hm[ ess > 0 ] + 1  # pos values
      
      if(REDUCT){
        rndTot <- rndTot + rndEff
      }
      
      if(TRAITS){
        yw     <- sweep(yp,1,rowSums(yp),'/')
        yw[yw <= 0]   <- 0
        yw[is.na(yw)] <- 0
        Ttrait <- .gjamPredictTraits(yw,specTrait[,colnames(yp)], traitTypes)
        tpred  <- tpred + Ttrait
        tpred2 <- tpred2 + Ttrait^2
      }
    }
  }     
  
  ################# end gibbs loop ####################
  
  
  otherpar$S <- S 
  otherpar$Q <- Q
  otherpar$snames <- snames
  otherpar$xnames <- xnames
  
  presence <- presence/ntot
  
  if(RICHNESS){
    missRows <- sort(unique(ymiss[,1]))
    richNonMiss <- richness/ntot            #only non-missing plots
    yr  <- as.matrix(ydata[,inRichness]) 
    yr[yr > 0] <- 1
    yr <- rowSums(yr,na.rm=T)
    vv  <- matrix(as.numeric(colnames(richNonMiss)),n,
                  ncol(richNonMiss),byrow=T)
    rmu <- rowSums( vv * richNonMiss )/rowSums(richNonMiss)
    
    rsd <- sqrt( rowSums( vv^2 * richNonMiss )/rowSums(richNonMiss) - rmu^2)
    
    vv  <- matrix(as.numeric(colnames(richFull)),n,ncol(richFull),byrow=T)
    rfull <- rowSums( vv * richFull )/rowSums(richFull)
    rfull[missRows] <- NA
    rmu <- rowSums(presence)
    
    shan <- sweep(y[,inRichness], 1, rowSums(y[,inRichness]), '/')
    shan[shan == 0] <- NA
    shanObs <- -rowSums(shan*log(shan),na.rm=T)
    
    richness <- cbind(yr, rmu, rsd, rfull, shanObs, shannon/ntot )
    colnames(richness) <- c('obs','predMu','predSd','predNotMissing',
                            'H_obs', 'H_pred')
    if(TIME)richness[timeZero,] <- NA
    
    ypredPresMu  <- ypredPres/ypredPresN   #predictive mean and se given presence
    ypredPresMu[ypredPresN == 0] <- 0
    yvv <- ypredPres2/ypredPresN - ypredPresMu^2
    yvv[!is.finite(yvv)] <- 0
    ypredPresSe <- sqrt(yvv)
  }
  
  if('OC' %in% typeNames){
    ordMatShift <- matrix(ordShift,n,length(ordCols),byrow=T)
    onames <- snames[ordCols]
    wb <- match(paste(onames,'intercept',sep='_'), colnames(bgibbs))
    bgibbs[,wb] <- bgibbs[,wb] + matrix(ordShift,ng,length(ordCols),byrow=T)
    bgibbsUn[,wb] <- bgibbsUn[,wb] + matrix(ordShift,ng,length(ordCols),byrow=T)
    y[,ordCols] <- y[,ordCols] + ordMatShift
  }
  
  if(mmiss > 0){
    ymissPred[ymiss]  <- ymissPred[ymiss]/ntot
    yd <- ymissPred2[ymiss]/ntot - ymissPred[ymiss]^2
    yd[!is.finite(yd)| yd < 0] <- 0
    ymissPred2[ymiss] <- sqrt(yd)
    
    if('OC' %in% typeNames){
      ymissPred[,ordCols] <- ymissPred[,ordCols] + ordMatShift
    }
  }
  
  xunstand    <- .getUnstandX(x, standRows, standMu[,1],
                              standMat[,1], interBeta$intMat)$xu
  
  rmspeBySpec <- sqrt( colSums(yerror)/ntot/n )
  rmspeAll    <- sqrt( sum(yerror)/ntot/n/S )
  
  sMean <- sMean/ntot
  
  if(TIME){
    
    xtime <- xpred*0
    xtime[,xnames] <- x
    xtime[,xlnames] <- xl
    
    xlunstand    <- .getUnstandX(xl, standRowsL, standMuL[,1],
                                 standMatL[,1], interLambda$intMat)$xu
    xtimeUn <- xtime*0
    xtimeUn[,xnames]  <- xunstand
    xtimeUn[,xlnames] <- xlunstand
    
    loL <- hiL <- lambdaMuUn <- lambdaSeUn <- lambda*0
    tmp1 <- colMeans(ggibbs[burnin:ng,])    #unstandardized
    tmp2 <- apply(ggibbs[burnin:ng,],2,sd)
    lambdaMuUn[ gindex[,c('rowG','colW')] ] <- tmp1
    lambdaSeUn[ gindex[,c('rowG','colW')] ] <- tmp2
    loL[gindex[,c('rowG','colW')] ] <- loLmat[wL]
    hiL[gindex[,c('rowG','colW')] ] <- hiLmat[wL]
    
    loA <- hiA <- alphaMu <- alphaSe <- matrix(0,S,S)
    tmp1 <- colMeans(alphaGibbs[burnin:ng,])    #unstandardized
    tmp2 <- apply(alphaGibbs[burnin:ng,],2,sd)
    alphaMu[ aindex[,c('toW','fromW')] ] <- tmp1
    alphaSe[ aindex[,c('toW','fromW')] ] <- tmp2
    loA[ aindex[,c('toW','fromW')] ] <- loAmat[wA]
    hiA[ aindex[,c('toW','fromW')] ] <- hiAmat[wA]
    
    gsensMu <- colMeans(gsensGibbs[burnin:ng,]) 
    gsensSd <- apply(gsensGibbs[burnin:ng,],2,sd)
    asensMu <- colMeans(asensGibbs[burnin:ng,])
    asensSd <- apply(asensGibbs[burnin:ng,],2,sd)
  }
  
  tmp <- .chain2tab(bgibbs[burnin:ng,], snames, xnames)
  betaStandXmu <- tmp$mu
  betaStandXTable <- tmp$tab
  
  tmp <- .chain2tab(bgibbsUn[burnin:ng,], snames, xnames)
  betaMu <- tmp$mu
  betaTable <- tmp$tab
  
  tmp <- .chain2tab(bFacGibbs[burnin:ng,], snames[notOther], rownames(agg))
  betaStandXWmu <- tmp$mu
  betaStandXWTable <- tmp$tab
  
  tmp <- .chain2tab(fSensGibbs[burnin:ng,,drop=F])
  sensTable <- tmp$tab[,1:4]
  
  yMu <- ypred/ntot
  y22 <- ypred2/ntot - yMu^2
  y22[y22 < 0] <- 0
  ySd <- sqrt(y22)
  
  cMu <- cuts
  cSe <- numeric(0)
  
  wMu <- wpred/ntot
  wpp <- pmax(0,wpred2/ntot - wMu^2)
  wSd <- sqrt(wpp)
  
  if('OC' %in% typeNames){
    yMu[,ordCols] <- yMu[,ordCols] + ordMatShift
    wMu[,ordCols] <- wMu[,ordCols] + ordMatShift
  }
  
  meanDev <- sumDev/ntot
  
  tmp <- .dMVN(wMu[,notOther],x%*%betaMu[,notOther],
               sMean[notOther,notOther], log=T)
  pd  <- meanDev - 2*sum(tmp )
  DIC <- pd + meanDev
  
  yscore <- colSums( .getScoreNorm(y[,notOther],yMu[,notOther],
                                   ySd[,notOther]^2),na.rm=T )  # gaussian w
  xscore <- xpredMu <- xpredSd <- NULL
  standX <- xmissMu <- xmissSe <- NULL
  
  if(RANDOM){
    ns <- 500
    simIndex <- sample(burnin:ng,ns,replace=T)
    tmp <- .expandSigmaChains(snames, alphaVarGibbs, otherpar, simIndex=simIndex,
                              sigErrGibbs, kgibbs, REDUCT=F)
    alphaRandGroupVarMu <- tmp$sMu
    alphaRandGroupVarSe <- tmp$sSe
    alphaRandByGroup <- alphaRanSums/ntot
    
  }
  
  if(PREDICTX){
    xpredMu <- predx/ntot
    xpredSd <- predx2/ntot - xpredMu^2
    xpredSd[xpredSd < 0] <- 0
    xpredSd <- sqrt(xpredSd)
    
    xrow <- standRows
    xmu  <- standMu[,1]
    xsd  <- standMat[,1]
    
    if(TIME){
      xrow <- c(standRows, standRowsL) 
      ww   <- !duplicated(names(xrow))
      xrow <- names(xrow)[ww]
      xmu  <- c(standMu[xrow,1], standMuL[xrow,1])
      xsd  <- c(standMat[xrow,1],standMatL[xrow,1])
      #  xrow <- names(xrow)[ww]
      #  xrow <- match(xrow,colnames(xpredMu))
      #  names(xrow) <- colnames(xpredMu)[xrow]
    }
    
    xpredMu <- .getUnstandX(xpredMu, xrow, xmu, xsd, intMat)$xu
    xpredSd[,xrow] <- xpredSd[,xrow]*matrix( xsd[xrow], n, length(xrow),
                                             byrow=T ) 
    
    if(TIME){
      if(Q == 2)xscore <- mean( .getScoreNorm(xtime[,2],
                                              xpredMu[,2],xpredSd[,2]^2) )
      if(Q > 2)xscore <- colMeans(.getScoreNorm(xtime[,-1],
                                                xpredMu[,-1],xpredSd[,-1]^2) )
    }else{
      if(Q == 2)xscore <- mean( .getScoreNorm(x[,2],
                                              xpredMu[,2],xpredSd[,2]^2) )
      if(Q > 2)xscore <- colMeans(.getScoreNorm(x[,-1],
                                                xpredMu[,-1],xpredSd[,-1]^2) )
    }
    
    if(TIME){
      wz <- wMu
      wz[wz < 0] <- 0
      Vmat[tindex[,2],] <- wz[tindex[,2], 
                              gindex[,'colW']]*xl[tindex[,2], gindex[,'colX']]
      Vmat[timeZero,]   <- wz[timeZero, 
                              gindex[,'colW']]*xl[timeZero, gindex[,'colX']]
      Umat <- wz[,uindex[,1]]*wz[,uindex[,2]] 
      
      Amat[ aindex[,c('rowA','fromW')] ] <- alphaMu[ aindex[,c('toW','fromW')] ]
      Lmat[ gindex[,c('rowL','colW')] ] <- lambdaMuUn[ gindex[,c('rowG','colW')] ]
      
      muw <- x%*%betaMu[,notOther] + Vmat%*%Lmat[,notOther] + Umat%*%Amat[,notOther]
      
      tmp <- .dMVN(wMu[,notOther],muw[,notOther],
                   sMean[notOther,notOther], log=T )
      pd  <- meanDev - 2*sum(tmp )
      DIC <- pd + meanDev
    }
  }
  
  if(nmiss > 0){
    xmissMu <- xmissSum/ntot
    xmissSe <- sqrt( xmissSum2/ntot - xmissMu^2 )
  }
  
  if(length(standRows) > 0){                #unstandardize
    standX <- cbind(standMu[,1],standMat[,1])
    colnames(standX) <- c('xmean','xsd')
    rownames(standX) <- rownames(standMat)
  }
  
  # betaSens, sigma and R
  
  ns <- 500
  simIndex <- sample(burnin:ng,ns,replace=T)
  
  tmp <- .expandSigmaChains(snames, sgibbs, otherpar, simIndex=simIndex,
                            sigErrGibbs, kgibbs, REDUCT)
  corMu <- tmp$rMu; corSe <- tmp$rSe
  sigMu  <- tmp$sMu; sigSe  <- tmp$sSe
  
  whichZero <- which(lo/ntot < ematAlpha & 
                       hi/ntot < ematAlpha,arr.ind=T) #not different from zero
  whConZero <- which(lm/ntot < ematAlpha & 
                       hm/ntot < ematAlpha,arr.ind=T)
  
  ematrix  <- emat/ntot
  fmatrix  <- fmat/ntot
  
  tMu <- tSd <- tMuOrd <- btMu <- btSe <- stMu <- stSe <- numeric(0)
  
  if(TRAITS){
    
    tMu <- tpred/ntot
    tSd <- sqrt(tpred2/ntot - tMu^2)
    wo  <- which(traitTypes == 'OC')    #predict ordinal scores
    M   <- ncol(tMu)
    
    if(length(wo) > 0){
      tMuOrd <- tMu*0
      for(j in wo)tMuOrd[,j] <- round(tMu[,j],0) - 1
      tMuOrd <- tMuOrd[,wo]
    }
    
    tmp <- .chain2tab(bTraitGibbs[burnin:ng,], tnames, xnames) #standardized
    betaTraitXMu <- tmp$mu
    betaTraitXTable <- tmp$tab
    
    tmp <- .chain2tab(mgibbs[burnin:ng,], tnames, tnames) 
    varTraitMu <- tmp$mu
    varTraitTable <- tmp$tab
    
    tmp <- .chain2tab(bTraitFacGibbs[burnin:ng,], tnames, rownames(tagg) )
    betaTraitXWmu <- tmp$mu
    betaTraitXWTable <- tmp$tab
  }
  
  if('OC' %in% typeNames){
    nk  <- length(ordCols)
    nc  <- ncut - 3
    
    os <- rep(ordShift,nc)
    
    cgibbs <- cgibbs + matrix(os,ng,length(os),byrow=T)
    
    tmp <- .processPars(cgibbs)$summary
    cMu <- matrix(tmp[,'estimate'],nk,nc)
    cSe <- matrix(tmp[,'se'],nk,ncut-3)
    cMu <- cbind(ordShift,cMu)
    cSe <- cbind(0,cSe)
    colnames(cMu) <- colnames(cSe) <- cnames[-c(1,ncut)]
    rownames(cMu) <- rownames(cSe) <- snames[ordCols]
    breakMat[ordCols,c(2:(2+(ncol(cMu))-1))] <- cMu
  }
  
  if('PA' %in% typeNames){
    zMu <- yMu
    zSd <- ySd
  }
  
  
  # outputs
  if(length(reductList) == 0)reductList <- list(N = 0, r = 0)
  reductList$otherpar <- otherpar
  
  modelList$effort    <- effort;      modelList$formula <- formula
  modelList$typeNames <- typeNames;    modelList$censor <- censor
  modelList$effort    <- effort; modelList$holdoutIndex <- holdoutIndex
  modelList$REDUCT    <- REDUCT;       modelList$TRAITS <- TRAITS
  modelList$ematAlpha <- ematAlpha; modelList$traitList <- traitList
  modelList$reductList <- reductList; modelList$ng <- ng
  modelList$burnin <- burnin
  
  inputs <- list(xdata = xdata, x = xunstand, standX = standX,
                 standMat = standMat, standRows = standRows, y = y, 
                 notOther = notOther, other = other, breakMat = breakMat, 
                 designTable = designTable, classBySpec = classBySpec, 
                 factorBeta = factorBeta, interBeta = interBeta,
                 linFactor = linFactor, intMat = intMat, RANDOM = RANDOM)
  missing <- list(xmiss = xmiss, xmissMu = xmissMu, xmissSe = xmissSe, 
                  ymiss = ymiss, ymissMu = ymissPred, ymissSe = ymissPred2)
  parameters <- list(betaMu = betaMu, betaTable = betaTable, 
                     betaStandXmu = betaStandXmu, 
                     betaStandXTable = betaStandXTable,
                     betaStandXWmu =  betaStandXWmu,
                     betaStandXWTable = betaStandXWTable,
                     corMu = corMu, corSe = corSe, 
                     sigMu = sigMu, sigSe = sigSe, 
                     ematrix = ematrix, fmatrix = fmatrix,
                     whichZero = whichZero, whConZero = whConZero,
                     wMu = wMu, wSd = wSd, sensTable = sensTable)
  prediction <- list(presence = presence, xpredMu = xpredMu, xpredSd = xpredSd,
                     ypredMu = yMu, ypredSd = ySd, richness = richness)
  chains <- list(sgibbs = sgibbs, bgibbs = bgibbs, bgibbsUn = bgibbsUn,
                 fSensGibbs = fSensGibbs, bFacGibbs = bFacGibbs) 
  fit <- list(DIC = DIC, yscore = yscore, 
              xscore = xscore, rmspeAll = rmspeAll,
              rmspeBySpec = rmspeBySpec)
  if(FULL)chains <- append(chains, list(ygibbs = ygibbs))
  if(RANDOM){
    parameters <- append(parameters,
                         list( randGroupVarMu = alphaRandGroupVarMu,
                               randGroupVarSe = alphaRandGroupVarSe,
                               randByGroup = alphaRandByGroup) )
  }
  if(RICHNESS){
    prediction <- append(prediction, 
                         list(yPresentMu = ypredPresMu, yPresentSe = ypredPresSe))
  }
  if(REDUCT) {
    parameters <- append(parameters, list(rndEff = rndTot/ntot))#, specRand = specRand))
    chains <- append(chains,list(kgibbs = kgibbs, sigErrGibbs = sigErrGibbs))
  }
  
  if('OC' %in% typeNames){
    parameters <- c(parameters,list(cutMu = cMu, cutSe = cSe))
    chains <- c(chains,list(cgibbs = cgibbs))
    modelList <- c(modelList,list(yordNames = yordNames))
  }
  
  if(TRAITS){
    parameters <- c(parameters,
                    list(betaTraitXMu = betaTraitXMu, 
                         betaTraitXTable = betaTraitXTable,
                         varTraitMu = varTraitMu, 
                         varTraitTable = varTraitTable,
                         betaTraitXWmu = betaTraitXWmu,
                         betaTraitXWTable = betaTraitXWTable))
    prediction <- c(prediction, list(tMuOrd = tMuOrd, tMu = tMu, tSe = tSd))
    chains <- append( chains,list(bTraitGibbs = bTraitGibbs,
                                  bTraitFacGibbs = bTraitFacGibbs,
                                  mgibbs = mgibbs) ) 
  }
  if(TIME){
    inputs <- c(inputs, list(xtime = xtime, timeZero = timeZero,
                             interLambda = interLambda, 
                             factorLambda = factorLambda))
    chains <- c(chains, list(ggibbs = ggibbs, alphaGibbs = alphaGibbs,
                             gsens = gsensGibbs, asens = asensGibbs))
    parameters <- c(parameters, 
                    list(lambdaMuUn = lambdaMuUn, lambdaSeUn = lambdaSeUn, 
                         lambdaLo = loL, lambdaHi = hiL,
                         alphaMu = alphaMu, alphaSe = alphaSe,
                         alphaLo = loA, alphaHi = hiA,
                         gsensMu = gsensMu, gsensSe = gsensSd, 
                         asensMu = asensMu, asensSe = asensSd,
                         aindex = aindex, wA = wA, unidex = uindex))
  }
  
  chains     <- chains[ sort( names(chains) )]
  fit        <- fit[ sort( names(fit) )]
  inputs     <- inputs[ sort( names(inputs) )]
  missing    <- missing[ sort( names(missing) )]
  modelList  <- modelList[ sort( names(modelList) )]
  parameters <- parameters[ sort( names(parameters) )]
  prediction <- prediction[ sort( names(prediction) )]
  
  all <- list(chains = chains, fit = fit, inputs = inputs, missing = missing,
              modelList = modelList, parameters = parameters,
              prediction = prediction)
  all$call <- match.call()
  all <- all[ sort(names(all)) ]
  class(all) <- "gjam"
  
  all
}


.getPars00 <- function(CLUST, x, N, r, Y, B, D, Z, sigmaerror, K, pvec,
                     alpha.DP, inSamples,...){
  
  # Y includes all terms but x%*%beta
  
  nn   <- length(inSamples)
  p    <- ncol(x)
  S    <- ncol(Y)
  ntot <- nrow(Y)
  nn   <- length(inSamples)
  
  covR <- solveRcpp( (1/sigmaerror)*crossprod(Z[K,]) + diag(r) ) # Sigma_W
  z1   <- crossprod( Z[K,]/sigmaerror,t(Y - x%*%t(B)) )
  RR   <- rmvnormRcpp(ntot, mu = rep(0,r), sigma = covR ) + t(crossprod( covR,z1))
  if(nn < ntot)RR[-inSamples,] <- rmvnormRcpp(ntot-nn,mu=rep(0,r), sigma=diag(r))
  rndEff <- RR%*%t(Z[K,])
  
  res        <- sum((Y[inSamples,] - x[inSamples,]%*%t(B) - rndEff[inSamples,] )^2)
  sigmaerror <- 1/rgamma(1,shape=(S*nn + 1)/2, rate=res/2)
  
  if(CLUST){   #only until convergence
    avec <- 1/rgamma(r, shape = (2 + r )/2,
                     rate = ((1/1000000) + 2*diag(solveRcpp(D)) ) )
    
    D    <- .riwish(df = (2 + r + N - 1), S = (crossprod(Z) + 2*2*diag(1/avec)))
    Z    <- fnZRcpp(kk=K, Yk=Y[inSamples,], Xk=x[inSamples,], Dk=D, Bk=B,
                    Wk=RR[inSamples,], sigmasqk=sigmaerror, Nz=N)
    
    pmat <- getPmatKRcpp(pveck = pvec,Yk = Y[inSamples,], Zk = Z,
                         Xk = x[inSamples,], Bk = B, Wk = RR[inSamples,],
                         sigmasqk = sigmaerror)
    K    <- unlist( apply(pmat, 1, function(x)sample(1:N, size=1, prob=x)) )
    
    #pvec <- .sampleP(N = N, avec = rep(alpha.DP/N,(N-1)),
                     #bvec = ((N-1):1)*alpha.DP/N, K = K)
    pvec <- .sampleP(N=N, avec=rep(1,(N-1)),
                                  bvec=rep(alpha.DP,(N-1)), K=K)
    # 
    # alpha.DP<-rgamma(1, shape=N+2-1, rate = 1/2-log(pvec[N]))
  }
  
  list(A = Z[K,], D = D, Z = Z, K = K, pvec = pvec,
       sigmaerror = sigmaerror, rndEff = rndEff,alpha.DP=alpha.DP)
}

.paramWrapper00 <- function(REDUCT, inSamples,SS){
  
  if(REDUCT){
    
    function(CLUST, x,beta,Y,otherpar){
      
      N  <- otherpar$N
      r  <- otherpar$r
      D  <- otherpar$D
      Z  <- otherpar$Z
      sigmaerror <- otherpar$sigmaerror
      K          <- otherpar$K
      pvec       <- otherpar$pvec
      alpha.DP   <- otherpar$alpha.DP
      tmp        <- .getPars00(CLUST, x = x, N = N, r = r, Y = Y, B = t(beta),
                             D = D, Z = Z, sigmaerror = sigmaerror,
                             K = K, pvec = pvec, alpha.DP = alpha.DP,
                             inSamples = inSamples, SELECT = F)
      
      sg <- with(tmp, .expandSigma(sigma = tmp$sigmaerror, SS, Z = tmp$Z,
                                   K = tmp$K, REDUCT=T))
      
      otherpar <- list(A = tmp$A, N = N, r = r, D = tmp$D, Z = tmp$Z,
                       sigmaerror = tmp$sigmaerror,
                       pvec = tmp$pvec, K = tmp$K, alpha.DP = alpha.DP)
      
      return(list(sg = sg, rndEff = tmp$rndEff, otherpar = otherpar))
    }
    
  } else {
    
    function(CLUST, x, beta,Y,otherpar){
      
      sigmaDf  <- otherpar$sigmaDf
      XX  <- crossprod(x[inSamples,])
      IXX <- solveRcpp(XX)
      WX  <- crossprod(x[inSamples,], Y[inSamples,])
      WIX <- IXX%*%WX
      
      sg <- .updateWishartNoPrior( x[inSamples,], Y[inSamples,], sigmaDf,
                                   beta = beta, IXX = IXX, WX = WX, WIX = WIX,
                                   TRYPRIOR = T)$sigma
      otherpar=list(Z = NA, K = NA, sigmaDf = sigmaDf)
      
      return(list(sg = sg, otherpar = otherpar))
    }
  }
}
