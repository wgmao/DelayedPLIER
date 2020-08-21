library(HDF5Array)
library(glmnet)

rowNorm=function(x){
  s=DelayedMatrixStats::rowSds(x)
  m=DelayedMatrixStats::rowMeans2(x)
  
  x=sweep(x,1,m,"-")
  x=sweep(x,1,s,"/")
  return(x)
}#rowNorm


commonRows=function(data1, data2){
  intersect(rownames(data1), rownames(data2))
}#commonRows




num.pc = function (data, method="elbow", B = 20, seed = NULL) 
{
  
  method=match.arg(method, c("elbow", "permutation"))
  if (!is.null(seed)) {
    set.seed(seed)
  }
  warn <- NULL
  if((class(data)!="list")&(class(data)!="rsvd")){
    message("Computing svd")
    n <- ncol(data)
    m <- nrow(data)
    data=rowNorm(data)
    if(n<500){
      k=n
    }
    else{
      k=max(200,n/4)
    }
    if(k==n){
      uu <- svd(data)
    }
    else{
      set.seed(123456);uu <- rsvd(data,k, q=3)
    }
  }
  else if (!is.null(data[["d"]])){
    if(method=="permutation"){
      message("Original data is needed for permutation method.\nSetting method to elbow")
      method="elbow"
    }
    
    uu=data
    
  }
  
  if(method=="permutation"){
    #nn = min(c(n, m))
    dstat <- uu$d[1:k]^2/sum(uu$d[1:k]^2)
    dstat0 <- matrix(0, nrow = B, ncol = k)
    for (i in 1:B) {
      dat0 <- t(apply(data, 1, sample, replace = FALSE))
      if(k==n){
        uu0 <- svd(dat0)
      }
      else{
        set.seed(123456);
        uu0 <- rsvd(dat0,k, q=3)
      }
      dstat0[i, ] <- uu0$d[1:k]^2/sum(uu0$d[1:k]^2)
    }
    psv <- rep(1, k)
    for (i in 1:k) {
      psv[i] <- mean(dstat0[, i] >= dstat[i])
    }
    for (i in 2:k) {
      psv[i] <- max(psv[(i - 1)], psv[i])
    }
    
    nsv <- sum(psv <= 0.1)
  }
  else if (method=="elbow"){
    x=smooth(xraw<-abs(diff(diff(uu$d))), twiceit = T)
    #plot(x)
    
    
    nsv=which(x<=quantile(x, 0.5))[1]+1
    
  }
  return(nsv)
}#num.pc




pinv.ridge=function (m, alpha = 0) 
{
  msvd = svd(m)
  if (length(msvd$d) == 0) {
    return(array(0, dim(m)[2:1]))
  }
  else {
    if (alpha > 0) {
      ss = (msvd$d^2) + alpha^2
      msvd$d = ss/msvd$d
    }
    out = msvd$v %*% (1/msvd$d * t(msvd$u))
    rownames(out) = rownames(m)
    colnames(out) = colnames(m)
    out
  }
}#pinv.ridge




solveU=function(Z,  Chat, priorMat, penalty.factor,pathwaySelection="fast", glm_alpha=0.9, maxPath=10, target.frac=0.7, L3=NULL){
  
  
  Ur=Chat%*%Z #get U by OLS
  #Ur=apply(-Ur,2,rank) #rank
  Ur = t(DelayedMatrixStats::colRanks(-Ur, ties.method = "first"))
  Urm= apply(Ur,1,min)
  
  U=matrix(0,nrow=ncol(priorMat), ncol=ncol(Z))
  if(is.null(L3)){
    
    lambdas=exp(seq(-4,-12,-0.125))
    results=list()
    lMat=matrix(nrow=length(lambdas), ncol=ncol(Z))
    for(i in 1:ncol(Z)){
      if(pathwaySelection=="fast"){
        iip=which(Ur[,i]<=maxPath)
      }else{
        iip=which(Urm<=maxPath)
      }#else
      gres=glmnet(y=Z[,i], x=priorMat[,iip], penalty.factor = penalty.factor[iip], alpha=glm_alpha, lower.limits=0, lambda = lambdas,intercept=T,  standardize=F )
      #plot(gres)
      gres$iip=iip
      lMat[,i]=colSums(as.matrix(gres$beta)>0)
      results[[i]]=gres
    }#for i
    fracs=rowMeans(lMat>0)
    iibest=which.min(abs(target.frac-fracs))
    iibest
    
    
    for(i in 1:ncol(Z)){
      U[results[[i]]$iip,i]=results[[i]]$beta[,iibest]
    }#for i
    rownames(U)=colnames(priorMat)
    colnames(U)=1:ncol(Z)
    
    Utmp=solveU(Z,  Chat, priorMat, penalty.factor,pathwaySelection="fast", glm_alpha=0.9, maxPath=10,  L3=lambdas[iibest])
    
    #stop()
    return(list(U=U, L3=lambdas[iibest]))  
  }
  else{ #do one fit with a given lambda
    for(i in 1:ncol(Z)){
      if(pathwaySelection=="fast"){
        iip=which(Ur[,i]<=maxPath)
      }else{
        iip=which(Urm<=maxPath)
      }#else
      gres=glmnet(y=Z[,i], x=priorMat[,iip], penalty.factor = penalty.factor[iip], alpha=glm_alpha, lower.limits=0, lambda = L3,intercept=T,  standardize=F )
      U[iip,i]=as.numeric(gres$beta)
    }
    
    return(U)
  }
}#solveU


############################################################################
copyMat=function(mat, zero=F){
  if (class(mat)=="matrix"){
    matnew=matrix(nrow=nrow(mat), ncol=ncol(mat))
    rownames(matnew)=rownames(mat)
    colnames(matnew)=colnames(mat)
    if(zero)
      matnew[]=0
    matnew
    
  }else if (class(mat)=="DelayedMatrix"){
    matnew <- mat
    if (zero){
      matnew[] <- 0
    }else{
      matnew[] <- NA
    }#else
    return(matnew)
  }#else if
}#copyMat


BH= function(pval){p.adjust(pval, method="BH")}


crossVal=function(plierRes, data, priorMat, priorMatcv){
  
  out=matrix(ncol=4, nrow=0)
  ii=which(colSums(plierRes$U)>0)
  Uauc=copyMat(plierRes$U,T)
  Up=copyMat(plierRes$U,T)
  Up[]=1
  for ( i in ii){
    
    iipath=which(plierRes$U[,i]>0)
    
    if (length(iipath) > 1){
      for(j in iipath){
        iiheldout=which((rowSums(priorMat[,iipath, drop=F])==0) |(priorMat[,j]>0&priorMatcv[,j]==0))
        aucres=AUC(priorMat[iiheldout,j], plierRes$Z[iiheldout,i])
        out=rbind(out,c(colnames(priorMat)[j], i, aucres$auc, aucres$pval))
        Uauc[j,i]=aucres$auc
        Up[j,i]=aucres$pval
      }}else{
        j <- iipath
        iiheldout=which((rowSums(matrix(priorMat[,iipath],ncol=1))==0) |(priorMat[,j]>0&priorMatcv[,j]==0))
        aucres=AUC(priorMat[iiheldout,j], plierRes$Z[iiheldout,i])
        out=rbind(out,c(colnames(priorMat)[j], i, aucres$auc, aucres$pval))
        Uauc[j,i]=aucres$auc
        Up[j,i]=aucres$pval
      }#else
  }
  out=data.frame(out,stringsAsFactors = F)
  out[,3]=as.numeric(out[,3])
  out[,4]=as.numeric(out[,4])
  out[,5]=BH(out[,4])
  colnames(out)=c("pathway", "LV index", "AUC", "p-value", "FDR") 
  return(list(Uauc=Uauc, Upval=Up, summary=out))
}#crossVal


getAUC=function(plierRes, data, priorMat){
  Y=data
  B=plierRes$B
  Z=plierRes$Z
  k=ncol(Z)
  L1=plierRes$L1
  L2=plierRes$L2
  
  Zcv=copyMat(Z)
  
  for (i in 1:5){
    ii=(0:(floor(nrow(data)/5)-1))*5+i
    ii=ii[ii<=nrow(Z)]
    
    #Bcv=solve(crossprod(Z[-ii,])+L2*diag(k))%*%t(Z[-ii,])%*%Y[-ii,]
    Bcv=solve(t(Z[-ii,])%*%Z[-ii,]+L2*DelayedArray(diag(k)))%*%t(Z[-ii,])%*%Y[-ii,]
    #Zcv[ii,]=Y[ii, ]%*%t(Bcv)%*%solve(tcrossprod(Bcv)+L1*diag(k))
    Zcv[ii,]=Y[ii, ]%*%t(Bcv)%*%solve(Bcv %*% t(Bcv)+L1*DelayedArray(diag(k)))
  }#for i
  
  out=matrix(ncol=4, nrow=0)
  ii=which(colSums(plierRes$U)>0)
  Uauc=copyMat(plierRes$U,T)
  Up=copyMat(plierRes$U,T)
  Up[]=1;
  
  for ( i in ii){
    iipath=which(plierRes$U[,i]>0)
    
    for(j in iipath){
      aucres=AUC(priorMat[,j], Zcv[,i])
      out=rbind(out,c(colnames(priorMat)[j], i, aucres$auc, aucres$pval))
      Uauc[j,i]=aucres$auc
      Up[j,i]=aucres$pval
    }#for j
  }#for i
  
  out=data.frame(out,stringsAsFactors = F)
  out[,3]=as.numeric(out[,3])
  out[,4]=as.numeric(out[,4])
  out[,5]=BH(out[,4])
  colnames(out)=c("pathway", "LV index", "AUC", "p-value", "FDR") 
  
  return(list(Uauc=Uauc, Upval=Up, summary=out))
}#getAUC



AUC<-function(labels, values){
  posii=which(labels>0)
  negii=which(labels<=0)
  posn=length(posii)
  negn=length(negii)
  posval=values[posii]
  negval=values[negii]
  myres=list()
  if(posn>0&negn>0){
    res=wilcox.test(posval, negval, alternative="greater", conf.int=TRUE);
    
    myres$low=res$conf.int[1]
    myres$high=res$conf.int[2]
    myres$auc=(res$statistic)/(posn*negn)
    myres$pval=res$p.value
  }
  else{
    myres$auc=0.5
    myres$pval=NA
  }
  return(myres)
}#AUC


nameB=function(plierRes, top=1, fdr.cutoff=0.01, use=c("coef", "AUC")){
  use=match.arg(use, c("coef", "AUC"))
  names=vector("character",ncol(plierRes$U))
  if(use=="coef"){
    Uuse=plierRes$U
  }
  else{
    Uuse=plierRes$Uauc
  }
  if(!is.null(plierRes[["Up"]])){
    pval.cutoff=max(plierRes$summary[plierRes$summary[,5]<fdr.cutoff,4])
    
    Uuse[plierRes$Up>pval.cutoff]=0
    
  }
  else{
    warning("No p-values in PLIER object: using coefficients only")
  }
  mm=apply(Uuse,2,max)
  for(i in 1:ncol(plierRes$U)){
    if(mm[i]>0){
      names[i]=paste(i,names(sort(Uuse[,i],T))[1:top], sep=",")
    }
    else if(max(plierRes$U[,i])>0){
      names[i]=paste(i,names(sort(plierRes$U[,i],T))[1:top], sep=",")
    }
    else{
      names[i]=paste("LV",i)
    }
  }
  
  names
}#nameB



############################################################################
PLIER=function(data, priorMat,svdres=NULL, k=NULL, L1=NULL, L2=NULL, L3=NULL,  frac=0.7,  max.iter=350, trace=F, scale=T, Chat=NULL, maxPath=10, doCrossval=T, penalty.factor=rep(1,ncol(priorMat)), glm_alpha=0.9, minGenes=10, tol=1e-6, seed=123456, allGenes=F, rseed=NULL, pathwaySelection=c("complete", "fast"), output_path = "output/"){
  
  pathwaySelection=match.arg(pathwaySelection, c("complete", "fast"))
  
  if(scale){
    Y=rowNorm(data)
  }
  else{
    Y=data
  }
  
  if(nrow(priorMat)!=nrow(data) || !all(rownames(priorMat)==rownames(data))){
    if(!allGenes){
      cm=commonRows(data, priorMat)
      message(paste("Selecting common genes:", length(cm)))
      priorMat=priorMat[cm,]
      Y=Y[cm,]
    }
    else{
      extra.genes=setdiff(rownames(data), rownames(priorMat))
      eMat=matrix(0, nrow=length(extra.genes), ncol=ncol(priorMat))
      rownames(eMat)=extra.genes
      priorMat=rbind(priorMat, eMat)
      priorMat=priorMat[rownames(data),]
    }
    
  }
  numGenes=colSums(priorMat)
  
  heldOutGenes=list()
  iibad=which(numGenes<minGenes)
  priorMat[, iibad]=0
  message(paste("Removing", length(iibad), "pathways with too few genes"))
  if(doCrossval){
    
    
    priorMatCV=priorMat
    if(!is.null(seed))
      set.seed(seed)
    for(j in 1:ncol(priorMatCV)){
      
      iipos=which(priorMatCV[,j]>0)
      iiposs=sample(iipos, length(iipos)/5)
      priorMatCV[iiposs,j]=0
      heldOutGenes[[colnames(priorMat)[j]]]=rownames(priorMat)[iiposs]
      
    }#for j
    C = priorMatCV
  }else{
    C=priorMat
  }#else
  
  nc=ncol(priorMat)
  ng=nrow(data)
  ns=ncol(data)
  
  Bdiff=-1
  BdiffTrace=double()
  BdiffCount=0
  if(is.null(Chat)){
    Cp=crossprod(C)
    Chat=pinv.ridge(crossprod(C), 5)%*%(t(C))
  }
  
  Yseq = Y^2
  YsqSum=sum(DelayedMatrixStats::rowSums2(Yseq))
  #compute svd and use that as the starting point
  
  if(!is.null(svdres) && nrow(svdres$v)!=ncol(Y)){
    message("SVD V has the wrong number of columns")
    svdres=NULL
  }
  
  if(is.null(svdres)){
    message("Computing SVD")
    if(ns>500){
      message("Using rsvd")
      set.seed(123456);
      #svdres=rsvd(Y, k=min(ns, max(200, ns/4)), q=3)
      svdres = BiocSingular::runRandomSVD(Y, k = min(ns, max(200, ns/4)), center = F, scale = F)
    }
    else{
      svdres=BiocSingular::runRandomSVD(Y, k = min(ng, ns))
    }
    message("Done")
  }
  if(is.null(k)){
    k=num.pc(svdres)*2
    k <- min(k, floor(ncol(Y)*0.9))
    message("k is set to ", k)
  }
  
  
  if(is.null(L2)){
    show(svdres$d[k])
    L2=svdres$d[k]
    print(paste0("L2 is set to ",L2))
  }
  
  if(is.null(L1)){
    L1=L2/2
    print(paste0("L1 is set to ",L1))
  }
  
  
  B=t(svdres$v[1:ncol(Y), 1:k]%*%diag(svdres$d[1:k]))
  Z=(Y%*%t(B))%*%solve(tcrossprod(B)+L1*diag(k))
  Z[Z<0]=0
  if(!is.null(rseed)){
    message("using random start")
    set.seed(rseed)
    B=t(apply(B, 1, sample))
    Z=apply(Z,2,sample)
  }
  
  
  
  
  
  B <- DelayedArray(B)
  U=matrix(0,nrow=ncol(C), ncol=k)
  
  
  round2=function(x){signif(x,4)}
  message(paste0("errorY (SVD based:best possible) = ", round2(mean((Y-Z%*%B)^2))))
  
  
  iter.full.start=iter.full=20
  
  curfrac=0
  nposlast=Inf
  npos=-Inf
  if(!is.null(L3)){
    L3.given=T
  }else{
    L3.given=F
  }#else
  
  
##############################################################
for ( i in 1:max.iter){

    if(i>=iter.full.start){
    
      if(i==iter.full & !L3.given){ #update L3 to the target fraction
        Ulist=solveU(Z, Chat, C, penalty.factor, pathwaySelection, glm_alpha, maxPath, target.frac = frac)
        U=Ulist$U
        L3=Ulist$L3
        message(paste("New L3 is", L3))
        iter.full=iter.full+iter.full.start
      }else{
        #HERE
        #solveU=function(Z,  Chat, priorMat, penalty.factor,pathwaySelection="fast", glm_alpha=0.9, maxPath=10, target.frac=0.7, L3=NULL)
        
        U=solveU(Z, Chat, C, penalty.factor, pathwaySelection, glm_alpha, maxPath, L3=L3)
      }#else
        curfrac=(npos<-sum(apply(U,2,max)>0))/k
        Z1=Y%*%t(B)
        Z2=DelayedArray(L1*C%*%U)
        ratio=median((Z2/Z1)[Z2>0&Z1>0])
        #Z=(Z1+Z2)%*%solve(tcrossprod(B)+L1*diag(k))
        Z=(Z1+Z2)%*%solve(B%*% t(B)+L1*DelayedArray(diag(k)))
    }else{
      #Z=(Y%*%t(B))%*%solve(tcrossprod(B)+L1*diag(k))
      Z=(Y%*%t(B))%*%solve(B%*% t(B)+L1*DelayedArray(diag(k)))
    }#else
    
    
    Z[Z<0]=0
    oldB=B
    B=solve(t(Z)%*%Z+DelayedArray(L2*diag(k)))%*%t(Z)%*%Y
    
    
    num <- sum(DelayedMatrixStats::rowSums2((B-DelayedArray(oldB))^2))
    den <- sum(DelayedMatrixStats::rowSums2(B^2))
    Bdiff = num/den
    BdiffTrace=c(BdiffTrace, Bdiff)
    
    
    err0=sum((Y-Z%*%B)^2)+sum((Z-DelayedArray(C%*%U))^2)*L1+sum(B^2)*L2
    if(trace & i >=iter.full.start){
      
      message(paste0("iter",i, " errorY= ",erry<-round2(mean((Y-Z%*%B)^2)), ", prior information ratio= ", round(ratio,2), ", Bdiff= ",round2(Bdiff), ", Bkappa=", round2(kappa(B))), ";pos. col. U=", sum(colSums(U)>0))
    }else if (trace){
      message(paste0("iter",i, " errorY= ",erry<-round2(mean((Y-Z%*%B)^2)), ", Bdiff= ",round2(Bdiff), ", Bkappa=", round2(kappa(B))))
    }#else if
    
    if(i>52&&Bdiff>BdiffTrace[i-50]){
      BdiffCount=BdiffCount+1
      message("Bdiff is not decreasing")
    }else if(BdiffCount>1){
      BdiffCount=BdiffCount-1
    }#else if
    
    if(Bdiff<tol){
      message(paste0("converged at  iteration ", i))
      break
    }#if
    if( BdiffCount>5){
      message(paste0("converged at  iteration ", i, " Bdiff is not decreasing"))
      break
    }#if
}#for i in 1:max.iter
##############################################################
  
  
  rownames(U)=colnames(priorMat)
  colnames(U)=rownames(B)=paste0("LV", 1:k)
  
  out=list(residual=(Y-Z%*%B), B=B, Z=Z, U=U, C=C, L1=L1, L2=L2, L3=L3, heldOutGenes=heldOutGenes)
  
  if(doCrossval){
    outAUC=crossVal(out, Y, priorMat, priorMatCV)
  }else{
    message("Not using cross-validation. AUCs and p-values may be over-optimistic")
    outAUC=getAUC(out, Y, priorMat)
  }#else
  
  out$withPrior=which(colSums(out$U)>0)
  out$Uauc=outAUC$Uauc
  out$Up=outAUC$Upval
  out$summary=outAUC$summary
  tt=apply(out$Uauc,2,max)
  message(paste("There are", sum(tt>0.70), " LVs with AUC>0.70"))
  
  rownames(out$B)=nameB(out)
  
  #residual, B, Z DelayedArray
  writeHDF5Array(out$residual, filepath = paste0(output_path, "residual.hdf5"), name = "count")
  writeHDF5Array(out$B, filepath = paste0(output_path, "B.hdf5"), name = "count")
  writeHDF5Array(out$Z, filepath = paste0(output_path, "Z.hdf5"), name = "count")
  
  out
}#PLIER
