

getLLK=function(xobs,yobs,x_est,y_est,xsig,sigma2){
  
  xllk <- pnorm(xobs,x_est,xsig,log=TRUE)
  yllk <- pnorm(yobs,y_est,sqrt(sigma2),log=TRUE)
  
  llike <- xllk+yllk
  
  return(llike)
}

PriorLLK=function(coefs,smoothParam,sigma2){
  
  priorSmooth=dunif(smoothParam,0,2,log=T)
  
  lpriorcoef = dlnorm(coefs[length(coefs)],1,100,log=T)
  for(i in (length(coefs)-1):1)
    lpriorcoef=lpriorcoef+dnorm(log(coefs[i]),coefs[i+1],smoothParam,log=T)
  
  lpriorerr=log(1/dgamma(sigma2,.01,.01))
  
  return(priorSmooth+lpriorcoef+lpriorerr)
}



findInitialSpline=function(x_est,bases,knotseq,yobs,designMatrix){

  min=999999
  for(i in seq(.5,7,by=.1)){
    coefs=seq(.1,i,length=ncol(designMatrix))
    y_est=coefs%*%t(designMatrix)
    if(sum(abs(yobs-y_est))<min){ min=sum(abs(yobs-y_est)); save=i}
  }
  coef=seq(.1,save,length=ncol(designMatrix))
  return(coef)
}


Ispline=function(intknots,lowept,upperept){
  k=3
  #determine knot sequence
  knotseq=c(rep(lowept,k+1),intknots,rep(upperept,k+1))
  numbases=length(knotseq)-k-2

  #create matrix of bases
  bases=matrix(NA,nrow=numbases,ncol=2)
  for(i in 1:numbases) bases[i,]=c(knotseq[i+1],knotseq[i+k+1])

  return(list(bases=bases,knotseq=knotseq))
}


getIsplineC=function(x_est,knotseq,bases){

  numBases=nrow(bases)
  lx_est=length(x_est)

  mat=rep(0,numBases*lx_est)

  storage.mode(mat) <- "double"
  storage.mode(knotseq) <- "double"
  storage.mode(x_est) <- "double"
  storage.mode(numBases) <- "integer"
  storage.mode(lx_est) <- "integer"
  temp=.C("getIspline",x_est,lx_est,knotseq,mat,numBases)
  designMatrix=matrix(temp[[4]],ncol=numBases,nrow=lx_est)
  return(designMatrix)

}


updateX=function(x_est,y_est,knotseq,bases,lowept,upperept,
                           coefs,xsig,sigma2,xobs,yobs,nobs){

    px_est = x_est + rnorm(nobs,0,.5)
    px_est[px_est<lowept] = lowept+.01
    px_est[px_est>upperept] = upperept-.01
    designMatrix=getIsplineC(px_est,knotseq,bases)
    py_est=as.numeric(coefs%*%t(designMatrix))
    llike=getLLK(xobs,yobs,px_est,y_est,xsig,sigma2)
    llike[is.na(llike)]=-99
    mhrat = runif(nobs)
    cond1=log(mhrat)<llike
    x_est[cond1]=px_est[cond1]
    designMatrix=getIsplineC(x_est,knotseq,bases)
    y_est=coefs%*%t(designMatrix)

  return(list(x_est=x_est,y_est=y_est))

}


updateCoefs=function(smoothParam,coefs,iter,bases,x_est,y_est,knotseq,xobs,yobs,xsig,sigma2,coefMat){

  if(iter<1000){
    Sigma=diag(.001,nrow=ncol(coefMat),ncol=ncol(coefMat))
  }else{
    Sigma=cov(coefMat[(iter-900):(iter-1),])
  }

  propCoef=tryCatch(as.numeric(exp(rmvnorm(n=1,as.numeric(log(coefs)),sigma=Sigma))),
          error=function(e){as.numeric(exp(rmvnorm(n=1,as.numeric(log(coefs)),
          sigma=diag(.001,nrow=nrow(bases),ncol=nrow(bases)))))})

  accept=0
  designMatrix=getIsplineC(x_est,knotseq,bases)
  newY=propCoef%*%t(designMatrix)

  #log likelihood and priors
  oldLLK=sum(getLLK(xobs,yobs,x_est,y_est,xsig,sigma2))
  oldPRIORLLK=sum(getPRIORLLK_spline(coefs,smoothParam,sigma2))
  newLLK=sum(getLLK(xobs,yobs,x_est,newY,xcens,ycens,xsig,sigma2))
  newPRIORLLK=sum(getPRIORLLK_spline(propCoef,smoothParam,sigma2)

  #accept/reject
  A=newLLK+newPRIORLLK-oldLLK-oldPRIORLLK
  if(is.nan(A)==FALSE & log(runif(1)) < A){
    y_est=newY
    coefs=propCoef
    accept=1
  }

  return(list(acceptCoef=accept,coefs=coefs,y_est=y_est))

}

updateSmoothParm=function(smoothParam,coefs){

  accept=0
  propSmooth=rnorm(1,smoothParam,.2)
  if(propSmooth<0 | propSmooth>2)
    return(list(smooth=smoothParam,accept=accept))
  oldPRIORLLK=sum(getPRIORLLK_spline(coefs,smoothParam,sigma2))
  newPRIORLLK=sum(getPRIORLLK_spline(coefs,propSmooth,sigma2))
  #accept/reject
  if(log(runif(1)) < newPRIORLLK-oldPRIORLLK){
    smoothParam=propSmooth
    accept=1
  }

  return(list(smooth=smoothParam,accept=accept))
}

updateVariance=function(xobs,yobs,x_est,y_est,xsig,sigma2,smoothParam,coefs){
  
  accept=0
  propSigma2=rnorm(1,sigma2,.2)
  if(propSigma2<0)
    return(sigma2)
  #log likelihood and priors
  oldLLK=sum(getLLK(xobs,yobs,x_est,y_est,xsig,sigma2))
  oldPRIORLLK=sum(getPRIORLLK_spline(coefs,smoothParam,sigma2))
  newLLK=sum(getLLK(xobs,yobs,x_est,y_est,xcens,ycens,xsig,propSigma2))
  newPRIORLLK=sum(getPRIORLLK_spline(propCoef,smoothParam,sigma2)
                  
  #accept/reject
  A=newLLK+newPRIORLLK-oldLLK-oldPRIORLLK
  if(is.nan(A)==FALSE & log(runif(1)) < A){
    sigma2=propSigma2
  }
  
  return(sigma2)
}
