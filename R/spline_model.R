fitSplineMon=function(xobs,yobs,coefs,x_est,y_est,xgrid,lowept,upperept,knotseq,bases,designMatrixGrid,
                                    xsig=1,numIter=22000,burnin=14000,thin=4){


  fitMat=matrix(nrow=numIter-burnin,ncol=length(xgrid))
  coefMat=matrix(nrow=numIter,ncol=length(coefs))
  sigma2_save=rep(NA,numIter)
  smoothParam_save=rep(NA,numIter)
  
  nobs=length(xobs)

  smoothParam=1


  for(iter in 1:numIter){


    ### Update x_est (and corresponding y_est)
    parms=updateX(x_est,y_est,knotseq,bases,lowept,upperept,coefs,xsig,ysig,xobs,yobs,nobs)
    x_est=parms$x_est; y_est=as.numeric(parms$y_est)

    
    ### Update Coefs
    parms=updateSplineCoefs(smoothParam,coefs,iter,bases,x_est,
                            y_est,knotseq,xobs,yobs,xcens,ycens,xsig,ysig,coefMat)
    coefs=parms$coefs; y_est=as.numeric(parms$y_est)


    ### Update smoothness parameter
    parms=updateSmoothParm(smoothParam,coefs)
    smoothParam=parms$smooth
    
    ### Update error variance
    sigma2=(xobs,yobs,x_est,y_est,xsig,sigma2,smoothParam,coefs)

    # smoothParam_sav[iter]=smoothParam
    coefMat[iter,]=log(coefs)
    sigma2_save[iter]=sigma2
    smoothParam_save[iter]=smoothParam
    if(iter>burnin){
      fitMat[iter-burnin,]=coefs%*%t(designMatrixGrid)
    }

  }

  ### Thin
  MICDens=MICDens[seq(1,nrow(MICDens),by=thin),]

  return(list(MICDens=MICDens,fitMat=fitMat))

}





