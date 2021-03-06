fitSplineMon=function(xobs,yobs,coefs,x_est,y_mu,xgrid,lowept,upperept,knotseq,bases,designMatrixGrid,
                                    xsig=1,numIter=22000,burnin=14000,thin=4){

  numIter=10000
  burnin=2000
  fitMat=matrix(nrow=numIter-burnin,ncol=length(xgrid))
  coefMat=matrix(nrow=numIter,ncol=length(coefs))
  sigma2_save=rep(NA,numIter)
  smoothParam_save=rep(NA,numIter)
  
  nobs=length(xobs)

  smoothParam=2
  
  begin=Sys.time()


  for(iter in 1:numIter){

    if(iter %% 500==0){
      cat('iter: ',iter,'\n')
      cat(Sys.time()-being,'\n')
    }
    
    ### Update Coefs
    parms=updateCoefs(smoothParam,coefs,iter,bases,y_mu,knotseq,xobs,yobs,sigma2,coefMat)
    coefs=parms$coefs; y_mu=as.numeric(parms$y_mu)


    ### Update smoothness parameter
    parms=updateSmoothParm(smoothParam,coefs)
    smoothParam=parms$smooth
    
    ### Update error variance
    sigma2=updateVariance(xobs,yobs,y_mu,sigma2,smoothParam,coefs)

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






