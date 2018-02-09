
run_spline_model=function(xobs,yobs,xcens,ycens){

  ### Initialize
  nobs=length(xobs)
  xgrid=seq(min(xobs)-1,max(xobs)+1,length=1000)
  parms=initialize_parms_spline(xobs,yobs,xcens,ycens,xgrid)
  xtrue=parms$xtrue
  coefs=parms$coefs
  ytrue=as.numeric(parms$ytrue)
  lowept=parms$lowept
  upperept=parms$upperept
  knotseq=parms$knotseq
  bases=parms$bases
  designMatrixGrid=parms$designMatrixGrid

  ### Run Model
  parms=bayesian_mon_errors_spline(xobs,yobs,xcens,ycens,coefs,xtrue,ytrue,xgrid,lowept,upperept,knotseq,bases,designMatrixGrid)
  MICDens=parms$MICDens; fitMat=parms$fitMat

  return(list(MICDens=MICDens,fitMat=fitMat))
}

