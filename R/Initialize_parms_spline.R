initialize_parms_spline=function(xobs,yobs,xcens,ycens,xgrid){

  ## Initialize xtrue
  xtrue=rep(NA,length(xobs))

  n=sum(xcens==0)
  n1=sum(xcens==-1)
  n2=sum(xcens==1)

  xtrue[xcens==0]=xobs[xcens==0]-runif(n,0,1)
  xtrue[xcens==-1]=xobs[xcens==-1]-runif(n1,.5,2)
  xtrue[xcens==1]=xobs[xcens==1]+runif(n2,.5,1.5)


  ### Initalize Isplines
  lowept=min(xobs)-.5
  upperept=max(xobs)+.5
  if(sum(xcens==1)>0 | sum(ycens==1)>0)
    upperept=max(xobs)+3
  if(sum(xcens==-1)>0 | sum(ycens==-1)>0)
    lowept=min(xobs)-3
  dist=1
  parms=Ispline(seq(min(xobs),max(xobs),by=dist),lowept,upperept)
  bases=parms$bases;
  knotseq=parms$knotseq

  #### Initial spline coefficients
  designMatrix=getIsplineC(xtrue,knotseq,bases)
  coefs=findInitialSpline(xtrue,bases,knotseq,yobs,designMatrix)
  ytrue=coefs%*%t(designMatrix)

  designMatrixGrid=getIsplineC(xgrid,knotseq,bases)



  return(list(coefs=coefs,xtrue=xtrue,bases=bases,knotseq=knotseq,ytrue=ytrue,
              lowept=lowept,upperept=upperept,designMatrixGrid=designMatrixGrid))
}
