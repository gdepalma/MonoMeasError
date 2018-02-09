spline_diagnostic=function(xgrid,xobs,yobs,xcens,ycens,MICDens,fitMat,acceptCoef,xtrue_sav,coefMat,
                           smoothAccept,smoothParam_sav){


  xobs1=xobs
  yobs1=yobs
  xobs[xcens==1 & xobs==max(xobs)]=max(xobs)+1
  xobs[xcens==-1 & xobs==min(xobs)]=min(xobs)-1
  yobs[ycens==1 & yobs==max(yobs)]=max(yobs)+1
  yobs[ycens==-1 & yobs==min(yobs)]=min(yobs)-1
  a1=data.frame(table(xobs,yobs))
  a1$xobs=as.numeric(as.character(a1$xobs))
  a1$yobs=as.numeric(as.character(a1$yobs))
  a1=a1[a1$Freq>0,]



  ### MIC Density
  densDat=data_frame(xgrid,y=apply(MICDens,2,mean))
  densDat$lower=apply(MICDens,2,function(x) quantile(x,probs=c(.05)))
  densDat$upper=apply(MICDens,2,function(x) quantile(x,probs=c(.95)))


  pltDens=ggplot(densDat,aes(x=xgrid,y))+geom_line(color='darkblue')+
    geom_ribbon(aes(ymin=lower, ymax=upper),alpha=.5,fill='cyan4')+
    scale_x_continuous(breaks = seq(min(xobs1)-1,max(xobs1)+1,by=1),
                       labels = c(paste("<",min(xobs1),sep=''),seq(min(xobs1),max(xobs1),by=1), paste(">",max(xobs1),sep='')),
                       limits = c(min(xobs1)-1,max(xobs1)+1))+
    theme_dbets()+
    labs(title='',y='',x=expression(MIC~(log["2"]~ug/mL)))

  ### MIC/DIA Relationship
  gxDat=data_frame(xgrid,y=apply(fitMat,2,mean))
  gxDat$lower=apply(fitMat,2,function(x) quantile(x,probs=c(.05)))
  gxDat$upper=apply(fitMat,2,function(x) quantile(x,probs=c(.95)))


  pltRel=ggplot(data=a1,aes(x=xobs,y=yobs,label=Freq))+geom_text(size=3.2,color='black')+
    geom_line(data=gxDat,aes(x=xgrid,y=y),color='darkblue',inherit.aes = FALSE)+
    geom_ribbon(data=gxDat,aes(x=xgrid,ymin=lower, ymax=upper),fill='cyan4',alpha=.5,inherit.aes = FALSE)+
    scale_x_continuous(breaks = seq(min(xobs1)-1,max(xobs1)+1,by=1),
                       labels = c(paste("<",min(xobs1),sep=''),seq(min(xobs1),max(xobs1),by=1), paste(">",max(xobs1),sep='')))+
    scale_y_continuous(breaks = seq(min(yobs1)-1,max(yobs1)+1,by=1),
                       labels = c(paste("<",min(yobs1),sep=''),seq(min(yobs1),max(yobs1),by=1), paste(">",max(yobs1),sep='')))+
    theme_dbets()+
    labs(title='Spline Model',y='DIA (mm)',x="")+
    coord_cartesian(ylim =c(min(yobs)-1,max(yobs)+1))


  plt1 <- ggplot_gtable(ggplot_build(pltRel))
  plt2 <- ggplot_gtable(ggplot_build(pltDens))
  maxWidth = unit.pmax(plt1$widths[2:3], plt2$widths[2:3])
  plt1$widths[2:3] <- maxWidth
  plt2$widths[2:3] <- maxWidth
  plot(grid.arrange(plt1, plt2, ncol=1, heights=c(5,2)))

  par(mfrow=c(2,1))
  xtrue_mean=apply(xtrue_sav,2,mean)
  hist(xtrue_mean,freq=FALSE)
  hist(xobs,freq=FALSE)

  cat("Mean Coefficients Acceptance: ",mean(acceptCoef),"\n")
  cat("Mean Smooth Acceptance: ",mean(smoothAccept),"\n")



  par(mfrow=c(1,1))

  tmp=gather(data.frame(coefMat),val)
  tmp$idx=rep(1:nrow(coefMat),ncol(coefMat))

  plt=ggplot(tmp,aes(x=idx,y=value))+geom_line()+facet_wrap(~val,scales='free_y')
  plot(plt)

  tmp=data.frame(idx=1:length(smoothParam_sav),value=smoothParam_sav)
  plt=ggplot(tmp,aes(x=idx,y=value))+geom_line()+labs(title='Smooth Parameter')
  plot(plt)


  return(invisible())

}
