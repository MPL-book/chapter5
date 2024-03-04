#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot.coxph_mpl=function(x,se="M2QM2",ask=TRUE,which=1:4,upper.quantile=.95,...){
  which.plot=rep(TRUE,4)
  if(!is.null(which)){which.plot[-which]=FALSE}
  if(sum(which.plot)==1){ask=FALSE}    
  if(ask){oask <- devAskNewPage(TRUE)
  on.exit(devAskNewPage(oask))
  }
  control = x$control
  knots   = x$knots
  pos     = x$coef$Theta<x$control$min.theta
  n.x       = 1000
  V_x_X     = seq(x$knots$Alpha[1],max(x$knots$Alpha),length=n.x)    
  colw      = terrain.colors(x$dim$m+1)
  prob      = upper.quantile
  #quant     = quantile(x$data$time,prob=prob)
  M_psi_Xm  = basis_mpl(V_x_X,knots,control$basis,control$order,which=1)    
  if(which.plot[1]){
    plot(1,1,pch="",xlim=range(V_x_X),ylim=max(M_psi_Xm)*c(-.05,1),axes=FALSE,
         xlab="Survival time",ylab=expression(psi[u]^{o}*(t)),
         main=paste(if(control$basis=="uniform"){"Uniform"}else{
           if(control$basis=="gaussian"){"Gaussian"}else{
             if(control$basis=="msplines"){"M-spline"}else{
               if(control$basis=="epanechikov"){"Epanechikov"}}}},
           " bases used to approximate the baseline hazard\n",
           "(",x$dim$m," bases)",sep=""))
    abline(v=knots$Alpha,col=gray(.9))
    #abline(v=quant,col=gray(.5),lty=2)
    for(u in 1:x$dim$m){lines(V_x_X,M_psi_Xm[,u],col=colw[u],lty=2*pos[u]+1)}
    axis(2,las=2)
    axis(1)
    if(any(pos)){
      legend("topright",legend=c("Observed","Censored","Knots",as.expression(bquote(widehat(F)^-1*(.(prob)))),
                                 expression(hat(theta)[u]==0),expression(hat(theta)[u]>0)),
             pch=c(1,1,NA,NA,NA,NA),lty=c(0,0,1,2,3,1),col=c("black","red",gray(.9),gray(.5),colw[1],colw[1]),ncol=1,bty="n",cex=.75)
    }else{
      legend("topright",legend=c("Observed","Censored","Knots",as.expression(bquote(widehat(F)^-1*(.(prob))))),
             pch=c(1,1,NA,NA),lty=c(0,0,1,2),col=c("black","red",gray(.9),gray(.5)),ncol=1,bty="n",cex=.75)            
    }
  }
  if(any(which.plot[2:4])){
    cov_ThetaTheta = x$covar[[se]][-c(1:x$dim$p),-c(1:x$dim$p)]    
    xlim = range(x$knots$Alpha)
    plot_bh <- function(j,V_x_X,Theta,covar,control,knots,pos,prob,xlim,...){
      M_Ppsi_Xm   = basis_mpl(V_x_X,knots,control$basis,control$order,which=as.numeric(j>1)+1)
      V_sd2.Hh0.X = diag(M_Ppsi_Xm[,]%*%covar[,]%*%t(M_Ppsi_Xm[,]))
      pos.var     = V_sd2.Hh0.X>0
      V_Hh0_X     = c(M_Ppsi_Xm%*%matrix(Theta,ncol=1))[pos.var]
      V_x_X       = V_x_X[pos.var]
      V_sd.Hh0.X  = sqrt(V_sd2.Hh0.X[pos.var])
      upper       = V_Hh0_X+2*V_sd.Hh0.X
      lower       = V_Hh0_X-2*V_sd.Hh0.X
      lower[lower<0] = 0
      if(j==3){
        V_Hh0_X = exp(-V_Hh0_X)
        lower   = exp(-lower)
        upper   = exp(-upper)
      }
      plot(1,1,pch="",xlim=xlim,axes=FALSE,
           ylim=c(ifelse(j<3,0,min(upper[V_x_X<xlim[2]])),ifelse(j<3,max(upper[V_x_X<xlim[2]]),max(lower[V_x_X<xlim[2]]))),
           main=paste("Estimate of the",c(" baseline hazard"," cumulative baseline hazard"," baseline survival")[j]," function",sep=""),
           xlab="Survival time",ylab=c(expression(h[0]*(t)),expression(H[0]*(t)),expression(S[0]*(t)))[j])
      #rect(quant,0,max(V_x_X),max(upper)*1.5,col=gray(.9),border = NA)
      #abline(v=quant,col=gray(.5),lty=2)            
      axis(2,pos=0,las=2)
      axis(1,pos=0)
      xx = c(V_x_X,V_x_X[length(V_x_X):1])
      yy = c(upper,lower[length(V_x_X):1])
      confcol = paste(substr(colw[length(Theta)],1,7),70,sep="")
      polygon(xx,yy,col=confcol,border = "gray")
      lines(V_x_X,V_Hh0_X,lwd=1.1,col=colw[1])
      legend(ifelse(j<3,"topleft","topright"),legend=c("Estimate","95% conf. interval",
                                                       as.expression(bquote(widehat(F)^-1*(.(prob))))),pch=c(NA,15,NA),lty=c(1,0,3),
             col=c(colw[1],confcol,gray(.5)),ncol=1,bty="n",cex=.75)        
    }
    if(which.plot[2]){plot_bh(1,V_x_X,x$coef$Theta,cov_ThetaTheta,control,knots,pos,prob,xlim)}
    if(which.plot[3]){plot_bh(2,V_x_X,x$coef$Theta,cov_ThetaTheta,control,knots,pos,prob,xlim,...)}    
    if(which.plot[4]){plot_bh(3,V_x_X,x$coef$Theta,cov_ThetaTheta,control,knots,pos,prob,xlim,...)}                
  }
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print.summary.coxph_mpl=function(x,se="M2QM2",...) {
  inf = x$inf
  cat("\n")
  print(inf$call)
  cat("\n-----\n\n")
  cat("Cox Proportional Hazards Model Fit Using MPL","\n\n\n")
  cat("Penalized log-likelihood  :  ",inf$ploglik,"\n",sep="")
  cat(ifelse(inf$control$max.iter[1]==1,
             "Fixed smoothing value     :  ",
             "Estimated smoothing value :  "),
      inf$control$smooth,"\n",sep="")
  cat("Convergence               : ",
      if(inf$control$max.iter[1]==1){
        ifelse(inf$iter[2]<inf$control$max.iter[2],
               paste(c("Yes (",inf$iter[2]," iter.)"),collapse=""),"NO")        
      }else{
        ifelse(inf$iter[1]<inf$control$max.iter[1]&inf$iter[2]<inf$control$max.iter[3],
               paste(c("Yes (",inf$iter[1]," + ",inf$iter[2]," iter.)"),collapse=""),"NO")
      },"\n\n")
  cat("Data             : ",inf$data,"\n",sep="")    
  number=format(c(inf$dim$n,inf$dim$n.obs,inf$dim$n-inf$dim$n.obs))
  percentage=format(c(0,inf$dim$n.obs,inf$dim$n-inf$dim$n.obs)/inf$dim$n*100)
  cat("Number of obs.   : ",number[1],"\n",sep="")
  cat("Number of events : ",number[2]," (",percentage[2],"%)\n",sep="")
  cat("Number of cens.  : ",number[3]," (",percentage[3],"%)\n\n",sep="")
  cat("Regression parameters : ",deparse(inf$call[[2]]),"\n",sep="")
  printCoefmat(x$Beta, P.values=TRUE, has.Pvalue=TRUE,...)    
  cat("\nBaseline hasard parameters approximated using",
      if(inf$control$basis=="uniform"){"a step function"}else{
        if(inf$control$basis=="gaussian"){"Gaussian splines"}else{
          if(inf$control$basis=="msplines"){"M-splines"}else{
            "Epanechikov splines"}}},":\n")
  if(inf$control$basis=="uniform"){
    cat(paste(" (",inf$dim$m," equal events bins)\n",sep=""))
  }else{
    cat(paste(" (2 (min/max) + ",inf$control$n.knots[1]," quantile knots + ",inf$control$n.knots[2]," equally spaced knots",
              if(inf$control$basis=="msplines"|inf$control$basis=="epanechikov"){
                paste(" +",inf$control$order,"(order) - 2")}," = ",inf$dim$m," parameters)",sep=""),"\n")}
  if(inf$full){
    printCoefmat(x$Theta, P.values=TRUE, has.Pvalue=TRUE,...)    
  }else{print(x$Theta,...)}
  cat("\n-----\n\n")
}




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print.coxph_mpl=function(x, ...) {
  cat("\n")
  print(x$call)
  cat("\nLog-likelihood : ",x$ploglik[1],"\n",sep="")    
  cat("\nRegression parameters :\n")
  vect=c(x$coef$Beta)
  names(vect)=dimnames(x$data$X)[[2]]
  print(vect, ...)
  cat("\nBaseline hasard parameters : \n")
  print(x$coef$Theta, ...)
  cat("\n")
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
summary.coxph_mpl=function(object,se="M2QM2",full=FALSE,...) {
  col.names = c("Estimate", "Std. Error", "z-value", "Pr(>|z|)")
  seB   = object$se$Beta[[se]]
  matxB = cbind(object$coef$Beta,seB,object$coef$Beta/seB,2*(1-pnorm(abs(object$coef$Beta/seB))))
  dimnames(matxB)=list(colnames(object$data$X),col.names)
  pos = object$coef$Theta>object$control$min.theta
  if(full){
    seT   = object$se$Theta[[se]]
    matxT = cbind(object$coef$Theta,seT,object$coef$Theta/seT,2*(1-pnorm(abs(object$coef$Theta/seT))))
    dimnames(matxT)=list(format(1:object$dim$m,just="right"),col.names)
    matxT = matxT[pos,]
  }else{    matxT = object$coef$Theta[pos]
  names(matxT) = seq(1,object$dim$m)[pos]}
  out = list(Beta = matxB, Theta = matxT, inf = list(iter = object$iter, call=object$call, data = object$data$name,
                                                     control = object$control, dim = object$dim, trunc = trunc, full = full, ploglik = object$ploglik[1]))
  class(out) = "summary.coxph_mpl"
  out
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
coef.summary.coxph_mpl=function(object, parameters = "Beta", ...) {
  object[[parameters]]
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
coef.coxph_mpl=function(object, parameters = "Beta", ...) {
  out = object$coef[[parameters]]
  if(parameters == "Beta") names(out) = colnames(object$data$X)
  else names(out) = 1:object$dim$m
  out
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
residuals.coxph_mpl=function(object,...) {
  control  = object$control
  out = as.data.frame(matrix(NA,object$dim$n,5))
  colnames(out)  = c("time1", "time2","censoring","coxsnell","martingale")
  out$time1      = object$data$time[,1L]
  out$time2      = object$data$time[,2L]
  out$censoring  = object$data$censoring
  i.r   = which(object$data$censoring==0)
  i.obs = which(object$data$censoring==1)
  i.l   = which(object$data$censoring==2)
  i.ic  = which(object$data$censoring==3)
  M_Psi_Xm1 = basis_mpl(out$time1,object$knots,control$basis,control$order,which=2)
  S1 = exp(-exp(object$data$X%*%object$coef$Beta)*M_Psi_Xm1%*%object$coef$Theta)
  r = rep(NA, object$dim$n)
  if(length(i.r)>0) r[i.r]   = 1 - log(S1[i.r])
  if(length(i.obs)>0) r[i.obs] = -log(S1[i.obs])
  if(length(i.l)>0) r[i.l]   = (1-S1[i.l]*(1-log(S1[i.l])))/(1-S1[i.l])
  if(length(i.ic)>0) {
    M_Psi_Xm2 = basis_mpl(out$time2[i.ic],object$knots,control$basis,control$order,which=2)
    S2 = exp(-exp(object$data$X[i.ic,,drop=F]%*%object$coef$Beta)*M_Psi_Xm2%*%object$coef$Theta)
    r[i.ic]  = (S1[i.ic]*(1-log(S1[i.ic])) - S2*(1-log(S2)))/(S1[i.ic] - S2)
  }

  out$coxsnell   = r 
  out$martingale = 1 - out$coxsnell
  class(out) =c("residuals.coxph_mpl","data.frame")
  out
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot.residuals.coxph_mpl=function(x,ask=TRUE,which=1:2,upper.quantile=.95,...){
  prob = upper.quantile
  which.plot=rep(TRUE,2)
  if(!is.null(which)){which.plot[-which]=FALSE}
  if(sum(which.plot)==1){ask=FALSE}
  if(ask){oask <- devAskNewPage(TRUE)
  on.exit(devAskNewPage(oask))
  }
  # martingale
  if(which.plot[1]){
    plot(1:nrow(x),x$martingale,col=(!(x$censoring==1))+1,ylab="",
         xlab="Index",main="Martingale Residuals",axes=FALSE,
         ylim=c(min(x$martingale),1+(1-min(x$martingale))*.075),...)
    mtext(expression(delta[i]-plain(e)^(x[i]^T*hat(beta))*widehat(H)[0](t[i])),2,padj=-2,...)
    abline(h=1,col="light gray")
    abline(h=0,col="light gray",lty=2)
    axis(1,...)
    axis(2,...)
    legend("top",ncol=2,legend=c("observed","censored"),col=c(1,2),pch=1,cex=.75)
  }
  # coxsnell
  if(which.plot[2]){
    plot(1:nrow(x),x$coxsnell,col=(!(x$censoring==1))+1,ylab="",
         xlab="Index",main="Adjusted Cox & Snell Residuals",axes=FALSE,
         ylim=c(0,max(x$coxsnell)*1.075),...)
    mtext(expression(plain(e)^(x[i]^T*hat(beta))*widehat(H)[0](t[i])),2,padj=-2,...)
    abline(h=1,col="light gray",lty=2)
    axis(1,pos=0,...)
    axis(2,...)
    legend("top",ncol=2,legend=c("observed","censored"),col=c(1,2),pch=1,cex=.75)
  }
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
predict.coxph_mpl=function(object,se="M2QM2",type="risk",i=NULL,time=NULL,upper.quantile=.95,...) {
  prob  = upper.quantile
  covar = object$covar[[se]]
  Beta  = object$coef$Beta
  Theta = object$coef$Theta
  p     = object$dim$p
  m     = object$dim$m
  if(length(i)>1){warning("only the first observation will be considered\n",call. = FALSE)}
  # time
  if(is.null(time)){
    n.x   = 1000
    V_x_X = seq(object$knots$Alpha[1],max(object$knots$Alpha),length=n.x)
  }else{	
    n.x   = length(time)
    V_x_X = time
  }
  # x
  xTB = if(is.null(i)){apply(object$data$X,2,mean)}else{object$data$X[i[1],,drop=FALSE]}
  Mu  = c(exp(xTB%*%Beta))
  # risk
  out = data.frame(time = V_x_X, mid = NA, se = NA, low = NA, high = NA)
  if(type=="risk"){
    M_psi_Xm = basis_mpl(V_x_X,object$knots,object$control$basis,object$control$order,which=1)
    out$mid  = Mu*M_psi_Xm%*%Theta
    # correction factor
    M_corr_mpm = matrix(c(rep(M_psi_Xm%*%Theta*Mu,p)*rep(Beta,each=n.x),Mu*M_psi_Xm),ncol=p+m)
    out$se     = sqrt(diag(M_corr_mpm%*%object$covar[[se]]%*%t(M_corr_mpm)))
    out$low    = out$mid - 2*out$se; out$low[out$low<0] = 0
    out$high   = out$mid + 2*out$se
    # survival
  }else{
    M_Psi_Xm = basis_mpl(V_x_X,object$knots,object$control$basis,object$control$order,which=2)
    out$mid  = exp(-Mu*M_Psi_Xm%*%Theta)
    # correction factor
    M_corr_mpm = rep(-out$mid*Mu,p+m)*matrix(c(rep(M_Psi_Xm%*%Theta,p)*rep(Beta,each=n.x),M_Psi_Xm),ncol=p+m)
    out$se     = sqrt(diag(M_corr_mpm%*%object$covar[[se]]%*%t(M_corr_mpm)))
    out$low    = out$mid - 2*out$se; out$low[out$low<0] = 0
    out$high   = out$mid + 2*out$se; out$high[out$high>1] = 1
  }
  # out
  times=c(object$data$time[,1L],object$data$time[which(object$data$censoring==3),2L])
  attributes(out)$inf = list(i=i[1], upper.quantile=prob, upper.value = quantile(times,prob),
                             max = max(times), user.time = !is.null(time), m=m, risk = type=="risk")
  colnames(out)[2]=type
  class(out) =c("predict.coxph_mpl","data.frame")
  out
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot.predict.coxph_mpl=function(x,...){
  inf     = attr(x,"inf")
  colw    = terrain.colors(3)[1:2]
  pos.var = x$se>0
  main    = paste(ifelse(inf$risk,"Predicted instantaneous risk at time t",
                         "Predicted probability of survival after time t"),"\n",
                  ifelse(is.null(inf$i),"for 'average' covariates",paste("for observation",inf$i)),sep="")
  plot(1,1,pch="",axes=FALSE,xlab="Time",ylim=c(0,max(x$high)),main=main,
       xlim=if(inf$user.time){c(0.5,nrow(x)+.5)}else{c(min(x$time),inf$upper.value)},ylab="",...)
  axis(2,...)
  abline(h=0,col="light gray")
  mtext(if(inf$risk){expression(hat(h)(t[i]))}else{expression(widehat(S)(t[i]))},2,padj=-2,...)
  if(inf$user.time){
    axis(1,at=1:nrow(x),labels=x$time,tick=TRUE,pos=0,...)
    for(tw in 1:nrow(x)){
      if(x$se[tw]>0){
        arrows(tw, y0=x$low[tw], y1=x$high[tw], angle=90, code=3, col=colw[2],...)
      }
    }
    points(1:nrow(x),x[,2],col=colw[1],...)
    legend("topleft",legend=c("Estimate","95% conf. interval"),pch=c(NA,15),lty=c(1,0),
           col=c(colw[1],confcol),ncol=1,bty="n",cex=.75)
  }else{
    axis(1,pos=0,...)
    rect(inf$upper.value,0,inf$max,max(x$high)*1.5,col=gray(.9),border = NA)
    if(inf$upper.quantile<1)abline(v=inf$upper.value,col=gray(.5),lty=2)
    xx = c(x$time[pos.var],x$time[pos.var][length(x$time[pos.var]):1])
    yy = c(x$high[pos.var],x$low[pos.var][length(x$time[pos.var]):1])
    confcol = paste(substr(colw[2],1,7),70,sep="")
    polygon(xx,yy,col=confcol,border = "gray",...)
    lines(x$time,x[,2],lwd=1.1,col=colw[1],...)
    legend(ifelse(inf$risk,"topleft","topright"),legend=c("Estimate","95% conf. interval",
                                                          as.expression(bquote(widehat(F)^-1*(.(inf$upper.quantile))))),pch=c(NA,15,NA),lty=c(1,0,3),
           col=c(colw[1],confcol,gray(.5)),ncol=1,bty="n",cex=.75)
  }
}



