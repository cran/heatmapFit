heatmap.fit<-function(form, fam, dat, reps=1000, span.l="aicc", color=F){
  
  
  # This function calculates the deviation between predictions from a GLM model and
  # non-parametrically smoothed empirical frequencies. Models that are good at prediction will tend
  # to have small deviation; for example, if a model predicts that Pr(y=1)=k%, about k% of observations
  # with this predicted probability should have y=1.
  
  dat<-subset(dat, select=all.vars(form))
  dat<-na.omit(dat)
  n<-dim(dat)[1]
  
  store<-glm(formula=form, family=fam, data=dat)                          # fit the model to the data
  pred<-predict(store, type="response")                                   # obtain p-hats from model
  YY<-eval(parse(text=paste("dat$",all.vars(form)[1], sep="")))           # call the y-variable YY
  
  dat.names<-c(names(dat),"index.counter")                                # add an index counter to the data
  dat<-data.frame(cbind(dat, 1:dim(dat)[1]))                              #
  names(dat)<-dat.names                                                   #
  
  # if an optimal bandwidth is specified, find it
  if(span.l=="aicc" | span.l=="gcv"){
    
    # use Michael Friendly's function for calculating AIC/GCV from a loess
    loess.aic <- function (x) {
      
      # Written by Michael Friendly 
      # http://tolstoy.newcastle.edu.au/R/help/05/11/15899.html
      
      if (!(inherits(x,"loess"))) stop("Error: argument must be a loess object")
      # extract values from loess object
      span <- x$pars$span
      n <- x$n
      traceL <- x$trace.hat
      sigma2 <- sum( x$residuals^2 ) / (n-1)
      delta1 <- x$one.delta
      delta2 <- x$two.delta
      enp <- x$enp
      
      aicc <- log(sigma2) + 1 + 2* (2*(traceL+1)) / (n-traceL-2)
      aicc1<- n*log(sigma2) + n* ( (delta1/delta2)*(n+enp)/(delta1^2/delta2)-2 )
      gcv  <- n*sigma2 / (n-traceL)^2
      
      result <- list(span=span, aicc=aicc, aicc1=aicc1, gcv=gcv)
      return(result)
      
    }
    
    
    # this is the optimization object; it just returns the AIC or GCV for a bandwidth argument
    smooth.err<-function(span.arg){
      assign("last.warning", NULL, envir = baseenv())
      ok<-T
      plot.model<-withCallingHandlers(tryCatch(loess(YY~pred, degree=1, span=span.arg)),  warning = function(w){ok<<-F; invokeRestart("muffleWarning")})
      if(ok==T){return(eval(parse(text=paste("loess.aic(plot.model)$", span.l, sep=""))))}
      if(ok==F){return(2e10)}
    }  
    
    # do the optimization, set the span argument to the optimal value
    span.l.name<-span.l
    span.l<-optimize(f=smooth.err, interval=c(0.01, 0.99))$minimum                          
    cat(span.l.name, "Chosen Span = ", span.l, "\n")
    
  }
  
  ok<-T
  plot.model<-withCallingHandlers(tryCatch(loess(YY~pred, degree=1, span=span.l)),  warning = function(w){ok<<-F; invokeRestart("muffleWarning")})
  # if a problem is detected with the GCV/AIC-chosen span, default to a 75% bandwidth
  if(ok==F){cat("Defaulting to span = 0.75", "\n"); span.l<-0.75; plot.model<-loess(YY~pred, degree=1, span=span.l)}
  y.obs<-predict(plot.model, newdata=pred)                             # determine observed y using loess smooth
  for(j in 1:length(y.obs)){y.obs[j]<-max(min(y.obs[j],1),0)}          # keep y.obs in bounds
  
  # prepare for heat map plot: predict empirical y at a bunch of points y-hat
  tick<-(max(pred)-min(pred))/500
  pr<-seq(from=min(pred),to=max(pred),by=tick)
  yo<-predict(plot.model,newdata=pr)
  for(j in 1:length(yo)){yo[j]<-max(min(yo[j],1),0)}   # keep yo in bounds
  
  # determine the distribution of the heat map line
  # using bootstrapping
  y.obs.boot.store<-matrix(data=NA, nrow=reps, ncol=length(pred))     # for calculating heat map at each observation
  y.obs.bs<-matrix(data=NA, nrow=reps, ncol=length(pr))   # for the heat map plot
  
  cat(c("Generating Bootstrap Predictions...","\n"))
  pb <- txtProgressBar(min = 0, max = reps, style = 3)                    # text progress bar for bootstrap replicates
  o<-order(pred)
  for(i in 1:reps){
    
    setTxtProgressBar(pb, i)
    
    boot.pred<-pred
    boot.y<-ifelse(runif(n, min=0, max=1)<boot.pred,1,0)                  # generate a bootstrapped y dataset from the bootstrapped data
    
    plot.model3<-withCallingHandlers(tryCatch(loess(boot.y~boot.pred, degree=1, span=span.l)),  warning = function(w){invokeRestart("muffleWarning")})  # calculate heat map line for each bootstrap; suppress minor warnings
    
    y.obs.boot<-predict(plot.model3, newdata=boot.pred)                        # determine observed y using loess smooth
    for(j in 1:length(y.obs.boot)){y.obs.boot[j]<-max(min(y.obs.boot[j],1),0)}    # keep y.obs in bounds
    y.obs.boot.store[i,]<-y.obs.boot[o]                              # save the bootstrapped y obs values
    
    y.obs.boot.two<-predict(plot.model3, newdata=pr)                        # determine observed y using loess smooth
    for(j in 1:length(y.obs.boot.two)){y.obs.boot.two[j]<-max(min(y.obs.boot.two[j],1),0)}    # keep y.obs in bounds
    y.obs.bs[i,]<-y.obs.boot.two                           # save the bootstrapped y obs values
    
  }
  close(pb)
  
  y.obs<-y.obs[o]
  y.obs.prob<-c()
  for(i in 1:length(pred)){
    
    y.obs.prob[i]<-((1/length(y.obs.boot.store[,i]))*(sum(y.obs[i]<=y.obs.boot.store[,i])+sum(y.obs[i]<y.obs.boot.store[,i]))/2)
    y.obs.prob[i]<-ifelse(y.obs.prob[i]<0.5,y.obs.prob[i],1-y.obs.prob[i])
    
  }
  
  y.obs.prob2<-c()
  for(i in 1:length(pr)){
    
    y.obs.prob2[i]<-((1/length(y.obs.bs[,i]))*(sum(yo[i]<=y.obs.bs[,i])+sum(yo[i]<y.obs.bs[,i]))/2)
    y.obs.prob2[i]<-ifelse(y.obs.prob2[i]<0.5,y.obs.prob2[i],1-y.obs.prob2[i])
    
  }
  
  # Construct the heat map plot
  
  def.par <- par(no.readonly = TRUE)
  o<-order(pred)
  strwidth
  nf <- layout(matrix(c(1,2),1,2,byrow=TRUE), widths=c(12,4), heights=15)
  par(oma=c(0,0,3,0))
  plot(y.obs~pred[o], type="n", ylim=c(min(c(pred,y.obs)), max(c(pred,y.obs))),ylab="Smoothed Empirical Pr(y=1)", xlab="Model Prediction, Pr(y=1)", main="Heat Map Plot")
  f<-cbind(yo,pr)
  if(color==T){
    for(i in 1:length(pr)){
      segments(f[i,2]-(tick/2),f[i,1],f[i,2]+(tick/2),f[i,1],col=rgb(red=2*255*(.5-y.obs.prob2[i]),green=0,blue=2*255*(y.obs.prob2[i]),maxColorValue = 255),lwd=5)
    }
  }else{
    for(i in 1:length(pr)){
      segments(f[i,2]-(tick/2),f[i,1],f[i,2]+(tick/2),f[i,1],col=gray((1/.6)*(y.obs.prob2[i])),lwd=5)
    }   
  }
  abline(0,1, lty=2)
  legend("topleft", lty=c(1,2), lwd=c(5,1), legend=c("heat map line","perfect fit"))
  rug(pred[o])
  par(mar=c(3,0,2,0))
  clr<-seq(.001,0.499,by=0.001)
  x.clr<-rep(5,length(clr))
  CLR<-cbind(clr,x.clr)
  if(color==T){
    plot(CLR[,1]~CLR[,2],bty="n",pch=15,xaxt="n",yaxt="n",xlab=" ",ylab=" ",main="p-Value\nLegend",xlim=c(4.9,5.1),col=rgb(red=2*255*(.5-CLR[,1]),green=0,blue=2*255*(CLR[,1]),maxColorValue = 255))
  }else{
    plot(CLR[,1]~CLR[,2],bty="n",pch=15,xaxt="n",yaxt="n",xlab=" ",ylab=" ",main="p-Value\nLegend",xlim=c(4.9,5.1),col=gray((1/.6)*(CLR[,1])))  
  } 
  axis(2,at=c(seq(.000,0.5,by=0.1)),labels=c(0.01,seq(.1,0.5,by=0.1)),line=-2,las=2)
  mtext("Predicted Probability Deviation \n Model Predictions vs. Empirical Frequency",outer=T,line=0,adj=.5)
  par(def.par)
  
  out1<-y.obs.prob
  heatmapstat<-sum(out1<=.1)/length(out1)
  cat("\n")
  cat("*******************************************", " \n", sep="")
  cat(heatmapstat*100, "% of Observations have one-tailed p-value <= 0.10", " \n", sep="")
  cat("Expected Maximum = 20%", " \n", sep="")
  cat("*******************************************", " \n", sep="")
  cat("\n")
  
  return(list(heatmap.obs.p=out1))
  
  
  
}
