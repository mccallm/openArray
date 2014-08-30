modelQPCR <- function(splitFrame,Cstr,isPlot=FALSE) {
  # splitFrame is a data.frame from split of full table into unique cycle vs value series
  
  modelVars = list(n=length(splitFrame$Value),error=NULL,optimResults=NULL,modelType='5-parameter logistic',convCheck=FALSE,modelFit=NULL);
  modelVars$SST = sum((splitFrame$Value-mean(splitFrame$Value))^2);
  modelVars$yRange = as.numeric(quantile(splitFrame$Value,c(0.1,0.9)))
  modelVars$init = c(modelVars$yRange[1],modelVars$yRange[2]-modelVars$yRange[1],which.min(abs(mean(modelVars$yRange)-splitFrame$Value)),0.5,1)
  
  modelVars$optimResults = constrOptim(modelVars$init,qpcrLogistic.robust.optim,grad=NULL,control=list(trace=0,maxit=10000000,fnscale=1),ui=Cstr$U,ci=Cstr$C,Data=list(x=splitFrame$Cycle,y=splitFrame$Value,SST=modelVars$SST),isOptim=TRUE);
  
  if (class(modelVars$optimResults)=='list') {
  	if ((modelVars$optimResults$convergence==0)|(modelVars$optimResults$convergence==10)) {
  		modelVars$modelFit = qpcrLogistic.robust.optim(modelVars$optimResults$par,Data=list(x=splitFrame$Cycle,y=splitFrame$Value,SST=modelVars$SST),FALSE);
  	} else {
  		modelVars$modelFit = list('parameters'=c(NA,NA,NA,NA,NA),'Rsq'=0);
  	}
  }

  if (isPlot) {
    plot(splitFrame$Cycle,splitFrame$Value,pch='o'); 
    lines(modelVars$modelFit$xPlot,modelVars$modelFit$yFitPlot);
  }
  
  return(modelVars);  
}
