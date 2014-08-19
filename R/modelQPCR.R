modelQPCR <- function(Data,Cstr,isPlot=FALSE) {
  # Data is a list with $x as the independent var, $y as the dependent var
  # $SST is sum of squares total (y), $init are the initial paramter estimates  
  
  localVars = list(n=length(Data$y),error=NULL,optimResults=NULL,modelType='5-parameter logistic',Data=Data,convCheck=FALSE,modelFit=NULL);
  localVars$optimResults = constrOptim(Data$init,qpcrLogistic.robust.optim,grad=NULL,control=list(trace=0,maxit=10000000,fnscale=1),ui=Cstr$U,ci=Cstr$C,Data=Data[c("x","y","SST")],isOptim=TRUE);
  if (class(localVars$optimResults)=='list') { if ((localVars$optimResults$convergence==0)|(localVars$optimResults$convergence==10)) { localVars$convCheck = TRUE; } }
  if (localVars$convCheck) {
    localVars$modelFit = qpcrLogistic.robust.optim(localVars$optimResults$par,Data,FALSE);
  } else {
    localVars$modelFit = list('parameters'=c(NA,NA,NA,NA,NA),'Rsq'=0);
  }    
  
  fitModel = localVars;
  
  if (isPlot) {
    plot(Data$x,Data$y,pch='o'); 
    lines(fitModel$modelFit$xPlot,fitModel$modelFit$yFitPlot);
  }
  
  return(fitModel);  
}

