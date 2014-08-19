qpcrLogistic.robust.optim <- function(P,Data,isOptim) {
  # P[1] = baseline background signal, P[2] = amplitude of signal
  # P[3] = inflection, P[4] = slope, P[5] = centrality factor
  # observedData to contain [X = cycle number (numeric), Y = raw signal level (numeric)    
  # isOptim = set to TRUE if R function optim is calling this function or FALSE for full output. 
  
  # calculate the fitted values (generalized logistic function)
  fit = list(x=Data$x,y=Data$y);
  fit$yFit = P[1]+(P[2]/((1+(P[5]*exp(-P[4]*(fit$x-P[3]))))^(1/P[5])));
  
  # an R2 statistic will be constructed, the follow data will be calculated
  fit$residuals = fit$yFit-fit$y;
  fit$SSE = sum(fit$residuals^2);
  fit$SST = Data$SST
  fit$Rsq = 1-(fit$SSE/fit$SST);
  
  # robust regression approach using a modified Welch function to "cap" the effect that outliers can have
  # the tuning parameter "tuneConstant" is in units of the dependent variable and is basically the "cap" level
  fit$tuneConstant = 500;
  fit$rhoError = (fit$tuneConstant^2)*(1-exp((-2/(fit$tuneConstant^2))*(fit$residuals^2)));
  
  fit$modelFunction <- function(x,p) {
    return(p[1]+(p[2]/((1+(p[5]*exp(-p[4]*(x-p[3]))))^(1/p[5]))));            
  }
    
  fit$modeldydx <- function(x,p) {
    return(p[2]*p[4]*exp(-p[4]*(x-p[3]))*((1+(p[5]*exp(-p[4]*(x-p[3]))))^-((1+p[5])/p[5])));      
  }
  
  # in addition to robust estimate, I have implemented a weighting function which focuses the optimization on the
  fit$weights = (1+((3-1)*fit$modeldydx(fit$x,P)/fit$modeldydx(P[3],P)));
  #fit$weights = rep(1,localVars$n);
  
  # the m-estimator - weighted and robust per point
  #fit$Mestimator = sum(fit$weights*fit$SSE)/sum(fit$weights);
  fit$Mestimator = sum(fit$weights*fit$rhoError)/sum(fit$weights);
  
  if (isOptim) {
    # return M-Estimator
    return(fit$Mestimator);   
  }
  else {
    # return full fit
    fit$parameters = P;   
    fit$xPlot = seq(0,max(fit$x)+7,0.1);
    fit$yFitPlot = fit$modelFunction(fit$xPlot,fit$parameters);
    
    return(fit);
  } 
}
