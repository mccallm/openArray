processQPCR <- function(filename,fileFormat="default") {

  # read in raw data table
  if (fileFormat=="default") {
    d = read.table(filename, header=TRUE, dec=".", sep=",", comment.char="");
  } else if (fileFormat=="LifeTech") {
    d = read.table(filename, skip=15, header=TRUE, dec=".", sep=",", comment.char="");
    names(d)[names(d)=="Barcode"] = "Chip.Id";
    names(d)[names(d)=="Well"] = "Chip.Well";
    d[,c("Sample.Id","Feature.Set")] = read.table(text=as.character(d$Sample.Name),sep='_');
    names(d)[names(d)=="Target.Name"] = "Feature.Id";
    names(d)[names(d)=="Cycle.Number"] = "Cycle";
    names(d)[names(d)=="Rn"] = "Value";
  } else {
    stop("Error: Input file format not recognized.");
  }

  d = d[,c("Chip.Id","Chip.Well","Sample.Id","Feature.Set","Feature.Id","Cycle","Value")];    
  
  #Some basic error check
  if(is.null(d$Value)){
    stop("Error: You must have a Value column for the raw fluorescence values.");
  }

  # find unique data series records  
  u = unique(d[,c("Chip.Id","Chip.Well","Sample.Id","Feature.Id","Feature.Set")])
  
  # modeling constraints
  Cstr = list();
  Cstr$U = matrix(NA,6,5);
  Cstr$C = matrix(NA,6,1);
  Cstr$U[1,] = c(0,0,0,-1,log(2)); Cstr$C[1,1] = 0; # max per cycle change is 2x
  Cstr$U[2,] = c(0,0,0,1,0); Cstr$C[2,1] = 0.001; # slope parameter must be > 0 
  Cstr$U[3,] = c(0,0,0,0,1); Cstr$C[3,1] = 0.001; # acentrality parameter must be > 0
  Cstr$U[4,] = c(0,0,1,0,0); Cstr$C[4,1] = 0; # calculated expression index is between 0 and 60 cycles
  Cstr$U[5,] = c(0,0,-1,0,0); Cstr$C[5,1] = -60;
  Cstr$U[6,] = c(0,1,0,0,0); Cstr$C[6,1] = 0.001; # amplitude of pcr reaction must be > 0

  # For all ChipID x ChipWell combinations
  bSize = 100;
  nRecords = nrow(u);
  nodes = makeForkCluster(nnodes=3);
  for (i in seq(1,nRecords,bSize)) {
    b = list();
    print(i/nRecords)
    for (j in 1:min(bSize,nRecords-i+1)) {
      b[[j]] = list();
      idx = which((d$Chip.Id==u[i+j-1,]$Chip.Id)&(d$Chip.Well==u[i+j-1,]$Chip.Well)&(d$Sample.Id==u[i+j-1,]$Sample.Id)&(d$Feature.Id==u[i+j-1,]$Feature.Id)&(d$Feature.Set==u[i+j-1,]$Feature.Set));
      b[[j]]$x = d$Cycle[idx];
      b[[j]]$y = d$Value[idx];
      b[[j]]$SST = sum((b[[j]]$y-mean(b[[j]]$y))^2);
      yRange = as.numeric(quantile(b[[j]]$y,c(0.1,0.9)));
      b[[j]]$init = c(yRange[1],yRange[2]-yRange[1],which.min(abs(mean(yRange)-b[[j]]$y)),0.5,1);
    }
    m = clusterApplyLB(cl=nodes,x=b,fun=modelQPCR,Cstr=Cstr);
    for (j in 1:min(bSize,nRecords-i+1)) {
      u[i+j-1,"p1"] = m[[j]]$modelFit$parameters[1];
      u[i+j-1,"p2"] = m[[j]]$modelFit$parameters[2];
      u[i+j-1,"p3"] = m[[j]]$modelFit$parameters[3];
      u[i+j-1,"p4"] = m[[j]]$modelFit$parameters[4];
      u[i+j-1,"p5"] = m[[j]]$modelFit$parameters[5];
      u[i+j-1,"r2"] = m[[j]]$modelFit$Rsq;
    }
  }
  stopCluster(nodes);
  
  return(u);
}
