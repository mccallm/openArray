processQPCR <- function(dat) {

  # modeling constraints
  Cstr = list()
  Cstr$U = matrix(NA,6,5)
  Cstr$C = matrix(NA,6,1)
  Cstr$U[1,] = c(0,0,0,-1,log(2)); Cstr$C[1,1] = 0 # max per cycle change is 2x
  Cstr$U[2,] = c(0,0,0,1,0); Cstr$C[2,1] = 0.001 # slope parameter must be > 0 
  Cstr$U[3,] = c(0,0,0,0,1); Cstr$C[3,1] = 0.001 # acentrality parameter must be > 0
  Cstr$U[4,] = c(0,0,1,0,0); Cstr$C[4,1] = 0 # calculated expression index is between 0 and 60 cycles
  Cstr$U[5,] = c(0,0,-1,0,0); Cstr$C[5,1] = -60
  Cstr$U[6,] = c(0,1,0,0,0); Cstr$C[6,1] = 0.001 # amplitude of pcr reaction must be > 0

  # For all ChipID x ChipWell combinations
  nodes = makeForkCluster(nnodes=6)
  m = clusterApplyLB(cl=nodes,x=dat,fun=modelQPCR,Cstr=Cstr)
  stopCluster(nodes)
  
  u = data.frame()
  for (j in 1:length(dat)) {
      u[j,"Chip.Id"] = dat[[j]]$Chip.Id[1]
      u[j,"Chip.Well"] = dat[[j]]$Chip.Well[1]
      u[j,"Sample.Id"] = dat[[j]]$Sample.Id[1]
      u[j,"Feature.Set"] = dat[[j]]$Feature.Set[1]
      u[j,"Feature.Id"] = dat[[j]]$Feature.Id[1]
      u[j,"p1"] = m[[j]]$modelFit$parameters[1]
      u[j,"p2"] = m[[j]]$modelFit$parameters[2]
      u[j,"p3"] = m[[j]]$modelFit$parameters[3]
      u[j,"p4"] = m[[j]]$modelFit$parameters[4]
      u[j,"p5"] = m[[j]]$modelFit$parameters[5]
      u[j,"r2"] = m[[j]]$modelFit$Rsq
  }
  
  ## tranform to standard miRNA x sample matrix format
  cols <- u$Sample.Id
  tmp <- split(u,cols)
  ## order each col by the Feature.Id and Feature.Set
  tmp <- lapply(tmp, function(x) x[order(paste(x$Feature.Set,x$Feature.Id,sep="::")),])
  ## merge cols to create data matrices
  e <- qc <- matrix(nrow=length(tmp[[1]]$Feature.Id),ncol=length(tmp))
  for(k in 1:ncol(e)){
    e[,k] <- tmp[[k]]$p3
    qc[,k] <- tmp[[k]]$r2
  }
  rownames(e) <- rownames(qc) <- paste(tmp[[1]]$Feature.Set,tmp[[1]]$Feature.Id,sep="::")
  colnames(e) <- colnames(qc) <- names(tmp)
  
  ## we may need to change this later to flag controls
  ft <- rep("Target",nrow(e))  
  
  fl <- matrix("Passed",nrow=nrow(e),ncol=ncol(e))
  fl[which(qc<0.99,arr.ind=TRUE)] <- "Flagged"
  colnames(fl) <- colnames(e)
  rownames(fl) <- rownames(e)

  fc <- matrix("OK",nrow=nrow(e),ncol=ncol(e))
  fc[which(is.na(e),arr.ind=TRUE)] <- "Undetermined"
  colnames(fc) <- colnames(e)
  rownames(fc) <- rownames(e)
  
  obj <- new("qPCRset", exprs=e, quality=qc, flag=fl)
  featureNames(obj) <- rownames(e)
  featureType(obj) <- ft
  featureCategory(obj) <- as.data.frame(fc)
  
  tab <- data.frame(sampleName=names(tmp), 
                    chipID=unlist(lapply(tmp,function(x) x$Chip.Id[1]))
  )
  phenoData(obj) <- AnnotatedDataFrame(data=tab)
  
  return(obj)
}
