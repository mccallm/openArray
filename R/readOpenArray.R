readOpenArray <- function(filename, fileFormat="default") {

  # read in raw data table
  if (fileFormat=="default") {
    d <- read.table(filename, header=TRUE, dec=".", sep=",", comment.char="")
  } else if (fileFormat=="LifeTech") {
    d <- read.table(filename, skip=15, header=TRUE, dec=".", sep=",", comment.char="")
    names(d)[names(d)=="Barcode"] <- "Chip.Id"
    names(d)[names(d)=="Well"] <- "Chip.Well"
    d[,c("Sample.Id","Feature.Set")] <- read.table(text=as.character(d$Sample.Name),sep='_')
    names(d)[names(d)=="Target.Name"] <- "Feature.Id"
    names(d)[names(d)=="Cycle.Number"] <- "Cycle"
    names(d)[names(d)=="Rn"] <- "Value"
  } else {
    stop("Error: Input file format not recognized.")
  }

  d = d[,c("Chip.Id","Chip.Well","Sample.Id","Feature.Set","Feature.Id","Cycle","Value")];    
##  d = d[order(d$Chip.Id,d$Chip.Well,d$Sample.Id,d$Feature.Set,d$Feature.Id),];

  #Some basic error check
  if(is.null(d$Value)){
    stop("Error: You must have a Value column for the raw fluorescence values.");
  }
  
  tmp <- split(d,paste(d$Chip.Id,d$Chip.Well,d$Sample.Id,d$Feature.Set,d$Feature.Id,sep="::"))
  dat <- lapply(tmp, function(x){
    obj <- list()
    ## we don't really need the following, it's all in the names of the list
    ## but for now keep it for error checking
    obj$Chip.Id <- as.character(x$Chip.Id[1])
    obj$Chip.Well <- as.character(x$Chip.Well[1])
    obj$Sample.Id <- as.character(x$Sample.Id[1])
    obj$Feature.Set <- as.character(x$Feature.Set[1])
    obj$Feature.Id <- as.character(x$Feature.Id[1])
    ## the following we will always need
    obj$x <- x$Cycle
    obj$y <- x$Value
  	obj$SST = sum((obj$y-mean(obj$y))^2)
	  yRange = as.numeric(quantile(obj$y,c(0.1,0.9)))
	  obj$init = c(yRange[1],yRange[2]-yRange[1],which.min(abs(mean(yRange)-obj$y)),0.5,1)
    return(obj)
  })

  return(dat)
}
