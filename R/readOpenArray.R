readOpenArray <- function(filename) {

  tmp <- readLines(filename, n=50)
  iStart <- grep("\"Experiment Name\"",tmp)
    
  d <- read.table(filename, skip=iStart-1, header=TRUE, dec=".", sep=",", comment.char="")
  names(d)[names(d)=="Barcode"] <- "Chip.Id"
  names(d)[names(d)=="Well"] <- "Chip.Well"
  snames <- as.character(d$Sample.Name)
  nsnames <- nchar(snames)
  d[,"Sample.Id"] <- gsub(".$","",snames)
  d[,"Feature.Set"] <- substr(snames,nsnames,nsnames)
  names(d)[names(d)=="Target.Name"] <- "Feature.Id"
  names(d)[names(d)=="Cycle.Number"] <- "Cycle"
  names(d)[names(d)=="Rn"] <- "Value"

  d = d[,c("Chip.Id","Chip.Well","Sample.Id","Feature.Set","Feature.Id","Cycle","Value")];    

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
