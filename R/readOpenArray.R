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
  d = d[order(d$Chip.Id,d$Chip.Well,d$Sample.Id,d$Feature.Set,d$Feature.Id),];

  #Some basic error check
  if(is.null(d$Value)){
    stop("Error: You must have a Value column for the raw fluorescence values.");
  }

  nRecords = nrow(d)
  dat = list()
  c = d[1,c("Chip.Id","Chip.Well","Sample.Id","Feature.Set","Feature.Id")]
  j = 1
  i = 1
  
  while (i <= nRecords) {
  	dat[[j]] = list("Chip.Id"=c$Chip.Id,"Chip.Well"=c$Chip.Well,"Sample.Id"=c$Sample.Id,"Feature.Set"=c$Feature.Set,"Feature.Id"=c$Feature.Id)
  	while ((i <= nRecords) & (prod(c == d[i,c("Chip.Id","Chip.Well","Sample.Id","Feature.Set","Feature.Id")])==1)) {
  		dat[[j]]$x = c(dat[[j]]$x,d$Cycle[i])
  		dat[[j]]$y = c(dat[[j]]$y,d$Value[i])
  		i = i+1
  	}
  	dat[[j]]$SST = sum((dat[[j]]$y-mean(dat[[j]]$y))^2)
	  yRange = as.numeric(quantile(dat[[j]]$y,c(0.1,0.9)))
	  dat[[j]]$init = c(yRange[1],yRange[2]-yRange[1],which.min(abs(mean(yRange)-dat[[j]]$y)),0.5,1)
	
  	c = d[i,c("Chip.Id","Chip.Well","Sample.Id","Feature.Set","Feature.Id")]
	  j = j+1
  }
  return(dat)
}
