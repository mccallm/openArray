readOpenArray <- function(filename) {

  tmp <- readLines(filename, n=50, encoding="C")
  iStart <- grep("Experiment Name", tmp, perl=TRUE)[1]
    
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
  
  dat = split(d,paste(d$Chip.Id,d$Chip.Well,d$Sample.Id,d$Feature.Set,d$Feature.Id,sep="::"))
  
  return(dat)
}
