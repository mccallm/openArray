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
  
  return(d)
}
