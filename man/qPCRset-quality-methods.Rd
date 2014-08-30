\name{qPCRset-quality-methods}
\alias{qPCRset-quality-methods}
\alias{quality,qPCRset-method}
\alias{quality<-,qPCRset-method}
\alias{quality}
\alias{quality<-}
\title{
  Accessor and replacement functions for the qPCRset quality field.
}
\description{These functions allow the user to access and replace an
  AssayData element called quality for a qPCRset object. The quality
  field is used by the processQPCR function to store a measure of the
  quality of each expression estimate.}

\examples{
  dat <- readOpenArray(
             system.file("extdata/rawData.csv", package = "openArray")) 
  obj <- processQPCR(dat)
  qc <- quality(obj)
}

\keyword{manip}