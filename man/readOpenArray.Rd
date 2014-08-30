\name{readOpenArray}
\alias{readOpenArray}
\title{Read in OpenArray Amplification Curves}
\description{
  This function reads in the raw amplification curve data exported by
  the OpenArray software.
}
\usage{
  readOpenArray(filename)
}
\arguments{
  \item{filename}{the name of the csv file containing the raw data exported from the LifeTech OpenArray software}
}
\value{
  The function returns an object in the format required by the
  processQPCR function. 
}
\author{Alexander Baras and Matthew N. McCall}
\examples{
  dat <- readOpenArray(
             system.file("extdata/rawData.csv", package = "openArray"))
}
\keyword{manip}

