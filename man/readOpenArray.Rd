\name{readOpenArray}
\alias{readOpenArray}
\title{Read in OpenArray Amplification Curves}
\description{
  This function reads in the raw amplification curve data exported by
  the OpenArray software.
}
\usage{
  readOpenArray(filename, fileFormat="default")
}
\arguments{
  \item{filename}{the name of the file containing the raw data}
  \item{fileFormat}{the format of the raw data file. Current options are "default" or "LifeTech"}
}
\value{
  The function returns a data.frame in the format required by the
  processQPCR function. 
}
\author{Alexander Baras and Matthew N. McCall}
\examples{
  
}
\keyword{manip}

