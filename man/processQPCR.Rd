\name{processQPCR}
\alias{processQPCR}
\title{Processing of OpenArray Amplification Curves}
\description{
  This function fits a 5 parameter curve to amplification curves from the OpenArray platform and returns expression estimates and quality scores.
}
\usage{
  processQPCR(filename, fileFormat="default")
}
\arguments{
  \item{filename}{the name of the file containing the raw data}
  \item{fileFormat}{the format of the raw data file. Current options are "default" or "LifeTech"}
}
\value{
  The function returns an HTqPCR object containing expression estimates and quality scores.
}
\author{Alexander Baras and Matthew N. McCall}
\examples{

}
\keyword{manip}

