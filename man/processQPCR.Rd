\name{processQPCR}
\alias{processQPCR}
\title{Processing of OpenArray Amplification Data}
\description{
  This function fits a 5 parameter curve to qPCR amplification measurements from the OpenArray platform and returns expression estimates and quality scores.
}
\usage{
  processQPCR(dat)
}
\arguments{
  \item{dat}{a data.frame containing the amplification measurements}
}
\value{
  The function returns an HTqPCR object containing expression estimates and quality scores.
}
\author{Alexander Baras and Matthew N. McCall}
\examples{

}
\keyword{manip}

