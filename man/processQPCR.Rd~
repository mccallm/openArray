\name{tnetfit}
\alias{tnetfit}
\title{Ternary Network Fitting}
\description{
  This function fits a ternary network based on perturbation experiments.
}
\usage{
tnetfit(steadyStateObj, perturbationObj, params=ternaryFitParameters(),
xSeed=NA)
}
\arguments{
  \item{steadyStateObj}{a matrix of steady gene expression observations
    from a perturbation experiment. Rows are genes and columns are experiments.}
  \item{perturbationObj}{a matrix of perturbation experiments. Rows are
    genes and columns are experiments.}
  \item{params}{a ternaryFitParameters object}  
  \item{xSeed}{an integer random seed. If NA, a random seed is generated.}
}
\value{
  The function returns a ternaryFit object.
}
\author{Matthew N. McCall and Anthony Almudevar}
\seealso{Almudevar A, McCall MN, McMurray H, Land H (2011). Fitting
  Boolean Networks from Steady State Perturbation Data, Statistical
  Applications in Genetics and Molecular Biology, 10(1): Article 47.}
\examples{
ssObj <- matrix(c(1,1,1,0,1,1,0,0,1),nrow=3)
pObj <- matrix(c(1,0,0,0,1,0,0,0,1),nrow=3)
rownames(ssObj) <- rownames(pObj) <- colnames(ssObj) <- colnames(pObj) <- c("Gene1","Gene2","Gene3")
tnfitObj <- tnetfit(ssObj, pObj)
}
\keyword{manip}

