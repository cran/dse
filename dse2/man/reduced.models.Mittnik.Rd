\name{reduced.models.Mittnik}
\alias{reduced.modeles.Mittnik}
\title{Reduced Models via Mittnik SVD balancing}
\description{Reduced Models via Mittnik SVD balancing.}
\usage{
	reduced.models.Mittnik(largeModel)
}
\arguments{
\item{ largeModel}{An SS TSmodel.}
}
\value{A list of state space TSmodels with smaller state dimensions.}
\details{
The largeModel is balanced by the SVD method promoted by Mittnik (see 
reduction.Mittnik) and then models for every state dimension up to the 
state dimension of the largeModel are return. Note that this procedure 
does not result in smaller models which are balanced.
}
\seealso{
\code{\link{reduction.Mittnik}}
}
\examples{
if(is.R()) data("eg1.DSE.data.diff", package="dse1")
z <-  reduced.models.Mittnik(to.SS(est.VARX.ls(eg1.DSE.data.diff)))
}
\keyword{DSE}
\keyword{ts}

