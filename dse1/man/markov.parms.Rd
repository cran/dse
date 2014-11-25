\name{markov.parms}
\alias{markov.parms}
\title{Markov Parameters}
\description{Construct a Matrix of the Markov Parameters}
\usage{
    markov.parms(model, blocks=NULL)
}
\arguments {
    \item{model}{An ARMA or SS TSmodel.}
    \item{blocks}{Number of blocks to calculate.}
}
\value{A matrix}
\details{
    Construct a matrix with partitions [M0|...|Mi] giving the Markov
    parameters Mi, i+1 = blocks
    where each Mi is a p by (m+p) matrix, (m is the dimension of the exogeneous 
    series and p is the dimension of endogeneous series)
    ie.  y(t) = e(t) + M [u'(t)|y'(t-1) | u'(t-1)|y'(t-2)]' 
    This requires that models be transformed so that lagged endogeneous variables
    are inputs.  See Mittnik p1190.
    If blocks=NULL (the default) then at least 3 blocks are generated, and
    up to n+1, but the series is truncated if the blocks are effectively zero.
    This will affect the size of the Hankel matrix.
}
\seealso{
    \code{\link{balance.Mittnik.svd}}
}
\references{
    S.Mittnik (1989), Multivariate Time Series Analysis With State 
    Space Models, Computers Math Appl. Vol 17, No 8/9, pp1189-1201.
   
    S.Mittnik (1990), Macroeconomic Forecasting Experience With 
    Balance State Space Models,  
    International Journal Of Forecasting, Vol 6, pp337-348.
  
    S.Mittnik (1990), Forecasting With Balanced State Space 
    Representations of Multivariate Distributed Lag Models. 
    J. of Forecasting, Vol.9, 207-218.  
}
\examples{
    if(is.R()) data("eg1.DSE.data.diff", package="dse1")
    model <- est.VARX.ls(eg1.DSE.data.diff)
    markov.parms(model)
}
%\keyword{DSE}
\keyword{ts}
\keyword{algebra}

