\name{project}
\alias{project}
\title{project}
\description{project}
\usage{
	project(c1, c2, signif = 0.05, eps=1e-5, fuzz=1e-14,
	    show.details=F, warn=T)
}
\arguments {
	\item{c1}{see details}
	\item{c2}{see details}
}
\value{xxx}
\details{ Under Development

	c1 and c2 should be curvature summary statistics as returned by curvature.
	It is assummed that c1 is a submodel of c2.
	The tangent space (parameter effects curvature) of the sub-model c1 is a 
	subspace of the tangent space of the larger model. The acceleration space 
	(first normal space, intrinsic curvature effects) of the sub-model is a 
	subspace of the direct sum of the tangent and acceleration spaces of the 
	larger model, so the intrinsic 
	Tangent and acceleration vectors of the submodel c1 can be projected onto the 
	tangent and acceration spaces of the larger model c2. 
	These are called T1inT2, N1inT2, and N1inN2 (T for tangent, N for normal).
	A second projection (restriction) of c2 onto the tangent and acceleration 
	spaces of c1 is less interesting but may be a useful check. The intrinsic 
	curvature of the larger model should also be intrinsic on the sub-model and
	the parameter effects of the larger model may be partly intrinsic on the 
	sub-model. These two projected onto the intrinsic curvature space of
	the submodel (N2andT2inN1) should compare with the intrinsic curvature of
	the submodel. 

}

\seealso{
	\code{\link{curvature}}
}

\keyword{DSE}
\keyword{ts}

