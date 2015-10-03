\name{kelleyskew}
\alias{kelleyskew}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function that compute the Kelley's skewness
}
\description{
The function compute the Kelley's skewness index on a set of observations.
}
\usage{
kelleyskew(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x }{
A vector of observations.
}}

\value{
A scalar which measures the Kelley's skewness}

\author{Christian Hennig, Cinzia Viroli}

\seealso{
See Also \code{\link{galtonskew}}
}
\examples{
data(ais)
kelleyskew(ais[,4])
}
\keyword{multivariate}
