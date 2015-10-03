\name{galtonskew}
\alias{galtonskew}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function that compute the Galton's skewness
}
\description{
The function compute the Galton's skewness index on a set of observations.
}
\usage{
galtonskew(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x }{
A vector of observations.
}}

\value{
A scalar which measures the Galton's skewness}


\author{Christian Hennig, Cinzia Viroli}

\seealso{
See Also \code{\link{kelleyskew}}
}
\examples{
data(ais)
galtonskew(ais[,4])
}
\keyword{multivariate}
