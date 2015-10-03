\name{skewness}
\alias{skewness}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function that compute the conventional skewness measure
}
\description{
A function that compute the conventional skewness measure according to the third standardized moment of x
}
\usage{
skewness(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x }{
A vector of observations.
}}

\value{
A scalar which measures the skewness}


\author{Christian Hennig, Cinzia Viroli}

\seealso{
See Also \code{\link{galtonskew}}
}
\examples{
data(ais)
skewness(ais[,4])
}
\keyword{multivariate}
