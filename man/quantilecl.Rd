\name{quantilecl}
\alias{quantilecl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function that applies the quantile classifier for a given set of quantile probabilities and selects the best quantile classifier in the training set.
}
\description{
The function applies the quantile classifier for a set of quantile probabilities and selects the optimal probability that minimize the misclassification rate in the training set. 
}
\usage{
quantilecl(train, test, cl, theta = NULL, 
cl.test = NULL, skew.correct="Galton")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{train }{A matrix of data (the training set) with observations in rows and variables in columns. It can be a matrix or a dataframe.}
  \item{test }{A matrix of data (the test set) with observations in rows and variables in columns. It can be a matrix or a dataframe.}
  \item{cl }{
A vector of class labels for each sample of the training set.  It can be factor or numerical.
}
  \item{theta }{ A vector of quantile probabilities (optional) }
  \item{cl.test }{ If available, a vector of class labels for each sample of the test set (optional)}
  \item{skew.correct }{ Skewness measures applied to correct the skewness direction of the variables. The possibile choices are: Galton's skewness (default), Kelley's skewness and the conventional skewness index based on the third standardized moment}
}
\details{
\code{quantile_cl} carries out the quantile classifier for a set of quantile probabilities and selects the optimal probability that minimize the misclassification rate in the training set. The values of the quantile probabilities can be given in input or automatically selected in a equispaced range of 49 values between 0 and 1. The data in the training and test samples are preprocessed so that the variables used for the quantile estimator all have the same (positive) direction of skewness according to different measures of skewness: Galton's skewness, Kelley's skewness or conventional skewness index.
}
\value{
A list with components
\item{train.rates }{Misclassification errors for each quantile probability in the training set}
\item{test.rates }{Misclassification errors for each quantile probability in the test set}
\item{thetas }{The list of optimal quantile probabilities for each variable}
\item{theta.choice }{The quantile probability that gives the less misclassification error in the training set}
\item{me.train }{Misclassification error in the training set}
\item{me.test }{Misclassification error in the test set (only if \code{cl.test} is available)}
\item{train }{The matrix of data (training set) with observations in rows and variables in columns}
\item{test }{The matrix of data (test set) with observations in rows and variables in columns}
\item{cl.train }{Predicted classification in the training set}
\item{cl.test }{Predicted classification in the test set}
\item{cl.train.0 }{The true classification labels in the training set}
\item{cl.test.0 }{The true classification labels in the test set (if available)} 
}

\author{Christian Hennig, Cinzia Viroli}

\seealso{
See Also \code{\link{quantilecl.vw}}
}
\examples{
data(ais)
x=ais[,3:13]
cl=as.double(ais[,1])
set.seed(22)
index=sample(1:202,152,replace=FALSE)
train=x[index,]
test=x[-index,]
cl.train=cl[index]
cl.test=cl[-index]
out.q=quantilecl(train,test,cl.train,cl.test=cl.test)
out.q$me.test
print(out.q)
plot(out.q)
}
\keyword{multivariate}
