\name{theta.cl}
\alias{theta.cl}

\title{
A function to perform the quantile classifier for a given quantile probability
}
\description{
Given a certain quantile probability, the function compute the quantile classifier on the training set and gives the predicted class labels in the training and test set.It also computes the training misclassification rate and the test misclassification rate, when the truth labels of the test set are available. When the quantile probability is 0.5 the function compute the median classifier.
}
\usage{
theta.cl(train, test, cl, theta, cl.test = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{train }{A matrix of data (the training set) with observations in rows and variables in columns. It can be a matrix or a dataframe.}
  \item{test }{A matrix of data (the test set) with observations in rows and variables in columns.  It can be a matrix or a dataframe.}
  \item{cl }{A vector of class labels for each sample of the training set. It can be factor or numerical.}
  \item{theta }{The quantile probability. If 0.5 the median classifier is applied}
  \item{cl.test }{If available, a vector of class labels for each sample of the test set (optional)}
}
\details{
\code{theta.cl} carries out quantile classifier for a given quantile probability.
}
\value{
A list with components
 \item{cl.train }{Predicted classification in the training set}
 \item{cl.test }{Predicted classification in the test set}
 \item{me.train }{Misclassification error in the training set}
 \item{me.test }{Misclassification error in the test set (only if \code{cl.test} is available)}
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{Christian Hennig, Cinzia Viroli}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See Also \code{\link{centroidcl}}
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
out.m=theta.cl(train,test,cl.train,0.5,cl.test)
out.m$me.test
misc(out.m$cl.test,cl.test)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
