\name{centroidcl}
\alias{centroidcl}
\title{
A function that performs the centroid classifier
}
\description{
Given a training and a test set, the function apply the centroid classifier and returns the classification labels of the observations in the training and in test set. It also gives the training misclassification rate and the test misclassification rate, if the truth class labels of the test set are provided in input.
}
\usage{
centroidcl(train, test, cl, cl.test = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{train}{ A matrix of data (the training set) with observations in rows and variables in column. It can be a matrix or a dataframe.
}
  \item{test}{ A matrix of data (the test set) with observations in rows and variables in columns. It can be a matrix or a dataframe.
}
  \item{cl}{ A vector of class labels for each sample of the training set. It can be factor or numerical.
}
  \item{cl.test}{ A vector of class labels for each sample of the test set (optional)
}
}
\details{
\code{centroidcl} carries out the centroid classifier and predicts classification.
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
See Also \code{\link{theta.cl}}
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
out.c=centroidcl(train,test,cl.train,cl.test)
out.c$me.test
misc(out.c$cl.test,cl.test)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
