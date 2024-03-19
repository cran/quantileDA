\name{quantilecl.vw}
\alias{quantilecl.vw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to apply the quantile classifier that uses a different optimal quantile probability for each variable}
\description{
A function to apply the quantile classifier that uses a different optimal quantile probability for each variable
}
\usage{
quantilecl.vw(train, test, cl, theta = NULL, cl.test = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{train }{
A matrix of data (the training set) with observations in rows and variables in columns. It can be a matrix or a dataframe.
}
  \item{test }{
A matrix of data (the test set) with observations in rows and variables in columns. It can be a matrix or a dataframe.
}
  \item{cl }{
A vector of class labels for each sample of the training set.  It can be factor or numerical.
}
  \item{theta }{Given $p$ variables, a vector of length $p$ of quantile probabilities (optional) }
  \item{cl.test }{ If available, a vector of class labels for each sample of the test set (optional)}
}
\details{
\code{quantilecl.vw} carries out the quantile classifier by using a different optimal quantile probability for each variable selected in the training set.
}
\value{
A list with components
\item{Vseq }{The value of the objective function at each iteration}
\item{thetas }{The vector of quantile probabilities}
\item{me.train }{Misclassification error for the best quantile probability in the training set}
\item{me.test }{Misclassification error for the best quantile probability in the test set (only if \code{cl.test} is available)}
\item{cl.train }{Predicted classification in the training set}
\item{cl.test }{Predicted classification in the test set}
\item{lambda }{The vector of estimated scale parameters}
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{Marco Berrettini, Christian Hennig, Cinzia Viroli}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See Also \code{\link{quantilecl}}
}
\examples{
data(ais)
x=ais[,3:7]
cl=as.double(ais[,1])
set.seed(22)
index=sample(1:202,152,replace=FALSE)
train=x[index,]
test=x[-index,]
cl.train=cl[index]
cl.test=cl[-index]
out.q=quantilecl.vw(train,test,cl.train,cl.test=cl.test)
out.q$me.test
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
