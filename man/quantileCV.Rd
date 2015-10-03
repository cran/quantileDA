\name{quantileCV}
\alias{quantileCV}

\title{
 A function to cross-validate the quantile classifier
}
\description{
Balanced cross-validation for the quantile classifier
}
\usage{
quantileCV(x, cl, nfold = min(table(cl)), 
folds = balanced.folds(cl, nfold), theta=NULL, seed = 1, varying = FALSE)
}
\arguments{
  \item{x}{ A matrix of data (the training set) with observations in rows and variables in columns (it can be a matrix or a dataframe)}
  \item{cl}{ A vector of class labels for each sample (factor or numerical)}
  \item{nfold}{ Number of cross-validation folds. Default is the smallest class size. Admitted values are from 1 to the smallest class size as maximum fold number.}
  \item{folds}{ A list with nfold components, each component a vector of indices of the samples in that fold. By default a (random) balanced cross-validation is used}
\item{theta }{ A vector of quantile probabilities (optional) }
  \item{seed}{ Fix the seed of the running. Default is 1}
  \item{varying}{ If TRUE a different quantile for each variable is selected in the training set. If FALSE (default) an unique quantile is used.}
}
\details{\code{quantileCV} carries out cross-validation for a quantile classifier.}

\value{
A list with components
 \item{test.rates }{Mean of misclassification errors in the cross-validation test sets for each quantile probability (available if \code{varying} is FALSE)}
 \item{train.rates }{Mean of misclassification errors in the cross-validation train sets for each quantile probability (available if \code{varying} is FALSE)}
 \item{thetas }{The fitted quantile probabilities}
 \item{theta.choice }{Value of the chosen quantile probability in the training set}
 \item{me.test }{Misclassification errors in the cross validation test sets for the best quantile probability}
 \item{me.train }{Misclassification errors in the cross validation training sets for the best quantile probability}
 \item{me.median }{Misclassification errors in the cross validation test sets of the median classifier}
 \item{me.centroid }{Misclassification errors in the cross validation test sets of the centroid classifier}
 \item{folds }{The cross-validation folds used}
}

\author{Christian Hennig, Cinzia Viroli}


\examples{
data(ais)
x=ais[,3:13]
cl=as.double(ais[,1])
out=quantileCV(x,cl,nfold=10)
}

\keyword{multivariate}
