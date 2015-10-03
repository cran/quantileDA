\name{ais}
\alias{ais}
\docType{data}
\title{Australian Institute of Sport data}
\description{Data on 102 male and 100 female athletes collected at the Australian Institute of Sport, courtesy of Richard Telford and Ross Cunningham.}
\usage{data(ais)}
\format{
  A data frame with 202 observations on the following 13 variables.
  \describe{
    \item{\code{sex}}{A factor with levels \code{female} \code{male}}
    \item{\code{sport}}{A factor with levels \code{B_Ball} \code{Field} \code{Gym} \code{Netball} \code{Row} \code{Swim} \code{T_400m} \code{T_Sprnt} \code{Tennis} \code{W_Polo}}
    \item{\code{rcc}}{A numeric vector: red cell count}
    \item{\code{wcc}}{A numeric vector: white cell count}
    \item{\code{Hc}}{A numeric vector: Hematocrit}
    \item{\code{Hg}}{A numeric vector: Hemoglobin}
    \item{\code{Fe}}{A numeric vector: plasma ferritin concentration}
    \item{\code{bmi}}{A numeric vector: body mass index}
    \item{\code{ssf}}{A numeric vector: sum of skin folds}
    \item{\code{Bfat}}{A numeric vector: body fat percentage}
    \item{\code{lbm}}{A numeric vector:  lean body mass}
    \item{\code{Ht}}{A numeric vector: height (cm)}
    \item{\code{Wt}}{A numeric vector: weight (kg)}
  }
}

\source{
Cook and Weisberg (1994), An Introduction to Regression Graphics. John Wiley & Sons, New York. 
}
\examples{
data(ais)
attach(ais)
pairs(ais[,c(3:4,10:13)], main = "AIS data")
plot(Wt~sport)
}
\keyword{datasets}
