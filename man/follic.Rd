\name{follic}
\alias{follic}
\docType{data}
\title{Follicular cell lymphoma study}
\description{
Competing risk data set involving follicular cell lymphoma from Pintilie (2007)
}
\usage{data(follic)}
\format{
A data frame containing:
  \describe{
    \item{\code{age}}{age}
    \item{\code{hgb}}{hemoglobin (g/l)}
    \item{\code{clinstg}}{clinical stage: 1=stage I, 2=stage II}
    \item{\code{ch}}{chemotherapy}
    \item{\code{rt}}{radiotherapy}
    \item{\code{time}}{	first failure time}
    \item{\code{status}}{Reason of failure: 1: Relapse, 2: Death, 0: No response}
  }
}

\references{
Pintilie M., (2006) Competing Risks: A Practical Perspective. West Sussex: John Wiley and Sons.
}

\keyword{datasets}
