\name{doPS}
\alias{doPS}

\title{
Generate the Inverse Probability Treatment Weights
}

\description{
\code{doPS} calculates the unstabilized and stabilized inverse probability treatment weights (IPW) for average treatment effect using propensity score. The propensity score is calculated by \code{twang} package using the boosted logistic regression.
}

\usage{
doPS(data,Trt,Trt.name,VARS.)
}

\arguments{

\item{data}{
The dataset, includes treatment assignment as well as covariates
}

\item{Trt}{
The name of the treatment variable in the dataset.
}

\item{Trt.name}{
The treated group name of the treatment variable in the dataset.
}

\item{VARS.}{
The vector of the name of potential confounding variables in the dataset.
}

}

\value{

doPS returns an object of class "PS".
An object of class "PS" is a list containing the following components:

\item{Data}{A new dataset which excludes all the missing value on the potential confounders from input data, add the propensity score and IPW into the new dataset.
 \describe{
   \item{ps_ate}{The estimated propensity scores with estimand of interest as ATE}
   \item{ipw_ate_unstab}{Unstabilized ipw calculated from ps_ate}
   \item{ipw_ate_stab}{Stabilized ipw calculated from ps_ate}
  }}

\item{ps}{an object of class \code{ps}, See the help for \link[twang]{ps} for details of the ps class.
}

}

\details{

The treatment variable should only contain 2 levels of treatment, and one should be viewed as treated group and another is control group.

For stabilized weights:

For the treated individuals, we assign the IPW: w = Pr(T=1)/Pr(T=1|X=x), for control individuals, the stabilized weight is: w = (1-Pr(T=1))/(1-Pr(T=1|X=x)).

}

\keyword{propensity score}
\keyword{inverse probability weighting}

\seealso{ \code{\link[twang]{ps}} }
