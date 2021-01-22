\name{weight_cause_cox}
\alias{weight_cause_cox}

\title{
Inverse probability weighted proportional cause-specific hazards model
}

\description{
\code{weight_cause_cox} fits the marginal structural proportional cause-specific hazards model using the inverse probability treatment weights.
}

\usage{
weight_cause_cox(data=,
                 time, time2 = NULL,
                 Event.var, Event,
                 weight.type,
                 ties = NULL)
}

\arguments{

\item{data}{
The dataset, output of \code{doPS}
}

\item{time}{
See also \code{\link[survival]{Surv}}, for right censored data, this is the follow up time. For interval data, the first argument is the starting time for the interval.
}

\item{time2}{
See also \code{\link[survival]{Surv}}, ending time of the interval for interval censored or counting process data only. Intervals are assumed to be open on the left and closed on the right, (start, end]. For counting process data, event indicates whether an event occurred at the end of the interval.
}

\item{Event.var}{
The variable name for the event indicator which typically has at least 3 levels.
}

\item{Event}{
Event of interest, the rest of the event are treating as competing event.
}

\item{weight.type}{
Type of inverse probability weights. Possible values are "Unstabilized" and "Stabilized".
}

\item{ties}{
See also \code{\link[survival]{coxph}}, a character string specifying the method for tie handling. If there are no tied death times all the methods are equivalent.
}

}

\value{

Returns a table containing the estimated coefficient of the treatment effect, the robust standard error of the coefficient, estimated hazard ratio and 95\% CI for the hazard ratio.


}

\details{

The marginal structural cause-specific Cox model for cause j usually has the form:

\deqn{
    \lambda^{a}_j (t) \equiv \lambda_{T^{a},J^{a}=j}(t) = \lambda_{0j}e^{\beta*a},
}

where \eqn{T^{a}}, \eqn{J^{a}} is the counterfactural survival time and cause for treatment \eqn{a (=0,1)}, \eqn{\lambda_{0j}} is the unspecified baseline cause-specific hazard for cause j, and \eqn{\beta} is the treatment effect.

}

\keyword{cause specific cox}
\keyword{inverse probability weighting}

\seealso{ \code{\link[survival]{coxph}} }
