\name{cif_est}
\alias{cif_est}

\title{
Estimated cumulative incidence function
}

\description{
\code{cif_est} estimates the cumulative incidence function (CIF, i.e.risk) based on the cause-specific regression results with 95\% confidence intervals, it also calculates the risk ratio and risk difference for the specific time point.
}

\usage{
cif_est(data=,
        time, time2 = NULL,
        Event.var, Events, cif.event,
        weight.type,
        ties = NULL,
        risktab = TRUE, risk.time = NULL)
}

\arguments{

\item{data}{
The dataset, output of \code{doPS}
}

\item{time}{
See \code{weight_cause_cox}.
}

\item{time2}{
See \code{weight_cause_cox}.
}

\item{Event.var}{
The variable name for the event indicator which typically has at least 3 levels.
}

\item{Events}{
The vector of all the event name, the rest of levels in the \code{Event.var} will be treated as loss to follow up (i.e. right censoring).
}

\item{cif.event}{
Value of event of interest for the CIF.
}

\item{weight.type}{
See \code{weight_cause_cox}.
}

\item{ties}{
See \code{weight_cause_cox}.
}

\item{risktab}{
Indicator whether the risk ratio and risk difference table should be returned.
}

\item{risk.time}{
If \code{risktab}, the specific time point for calculating the risk ratio and risk difference.
}

}

\value{

Returns a table containing the estimated CIF for the event of interest for control and treated group and their 95\% confidence intervals.

If \code{risktab}, will return the risk ratio and risk difference at time \code{risk.time}, and their 95\% confidence intervals.


}

\details{

After estimating the parameters in the cause-specific hazard \eqn{\lambda_{j}^a} using IPW, we could estimate the corresponding CIF:

\deqn{
    \hat{P}(T^a<t,J^a=j) = \int_{0}^{t} \hat{S}^a(u) d\hat{\Lambda}_{j}^a(u),
}

where \eqn{\hat{S}^a} is the estimated overall survial function for \eqn{T^a}, \eqn{\hat{S}^a(u) = e^{-\sum_j\hat{\Lambda}_{j}^a(u)}}, \eqn{\hat\Lambda_{j}^a(u) = \hat\Lambda_{0j}(u)e^{\hat\beta*a}}, and \eqn{\hat\Lambda_{0j}(u)} is a Breslow-type estimator of the baseline cumulative hazard.



}

\keyword{cumulative incidence function}

\references{
Hou, J., Paravati, A., Hou, J., Xu, R., & Murphy, J. (2018). \dQuote{High-dimensional variable selection and prediction under competing risks with application to SEER-Medicare linked data,}
\emph{Statistics in Medicine} 37(24), 3486-3502.
}
