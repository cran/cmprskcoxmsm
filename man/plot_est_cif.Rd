\name{plot_est_cif}
\alias{plot_est_cif}

\title{
ggplot method for cif_est objects
}

\description{
This function produces a CIF plots for \code{cif_est} objects.
}

\usage{
plot_est_cif(cif.data, color = color, ci.cif = FALSE)
}

\arguments{

\item{cif.data}{
The dataset, output of \code{cif_est}
}

\item{color}{
Color for control and treatment group, has to be a vector of length 2.
}

\item{ci.cif}{
Indicator whether to plot the 95\% confidence interval area for the CIF.
}

}

\value{
A cumulative incidence function plot (with 95\% confidence interval area) corresponding to the output from \code{cif_est}.
}



\keyword{cumulative incidence function}


