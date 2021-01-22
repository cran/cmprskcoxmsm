\name{pvalFormat}
\alias{pvalFormat}

\title{
Adjust p-value format
}

\description{
Formats p-values for reports, can report adjusted pvalues
}

\usage{
pvalFormat(p.values, method = 'none', replace = FALSE)
}

\arguments{

\item{p.values}{
p-value
}

\item{method}{
pvalue adjustment, passed to p.adjust.methods
}

\item{replace}{
if TRUE, replaces p-values with their adjusted value
}

}

\value{

Return the formatted p-value: round the p-value, assign the significance sign to the p-value based on the significant level. Can be used directly to the LaTex report.}

\examples{
p <- 0.002354
print(p.1 <- pvalFormat(p))
}
\keyword{p-value}
