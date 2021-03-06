\name{circuit}
\alias{circuit}
\docType{data}

\title{Circuit boards data}
\description{
Number of nonconformities observed in 26 successive samples of 100 printed circuit boards.
Sample 6 and 20 are outside the control limits. Sample 6 was examined by a new inspector and he did not recognize several type of nonconformities that could have been present. Furthermore, the unusually large number of nonconformities in sample 20 resulted from a temperature control problem in the wave soldering machine, which was subsequently repaired.
The last 20 samples are further samples collected on inspection units (each formed by 100 boards).
}
\usage{data(circuit)}
\format{
  A data frame with 46 observations on the following 4 variables.
  \describe{
    \item{sample}{sample number}
    \item{x}{number of defectives in 100 printed circuit boards (inspection unit)}
    \item{size}{sample size}
    \item{trial}{trial sample indicator (TRUE/FALSE)}
  }
}
%\details{}
%\source{}
\references{
Montgomery, D.C. (1991) \emph{Introduction to Statistical Quality Control}, 2nd ed, New York, John Wiley & Sons, pp. 173--175
}
\examples{
data(circuit)
describe(circuit, by = trial)
boxplot(x/size ~ trial, data = circuit)
plot(x/size ~ sample, data = circuit, type="b")
}
\keyword{datasets}
