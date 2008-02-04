
%%% $Id: utilsForExamples.Rd,v 1.1 2008/02/04 19:58:52 goswami Exp $

\name{utilsForExamples}

\alias{MarkovSwitchingFuncGenerator}

\title{The utility function(s) for examples}

\description{
  The utility function(s) that are used in the example sections of the
  exported functions in this package.
}

\usage{
MarkovSwitchingFuncGenerator(seed = -975313579)
}

\arguments{
  \item{seed}{the seed for random number generation.}
}

\value{
  A list containing the objects to be used as arguments to the exported
  functions in the respective example sections of this package. 
}

\author{Gopi Goswami \email{goswami@stat.harvard.edu}}

\seealso{\code{\link{particleFilter}},
  \code{\link{auxiliaryParticleFilter}},
  \code{\link{sequentialMonteCarlo}}}

\keyword{datagen}

