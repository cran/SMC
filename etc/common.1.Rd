
\name{common.1}

\alias{common.1}

\title{The Common stuff for particleFilter, auxiliaryParticleFilter and
sequentialMonteCarlo}

\arguments{
  See the sections \emph{Details}, \emph{Required Functions} and
  \emph{Optional Functions} for explanation on the arguments and the
  return values of the \emph{following arguments that are themselves
  functions}.
  
  \item{nStreams}{\code{integer} \eqn{>} 0.}
  
  \item{nPeriods}{\code{integer} \eqn{>} 0.}
  
  \item{dimPerPeriod}{\code{integer} \eqn{>} 0.}
  
  \item{resampCriterionFunc}{\code{function} of four arguments
    \code{(currentPeriod, currentStreams, currentLogWeights, \dots)}.}
  
  \item{resampFunc}{\code{function} of four arguments
    \code{(currentPeriod, currentStreams, currentLogWeights, \dots)}.}
  
  \item{summaryFunc}{\code{function} of four arguments
    \code{(currentPeriod, currentStreams, currentLogWeights, \dots)}.}

  \item{nMHSteps}{\code{integer} \eqn{\ge}{>=} 0.}

  \item{MHUpdateFunc}{\code{function} of six arguments
    \code{(currentPeriod, nMHSteps, currentStreams, lag1Streams,
      lag1LogWeights, \dots)}.}
    
  \item{nStreamsPreResamp}{\code{integer} \eqn{>} 0.}
  
  \item{returnStreams}{\code{logical}.}

  \item{returnLogWeights}{\code{logical}.}
  
  \item{verboseLevel}{\code{integer}, a value \eqn{\ge}{>=} 2 produces a
    lot of output.}

}

\details{
  We introduce the following terms, which will be used in the sections
  \emph{Required Function} and \emph{Optional Function} below:

  \describe{

    \item{\code{stream}}{the state vector also called the particle, the
      hidden state or the latent variable. Below we will use the terms
      stream and state vector interchangeably.}
    
    \item{\code{dimPerPeriod}}{the dimension of the space, the state
      vectors live in.}

  }
}

\section{Optional function: resampCriterionFunc}{

  \describe{

    \item{Arguments:}{The following argument(s) require some
      explanation:}

  }
  
  \describe{

    \item{\code{currentStreams}}{a matrix with \code{dimPerPeriod}
      columns, the rows containing the updated streams for
      \code{currentPeriod}.}

    \item{\code{currentLogWeights}}{a vector of log weights
      corresponding to the streams in the argument matrix
      \code{currentStreams}.}

  }

  \describe{      

    \item{Return value:}{\code{TRUE} or \code{FALSE} reflecting the
      decision of the resampling scheme implemented by this function.}

  }

  \describe{

    \item{Note:}{The following points are in order:}

  }
  
  \describe{

    \item{--}{resampling schemes manily depend on
      \code{currentLogWeights}, the other two arguments might come in
      handy for implementing period or stream specific resampling
      schemes.}


    \item{--}{if \code{nStreamsPreResamp} > \code{nStreams}, then this
      function should always return \code{TRUE}.}

  }
}

\section{Optional function: resampFunc}{
  \describe{

    \item{Arguments:}{see the sub-section \emph{Arguments:} for section
      \emph{Optional function: resampCriterionFunc}.}

  }

  \describe{  

    \item{Return value:}{a \emph{named} list with the following
      components:}

  }
  
  \describe{

    \item{\code{currentStreams}}{a matrix of dimension \code{nStreams}
      \eqn{\times}{x} \code{dimPerPeriod}. The rows of this matrix
      contain the streams for period \code{currentPeriod + 1} that were
      resampled from those of the argument \code{currentStreams} matrix,
      which may contain \eqn{\ge}{>=} \code{nStreams} rows.}

    \item{\code{currentLogWeights}}{The log weights vector of length
      \code{nStreams}, associated with the streams that were resampled
      in the returned \code{currentStreams} matrix. Note, after the
      resampling step, usually all the log weights are set to 0.}

  }

  \describe{

    \item{Note:}{the components of the list returned by this function
      and the arguments to this function have two common names, namely,
      \code{currentStreams} and \code{currentLogWeights}. These entities
      have different meanings, as explained above. For example, the
      argument matrix \code{currentStreams} could possibly have
      \eqn{\ge}{>=} \code{nStreams} rows, whereas the returned
      \code{currentStreams} has exactly \code{nStreams} number of
      (resampled) streams in its rows.}

  }
}

\section{Optional function: summaryFunc}{
  \describe{

    \item{Arguments:}{The following argument(s) require some
      explanation:}

  }
  
  \describe{

    \item{\code{currentStreams}}{a matrix of dimension \code{nStreams}
      \eqn{\times}{x} \code{dimPerPeriod} of streams for
      \code{currentPeriod}.}

    \item{\code{currentLogWeights}}{a vector of log weights
      corresponding to the streams in the argument matrix
      \code{currentStreams}.}

  }

  \describe{

    \item{Return value:}{a vector of length of \code{dimSummPerPeriod}
      of summaries for \code{currentPeriod} given the
      \code{currentStreams} and the \code{currentLogWeights}.}

  }
}

\section{Optional function: MHUpdateFunc}{
  \describe{
    
    \item{Arguments:}{The following argument(s) require some
      explanation:}

  }
  
  \describe{

    \item{\code{nMHSteps}}{the number of Metropolis Hastings (MH) steps
      (iterations) to be performed.}
    
    \item{\code{currentStreams}}{a matrix of dimension \code{nStreams}
      \eqn{\times}{x} \code{dimPerPeriod} of streams for
      \code{currentPeriod}.}

    \item{\code{lag1Streams}}{a matrix of dimension \code{nStreams}
      \eqn{\times}{x} \code{dimPerPeriod} of streams for
      \code{currentPeriod - 1}.}

    \item{\code{lag1LogWeights}}{a vector of length \code{nStreams} of
      log weights corresponding to the streams in the argument matrix
      \code{lag1Streams}.}

  }

  \describe{

    \item{Return value:}{a \emph{named} list with the following
      components:}

  }
  
  \describe{        

    \item{\code{currentStreams}}{a matrix of dimension \code{nStreams}
      \eqn{\times}{x} \code{dimPerPeriod}. The rows of this matrix
      contain the streams for period \code{currentPeriod} that are
      (possibly) MH-updated versions of the rows of the argument
      \code{currentStreams} matrix.}

    \item{\code{acceptanceRates}}{a vector of length \code{nStreams},
      representing the acceptance rates of the \code{nMHSteps}-many MH
      steps for each of the streams in the rows of the argument
      \code{currentStreams} matrix.}

  }

  \describe{

    \item{Note:}{a positive value of \code{nMHSteps} performs as many MH
      steps on the rows of the argument \code{currentStreams}
      matrix. This is done to reduce the possible degeneracy after the
      resampling.}

  }
}

\value{
  This function returns a list with the following components:
  
  \item{draws}{a list with the following components: \code{summary},
    \code{propUniqueStreamIds}, \code{streams}, \code{logWeights},
    \code{acceptanceRates}. See the section \emph{Note} for more
    details.}

  \item{nStreams}{the \code{nStreams} argument.}

  \item{nPeriods}{the \code{nPeriods} argument.}

  \item{dimPerPeriod}{the \code{dimPerPeriod} argument.}

  \item{nStreamsPreResamp}{the \code{nStreamsPreResamp} argument.}

  \item{nMHSteps}{the \code{nMHSteps} argument.}
  
  \item{time}{the time taken by the run.}  

}

\note{
  The effect of leaving the default value \code{NULL} for some of the
  arguments above are as follows:

  \describe{

    \item{\code{resampCriterionFunc}}{the builtin resampling criterion,
      namely, resample when square of the coefficient of variation of
      the weights \eqn{\ge}{>=} 1, is used.}

    \item{\code{resampFunc}}{the builtin resampling function, which
      resamples streams with probability proportional to their weights,
      is used.}

    \item{\code{summaryFunc}}{the builtin summary function, which
      returns the weighted average of each of the \code{dimPerPeriod}
      dimensions, is used.}

    \item{\code{nStreamsPreResamp}}{it is set to \code{nStreams}.}

  }

  This function returns a list with component called \code{draw}. The
  detailed description of this component, as promised in section
  \emph{Value}, is as follows. It is a list itself with the following
  components:

  \describe{

    \item{\code{summary}}{a matrix of dimension \code{nPeriods}
      \eqn{\times}{x} \code{dimSummPerPeriod}.}

    \item{\code{propUniqueStreamIds}}{a vector of length
      \code{nPeriods}. The values are either proportions of unique
      stream ids accpeted (at each period) if resampling was done or
      \code{NA}.}
      
    \item{\code{streams}}{an array of dimension \code{nStreams}
      \eqn{\times}{x} \code{dimPerPeriod} \eqn{\times}{x}
      \code{nPeriods}. This is returned if \code{returnStreams = TRUE}.}

    \item{\code{logWeights}}{a matrix of dimension \code{nStreams}
      \eqn{\times}{x} \code{nPeriods}. This is returned if
      \code{returnLogWeights = TRUE}.}

    \item{\code{acceptanceRates}}{a matrix of dimension \code{nStreams}
      \eqn{\times}{x} \code{nPeriods}. This is returned if
      \code{nMHSteps > 0}.}

  }
}

\section{Warning}{
  Using very small values (\eqn{\le}{<=} \code{1e3}) for \code{nStreams}
  might not give reliable results.
}

\author{Gopi Goswami \email{goswami@stat.harvard.edu}}


