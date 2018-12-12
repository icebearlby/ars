\name{ars}
\alias{ars}

\title{Perform adaptive rejection sampling for log-concave density functions}

\description{
  \sQuote{ars} provides an easy interface to perform adaptive rejection
  sampling using the algorithm proposed in Gilks et al (1992). The
  function integrates the four stages of the algorithm performed by the
  package \pkg{ars} to perform the complete adaptive rejection sampling
  algorithm. It allows the user to specify the type of concave density
  function to be sampled from, and if the user wants he/she can set the
  bounds of the domain. 
}

\usage{
ars(M,
    lb = -Inf,
    ub = Inf,
    f,
    width = 0.5,
    mod_val = 0)
}

\arguments{
  \item{M}{
    Number of observations to be sampled. The object should be of the
  type \sQuote{integer}
  }

  \item{lb}{
    Lower bound for the domain (by default -Inf). See \sQuote{Details}.
  }

  \item{ub}{
    Upper bound for the domain (by default -Inf). See \sQuote{Details}.
  }

  \item{f}{
    Argument of the class \sQuote{function} specifying a univariate
    log-concave density function. See \sQuote{Details}.
  }

  \item{width}{
    Step value around highest probability density point (by default 0.5). The object
    should be of the type \sQuote{integer}.
  }

  \item{mod_val}{
    Highest probability density point (by default 0). The object
    should be of the type \sQuote{integer}.
  }
}

\details{

  This function is used to conduct adaptive rejection sampling acting
  as a wrapper over all the stages of the algorithm provided in the \pkg{ars}
  package. It provides an interface to select and control the process of
  adaptive rejection sampling. It includes choice of the concave density
  function and the bounds for the domain of the function. The user can
  choose to have infinite bounds. The algorithm for adaptive rejection
  sampling is detailed below:

  \enumerate{
    \item{M: is defined as the number of observations to be sampled from
    the density function.}

    \item{lb: is the lower bound for the domain of the density
      function. If the user decided to have infinite lower bound, the
      function \sQuote{initialize_sample()} initializes the sample by
      creating a new lower bound using the point of highest probabily
      density and step function around it.}

    \item{ub: is the upper bound for the domain of the density
      function. If the user decided to have infinite upper bound, the
      function \sQuote{initialize_sample()} initializes the sample by
      creating a new upper bound using the point of highest probabily
      density and step function around it.}

    \item{f: is the univariate density function provided by the
      user. The function should be a continuous, differentiable, concave
      function. Common log-concave density functions are provided in the
      package. These include: Standard Normal, Gamma (chi-squared),
      Uniform, Beta, Logistic, T-Distribution, F-Distribution,
      Exponential, Log-normal.}
  }
  
  The function divides the sampling algorithm into four stages. After
  the first stage of initialization, the remaining three stages are
  repeated until \sQuote{M} sample points have been generated:
  \enumerate{
    \item{Initialization: See \code{\link{initialize_sample}} for more details.}

    \item{Sampling: See \code{\link{create_samples}} for more details.}

    \item{Squeezing and Rejection tests: See \code{\link{sr_test}} for more details.}

    \item{Updation: See \code{\link{update_hfamily}} and
      \code{\link{update_z}} for more details.}
  }
}	
  
\author{Baoyue Liang, Sargam Jain, Dandan Ru}

\seealso{
  \code{\link{initialize_sample}},
  \code{\link{initialize_hfamily}},	
  \code{\link{initialize_z}},
  \code{\link{create_samples}},
  \code{\link{calculate_scdf}},
  \code{\link{u}},
  \code{\link{l}},
  \code{\link{update_hfamily}},
  \code{\link{update_z}},
}

\examples{ 
## Adaptive rejection sampling for standard normal
es <- ars(M = 10000,
          lb = -1,
          ub = 1,
          f = dnor)

## Adaptive rejection sampling for logistic
es <- ars(M = 10000,
          lb = -1,
          ub = 1,
          f = dlogit)
}

\keyword{ars}
