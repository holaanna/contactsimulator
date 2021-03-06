% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/contact.R
\name{contact}
\alias{contact}
\title{Constructor of the basic contact object}
\usage{
contact(data, rast, initializer, times, t0, ..., grid.lines, pop.grid,
  params, grid.size, age.level, nkernls, age.dist, t.max, t.intervention,
  kernelmodel, obsnames, statenames, paramnames, covarnames, zeronames,
  rmodel, r.model, PACKAGE, d.model)
}
\arguments{
\item{data}{either a data frame holding the time series data,
or an object of class \sQuote{contact},
i.e., the output of another \pkg{contact} calculation.}

\item{rast}{an object of class rasterlayer containing the landscape informationn.}

\item{initializer}{the initial state contaning
\describe{
  \item{row}{The row at which the cell containing the host lies in}
  \item{col}{The column at which the cell containing the host lies in}
  \item{typ}{The type of the host (e.g. backyard 1 /plantation 0 hosts)}
  \item{x}{x-coordinate}
  \item{y}{y-coordinate}
  \item{t_e}{exposure time}
  \item{t_i}{infection time}
 }}

\item{times}{the times at which observations are made.
\code{times} must indicate the column of observation times by name or index.
The time vector must be numeric and non-decreasing.
Internally, \code{data} will be internally coerced to an array with storage-mode \code{double}.}

\item{t0}{The zero-time, i.e., the time of the initial state.
This must be no later than the time of the first observation, i.e., \code{t0 <= times[1]}.}

\item{\dots}{additional arguments supply new or modify existing model characteristics or components.
See \code{\link{contact}} for a full list of recognized arguments.

When named arguments not recognized by \code{\link{contact}} are provided, these are made available to all basic components via the so-called \dfn{userdata} facility.
This allows the user to pass information to the basic components outside of the usual routes of covariates (\code{covar}) and model parameters (\code{params}).}

\item{grid.lines}{A 6 columns data frame with columns names as coor_x_1, coor_y_1, coor_x_2, coor_y_2, orient_line.
\describe{
    \item{coor_x_1, coor_y_1}{Coordinates of the left end point of the grid line }
    \item{coor_x_2, coor_y_2}{Coordinates of the right end point of the grid line }
    \item{orient_line}{Line orientation}
    \enumerate{
       \item indicates horizontal orientation
       \item indicates vetical orientation
    }
    \item{k_line}{Line numbering: bottom to top, then left to right}
}}

\item{pop.grid}{Population density of the grid a case resides. This is filled from bottom to top, then left to right.}

\item{params}{optional; named numeric vector of parameters.
This will be coerced internally to storage mode \code{double}.}

\item{age.level, age.dist}{Vectors of age level and the propportion of each age group respectively. See details.}

\item{obsnames}{optional character vector;
names of the observables.
It is not usually necessary to specify \code{obsnames} since, by default,
these are read from the names of the data variables.}

\item{statenames}{optional character vector;
names of the latent state variables.
It is typically only necessary to supply \code{statenames} when C snippets are in use.}

\item{paramnames}{optional character vector;
names of model parameters.
It is typically only necessary to supply \code{paramnames} when C snippets are in use.}

\item{covarnames}{optional character vector;
names of the covariates.
It is not usually necessary to specify \code{covarnames} since, by default,
these are read from the names of the covariates.}

\item{t_intervention}{Start of the intervention if any.}

\item{EI_model}{Take integer values to specify the type of model used for the latent period. See \code{\link{E_to_I}}}

\item{kern.}{model Take integer values to specify the type of dispersal kernel used. See \code{\link{Samp_dis}}}

\item{r.process}{simulator of the latent state process.
Setting \code{r.process=NULL} removes the latent-state simulator.}

\item{d.process}{optional;
specification of the probability density evaluation function of the unobserved state process.
Setting \code{dprocess=NULL} removes the latent-state density evaluator.}

\item{r.measure}{simulator of the measurement model, specified either as a C snippet, an \R function, or the name of a pre-compiled native routine available in a dynamically loaded library.
Setting \code{rmeasure=NULL} removes the measurement model simulator.}

\item{d.measure}{evaluator of the measurement model density, specified either as a C snippet, an \R function, or the name of a pre-compiled native routine available in a dynamically loaded library.
Setting \code{dmeasure=NULL} removes the measurement density evaluator.}

\item{r.prior}{optional; prior distribution sampler, specified either as a C snippet, an \R function, or the name of a pre-compiled native routine available in a dynamically loaded library.
Setting \code{rprior=NULL} removes the prior distribution sampler.}

\item{d.prior}{optional; prior distribution density evaluator, specified either as a C snippet, an \R function, or the name of a pre-compiled native routine available in a dynamically loaded library.
Setting \code{dprior=NULL} resets the prior distribution to its default, which is a flat improper prior.}

\item{partrans}{optional parameter transformations, constructed using.


If a covariate table is supplied, then the value of each of the covariates is interpolated as needed.
The resulting interpolated values are made available to the appropriate basic components.
See the documentation for \code{\link{covariate_table}} for details.}

\item{verbose}{logical; if \code{TRUE}, diagnostic messages will be printed to the console.}
}
\value{
The \code{contact} constructor function returns an object, call it \code{C}, of class \sQuote{contact}.
\code{C} contains, in addition to the data, any elements of the model that have been specified as arguments to the \code{contact} constructor function.
One can add or modify elements of \code{C} by means of further calls to \code{contact}, using \code{C} as the first argument in such calls.
One can pass \code{C} to most of the \pkg{contact} package methods via their \code{data} argument.
}
\description{
This function constructs a \sQuote{contact} object, encoding a partially-observed process model together with a spatio-temporal date.
As such, it is central to all the package's functionality.
One implements the \acronym{contact type} model by specifying some or all of its \emph{basic components}.
These comprise:
\describe{
\item{kernelmodel,}{which indicate the kenel to sample from;}
\item{r.process,}{the simulator of the unobserved Markov state process;}
\item{d.process,}{the evaluator of the probability density function for transitions of the unobserved Markov state process;}
\item{r.measure,}{the simulator of the observed process, conditional on the unobserved state;}
\item{d.measure,}{the evaluator of the measurement model probability density function;}
\item{r.prior,}{which samples from a prior probability distribution on the parameters;}
\item{d.prior,}{which evaluates the prior probability density function;}
\item{partrans,}{which performs parameter transformations.}
}
}
\details{
Each basic component is supplied via an argument of the same name.
These can be given in the call to \code{contact}, or to many of the package's other functions.
In any case, the effect is the same: to add, remove, or modify the basic component.

Each basic component can be furnished using C snippets, \R functions, or pre-compiled native routine available in user-provided dynamically loaded libraries.
}
\section{Note}{


\strong{It is not typically necessary (or indeed often feasible) to define all of the basic components for any given purpose.
Each \pkg{contact} algorithm makes use of only a subset of these components.
Any algorithm requiring a component that is not present will generate an error letting you know that you have not provided a needed component.
FIXME }
}

\author{
Hola K. Adrakey
}
