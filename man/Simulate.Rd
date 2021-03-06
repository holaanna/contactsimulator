% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nsimulations.R
\name{Simulate}
\alias{Simulate}
\title{Generate a certain number of temporal observation from the process}
\usage{
Simulate(param, grid_lines, pop_grid, grid_size = 5000,
  age_level = c(1, 1), age_dist = c(1, 0), m_start = 1,
  t_max = 118, t_intervention = 365, EI_model = 1, kern_model = 4,
  times, obs, nsim)
}
\arguments{
\item{param}{Indicating a data frame containing a vector of parameters including:
\describe{
  \item{epsion}{The primary infection rate. See \code{\link{func_time_beta}}}
  \item{beta_0}{Baseline or average transmission rate. See \code{\link{func_time_beta}}}
  \item{beta_1}{Amplitude of the seasonality. See \code{\link{func_time_beta}}}
  \item{alpha1,alpha2}{The dispersal kernel parameters.}
  \item{mu_lat,var_lat}{mean and variance of the latent period. See \code{\link{E_to_I}} for details.}
  \item{t0}{Time at which the primary source became active}.
  \item{omega}{Period of the forcing. See \code{\link{func_time_beta}}}
  \item{gama}{The mean proportion of short range dispersal events.}.
 }}

\item{pop_grid}{Population density of the grid a case resides. This is filled from bottom to top, then left to right.}

\item{grid_size}{Grid resolution
//@inheritParams circle_line_intersections}

\item{age_level}{Vectors of age level and the propportion of each age group respectively. See details.}

\item{age_dist}{Vectors of age level and the propportion of each age group respectively. See details.}

\item{m_start}{The size of initial cases. Default is 1.}

\item{t_max}{Final observation time.}

\item{t_intervention}{Start of the intervention if any.}

\item{EI_model}{Take integer values to specify the type of model used for the latent period. See \code{\link{E_to_I}}}

\item{kern_model}{Take integer values to specify the type of dispersal kernel used. See \code{\link{Samp_dis}}}

\item{obs}{The obervation requiered from during the period time}

\item{nsim}{The total size of the simulation requierd.}

\item{time}{The sequence of time at which the observation will occure.}
}
\value{
A data frame with the time, observation made and the nth simulation.
}
\description{
Temporal oservation using the contact type model.
}
\details{
\code{Simulate} Count the number of observation during a certain period.
}
