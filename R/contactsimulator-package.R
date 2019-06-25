##' Inference for spatio-temporal contact type models
##'
##' The \pkg{contactsimulator} package provides facilities for inference on spatial-temporal
##' data using contact type distribution models.
##' These models are also known as state-space models, hidden Markov models, or
##' nonlinear stochastic dynamical systems.  One can use \pkg{contactsimulator} to fit
##' nonlinear, non-Gaussian dynamic models to time-series data.  The package is
##' both a set of tools for data analysis and a platform upon which statistical
##' inference methods for contacts models can be implemented.
##'
##' @name contactsimulator-package
##' @aliases contactsimulator,package contactsimulator-package
##' @docType package
##' @author Hola Adrakey
##' @family information on model implementation
##' @family contactsimulator parameter estimation methods
##' @family elementary contact methods
##'
##' @section Data analysis using \pkg{contactsimulator}:
##' \pkg{contactsimulator} provides algorithms for
##' \enumerate{
##' \item simulation of stochastic
##' dynamical systems; see \code{\link[=simulate-contactsimulator]{simulate}}
##' }
##' The package
##' also provides various tools for plotting and extracting information on
##' models and data.


#' @importFrom Rdpack reprompt
#'
#' @useDynLib contactsimulator
#' @importFrom Rcpp sourceCpp
NULL
