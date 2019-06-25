## This file defines "contact model", the basic constructor of the contact class
#' Constructor of the basic contact object
#'
#' This function constructs a \sQuote{contact} object, encoding a partially-observed process model together with a spatio-temporal date.
#' As such, it is central to all the package's functionality.
#' One implements the \acronym{contact type} model by specifying some or all of its \emph{basic components}.
#' These comprise:
#' \describe{
#' \item{kernelmodel,}{which indicate the kenel to sample from;}
#' \item{r.process,}{the simulator of the unobserved Markov state process;}
#' \item{d.process,}{the evaluator of the probability density function for transitions of the unobserved Markov state process;}
#' \item{r.measure,}{the simulator of the observed process, conditional on the unobserved state;}
#' \item{d.measure,}{the evaluator of the measurement model probability density function;}
#' \item{r.prior,}{which samples from a prior probability distribution on the parameters;}
#' \item{d.prior,}{which evaluates the prior probability density function;}
#' \item{partrans,}{which performs parameter transformations.}
#' }
#'
#' Each basic component is supplied via an argument of the same name.
#' These can be given in the call to \code{contact}, or to many of the package's other functions.
#' In any case, the effect is the same: to add, remove, or modify the basic component.
#'
#' Each basic component can be furnished using C snippets, \R functions, or pre-compiled native routine available in user-provided dynamically loaded libraries.
#'
#' @name contact
#' @rdname contact
#' @include contact_class.R contact_fun.R cppsnippet.R  builder.R
#' @importFrom stats setNames
#'
#'
#' @param data either a data frame holding the time series data,
#' or an object of class \sQuote{contact},
#' i.e., the output of another \pkg{contact} calculation.
#'
#' @param rast an object of class rasterlayer containing the landscape informationn.
#'
#' @param initializer the initial state contaning
#'       \describe{
#'         \item{row}{The row at which the cell containing the host lies in}
#'         \item{col}{The column at which the cell containing the host lies in}
#'         \item{typ}{The type of the host (e.g. backyard 1 /plantation 0 hosts)}
#'         \item{x}{x-coordinate}
#'         \item{y}{y-coordinate}
#'         \item{t_e}{exposure time}
#'         \item{t_i}{infection time}
#'        }
#' @inheritParams Simulate_contact_model
#' @param times the times at which observations are made.
#' \code{times} must indicate the column of observation times by name or index.
#' The time vector must be numeric and non-decreasing.
#' Internally, \code{data} will be internally coerced to an array with storage-mode \code{double}.
#'
#' @param t0 The zero-time, i.e., the time of the initial state.
#' This must be no later than the time of the first observation, i.e., \code{t0 <= times[1]}.
#'
#' @param grid.lines A 6 columns data frame with columns names as coor_x_1, coor_y_1, coor_x_2, coor_y_2, orient_line.
#' \describe{
#'     \item{coor_x_1, coor_y_1}{Coordinates of the left end point of the grid line }
#'     \item{coor_x_2, coor_y_2}{Coordinates of the right end point of the grid line }
#'     \item{orient_line}{Line orientation}
#'     \enumerate{
#'        \item indicates horizontal orientation
#'        \item indicates vetical orientation
#'     }
#'     \item{k_line}{Line numbering: bottom to top, then left to right}
#' }
#'
#' @param pop.grid  Population density of the grid a case resides. This is filled from bottom to top, then left to right.
#'
#' @param age.level,age.dist Vectors of age level and the propportion of each age group respectively. See details.
#'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          max  Final observation time.
#' @param t_intervention Start of the intervention if any.
#'
#' @param EI_model Take integer values to specify the type of model used for the latent period. See \code{\link{E_to_I}}
#' @param kern.  model Take integer values to specify the type of dispersal kernel used. See \code{\link{Samp_dis}}

#'
#' @param r.process simulator of the latent state process.
#' Setting \code{r.process=NULL} removes the latent-state simulator.
#'
#' @param d.process optional;
#' specification of the probability density evaluation function of the unobserved state process.
#' Setting \code{dprocess=NULL} removes the latent-state density evaluator.
#'
#' @param r.measure simulator of the measurement model, specified either as a C snippet, an \R function, or the name of a pre-compiled native routine available in a dynamically loaded library.
#' Setting \code{rmeasure=NULL} removes the measurement model simulator.
#'
#' @param d.measure evaluator of the measurement model density, specified either as a C snippet, an \R function, or the name of a pre-compiled native routine available in a dynamically loaded library.
#' Setting \code{dmeasure=NULL} removes the measurement density evaluator.
#'
#'
#' @param r.prior optional; prior distribution sampler, specified either as a C snippet, an \R function, or the name of a pre-compiled native routine available in a dynamically loaded library.
#' Setting \code{rprior=NULL} removes the prior distribution sampler.
#'
#' @param d.prior optional; prior distribution density evaluator, specified either as a C snippet, an \R function, or the name of a pre-compiled native routine available in a dynamically loaded library.
#' Setting \code{dprior=NULL} resets the prior distribution to its default, which is a flat improper prior.
#'
#' @param partrans optional parameter transformations, constructed using.
#'
#'
#' If a covariate table is supplied, then the value of each of the covariates is interpolated as needed.
#' The resulting interpolated values are made available to the appropriate basic components.
#' See the documentation for \code{\link{covariate_table}} for details.
#'
#' @param params optional; named numeric vector of parameters.
#' This will be coerced internally to storage mode \code{double}.
#'
#' @param obsnames optional character vector;
#' names of the observables.
#' It is not usually necessary to specify \code{obsnames} since, by default,
#' these are read from the names of the data variables.
#'
#' @param statenames optional character vector;
#' names of the latent state variables.
#' It is typically only necessary to supply \code{statenames} when C snippets are in use.
#'
#' @param paramnames optional character vector;
#' names of model parameters.
#' It is typically only necessary to supply \code{paramnames} when C snippets are in use.
#'
#' @param covarnames optional character vector;
#' names of the covariates.
#' It is not usually necessary to specify \code{covarnames} since, by default,
#' these are read from the names of the covariates.
#'
#' @param \dots additional arguments supply new or modify existing model characteristics or components.
#' See \code{\link{contact}} for a full list of recognized arguments.
#'
#' When named arguments not recognized by \code{\link{contact}} are provided, these are made available to all basic components via the so-called \dfn{userdata} facility.
#' This allows the user to pass information to the basic components outside of the usual routes of covariates (\code{covar}) and model parameters (\code{params}).
#'
#' @param verbose logical; if \code{TRUE}, diagnostic messages will be printed to the console.
#'
#' @return
#' The \code{contact} constructor function returns an object, call it \code{C}, of class \sQuote{contact}.
#' \code{C} contains, in addition to the data, any elements of the model that have been specified as arguments to the \code{contact} constructor function.
#' One can add or modify elements of \code{C} by means of further calls to \code{contact}, using \code{C} as the first argument in such calls.
#' One can pass \code{C} to most of the \pkg{contact} package methods via their \code{data} argument.
#'
#' @section Note:
#'
#' \strong{It is not typically necessary (or indeed often feasible) to define all of the basic components for any given purpose.
#' Each \pkg{contact} algorithm makes use of only a subset of these components.
#' Any algorithm requiring a component that is not present will generate an error letting you know that you have not provided a needed component.
#' FIXME }
#'
#' @author Hola K. Adrakey
#'
NULL

#' @rdname contact
#' @export
contact <- function (data, rast, initializer, times, t0, ..., grid.lines,
                     pop.grid, params, grid.size, age.level, nkernls,
                     age.dist, t.max, t.intervention, kernelmodel,
                     obsnames, statenames, paramnames, covarnames, zeronames,
                     rmodel, r.model, PACKAGE,
                     d.model){

  ep <- paste0("in ", sQuote("contact"),": ")

  if (missing(data))
    stop(ep,sQuote("data")," is a required argument",call.=FALSE)

  if (is(data,"contact")){
    ## data is a contact object:
    ## extract missing arguments from it

    if (nargs()==1) return(data)

    if(missing(times)) times <- data@times
    if (missing(t0))   t0 <- data@t0
    if (missing(grid.lines)) grid.lines <- data@grid.lines
    if (missing(pop.grid))  pop.grid <- data@pop.grid
    if (missing(params))  params <- data@params
    if (missing(grid.size))  grid.size <- data@grid.size
    if (missing(age.level))  age.level <- data@age.level
    if (missing(age.dist))  age.dist <- data@age.dist
    if(missing(t.max)) t.max <- data@t.max
    if (missing(t.intervention))   t.intervention <- data@t.intervention
    if (missing(kernelmodel)) kernelmodel <- data@kernelmodel
    if (missing(rmodel))  rmodel <- data@rmodel
    if (missing(r.model))  r.model <- data@r.model
    if (missing(d.model))  d.model <- data@d.model
    if (missing(nkernls)) nkernls <- data@nkernls
    if (missing(zeronames)) zeronames <- data@zeronames
    if (missing(obsnames)) obsnames <- character(0)
    if (missing(statenames)) statenames <- character(0)
    if (missing(paramnames)) paramnames <- character(0)
    if (missing(covarnames)) covarnames <- character(0)
    if (missing(PACKAGE)) PACKAGE <- character(0)

    tryCatch(
      contact.internal(
        data=data@data,
        initializer = initializer,
        times = times,
        t0 = t0,
        grid.lines = grid.lines,
        pop.grid = pop.grid,
        params = params,
        grid.size = grid.size,
        age.level = age.level,
        age.dist = age.dist,
        t.max = t.max,
        t.intervention = t.intervention,
        kernelmodel = kernelmodel,
        rmodel = rmodel,
        r.model = r.model,
        PACKAGE = PACKAGE,
        d.model = d.model,
        nkernls = nkernls,
        userdata=data@userdata,
        statenames=statenames,
        paramnames=paramnames,
        covarnames=covarnames,
        obsnames=obsnames,
        zeronames=zeronames,
        ...
      ),
      error = function (e) {
        stop(ep,conditionMessage(e),call.=FALSE)
      }
    )
  }
  else{
  ## construct a contact object de nouveau

    if (is.data.frame(data)){

      ## data is a data frame
      if ((is.numeric(times) && (times<1 || times>ncol(data) ||times!=as.integer(times))) ||
          (is.character(times) && (!(times%in%names(data)))) ||
          (!is.numeric(times) && !is.character(times)) ||
          length(times)!=1) {
        stop(ep,"when ",sQuote("data")," is a data frame, ",sQuote("times"),
             " must identify a single column of ",sQuote("data"),
             " either by name or by index.",call.=FALSE)
      }

      if (is.numeric(times)) {
        tpos <- as.integer(times)
      } else if (is.character(times)) {
        tpos <- match(times,names(data))
      }
      times <- data[[tpos]]
      data <- do.call(rbind, lapply(data[-tpos], as.numeric))
      data <- as.data.frame(data)

    }
    else {
      stop(ep,sQuote("data"),
           " must be a data frame or an object of class ",sQuote("contact"),
           call.=FALSE)
    }


    if (missing(rast) | all(is.na(raster::values(rast)))) {
      if (missing(grid.lines)) stop(ep,"when ", sQuote("rast"), " is missing or null values",sQuote("grid.lines")," must be specified",call. = FALSE)
      if (missing(pop.grid)) stop(ep,"when ", sQuote("rast"), " is missing or null values",sQuote("pop.grid")," must be specified",call. = FALSE)
      if (missing(grid.size)) stop(ep,"when ", sQuote("rast"), " is missing or null values",sQuote("grid.size")," must be specified",call. = FALSE)
    }
    else{
      size<- raster::res(rast)[1]
      n_row_grid <- nrow_grid <- raster::nrow(rast)
      n_col_grid <- ncol_grid <- raster::ncol(rast)
      grid.size <- raster::res(rast)[1]     # Resolution

      n_line <- (nrow_grid+1) + (ncol_grid +1)  # Number of grid  lines

      x_min <- raster::xmin(rast)  # min max of the bounding box
      x_max <- raster::xmax(rast)

      y_min <- raster::ymin(rast)
      y_max <- raster::ymax(rast)

      pop_per_grid <- round(raster::values(rast)*size^2)
      pop_per_grid[is.na(pop_per_grid)] <- 0
      mat <- matrix(pop_per_grid,nrow  <-  nrow_grid, byrow = TRUE)
      pop.grid <- flipdim(mat)     # population per grid

      # Structure of the grid
      x <- seq(x_min,x_max,grid.size)
      y <- seq(y_min,y_max,grid.size)

      grid.lines <- array(0,c(n_line,6))
      for(i in 1:n_line){
        if(i<=(nrow_grid +1)){
          grid.lines[i,] <- c(i,1,x[1],y[i],x[length(x)],y[i])
        }
        else{
          grid.lines[i,] <- c(i,2,x[i-length(y)],y[1],x[i-length(y)],y[length(y)])
        }
      }

      grid.lines <- as.data.frame(grid.lines)
      colnames(grid.lines)<- c("indx","orient_line","coor_x_1","coor_y_1","coor_x_2","coor_y_2")
    }

   if(missing(t.max)) t.max <- max(times) +1


    #if(all(raster::values(rast)))

    tryCatch(
      contact.internal(
        data=data,
        initializer = initializer,
        times = times,
        t0 = t0,
        grid.lines = grid.lines,
        pop.grid = pop.grid,
        params = params,
        grid.size = grid.size,
        age.level = age.level,
        age.dist = age.dist,
        t.max = t.max,
        t.intervention = t.intervention,
        kernelmodel = kernelmodel,
        rmodel = rmodel,
        r.model = r.model,
        PACKAGE = PACKAGE,
        d.model = d.model,
        nkernls = nkernls,
        statenames=statenames,
        paramnames=paramnames,
        covarnames=covarnames,
        obsnames=obsnames,
        zeronames=zeronames,
        ...
      ),
      error = function (e) {
        stop(ep,conditionMessage(e),call.=FALSE)
      }
    )

  }

}




contact.internal <- function(data, rast, initializer, times, t0, grid.lines,
                             pop.grid, params, grid.size, age.level,
                             age.dist, t.max, t.intervention, kernelmodel,
                             rmodel, r.model, PACKAGE, r.prior,zeronames,
                             obsnames, statenames, paramnames, covarnames,
                             d.prior, d.model, fromEstimationScale, toEstimationScale,
                             covar, tcovar, nkernls, userdata,...,
                             verbose = getOption("verbose",FALSE)){


  ep <- paste0("in ",sQuote("contact"),": ")

  ##preliminary error checking
  if (missing(times)) stop(sQuote("times")," is a requiered argument",call. = FALSE)
  if (missing(rast)) {
    if (missing(grid.lines)) stop(sQuote("grid.lines")," is a requiered argument",call. = FALSE)
    if (missing(pop.grid)) stop(sQuote("pop.grid")," is a requiered argument",call. = FALSE)
    if (missing(grid.size)) stop(sQuote("grid.size")," is a requiered argument",call. = FALSE)
  }

  if (missing(t0))  stop(sQuote("t0")," is a requiered argument",call. = FALSE)
  if (missing(params)) params <- numeric(0)
  if (missing(userdata)) userdata <- list()
  added.userdata <- list(...)

  if (length(added.userdata)>0) {
    message("In ", sQuote("contact"),
            ": the following unrecongnized argument(s) ",
            "will be sotored for use by user-defined functions: ",
            paste(sQuote(names(added.userdata)),collapse = ","))
    userdata[names(added.userdata)] <- added.userdata
  }

  if(missing(age.level)) age.level<- c(1,1)
  if(missing(age.dist)) age.dist <- c(1,0)
  if(missing(t.intervention)) t.intervention <- 10000000
  ## name of shared object library
  if (missing(PACKAGE)) PACKAGE <- NULL
  PACKAGE <- as.character(PACKAGE)

  ## deal with mising components
  if (missing(rmodel)) rmodel <- NULL
  if (missing(r.model)) r.model <- NULL
  if (missing(d.model)) d.model <- NULL
  if (missing(r.prior)) r.prior <- NULL
  if (missing(d.prior)) d.prior <- NULL
  if (missing(nkernls)) nkernls <- 1
  if (missing(kernelmodel)) kernelmodel <- 5  # Exponential
  if (missing(initializer))  initializer <- as.data.frame(array(data=0,dim = c(0,0)))
  if (missing(fromEstimationScale)) fromEstimationScale <- NULL
  if (missing(toEstimationScale)) toEstimationScale <- NULL

  ## defaults for names of states, parameters, and accumulator variables
  ## defaults for names of states, parameters, and accumulator variables
  if (missing(statenames)) statenames <- character(0)
  if (missing(paramnames)) paramnames <- character(0)
  if (missing(zeronames)) zeronames <- character(0)
  statenames <- as.character(statenames)
  paramnames <- as.character(paramnames)
  zeronames <- as.character(zeronames)

  ## check for duplicate names
  if (anyDuplicated(statenames)) {
    stop("all ",sQuote("statenames")," must be unique", call.=FALSE)
  }
  if (anyDuplicated(paramnames)) {
    stop("all ",sQuote("paramnames")," must be unique", call.=FALSE)
  }
  if (anyDuplicated(zeronames)) {
    stop("all ",sQuote("zeronames")," must be unique", call.=FALSE)
  }

  ## check the parameters and force them to be double-precision
  if (length(params)>0) {
    if (is.null(names(params)) || !is.numeric(params))
      stop(sQuote("params")," must be a named numeric vector",
           call.=FALSE)
  }
  storage.mode(params) <- 'double'

  ## check the data and store it as double-precision matrix
  if (!is.data.frame(data))
    stop(sQuote("data")," must be data frame",call.=FALSE)
  storage.mode(data) <- 'double'
  if (missing(obsnames) || length(obsnames)==0) obsnames <- rownames(data)
  obsnames <- as.character(obsnames)

  ## check times
  if (!is.numeric(times) || any(is.na(times)) || !all(diff(times)>0))
    stop(sQuote("times")," must be an increasing numeric vector",call.=FALSE)
  storage.mode(times) <- 'double'

  ## check initializer
  if (!is.data.frame(initializer))
    stop(sQuote("initializer")," must be data frame",call.=FALSE)
  storage.mode(initializer) <- 'double'

  ## check t0
  if (!is.numeric(t0) || length(t0) > 1)
    stop("the zero-time ",sQuote("t0")," must be a single number",call.=FALSE)
  storage.mode(t0) <- 'double'

  ## check the kernels
  if(nkernls==1 & is.null(kernelmodel)){
    stop("the kernel model ",sQuote("kernelmodel")," must be specified",call.=FALSE)
  }



  ## check the statenames

  ## check and arrange covariates
  if (missing(covar)) {
    covar <- matrix(data=0,nrow=0,ncol=0)
    tcovar <- numeric(0)
  } else if (missing(tcovar)) {
    stop("if ",sQuote("covar")," is supplied, ",
         sQuote("tcovar")," must also be supplied",call.=FALSE)
  } else if (is.data.frame(covar)) {
    if ((is.numeric(tcovar) && (tcovar<1 || tcovar>length(covar))) ||
        (is.character(tcovar) && (!(tcovar%in%names(covar)))) ||
        (!is.numeric(tcovar) && !is.character(tcovar))) {
      stop("if ",sQuote("covar")," is a data frame, ",
           sQuote("tcovar")," should indicate the time variable",call.=FALSE)
    } else if (is.numeric(tcovar)) {
      tpos <- tcovar
      tcovar <- covar[[tpos]]
      covar <- as.matrix(covar[-tpos])
    } else if (is.character(tcovar)) {
      tpos <- match(tcovar,names(covar))
      tcovar <- covar[[tpos]]
      covar <- as.matrix(covar[-tpos])
    }
  } else {
    covar <- as.matrix(covar)
  }
  if (missing(covarnames) || length(covarnames)==0)
    covarnames <- as.character(colnames(covar))
  if (!all(covarnames%in%colnames(covar))) {
    missing <- covarnames[!(covarnames%in%colnames(covar))]
    stop("covariate(s) ",
         paste(sapply(missing,sQuote),collapse=","),
         " are not among the columns of ",sQuote("covar"),call.=FALSE)
  }
  storage.mode(tcovar) <- "double"
  storage.mode(covar) <- "double"


  snips <- list()
  if (is(rmodel,"Cppsnippet"))
    snips <- c(snips,rmodel=rmodel@text)

  if (is(kernelmodel,"Cppsnippet"))
    snips <- c(snips,kernelmodel=kernelmodel@text)

  if (is(r.prior,"Cppsnippet"))
    snips <- c(snips,r.prior=r.prior@text)

  if (is(d.prior,"Cppsnippet"))
    snips <- c(snips,d.prior=d.prior@text)

  rsnips <- list()  # class Rsinps
  if (is(rmodel,"Rsnippet"))
    rsnips <- c(rsnips,rmodel=rmodel@text)

  if (is(kernelmodel,"Rsnippet"))
    rsnips <- c(rsnips,kernelmodel=kernelmodel@text)

  if (is(r.prior,"Rsnippet"))
    rsnips <- c(rsnips,r.prior=r.prior@text)

  if (is(d.prior,"Rsnippet"))
    rsnips <- c(rsnips,d.prior=d.prior@text)


  # build of modules. Note that we consider them here as libraries
  if (length(snips)>0) {
    libname <- tryCatch(
      do.call(
        contactCBuilder,
        c(
          list(
            obsnames=obsnames,
            statenames=statenames,
            paramnames=paramnames,
            covarnames=covarnames,
            verbose=verbose ,
            snips=snips
          )
        )
      ),
      error = function (e) {
        stop("error in building shared-object library from C++ snippets: ",
             conditionMessage(e),call.=FALSE)
      }
    )
    libname <- libname@.xData$moduleName
  }
  # else {
  #   libname <- ''
  # }

  if (length(rsnips)>0) {
    libname <- tryCatch(
      do.call(
        contactCBuilder,
        c(
          list(
            obsnames=obsnames,
            statenames=statenames,
            paramnames=paramnames,
            covarnames=covarnames,
            verbose=verbose,
            snips=snips
          )

        )
      ),
      error = function (e) {
        stop("error in building shared-object library from C++ snippets: ",
             conditionMessage(e),call.=FALSE)
      }
    )
    libname <- libname$name
  }
  # else {
  #   libname <- ''
  # }

  ## handle rmodel
  rmodel <- contact_fun(
    f=rmodel,
    PACKAGE=PACKAGE,
    libname=libname,
    proto=quote(rmodel(x,t,params,...)),
    statenames=statenames,
    paramnames=paramnames,
    obsnames=obsnames,
    covarnames=covarnames,
    slotname="rmodel"
  )

  ## handle r.model
  r.model <- contact_fun(
    f=r.model,
    PACKAGE=PACKAGE,
    libname=libname,
    proto=quote(r.model(x,t,params,...)),
    statenames=statenames,
    paramnames=paramnames,
    obsnames=obsnames,
    covarnames=covarnames,
    slotname="r.model"
  )

  ## handle d.model
  d.model <- contact_fun(
    f=d.model,
    PACKAGE=PACKAGE,
    libname=libname,
    proto=quote(d.model(x,t1,t2,params,...)),
    statenames=statenames,
    paramnames=paramnames,
    obsnames=obsnames,
    covarnames=covarnames,
    slotname="d.model"
  )


  ## handle r.model
  r.prior <- contact_fun(
    f=r.prior,
    PACKAGE=PACKAGE,
    libname=libname,
    proto=quote(r.model(params,...)),
    statenames=statenames,
    paramnames=paramnames,
    obsnames=obsnames,
    covarnames=covarnames,
    slotname="r.prior"
  )

  ## handle d.model
  d.prior <- contact_fun(
    f=d.prior,
    PACKAGE=PACKAGE,
    libname=libname,
    proto=quote(d.prior(params,...)),
    statenames=statenames,
    paramnames=paramnames,
    obsnames=obsnames,
    covarnames=covarnames,
    slotname="d.prior"
  )

  ## handle kernel model
  # kernelmodel <- contact.fun(
  #   f=kernelmodel,
  #   PACKAGE=PACKAGE,
  #   libname=libname,
  #   proto=quote(kernelmodel(r,params,...)),
  #   statenames=statenames,
  #   paramnames=paramnames,
  #   obsnames=obsnames,
  #   covarnames=covarnames,
  #   slotname="kernelmodel"
  # )



  new(
    'contact',
    data = data,
    initializer = initializer,
    times = times,
    t0 = t0,
    grid.lines = grid.lines,
    pop.grid = pop.grid,
    params = params,
    grid.size = grid.size,
    age.level = age.level,
    age.dist = age.dist,
    t.max = t.max,
    t.intervention = t.intervention,
    kernelmodel = kernelmodel,
    rmodel = rmodel,
    r.model = r.model,
    r.prior = r.prior,
    d.prior = d.prior,
    d.model = d.model,
    nkernls = nkernls,
    zeronames = zeronames,
    userdata = userdata
  )
}

