## This file defines "contact model", the basic constructor of the contact class
#' @include contact_class.R
NULL
contact.internal <- function(data, initializer, times, t0, grid.lines,
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
  if (missing(grid.lines)) stop(sQuote("grid.lines")," is a requiered argument",call. = FALSE)
  if (missing(pop.grid)) stop(sQuote("pop.grid")," is a requiered argument",call. = FALSE)
  if (missing(grid.size)) stop(sQuote("grid.size")," is a requiered argument",call. = FALSE)
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
  if (missing(kernelmodel)) kernelmodel <- NULL
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
  rmodel <- contact.fun(
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
  r.model <- contact.fun(
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
  d.model <- contact.fun(
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
  r.prior <- contact.fun(
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
  d.prior <- contact.fun(
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
  kernelmodel <- contact.fun(
    f=kernelmodel,
    PACKAGE=PACKAGE,
    libname=libname,
    proto=quote(kernelmodel(r,params,...)),
    statenames=statenames,
    paramnames=paramnames,
    obsnames=obsnames,
    covarnames=covarnames,
    slotname="kernelmodel"
  )



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

contact <- function (data, initializer, times, t0, ..., grid.lines,
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

    #if()

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
