#============================================================================
# Define the contact model class
#============================================================================
#' An S4 class to represent a contact model
#' @include builder.R
NULL
setClass(
  "contact",
          slots = c(
            data = "data.frame",
            initializer = "data.frame",
            times= "numeric",
            t0 = 'numeric',
            grid.lines= "data.frame",
            pop.grid = "matrix",
            params = "numeric",
            covar = "matrix",
            tcovar = "numeric",
            grid.size = "numeric",
            age.level = "numeric",
            age.dist  ="numeric",
            t.max = "numeric",
            t.intervention = "numeric",
            kernelmodel = "contact.fun",
            rmodel = "contact.fun",
            r.model = "contact.fun",
            d.model = "contact.fun",
            d.prior = 'contact.fun',
            r.prior = 'contact.fun',
            states = "array",
            nkernls = "numeric",
            has.trans = 'logical',
            from.trans = 'contact.fun',
            to.trans = 'contact.fun',
            zeronames = "character",
            userdata = 'list'
          ),

  prototype = prototype(
    data = as.data.frame(array(data=numeric(0),dim=c(0,0))),
    initializer = as.data.frame(array(data=numeric(0),dim=c(0,0))),
    times = numeric(0),
    t0=numeric(0),
    grid.lines = as.data.frame(array(data=numeric(0),dim=c(0,0))),
    pop.grid = array(data=numeric(0),dim=c(0,0)),
    params =numeric(0),
    grid.size = numeric(0),
    age.level = numeric(0),
    age.dist = numeric(0),
    t.max = numeric(0),
    t.intervention = numeric(0),
    covar=array(data=numeric(0),dim=c(0,0)),
    tcovar=numeric(0),
    kernelmodel = contact.fun(slotname="kernelmodel"),
    rmodel = contact.fun(slotname="rmodel"),
    r.model = contact.fun(slotname="r.model"),
    d.model = contact.fun(slotname="d.model"),
    r.prior = contact.fun(slotname="r.prior"),
    d.prior = contact.fun(slotname="d.prior"),
    states = array(data=numeric(0),dim=c(0,0)),
    nkernls = numeric(0),
    from.trans=contact.fun(slotname="fromEstimationScale"),
    to.trans=contact.fun(slotname="toEstimationScale"),
    userdata=list()

  ),
  validity= function (object){
    retval <- character(0)
    if (length(object@data)<1)
      retval <- append(retval,paste(sQuote("data"),"is a required argument"))
    if (length(object@times)<1)
      append(retval,paste(sQuote("data"),"is a required argument"))
    if (length(object@grid.lines)<1)
      append(retval,paste(sQuote("grid.lines"),"is a required argument"))
    if (length(object@pop.grid)<1)
      append(retval,paste(sQuote("pop.grid"),"is a required argument"))
    if ((ncol(object@pop.grid)+nrow(object@pop.grid) + 2)!=nrow(object@grid.lines))
      retval <- append(retval,paste("the number of colums + the number of rows +2 should match the number of lines"))
    if (!is.numeric(object@params) || (length(object@params)>0 && is.null(names(object@params))))
      retval <- append(retval,paste(sQuote("params"),"must be a named numeric vector"))
    # if (ncol(object@data)!=length(object@times))
    #   retval <- append(retval,paste("the length of",sQuote("times"),"should match the number of observations"))
    if (length(object@t0)<1)
      retval <- append(retval,paste(sQuote("t0"),"is a required argument"))
    if (length(object@t0)>1)
      retval <- append(retval,paste(sQuote("t0"),"must be a single number"))
    # if (object@t0 > object@times[1])
    #   retval <- append(retval,paste("the zero-time",sQuote("t0"),
    #                                 "must occur no later than the first observation"))
    if(length(object@age.level)<1)
      append(retval,paste(sQuote("age.level"),"is a required argument"))
    if(length(object@age.dist)<1)
      append(retval,paste(sQuote("age.dist"),"is a required argument"))
    if(length(object@grid.size)<1)
      append(retval,paste(sQuote("grid.size"),"is a required argument"))
    if(length(object@age.dist)!=length(object@age.level))
      retval <- append(retval,paste("the length of",sQuote("age.level"),"should match the number of age group"))
    if (length(object@kernelmodel)==0)
      retval <- append(retval,paste(sQuote("kernelmodel"),"is requiered"))
    if (length(object@rmodel)==0)
      retval <- append(retval,paste(sQuote("rmodel"),"is requiered"))
    # if(object@nkernls<=1 && missing(object@kernel.model1))
    #   retval <- append(retval,paste(sQuote("kernel.model1"),"is a required argument"))
    # if(object@nkernls==2 && (missing(object@kernel.model1) || missing(object@kernel.model2)) )
    #   retval <- append(retval,paste(sQuote("kernel.model1 and kernel.mdoel2"),"are required arguments"))
    if (length(object@tcovar)!=nrow(object@covar)) {
      retval <- append(
        retval,
        paste(
          "the length of",sQuote("tcovar"),
          "should match the number of rows of",sQuote("covar")
        )
      )
    }
    if (!is.numeric(object@tcovar))
      retval <- append(
        retval,
        paste(
          sQuote("tcovar"),
          "must either be a numeric vector or must name a numeric vector in the data frame",
          sQuote("covar")
        )
      )
    if (length(retval)==0) TRUE else retval
  }

)

