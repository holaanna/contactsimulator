#============================================================================
# A class for functions that may be defined in R, using pre-written native routines of Cpp functions
#============================================================================
#' A class for functions that may be defined in R, using pre-written native routines of Cpp functions

contactfunmode <- list(undef=-1L, Rfun=0L, native=1L, regNative=2L)

#' @include cppsnippet.R rsnippet.R generics.R
NULL

setClass(
  "contact.fun",
  slots = c(
    native.fun = "character",
    PACKAGE = "character",
    mode = "integer",
    address = "externalptr",
    obsnames = 'character',
    statenames = 'character',
    paramnames = 'character',
    covarnames = 'character',
    purpose = "character"
  ),
  prototype = prototype(
    native.fun=character(0),
    PACKAGE = character(0),
    mode =  contactfunmode$undef,
    obsnames = character(0),
    statenames = character(0),
    paramnames = character(0),
    covarnames = character(0),
    purpose = "a needed function"
  )

)

setGeneric("contact.fun", function(f, ...) standardGeneric("contact.fun"))

setMethod(
  "contact.fun",
  signature = signature(f="missing"),
  definition = function(f, slotname = NULL, ...){
    new("contact.fun", purpose=as.character(slotname))
  }
  )

setMethod(
  "contact.fun",
  signature = signature(f="NULL"),
  definition = function(f, slotname = NULL, ...){
    new("contact.fun", purpose=as.character(slotname))
  }
)

setMethod(
  "contact.fun",
  signature = signature(f="ANY"),
  definition=function (f, slotname = NULL, ...) {
    stop("bad option for ",sQuote(slotname)," argument",call.=FALSE)
  }
)

setMethod(
  "contact.fun",
  signature=signature(f="Rsnippet"),
  definition=function (f, PACKAGE, slotname = NULL, libname = NULL,
                       obsnames = character(0), statenames = character(0),
                       paramnames = character(0), covarnames = character(0), ...) {
    if (is.null(slotname))
      stop("in ",sQuote("contact.fun"),": unspecified ",
           sQuote("slotname"),call.=FALSE)
    if (is.null(libname))
      stop("in ",sQuote("contact.fun"),": unspecified ",
           sQuote("libname"),call.=FALSE)
    slotname <- as.character(slotname)
    libname <- as.character(libname)
    new(
      "contact.fun",
      native.fun=render(slotname,name=libname),
      PACKAGE=as.character(libname),
      mode=contactfunmode$Rfun,
      obsnames=obsnames,
      statenames=statenames,
      paramnames=paramnames,
      covarnames=covarnames,
      purpose=slotname
    )
  }
)

setMethod(
  "contact.fun",
  signature=signature(f="Cppsnippet"),
  definition=function (f, PACKAGE, slotname = NULL, libname = NULL,
                       obsnames = character(0), statenames = character(0),
                       paramnames = character(0), covarnames = character(0), ...) {
    if (is.null(slotname))
      stop("in ",sQuote("contact.fun"),": unspecified ",
           sQuote("slotname"),call.=FALSE)
    if (is.null(libname))
      stop("in ",sQuote("contact.fun"),": unspecified ",
           sQuote("libname"),call.=FALSE)
    slotname <- as.character(slotname)
    libname <- as.character(libname)
    new(
      "contact.fun",
      native.fun=render(slotname,name=libname),
      PACKAGE=as.character(libname),
      mode=contactfunmode$regNative,
      obsnames=obsnames,
      statenames=statenames,
      paramnames=paramnames,
      covarnames=covarnames,
      purpose=slotname
    )
  }
)



setMethod(
  "contact.fun",
  signature = signature(f="character"),
  definition = function(f, PACKAGE = NULL,
                        slotname = NULL, ...){
    new(
      "contact.fun",
      native.fun = f,
      PACKAGE = as.character(PACKAGE),
      mode = contactfunmode$native,
      purpose=as.character(slotname)
    )
  }
)

setMethod(
  "contact.fun",
  signature=signature(f="contact.fun"),
  definition=function (f, ...) f
)

setMethod(
  'show',
  signature=signature('contact.fun'),
  definition=function (object) {
    mode <- object@mode
    if (mode==contactfunmode$Rfun) { # R function
      show(object@R.fun)
    } else if (mode==contactfunmode$native) { # user supplied native code
      cat("native function ",sQuote(object@native.fun),sep="")
      if (length(object@PACKAGE)>0)
        cat(", dynamically loaded from ",sQuote(object@PACKAGE),sep="")
      cat ("\n")
    } else if (mode==contactfunmode$regNative) { # built from C snippets
      cat("native function ",sQuote(object@native.fun),sep="")
      if (length(object@PACKAGE)>0)
        cat(", defined by a Cppsnippet",sep="")
      cat ("\n")
    } else {
      cat("not specified\n")
    }
  }
)

setMethod(
  'print',
  'contact.fun',
  function (x, ...) show(x)
)
