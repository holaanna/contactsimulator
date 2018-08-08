#Function to create the cpp module
#' @include contact_fun.R
NULL

contactCBuilder <- function (name = NULL, snips, statenames, paramnames, covarnames, obsnames, params,
                          verbose = getOption("verbose",FALSE)){


  if (!is.null(name)) name <- cleanForC(name)
  id <- randomName(4)

  ep <- paste0("in ",sQuote("contact"),": ")

  # if(length(paramnames)!=length(params))
  #   stop("length of ",sQuote("paramnames")," should much with the length of ","params",call.=FALSE)

  statenames <- cleanForC(statenames)
  paramnames <- cleanForC(paramnames)
  covarnames <- cleanForC(covarnames)
  obsnames <- cleanForC(obsnames)
  parm <- list()
  nam <- list()
  ## variable/parameter/observations definitions
  for (v in seq_along(paramnames)) {
    #parm<- c(parm,render(define$var,variable=paramnames[v], value=params[v]))
    nam <- c(nam,render(define$parmcreate, variable=paramnames[v]))
  }
  nam <- c(nam,define$closfun)
  #print("here")
  md <- list()
  for (v in 1:(length(statenames)-1)) {
    md <- c(md ,paste(render(define$modulercpp,variable=paste0(statenames[v],statenames[v+1]))),define$listcreat,nam)

  }

  md <- c(md,paste(render(define$modulercpp,variable="kernelmodel")),define$listcreat,nam)

  if (is.null(name))
    name <- paste0("contact_",digest(" ",serialize=FALSE))

inc <- paste(define$NameSpace,paste(snips,collapse = "\n"),define$modrcpp,paste(md,collapse = "\n"),define$foot)
fx <- cxxfunction(signature(), plugin="Rcpp", include=inc)
csrc <- Module("yada", getDynLib(fx))
  #invisible(list(name = name, src = csrc))
return(csrc)
}

#Function to create the R module
contactRBuilder <- function (name = NULL, snips, statenames, paramnames, covarnames, obsnames, params,
                              verbose = getOption("verbose",FALSE)){

  if (!is.null(name)) name <- cleanForC(name)
  id <- randomName(4)

  ep <- paste0("in ",sQuote("contact"),": ")

  csrc <- paste(paste(snips,collapse = "\n"))
  fileName <- tempfile(fileext = ".R")
  writeLines(csrc, fileName)
  invisible(name = name, src = use(fileName))
}


cleanForC <- function (text) {
  text <- as.character(text)
  text <- gsub("\\.","_",text)
  text <- gsub("-","_",text)
  text
}

define <- list(
  NameSpace = "\n using namespace Rcpp;\n",
  var="double {%variable%} ={%value%};\n",
  varR=" {%variable%} <- {%value%}\n",
  modrcpp = "\nRCPP_MODULE(yada){\n",
  modulercpp= "function( \"{%variable%}\", &{%variable%},\n",
  listcreat = "List::create(",
  parmcreate = " _[\"{%variable%}\"],",
  closfun = " _[\"...\"]));\n",
  foot= "\n}\n\n"
)

## define the rcpp module
header <- list(

  file="/* pomp model file: {%name%} */\n/* {%id%} */\n\n#include <{%pompheader%}>\n#include <R_ext/Rdynload.h>\n\n",
  registration="\nvoid R_init_{%name%} (DllInfo *info)\n{\n",
  rmeasure="\nvoid __pomp_rmeasure (double *__y, const double *__x, const double *__p, const int *__obsindex, const int *__stateindex, const int *__parindex, const int *__covindex, int __ncovars, const double *__covars, double t)\n{\n",
  dmeasure= "\nvoid __pomp_dmeasure (double *__lik, const double *__y, const double *__x, const double *__p, int give_log, const int *__obsindex, const int *__stateindex, const int *__parindex, const int *__covindex, int __ncovars, const double *__covars, double t)\n{\n",
  step.fn="\nvoid __pomp_stepfn (double *__x, const double *__p, const int *__stateindex, const int *__parindex, const int *__covindex, int __covdim, const double *__covars, double t, double dt)\n{\n",
  rate.fn="\ndouble __pomp_ratefn (int j, double t, double *__x, const double *__p, const int *__stateindex, const int *__parindex, const int *__covindex, int __covdim, const double *__covars){\n  double rate = 0.0;  \n",
  skeleton="\nvoid __pomp_skelfn (double *__f, const double *__x, const double *__p, const int *__stateindex, const int *__parindex, const int *__covindex, int __ncovars, const double *__covars, double t)\n{\n",
  initializer="\nvoid __pomp_initializer (double *__x, const double *__p, double t, const int *__stateindex, const int *__parindex, const int *__covindex, const double *__covars)\n{\n",
  fromEstimationScale="\nvoid __pomp_par_trans (double *__pt, const double *__p, const int *__parindex)\n{\n",
  toEstimationScale="\nvoid __pomp_par_untrans (double *__pt, const double *__p, const int *__parindex)\n{\n",
  rprior="\nvoid __pomp_rprior (double *__p, const int *__parindex)\n{\n",
  dprior="\nvoid __pomp_dprior (double *__lik, const double *__p, int give_log, const int *__parindex)\n{\n"
)

randomName <- function (size = 4, stem = "") {
  paste0(stem,
         " Time: ",format(Sys.time(),"%Y-%m-%d %H:%M:%OS3 %z"),
         " Salt: ",
         paste(
           format(
             as.hexmode(ceiling(runif(n=size,max=2^24))),
             upper.case=TRUE
           ),
           collapse=""
         )
  )
}

rcppmod<- list()
render <- function (template, ...) {
  vars=list(...)
  if (length(vars)==0) return(template)
  n <- sapply(vars,length)
  if (!all((n==max(n))|(n==1)))
    stop("in ",sQuote("render"),"incommensurate lengths of replacements",call.=FALSE)
  short <- which(n==1)
  n <- max(n)
  for (i in short) vars[[i]] <- rep(vars[[i]],n)

  retval <- vector(mode="list",length=n)
  for (i in seq_len(n)) {
    tpl <- template
    for (v in names(vars)) {
      src <- sprintf("\\{%%%s%%\\}",v)
      tgt <- vars[[v]][i]
      tpl <- gsub(src,tgt,tpl)
    }
    retval[[i]] <- tpl
  }
  do.call(paste0,retval)
}
