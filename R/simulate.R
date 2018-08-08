simulate.internal <- function(object, nsim = 1L, seed = NULL, params,
                               states = FALSE, obs = FALSE,
                               times, t0, as.data.frame = FALSE,
                               include.data = FALSE,
                               .getnativesymbolinfo = TRUE,
                               verbose = getOption("verbose", FALSE), ...){

  ep <- paste0("in ",sQuote("simulate"),": ")

  if (missing(times))
    times <- time(object,t0=FALSE)
  else
    times <- as.numeric(times)

  if (missing(t0))
    t0 <- timezero(object)
  else
    t0 <- as.numeric(t0)

  obs <- as.logical(obs)
  states <- as.logical(states)
  as.data.frame <- as.logical(as.data.frame)
  include.data <- as.logical(include.data)

  if (missing(params))
    params <- coef(object)

  if (length(params)==0)
    stop(ep,"no ",sQuote("params")," specified",call.=FALSE)

  params <- as.matrix(params)

  ## set the random seed (be very careful about this)
  seed <- as.integer(seed)
  if (length(seed)>0) {
    if (!exists('.Random.seed',envir=.GlobalEnv)) set.seed(NULL)
    save.seed <- get('.Random.seed',envir=.GlobalEnv)
    set.seed(seed)
  }

  if (!obs && !states)
    object <- as(object,"contact")
  #Load the object
  data=data@data
  initializer = initializer
  grid.lines = grid.lines
  pop.grid = pop.grid
  params = params
  grid.size = grid.size
  age.level = age.level
  age.dist = age.dist
  t.max = t.max
  t.intervention = t.intervention
  kernel.model1 = kernel.model1
  kernel.model2 = kernel.model2
  EI.model = EI.model
  r.model = r.model
  PACKAGE = PACKAGE
  d.model = d.model
  nkernls = nkernls
  userdata=data@userdata

}
