#-----------------------------------------------------------------------------------------------------------------
#                             Temporal trajectory
#-----------------------------------------------------------------------------------------------------------------
#' Temporal trajectory to computing the confidence band
#'
#' Trajectory.
#'
#'\code{Trajectory} Provide a the trjaectory of the simulated epidemic.
#'
#' @param epi A matrix specifying the simulated epidemic. Note that the object should be of the form of the output of \code{\link{Simulate_contact_control}} or \code{\link{Simulate_contact_model}}
#' @param Dat_obs A 2 columns data frame wiht compoents:
#'          \describe{
#'            \item{times}{ The sequence of times at which the observations are made}
#'            \item{obs}{the number of observation recorded at a given time}
#'          }
#' @param  var The observation variable. t_r is the default.
#' @return A data frame of observation. By default it considers the removals.
#'
#' @examples
#' @export
Trajectory<-function(epi=NULL, Dat_obs=NULL, var="t_r"){

  # Control the input
  if(is.null(epi) | is.null(Dat_obs)){
    stop("Invalid input")
  }

  if(!(all(colnames(epi)%in%c("k", "coor_x", "coor_y", "t_e", "t_i", "t_r","age", "infected_source", "row", "col", "typ"))) | !(all(c("times","obs")%in%colnames(Dat_obs)))){
    stop("Check the colnames of your inputs ")
  }

  tra<- traj(epi[,var], Dat_obs$times)

  return(data.frame(Dat_obs,sim=tra))

}


#-----------------------------------------------------------------------------------------------------------------
#                             Temporal trajectory for all simulations
#-----------------------------------------------------------------------------------------------------------------
#' Temporal trajectory to computing the confidence band
#'
#' Trajectories.
#'
#'\code{Trajectory_all} Provide a the trjaectories of the simulated epidemics.
#'
#' @param epi A list of matrix specifying the simulated epidemics. Note that the object should be of the form of the output of \code{\link{Simulate_contact_control}} or \code{\link{Simulate_contact_model}}
#' @param Dat_obs A 2 columns data frame wiht compoents:
#'          \describe{
#'            \item{times}{ The sequence of times at which the observations are made}
#'            \item{obs}{the number of observation recorded at a given time}
#'          }
#' @param  var The observation variable. t_r is the default.
#' @return A data frame consting of the the time points and the observation with the simulation. By default it considers the removals.
#'
#' @examples
#' @export
Trajectory_all<-function(epi=NULL, Dat_obs=NULL, var="t_r"){

  # Control the input
  if(is.null(epi) | is.null(Dat_obs)){
    stop("Invalid input")
  }

  if(!(all(colnames(epi[[1]])%in%c("k", "coor_x", "coor_y", "t_e", "t_i", "t_r","age", "infected_source", "row", "col", "typ"))) | !(all(c("times","obs")%in%colnames(Dat_obs)))){
    stop("Check the colnames of your inputs ")
  }
  xx<- lapply(epi,FUN = function(x){
    tra<- traj(x[,var], Dat_obs$times)
  })
  sim=do.call(cbind,xx)

  return(data.frame(times=Dat_obs$times,obs=Dat_obs$obs,sim))

}


