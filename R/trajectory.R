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

  if(!(all(colnames(epi)%in%c("k", "coor_x", "coor_y", "t_e", "t_i",  "t_d","t_r","age", "infected_source", "row", "col", "typ"))) | !(all(c("times","obs")%in%colnames(Dat_obs)))){
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

  if(!(all(colnames(epi[[1]])%in%c("k", "coor_x", "coor_y", "t_e", "t_i", "t_d","t_r","age", "infected_source", "row", "col", "typ"))) | !(all(c("times","obs")%in%colnames(Dat_obs)))){
    stop("Check the colnames of your inputs ")
  }
  xx<- lapply(epi,FUN = function(x){
    tra<- traj(x[,var], Dat_obs$times)
  })
  sim=do.call(cbind,xx)

  return(data.frame(times=Dat_obs$times,obs=Dat_obs$obs,sim))

}



#-----------------------------------------------------------------------------------------------------------------
#                             Temporal trajectory for all simulations for Ripley index
#-----------------------------------------------------------------------------------------------------------------
#' Temporal trajectory to computing the confidence band for Ripley index
#'
#' Trajectories.
#'
#'\code{Trajectory_ripley} Provide a the trjaectories of the simulated epidemics.
#'
#' @param epi A list of matrix specifying the simulated epidemics. Note that the object should be of the form of the output of \code{\link{Simulate_contact_control}} or \code{\link{Simulate_contact_model}}
#' @param Dat_obs A 2 columns data frame wiht compoents:
#'          \describe{
#'            \item{times}{ The sequence of times at which the observations are made}
#'            \item{obs}{the number of observation recorded at a given time}
#'          }
#' @param  var The observation variable. t_r is the default.
#' @param  rast The raster file of the region.
#' @return A data frame consting of the the time points and the observation with the simulation. By default it considers the removals.
#'
#' @examples
#' @export
Trajectory_ripley<-function(epi=NULL, Dat_obs=NULL, rast=NULL){

  # Control the input
  if(is.null(epi) | is.null(Dat_obs)){
    stop("Invalid input")
  }

  if(!(all(c("coor_x","coor_y")%in%colnames(epi[[1]]))) | !(all(c("coor_x","coor_y")%in%colnames(Dat_obs)))){
    stop("Check the colnames of your inputs ")
  }

  if(is.null(rast)){
    message("You have not provided the raster file - the extent of the observed will be used")
    min_x<- min(Dat_obs$coor_x)
    min_y<- min(Dat_obs$coor_y)
    max_x<- max(Dat_obs$coor_x)
    max_y<- max(Dat_obs$coor_y)
    win_x<- c(min_x-1,max_x+1)
    win_y<- c(min_y-1,max_y+1)
  }
  else{
    min_x<- xmin(rast)
    min_y<- ymin(rast)
    max_x<- xmax(rast)
    max_y<- ymax(rast)
    win_x<- c(min_x,max_x)
    win_y<- c(min_y,max_y)
  }
  options(warn = -1)
  u<- spatstat::ppp(Dat_obs$coor_x,Dat_obs$coor_y,win_x,win_y)
  L<- as.data.frame(suppressWarnings(spatstat::Kest(u,correction="iso")))
  DF<- array(c(L$r,L$iso),c(nrow(L),2))

  xx<- lapply(epi,FUN = function(x){
    u<- spatstat::ppp(x$coor_x,x$coor_y,win_x,win_y)
    L<- as.data.frame(suppressWarnings(spatstat::Kest(u,correction="iso")))
    return(L$iso)
  })
  sim=do.call(cbind,xx)

  return(data.frame(r=DF[,1],obs=DF[,2],sim))

}


#-----------------------------------------------------------------------------------------------------------------
#                             Temporal evolution of the max distance
#-----------------------------------------------------------------------------------------------------------------
#' Temporal trajectory to computing the wave speed
#'
#' Trajectories.
#'
#'\code{Trajectory_dis} Provide a the trjaectories of the maximum distance the inoculum has travelled.
#'
#' @param epi A list of matrix specifying the simulated epidemics. Note that the object should be of the form of the output of \code{\link{Simulate_contact_control}} or \code{\link{Simulate_contact_model}}
#' @param Dat_obs A 2 columns data frame with components:
#'          \describe{
#'            \item{times}{ The sequence of times at which the observations are made}
#'            \item{obs}{the number of observation recorded at a given time}
#'          }
#' @param  var The observation variable. t_r is the default.
#' @return A data frame consting of the the time points and the max distance the wave has traveled. By default it considers the removals.
#'
#' @examples
#' @export
Trajectory_dis<-function(epi=NULL, Dat_obs=NULL, var="t_r"){

  # Control the input
  if(is.null(epi) | is.null(Dat_obs)){
    stop("Invalid input")
  }

  if(!(all(colnames(epi[[1]])%in%c("k", "coor_x", "coor_y", "t_e", "t_i", "t_d","t_r","age", "infected_source", "row", "col", "typ"))) | !(all(c("times","obs")%in%colnames(Dat_obs)))){
    stop("Check the colnames of your inputs ")
  }
  xx<- lapply(epi,FUN = function(x){
    first_inf<- as.numeric(x[1,c("coor_x", "coor_y")])
    df<- x%>%dplyr::mutate(dis=sqrt((first_inf[1]-coor_x)^2+(first_inf[2]-coor_y)^2))
    tra<- sapply(Dat_obs$times,FUN = function(tt){
      df_filt<- df%>%dplyr::filter(floor(t_r)<=tt)
      if(nrow(df_filt)==0){
        mx_d<- 0
      }
      else{
        mx_d<- max(df%>%dplyr::filter(floor(t_r)<=tt)%>%dplyr::select(dis))/tt
      }

    })
  })
  sim=do.call(cbind,xx)

  return(data.frame(times=Dat_obs$times,obs=Dat_obs$obs,sim))

}


