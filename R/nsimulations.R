# ==============================================================================
# Generate nsim simulation from the contact model
#==============================================================================

#' Generate a certain number of temporal observation from the process
#'
#' Temporal oservation using the contact type model.
#'
#'\code{Simulate} Count the number of observation during a certain period.
#'
#' @inheritParams Simulate_contact_model
#' @param nsim The total size of the simulation requierd.
#' @param time The sequence of time at which the observation will occure.
#' @param obs  The obervation requiered from during the period time
#'
#' @return A data frame with the time, observation made and the nth simulation.
#'
#' @example
#' @export
Simulate<- function(param, grid_lines, pop_grid, grid_size=5000, age_level=c(1,1),age_dist=c(1,0), m_start=1, t_max=118, t_intervention=365, EI_model=1, kern_model=4, times, obs, nsim){


  OBS <- data.frame(times=times[-1], Obs=obs, nsim="data")
  for(i in 1:nsim){
    print(i)
    epi<- Simulate_contact_model(param,grid_lines=grid_lines, pop_grid=pop_grid,grid_size = size,t_max = t_max,t_intervention = t_intervention, EI_model = EI_model, kern_model = kern_model)

    df<- data.frame(times=times[-1], Obs=Sub_set(epi$t_r,times), nsim=as.character(i))
    OBS<- rbind(OBS,df)
  }
  return(OBS)
  }
