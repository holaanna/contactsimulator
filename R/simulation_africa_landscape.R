#' Simulation of banana landscape to address the question of clean seed in Africa.
#' Clean seed are recovered at the end of the season and then used to generate new plantation.
#'
#' Epidemic simulation using the contact type model with dynamical host densites.
#'
#'\code{Simulate_contact_control_LER_africa} provide a simulation of the epidemic process applying the Malawi management strategy.
#'
#' @inheritParams Simulate_contact_control
#' @param rast A raster file containing info related to grid_lines etc. Need to construnct the object contact to handle all of these
#' @param plant_proc The cells ids guiding the planting protoccol
#' @param prev The disease prevalence in the region. Note that this is redundund when suckers are not imported from external source
#' @param nb_p_grid Cell population size
#' @seealso  \code{\link{Simulate_contact_control_LER_farm}}
#' @references
#' \insertRef{KR08}{contactsimulator}
#' \insertRef{Mee11}{contactsimulator}
#' @example examples/africa_landscape_example.R
#' @export
Simulate_contact_control_LER_africa <- function(rast, param, grid_lines, pop_grid, grid_size=70, age_level=c(1,1),age_dist=c(1,0), m_start=1,t_b=100000, t_max=1000, t_intervention=100000, t_obs=seq(0,1000,100), EI_model=1, kern_model=4,rad=1000,prev=0.003,nb=3,leav=c(3,6),ini=NULL, plant_proc, type_cont=TRUE, nb_p_grid=140){

  #Set parameters
  epsilon <- param$epsilon
  #beta <- param$beta
  b1 <- param$b1
  alpha1 <- param$alpha1
  alpha2 <- param$alpha2
  mu_lat <- param$mu_lat
  var_lat <- param$var_lat
  omega <- param$omega
  c<- param$c
  delta<- param$delta
  t0<- param$t0

  beta_1<- param$beta_1;
  beta_2<- param$beta_2;

  ru <- param$gama

  min_coor_x <- min(c(grid_lines$coor_x_1,grid_lines$coor_x_2)) # bounds of the region
  max_coor_x <- max(c(grid_lines$coor_x_1,grid_lines$coor_x_2))
  min_coor_y <- min(c(grid_lines$coor_y_1,grid_lines$coor_y_2))
  max_coor_y <- max(c(grid_lines$coor_y_1,grid_lines$coor_y_2))

  x_intervals <- seq(min_coor_x,max_coor_x,grid_size)
  y_intervals <- seq(min_coor_y,max_coor_y,grid_size)

  n_line <- max(grid_lines$indx) + 1 #+ 1 # number of lines
  n_row_grid <- nrow(pop_grid) # number of rows of grids
  n_col_grid <- ncol(pop_grid)  # number of cols of grids



  raster::values(rast)[plant_proc[1]]<- nb_p_grid # Initial plantation
  pop_grid<- flipdim(matrix(raster::values(rast),nrow = nrow_grid, byrow = TRUE))
  pop_grid_old = pop_grid

  # else if(is.null(f_rast)){
  #   farm_pos_cat<- data.frame(ro=1,co=1,cat="A",vis_int=360,tim_lst_pos=-10000,nb_round=1,sweep=1)
  # }``



  simulated_epi <- data.frame(k=numeric(0), coor_x=numeric(0), coor_y=numeric(0), t_e=numeric(0), t_i=numeric(0), t_d=numeric(0), t_r=numeric(0),age=numeric(0), infected_source=numeric(0), row=numeric(0), col=numeric(0), typ=numeric(0), season=numeric(0))
  simulated_pro<- data.frame(seas=numeric(0), n_suck_inrn=numeric(0),n_suck_imp=numeric(0), n_plantation=numeric(0), n_remv=numeric(0)) # simulation protocole
  k_11<- 0
  plant_age<- numeric(length(plant_proc))
  plant_age[1]<- 0


  # Extract data from farm and bakckyard
  #print(c(b_rast[355,140],f_rast[355,140]))

#show(param)
  #print(c(b_rast[355,140],f_rast[355,140]))
  if(t_max<=365){
    vist_int<- 2*t_max
  }
  else{
    vist_int<- seq(365,t_max,365)
  }
 # Season
 t_seas <- seq(365,5000,365)
  ## initialize the index cases ##

  if(!is.null(ini)){
    if(!is.data.frame(ini)){
      stop("The initial info ini must be a data frame")
    }
    if(!(all(colnames(ini)%in%c("x", "y", "t_e", "t_i", "row", "col", "typ", "age")))){
      stop("Check the colnames of the initial state ini ")
    }
    m_start<- nrow(ini)
  }

  index_k <- 0:(m_start-1)
  tp_indx<- n_indx<- m_indx<- index_coor_x <- index_coor_y<- numeric(m_start)
  #Assume that the initial infection(s) are located in the centroid
  cent_coor<- sp::coordinates(rast)

  # m_indx<- cent
# show(param)
  if(is.null(ini)){

  for(i in 1:m_start){

     m_grid <-  raster::rowFromCell(rast,plant_proc[1])# the mth row of the grids
     n_grid <-  raster::colFromCell(rast,plant_proc[1])# the mth row of the grids

    index_coor_x[i] <- xFromCell(rast,plant_proc[1])  + sample(c(-1,1),1)*runif(1,0,grid_size/2) #stats::runif(1,min=x_intervals[n_grid],max=x_intervals[n_grid+1]) # random lands at a point in the grid selected
    index_coor_y[i] <- yFromCell(rast,plant_proc[1])  + sample(c(-1,1),1)*runif(1,0,grid_size/2)

    n_indx[i]=n_grid
    m_indx[i]=nrow(pop_grid) - m_grid + 1

  }
#show(param)
  # index_coor_x <- stats::runif(m_start,min_coor_x,max_coor_x)
  # index_coor_y <- stats::runif(m_start,min_coor_y,max_coor_y)
  dt<- stats::rexp(1)/epsilon    # Sellke
  #dt<- 0

  if(is.infinite(dt)){  # When running the simulation from the start e.g for parameter estimation
    index_t_e <- rep(t0,m_start)


    # index_t_i <- index_t_e + E_to_I(EI_model,t0 , mu_lat, var_lat)
  }
  else{
     index_t_e <- rep(t0+dt,m_start)
  }
  # index_t_e <- rep(t0+dt,m_start) # all assumed to infected at time=t0
  #index_t_e <- rep(t0,m_start) # all assumed to infected at time=t0

  # Latent period depending on the model

  index_t_i <- index_t_e + rBTFinv3(EI_model,index_t_e, mu_lat, var_lat, leav[1])

index_age <- sample(age_level,size=m_start,prob=age_dist,replace=T)

  }
  else{  # If initial starting points are provided
    index_coor_x<- ini$x
    index_coor_y<- ini$y
    n_indx<- ini$col+1
    m_indx<- ini$row+1
    tp_indx<- ini$typ
    index_t_e<- ini$t_e
    index_t_i<- ini$t_i
    index_age <- ini$age
    }#####
   # index_t_r <- 2*t_max

   # index_t_r <- index_t_i +  stats::rexp(m_start,rate=1/c)

 index_t_d <- index_t_r <- index_t_i
  for(i in 1:m_start){
    if(index_coor_x[i]<min_coor_x | index_coor_x[i]>max_coor_x | index_coor_y[i]<min_coor_y | index_coor_y[i]>max_coor_y){
       index_t_d[i] <- 2*t_max
     }
     else{
        index_t_d[i] <- index_t_i[i] +  rBTFinv3(EI_model,index_t_i[i], mu_lat, var_lat,nb)
     }
    # if(index_t_r[i]>t_obs){
      index_t_r[i] <- 2*t_max
    # }
  }

  rast_list<- rast

  index_source <- rep(9999,m_start) # all assumed to be infected by the background
  # show(length(index_k))
  # show(length(index_coor_x))
  # show(length(index_coor_y))
  # show(length(index_t_e))
  # show(length(index_t_i))
  # show(length(index_t_d))
  # show(length(index_t_r))
  # show(length(index_age))
  # show(length(index_source))
  # show(length(m_indx))
  # show(length(n_indx))
  # show(length(tp_indx))
  simulated_epi[1:m_start,] <- cbind(index_k, index_coor_x,index_coor_y,index_t_e,index_t_i, index_t_d, index_t_r,index_age,index_source,m_indx,n_indx,tp_indx, plant_age[1])

  #show(param)
  t_now <- min(simulated_epi$t_i)

  simulated_epi_sub <- subset(simulated_epi, simulated_epi$t_i<=t_now & simulated_epi$t_r>t_now) # those are currently infectious

  num_infection <- nrow(simulated_epi)

  t_next <- t_now # to start the while loop

  #show(min(f_rast))
  t_i_new<- index_t_i
  while(t_next<(t_max)){
   # show(t_next)
    #show(min(pop_grid))
    ### simulate the timings, and the source, of next infection ###

    simulated_epi_sub <- subset(simulated_epi, simulated_epi$t_i<=t_now & simulated_epi$t_r>t_now) # those are currently infectious


    #Risk due to secondary infection
    if(nrow(simulated_epi_sub)>=1){

      beta_infectious <- sapply(simulated_epi_sub$age, FUN=beta_by_age, c(beta_1,beta_2))
      total_beta <- sum(beta_infectious)
      joint_I_R <- c(simulated_epi_sub$t_i,simulated_epi_sub$t_r)
      min_I_R <- min(joint_I_R[which(joint_I_R>t_now)])
    }

    #Risk due to primary infection
    if(nrow(simulated_epi_sub)<1){
      total_beta <- 0
      min_I_R <- t_now
    }

    #show(max(pop_grid))
    # show(c(t_now, t_intervention, total_beta, epsilon, omega, b1, t_max) )
    t_next <- simulate_NHPP_next_event (t_now=t_now, t_intervention=t_intervention, sum_beta=total_beta, epsilon=epsilon, omega=omega, b1=b1, t_max=t_max) # simulate the next infection time using thinning algorithm
    # show(t_next)
    min_tim_cont<- t_intervention

    while(t_next>=min(c(min_I_R,min_tim_cont,min(t_obs), t_seas[1])) & t_next!=Inf & t_max>min(c(min_I_R,min_tim_cont,min(t_obs), t_seas[1]))){  # If next event is a removal

      kk<- (min(min_I_R,min_tim_cont)==min_tim_cont) +1
     #show(c(kk))
      tim_next_event<- min(c(min_I_R,min_tim_cont,min(t_obs), t_seas[1]))

      if(tim_next_event==min_I_R){ # next event is the baseline process
                #indx_inf<- which(simulated_epi$t_i==min_I_R)  # Increase number of plants infected in the cell
                #farm_inf[simulated_epi$row[indx_inf],simulated_epi$col[indx_inf]]<- farm_inf[simulated_epi$row[indx_inf],simulated_epi$col[indx_inf]] +1
                t_now <- min_I_R

                ll<- which(simulated_epi_sub$t_r==min_I_R)

      }

      else if(tim_next_event==min(t_obs)){ # Removal during the survey
        #show(t_max)
         plan_rem<- simulated_epi%>%dplyr::filter(t_d<t_now, t_r>t_now)%>%dplyr::arrange(season)
         plan_rem<- as.numeric(plan_rem[,1])
         #show(plan_rem)
         if(length(plan_rem)>0){
           samp<- ceiling(length(plan_rem)*delta)
           simulated_epi[plan_rem[1:samp] + 1,"t_r"]<- min(t_obs)  # Remove all
         }
         t_obs[which(t_obs==min(t_obs))]<- 2*t_max
      }
      # End of season. Generate new plantationa/cell for the new season
      else{
      #show(t_max)
       # show(c(1,min(values(rast))))
        indx_pop <- which(raster::values(rast)>0)
        indx_npop<- which(raster::values(rast)==0)
        #print(length(indx_pop))
        rast_inf<- rast
        raster::values(rast_inf)<- 0
        # Plant remove during the previous season
        plant_rem<- subset(simulated_epi, simulated_epi$t_r<=t_now & simulated_epi$t_r>t_seas[1]-365)
        #show(plant_rem)
        # Cells where removals occured
        count_rm_cell<- 0
        if(nrow(plant_rem)>0){
          cell_rm <- cellFromXY(rast, cbind(plant_rem$coor_x,plant_rem$coor_y))
          count_rm_cell<- table(cell_rm)
          raster::values(rast)[as.numeric(names(count_rm_cell))]<- raster::values(rast)[as.numeric(names(count_rm_cell))] - as.numeric(count_rm_cell)

          #show(cell_rm)

        }
        #show(c(2,min(raster::values(rast))))
         # show(count_rm_cell)
         # show(length(plant_rem))
                   #show(count_rm_cell)
        # Nb of suckers generated by infected plants but not detected
        plant_inf<- subset(simulated_epi, simulated_epi$t_r>=t_now & simulated_epi$t_i>t_seas[1]-365 & typ==0)

        if(nrow(plant_inf)==0){
          count_inf_cell<- 0
        }
        else{
          cell_inf <- cellFromXY(rast, cbind(plant_inf$coor_x,plant_inf$coor_y))
          count_inf_cell<- table(cell_inf)

        }
        raster::values(rast_inf)[as.numeric(names(count_inf_cell))]<- raster::values(rast_inf)[as.numeric(names(count_inf_cell))] + as.numeric(count_inf_cell)

        raster::values(rast)<- raster::values(rast) - raster::values(rast_inf)
         #N_suckers_non_inf<- raster::values(rast)
        #
        # if(any(N_suckers_non_inf>0)){
        #   N_suckers_non_inf[which(N_suckers_non_inf>0)]<- sapply(N_suckers_non_inf[which(N_suckers_non_inf>0)],function(x){
        #     return(sum(replicate(x,sample(2:5,1))))
        #   })
        # }
        # N_suckers_non_inf<- sapply(raster::values(rast),function(x){
        #      return(sum(replicate(x,sample(2:5,1))))
        #    })

        N_suckers_non_inf<- raster::values(rast)
        N_suckers_non_inf[N_suckers_non_inf>0]<- sapply(N_suckers_non_inf[N_suckers_non_inf>0], function(x){
          return(sum(sample(2:5,x,replace = T)))
        })
        N_suckers_non_inf<- N_suckers_non_inf + raster::values(rast)
        # N_suckers_inf<- raster::values(rast_inf)
        # if(any(N_suckers_inf>0)){
        #   N_suckers_inf[which(N_suckers_inf>0)]<- sapply(N_suckers_inf[which(N_suckers_inf>0)],function(x){
        #     return(sum(replicate(x,sample(2:5,1))))
        #   })
        # }
        N_suckers_inf<- raster::values(rast_inf)
        if(!all(N_suckers_inf==0)){
          N_suckers_inf[N_suckers_inf>0]<- sapply(N_suckers_inf[N_suckers_inf>0], function(x){
          return(sum(sample(2:5,x,replace = T)))
        })
        N_suckers_inf<- N_suckers_inf + raster::values(rast_inf)
        }

        #N_suckers_inf<- raster::values(rast_inf)*sample(2:5,length(raster::values(rast_inf)), replace = T) + raster::values(rast_inf)

        N_suckers <- N_suckers_non_inf + N_suckers_inf
        N_new_plantation <- sum(N_suckers)%/%nb_p_grid
        N_suckers_remain <- sum(N_suckers)%%nb_p_grid

         #Proportion of infected suckers

        p_suck_inf<- (sum(N_suckers_inf) - sum(raster::values(rast_inf)))*runif(1,0.7,0.9)/(sum(N_suckers))
        # Cell to fill in with infected suckers from within the plantation
        if(nrow(plant_rem)>0){
          rm_cell<- floor(as.numeric(count_rm_cell)*p_suck_inf)
          for(i in 1:length(count_rm_cell)){
                if(rm_cell[i]>0){
                   m_indx <-  raster::rowFromCell(rast,names(count_rm_cell)[i])# the mth row of the grids
                   n_indx <-  raster::colFromCell(rast,names(count_rm_cell)[i])# the mth row of the grids

                  index_coor_x <- xFromCell(rast,names(count_rm_cell)[i])  + sample(c(-1,1),rm_cell[i],replace = T)*runif(rm_cell[i],0,grid_size/2) #stats::runif(1,min=x_intervals[n_grid],max=x_intervals[n_grid+1]) # random lands at a point in the grid selected
                  index_coor_y <- yFromCell(rast,names(count_rm_cell)[i])  + sample(c(-1,1),rm_cell[i],replace = T)*runif(rm_cell[i],0,grid_size/2) #stats::runif(1,min=x_intervals[n_grid],max=x_intervals[n_grid+1]) # random lands at a point in the grid selected

                  k <- (num_infection+1):(num_infection + rm_cell[i]) # k=0,1,2...
                  index_t_e<- rep(t_seas[1],rm_cell[i])
                  index_t_i <- index_t_e +  rBTFinv3(EI_model,index_t_e,mu_lat,var_lat,leav[1]) #E_to_I(EI_model,t_now , mu_lat, var_lat)
                  index_t_d <- index_t_i + rBTFinv3(EI_model,index_t_i, mu_lat, var_lat,nb)
                  index_t_r <- rep(2*t_max,rm_cell[i])

                  simulated_epi[k,] <- cbind(k, index_coor_x,index_coor_y,index_t_e,index_t_i, index_t_d, index_t_r,1,9999,nrow(pop_grid)-m_indx+1,n_indx,0,plant_age[names(count_rm_cell)[i]])

                  num_infection<- num_infection + rm_cell[i]
                }


              }
        }
       #show(t_max)
        # Cells with surplus suckers and non-surplus
        cell_suplus<- which(N_suckers>nb_p_grid)
        cell_non_suplus<- which(N_suckers<nb_p_grid & N_suckers>0)

        N_suckers_surplus<- sum(N_suckers) - nb_p_grid*length(cell_suplus)

        if(N_new_plantation>length(plant_proc)){
          N_new_plantation<- length(plant_proc)
          N_suckers_remain<- 0
        }
        #show(c(t_seas[1],length(cell_suplus),length(cell_non_suplus),N_new_plantation,N_suckers_remain))
        if(length(cell_non_suplus)==0){
           if(length(cell_suplus)==0){
            stop("No suckers left to continue.")
           }
          else{
            #if(length(cell_suplus)<length(indx_pop)){
              if(N_new_plantation>length(indx_pop)){

                raster::values(rast)[indx_pop]<- nb_p_grid
                raster::values(rast)[plant_proc[(length(indx_pop)+1):N_new_plantation]]<- nb_p_grid
                plant_age[(length(indx_pop)+1):N_new_plantation]<- t_seas[1]
                if(N_suckers_remain>0){
                  if(type_cont){
                    raster::values(rast)[plant_proc[N_new_plantation+1]] <- nb_p_grid

                    Nb_inf_imp<- ceiling(prev*(nb_p_grid-N_suckers_remain))
                    if(Nb_inf_imp>0){
                      plant_age[N_new_plantation+1]<- t_seas[1]
                      #show(c(N_new_plantation,num_infection,Nb_inf_imp))
                      #show(Nb_inf_imp)
                      #randomly smaple cell in wich imported infected suckers lie in
                      rand_cel_suck<- sample(plant_proc[(length(indx_pop)+1):(N_new_plantation+1)],1)
                      m_indx <-  raster::rowFromCell(rast,rand_cel_suck)# the mth row of the grids
                      n_indx <-  raster::colFromCell(rast,rand_cel_suck)# the mth row of the grids

                      index_coor_x <- xFromCell(rast,rand_cel_suck)  + sample(c(-1,1),Nb_inf_imp,replace = T)*runif(Nb_inf_imp,0,grid_size/2) #stats::runif(1,min=x_intervals[n_grid],max=x_intervals[n_grid+1]) # random lands at a point in the grid selected
                      index_coor_y <- yFromCell(rast,rand_cel_suck)  + sample(c(-1,1),Nb_inf_imp,replace = T)*runif(Nb_inf_imp,0,grid_size/2) #stats::runif(1,min=x_intervals[n_grid],max=x_intervals[n_grid+1]) # random lands at a point in the grid selected

                      k <- (num_infection+1):(num_infection + Nb_inf_imp) # k=0,1,2...
                      index_t_e<- rep(t_seas[1],Nb_inf_imp)
                      index_t_i <- index_t_e +  rBTFinv3(EI_model,index_t_e,mu_lat,var_lat,leav[1]) #E_to_I(EI_model,t_now , mu_lat, var_lat)
                      index_t_d <- index_t_i +  rBTFinv3(EI_model,index_t_i, mu_lat, var_lat,nb)
                      index_t_r <- rep(2*t_max,Nb_inf_imp)

                      simulated_epi[k,] <- cbind(k, index_coor_x,index_coor_y,index_t_e,index_t_i, index_t_d, index_t_r,1,9999,nrow(pop_grid)-m_indx+1,n_indx,0,plant_age[N_new_plantation+1])
                       # show(c(t_max,1))
                      num_infection<- num_infection + Nb_inf_imp
                    }


                  }
                  else{
                    raster::values(rast)[plant_proc[N_new_plantation+1]] <- N_suckers_remain

                  }
                }
               # show(min(values(rast)))
              }
              else{
                raster::values(rast)[cell_suplus]<- nb_p_grid
                raster::values(rast)[plant_proc[1:N_new_plantation][-cell_suplus]]<- nb_p_grid

                 if(N_suckers_remain>0){
                  if(type_cont){
                    raster::values(rast)[plant_proc[N_new_plantation+1]] <- nb_p_grid # Add proportion infected here
                    plant_age[N_new_plantation+1]<- t_seas[1]
                    rand_cel_suck<- sample(plant_proc[1:(N_new_plantation+1)][-cell_suplus],1)
                    # Nber of infected suckers imported
                    Nb_inf_imp<- ceiling(prev*(nb_p_grid-N_suckers_remain))
                    if(Nb_inf_imp>0){
                      m_indx <-  raster::rowFromCell(rast,rand_cel_suck)# the mth row of the grids
                      n_indx <-  raster::colFromCell(rast,rand_cel_suck)# the mth row of the grids

                      index_coor_x <- xFromCell(rast,rand_cel_suck)  + sample(c(-1,1),Nb_inf_imp,replace = T)*runif(Nb_inf_imp,0,grid_size/2) #stats::runif(1,min=x_intervals[n_grid],max=x_intervals[n_grid+1]) # random lands at a point in the grid selected
                      index_coor_y <- yFromCell(rast,rand_cel_suck)  + sample(c(-1,1),Nb_inf_imp,replace = T)*runif(Nb_inf_imp,0,grid_size/2) #stats::runif(1,min=x_intervals[n_grid],max=x_intervals[n_grid+1]) # random lands at a point in the grid selected

                      k <- (num_infection+1):(num_infection + Nb_inf_imp) # k=0,1,2...
                      index_t_e<- rep(t_seas[1],Nb_inf_imp)
                      index_t_i <- index_t_e +  rBTFinv3(EI_model,index_t_e,mu_lat,var_lat,leav[1]) #E_to_I(EI_model,t_now , mu_lat, var_lat)
                      index_t_d <- index_t_i +  rBTFinv3(EI_model,index_t_i, mu_lat, var_lat,nb)
                      index_t_r <- rep(2*t_max,Nb_inf_imp)

                      simulated_epi[k,] <- cbind(k, index_coor_x,index_coor_y,index_t_e,index_t_i, index_t_d, index_t_r,1,9999,nrow(pop_grid)-m_indx+1,n_indx,0,plant_age[N_new_plantation+1])

                      num_infection<- num_infection + Nb_inf_imp
                    }

                  }
                  else{
                    raster::values(rast)[plant_proc[N_new_plantation+1]] <- N_suckers_remain

                  }
                }
              }
          #}
          }


        }
        else{
          if(length(cell_suplus)==0){
            if(type_cont){

              # Nber of infected suckers imported
              for(i in 1:length(cell_non_suplus)){
                Nb_inf_imp<- ceiling(prev*(nb_p_grid-raster::values(rast)[cell_non_suplus[i]]))
                if(Nb_inf_imp>0){
                  m_indx <-  raster::rowFromCell(rast,cell_non_suplus[i])# the mth row of the grids
                  n_indx <-  raster::colFromCell(rast,cell_non_suplus[i])# the mth row of the grids
                  #plant_age[cell_non_suplus[i]]<- t_seas[1]

                  index_coor_x <- xFromCell(rast,cell_non_suplus[i])  + sample(c(-1,1),Nb_inf_imp,replace = T)*runif(Nb_inf_imp,0,grid_size/2) #stats::runif(1,min=x_intervals[n_grid],max=x_intervals[n_grid+1]) # random lands at a point in the grid selected
                  index_coor_y <- yFromCell(rast,cell_non_suplus[i])  + sample(c(-1,1),Nb_inf_imp,replace = T)*runif(Nb_inf_imp,0,grid_size/2) #stats::runif(1,min=x_intervals[n_grid],max=x_intervals[n_grid+1]) # random lands at a point in the grid selected

                  k <- (num_infection+1):(num_infection + Nb_inf_imp) # k=0,1,2...
                  index_t_e<- rep(t_seas[1],Nb_inf_imp)
                  index_t_i <- index_t_e +  rBTFinv3(EI_model,index_t_e,mu_lat,var_lat,leav[1]) #E_to_I(EI_model,t_now , mu_lat, var_lat)
                  index_t_d <- index_t_i + rBTFinv3(EI_model,index_t_i, mu_lat, var_lat,nb)
                  index_t_r <- rep(2*t_max,Nb_inf_imp)

                  simulated_epi[k,] <- cbind(k, index_coor_x,index_coor_y,index_t_e,index_t_i, index_t_d, index_t_r,1,9999,nrow(pop_grid)-m_indx+1,n_indx,0,plant_age[cell_non_suplus[i]])

                  num_infection<- num_infection + Nb_inf_imp
                }


              }
              raster::values(rast)[cell_non_suplus] <- nb_p_grid # Add proportion infected here
            }

          }
          else{
            #show(t_seas[1])
            if(N_suckers_surplus>sum((nb_p_grid - cell_non_suplus))){
              raster::values(rast)[cell_non_suplus]<- nb_p_grid

              raster::values(rast)[indx_pop]<- nb_p_grid
              raster::values(rast)[plant_proc[(length(indx_pop)+1):N_new_plantation]]<- nb_p_grid
              plant_age[(length(indx_pop)+1):N_new_plantation]<- t_seas[1]
                if(N_suckers_remain>0){
                  if(type_cont){
                    raster::values(rast)[plant_proc[N_new_plantation+1]] <- nb_p_grid
                    Nb_inf_imp<- ceiling(prev*(nb_p_grid-N_suckers_remain))
                    if(Nb_inf_imp>0){
                      m_indx <-  raster::rowFromCell(rast,plant_proc[N_new_plantation+1])# the mth row of the grids
                      n_indx <-  raster::colFromCell(rast,plant_proc[N_new_plantation+1])# the mth row of the grids
                      plant_age[N_new_plantation+1]<- t_seas[1]
                      rand_cel_suck<- sample(plant_proc[(length(indx_pop)+1):(N_new_plantation+1)],1)

                      index_coor_x <- xFromCell(rast,rand_cel_suck)  + sample(c(-1,1),Nb_inf_imp,replace = T)*runif(Nb_inf_imp,0,grid_size/2) #stats::runif(1,min=x_intervals[n_grid],max=x_intervals[n_grid+1]) # random lands at a point in the grid selected
                      index_coor_y <- yFromCell(rast,rand_cel_suck)  + sample(c(-1,1),Nb_inf_imp,replace = T)*runif(Nb_inf_imp,0,grid_size/2) #stats::runif(1,min=x_intervals[n_grid],max=x_intervals[n_grid+1]) # random lands at a point in the grid selected

                      k <- (num_infection+1):(num_infection + Nb_inf_imp) # k=0,1,2...
                      index_t_e<- rep(t_seas[1],Nb_inf_imp)
                      index_t_i <- index_t_e +  rBTFinv3(EI_model,index_t_e,mu_lat,var_lat,leav[1]) #E_to_I(EI_model,t_now , mu_lat, var_lat)
                      index_t_d <- index_t_i + rBTFinv3(EI_model,index_t_i, mu_lat, var_lat,nb)
                      index_t_r <- rep(2*t_max,Nb_inf_imp)

                      simulated_epi[k,] <- cbind(k, index_coor_x,index_coor_y,index_t_e,index_t_i, index_t_d, index_t_r,1,9999,nrow(pop_grid)-m_indx+1,n_indx,0,plant_age[N_new_plantation+1])

                      num_infection<- num_infection + Nb_inf_imp
                    }


                  }
                  else{
                    raster::values(rast)[plant_proc[N_new_plantation+1]] <- N_suckers_remain
                  }
                }
            }
            else{
              # Randomly select from the suckers generated and fill each remaining cell proportional to their deficit
              samp_suck<- rmultinom(1,N_suckers_surplus,(nb_p_grid - cell_non_suplus)/sum((nb_p_grid - cell_non_suplus)))  #(700 - cell_non_suplus)/sum((700 - cell_non_suplus))*N_suckers_surplus
              raster::values(rast)[cell_non_suplus] <- raster::values(rast)[cell_non_suplus] + samp_suck
              if(type_cont){
                for(i in 1:length(cell_non_suplus)){
                Nb_inf_imp<- ceiling(prev*(nb_p_grid-raster::values(rast)[cell_non_suplus[i]]))
                if(Nb_inf_imp>0){
                  m_indx <-  raster::rowFromCell(rast,cell_non_suplus[i])# the mth row of the grids
                  n_indx <-  raster::colFromCell(rast,cell_non_suplus[i])# the mth row of the grids

                  index_coor_x <- xFromCell(rast,cell_non_suplus[i])  + sample(c(-1,1),Nb_inf_imp,replace = T)*runif(Nb_inf_imp,0,grid_size/2) #stats::runif(1,min=x_intervals[n_grid],max=x_intervals[n_grid+1]) # random lands at a point in the grid selected
                  index_coor_y <- yFromCell(rast,cell_non_suplus[i])  + sample(c(-1,1),Nb_inf_imp,replace = T)*runif(Nb_inf_imp,0,grid_size/2) #stats::runif(1,min=x_intervals[n_grid],max=x_intervals[n_grid+1]) # random lands at a point in the grid selected

                  k <- (num_infection+1):(num_infection + Nb_inf_imp) # k=0,1,2...
                  index_t_e<- rep(t_seas[1],Nb_inf_imp)
                  index_t_i <- index_t_e +  rBTFinv3(EI_model,index_t_e,mu_lat,var_lat,leav[1]) #E_to_I(EI_model,t_now , mu_lat, var_lat)
                  index_t_d <- index_t_i + rBTFinv3(EI_model,index_t_i, mu_lat, var_lat,nb)
                  index_t_r <- rep(2*t_max,Nb_inf_imp)

                  simulated_epi[k,] <- cbind(k, index_coor_x,index_coor_y,index_t_e,index_t_i, index_t_d, index_t_r,1,9999,nrow(pop_grid)-m_indx+1,n_indx,0,plant_age[cell_non_suplus[i]])

                  num_infection<- num_infection + Nb_inf_imp
                }


              }
              raster::values(rast)[cell_non_suplus] <- nb_p_grid # Add proportion infected here
              }
              }

          }
        }


     k_11<- k_11 + 1
    if(type_cont){
      if(N_suckers_remain>0){
        simulated_pro[k_11,]<- c(t_seas[1], sum(N_suckers), nb_p_grid-N_suckers_remain, N_new_plantation,sum(as.numeric(count_rm_cell)))
      }
      else{
         simulated_pro[k_11,]<- c(t_seas[1], sum(N_suckers), 0, N_new_plantation,sum(as.numeric(count_rm_cell)))
      }

    }
     else{
       simulated_pro[k_11,]<- c(t_seas[1], sum(N_suckers), 0, N_new_plantation,sum(as.numeric(count_rm_cell)))
     }


     rast_list<- raster::stack(rast_list,rast)
     t_seas<- t_seas[-1]

     pop_grid<- flipdim(matrix(raster::values(rast), nrow = nrow_grid, byrow = TRUE))
     #show(min(pop_grid))
    # Number of cryptic suckers
      }
      #show(c(t_next,min(pop_grid)))
      simulated_epi_sub <- subset(simulated_epi, simulated_epi$t_i<=t_now & simulated_epi$t_r>t_now) # those are currently infectious
#   show(nrow(simulated_epi_sub))

      if(nrow(simulated_epi_sub)>=1){

        beta_infectious <- sapply(simulated_epi_sub$age, FUN=beta_by_age, c(beta_1,beta_2))

        total_beta <- sum(beta_infectious)
        joint_I_R <- c(simulated_epi_sub$t_i,simulated_epi_sub$t_r)
        min_I_R <- min(joint_I_R[which(joint_I_R>t_now)])

        t_next <- simulate_NHPP_next_event (t_now=t_now,  t_intervention=t_intervention, sum_beta=total_beta, epsilon=epsilon, omega=omega, b1=b1,t_max=t_max)

      }
      # show(c(min_tim_cont,min_I_R,t_now,nrow(simulated_epi_sub),t_next))
      if(nrow(simulated_epi_sub)<1){
        total_beta <- 0
        t_next <- simulate_NHPP_next_event (t_now=t_now,  t_intervention=t_intervention, sum_beta=total_beta, epsilon=epsilon, omega=omega, b1=b1, t_max=t_max)
        #print(c(t_next,t_now,3))
        break # needed when consider background infection
      }

      # show(farm_pos_cat)
      if(t_max<min(min_I_R,min_tim_cont,t_next)){
        # show(c(min_tim_cont,min_I_R,t_next))
        break;
      }




    } # end of while(t_next>=min_I_R)

    #show(c(t_next,t_now,k,nrow(simulated_epi_sub)))
    k <- num_infection + 1 - 1 # k=0,1,2...
    t_old<- t_now
    t_now <- t_next
    t_i_new <- t_now +  rBTFinv3(EI_model,t_now,mu_lat,var_lat,leav[1]) #E_to_I(EI_model,t_now , mu_lat, var_lat)
    t_d_new <- t_i_new + rBTFinv3(EI_model,t_i_new, mu_lat, var_lat,nb)
    t_r_new <- 2*t_max
    # if(t_now<t_obs){
    #   t_r_new <- t_i_new + stats::rexp(1,rate=1/c)
    #   if(t_r_new>t_obs){
    #     t_r_new <- 2*t_max
    #   }
    # }
    # else{
    #   t_r_new <- 2*t_max
    # }




#show(c(t_next,nrow(simulated_epi_sub),beta_infectious))
    # if(nrow(simulated_epi_sub)>=1) source <- sample(c(9999,simulated_epi_sub$k),size=1, prob=c(epsilon/365,beta_infectious/365*(1+b1*cos(omega*t_old)))) # 9999 = background
if(nrow(simulated_epi_sub)>=1) source <- sample(c(9999,simulated_epi_sub$k),size=1, prob=c(epsilon,beta_infectious)) # 9999 = background
    #show(c(source,nrow(simulated_epi)))
#show(c(t_next,min(pop_grid)))

    if(nrow(simulated_epi_sub)<1) source <- 9999 # 9999 = background
#show(nrow(simulated_epi_sub))
    age <- sample(age_level,size=1,prob=age_dist)
    #print(c(t_next,t_now,2))
    ### simulate the coordinates of the new infection (above) ###

    x_new = min_coor_x - 5 # to start the while loop
    y_new = min_coor_y - 5 # to start the while loop

    m_1=n_1=1000
    #pop_grid = f_rast + b_rast
    #show(min(pop_grid))
    siz<- -1
    tt<- 0

    while(x_new<min_coor_x | x_new>max_coor_x | y_new<min_coor_y | y_new>max_coor_y | siz<0){
      if(is.na(simulated_epi$coor_x[source+1]) & source<9999){
        show(c(1,x_new,y_new,siz,r,source,nrow(simulated_epi),simulated_epi$coor_x[source+1],simulated_epi$coor_y[source+1]))
        show(simulated_epi[source+1,])
      }
            # show(c(t_next,min(pop_grid),max(pop_grid),source))
      if (source!=9999){
        if(simulated_epi[source+1,"typ"]==1){
          r <- abs(Samp_dis (kern_model,ru, alpha2, alpha2))
        }
        else{
          r <- abs(Samp_dis (kern_model,ru, alpha1, alpha1))
        }


        set_points <- circle_line_intersections (circle_x=simulated_epi$coor_x[source+1],circle_y=simulated_epi$coor_y[source+1], r,  n_line=n_line, grid_lines=grid_lines)

        n_set_points = nrow(set_points)
        # show(c(r,n_set_points))

        if (n_set_points>=1) {

          arcs <- func_arcs_attributes(set_points, pop_grid, r, min_coor_x, min_coor_y, grid_size, n_row_grid, n_col_grid)
          #arcs<- na.omit(arcs)
          #show(arcs)
          arcs[is.na(arcs)] <- 0
          arcs$mass <- arcs$dens*arcs$len_arc

          arcs$theta=set_points$theta

          #arcs <- arcs[order(arcs$theta_abs),]

          sum_arcs_den <- sum(arcs$dens)
          #show(c(r,simulated_epi$coor_x[source+1],simulated_epi$coor_y[source+1],min(arcs$mass)))
          #show(min(arcs$mass))
          #show(c(sum_arcs_den,min(arcs$mass),min(arcs$dens),arcs$len_arc))
          if (sum_arcs_den>0){
            k_segment <- sample(1:nrow(arcs),size=1,prob=arcs$mass) # decide which segment the new infection would lie
            #print(k_segment)
            theta_within_segment <- stats::runif(1, min=0, max=arcs$theta_abs[k_segment]) # uniformly draws a point within the segment chosen above

            if (k_segment==1) theta_from_y_eq_0 <- arcs$theta[nrow(arcs)] + theta_within_segment # the measure of theta from y=0
            if (k_segment!=1) theta_from_y_eq_0 <- arcs$theta[k_segment-1] + theta_within_segment # the measure of theta from y=0

            x_new <- simulated_epi$coor_x[source+1] + r*cos(theta_from_y_eq_0) # the coor_x of the new infection
            y_new <- simulated_epi$coor_y[source+1] + r*sin(theta_from_y_eq_0) # the coor_x of the new infection
            # print(c(r*cos(theta_from_y_eq_0)))
            n=ceiling((x_new-min_coor_x)/grid_size)
            m=ceiling((y_new-min_coor_y)/grid_size)

            m=m_1=arcs$m_th_row[k_segment]
            n=n_1=arcs$n_th_col[k_segment]

          }

          ###
          if (sum_arcs_den==0){
            tt<- 1
            break
            k_grid <- sample(1:length(as.numeric(pop_grid)),size=1, prob=as.numeric(pop_grid))
            m_grid <- k_grid%%nrow(pop_grid) # the mth row of the grids

            if(m_grid==0) m_grid <- nrow(pop_grid)
            n_grid <- ceiling(k_grid/nrow(pop_grid))  # nth column ..
            x_new <- stats::runif(1,min=x_intervals[n_grid],max=x_intervals[n_grid+1]) # random lands at a point in the grid selected
            y_new <- stats::runif(1,min=y_intervals[m_grid],max=y_intervals[m_grid+1])

            n=n_grid
            m=m_grid
            m_1 = m
            n_1 = n

          }
          ###

        #show(arcs)
        }

        #print(c(r,n_set_points))
        if (n_set_points<1) {

          theta_from_y_eq_0 <-  stats::runif(1, min=0,max=(2*pi)) # draw theta uniformly between [0,2pi]

          x_new <- simulated_epi$coor_x[source+1] + r*cos(theta_from_y_eq_0) # the coor_x of the new infection
          y_new <- simulated_epi$coor_y[source+1] + r*sin(theta_from_y_eq_0) # the coor_x of the new infection

          n=ceiling((x_new-min_coor_x)/grid_size)
          m=ceiling((y_new-min_coor_y)/grid_size)
          m_1 = m
          n_1 = n

          if(x_new<min_coor_x | x_new>max_coor_x | y_new<min_coor_y | y_new>max_coor_y){
            tt<- 1
            break
          }

          #message(c(r))
        }

        # print(c(r,n_set_points))
      }

      if(source==9999){
        # r <- abs(Samp_dis (kern_model,ru, alpha1, alpha2))
        # if(t_next<365){
        #
        # }
        # if(r<)
        #print(pop_grid)
        k_grid <- sample(1:length(as.numeric(pop_grid)),size=1, prob=as.numeric(pop_grid))
        m_grid <- k_grid%%nrow(pop_grid) # the mth row of the grids
        if(m_grid==0) m_grid <- nrow(pop_grid)
        n_grid <- ceiling(k_grid/nrow(pop_grid))  # nth column ..
        x_new <- stats::runif(1,min=x_intervals[n_grid],max=x_intervals[n_grid+1]) # random lands at a point in the grid selected
        y_new <- stats::runif(1,min=y_intervals[m_grid],max=y_intervals[m_grid+1])

        n=n_grid
        m=m_grid
        m_1 = m
        n_1 = n


      }
    # show(c(t_next,max(pop_grid)))
      #show(c(m,n,source))
      siz<- pop_grid[m,n]-1
      #show(c(2,x_new,y_new,siz,m,n))
      # show(c(x_new, y_new, siz))
    } # end while(x_new<min_coor_x | x_new>max_coor_x | y_new<min_coor_y | y_new>max_coor_y){
    #show(c(t_next,t_now,3))
    if(tt==0){
      if(is.infinite(t_next)){
        break
      }
      else{

        pop_grid[m,n]<- pop_grid[m,n] - 1
        simulated_epi[k+1,] <- c(k, x_new, y_new, t_now, t_i_new, t_d_new, t_r_new, age, source, m,n,0, plant_age[cellFromRowCol(rast,nrow(pop_grid)-m+1,n)])

       }
       num_infection <- num_infection + 1
    }

    # Control proccedures; what if scenarios


  } # end of while(t_next<t_max)
  # show(farm_pos_cat)

  return(list(simulated_epi, simulated_pro,rast_list))
}

