#' Simulation of the epidemic control process considering a distribution for the leaf
#'
#' Epidemic simulation using the contact type model with the Australian control strategies.
#'
#'\code{Simulate_contact_control_LER_farm} provide the simulation of the epidemic process with the Australian BBTV management plan uisnt the obserbation times.
#' @param t_obs The sequence of time at which the observations are made
#' @param ini A data frame of the intial state:
#'       \describe{
#'         \item{row}{The row at which the cell containing the host lies in}
#'         \item{col}{The column at which the cell containing the host lies in}
#'         \item{typ}{The type of the host (e.g. backyard 1 /plantation 0 hosts)}
#'         \item{x}{x-coordinate}
#'         \item{y}{y-coordinate}
#'         \item{t_e}{exposure time}
#'         \item{t_i}{infection time}
#'        }
#'
#' @inheritParams Simulate_contact_control
#' @seealso  \code{\link{Simulate_contact_control}}
#' @references
#' \insertRef{KR08}{contactsimulator}
#' \insertRef{Mee11}{contactsimulator}
#'
#' @example examples/plantation_example.R
#' @export
Simulate_contact_control_LER_farm <- function(f_rast=NULL, b_rast=NULL, farm_pos_cat=NULL, vis_int_per_cat=NULL, param, grid_lines, pop_grid, grid_size=500, age_level=c(1,1),age_dist=c(1,0), m_start=1,t_b=100000, t_max=1000, t_intervention=100000, t_obs=seq(0,1000,100), EI_model=1, kern_model=4,rad=1000,sweep_prop=c(.5,.5),back_p=c(.7,.5),rate_det=c(0.3,1),int_det=c(30,90,180),nb_in_b=1,nb=3,leav=c(3,6),ini=NULL){

  #Set parameters
  epsilon <- param$epsilon
  beta <- param$beta
  b1 <- param$b1
  alpha1 <- param$alpha1
  alpha2 <- param$alpha2
  mu_lat <- param$mu_lat
  var_lat <- param$var_lat
  omega <- param$omega
  c<- param$c
  delta<- param$delta
  t0<- param$t0

  beta_1<- beta;
  beta_2<- beta;

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

  pop_grid_old = pop_grid
  #Plantations
  if(is.null(f_rast) & is.null(b_rast) & is.null(farm_pos_cat) & is.null(farm_pos_cat)){
    f_rast<- pop_grid
    b_rast<- pop_grid*0
    nb_in_b<- 0

		k_grid <- sample(1:length(as.numeric(pop_grid)),size=1, prob=as.numeric(pop_grid))
		m_grid <- k_grid%%nrow(pop_grid) # the mth row of the grids
		if(m_grid==0) m_grid <- nrow(pop_grid)
		n_grid <- ceiling(k_grid/nrow(pop_grid))  # nth column .
    farm_pos_cat<- data.frame(ro=m_grid,co=n_grid,cat="A",vis_int=360,tim_lst_pos=-10000,nb_round=1,sweep=1)

  }
  # else if(is.null(f_rast)){
  #   farm_pos_cat<- data.frame(ro=1,co=1,cat="A",vis_int=360,tim_lst_pos=-10000,nb_round=1,sweep=1)
  # }``

if(max(t_obs)<=t0){
  farm_pos_cat$cat<- "A"
  farm_pos_cat$tim_lst_pos<- -10000
  farm_pos_cat$vis_int<- 360
}


  simulated_epi <- data.frame(k=numeric(0), coor_x=numeric(0), coor_y=numeric(0), t_e=numeric(0), t_i=numeric(0), t_d=numeric(0), t_r=numeric(0),age=numeric(0), infected_source=numeric(0), row=numeric(0), col=numeric(0), typ=numeric(0))

  # Extract data from farm and bakckyard
  #print(c(b_rast[355,140],f_rast[355,140]))
  b_rast<- b_rast*nb
  f_rast_p<- f_rast/(b_rast+f_rast)              # proportion occupied by farm and backyard plant resp
  b_rast_p<- b_rast/(b_rast+f_rast)
  f_rast_p[is.na(f_rast_p)]<- 0
  b_rast_p[is.na(b_rast_p)]<- 0
  farm_pos_cat$vis<- max(t_obs) + farm_pos_cat$vis_int
  farm_pos_cat$cat<- as.character(farm_pos_cat$cat)
  farm_inf<- array(0,c(nrow(f_rast),ncol(f_rast)))
  time_sinc_last_pos<- array(0,c(nrow(f_rast),ncol(f_rast)))
  time_sinc_last_pos1<- array(0,c(nrow(f_rast),ncol(f_rast)))
#show(param)
  #print(c(b_rast[355,140],f_rast[355,140]))
  if(t_max<=365){
    vist_int<- 2*t_max
  }
  else{
    vist_int<- seq(365,t_max,365)
  }

  ## initialize the index cases ##

  if(!is.null(ini)){
    if(!is.data.frame(ini)){
      stop("The initial info ini must be a data frame")
    }
    if(!(all(colnames(ini)%in%c("x", "y", "t_e", "t_i", "row", "col", "typ")))){
      stop("Check the colnames of the initial state ini ")
    }
    m_start<- nrow(ini)
  }

  index_k <- 0:(m_start-1)
  tp_indx<- n_indx<- m_indx<- index_coor_x <- index_coor_y<- numeric(m_start)
# show(param)
  if(is.null(ini)){

  for(i in 1:m_start){
    if(t0<=max(t_obs)){
      ss<- 1:nrow(farm_pos_cat)
      k_grid <- sample(ss,size=1)
    }
    else{
      ss<- which(farm_pos_cat$cat=="E" | farm_pos_cat$cat=="D" | farm_pos_cat$cat=="C" )
      if(length(ss)==0 ){
        ss<- which(farm_pos_cat$cat=="B" )
        if(length(ss)==0){
          ss<- which(farm_pos_cat$cat=="A" )
          if(length(ss)==0){
            stop("Farm need to be in one of the categories: A, B, C or D")

          }
        }

        k_grid <- sample(ss,size=1)
        #show(i)
        # k_grid <- sample(1:length(as.numeric(pop_grid)),size=1, prob=as.numeric(pop_grid))
        # m_grid <- k_grid%%nrow(pop_grid) # the mth row of the grids
        # if(m_grid==0) m_grid <- nrow(pop_grid)
        # n_grid <- ceiling(k_grid/nrow(pop_grid))  # nth column ..
      }

      else{
        k_grid <- sample(ss,size=1)


      }
    }
#show(param)
    m_grid <-  farm_pos_cat$ro[k_grid]# the mth row of the grids
    n_grid <-  farm_pos_cat$co[k_grid]# the mth row of the grids

    index_coor_x[i] <- stats::runif(1,min=x_intervals[n_grid],max=x_intervals[n_grid+1]) # random lands at a point in the grid selected
    index_coor_y[i] <- stats::runif(1,min=y_intervals[m_grid],max=y_intervals[m_grid+1])

    n_indx[i]=n_grid
    m_indx[i]=m_grid
     #show(c(m_grid,n_grid,dim(f_rast_p)))
    u1<- f_rast_p[m_grid,n_grid]
    u2<- b_rast_p[m_grid,n_grid]

#show(c(u1,u2))
    u<- max(u1,u2)
    if(length(ss)==0 ){
      if(u2>0){
        tp_indx[i]<- 1
        b_rast[m_grid,n_grid]<- b_rast[m_grid,n_grid] - 1
      }
      else{
        tp_indx[i]<- 0
        f_rast[m_grid,n_grid]<- f_rast[m_grid,n_grid] - 1
        farm_inf[m_grid,n_grid]<- farm_inf[m_grid,n_grid] + 1

      }
    }

    else{
      tp_indx[i]<- 0
      f_rast[m_grid,n_grid]<- f_rast[m_grid,n_grid] - 1
      farm_inf[m_grid,n_grid]<- farm_inf[m_grid,n_grid] + 1

    }



  }
  #show(c(m_grid,n_grid))
  # index_coor_x <- stats::runif(m_start,min_coor_x,max_coor_x)
  # index_coor_y <- stats::runif(m_start,min_coor_y,max_coor_y)
  dt<- stats::rexp(1)/(epsilon/365.0)    # Sellke
  #dt<- 0

  if(dt>0){  # When running the simulation from the start e.g for parameter estimation
    index_t_e <- rep(t0+dt,m_start)
    # index_t_i <- index_t_e + E_to_I(EI_model,t0 , mu_lat, var_lat)
  }
  else{
    index_t_e <- rep(t0,m_start)

  }
  # index_t_e <- rep(t0+dt,m_start) # all assumed to infected at time=t0
  #index_t_e <- rep(t0,m_start) # all assumed to infected at time=t0

  # Latent period depending on the model

  index_t_i <- index_t_e + rBTFinv3(EI_model,index_t_e, mu_lat, var_lat, leav[1])



  }
  else{  # If initial starting points are provided
    index_coor_x<- ini$x
    index_coor_y<- ini$y
    n_indx<- ini$col+1
    m_indx<- ini$row+1
    tp_indx<- ini$typ
    index_t_e<- ini$t_e
    index_t_i<- ini$t_i
    # b_in<- which(tp_indx==1)
    # f_in<- which(tp_indx==0)
    # if(length(b_in)>0){
    #   b_rast[m_grid,n_grid]<- b_rast[m_grid,n_grid] - 1
    # }
    for(i in 1:m_start){
      if(tp_indx[i]==1){
        b_rast[m_indx[i],n_indx[i]]<- b_rast[m_indx[i],n_indx[i]] - 1
      }
      else{
        #print(c(f_rast[m_indx[i],n_indx[i]],m_indx[i],n_indx[i]))
        f_rast[m_indx[i],n_indx[i]]<- f_rast[m_indx[i],n_indx[i]] - 1
        farm_inf[m_indx[i],n_indx[i]]<- farm_inf[m_indx[i],n_indx[i]] + 1
      }
      show(c(m_indx[i],n_indx[i]))
    }

    }#####
   # index_t_r <- 2*t_max

   # index_t_r <- index_t_i +  stats::rexp(m_start,rate=1/c)
index_t_d <- index_t_i +  rBTFinv3(EI_model,index_t_i, mu_lat, var_lat,nb)
 index_t_r <- index_t_i
  for(i in 1:m_start){
    # if(index_t_r[i]>t_obs){
      index_t_r[i] <- 2*t_max
    # }
  }

  index_age <- sample(age_level,size=m_start,prob=age_dist,replace=T)
  index_source <- rep(9999,m_start) # all assumed to be infected by the background

# show(length(index_k))
# show(length(index_coor_x))
# show(length(index_coor_y))
# show(length(index_t_e))
# show(length(index_t_i))
# show(length(index_t_r))
# show(length(index_age))
# show(length(index_source))
# show(length(m_indx))
# show(length(n_indx))
# show(length(tp_indx))
  simulated_epi[1:m_start,] <- c(index_k, index_coor_x,index_coor_y,index_t_e,index_t_i, index_t_d, index_t_r,index_age,index_source,m_indx,n_indx,tp_indx)

#show(param)
  t_now <- min(simulated_epi$t_i)

  simulated_epi_sub <- subset(simulated_epi, simulated_epi$t_i<=t_now & simulated_epi$t_r>t_now) # those are currently infectious

  num_infection <- nrow(simulated_epi)
 if(t_now>min(farm_pos_cat$vis)){ # To ensure that the initial visit occurs beyond the first infection
   farm_pos_cat$vis<- t_now + farm_pos_cat$vis_int
 }

  t_next <- t_now # to start the while loop

  #show(min(f_rast))
  t_i_new<- index_t_i
  while(t_next<(t_max)){
    #show(t_next)
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


    t_next <- simulate_NHPP_next_event (t_now=t_now, t_intervention=t_intervention, sum_beta=total_beta, epsilon=epsilon, omega=omega, b1=b1, t_max=t_max) # simulate the next infection time using thinning algorithm
    #show(c(t_next,t_now,1))
    #print(c(t_next,total_beta,omega,t_max,t_now,t_intervention))
    min_tim_cont<- min(farm_pos_cat$vis)
    #show(min(min_tim_cont))
    #show(min(c(min_I_R,min_tim_cont,min(t_obs))))
    while(t_next>=min(c(min_I_R,min_tim_cont,min(t_obs))) & t_next!=Inf & t_max>min(c(min_I_R,min_tim_cont,min(t_obs)))){  # If next event is a removal
#show(min(min_tim_cont))
      kk<- (min(min_I_R,min_tim_cont)==min_tim_cont) +1
     #show(c(kk))
      tim_next_event<- min(c(min_I_R,min_tim_cont,min(t_obs)))

      if(tim_next_event==min_I_R){ # next event is the baseline process
                #indx_inf<- which(simulated_epi$t_i==min_I_R)  # Increase number of plants infected in the cell
                #farm_inf[simulated_epi$row[indx_inf],simulated_epi$col[indx_inf]]<- farm_inf[simulated_epi$row[indx_inf],simulated_epi$col[indx_inf]] +1
                t_now <- min_I_R

                ll<- which(simulated_epi_sub$t_r==min_I_R)
                #show(ll)
                if(length(ll)>0){
                  if(simulated_epi[ll,"typ"]==0){

                    farm_inf[simulated_epi[ll,"row"],simulated_epi[ll,"col"]]<- farm_inf[simulated_epi[ll,"row"],simulated_epi[ll,"col"]] -1

                    f_rast[simulated_epi[ll,"row"],simulated_epi[ll,"col"]]<- f_rast[simulated_epi[ll,"row"],simulated_epi[ll,"col"]] +1
                  }
                  else{
                    b_rast[simulated_epi[ll,"row"],simulated_epi[ll,"col"]]<- b_rast[simulated_epi[ll,"row"],simulated_epi[ll,"col"]] +1

                  }
                  #pop_grid[simulated_epi[ll,"row"],simulated_epi[ll,"col"]]<- pop_grid[simulated_epi[ll,"row"],simulated_epi[ll,"col"]] +1

                }





      }
      else if(tim_next_event==min_tim_cont){  # control sweep
                indx_con<- which(farm_pos_cat$vis==min_tim_cont)
                ll<- (min_tim_cont>t_b) + 1  # baseline or alternative
                for(i in 1:length(indx_con)){
                   #show(c(i,length(indx_con),t_now))

                  nn_1<- farm_pos_cat[indx_con[i],c("co")]
                  mm_1<- farm_pos_cat[indx_con[i],c("ro")]
                  #show(c(mm_1,nn_1,farm_inf[mm_1,nn_1]))
                  # if(farm_inf[mm_1,nn_1]>0){
                  farm_inf[mm_1,nn_1]<- 0
                  plan_in_farm<- as.numeric(subset(simulated_epi, row==mm_1 & col==nn_1 & typ==0 & t_i<t_now &t_r>t_now)[,1])
                  plan_det<- integer(0)
                  if(length(plan_in_farm)>0 ){

                    if(ll==2 & (farm_pos_cat$cat[indx_con[i]]=="A" | farm_pos_cat$cat[indx_con[i]]=="B")){ # Growers monitoring
                      plan_in_farm_det_tim<- simulated_epi[plan_in_farm+1,"t_i"] + rBTFinv3(1,simulated_epi[plan_in_farm+1,"t_i"],mu_lat,var_lat,leav[2]-4)
                      samp1<- plan_in_farm[which(plan_in_farm_det_tim<t_now)]
                      if(length(samp1)>0){
                        samp<- runif(length(samp1))
                        indx_det<- samp1[which(samp<=rate_det[1])]
                        if(length(indx_det)>0){
                          plan_det<- indx_det
                        }

                      }

                    }
                    else{
                      # print(farm_pos_cat$cat[indx_con[i]])

                      expert_det<- simulated_epi[plan_in_farm+1,"t_i"] #+ rBTFinv3(1,simulated_epi[plan_in_farm+1,"t_i"],mu_lat,var_lat,leav[1])

                      # print(plan_in_farm)
                      if(length(expert_det)>0){
                        samp1<- plan_in_farm[which(expert_det<t_now)]

                        if(length(samp1)>0){
                          samp<- runif(length(samp1))
                          indx_det<- samp1[which(samp<=rate_det[2])]
                          if(length(indx_det)>0){
                            plan_det<- indx_det
                          }

                        }

                      }



                    }

                    if(length(plan_det)!=0){
                      # show(plan_in_farm)
                      simulated_epi[plan_det+1,"t_r"]<- min_tim_cont  # Remove all
                      f_rast[simulated_epi[plan_det+1,"row"],simulated_epi[plan_det+1,"col"]]<- f_rast[simulated_epi[plan_det+1,"row"],simulated_epi[plan_det+1,"col"]] +1
                      farm_inf[mm_1,nn_1]<- length(plan_det)

                      if(runif(1)<sweep_prop[ll]){ # Check if the plantation is considered for sweep
                        # Removal of backyard plant in a certain radius
                        x_intervals
                        # x_center<- min_coor_x + (nn_1-1)*grid_size + grid_size/2
                        # y_center<-  min_coor_y + (mm_1-1)*grid_size + grid_size/2

                        x_center<- x_intervals[nn_1] + grid_size/2
                        y_center<-  y_intervals[mm_1] + grid_size/2

                        back_in_sweep<- as.numeric(subset(simulated_epi, typ==1 & t_i<=t_now &t_r>t_max)[,1])
                        #show(length(back_in_sweep))
                        if(length(back_in_sweep)!=0){
                          #Select proportion in the sweep
                          expert_det<- simulated_epi[ back_in_sweep+1,"t_i"] #+ rBTFinv3(1,simulated_epi[ back_in_sweep+1,"t_i"],mu_lat,var_lat,leav[1])


                          back_in_sweep_det<- back_in_sweep[which(expert_det<t_now)]
                          if(length(back_in_sweep_det)>0){
                            prob<- round(length(back_in_sweep_det)*back_p[ll])
                            samp<- sample(back_in_sweep_det,prob)

                            dist_pot<- distanc(as.matrix(simulated_epi[samp+1,c("coor_x","coor_y")]),c(x_center,y_center))
                            indx_rem<- which(dist_pot<rad)
                            # show(dist_pot)
                            # show(c(x_center,y_center))
                            # show(c(mm_1,nn_1,prob))
                           # show(c(length(back_in_sweep),prob,min(dist_pot),nn_1,mm_1,length(indx_rem),length(back_in_sweep_det),x_center,y_center,t_now))
                            # show(samp)
                            if(length(indx_rem)>0){
                              simulated_epi[samp[indx_rem]+1,"t_r"]<- min_tim_cont  # Remove all
                              b_rast[simulated_epi[samp[indx_rem]+1,"row"],simulated_epi[samp[indx_rem]+1,"col"]]<- b_rast[simulated_epi[samp[indx_rem]+1,"row"],simulated_epi[samp[indx_rem]+1,"col"]] +1

                            }
                          }


                        }
                      }
                    }



                  }


                  # recorded in the last year

                  rem_last_year<- subset(simulated_epi, row==mm_1 & col==nn_1 & typ==0 & t_r>t_now-360)
                  farm_inf_on<- farm_inf
                  if(length(rem_last_year)>0){
                    farm_inf[mm_1,nn_1]<-  farm_inf[mm_1,nn_1] + nrow(rem_last_year)
                  }



                  # }

                  if(farm_pos_cat$vis[indx_con[i]]<farm_pos_cat$tim_lst_pos[indx_con[i]]){
                    # show(c(farm_inf[mm_1,nn_1],farm_inf_on[mm_1,nn_1]))
                    # show(farm_pos_cat[indx_con[i],])
                    # show(t_now)
                    #show(c(t_now,min_tim_cont-farm_pos_cat$tim_lst_pos[indx_con[i]]))
                  }

                  if(ll==1){  # baseline
                    baseline_con(farm_pos_cat,farm_inf,farm_inf_on,min_tim_cont,indx_rem,indx_con[i],mm_1,nn_1,int_det)
                  }
                  else{ # alternative
                    baseline_alt(farm_pos_cat,vis_int_per_cat,farm_inf,farm_inf_on,min_tim_cont,indx_rem,indx_con[i],mm_1,nn_1)
                  }
                  # if(farm_pos_cat$vis[indx_con[i]]<farm_pos_cat$tim_lst_pos[indx_con[i]]){
                  #   show(c(farm_inf[mm_1,nn_1],farm_inf_on[mm_1,nn_1]))
                  #   show(farm_pos_cat[indx_con[i],])
                  #   show(c(t_now,min_tim_cont-farm_pos_cat$tim_lst_pos[indx_con[i]]))
                  # }
                  # show(farm_inf[mm_1,nn_1])
                  # show(farm_pos_cat$vst)
                }
                t_now <- min_tim_cont
                # simulated_epi_sub <- subset(simulated_epi, simulated_epi$t_i<=t_now & simulated_epi$t_r>t_now)# those are currently infectious

                min_tim_cont<- min(farm_pos_cat$vis)
                if(ll==1){
                  # show(c(min_tim_cont,min_I_R,t_next))
                  #show(farm_pos_cat)
                }
                # show(farm_pos_cat)
                # show(c(min_tim_cont,ll))
                # show(c(min_tim_cont,min_I_R,t_next))
              }

      else{ # Removal during the survey
         plan_rem<- as.numeric(subset(simulated_epi, t_d<t_now &t_r>t_now)[,1])
         if(length(plan_rem)>0){
           samp<- round(length(plan_rem)*delta)
           simulated_epi[sample(plan_rem,samp) + 1,"t_r"]<- min(t_obs)  # Remove all
         }
         t_obs[which(t_obs==min(t_obs))]<- 2*t_max
      }
      simulated_epi_sub <- subset(simulated_epi, simulated_epi$t_i<=t_now & simulated_epi$t_r>t_now) # those are currently infectious
# show(nrow(simulated_epi_sub))

      if(nrow(simulated_epi_sub)>=1){

        beta_infectious <- sapply(simulated_epi_sub$age, FUN=beta_by_age, c(beta_1,beta_2))

        total_beta <- sum(beta_infectious)
        joint_I_R <- c(simulated_epi_sub$t_i,simulated_epi_sub$t_r)
        min_I_R <- min(joint_I_R[which(joint_I_R>t_now)])

        t_next <- simulate_NHPP_next_event (t_now=t_now,  t_intervention=t_intervention, sum_beta=total_beta, epsilon=epsilon, omega=omega, b1=b1,t_max=t_max)


        for(j in 1:nrow(simulated_epi_sub)){
          if(simulated_epi_sub$typ[j]==0){
            farm_inf[simulated_epi_sub$row[j],simulated_epi_sub$col[j]]<- farm_inf[simulated_epi_sub$row[j],simulated_epi_sub$col[j]] + 1
          }
        }
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




#show(sample(c(9999,simulated_epi_sub$k),size=1, prob=c(epsilon,beta_infectious)))
    # if(nrow(simulated_epi_sub)>=1) source <- sample(c(9999,simulated_epi_sub$k),size=1, prob=c(epsilon/365,beta_infectious/365*(1+b1*cos(omega*t_old)))) # 9999 = background
if(nrow(simulated_epi_sub)>=1) source <- sample(c(9999,simulated_epi_sub$k),size=1, prob=c(epsilon,beta_infectious)) # 9999 = background

    if(nrow(simulated_epi_sub)<1) source <- 9999 # 9999 = background
#show(nrow(simulated_epi_sub))
    age <- sample(age_level,size=1,prob=age_dist)
    #print(c(t_next,t_now,2))
    ### simulate the coordinates of the new infection (above) ###

    x_new = min_coor_x - 5 # to start the while loop
    y_new = min_coor_y - 5 # to start the while loop
    #print(t_now)
    # if(source!=9999){
    #   #print(c(simulated_epi$row[source+1],simulated_epi$col[source+1],source,k_grid))
    #   ru<- pop_grid[simulated_epi$row[source+1],simulated_epi$col[source+1]]/pop_grid_old[simulated_epi$row[source+1],simulated_epi$col[source+1]]
    #
    # }
    # ru <- .47
    # show(source)
    m_1=n_1=1000
    pop_grid = f_rast + b_rast
    # show(min(pop_grid))
    siz<- -1

    while(x_new<min_coor_x | x_new>max_coor_x | y_new<min_coor_y | y_new>max_coor_y | siz<0){
    #   show(x_new)
    # show(y_new)
    # show(min_coor_x)
    # show(max_coor_x )
    # show(min_coor_y)
    # show(max_coor_y)
      # show((x_new<min_coor_x) +(x_new>max_coor_x) + (y_new<min_coor_y) + (y_new>max_coor_y))
      # show((siz))
      # show(t_next)
      #pop_grid = pop_grid_before
      #show(c(min(as.numeric(pop_grid))))
      if (source!=9999){
        # show(kern_model)
        # show(ru)
        # show(alpha1)
        # show(alpha2)
        #show(c(kern_model,ru,alpha1,alpha2))
        r <- abs(Samp_dis (kern_model,ru, alpha1, alpha2))

        set_points <- circle_line_intersections (circle_x=simulated_epi$coor_x[source+1],circle_y=simulated_epi$coor_y[source+1], r,  n_line=n_line, grid_lines=grid_lines)

        n_set_points = nrow(set_points)
        # show(c(r,n_set_points,min(as.numeric(pop_grid))))
        if (n_set_points>=1) {
          # show(set_points)
          # show(source)
          # show(simulated_epi)
          arcs <- func_arcs_attributes(set_points, pop_grid, r, min_coor_x, min_coor_y, grid_size, n_row_grid, n_col_grid)

          arcs$mass <- arcs$dens*arcs$len_arc

          arcs$theta=set_points$theta

          #arcs <- arcs[order(arcs$theta_abs),]

          sum_arcs_den <- sum(arcs$dens)
           # print(c(sum_arcs_den,n_line,min(as.numeric(pop_grid))))
          # print(c(simulated_epi$coor_x[source+1],simulated_epi$coor_y[source+1]))
          # print(as.data.frame(arcs$n_th_col,arcs$m_th_row))
          # show(c(sum(arcs$mass),0))
          if (sum_arcs_den>0){
            # print(nrow(arcs))
            #print(arcs$mass)
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

            #print(c(arcs$n_th_col[k_segment],arcs$m_th_row[k_segment]))
            # print(c(m,n,k_segment,n_set_points))
            # print(c(pop_grid[m,n],r))
            # print(arcs)
            # print(set_points)
            # #print(c(simulated_epi$coor_x[source+1],simulated_epi$coor_y[source+1],r))

            m_1=arcs$m_th_row[k_segment]
            n_1=arcs$n_th_col[k_segment]
            # if(any(c(m,n)!=c(arcs$m_th_row[k_segment],arcs$n_th_col[k_segment]))){
            #   print(c(m,n,k_segment,n_set_points))
            #   print(c(pop_grid[m,n],r))
            #   print(c(simulated_epi$coor_x[source+1],simulated_epi$coor_y[source+1],r))
            #   print(arcs)
            #   print(set_points)
            #   print(c(x_new,y_new))
            # }

            # print(c(m,n,1))
          }

          ###
          if (sum_arcs_den==0){
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

          #message(c(r))
        }

        # print(c(r,n_set_points))
      }

      if(source==9999){

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


      siz<- pop_grid[m,n]-1
    } # end while(x_new<min_coor_x | x_new>max_coor_x | y_new<min_coor_y | y_new>max_coor_y){
    #show(c(t_next,t_now,3))
    if(is.infinite(t_next)){
      break
    }
    else{
      # label the plant
      # show(c(m,n))
      f_rast_p<- f_rast/(b_rast+f_rast)              # proportion occupied by farm and backyard plant resp
      b_rast_p<- b_rast/(b_rast+f_rast)
      f_rast_p[is.na(f_rast_p)]<- 0
      b_rast_p[is.na(b_rast_p)]<- 0

      u1<- f_rast_p[m,n]
      u2<- b_rast_p[m,n]

      if(all(c(u1,u2)==0)){
        tp<- 1
        #show(c(k,u1,u2,source,n_set_points,sum_arcs_den,r,m,n))
      }
      else{
        u<- max(u1,u2)
        #print(c(u1,u2,m,n))
        if(u==u1){
          if(u>runif(1)){
            tp<- 0
            f_rast[m,n]<- f_rast[m,n] - 1
          }
          else{
           # tp<- 1
            tp<- 0
            b_rast[m,n]<- b_rast[m,n] - 1
          }

        }
        else{
          if(u>runif(1)){
            #tp<- 1
            tp<- 0
            b_rast[m,n]<- b_rast[m,n] - 1
          }
          else{
            tp<- 0
            # f_rast[m,n]<- f_rast[m,n] - 1
          }
        }





      }

      # if(tp==0){
      #   farm_inf[m,n]<- farm_inf[m,n] + 1
      # }
      #pop_grid[m,n]=pop_grid[m,n]-1
      simulated_epi[k+1,] <- c(k, x_new, y_new, t_now, t_i_new, t_d_new, t_r_new, age, source, m,n,tp)

    }
    ####

    num_infection <- num_infection + 1

    # Control proccedures; what if scenarios


  } # end of while(t_next<t_max)
  # show(farm_pos_cat)

  return(simulated_epi)
}

