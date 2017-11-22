#' Simulation of the epidemic process
#'
#' Epidemic simulation using the contact type model.
#'
#'\code{Simulate_contact_model} provide the simulation of the epidemic process and the
#'       maximum distance the wave can travel at each particular obaservation date.
#'
#' @param param Indicating a data frame contaning a vector of parameters including:
#'       \describe{
#'         \item{epsion}{The primary infection rate. See \code{\link{func_time_beta}}}
#'         \item{beta}{Rate of infection. See \code{\link{func_time_beta}}}
#'         \item{alpha}{The dispersal kernel parameter.}
#'         \item{mu_lat,var_lat}{mean and variance of the latent period. See \code{\link{E_to_I}} for details.}
#'         \item{omega}{Parameter characterising the seasality. See \code{\link{func_time_beta}}}
#'        }
#' @inheritParams func_arcs_attributes
#' @param age_level,age_dist Vectors of age level and the propportion of each age group respectively. See details.
#' @param m_start The size of initial cases. Default is 1.
#' @param t_max  Final observation time.
#' @param t_intervention Start of the intervention if any.
#' @param EI_model Take integer values to specify the type of model used for the latent period. See \code{\link{E_to_I}}
#' @return A data frame with components:
#'       \describe{
#'         \item{k}{Index of the case}
#'         \item{coor_x,coor_y}{Location of the infection premisse}
#'         \item{t_e,ti,tr}{Exposure, infection and removal times respectively of the new premisses.}
#'         \item{age}{Age group category where the new premisse belongs to.}
#'         \item{infected_source}{Source of infection. 9999 for primary infection.}
#'         \item{value_indx}{Grid index in which lies the infection premisse in the raster.}
#'        }
#' @details The contact model developped here only consider a short range dispersal kernal as BBTV (the pathogen considered) only spead
#'          at short range from the source (ref). Consequentely, we use an expoential dispersal parameter \eqn{\alpha} of the form:
#'          \deqn{f(r,\alpha)=\alpha exp(-\alpha r)}
#' @seealso  \code{\link{E_to_I}}, \code{\link{func_time_beta}}.
#' @examples
#' data(bbtv)
#' attach(bbtv)
#' Dat<- bbtv[,c(2:6,8,10)]     # Coonsidering the essential part of the data
#' Dat1<-subset(Dat,Dat$latitude> -27.3 & Dat$processedbananas%in%c("P&I","P", "NI") )  # data up in queensland (noth of brisbane)
#' Dat1$treatmentdate[is.na(Dat1$treatmentdate)]<- Dat1$date[is.na(Dat1$treatmentdate)] # When NA, consider removal date as
#' Dat1$detection<-as.numeric(difftime(as.Date(Dat1$date), as.Date("2011/01/01"), unit="days")) + runif(1)
#' Dat1$removal<-as.numeric(difftime(as.Date(Dat1$treatmentdate), as.Date("2011/01/01"), unit="days"))
#' Datt<-Dat1[,c(3,4,7,8,2)]
#' Datt$detection=Datt$detection + runif(nrow(Datt))
#' Datt=Datt[with(Datt,order(Datt$detection)),]
#'
#'
#' # Australian reference system
#' sp::coordinates(Datt) <- c("longitude", "latitude")
#' sp::proj4string(Datt) <- sp::CRS("+init=epsg:4326")
#' australianCRS <- sp::CRS("+init=epsg:3577")
#'
#' pointsinaustraliangrid = sp::spTransform(Datt,australianCRS)
#'
#' # Raster
#' rast <- raster::raster()
#' raster::extent(rast) <- raster::extent(pointsinaustraliangrid) # Set same extent
#'
#' raster::res(rast)=5000 # Set resolution
#'
#' # And then ... rasterize it! This creates a grid version
#' # of your points using the cells of rast,
#'
#' rast2 <- raster::rasterize(pointsinaustraliangrid, rast, 1, fun=sum)
#'
#' # Extract infos on the grid
#'
#'
#' n_row_grid=nrow_grid=raster::nrow(rast)
#' n_col_grid=ncol_grid=raster::ncol(rast)
#' grid_size=raster::res(rast)[1]     # Resolution
#'
#' n_line=(nrow_grid+1) + (ncol_grid +1)  # Number of grid  lines
#'
#' x_min=raster::xmin(rast)  # min max of the bounding box
#' x_max=raster::xmax(rast)
#'
#' y_min=raster::ymin(rast)
#' y_max=raster::ymax(rast)
#'
#' da=as.data.frame(pointsinaustraliangrid)
#'
#' pop_per_grid=raster::values(rast2)
#' pop_per_grid[is.na(pop_per_grid)]=0
#' mat=matrix(pop_per_grid,nrow = nrow_grid, byrow = TRUE)
#' pop_grid=apply(mat,2,rev)     # population per grid
#'
#' # Structure of the grid
#' x=seq(x_min,x_max,grid_size)
#' y=seq(y_min,y_max,grid_size)
#'
#' grid_lines=array(0,c(n_line,6))
#' for(i in 1:n_line){
#'   if(i<=(nrow_grid +1)){
#'     grid_lines[i,]=c(i,1,x[1],y[i],x[length(x)],y[i])
#'   }
#'   else{
#'     grid_lines[i,]=c(i,2,x[i-length(y)],y[1],x[i-length(y)],y[length(y)])
#'   }
#' }
#'
#' grid_lines=as.data.frame(grid_lines)
#' colnames(grid_lines)<- c("indx","orient_line","coor_x_1","coor_y_1","coor_x_2","coor_y_2")
#' circle_x=2022230
#' circle_y=-3123109
#' r=10000
#'# Simulation with exponential kernel
#' alpha<- 30; beta<- 0.012; epsilon<- 0.02; omega<- 0.12; mu_lat<- 30; var_lat<- 20; t0<- 0; c<- 20;
#' param=data.frame(alpha=alpha, beta=beta, epsilon=epsilon, omega=omega, mu_lat=mu_lat, var_lat=var_lat, t0=t0, c=c)
#'
#' Simulate_contact_model(param, grid_lines, pop_grid)
#'
#' detach(bbtv)
#'
#' @export
Simulate_contact_model<- function(param, grid_lines, pop_grid, grid_size=5000, age_level=c(1,1),age_dist=c(1,0), m_start=1, t_max=118, t_intervention=365, EI_model=1){



  #Set parameters
  epsilon <- param$epsilon
  beta <- param$beta
  alpha <- param$alpha
  mu_lat <- param$mu_lat
  var_lat <- param$var_lat
  omega <- param$omega
  c<- param$c
  t0<- param$t0

  beta_1<- beta;
  beta_2<- beta;


  min_coor_x <- min(c(grid_lines$coor_x_1,grid_lines$coor_x_2)) # bounds of the region
  max_coor_x <- max(c(grid_lines$coor_x_1,grid_lines$coor_x_2))
  min_coor_y <- min(c(grid_lines$coor_y_1,grid_lines$coor_y_2))
  max_coor_y <- max(c(grid_lines$coor_y_1,grid_lines$coor_y_2))

  x_intervals <- seq(min_coor_x,max_coor_x,grid_size)
  y_intervals <- seq(min_coor_y,max_coor_y,grid_size)

  n_line <- max(grid_lines$indx) #+ 1 # number of lines
  n_row_grid <- nrow(pop_grid) # number of rows of grids
  n_col_grid <- ncol(pop_grid)  # number of cols of grids


  simulated_epi <- data.frame(k=numeric(0), coor_x=numeric(0), coor_y=numeric(0), t_e=numeric(0), t_i=numeric(0), t_r=numeric(0), age=numeric(0), infected_source=numeric(0), value_indx=numeric(0))

  ## initialize the index cases ##

  index_k <- 0:(m_start-1)
  index_coor_x <- stats::runif(m_start,min_coor_x,max_coor_x)
  index_coor_y <- stats::runif(m_start,min_coor_y,max_coor_y)
  dt<- stats::rexp(1)/epsilon    # Sellke

  index_t_e <- rep(t0+dt,m_start) # all assumed to infected at time=t0

  # Latent period depending on the model

  index_t_i <- index_t_e + E_to_I(EI_model,t0 , mu_lat, var_lat)
  index_t_r <- index_t_i + stats::rexp(m_start,rate=1/c)
  index_age <- sample(age_level,size=m_start,prob=age_dist,replace=T)
  index_source <- rep(9999,m_start) # all assumed to be infected by the background


  m=ceiling((index_coor_x-min_coor_x)/grid_size)-1
  n=ceiling((index_coor_y-min_coor_y)/grid_size)-1

  simulated_epi[1:m_start,] <- c(index_k, index_coor_x,index_coor_y,index_t_e,index_t_i,index_t_r,index_age,index_source,(n_row_grid-m-2)*n_col_grid + n)
  #####

  t_now <- min(simulated_epi$t_i)

  simulated_epi_sub <- subset(simulated_epi, simulated_epi$t_i<=t_now & simulated_epi$t_r>t_now) # those are currently infectious

  num_infection <- nrow(simulated_epi)


  t_next <- t_now # to start the while loop
  while(t_next<(t_max) | max(pop_grid)==0){

    ### simulate the timings, and the source, of next infection ###

    simulated_epi_sub <- subset(simulated_epi, simulated_epi$t_i<=t_now & simulated_epi$t_r>t_now) # those are currently infectious


  #Risk due to secondary infection
    if(nrow(simulated_epi_sub)>=1){
      beta_infectious <- sapply(simulated_epi_sub$age, FUN=beta_by_age, c(beta_1,beta_2))
      total_beta <- sum(beta_infectious)
      joint_I_R <- c(simulated_epi_sub$t_i,simulated_epi_sub$t_r)
      min_I_R <- min(joint_I_R[which(joint_I_R>t_now)])
     # print(total_beta)
    }

    #Risk due to primary infection
    if(nrow(simulated_epi_sub)<1){
      total_beta <- 0
    }


    t_next <- simulate_NHPP_next_event (t_now=t_now, t_intervention=t_intervention, sum_beta=total_beta, epsilon=epsilon, omega=omega, t_max=t_max) # simulate the next infection time using thinning algorithm

    #print(t_next)

    while(t_next>=min_I_R & t_next!=Inf){  # If next event is a removal

      t_now <- min_I_R
      simulated_epi_sub <- subset(simulated_epi, simulated_epi$t_i<=t_now & simulated_epi$t_r>t_now) # those are currently infectious


      if(nrow(simulated_epi_sub)>=1){

        beta_infectious <- sapply(simulated_epi_sub$age, FUN=beta_by_age, c(beta_1,beta_2))

        total_beta <- sum(beta_infectious)

        joint_I_R <- c(simulated_epi_sub$t_i,simulated_epi_sub$t_r)
        min_I_R <- min(joint_I_R[which(joint_I_R>t_now)])
        t_next <- simulate_NHPP_next_event (t_now=t_now,  t_intervention=t_intervention, sum_beta=total_beta, epsilon=epsilon, omega=omega, t_max=t_max)

      }

      if(nrow(simulated_epi_sub)<1){
        total_beta <- 0
        t_next <- simulate_NHPP_next_event (t_now=t_now,  t_intervention=t_intervention, sum_beta=total_beta, epsilon=epsilon, omega=omega, t_max=t_max)
        break # needed when consider background infection
      }

    } # end of while(t_next>=min_I_R)


    k <- num_infection + 1 - 1 # k=0,1,2...
    t_now <- t_next
    t_i_new <- t_now +  E_to_I(EI_model,t_now , mu_lat, var_lat)
    t_r_new <- t_i_new + stats::rexp(1,rate=1/c)


    if(nrow(simulated_epi_sub)>=1) source <- sample(c(9999,simulated_epi_sub$k),size=1, prob=c(epsilon,beta_infectious)) # 9999 = background

    if(nrow(simulated_epi_sub)<1) source <- 9999 # 9999 = background

    age <- sample(age_level,size=1,prob=age_dist)

    ### simulate the coordinates of the new infection (above) ###

    x_new = min_coor_x - 5 # to start the while loop
    y_new = min_coor_y - 5 # to start the while loop
    while(x_new<min_coor_x | x_new>max_coor_x | y_new<min_coor_y | y_new>max_coor_y){

      if (source!=9999){

        ru<- stats::runif(1)
        if(ru<.5)
        r <- stats::rexp(1,rate=alpha)

        set_points <- circle_line_intersections (circle_x=simulated_epi$coor_x[source+1],circle_y=simulated_epi$coor_y[source+1], r,  n_line=n_line, grid_lines=grid_lines)

        n_set_points = nrow(set_points)

        if (n_set_points>=1) {
          # show(set_points)
          # show(source)
          # show(simulated_epi)
          arcs <- func_arcs_attributes(set_points, pop_grid, r, min_coor_x, min_coor_y, grid_size, n_row_grid, n_col_grid)
          arcs$mass <- arcs$dens*arcs$len_arc
          arcs <- arcs[order(arcs$theta_abs),]

          sum_arcs_den <- sum(arcs$dens)

          if (sum_arcs_den>0){
            k_segment <- sample(1:nrow(arcs),size=1,prob=arcs$mass) # decide which segment the new infection would lie
            theta_within_segment <- stats::runif(1, min=0, max=arcs$theta_abs[k_segment]) # uniformly draws a point within the segment chosen above

            if (k_segment==1) theta_from_y_eq_0 <- theta_within_segment # the measure of theta from y=0
            if (k_segment!=1) theta_from_y_eq_0 <- sum(arcs$theta_abs[1:(k_segment-1)]) + theta_within_segment # the measure of theta from y=0

            x_new <- simulated_epi$coor_x[source+1] + r*cos(theta_from_y_eq_0) # the coor_x of the new infection
            y_new <- simulated_epi$coor_y[source+1] + r*sin(theta_from_y_eq_0) # the coor_x of the new infection

            n<- arcs$n_th_col
            m<- arcs$m_th_row
            if(max(pop_grid)>1){  # Density instead of sum
              pop_grid[m+1,n+1]<- pop_grid[m+1,n+1] -1/grid_size # Rememeber to divide by the resolution for the density
            }
            else{
              pop_grid[m+1,n+1]<- pop_grid[m+1,n+1] -1
            }


          }

          ###
          if (sum_arcs_den==0){
            theta_from_y_eq_0 <-  stats::runif(1, min=0,max=(2*pi)) # draw theta uniformly between [0,2pi]
            x_new <- simulated_epi$coor_x[source+1] + r*cos(theta_from_y_eq_0) # the coor_x of the new infection
            y_new <- simulated_epi$coor_y[source+1] + r*sin(theta_from_y_eq_0) # the coor_x of the new infection

            m=ceiling((x_new-min_coor_x)/grid_size)-1
            n=ceiling((y_new-min_coor_y)/grid_size)-1
          }
          ###


        }


        if (n_set_points<1) {
          theta_from_y_eq_0 <-  stats::runif(1, min=0,max=(2*pi)) # draw theta uniformly between [0,2pi]

          x_new <- simulated_epi$coor_x[source+1] + r*cos(theta_from_y_eq_0) # the coor_x of the new infection
          y_new <- simulated_epi$coor_y[source+1] + r*sin(theta_from_y_eq_0) # the coor_x of the new infection

          m=ceiling((x_new-min_coor_x)/grid_size)-1
          n=ceiling((y_new-min_coor_y)/grid_size)-1

          if(pop_grid[n+1,m+1]>0){
            if(max(pop_grid)>1){  # Density instead of sum
              pop_grid[n+1,m+1]<- pop_grid[n+1,m+1] -1/grid_size # Rememeber to divide by the resolution for the density
            }
            else{
              pop_grid[n+1,m+1]<- pop_grid[n+1,m+1] -1
            }
          }

        }


      }

      if(source==9999){

        k_grid <- sample(1:length(as.numeric(pop_grid)),size=1, prob=as.numeric(pop_grid))
        m_grid <- k_grid%%nrow(pop_grid) # the mth row of the grids
        if(m_grid==0) m_grid <- nrow(pop_grid)
        n_grid <- ceiling(k_grid/nrow(pop_grid))  # nth column ..
        x_new <- stats::runif(1,min=x_intervals[n_grid],max=x_intervals[n_grid+1]) # random lands at a point in the grid selected
        y_new <- stats::runif(1,min=y_intervals[m_grid],max=y_intervals[m_grid+1])

        n=n_grid-1
        m=m_grid-1
        if(pop_grid[m_grid,n_grid]>0){
        if(max(pop_grid)>1){  # Density instead of sum
          pop_grid[m_grid,n_grid]<- pop_grid[m_grid,n_grid] -1/grid_size # Rememeber to divide by the resolution for the density
        }
        else{
          pop_grid[m_grid,n_grid]<- pop_grid[m_grid,n_grid] -1
        }
        }

      }



    } # end while(x_new<min_coor_x | x_new>max_coor_x | y_new<min_coor_y | y_new>max_coor_y){


    simulated_epi[k+1,] <- c(k, x_new, y_new, t_now, t_i_new, t_r_new, age, source, (n_row_grid-m-2)*n_col_grid + n)
    ####

    num_infection <- num_infection + 1

  } # end of while(t_next<t_max)


return(simulated_epi)
}
