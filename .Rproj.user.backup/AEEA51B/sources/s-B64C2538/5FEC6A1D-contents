#' Simulation of the epidemic process
#'
#' Epidemic simulation using the contact type model.
#'
#'\code{Simulate_contact_model} provide the simulation of the epidemic process and the
#'       maximum distance the wave can travel at each particular obaservation date.
#'
#' @param param Indicating a data frame containing a vector of parameters including:
#'       \describe{
#'         \item{epsion}{The primary infection rate. See \code{\link{func_time_beta}}}
#'         \item{beta_0}{Baseline or average transmission rate. See \code{\link{func_time_beta}}}
#'         \item{beta_1}{Amplitude of the seasonality. See \code{\link{func_time_beta}}}
#'         \item{alpha1,alpha2}{The dispersal kernel parameters.}
#'         \item{mu_lat,var_lat}{mean and variance of the latent period. See \code{\link{E_to_I}} for details.}
#'         \item{t0}{Time at which the primary source became active}.
#'         \item{omega}{Period of the forcing. See \code{\link{func_time_beta}}}
#'         \item{gama}{The mean proportion of short range dispersal events.}.
#'        }
#' @inheritParams func_arcs_attributes
#' @param age_level,age_dist Vectors of age level and the propportion of each age group respectively. See details.
#' @param m_start The size of initial cases. Default is 1.
#' @param t_max  Final observation time.
#' @param t_intervention Start of the intervention if any.
#' @param EI_model Take integer values to specify the type of model used for the latent period. See \code{\link{E_to_I}}
#' @param kern_model Take integer values to specify the type of dispersal kernel used. See \code{\link{Samp_dis}}
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
#' @references
#' \insertRef{KR08}{contactsimulator}
#' \insertRef{Mee11}{contactsimulator}
#' @examples
#' data(bbtv)
#' attach(bbtv)
#' Dat<- bbtv[,c("longitude","latitude","BBTV","inspectiondate","leavesinfected","treatmentdate","location")]
#' Dat1<-subset(Dat,Dat$latitude> -27.4698 & Dat$BBTV%in%c("P&I","P", "NI") & difftime(as.Date(Dat$inspectiondate), as.Date("2010/01/01"), unit="days")>=0)  # data up in queensland
#' Dat1$treatmentdate[is.na(Dat1$treatmentdate)]<- Dat1$inspectiondate[is.na(Dat1$treatmentdate)]
#' Dat1$detection<-as.numeric(difftime(as.Date(Dat1$inspectiondate), as.Date("2010/01/01"), unit="days"))
#' Dat1$removal<-as.numeric(difftime(as.Date(Dat1$treatmentdate), as.Date("2010/01/01"), unit="days"))
#' Dat1$removal[which(Dat1$removal<0)]<- Dat1$detection[which(Dat1$removal<0)]
#' Datt<-Dat1[,c("longitude","latitude","BBTV","leavesinfected","detection","removal")]
#'
#'
#' Datt[which(Datt$leavesinfected=="LOTS"),"leavesinfected"]<- 45
#' Datt[which(Datt$leavesinfected=="1,2,4"),"leavesinfected"]<- 2.3
#' Datt[which(Datt$leavesinfected=="'3"),"leavesinfected"]<- 3
#' Datt[which(Datt$leavesinfected=="2 +bunch"),"leavesinfected"]<- 2
#' Datt[which(Datt$leavesinfected=="3 +bunch"),"leavesinfected"]<- 3
#' Datt[which(Datt$leavesinfected=="4+BUNCH"),"leavesinfected"]<- 4
#' Datt[which(Datt$leavesinfected=="avg 3.2"),"leavesinfected"]<- 3.2
#' Datt[which(Datt$leavesinfected=="1-6, avg 3.5"),"leavesinfected"]<- 3.5
#' Datt[which(Datt$leavesinfected=="all"),"leavesinfected"]<- 45
#'
#'
#' leav=sapply(Datt[,"leavesinfected"],function(x){
#'   gsub("all/","",x)
#' })
#'
#' leav=sapply(leav,function(x){
#'   gsub("/all","",x)
#' })
#'
#' leav[grepl("[+]",leav)]<- 45  # Assuming 45 leaves on a plant
#'
#' Datt$leavesinfected<- leav
#'
#' Datt=Datt[with(Datt,order(Datt$detection)),]
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
#'  size<- raster::res(rast)[1]
# Adding column at the top or bottom of the grid if raster leaves points out
#' dif=(raster::xmax(pointsinaustraliangrid)-raster::xmin(pointsinaustraliangrid))/size
#' cei= abs(ceiling(dif))
#'
#' if(cei!=dif){
#'   if(raster::xmax(rast)!=raster::xmax(pointsinaustraliangrid)){
#'     raster::xmax(rast)<- raster::xmin(rast) + size*cei
#'   }
#'   if(raster::xmin(rast)!=raster::xmin(pointsinaustraliangrid)){
#'     raster::xmin(rast)<- raster::xmax(rast) - size*cei
#'   }
#'
#' }
#'
#' # Adding row at the top or bottom of the grid if raster leaves points out
#'
#' dif1=(raster::ymax(pointsinaustraliangrid)-raster::ymin(pointsinaustraliangrid))/size
#' cei1= abs(ceiling(dif1))
#'
#'
#' if(cei1!=dif1){
#'   if(raster::ymax(rast)!=raster::ymax(pointsinaustraliangrid)){
#'     raster::ymax(rast)<- raster::ymin(rast) + size*cei1
#'   }
#'   if(raster::ymin(rast)!=raster::ymin(pointsinaustraliangrid)){
#'     raster::ymin(rast)<- raster::ymax(rast) - size*cei1
#'   }
#'
#' }
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
#' alpha<- 30; beta<- 0.012; epsilon<- 0.02; omega<- 0.12; mu_lat<- 30; var_lat<- 20; t0<- 0; c<- 20; b1<-0; gama<- 0.5
#' param=data.frame(alpha1=alpha, alpha2=alpha, beta=beta, epsilon=epsilon, omega=omega, mu_lat=mu_lat, var_lat=var_lat, t0=t0, c=c, b1=b1, gama=gama)
#'
#' Simulate_contact_model(param, grid_lines, pop_grid)
#'
#' detach(bbtv)
#'
#' @export
Simulate_contact_model<- function(param, grid_lines, pop_grid, grid_size=5000, age_level=c(1,1),age_dist=c(1,0), m_start=1, t_max=118, t_intervention=365, EI_model=1, kern_model=4){

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

  simulated_epi <- data.frame(k=numeric(0), coor_x=numeric(0), coor_y=numeric(0), t_e=numeric(0), t_i=numeric(0), t_r=numeric(0), age=numeric(0), infected_source=numeric(0), row=numeric(0), col=numeric(0))

  ## initialize the index cases ##

  index_k <- 0:(m_start-1)
  n_indx<- m_indx<- index_coor_x <- index_coor_y<- numeric(m_start)
  for(i in 1:m_start){
    k_grid <- sample(1:length(as.numeric(pop_grid)),size=1, prob=as.numeric(pop_grid))
    m_grid <- k_grid%%nrow(pop_grid) # the mth row of the grids
    if(m_grid==0) m_grid <- nrow(pop_grid)
    n_grid <- ceiling(k_grid/nrow(pop_grid))  # nth column ..
    index_coor_x[i] <- stats::runif(1,min=x_intervals[n_grid],max=x_intervals[n_grid+1]) # random lands at a point in the grid selected
    index_coor_y[i] <- stats::runif(1,min=y_intervals[m_grid],max=y_intervals[m_grid+1])

    n_indx[i]=n_grid
    m_indx[i]=m_grid

  }

  # index_coor_x <- stats::runif(m_start,min_coor_x,max_coor_x)
  # index_coor_y <- stats::runif(m_start,min_coor_y,max_coor_y)
  dt<- stats::rexp(1)/epsilon    # Sellke

  index_t_e <- rep(t0,m_start) # all assumed to infected at time=t0

  # Latent period depending on the model

  index_t_i <- index_t_e + E_to_I(EI_model,t0 , mu_lat, var_lat)
  index_t_r <- index_t_i + stats::rexp(m_start,rate=1/c)
  index_age <- sample(age_level,size=m_start,prob=age_dist,replace=T)
  index_source <- rep(9999,m_start) # all assumed to be infected by the background


  # m=ceiling((index_coor_x-min_coor_x)/grid_size)-1
  # n=ceiling((index_coor_y-min_coor_y)/grid_size)-1

  simulated_epi[1:m_start,] <- c(index_k, index_coor_x,index_coor_y,index_t_e,index_t_i,index_t_r,index_age,index_source,m_indx,n_indx)
  #####

  t_now <- min(simulated_epi$t_i)

  simulated_epi_sub <- subset(simulated_epi, simulated_epi$t_i<=t_now & simulated_epi$t_r>t_now) # those are currently infectious

  num_infection <- nrow(simulated_epi)


  t_next <- t_now # to start the while loop

  while(t_next<(t_max) & max(pop_grid)>0){
# show(t_next)
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
      min_I_R <- t_now
    }


    t_next <- simulate_NHPP_next_event (t_now=t_now, t_intervention=t_intervention, sum_beta=total_beta, epsilon=epsilon, omega=omega, b1=b1, t_max=t_max) # simulate the next infection time using thinning algorithm
    #print(c(t_next,t_now,1))
    #print(c(t_next,total_beta,omega,t_max,t_now,t_intervention))

    while(t_next>=min_I_R & t_next!=Inf){  # If next event is a removal

      t_now <- min_I_R
      simulated_epi_sub <- subset(simulated_epi, simulated_epi$t_i<=t_now & simulated_epi$t_r>t_now) # those are currently infectious


      if(nrow(simulated_epi_sub)>=1){

        beta_infectious <- sapply(simulated_epi_sub$age, FUN=beta_by_age, c(beta_1,beta_2))

        total_beta <- sum(beta_infectious)

        joint_I_R <- c(simulated_epi_sub$t_i,simulated_epi_sub$t_r)
        min_I_R <- min(joint_I_R[which(joint_I_R>t_now)])
        t_next <- simulate_NHPP_next_event (t_now=t_now,  t_intervention=t_intervention, sum_beta=total_beta, epsilon=epsilon, omega=omega, b1=b1,t_max=t_max)
        #print(c(t_next,t_now,2))
      }

      if(nrow(simulated_epi_sub)<1){
        total_beta <- 0
        t_next <- simulate_NHPP_next_event (t_now=t_now,  t_intervention=t_intervention, sum_beta=total_beta, epsilon=epsilon, omega=omega, b1=b1, t_max=t_max)
        #print(c(t_next,t_now,3))
        break # needed when consider background infection
      }

    } # end of while(t_next>=min_I_R)

    #print(c(t_next,t_now,k,nrow(simulated_epi_sub)))

    k <- num_infection + 1 - 1 # k=0,1,2...
    t_now <- t_next
    t_i_new <- t_now +  E_to_I(EI_model,t_now , mu_lat, var_lat)
    t_r_new <- t_i_new + stats::rexp(1,rate=1/c)


    if(nrow(simulated_epi_sub)>=1) source <- sample(c(9999,simulated_epi_sub$k),size=1, prob=c(epsilon,beta_infectious)) # 9999 = background

    if(nrow(simulated_epi_sub)<1) source <- 9999 # 9999 = background

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
    m_1=n_1=1000
    pop_grid_before = pop_grid
    while(x_new<min_coor_x | x_new>max_coor_x | y_new<min_coor_y | y_new>max_coor_y ){
      pop_grid = pop_grid_before
      if (source!=9999){
      # show(kern_model)
      # show(ru)
      # show(alpha1)
      # show(alpha2)

        r <- abs(Samp_dis (kern_model,ru, alpha1, alpha2))

        set_points <- circle_line_intersections (circle_x=simulated_epi$coor_x[source+1],circle_y=simulated_epi$coor_y[source+1], r,  n_line=n_line, grid_lines=grid_lines)

        n_set_points = nrow(set_points)
        #print(c(r,n_set_points))
        if (n_set_points>=1) {
          # show(set_points)
          # show(source)
          # show(simulated_epi)
          arcs <- func_arcs_attributes(set_points, pop_grid, r, min_coor_x, max_coor_y, grid_size, n_row_grid, n_col_grid)

          arcs$mass <- arcs$dens*arcs$len_arc

          arcs$theta=set_points$theta

          #arcs <- arcs[order(arcs$theta_abs),]

          sum_arcs_den <- sum(arcs$dens)
             # print(c(sum_arcs_den,n_line))
             # print(c(simulated_epi$coor_x[source+1],simulated_epi$coor_y[source+1]))
            # print(as.data.frame(arcs$n_th_col,arcs$m_th_row))
            # print(sum_arcs_den)
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
            theta_from_y_eq_0 <-  stats::runif(1, min=0,max=(2*pi)) # draw theta uniformly between [0,2pi]
            x_new <- simulated_epi$coor_x[source+1] + r*cos(theta_from_y_eq_0) # the coor_x of the new infection
            y_new <- simulated_epi$coor_y[source+1] + r*sin(theta_from_y_eq_0) # the coor_x of the new infection

            n=ceiling((x_new-min_coor_x)/grid_size)
            m=ceiling((-y_new+max_coor_y)/grid_size)

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



    } # end while(x_new<min_coor_x | x_new>max_coor_x | y_new<min_coor_y | y_new>max_coor_y){
    #print(c(t_next,t_now,3))
    if(is.infinite(t_next)){
      break
    }
    else{
      simulated_epi[k+1,] <- c(k, x_new, y_new, t_now, t_i_new, t_r_new, age, source, m,n)

    }
    ####

    num_infection <- num_infection + 1

  } # end of while(t_next<t_max)


return(simulated_epi)
}

