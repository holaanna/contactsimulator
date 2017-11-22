#' Visualisation of the dynamic of the epidemic
#'
#' Visualisation of the process.
#'
#'\code{Plot_dyn} Provide a basic mode of visualising the spread of the epidemic
#' @param rast The raser object provided.
#' @inheritParams Simulate_contact_model
#' @seealso \code{\link[raster]{raster}}
#' @param k Indicate the rank of the simulation
#'
#' @return A raster plot with points representing infections premisses:
#'
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
#'# Simulation with exponential kernel
#' alpha<- 30; beta<- 0.012; epsilon<- 0.02; omega<- 0.12; mu_lat<- 30; var_lat<- 20; t0<- 0; c<- 20;
#' param=data.frame(alpha=alpha, beta=beta, epsilon=epsilon, omega=omega, mu_lat=mu_lat, var_lat=var_lat, t0=t0, c=c)
#' Plot_dyn(rast2,param,grid_lines, pop_grid,k=1)
#' @export
Plot_dyn<-function(rast,param, grid_lines, pop_grid, age_level=c(1,1),age_dist=c(1,0), m_start=1, t_max=118, t_intervention=365, EI_model=3,k){

    # Extract infos on the grid


  # n_row_grid=nrow_grid=raster::nrow(rast)
  # n_col_grid=ncol_grid=raster::ncol(rast)
   grid_size=raster::res(rast)[1]     # Resolution
  #
  # n_line=(nrow_grid+1) + (ncol_grid +1)  # Number of grid  lines
  #
  # x_min=raster::xmin(rast)  # min max of the bounding box
  # x_max=raster::xmax(rast)
  #
  # y_min=raster::ymin(rast)
  # y_max=raster::ymax(rast)
  #
  # da=as.data.frame(pointsinaustraliangrid)
  #
  # pop_per_grid=raster::values(rast)
  # pop_per_grid[is.na(pop_per_grid)]=0
  # mat=matrix(pop_per_grid,nrow = nrow_grid, byrow = TRUE )
  # pop_grid=apply(mat,2,rev)     # population per grid
  #
  # # Structure of the grid
  # x=seq(x_min,x_max,grid_size)
  # y=seq(y_min,y_max,grid_size)
  #
  # grid_lines=array(0,c(n_line,6))
  # for(i in 1:n_line){
  #   if(i<=(nrow_grid +1)){
  #     grid_lines[i,]=c(i,1,x[1],y[i],x[length(x)],y[i])
  #   }
  #   else{
  #     grid_lines[i,]=c(i,2,x[i-length(y)],y[1],x[i-length(y)],y[length(y)])
  #   }
  # }
  #
  # grid_lines=as.data.frame(grid_lines)
  # colnames(grid_lines)<- c("indx","orient_line","coor_x_1","coor_y_1","coor_x_2","coor_y_2")

  #simulation and plot
  simulated_epi<- Simulate_contact_model(param, grid_lines, pop_grid, grid_size, age_level,age_dist, m_start, t_max, t_intervention, EI_model)
  raster::plot(rast,asp=1,xlab=" ",ylab=" ",main=bquote(paste('Simulation'==.(k))), xaxt="n",yaxt="n")

  # for(x in seq(0,t_max,by=25)){
  #   datas=simulated_epi[simulated_epi[,"t_e"]<=x,2:3,drop=F]


    graphics::points(simulated_epi$coor_x,simulated_epi$coor_y,cex=.2,pch=19)
    # Sys.sleep(0.5)
  # }


}
