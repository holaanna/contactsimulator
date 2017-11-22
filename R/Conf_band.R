#' Visualisation of the crdible band
#'
#' Posterior credible band.
#'
#'\code{Conf_Band} Provide a 95% credible band.
#'
#' @param rast Raster object.
#' @inheritParams Plot_dyn
#' @param Param The 8 colums data frame with columns names:
#'      \describe{
#'         \item{epsion}{The primary infection rate. See \code{\link{func_time_beta}}}
#'         \item{beta}{Rate of infection. See \code{\link{func_time_beta}}}
#'         \item{alpha}{The dispersal kernel parameter.}
#'         \item{mu_lat,var_lat}{mean and variance of the latent period. See \code{\link{E_to_I}} for details.}
#'         \item{omega}{Parameter characterising the seasality. See \code{\link{func_time_beta}}}
#'        }
#' @param Dat_obs A 2 columns data frame wiht compoents:
#'          \describe{
#'            \item{times}{ The sequence of times at which the observations are made}
#'            \item{remov}{the number of reovals recorded at a given time}
#'          }
#' @seealso \code{\link[raster]{raster}}, \code{\link{Plot_dyn}}, \code{\link{Simulate_contact_model}}, \code{\link{func_time_beta}}, \code{\link{E_to_I}}
#' @note The number of rows in \code{Param} corresponds to the size of the simulation expected.
#' @return A 95% credible band for the obsrvation. By default it considers the removals.
#'
#' @examples
#' data(rast2)
#' data(Dat_obs)
#' data(Param)
#' #Conf_Band(rast2,Param,Dat_obs=Dat_obs)
#' @export
Conf_Band<-function(rast,Param, age_level=c(1,1),age_dist=c(1,0), m_start=1, t_max=118, t_intervention=365, EI_model=3, Dat_obs){
  #Extrct infos from the raster
  n_row_grid=nrow_grid=raster::nrow(rast)  # number of rows of grids
  n_col_grid=ncol_grid=raster::ncol(rast)  # number of cols of grids
  grid_size=raster::res(rast)[1]     # Resolution

  n_line=(nrow_grid+1) + (ncol_grid +1)  # Number of grid  lines

  min_coor_x<- x_min<- raster::xmin(rast)  # min max of the bounding box
  max_coor_x <- x_max<- raster::xmax(rast)

  min_coor_y <- y_min<- raster::ymin(rast)
  max_coor_y <- y_max<- raster::ymax(rast)

  da=as.data.frame(pointsinaustraliangrid)

  pop_per_grid=raster::values(rast)
  pop_per_grid[is.na(pop_per_grid)]=0
  mat=matrix(pop_per_grid,nrow = nrow_grid, byrow = TRUE )
  pop_grid=apply(mat,2,rev)     # population per grid

  # Structure of the grid
  x_intervals <- x<- seq(x_min,x_max,grid_size)
  y_intervals <- y<- seq(y_min,y_max,grid_size)

  grid_lines=array(0,c(n_line,6))
  for(i in 1:n_line){
    if(i<=(nrow_grid +1)){
      grid_lines[i,]=c(i,1,x[1],y[i],x[length(x)],y[i])
    }
    else{
      grid_lines[i,]=c(i,2,x[i-length(y)],y[1],x[i-length(y)],y[length(y)])
    }
  }

  grid_lines=as.data.frame(grid_lines)
  colnames(grid_lines)<- c("indx","orient_line","coor_x_1","coor_y_1","coor_x_2","coor_y_2")

  #Set parameters
  epsilon <- Param$epsilon
  beta <- Param$beta
  alpha <- Param$alpha
  mu_lat <- Param$mu_lat
  var_lat <- Param$var_lat
  omega <- Param$omega
  c<- Param$c
  t0<- Param$t0


  N_sim<- nrow(Param)

  times<- Dat_obs[,"times"]
  obs<- Dat_obs[,"remov"]
  rem_siz=array(0,c(N_sim,(length(times)-1)))
  for(i in 1:N_sim){
    times[1]<- t0[i]
    param=data.frame(alpha=1/alpha[i], beta=beta[i], epsilon=epsilon[i], omega=omega[i], mu_lat=mu_lat[i], var_lat=var_lat[i], t0=t0[i], c=c[i])
    simulated_epi<- Simulate_contact_model(param, grid_lines, pop_grid, age_level, age_dist, m_start, t_max, t_intervention, EI_model)
    for(j in 2:(length(times))){
       rem_siz[i,j-1]<- length(which(simulated_epi[,"t_r"]>times[j-1] & simulated_epi[,"t_r"]<=times[j]))
    }

  }


  yv=times[-1]
  mat=rem_siz
  qv=c(0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45)
  ny=length(yv)
  yy=c(yv,yv[ny:1])
  graphics::plot(yv,obs,type='n',xlim=c(0,times[length(times)]),ylim=c(0,2*max(mat)),xlab='Time (days)',ylab='Removed')
  graphics::points((Dat_obs[,"times"])[-1],(Dat_obs[,"remov"])[-1],type='p')

  for(i in 1:length(qv)){
    v1=yv*0
    v2=v1
    for(j in 1:ny){
      v1[j]=stats::quantile(mat[j,],qv[i])
      v2[j]=stats::quantile(mat[j,],1-qv[i])
    }
    vv=c(v1,v2[ny:1])
    graphics::polygon(yy,vv,border=NA,col=grDevices::rgb(0,(0.5-qv[i])*2,(0.5-qv[i])*2))
  }


}
