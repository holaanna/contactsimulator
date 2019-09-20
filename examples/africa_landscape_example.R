library(raster)
library(tidyverse)
library(parallel)
library(pracma)
library(contactsimulator)

# Create a landscape
malawi<- raster(xmn=0, xmx=10000, ymn=0, ymx=10000)
res(malawi)<- 70
values(malawi)<- 0

rast<- malawi
size<- grid_size <- raster::res(rast)[1]
n_row_grid=nrow_grid=raster::nrow(rast)
n_col_grid=ncol_grid=raster::ncol(rast)
grid_size=raster::res(rast)[1]     # Resolution

n_line=(nrow_grid+1) + (ncol_grid +1)  # Number of grid  lines

x_min=raster::xmin(rast)  # min max of the bounding box
x_max=raster::xmax(rast)

y_min=raster::ymin(rast)
y_max=raster::ymax(rast)

pop_per_grid=round(raster::values(rast)*size^2)
pop_per_grid[is.na(pop_per_grid)]=0
mat=matrix(pop_per_grid,nrow = nrow_grid, byrow = TRUE)
pop_grid=flipdim(mat)     # population per grid

# Structure of the grid
x=seq(x_min,x_max,grid_size)
y=seq(y_min,y_max,grid_size)

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

t_obs <- seq(30,4000,30)
row<- sample(1:length(x),1)
col<- sample(1:length(y),1)
ini<- data.frame(x=runif(1,x[col],x[col+1]), y=runif(1,y[row],y[row+1]),t_e=0,t_i=0,typ=0,row=row-1,col=col-1)
ini<- data.frame(x=c(35,35), y=c(10005,10007),t_e=0,t_i=0,typ=0,row=142,col=0)

xx<- xFromCol(rast,1:nrow(pop_grid))
yy<- yFromRow(rast,1:nrow(pop_grid))
idx<- data.frame(x=xx[1]-res(rast)[1]/2-0.5,y=yy)
idx<- rbind(idx,data.frame(x=xx,y=rep(yy[1]+res(rast)[1]/2+0.5,length(xx))))
idx<- rbind(idx,data.frame(x=xx,y=rep(yy[length(xx)]-res(rast)[1]/2-0.5,length(xx))))
idx<- rbind(idx,data.frame(x=rep(xx[length(xx)]+res(rast)[1]/2+0.5,length(xx)),y=yy))
row<- c(1:nrow(pop_grid),rep(0,length(xx)),rep(nrow(pop_grid),length(xx)),1:nrow(pop_grid))
col<- c(rep(0,length(xx)),1:nrow(pop_grid),1:nrow(pop_grid),rep(nrow(pop_grid),length(xx)))
ini<- data.frame(idx,t_e=0,t_i=0,typ=1, age=0,row=nrow(pop_grid)-row+1,col=col)

t_max<- max(t_obs)+1
pop_grid[row,col]<- 700
param <- data.frame(epsilon=0,beta_1=0.42846, beta_2=0,c=30, delta=0.713969, b1=0.798782,alpha1=1/0.00189890,alpha2=100,t0=0,omega=0.0172142,mu_lat=0.062,var_lat=0.903)
param$gama=0.5
xx=Simulate_contact_control_LER_africa(rast,param=param, grid_lines=grid_lines, pop_grid=pop_grid,t_max = t_max, EI_model = 1,kern_model = 5,t_obs = t_obs,grid_size = grid_size,m_start = nrow(ini),nb=2,ini = ini,leav = c(2,6),plant_proc=1:length(values(rast)))

xx[[2]]%>%ggplot(aes(seas,n_suck)) + geom_point() + geom_line() +
  scale_x_continuous(breaks=xx[[2]]$seas,labels=paste0("season ", 1:length(xx[[2]]$seas)))+
  xlab("season") + ylab("Nomber of suckers")

xx[[2]]%>%ggplot(aes(seas,n_remv)) + geom_point() + geom_line()+
  scale_x_continuous(breaks=xx[[2]]$seas,labels=paste0("season ", 1:length(xx[[2]]$seas)))+
  xlab("season") + ylab("Nomber removed")

xx[[2]]%>%ggplot(aes(seas,n_plantation)) + geom_point() + geom_line()+
  scale_x_continuous(breaks=xx[[2]]$seas,labels=paste0("season ", 1:length(xx[[2]]$seas)))+
  xlab("season") + ylab("Nomber of plantation")



