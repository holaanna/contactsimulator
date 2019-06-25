\dontrun{
library(raster)
library(tidyverse)
library(parallel)
library(pracma)
library(contactsimulator)
  size=500
f<- system.file("external/rast_newr_NSW.tif", package="contactsimulator")
f<- system.file(paste0("external/rast_NSW_",size,".tif"), package="contactsimulator")
rast<- raster(f)
names(rast)<- "layer"
size<- raster::res(rast)[1]
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

# Individual raster
f<- system.file(paste0("external/rast_farms_NSW_",size,".tif"), package="contactsimulator")
f_rast<- raster(f)
pop=round(raster::values(f_rast)*size^2)
pop[is.na(pop)]=0
mat=matrix(pop,nrow = nrow_grid, byrow = TRUE)
f_rast=flipdim(mat)
y<- system.file(paste0("external/rast_backy_NSW_",size,".tif"), package="contactsimulator")
b_rast<- raster(y)
pop=round(raster::values(b_rast)*size^2)
pop[is.na(pop)]=0
mat=matrix(pop,nrow = nrow_grid, byrow = TRUE)
b_rast=flipdim(mat)

param<- data.frame(epsilon=10,beta=15,c=30,b1=0.3,alpha1=18,alpha2=1000,t0=0,omega=2*pi,mu_lat=0.062,var_lat=0.903,gama=0.5)
path="/run/user/1000/gvfs/sftp:host=login-cpu.hpc.cam.ac.uk/home/ha411/big_space/BBTV/Newrybar_partial/latent-leaf/"
path1="/run/user/1000/gvfs/sftp:host=login-cpu.hpc.cam.ac.uk/home/ha411/big_space/BBTV/"
dat=as.data.frame(read.table(paste0(path1,"NSW_8_1_19/expo-kernel/reso_500/Dout")))
dat=as.data.frame(read.table(paste0(path,"no-initial-infection/Dout")))

farm_pos_cat<- read.table(paste0(path1,"NSW_8_1_19/expo-kernel/farm_pos_cat.txt"),fill=T,header = T)
farm_pos_cat<- farm_pos_cat[with(farm_pos_cat,order(farm_pos_cat$cat)),]


M=read.table(paste0(path1,"NSW_8_1_19/expo-kernel/reso_500/cyclic_latent_omegapi/parameters_current.txt"),head=T)
M=read.table(paste0(path,"test-vague-prior/cyclic_latent_omegapi/parameters_current.txt"),head=T)
samp=sample(150000:nrow(M),1000);Param=data.frame(epsilon=M[samp,1],beta=M[samp,2],c=M[samp,5], delta=M[samp,6], b1=M[samp,7],alpha1=1/M[samp,8],alpha2=M[samp,9],t0=M[samp,10],omega=M[samp,11],mu_lat=M[samp,3],var_lat=M[samp,4])
Param$gama=0.5
ini<- data.frame(x=2069503,y=-3293246,t_e=0,t_i=0,typ=0,row=9,col=15)
ini=dat%>%filter(V5==min(V5))%>%dplyr::select(V2,V3,V4,V5,V10,V12,V13)%>%rename(x=V2,y=V3,t_e=V4,t_i=V5,typ=V10,row=V12,col=V13)
#ini<- ini%>%mutate(row=row+1,col=col+1)
t_obs<- read.table(paste0(path1,"NSW_8_1_19/expo-kernel/reso_500/obs_time1"),fill = T)[,1]
t_max=max(t_obs)+1
L=list()
sim=numeric(1000)
L<- lapply(1:1000,function(x){
  print(x)
  param=Param[x,]
  # param$delta=.6
  xx=Simulate_contact_control_LER_farm(param=param, grid_lines=grid_lines, pop_grid=pop_grid,t_max = t_max, EI_model = 1,kern_model = 5,t_obs = t_obs,grid_size = grid_size,m_start = 5,nb=2,ini = ini,leav = c(2,6))
  # xx=Simulate_contact_control_LER(param=param, grid_lines=grid_lines, pop_grid=pop_grid,t_max = t_max, EI_model = 1,kern_model = 5,t_obs = 2000,grid_size = grid_size,m_start = 5,nb=2,ini = ini)
  return(xx)
})


sim<- unlist(lapply(L,nrow))
hist(sim,breaks=100)
quantile(sim,c(0.025,0.5,0.975))

df<- L[[18]]
midx<- matrix(c(df$row,df$col),nrow = nrow(df))
df1=df%>%mutate(t_up=t_i,age=age-1,row=row-1,col=col-1,leaf=3,t_e=t_e,t_i=t_i,t_r=t_r)%>%
  dplyr::select(k,coor_x,coor_y,t_e,t_i,t_r,t_up,leaf,age,typ,infected_source,row, col)
df1$dens<- pop_grid[midx]
write.table(df1,"/run/user/1000/gvfs/sftp:host=login-cpu.hpc.cam.ac.uk/home/ha411/big_space/BBTV/Newrybar_partial/Simulation/Dout",row.names = FALSE)

df1<- read.table("/run/user/1000/gvfs/sftp:host=login-cpu.hpc.cam.ac.uk/home/ha411/big_space/BBTV/Newrybar_partial/Simulation/Dout",fill = T,header = T)
df1<- df1%>%mutate(t_e=365*t_e,t_i=365*t_i,t_r=365*t_r,t_up=365*t_up)

}
