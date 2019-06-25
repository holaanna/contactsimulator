library(pracma)
library(raster)
library(plotrix)
library(dplyr)
raster1<- raster(ncol=10,nrow=10, xmx=10000, xmn=0, ymn=0, ymx=10000)
ind<- sample(1:ncell(raster1),10)
values(raster1)[ind]<- round(runif(length(ind),0,100))

raster2<- raster(ncol=10,nrow=10, xmx=10000, xmn=0, ymn=0, ymx=10000)
ind1<- sample(1:ncell(raster2),20)
ind3<- which(ind1%in%ind)
ind1<- c(ind1,ind[-ind3][1:5])
values(raster2)[ind1]<- round(runif(length(ind1),0,10))

raster3<- mosaic(raster1,raster2,fun=sum)

size<- raster::res(raster3)[1]
# Extract infos on the grid


n_row_grid=nrow_grid=raster::nrow(raster3)
n_col_grid=ncol_grid=raster::ncol(raster3)
grid_size=raster::res(raster3)[1]     # Resolution

n_line=(nrow_grid+1) + (ncol_grid +1)  # Number of grid  lines

x_min=raster::xmin(raster3)  # min max of the bounding box
x_max=raster::xmax(raster3)

y_min=raster::ymin(raster3)
y_max=raster::ymax(raster3)

pop_per_grid=round(raster::values(raster3))
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



f_rast<- raster1
b_rast<- raster2

pop=round(raster::values(f_rast))
pop[is.na(pop)]=0
mat=matrix(pop,nrow = nrow_grid, byrow = TRUE)
f_rast=flipdim(mat)

pop=round(raster::values(b_rast))
pop[is.na(pop)]=0
mat=matrix(pop,nrow = nrow_grid, byrow = TRUE)
b_rast=flipdim(mat)

indx<- which(f_rast>0,arr.ind = T)          # coordinates of cells containing farms
t_max<- 1000
farm_pos_cat<- data.frame(ro=indx[,1],co=indx[,2],cat="A",vis_int=360,tim_lst_pos=2*t_max,nb_round=1,sweep=1)  # Farm position and category in the grid
vis_int_per_cat<- data.frame(cat=LETTERS[1:5],vis_int=c(360,360,30,30,30))
# Simulation with exponential kernel
alpha<- 30; beta<- 0.012; epsilon<- 0.02; omega<- 0.12; mu_lat<- 30; var_lat<- 20; t0<- 0; c<- 20; b1<-0; gama<- 0.5
param=data.frame(alpha1=alpha, alpha2=alpha, beta=beta, epsilon=epsilon, omega=omega, mu_lat=mu_lat, var_lat=var_lat, t0=t0, c=c, b1=b1, gama=gama)




set_points<- circle_line_intersections(epi[735,2],epi[735,3],258.9893,n_line,grid_lines)

atr<- func_arcs_attributes(set_points,pop_grid,258.9893,x_min,y_min,grid_size,n_row_grid,n_col_grid)


plot(raster3)
draw.circle(epi[734,2],epi[734,3],258.9893)
points(4050,3510,pch=19)

M=Simulate_contact_control(f_rast,b_rast,farm_pos_cat,vis_int_per_cat,param, grid_lines, pop_grid,t_b=1500,t_max = 1000,rad = 1000)
plot(raster3)
draw.circle(4050,3510,2000)
id<- which(M[,11]==0)
points(M[id,c("coor_x","coor_y")],cex=.5,pch=19,col="red")
points(M[-id,c("coor_x","coor_y")],cex=.5,pch=19,col="blue")





# Combined raster
# f<- system.file("external/rast_SEQ.tif", package="contactsimulator")
size=500
#f<- paste0("/media/sf_D_DRIVE/BBTV_PROJECT/Data/rast_SEQ_",size,".tif")
f<- paste0("/media/sf_D_DRIVE/BBTV_PROJECT/Data/rast_NSW_",size,".tif")
rast<- raster(f)
size<- raster::res(rast)[1]
# Extract infos on the grid


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
# f<- paste0("/media/sf_D_DRIVE/BBTV_PROJECT/Data/rast_farms_SEQ_",size,".tif")
f<- paste0("/media/sf_D_DRIVE/BBTV_PROJECT/Data/rast_farms_NSW_",size,".tif")
f_rast<- raster(f)
pop=round(raster::values(f_rast)*size^2)
pop[is.na(pop)]=0
mat=matrix(pop,nrow = nrow_grid, byrow = TRUE)
f_rast=apply(mat,2,rev)
# y<- paste0("/media/sf_D_DRIVE/BBTV_PROJECT/Data/rast_backy_SEQ_",size,".tif")
y<- paste0("/media/sf_D_DRIVE/BBTV_PROJECT/Data/rast_backy_NSW_",size,".tif")
b_rast<- raster(y)
pop=round(raster::values(b_rast)*size^2)
pop[is.na(pop)]=0
mat=matrix(pop,nrow = nrow_grid, byrow = TRUE)
b_rast=apply(mat,2,rev)


# farm_pos_cat<- read.table("/run/user/1000/gvfs/sftp:host=login-cpu.hpc.cam.ac.uk/home/ha411/big_space/BBTV/SEQ_26_11_18/Cauchy-kernel/farm_pos_cat.txt",fill=T,header = T)
farm_pos_cat<- read.table("/run/user/1000/gvfs/sftp:host=login-cpu.hpc.cam.ac.uk/home/ha411/big_space/BBTV/NSW_26_11_18/Cauchy-kernel/farm_pos_cat.txt",fill=T,header = T)


vis_int_per_cat<- data.frame(cat=LETTERS[1:5],vis_int=c(360,360,30,30,30))

#============================================================================================#
#                       Generate plot for simulation using draws from posterior
#============================================================================================#


times<- seq(1334,4388,30)
#times<- c(times, 365+times)
#obs<- c(212,90,57)
times1<- c(0,times)
OBS<- array(0,c(1000,length(times)))

#Risk of each cell
#indx<- which(pop_grid>0,arr.ind = T)
Risk<- array(0,c(nrow(farm_pos_cat),length(times)+2))
Risk<- as.data.frame(Risk)
Risk[,1:2]<- farm_pos_cat[,1:2]

Dat1<- read.table("/run/user/1000/gvfs/sftp:host=login-cpu.hpc.cam.ac.uk/home/ha411/big_space/BBTV/Newrybar_partial/Dout",sep=" ",header =T,stringsAsFactors = F)
dat=as.data.frame(read.table("/run/user/1000/gvfs/sftp:host=login-cpu.hpc.cam.ac.uk/home/ha411/big_space/BBTV/Newrybar_partial/Dout"))
# pop_grid <- as.matrix(read.table("pop_grid",sep=" ",header =F,fill = T))
# grid_lines <- read.table("grid_line.txt",sep=" ",header =T,stringsAsFactors = F)


Sim1=Sim=array(0,c(1000,9))
Li=vector("list",length = 9)
M=read.table("/run/user/1000/gvfs/sftp:host=login-cpu.hpc.cam.ac.uk/home/ha411/big_space/BBTV/Newrybar_partial/parameters_current.txt",fill=T)
sour=read.table("/run/user/1000/gvfs/sftp:host=login-cpu.hpc.cam.ac.uk/home/ha411/big_space/BBTV/Newrybar_partial/infected_source_current.txt",fill=T)
prop=numeric(nrow(sour)-1)
for(k in 1:(nrow(sour)-1)){
  kk=ss=0
  for( i in 1:ncol(sour)){
    if(sour[k,i]!=9999){
      kk = kk +1
      if(all(as.numeric(dat[i,c(10,11)])==as.numeric(dat[sour[k,i]+1,c(10,11)]))){
        ss = ss +1
      }
    }
  }
  prop[k]=ss/(ss+kk)
}
samp=sample(100000:nrow(M),1);Param=data.frame(epsilon=M[samp,1],beta=M[samp,2],c=M[samp,5],b1=M[samp,6],alpha1=M[samp,7],alpha2=M[samp,8],t0=M[samp,9],omega=M[samp,10],mu_lat=M[samp,3],var_lat=M[samp,4])

epsilon <- Param$epsilon
beta <- Param$beta
alpha1 <- Param$alpha1
alpha2 <- Param$alpha2
mu_lat <- Param$mu_lat
var_lat <- Param$var_lat
omega <- Param$omega
c<- Param$c
t0<- Param$t0
b1 <- Param$b1
prop <- mean(prop)

#define the control structure
A_B<- c(360,540,720)
C_D<- c(30,90,180)
Rat<- c(0.3,1)
radius<- c(0.25,0.5)

#Input

t_b<-max(times) + 2
t_max <- min(times)
#t_max<- 365
# sweep proportion
grid_size<- size
t_obs<- min(times)+1
ini=dat%>%filter(V5>min(V5))%>%filter(V4<=min(V5))%>%select(V2,V3,V4,V5,V8,V10,V11)%>%rename(x=V2,y=V3,t_e=V4,t_i=V5,typ=V8,row=V10,col=V11)
for(i in 1:1){
  param=data.frame(alpha1=1/alpha1[i],alpha2=alpha2[i], beta=beta[i], b1=b1[i], epsilon=epsilon[i], omega=omega[i], mu_lat=mu_lat[i], var_lat=var_lat[i], t0=t0[i], c=c[i], gama=prop)
  epi<- Simulate_contact_control(f_rast,b_rast,farm_pos_cat,vis_int_per_cat,param, grid_lines, pop_grid,grid_size = 500,t_b=max(times)+2,t_max = min(times), EI_model = 3,kern_model = 4,rad = 1000,sweep_prop = c(.5,.5))
  # sim[i]<- length(which(epi$t_i<=118))
  # sim1[i]<- length(which(epi$t_r<=118))
  # Sim[i,3*(k-1)+l]<- length(which(epi$t_r<=118))
  # Sim1[i,3*(k-1)+l]<- length(which(epi$t_i<=118))
  for(j in 1:(length(times))){
    indx1<- which(epi[,"t_e"]<=times1[j+1] & epi[,"t_e"]>times1[j])
    OBS[i,j]<- length(indx1)
    # for (ll in 1:nrow(Risk)) {
    #   id<- nrow(subset(epi,epi[,"t_e"]<=times1[j+1] & epi[,"t_e"]>times[j] & row==Risk[ll,1] & col==Risk[ll,2]))
    #   if(length(id)==0){
    #     id=0
    #   }
    #   Risk[ll,j+2]<- id
    # }

  }
}

k<- 472 + 1
r<- 70.88975

set_points<- circle_line_intersections(epi[k,2],epi[k,3],r,n_line,grid_lines)
atr<- func_arcs_attributes(set_points,pop_grid,82.52887,x_min,y_min,grid_size,n_row_grid,n_col_grid)
plot(rast,xlim=c(epi[k,2]-5000,epi[k,2]+5000),ylim=c(epi[k,3]-5000,epi[k,3]+5000))
draw.circle(epi[k,2],epi[k,3],r)
