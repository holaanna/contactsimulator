\dontrun{
library(raster)
library(tidyverse)
library(parallel)
library(pracma)
library(contactsimulator)
f<- system.file("external/rast_newr_NSW.tif", package="contactsimulator")
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
}
