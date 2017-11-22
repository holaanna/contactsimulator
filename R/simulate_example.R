data(bbtv)
attach(bbtv)
Dat<- bbtv[,c(2:6,8,10)]     # Coonsidering the essential part of the data
Dat1<-subset(Dat,Dat$latitude> -27.3 & Dat$processedbananas%in%c("P&I","P", "NI") )  # data up in queensland (noth of brisbane)
Dat1$treatmentdate[is.na(Dat1$treatmentdate)]<- Dat1$date[is.na(Dat1$treatmentdate)] # When NA, consider removal date as
Dat1$detection<-as.numeric(difftime(as.Date(Dat1$date), as.Date("2011/01/01"), unit="days")) + runif(1)
Dat1$removal<-as.numeric(difftime(as.Date(Dat1$treatmentdate), as.Date("2011/01/01"), unit="days"))
Datt<-Dat1[,c(3,4,7,8,2)]
Datt$detection=Datt$detection + runif(nrow(Datt))
Datt=Datt[with(Datt,order(Datt$detection)),]


# Australian reference system
sp::coordinates(Datt) <- c("longitude", "latitude")
sp::proj4string(Datt) <- sp::CRS("+init=epsg:4326")
australianCRS <- sp::CRS("+init=epsg:3577")

pointsinaustraliangrid = sp::spTransform(Datt,australianCRS)

pointsinaustraliangridextent = raster::extent(pointsinaustraliangrid)

pointsinaustraliangridextent@xmin=floor(as.numeric(pointsinaustraliangridextent@xmin)/5000)*5000
pointsinaustraliangridextent@ymin=floor(as.numeric(pointsinaustraliangridextent@ymin)/5000)*5000
pointsinaustraliangridextent@xmax=ceiling(as.numeric(pointsinaustraliangridextent@xmax)/5000)*5000
pointsinaustraliangridextent@ymax=ceiling(as.numeric(pointsinaustraliangridextent@ymax)/5000)*5000


#Raster
rast <- raster::raster()
raster::extent(rast) <- raster::extent(pointsinaustraliangrid) # Set same extent

raster::res(rast)=5000 #Set resolution

# And then ... rasterize it! This creates a grid version
# of your points using the cells of rast,

rast2 <- raster::rasterize(pointsinaustraliangrid, rast, 1, fun=sum)

# Extract infos on the grid


n_row_grid=nrow_grid=raster::nrow(rast)
n_col_grid=ncol_grid=raster::ncol(rast)
grid_size=raster::res(rast)[1]     # Resolution

n_line=(nrow_grid+1) + (ncol_grid +1)  # Number of grid  lines

x_min=raster::xmin(rast)  # min max of the bounding box
x_max=raster::xmax(rast)

y_min=raster::ymin(rast)
y_max=raster::ymax(rast)

da=as.data.frame(pointsinaustraliangrid)

pop_per_grid=raster::values(rast2)
pop_per_grid[is.na(pop_per_grid)]=0
mat=matrix(pop_per_grid,nrow = nrow_grid, byrow = TRUE)
pop_grid=apply(mat,2,rev)     # population per grid

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
circle_x=2022230
circle_y=-3123109
r=10000

detach(bbtv)

