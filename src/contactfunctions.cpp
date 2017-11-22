#include "header.h"
//// [[Rcpp::interfaces(r, cpp)]]
/*======================================================================================================================================
                                                     Intensity of infection
 ======================================================================================================================================*/
//' Functional form of the intensity (beta) as a function of t.
//'
//'\code{func_time_beta} returns the value of the intensity as a function of t.
//'
//' @param t,t_intervention  Time of the contact and the introduction of the control respectively espressed in days.
//' @param sum_beta,epsilon   Total risk from infected premisses and the primary infection rate respectively.
//' @param omega             Introduses the effect of seasonality in the infection rate.
//'
//' @details The particular case of transmission rate used here is given by \deqn{\beta(t)= n(t)\beta + \epsilon)cos^2(\omega t) if t<t_interventioin}
//'          where n(t) is the size of potential sources.
//'
//' @return It returns the rate of infection in the non-homogeneous poisson process.
//'
//' @examples
//' func_time_beta(20,50,0.08,0.2,0.3)
//' @export
// [[Rcpp::export]]
double func_time_beta (const double& t, const double& t_intervention, const double& sum_beta, const double& epsilon,  const double& omega){


  double beta_t;
  int l=(t<t_intervention);
  switch(l){
  case 1:{

    beta_t =  (sum_beta + epsilon)*pow(cos(omega*t/365),2);
    //cout<<beta<<"\t"<<beta_t<<"\n";
    break;
  }

  case 0:{ // To be edited acccordinly with presence of control
    beta_t = beta_t;
    break;
  }

  }

  return(beta_t);

}

/*======================================================================================================================================
                                                  Simulation of the time of the next event of a non-homogeneous oisson process
 ======================================================================================================================================*/

//' Simulation of the time of the next event of a non-homogeneous oisson process .
//'
//'\code{simulate_NHPP_next_event} computes the next-event time given the current state at t.
//'
//' @param t_now,t_intervention,t_max  current time, the introduction of the control and the final observation time (usually fixed)
//'        respectively espressed in days.
//' @param sum_beta,epsilon total risk from infected premisses and the primary infection rate respectively.
//' @param omega, introduses the effect of seasonality in the infection rate.
//'
//' @return It returns the next-event time in the non-homogeneous poisson process.
//'
//' @examples
//'simulate_NHPP_next_event(2,50,0.08,0.2,0.3,300)
//' @export
// [[Rcpp::export]]
double simulate_NHPP_next_event(const double t_now,const double t_intervention,const double sum_beta,const double epsilon,const double omega, double t_max){
  double u=1;
  double acp_pr=0;
  double t_sim,t_next;
  double total_rate_now = func_time_beta(t_now, t_intervention, sum_beta, epsilon, omega);

  if ((total_rate_now)>0) {
    // if ((sum_beta+alpha)>0) {

    while(u>acp_pr){

      t_sim = t_now + Rcpp::rexp(2,total_rate_now)[0];   // time of next event
      acp_pr = func_time_beta(t_sim, t_intervention, sum_beta, epsilon, omega)/(total_rate_now);  //acceptance probability

      u = R::runif(0,1);

      if(u<=acp_pr){
        t_next = t_sim;
      }

    }


    if (t_sim<=t_max){
      t_next = t_sim;
    } else {
      t_next=INFINITY;
    }


  }

  if ((total_rate_now)<=0) {
    t_next =INFINITY;
  }

  return(t_next);

}

/*======================================================================================================================================
                                                       Circcle line intersection in cpp-
 ======================================================================================================================================*/
vector<set_points_struct> circle_line_intersections_C (double circle_x, double circle_y, double r, int& n_line, grid_lines_struct& grid_lines){ // return a set of intersection points between the cirlce and the grid lines

  // note that formula assume the circle sits at (0,0), so need a shift of coordinates before and after//

  vector<set_points_struct>  set_points;
  set_points.reserve(n_line*2);

  double x_lb = -r;
  double x_ub = r;
  double y_lb = -r;
  double y_ub = r;


  for (int i=0;i<=(n_line-1);i++){

    if ( (grid_lines.orient_line[i]==1 & (grid_lines.coor_y_1[i]- circle_y)>=y_lb & (grid_lines.coor_y_1[i]- circle_y)<=y_ub) | (grid_lines.orient_line[i]==2 & (grid_lines.coor_x_1[i]-circle_x)>=x_lb & (grid_lines.coor_x_1[i]-circle_x)<=x_ub) ) { // only need to consider those grids line intersect with the smallest bounding box that contains the circle

      double dx = grid_lines.coor_x_2[i] - grid_lines.coor_x_1[i];
      double dy = grid_lines.coor_y_2[i] - grid_lines.coor_y_1[i];
      double dr = sqrt(pow(dx,2.0)+pow(dy,2.0));
      double D = (grid_lines.coor_x_1[i]-circle_x)*(grid_lines.coor_y_2[i]-circle_y) - (grid_lines.coor_x_2[i]-circle_x)*(grid_lines.coor_y_1[i]-circle_y);

      int sgn_dy=1;
      if (dy<0) sgn_dy = -1;

      double Delta = pow(r,2.0)*pow(dr,2.0) - pow(D,2.0);

      int l=(Delta>=0);
      switch(l){
      case 1:{ // tangent (one intersection) or two intersection



        if (Delta>0){

        double x_1 = (D*dy + sgn_dy*dx*sqrt(Delta))/pow(dr,2.0) + circle_x;
        double x_2 = (D*dy - sgn_dy*dx*sqrt(Delta))/pow(dr,2.0) + circle_x;

        double y_1 = (-D*dx + abs(dy)*sqrt(Delta))/pow(dr,2.0) + circle_y;
        double y_2 = (-D*dx - abs(dy)*sqrt(Delta))/pow(dr,2.0) + circle_y;

        double theta_1 = atan2((y_1-circle_y),(x_1-circle_x));
        double theta_2 = atan2((y_2-circle_y),(x_2-circle_x));  // Angle from the x-axis

        if (theta_1<0) theta_1 = 2*M_PI + theta_1;
        if (theta_2<0) theta_2 = 2*M_PI + theta_2;


        set_points_struct tmp_1, tmp_2;

        tmp_1.coor_x= x_1;
        tmp_1.coor_y=y_1;

        tmp_2.coor_x=x_2;
        tmp_2.coor_y=y_2;

        tmp_1.theta=theta_1;
        tmp_2.theta=theta_2;

        set_points.push_back(tmp_1);
        set_points.push_back(tmp_2);



      }


        if (Delta==0){
          double x = D*dy/pow(dr,2.0)+ circle_x;

          double y = -D*dx/pow(dr,2.0) + circle_y;

          double theta = atan2((y-circle_y),(x-circle_x));


          if (theta<0) theta = 2*M_PI + theta;

          set_points_struct tmp;

          tmp.coor_x=(x);
          tmp.coor_y=(y);

          tmp.theta=(theta);

          set_points.push_back(tmp);

        }

        break;
      }

      case 0:{ // no intersection
        // do nothing
        break;
      }

      }


    }



  }

  std::sort(set_points.begin(), set_points.end(), by_theta());

  return(set_points);
}

/*======================================================================================================================================
                                                        Intersection of a circle with cells-
 ======================================================================================================================================*/

//' Generates a set of intersection points between the cirlce and the grid lines.
//'
//'\code{circle_line_intersections} computes the intersections points of a given circle with the grid lines along with
//' the angle formed with the x-axis.
//'
//' @param circle_x,circle_y The the euclidean coordinates of the center of the circle.
//'
//' @param r The radius of the given circle.
//' @param n_line The number of grid lines.
//' @param grid_lines A 6 columns data frame with columns names as coor_x_1, coor_y_1, coor_x_2, coor_y_2, orient_line.
//' \describe{
//'     \item{coor_x_1, coor_y_1}{Coordinates of the left end point of the grid line }
//'     \item{coor_x_2, coor_y_2}{Coordinates of the right end point of the grid line }
//'     \item{orient_line}{Line orientation}
//'     \enumerate{
//'        \item indicates horizontal orientation
//'        \item indicates vetical orientation
//'     }
//'     \item{k_line}{Line numbering: bottom to top, then left to right}
//' }
//'
//' @return It returns a three columns data frame containing x-coordinate, y-coordanate of the intersection of the circle with the
//' grid, and the value of the angle betweem the x-axis and the line joining the center of the circle to the corresponding
//' intersection point.
//'
//' @examples
//' data(grid_line)
//' attach(grid_line)
//' circle_line_intersections(2022230,-3123109,10000,39,grid_line)
//' detach(grid_line)
//' @export
// [[Rcpp::export]]
DataFrame circle_line_intersections(double circle_x, double circle_y, double r, int n_line, DataFrame grid_lines){
  vector <set_points_struct> set_points;
  grid_lines_struct grid_lines1;

    grid_lines1.coor_x_1=grid_lines["coor_x_1"];
    grid_lines1.coor_y_1=grid_lines["coor_y_1"];
    grid_lines1.coor_x_2=grid_lines["coor_x_2"];
    grid_lines1.coor_y_2=grid_lines["coor_y_2"];
    grid_lines1.orient_line=grid_lines["orient_line"];

    set_points= circle_line_intersections_C(circle_x, circle_y, r, n_line, grid_lines1);

    int n_seqments = set_points.size();
    NumericVector x(n_seqments), y(n_seqments), theta(n_seqments);
    for(int i=0; i<n_seqments; i++){
      x[i]=set_points[i].coor_x;
      y[i]=set_points[i].coor_y;
      theta[i]=set_points[i].theta;
    }
    DataFrame set_points_R= DataFrame::create(Named("x")=x, Named("y")=y, Named("theta")=theta);

    return (set_points_R);
}

/*======================================================================================================================================
  Return a set of segments (with length, density, and absolute angle) correspond to a set_points  not exported version
======================================================================================================================================*/

vector<segments_struct> func_segments_attributes (DataFrame& set_points, NumericMatrix& pop_grid, double& r, para_aux& para_other){ // return a set of segments (with length, density, and absolute angle) correspond to a set_points

  vector<segments_struct> segments;

  int n_segments = set_points.nrows();

  segments.resize(n_segments);
  NumericVector theta = set_points["theta"];
  NumericVector x= set_points["x"];
  NumericVector y= set_points["y"];

  for (int i=0; i<=(n_segments-1); i++){
    if(i==0) {

     // segments[i].theta_abs = set_points[0].theta + (2*M_PI - set_points[set_points.size()-1].theta);
     segments[i].theta_abs = theta(0)+ (2*M_PI - theta(set_points.nrow()-1));


      // double x_1= set_points[set_points.size()-1].coor_x;
      // double y_1= set_points[set_points.size()-1].coor_y;
      // double x_2= set_points[0].coor_x;
      // double y_2= set_points[0].coor_y;
      double x_1= x(set_points.nrow()-1);
      double y_1= y(set_points.nrow()-1);
      double x_2= x(0);
      double y_2= y(0);
      double midpoint_x = (x_1+x_2)/2.0;
      double midpoint_y = (y_1+y_2)/2.0;
      int m = ceil((midpoint_y - para_other.y_min)/para_other.grid_size); // at mth row of the grid
      int n = ceil((midpoint_x - para_other.x_min)/para_other.grid_size); // at nth col..

      segments[i].m=m;
      segments[i].n=n;


      // switch(midpoint_x>para_other.x_max|midpoint_x<para_other.x_min|midpoint_y>para_other.y_max|midpoint_y<para_other.y_min){
      int l=(m>para_other.n_row_grid | n>para_other.n_col_grid | m<=0| n<=0);
      switch(l){
      case 1:{
        segments[i].den = 0.0;
        break;
      }
      case 0:{
        segments[i].den = pop_grid(m-1,n-1);
        break;
      }
      }


    }

    if(i!=0) {

// segments[i].theta_abs = set_points[i].theta - set_points[i-1].theta;
    segments[i].theta_abs = theta[i] - theta[i-1];

      // double x_1= set_points[i-1].coor_x;
      // double y_1= set_points[i-1].coor_y;
      // double x_2= set_points[i].coor_x;
      // double y_2= set_points[i].coor_y;
      double x_1= x(i-1);
      double y_1= y(i-1);
      double x_2= x(i);
      double y_2= y(i);
      double midpoint_x = (x_1+x_2)/2.0;
      double midpoint_y = (y_1+y_2)/2.0;
      int m = ceil((midpoint_y - para_other.y_min)/para_other.grid_size); // at mth row of the grid
      int n = ceil((midpoint_x - para_other.x_min)/para_other.grid_size); // at nth col..


      int l=(m>para_other.n_row_grid | n>para_other.n_col_grid | m<=0| n<=0);
      switch(l){
      case 1:{
        segments[i].den = 0.0;
        break;
      }
      case 0:{
        segments[i].den = pop_grid(m-1,n-1);
        break;
      }
      }
      segments[i].m=m;
      segments[i].n=n;

    }

    segments[i].len = r*(segments[i].theta_abs);



  }



  return(segments);

}

/*======================================================================================================================================
                         Return a set of segments (with length, density, and absolute angle) correspond to a set_points
 ======================================================================================================================================*/

//' Generates a set of segments (with length, density, and absolute angle) corresponding to a set of intersection points.
//'
//'\code{func_arcs_attributes} computes the intersectioins points of a given circle with the grid lines along with
//' the angle formed with the x-axis.
//'
//' @param set_points A data frame of intersection points between the cirlce and the grid lines.
//' @seealso{\code{\link{circle_line_intersections}}}
//' @param pop_grid  Population density of the grid a case resides. This is filled from bottom to top, then left to right.
//' @param r The travelling distance of the inoculum.
//' @param x_min,y_min  x/y min of the left 2 corners of the box.
//' @param n_row_grid,n_col_grid Number of rors and columns of the grid.
//' @param grid_size Grid resolution
//' //@inheritParams circle_line_intersections
//'
//' @return It returns a five columns data frame containing:
//'        \describe{
//'          \item{len_arc}{Length of the subtending arc delimited by the intersection points of the circle with center \code{circle_x} and \code{circle_y} with a grid}
//'          \item{dens}{Density of the grid where the new infection premisse resides}
//'          \item{theta}{Angle (between the source and the intersection points) specifying the direction of the inoculum}
//'        }
//'
//' @examples
//' data(bbtv)
//'   attach(bbtv)
//'   Dat<- bbtv[,c(2:6,8,10)]     # Coonsidering the essential part of the data
//'   Dat1<-subset(Dat,Dat$latitude> -27.3 & Dat$processedbananas%in%c("P&I","P", "NI") )  # data up in queensland (noth of brisbane)
//'   Dat1$treatmentdate[is.na(Dat1$treatmentdate)]<- Dat1$date[is.na(Dat1$treatmentdate)] # When NA, consider removal date as
//'   Dat1$detection<-as.numeric(difftime(as.Date(Dat1$date), as.Date("2011/01/01"), unit="days")) + runif(1)
//'   Dat1$removal<-as.numeric(difftime(as.Date(Dat1$treatmentdate), as.Date("2011/01/01"), unit="days"))
//'   Datt<-Dat1[,c(3,4,7,8,2)]
//' Datt$detection=Datt$detection + runif(nrow(Datt))
//'   Datt=Datt[with(Datt,order(Datt$detection)),]
//'
//' # Australian reference system
//' sp::coordinates(Datt) <- c("longitude", "latitude")
//'   sp::proj4string(Datt) <- sp::CRS("+init=epsg:4326")
//'   australianCRS <- sp::CRS("+init=epsg:3577")
//'
//'   pointsinaustraliangrid = sp::spTransform(Datt,australianCRS)
//'
//' #Raster
//'   rast <- raster::raster()
//'     raster::extent(rast) <- raster::extent(pointsinaustraliangrid) # Set same extent
//'
//'     raster::res(rast)=5000 #Set resolution
//'
//' # And then ... rasterize it! This creates a grid version
//' # of your points using the cells of rast,
//'
//'     rast2 <- raster::rasterize(pointsinaustraliangrid, rast, 1, fun=sum)
//'
//' # Extract infos on the grid
//'
//'
//'       n_row_grid=nrow_grid=raster::nrow(rast)
//'         n_col_grid=ncol_grid=raster::ncol(rast)
//'         grid_size=raster::res(rast)[1]     # Resolution
//'
//'         n_line=(nrow_grid+1) + (ncol_grid +1)  # Number of grid  lines
//'
//'         x_min=raster::xmin(rast)  # min max of the bounding box
//'         x_max=raster::xmax(rast)
//'
//'         y_min=raster::ymin(rast)
//'         y_max=raster::ymax(rast)
//'
//'         da=as.data.frame(pointsinaustraliangrid)
//'
//'         pop_per_grid=raster::values(rast2)
//'         pop_per_grid[is.na(pop_per_grid)]=0
//'       mat=matrix(pop_per_grid,nrow = nrow_grid, byrow = TRUE )
//'         pop_grid=apply(mat,2,rev)     # population per grid
//'
//' # Structure of the grid
//'         x=seq(x_min,x_max,grid_size)
//'           y=seq(y_min,y_max,grid_size)
//'
//'           grid_lines=array(0,c(n_line,6))
//'           for(i in 1:n_line){
//'             if(i<=(nrow_grid +1)){
//'               grid_lines[i,]=c(i,1,x[1],y[i],x[length(x)],y[i])
//'             }
//'             else{
//'               grid_lines[i,]=c(i,2,x[i-length(y)],y[1],x[i-length(y)],y[length(y)])
//'             }
//'           }
//'
//'           grid_lines=as.data.frame(grid_lines)
//'             colnames(grid_lines)<- c("indx","orient_line","coor_x_1","coor_y_1","coor_x_2","coor_y_2")
//'             circle_x=2022230
//'           circle_y=-3123109
//'           r=10000
//'
//'           set_points=circle_line_intersections(2022230,-3123109,10000,39,grid_lines)
//'             func_arcs_attributes(set_points, pop_grid, r, x_min, y_min, grid_size, n_row_grid, n_col_grid);
//'         detach(bbtv)
//'
//' @export
// [[Rcpp::export]]
DataFrame func_arcs_attributes(DataFrame set_points, NumericMatrix& pop_grid, double r, double x_min, double y_min, double grid_size, double n_row_grid, double n_col_grid){
  //DataFrame set_points;
  grid_lines_struct grid_lines1;
  para_aux para_other;

  vector<segments_struct> segments;

  // // Read in grid lines structure
  // grid_lines1.coor_x_1=grid_lines["coor_x_1"];
  // grid_lines1.coor_y_1=grid_lines["coor_y_1"];
  // grid_lines1.coor_x_2=grid_lines["coor_x_2"];
  // grid_lines1.coor_y_2=grid_lines["coor_y_2"];
  // grid_lines1.orient_line=grid_lines["orient_line"];

  // Read in grid parameters: size, density, boxes coordinate etc..
  para_other.dimen_x = n_row_grid;
  para_other.dimen_y = n_col_grid;
  para_other.grid_size = grid_size;
 // para_other.n_line = n_line;
 // para_other.total_pop_of_grid = pop_grid;
  para_other.x_min = x_min;
  para_other.y_min =y_min;

  //set_points= circle_line_intersections(circle_x, circle_y, r, n_line, grid_lines);  // Intersection points with angles
  //Rcout<<r<<"\n";
  segments = func_segments_attributes (set_points, pop_grid, r, para_other);  // angle and length of arcs
  //Rcout<<r<<"\n";


 //  Rcout<<segments.size()<<"\n";
  int n_seqments = segments.size();
  NumericVector len(n_seqments), theta(n_seqments), dens(n_seqments), n(n_seqments), m(n_seqments);
  for(int i=0; i<n_seqments; i++){
    len[i]=segments[i].len;
    theta[i]=segments[i].theta_abs;
    dens[i]=segments[i].den;
    n[i] = segments[i].n;
    m[i] = segments[i].m;
  }
  DataFrame segments_R= DataFrame::create(Named("len_arc")=len, Named("dens")=dens,Named("theta_abs")=theta, Named("n_th_col")=n, Named("m_th_row")=m);

  return (segments_R);
}

/*======================================================================================================================================
                                        Cyclic latent period using the leaves emergence rate (LER)
 ======================================================================================================================================*/

//' @export
// [[Rcpp::export]]
double f(double x0, double E,double A){
  double x, a, b, c, d, rat;
  double c1=2*M_PI/365;
  c=E;
  a=A;
  b=A/c1;
  d=-b*sin(c1*(c-15));
  rat=A*cos(c1*(c+x0-15)) + 0.062;
  return (rat*exp(-(d+b*sin(c1*(x0+c-15))+0.062*x0)));
}



/*======================================================================================================================================
                                Sampling from the cyclic distribution using M-H algorithm
 ======================================================================================================================================*/
//' @export
// [[Rcpp::export]]
NumericVector BTFinv1 (double E, double A, double t0)
{
  double q, x, pacc;
  NumericVector par(1000);
  par[0]=t0;
  for( int i=1;i<1000;i++){
    q=3*Rcpp::rnorm(2,0,1)[0];
    x=par[i-1]*exp(q);
    pacc=f(x,E,A)/f(par[i-1],E,A)*exp(q);
    if(pacc>Rcpp::runif(2,0,1)[0]){
      par[i]=x;
    }
    else{
      par[i]=par[i-1];
    }
  }



  return (par);
}

/*======================================================================================================================================
                            Latent period depending on the distribution used
 ======================================================================================================================================*/
//' Sample the latent period corresponding to each specified model.
//'
//'\code{E_to_I} Geratee a random draw from the distribution specified for the latent period.
//'
//' @param EI_model A given integer characterising the type pf distribution for the latent period.
//' \enumerate{
//'          \item Cyclic distribution using the LER (leaf emergence rate) derived by R Allen \url{http://www.publish.csiro.au/ar/pdf/ar9780535}:
//'          \deqn{\frac{dL}{dt}=0.056cos(t-15) + 0.062}
//'          \item Gamma distribution with given rate and shape parameters: \code{\link{rgamma}}
//'          \item Expoential distribution with given rate: \code{\link{rexp}}, the default.
//'
//'        }
//' @param E The time at which the contact occurs, or strictly speaking the exposure time.
//' @param mu_lat The mean latent period.
//' @param var_lat The variance of the latent period for the gamma distribution.
//'
//' @return It returns a random draw from the latent period given the time of exposure.
//'
//' @examples
//' E_to_I(1,-20,30,10)
//' E_to_I(2,-20,30,10)
//' E_to_I(3,-20,30,10)
//' @export
 // [[Rcpp::export]]
double E_to_I (int EI_model, double E, double mu_lat, double var_lat)
{
  double dt, a, b, ru;
  int k;
  switch (EI_model) {

  case 1:
    ru=runif(1,0,1)[0];
    k=ru*500;
    mu_lat=0.056;
    //dt=mean(BTFinv1(E,mu_lat,10.0));
    dt=BTFinv1(E,mu_lat,10.0)[k+450];
    break;

  case 2:
    a=mu_lat*mu_lat/var_lat;
    b=var_lat/mu_lat;
    dt=Rcpp::rgamma(1,a,1/b)[0];
    break;

  default:
    dt=Rcpp::rexp(1,1/mu_lat)[0];

  }
  return dt;
}

/*======================================================================================================================================
                                               Assign the proper beta for an individual by age
 ======================================================================================================================================*/
// [[Rcpp::export]]
double beta_by_age(int age, NumericVector beta_by_age_vector){

 return(beta_by_age_vector[age]);

}

/*======================================================================================================================================
                                      Sampling from the appropriate kernel
======================================================================================================================================*/
//' Sample the distance at which the inoculum will travel.
//'
//'\code{Samp_dis} Generate a random draw from the distribution specified for the latent period.
//'
//' @param kern_model A given integer characterising the type of distribution for the kernel for both short and long range interaction.
//' \enumerate{
//'          \item exponential-exponential
//'          \item cauchy-cauchy
//'          \item linear combination of both (see paper by Montemayer et. al), the default.
//'
//'        }
//' @param alpha1 Dispersal scale parameter for the local spread kernel.
//' @param alpha2 Dispersal scale parameter for the long range interaction.
//'
//' @return It returns a random distance the inoculum will travel to.
//'
//' @examples
//' Samp_dis (1, 0.2, 0.3)
//' @export
// [[Rcpp::export]]
double Samp_dis (int kern_model, double alpha1, double alpha2)
{
  double r;
  double gama = .5;
  switch (kern_model) {

  case 1:   // expo-expo
    if(runif(1,0,1)[0]< gama){
      r=Rcpp::rexp(1,1/alpha1)[0];
    }
    else{
      r=Rcpp::rexp(1,1/alpha2)[0];
    }

    break;

  case 2:  // Cauchy-cauchy
    if(runif(1,0,1)[0]< gama){
      r=Rcpp::rcauchy(1,0,alpha1)[0];
    }
    else{
      r=Rcpp::rcauchy(1,0,alpha2)[0];
    }

    break;

  default:
    r=Rcpp::rexp(1,1/alpha1)[0];
    // to be changed accordling

  }
  return r;
}

