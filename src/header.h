

 #include <RcppGSL.h>
#include <Rcpp.h>
using namespace Rcpp;

#ifndef _FUNCTIONS_H_INCLUDED_ //include guard
#define _FUNCTIONS_H_INCLUDED_

#include<iostream>
#include<vector>
#include<algorithm>
#include <fstream>
#include <string>
#include<sstream>
#include <ctime>
#include <cstdlib>
#include<time.h>
#include <cmath>
#include <stdio.h>
#include <numeric>
#include <set>
#include <unistd.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <boost/date_time/gregorian/gregorian.hpp>

#define path1 "/Users/sylauadmin/Dropbox/PU/scripts/ebola/contact_model_extension/simulation_inference/cplus_script/" /*the script directory*/
#define path3 "/Users/sylauadmin/PU/results/ebola/contact_model_extension/output/figures/"      /*as above*/

#define path2 "/Users/sylauadmin/PU/results/ebola/contact_model_extension/simulation/permutated_pop_grid/sim_10/"      /*the simulated data*/
#define path4 "/Users/sylauadmin/PU/results/ebola/contact_model_extension/simulation_inference/permutated_pop_grid/sim_10/mcmc samples 1/"   /*the mcmc samples */

//#define path4 "/home/sl451/HWU/spatial_seir_inference/output/base case/mcmc samples 3/"      /*as above*/
//#define path4 "/home/sl451/HWU/spatial_seir_inference/output/scenarios/scenario C/mcmc samples 1/"      /*as above*/

// #define path1 "/home/max/profile/HWU/spatial_seir_inference/script/cplus-script/" /*the working directory*/
// #define path2 "/home/max/profile/HWU/spatial_seir_simulation/output/simulated data/"      /*as above*/
// #define path3 "/home/max/profile/HWU/spatial_seir_inference/output/figures/"      /*as above*/
// #define path4 "/home/max/profile/HWU/spatial_seir_inference/output/mcmc samples/"      /*as above*/

// #define path1 "/home/max/HWU/spatial_seir_inference/script/cplus-script/" /*the working directory*/
// #define path2 "/home/max/HWU/spatial_seir_simulation/output/simulated data/"      /*as above*/
// #define path3 "/home/max/HWU/spatial_seir_inference/output/figures/"      /*as above*/
// #define path4 "/home/max/HWU/spatial_seir_inference/output/mcmc samples/"      /*as above*/

//#define EPSILON 0.0000001   // Define the tolerance when comparing doubles
//#define DOUBLE_EQ(x,v) (((v - EPSILON) < x) && (x <( v + EPSILON)))

using namespace std;
using namespace boost::gregorian;

struct para_key {
  double alpha,a,b,mu_lat,var_lat,c,d,k_1,k_2,beta,t_0;


  //double beta_1, beta_2, beta_3, beta_4;
  vector<double> beta_age;

  double omega;
  //
  double eta;
  //
  double stb_1, stb_2;
};

//----------------

struct para_aux{
  int n;
  // string kernel_type;
  double dimen_x,dimen_y;
  double t_max,unassigned_time;
  // int seed;
  int model;
  int n_row_grid, n_col_grid; // number of rows and colmumns of the grids
  double grid_size; // width of the square grid

  int n_line; // number of grid lines


  double x_min, y_min; //x/y min/max of the 4 corners of the bounding box
  double x_max,y_max;

  double total_pop_of_grid; // sum of pop density of all grids being considered

  double T; // the time where intervention starts

};

//-----------------

struct epi_struct {
  vector<int> k;
  vector<double> t_e, t_i,t_r;
  // vector<double> q, t_e, t_i,t_r;
  // vector<int> status;

  // vector<double> coor_x,coor_y;

  // vector<int> gp_stb;
  vector<double> stb;
  vector<int> onset_type; // 1=trustworthy onset, 0-onset requires imputations
  vector<int> lab_result; // 1=positive, 0=negative
  vector<int> age_gp; // age group =1,2,3,4
  vector<int> sex;


  vector<int> infected_source;

  // vector<int> grid_num;

  vector<int> grid_num_row; // indicates which row of the grids a case is in
  vector<int> grid_num_col; // indicates which column of the grids a case is in

  vector<double> pop; // the pop density of the grid a case resides

};


//----

// struct grid_struct{ // pop density of the grids and coordinates of center, and 4 corners of the grids in the study region (numbering: start from right-bottom and then anti-clockwise)

// vector<int> k_grid; //index for grids

// vector<double> coor_x,coor_y; // centers

// vector<double> coor_x_1,coor_y_1;
// vector<double> coor_x_2,coor_y_2;
// vector<double> coor_x_3,coor_y_3;
// vector<double> coor_x_4,coor_y_4;

// vector<double> pop;


// };
//-----

// struct grid_bound{ // coordinates of the 4 corners of the bounding box

// double x_min, y_min;
// double x_min,y_max;
// double x_max, y_max;
// double x_max,y_min;

// };
// //--

struct grid_lines_struct{ //

  vector<int> k_line; //index for grids, numbering: bottom to left, then left to right

  vector<int> orient_line; // orientation of line, 1=horizontial, 2= vertical

  vector<double> coor_x_1,coor_y_1; // first end point
  vector<double> coor_x_2,coor_y_2; // second end point

};
//----

struct set_points_struct{

  // vector<double> coor_x;
  // vector<double> coor_y;
  // vector<double> theta; // the theta measured from the first quadrant centered at the circle center (need change of coordinate)

  double coor_x;
  double coor_y;
  double theta; // note: the theta measured from the x-axis at first quadrant centered at the circle center (need change of coordinate)


};

struct by_theta {
  bool operator()(set_points_struct const &a, set_points_struct const &b) {
    return a.theta < b.theta;
  }
};

//---

struct segments_struct{

  double theta_abs ; // the absoulte angle between two consecutibe points (not the theta that measured from the x-axis used in set_point_struct)
  double len; // length of the arc segment between two consecutive intersection points
  double den; // density of grid that the arc segment ... resides
  int m, n;
};



//---

// for sorting according to t_e in imputing source //

struct paired_to_t_e {
  int source;
  double distance;
  double t_e;
};

struct by_t_e {
  bool operator()(paired_to_t_e const &a, paired_to_t_e const &b) {
    return a.t_e < b.t_e;
  }
};
//------------------

struct lh_SQUARE { // the contribution to likelihood function from each individual is subdivided into f_U, f_E, f_I,f_R
  vector<long double> f_U, q_T, kt_sum_U;
  vector<long double> f_E, g_E,k_sum_E, h_E,q_E,kt_sum_E;
  vector<long double> f_I, f_EnI;
  vector<long double> f_R, f_InR;

  vector<long double> f_Grid, f_Arc; // f_Grid=the likelihood for landing at a grid | distance; f_Arc = likelihood of landing at a point of an arc within the region (uniformly)

  double f_Surv; // survival prob for the whole time period

};

//------------------------------------------
struct quadratic_params
{
  double t2, l, b1, a1;
};



//----

double func_time_beta (const double& ,const double& ,const double& , const double& );

double func_time_alpha (const double& ,const double& ,const double& , const double& );

double quadratic (double x, void *params);
double quadratic_deriv (double x, void *params);
void quadratic_fdf (double x, void *params, double *y, double *dy);

double quadratic1 (double x, void *params);
double quadratic1_deriv (double x, void *params);
void quadratic1_fdf (double x, void *params, double *y, double *dy);

//----

vector<set_points_struct> circle_line_intersections (double , double , double , int& , grid_lines_struct& ); // return a set of intersection points between the cirlce and the grid lines

//--

vector<segments_struct> func_segments_attributes (vector<set_points_struct>& , vector < vector<double> >&, double&,para_aux&);


//---------------------
double dtnorm(double, double, double, double);

//----------------------------

long double func_distance_for_main (double , double , double , double );

inline long double func_kernel (double, double, double, double, double, double, const string&); // function prototype for calculating kernel distance

double log_lh_func (lh_SQUARE, int);

//--------------------------


//-----
class FUNC {

private:

  double alpha_Clh;
  // double beta_Clh;
  double a_Clh;
  double b_Clh;
  double mu_lat_Clh;
  double var_lat_Clh;
  double c_Clh;
  double d_Clh;
  double k_1_Clh;
  double k_2_Clh;
  double t_0_Clh;
  double omega_Clh;

  // double beta_1_Clh, beta_2_Clh, beta_3_Clh, beta_4_Clh;
  vector<double> beta_age_Clh;

  double eta_Clh;

  int n_Clh;
  // string kernel_type_Clh;
  double dimen_x_Clh,dimen_y_Clh;
  double t_max_Clh,unassigned_time_Clh;
  int    model_Clh;
  // int seed_Clh;

  int n_row_grid_Clh, n_col_grid_Clh, n_line_Clh;
  double grid_size_Clh;


  vector< vector<double> > coordinate_Clh;

  vector< vector<double> > pop_grid_Clh;

  // grid_struct grid_data_Clh;

  grid_lines_struct grid_lines_Clh;

  vector<int> xi_U_Clh, xi_E_Clh, xi_E_minus_Clh, xi_I_Clh, xi_R_Clh, xi_EnI_Clh, xi_InR_Clh;
  vector<double> t_e_Clh, t_i_Clh, t_r_Clh;

  vector<double> pop_Clh; // the pop density vector of individuals

  vector<int> age_gp_Clh;

  vector<int> index_Clh;

  vector<double> stb_Clh;

  vector <int> infected_source_Clh;

  para_aux para_other_Clh;

public:

  void set_para (para_key& para_current_arg, para_aux& para_other_arg, vector< vector<double> >& coordinate_arg,
                 vector<int> xi_U_arg, vector<int> xi_E_arg, vector<int> xi_E_minus_arg ,vector<int> xi_I_arg,
                 vector<int> xi_R_arg, vector<int> xi_EnI_arg, vector<int> xi_InR_arg, vector<double> t_e_arg,
                 vector<double> t_i_arg, vector<double> t_r_arg, vector<double>& pop_arg, vector<double>& stb_arg,
                 vector<int> index_arg, vector<int>& age_gp_arg, vector<int> infected_source_arg,
                 vector< vector<double> >& pop_grid_arg, grid_lines_struct& grid_lines_arg) {

    alpha_Clh=para_current_arg.alpha;
    // beta_Clh=para_current_arg.beta;
    a_Clh=para_current_arg.a;
    b_Clh=para_current_arg.b;
    mu_lat_Clh=para_current_arg.mu_lat;
    var_lat_Clh=para_current_arg.var_lat;
    c_Clh=para_current_arg.c;
    d_Clh=para_current_arg.d;
    k_1_Clh=para_current_arg.k_1;
    k_2_Clh=para_current_arg.k_2;
    t_0_Clh=para_current_arg.t_0;


    omega_Clh=para_current_arg.omega;


    // beta_1_Clh=para_current_arg.beta_1;
    // beta_2_Clh=para_current_arg.beta_2;
    // beta_3_Clh=para_current_arg.beta_3;
    // beta_4_Clh=para_current_arg.beta_4;

    beta_age_Clh = para_current_arg.beta_age;

    eta_Clh = para_current_arg.eta;

    n_Clh = para_other_arg.n;
    // kernel_type_Clh = para_other_arg.kernel_type;
    dimen_x_Clh = para_other_arg.dimen_x;
    dimen_y_Clh = para_other_arg.dimen_y;
    t_max_Clh = para_other_arg.t_max;
    model_Clh=para_other_arg.model;
    unassigned_time_Clh = para_other_arg.unassigned_time;
    // seed_Clh = para_other_arg.seed;

    n_row_grid_Clh = para_other_arg.n_row_grid;
    n_col_grid_Clh = para_other_arg.n_col_grid;

    grid_size_Clh = para_other_arg.grid_size;
    n_line_Clh = para_other_arg.n_line;

    para_other_Clh = para_other_arg;

    coordinate_Clh = coordinate_arg;

    pop_grid_Clh = pop_grid_arg;

    // grid_data_Clh =  grid_data_arg;
    grid_lines_Clh =  grid_lines_arg;

    xi_U_Clh = xi_U_arg;
    xi_E_Clh = xi_E_arg;
    xi_E_minus_Clh = xi_E_minus_arg;
    xi_I_Clh = xi_I_arg;
    xi_R_Clh = xi_R_arg;
    xi_EnI_Clh = xi_EnI_arg;
    xi_InR_Clh = xi_InR_arg;


    t_e_Clh = t_e_arg;
    t_i_Clh = t_i_arg;
    t_r_Clh = t_r_arg;

    pop_Clh = pop_arg;

    index_Clh = index_arg;

    stb_Clh = stb_arg;

    age_gp_Clh= age_gp_arg;

    infected_source_Clh = infected_source_arg;

  }


  void initialize_kernel_mat (vector< vector<double> >&, vector<double>&); // function prototype for initializing kernel distance

  void initialize_delta_mat (vector< vector<double> >&); // function prototype for initializing length of exposure time

  void initialize_lh_square (lh_SQUARE&, vector< vector<double> >, vector< vector<double> >, vector<double>&, const int);

};

//-----------------------------

class mcmc_UPDATE {

private:

  // double alpha_CUPDATE;
  // double beta_CUPDATE;
  // double a_CUPDATE;
  // double b_CUPDATE;
  // double c_CUPDATE;
  // double d_CUPDATE;
  // double k_1_CUPDATE;
  // double k_2_CUPDATE;

  int n_CUPDATE;

  int n_row_grid_CUPDATE, n_col_grid_CUPDATE;
  int n_line_CUPDATE;
  double grid_size_CUPDATE;

  // string kernel_type_CUPDATE;
  double dimen_x_CUPDATE,dimen_y_CUPDATE;
  double t_max_CUPDATE,unassigned_time_CUPDATE;
  int    model_CUPDATE;
  // int seed_CUPDATE;

  para_aux para_other_CUPDATE;

  vector< vector<double> > coordinate_CUPDATE;

  vector< vector<double> > pop_grid_CUPDATE;

  // grid_struct grid_data_CUPDATE;
  grid_lines_struct grid_lines_CUPDATE;

  vector < vector<double> > distance_mat_CUPDATE;

  // vector<int> xi_U_CUPDATE, xi_E_CUPDATE, xi_E_minus_CUPDATE, xi_I_CUPDATE, xi_R_CUPDATE;
  // vector<double> t_e_CUPDATE, t_i_CUPDATE, t_r_CUPDATE;

  vector<int> index_CUPDATE;

  vector<int> age_gp_CUPDATE;

  vector<double> pop_CUPDATE; // pop density vector of all individuals

public:

  void set_para (para_aux& para_other_arg, vector< vector<double> >& coordinate_arg, vector<int>& age_gp_arg, vector<double>& pop_arg, vector< vector<double> >&  pop_grid_arg, grid_lines_struct& grid_lines_arg) {

    // alpha_CUPDATE=para_current_arg.alpha;
    // beta_CUPDATE=para_current_arg.beta;
    // a_CUPDATE=para_current_arg.a;
    // b_CUPDATE=para_current_arg.b;
    // c_CUPDATE=para_current_arg.c;
    // d_CUPDATE=para_current_arg.d;
    // k_1_CUPDATE=para_current_arg.k_1;
    // k_2_CUPDATE=para_current_arg.k_2;

    n_CUPDATE = para_other_arg.n;
    // kernel_type_CUPDATE = para_other_arg.kernel_type;
    dimen_x_CUPDATE = para_other_arg.dimen_x;
    dimen_y_CUPDATE = para_other_arg.dimen_y;
    t_max_CUPDATE = para_other_arg.t_max;
    model_CUPDATE=para_other_arg.model;
    unassigned_time_CUPDATE = para_other_arg.unassigned_time;
    // seed_CUPDATE = para_other_arg.seed;

    n_row_grid_CUPDATE = para_other_arg.n_row_grid;
    n_col_grid_CUPDATE = para_other_arg.n_col_grid;

    grid_size_CUPDATE = para_other_arg.grid_size;
    n_line_CUPDATE = para_other_arg.n_line;

    para_other_CUPDATE = para_other_arg;

    coordinate_CUPDATE = coordinate_arg;
    pop_grid_CUPDATE = pop_grid_arg;

    age_gp_CUPDATE = age_gp_arg;

    pop_CUPDATE = pop_arg;

    // grid_data_CUPDATE = grid_data_arg;
    grid_lines_CUPDATE = grid_lines_arg;

    // distance_mat_CUPDATE = distance_mat_arg;

    /*
    xi_U_CUPDATE = xi_U_arg;
    xi_E_CUPDATE = xi_E_arg;
    xi_E_minus_CUPDATE = xi_E_minus_arg;
    xi_I_CUPDATE = xi_I_arg;
    xi_R_CUPDATE = xi_R_arg;

    t_e_CUPDATE = t_e_arg;
    t_i_CUPDATE = t_i_arg;
    t_r_CUPDATE = t_r_arg;

    index_CUPDATE = index_arg;*/
  }

public:

  // void eta_update(lh_SQUARE&, double&, const vector<int>&, const vector<int>&, const vector <double>&,const vector <double>&, const vector<int>&, const vector<double>& , para_key&, const vector<int>&, int);

  void alpha_update(lh_SQUARE& , double& , vector< vector<double> >&, const  vector< vector<double> >&, const vector<int>&, const vector<int>&,const vector<int>&, const vector<double>&, const vector<double>&, const vector<double>&, const  vector<int>&, para_key&,   const vector<double>& , vector<double>&, const vector<int>&,int, int );

  void omega_update(lh_SQUARE& , double& , vector< vector<double> >&, const  vector< vector<double> >&, const vector<int>&, const vector<int>&,const vector<int>&, const vector<double>&, const vector<double>&, const vector<double>&, const  vector<int>&, para_key&,   const vector<double>& , vector<double>&, const vector<int>&,int, int );


  // void beta_update(lh_SQUARE&, double&, const vector<int>&, const vector<int>&, const vector <double>&,const vector <double>&, const vector<int>&, const vector<double>& , para_key&, int);

  void beta_1_update(lh_SQUARE& , double& , vector< vector<double> >&, const  vector< vector<double> >&, const vector<int>&, const vector<int>&,const vector<int>&, const vector<double>&, const vector<double>&, const vector<double>&, const  vector<int>&, para_key&,   const vector<double>& , vector<double>&, const vector<int>&,int, int );
  void beta_2_update(lh_SQUARE& , double& , vector< vector<double> >&, const  vector< vector<double> >&, const vector<int>&, const vector<int>&,const vector<int>&, const vector<double>&, const vector<double>&, const vector<double>&, const  vector<int>&, para_key&,   const vector<double>& , vector<double>&,const vector<int>&, int , int);
  void beta_3_update(lh_SQUARE& , double& , vector< vector<double> >&, const  vector< vector<double> >&, const vector<int>&, const vector<int>&,const vector<int>&, const vector<double>&, const vector<double>&, const vector<double>&, const  vector<int>&, para_key&,   const vector<double>& , vector<double>&,const vector<int>&, int, int );
  void beta_4_update(lh_SQUARE& , double& , vector< vector<double> >&, const  vector< vector<double> >&, const vector<int>&, const vector<int>&,const vector<int>&, const vector<double>&, const vector<double>&, const vector<double>&, const  vector<int>&, para_key&,   const vector<double>& , vector<double>&, const vector<int>&,int, int );
  void beta_5_update(lh_SQUARE& , double& , vector< vector<double> >&, const  vector< vector<double> >&, const vector<int>&, const vector<int>&,const vector<int>&, const vector<double>&, const vector<double>&, const vector<double>&, const  vector<int>&, para_key&,   const vector<double>& , vector<double>&, const vector<int>&,int, int );

  void c_update(lh_SQUARE&, double&, const vector<int>&,const vector<int>&, const vector<int>&, const vector<double>&, const vector<double>&,const vector<int>&, para_key& ,const vector<int>&,int, int);
  void d_update(lh_SQUARE&, double&, const vector<int>&,const vector<int>&, const vector<int>&, const vector<double>&, const vector<double>&, const vector<int>&, para_key& ,const vector<int>&,int, int);
  void k_1_update(lh_SQUARE& , double& , vector< vector<double> >&, const  vector< vector<double> >&, const vector<int>&, const vector<int>&,const vector<int>&, const vector<double>&, const vector<double>&, const vector<double>&, const  vector<int>&, para_key&,   const vector<double>& , vector<double>&, const vector<int>&,int, int );
  void k_2_update(lh_SQUARE& , double& , vector< vector<double> >&, const  vector< vector<double> >&, const vector<int>&, const vector<int>&,const vector<int>&, const vector<double>&, const vector<double>&, const vector<double>&, const  vector<int>&, para_key&,   const vector<double>& ,vector<double>&, const vector<int>&,int, int );


  void mu_lat_update(lh_SQUARE&, double&, const vector<int>&,const vector<int>&,  const vector<int>&,  const vector<double>&, const vector<double>&,const vector<double>&, const vector<int>&, para_key&, const vector<int>&,int, int);
  void var_lat_update(lh_SQUARE&, double&, const vector<int>&,const vector<int>&, const vector<int>&, const vector<double>&, const vector<double>&,const vector<double>&, const vector<int>&, para_key&, const vector<int>&,int, int);

  void source_t_e_update(lh_SQUARE&, double&, const vector< vector<double> >&, vector< vector<double> >&, vector<int>&, vector<int>&, vector<int>& , const vector<int>&, vector<int>&, const vector<double>&, const vector<double>&, vector<double>& , vector<int>&, const para_key&,  const vector<double>& , const vector<double>&,  vector<int>&, vector<double>&, int, int, int);

  void t_e_replct(lh_SQUARE&, double&, const vector< vector<double> >&, vector< vector<double> >&, vector<int>&, vector<int>&, vector<int>& , const vector<int>&, vector<int>&, const vector<double>&, const vector<double>&, vector<double>& , vector<int>&, const para_key&,  const vector<double>& , const vector<double>&, const vector<int>&,int, int);

  void t_e_add_del(lh_SQUARE&, double&, const vector< vector<double> >&, vector< vector<double> >&, vector<int>&, vector<int>&, vector<int>& , const vector<int>&, vector<int>&, const vector<double>&, const vector<double>&, vector<double>& , vector<int>&, const para_key&,  const vector<double>& ,const vector<double>&,const vector<int>&,int, int);
  void t_i_update(lh_SQUARE& , double& , const vector< vector<double> >&, vector< vector<double> >& , const vector<int>& , const vector<int>& , const vector<int>& , const vector<int>& , const vector<int>& , const vector<int>& , const vector<int>& ,const vector<double>& , vector<double>& , const vector<double>& , const vector<int>& , const para_key& , const vector<double>& , const vector<double>& , const vector<int>&,int , int);

  void t_i_update_2(lh_SQUARE& , double& , const vector< vector<double> >&, vector< vector<double> >& , const vector<int>& , const vector<int>& , const vector<int>& , const vector<int>& , const vector<int>& , const vector<int>& , const vector<int>& ,const vector<double>& , vector<double>& , const vector<double>& , const vector<int>& , const para_key& , const vector<double>& , const vector<double>& , const vector<int>&,int, int );


  void t_i_update_with_obs(lh_SQUARE& , double& , const vector< vector<double> >&, vector< vector<double> >& , const vector<int>& , const vector<int>& , const vector<int>& , const vector<int>& , const vector<int>& , const vector<int>& , const vector<int>& ,const vector<double>& , vector<double>& , vector<double>& , const vector<double>& , const vector<int>& , const para_key& , const vector<double>& , const vector<double>& ,const vector<int>&, int, int );
  void t_0_update(lh_SQUARE&, double&, const vector< vector<double> >&, vector< vector<double> >&, vector<int>&, vector<int>&, vector<int>& , const vector<int>&, vector<int>&, const vector<double>&, const vector<double>&, vector<double>& , vector<int>&, para_key&,  const vector<double>& , const vector<double>&,  vector<int>&, vector<double>&, int, int);

  // void sample_source_distance (const vector< vector<double> >&, const vector<int>&, const vector<int>&,  const vector<double>&, const vector<double>&, const vector<double>&, const para_key&, const para_aux&, const vector<double>& , const vector<double>&, const vector<int>&,int); // sample sources and dispersal distance


  void residual_gnl(const vector< vector<double> >&, const vector<int>&, const vector<int>&, const vector<int>&,  const vector<double>&, const vector<double>&, const vector<double>&, const para_key&, const para_aux&,  const vector<double>& , const vector<double>&, const vector<int>&,int, int);

  void residual_kernel(const vector< vector<double> >&, const vector<int>&, const vector<int>&,  const vector<double>&, const vector<double>&, const vector<double>&, const para_key&, const para_aux&, const vector<double>& , const vector<double>&, const vector<int>&,int, int); // output the defined residuals

  void residual_lat( const vector<int>&, const vector<int>&,  const vector<int>& , const vector<double>&, const vector<double>&, const para_key&, const para_aux&,  const vector<int>&,int, int);



};


//-----------------------------


//-----------------------------

void IO_simpara(para_key&, para_aux&);
void IO_simdata(para_key, para_aux ,vector < vector<double> >& , epi_struct& ,vector<int>&, vector< vector<double> >&, grid_lines_struct&);
//void lh_plot (para_key, para_aux , vector< vector<double> >, vector<int>, vector<int> , vector<int> , vector<int> ,vector<int>, vector<int> ,vector<int>, vector<double> , vector<double> , vector<double>, vector<int>);


#endif
