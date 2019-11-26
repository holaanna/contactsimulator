#include "header.h"

// [[Rcpp::interfaces(r, cpp)]]
/*======================================================================================================================================
                              Intensity of the distributin of the number of Leaf that emergences between t1 and t2
 ======================================================================================================================================*/
//' Functional form of the intensity of the non-homogenous Poisson process.
//'
//'\code{Sum_Leaf_emergence_rate} returns the value of the intensity as a function of t1 and t2.
//'
//' @param t1,t2  Time of the start and end of observation respectively espressed in days.
//' @param a,b    Parameters of the model.
//'
//' @return It returns the intensity in the non-homogeneous poisson process.
//' @references
//' @examples
//' Leaf_emrgence_rate(20,50,0.062,0.056)
//' @export
// [[Rcpp::export]]
NumericVector Sum_Leaf_emergence_rate(NumericVector t1, NumericVector t2, double a, double b){
  double omega = 2*M_PI/365.0;
  NumericVector h_t = a*((t2-t1) + b/omega*(sin(omega*(t2-15)) - sin(omega*(t1-15))));
  return(h_t);
}

/*======================================================================================================================================
                              Log-likelihood
 ======================================================================================================================================*/
//' The log-likelihood of the non-homogenous Poisson process.
//'
//'\code{Log_likelihood} returns the value of the log-likelihood.
//'
//' @param data  Data frame with all records.
//' @param a,b    Parameters of the model.
//'
//' @return It returns the log likelihood.
//' @references
//' @examples
//' df<- data.frame(t_1=c(10,60,100),t_2=c(20,80,130),obs=c(2,3,4))
//' Leaf_emrgence_rate(df,0.062,0.056)
//' @export
// [[Rcpp::export]]
double Log_likelihood(DataFrame data, double a, double b){
  double lik=0;
  NumericVector obs = data["obs"];
  NumericVector t1 = data["t_1"];
  NumericVector t2 = data["t_2"];
  NumericVector h_t, r_t;
    h_t = Sum_Leaf_emergence_rate(t1,t2,a,b);
    r_t = a*(1 + b*cos(2*M_PI/365.0*(t2-15)));
    lik = -sum(h_t) + sum((obs-1)*log(h_t)) + sum(log(r_t))- sum(log(factorial(obs-1)));

  return(lik);
}

//' @export
// [[Rcpp::export]]
double Prior(double a, double a_p, double b, double b_p){
  return((-a/a_p - b/b_p)/0.5);
}
/*======================================================================================================================================
                             MCMC
 ======================================================================================================================================*/
//' The mcmc sampling of the non-homogenous Poisson process.
//'
//'\code{mcmc_leaf} returns the sample from the posterior distribution.
//'
//' @param data  Data frame with all records.
//' @param a0,b0    Innitial parameters of the model.
//'
//' @return It returns sample from the posterior distribution.
//' @references
//' @examples examples/africa_landscape_example.R
//' @export
// [[Rcpp::export]]
DataFrame mcmc_leaf(DataFrame data, double a0, double b0, int samp=1000){
  double q, x, pacc, tmp, qr, lp,lp_n;
  NumericVector a(samp);
  NumericVector b(samp);
  a[0]=a0;
  b[0]=b0;
  lp = Prior(a0,b0,a0,b0);
  for( int i=1;i<samp;i++){
// Update a
    q=.01*Rcpp::rnorm(2,0,1)[0];
    x=a[i-1]*exp(q);
    pacc=Log_likelihood(data, x, b[i-1]) - Log_likelihood(data, a[i-1], b[i-1]) + q;
    if(pacc>log(Rcpp::runif(2,0,1)[0])){
      a[i]=x;
    }
    else{
      a[i]=a[i-1];
    }

    // Update b
    q=.1*Rcpp::rnorm(2,0,1)[0];
    tmp =  exp(q)*(b[i-1]-0.)/(1-b[i-1]);
    x = 0. + (1-0.)*tmp/(1+tmp);
    qr = 2*log((x-0.)/b[i-1] - 0.) -q;
    pacc=Log_likelihood(data, a[i], x) - Log_likelihood(data, a[i], b[i-1]) + qr;
    // if(x<1){
    //   pacc=Log_likelihood(data, a[i], x) - Log_likelihood(data, a[i], b[i-1]) + qr;
    // }
    // else{
    //   pacc = -100000000000000000000000000000000000000000000000000000000000000.0;
    // }

    if(pacc>log(Rcpp::runif(2,0,1)[0])){
      b[i]=x;
    }
    else{
      b[i]=b[i-1];
    }
  }
  return(DataFrame::create(Named("a")=a, Named("b")=b));
}
