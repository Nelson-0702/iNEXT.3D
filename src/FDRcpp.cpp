#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List rcpp_hello_world() {

    CharacterVector x = CharacterVector::create( "foo", "bar" )  ;
    NumericVector y   = NumericVector::create( 0.0, 1.0 ) ;
    List z            = List::create( x, y ) ;

    return z ;
}

// [[Rcpp::export]]
NumericVector RFD(NumericMatrix x , int n  , double m , NumericVector q,double V_bar) {
  // x is x*i_vi; return length q
  int nrow = x.nrow();
  int qlength = q.length();
  NumericVector hhat(m);
  NumericVector out(qlength);
  //Rcout << "out in cpp: " << out << std::endl;
  for (int k = 0; k < m ; k++) {
    for (int i = 0; i < nrow; i++) {
      if ( x(i,0) >= k+1 && x(i,0) <= n-m+k+1 )
      {
        hhat[k] +=  x(i,1)*exp(Rf_lchoose(x(i,0), k+1)+Rf_lchoose(n-x(i,0), m-k-1)-Rf_lchoose(n, m)) ;
      }
      else
      {
        hhat[k] += 0 ;
      }
    }
  }
  //Rcpp::Rcout << "hhat in cpp: " << hhat << std::endl;
  for (int j = 0; j < qlength; j++ ){  
    for(int k = 0; k < m; k++){
      if(q[j] == 0){
        out[j] = hhat[k] + out[j];
      }else if(q[j] == 1){
        //Rcout << "q1 in cpp: " <<log ( (k+1) )<< std::endl;
        out[j] = -( (k+1) / (m*V_bar) ) * log ( (k+1) / (m*V_bar) ) * hhat[k] + out[j];
      }else if(q[j] == 2){
        out[j] = pow( ( (k+1) / (m*V_bar) ),2) * hhat[k] + out[j];
      }else{
        out[j] = pow( ( (k+1) / (m*V_bar) ),q[j]) * hhat[k] + out[j];
      }
    }
  }
  //Rcout << "out in cpp: " << out << std::endl;
  for(int j = 0; j < qlength; j++ ){
    if(q[j] == 0){
      out[j] = out[j] ;
    }else if(q[j] == 1){
      out[j] = exp(out[j]);
    }else if(q[j] == 2){
      out[j] = 1 / out[j];
    }else{
      out[j] = pow( (out[j]) , 1/(1-q[j]) );
    }
  }
  return out;
}

// [[Rcpp::export]]
double FDq0(double n, double f1, double f2, double h1,  double h2, double A) {
  double ans;
  if(h1+h2 == 0) {ans = 0;
  } else if (h2>0){
    ans = ((n-1)/n)*(pow(h1,2)/(2*h2));
  } else{
    ans = ((n-1)/n)*(h1*h1/2);
  }
  
  return(ans);
}
// [[Rcpp::export]]
double FDq1_1(double n, double h1, double A) {
    double q1 = 0;
    double h_est_2 = 0;
  
  if(A==1){
    h_est_2 = 0;
  }else{
    for(int r = 1; r < n; r++){
      q1 = q1 + pow((1-A),r)/r;
    }
    

  h_est_2 = (h1 / n) * pow(1 - A, (-n + 1)) * round((-log(A) - q1) * pow(10, 12)) / pow(10, 12);
    
    
  }
  return(h_est_2);
}
// [[Rcpp::export]]
double FDq2(NumericMatrix tmpxv, double n) {
  double ans = 0;
  for(int i = 0; i < tmpxv.nrow(); i++){
    ans = ans + (tmpxv(i,2)*tmpxv(i,1)*tmpxv(i,0)*(tmpxv(i,0)-1)/n/(n-1));
  }
  //Rcout << "ans: " << ans;
  ans = 1/ans;
  return(ans);
}
