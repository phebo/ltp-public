// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppNumerical.h>
#define PI 3.141592653589793238462643383280
using namespace Rcpp;
using namespace Numer;

class integrand: public Func
{
private:
  double lz;
  double s;
  double mux;
  double muy;
  double sig;
  double rho;
public:
  integrand(double lz_, double s_, double mux_, double muy_, double sig_, double rho_) :
    lz(lz_), s(s_), mux(mux_), muy(muy_), sig(sig_), rho(rho_) {}
  
  double operator()(const double& v) const
  {
    double x = lz - std::log(v) - mux;
    double y = s * v - muy;
    return std::exp(-(x*x - 2*rho*sig*x*y + sig*sig*y*y) /
                    (2*(1-rho*rho)*sig*sig));
  }
};

// [[Rcpp::export]]
NumericVector dnln(NumericVector x, double mux, double muy, double sig, double rho)
{
  const double lower = 0;
  const double upper = R_PosInf;
  double err_est;
  int err_code;
  int n = x.size();
  NumericVector lz = Rcpp::log(Rcpp::abs(x));
  NumericVector s = ifelse(x >= 0.0, 1.0, -1.0);
  NumericVector out = NumericVector(n);
  for(int i =0; i<n; i++) {
    integrand f(lz[i], s[i], mux, muy, sig, rho);
    out[i] = integrate(f, lower, upper, err_est, err_code) /
      (2 * PI * sig * std::sqrt(1-rho*rho));
  }
  out = ifelse(x != 0.0, out / Rcpp::abs(x),
               exp(-0.5*muy*muy+rho*sig*muy+0.5*(1-rho*rho)*sig*sig-mux) /
                 sqrt(2*PI));
  return out; 
}

