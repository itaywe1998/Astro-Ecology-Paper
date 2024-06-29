/*
 Copyright (C) 2024 Itay Weintraub ----
 This program comes with ABSOLUTELY NO WARRANTY. This is free software, and
 you are welcome to redistribute it under certain conditions. for details,
 see the GNU General Public License Agreement (in the file COPYING.txt).
 */

#include <Rcpp.h>
#include <math.h>
#include <iostream>
#include <string>

using namespace std;
using namespace Rcpp;
/* Apply twice continuously differentiable smoothed step function to a number x
 Input:
 - x: Distance from pole, measured in units of the pole-to-equator distance
 Output:
 - 0 if x < 0; 10*x^3-15*x^4+6*x^5 if 0 <= x <= 1; otherwise 1 */
// Do not erase Rcpp export statement so the functions will callable also from R
// [[Rcpp::export]]
double smoothstep(double x) {
  double y;
  if (x<0.0) y=0.0;
  else if (x>1.0) y=1.0;
  else y=x*x*x*(10.0+x*(-15.0+6.0*x));
  return(y);
}
/* Temperature as a function of space (x), time (t), and some climate parameters
 Input:
 - x: Vector of distances from pole, in units of the pole-to-equator distance
 - t: Time at which temperatures are evaluated (climate change starts at t = 0)
 - tE: Time at which climate change ends (so it lasts from t = 0 to t = tE)
 - C: Projected temperature increase 
 - T0: Initial  temperature 
 Output:
 - Vector of temperatures at each location x */
// [[Rcpp::export]]
double Temp(double t,double Tmin,double C, double tE) {
   return Tmin+C*smoothstep(t/tE);
}
/* Right-hand side of dynamical equations
 Input:
 - time: Time at which function is evaluated (explicit time-dependence)
 - state: Vector of state variables, with 2 entries. The first is population
 density n, the second is the optimal temperature trait m.
 - pars: Model parameters, given as members of a list
 Output:
 - The derivatives of the densities and trait means, as a vector in a list */
// [[Rcpp::export]]
List eqs(double time, NumericVector state, List pars) {
  // Parameters
  double nmin=pars["nmin"];
  double sigma=pars["sigma"], kappa=pars["kappa"];
  double T0=pars["T0"],tE=pars["tE"], C=pars["C"];
  double rho=pars["rho"],v = pars["v"];
  // Variables
  double ef,g, h2;
  double n,m,T;
  NumericVector dvdt(2);
  n=state[0]; 
  m=state[1]; 
  if(n<nmin){
    cout<<"Population died"<<endl;
    cout<<n<<endl;
    throw range_error(to_string(time));
  }
  T=Temp(time,  T0, C, tE); // Vector of temperatures
  ef=rho*exp(-(T-m)*(T-m)/(2.0*sigma*sigma))/sigma;
  g=ef*v*(T-m)/sigma;
  h2=0.5; // Heritability
  // Assign calculated rates to vector of derivatives for output
  dvdt[0]=n*(ef-kappa);
  dvdt[1]=h2*g;
  return(List::create(dvdt));
}
