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
 - Cmax: Projected temperature increase at North Pole
 - Cmin: Projected temperature increase at equator
 - Tmax: Initial mean temperature at equator
 - Tmin: Initial mean temperature at North Pole
 Output:
 - Vector of temperatures at each location x */
// [[Rcpp::export]]
NumericVector Temp(NumericVector x,double Tmax, double Tmin) {
    return((Tmax-Tmin)*x+Tmin);
}
/* Right-hand side of dynamical equations
 Input:
 - time: Time at which function is evaluated (explicit time-dependence)
 - state: Vector of state variables, with 2*S*L entries, where S is the number
 of species and L the number of patches. The first S*L entries are the
 densities, the second S*L entries are the trait means.
 - pars: Model parameters, given as members of a list
 Output:
 - The derivatives of the densities and trait means, as a vector in a list */
// [[Rcpp::export]]
List eqs(double time, NumericVector state, List pars) {
  // Parameters
  int S=pars["S"], L=pars["L"], lT=pars["lT"];
  double eta=pars["eta"], nmin=pars["nmin"], vbar=pars["vbar"];
  double aw=pars["aw"], bw=pars["bw"], kappa=pars["kappa"];
  NumericVector d=pars["d"], v=pars["v"], rho=pars["rho"];
  NumericMatrix mig=pars["mig"];
  NumericMatrix T_kozai=pars["T_kozai"];
  // Variables
  int i, j, k, l,ki, Tmin = 0, Tmax = 0;
  double sumgr, summig, sigma, ef, b, bsumgr, bsummig, g, q, dm, h2;
  NumericMatrix n(S,L), m(S,L), alpha(S,S), beta(S,S);
  NumericVector dvdt(2*S*L), x(L), T(L);
  // Assign state variables into matrices n and m; calculate local temperatures
  for(i = 0; i < S; i++){
    for (k=0; k<L; k++) {
        x[k]=k/((double)L-1.0); // Patch k's distance from pole
          n(i,k)=state[i+k*S]; // Density of species i in patch k
          if (n(i,k)<1.0e-10) n(i,k)=0.0; // Extinction threshold
          m(i,k)=state[S*L+i+k*S]; // Trait mean of species i in patch k

      }
  }
  
  double maxn= max(n);
  double meann = mean(n);
  double threshold = nmin;
  
  if(maxn<threshold || meann<0){
    cout<<maxn<<endl;
    if(maxn< threshold) cout<<"Population died"<<endl;
    throw range_error(to_string(time));
  }


  int prev, next,ns=5; // interpolating the given T table
  double fraction;
  for(ki=0; ki<lT; ki++){
    prev = T_kozai(ki,0);
    if (ki == lT-1){
      next = prev;
    }
    else{
      next = T_kozai(ki+1,0);
    }
    if(time >= prev && time<=next){
      fraction = ceil((time - prev) / (next-prev)*ns)/ns; //interpolating only ns points between table steps
      Tmin = (T_kozai(ki+1,1)-T_kozai(ki,1))*fraction + T_kozai(ki,1); //linear interpolation
      Tmax = (T_kozai(ki+1,2)-T_kozai(ki,2))*fraction + T_kozai(ki,2);
      break;
    }
  }
  T=Temp(x,  Tmax, Tmin); // Vector of temperatures

  for (k=0; k<L; k++) {
    // If we have temperature-dependent competition:
    for (i=0; i<(S-1); i++) {
      alpha(i,i)=1;
      for (j=i+1; j<S; j++) {
        dm=m(j,k)-m(i,k);
        alpha(i,j)=exp(-dm*dm/(eta*eta));
        alpha(j,i)=alpha(i,j);
        beta(i,j)=2.0*v[i]*alpha(i,j)*dm/(eta*eta);
        beta(j,i)=-beta(i,j)*v[j]/v[i];
      }
      alpha(S-1,S-1)=1;
    }
    for (i=0; i<S; i++) {
      sumgr=0.0;
      bsumgr=0.0;
      // Species interaction terms in density and then trait evolution equations
      for (j=0; j<S; j++) {
        sumgr+=-n(i,k)*alpha(i,j)*n(j,k);
        bsumgr+=beta(i,j)*n(j,k);
      }
      summig=0.0;
      bsummig=0.0;
      // Dispersal terms in density and then trait evolution equations
      for (l=0; l<L; l++) {
        summig+=mig(k,l)*n(i,l)-n(i,k)*mig(l,k);
        bsummig+=mig(k,l)*n(i,l)*(m(i,l)-m(i,k))/(n(i,k)+1.0e-10);
      }
      // Growth terms in the equations
      summig*=d[i];
      bsummig*=d[i];
      sigma=bw-aw*m(i,k);
      ef=rho[i]*exp(-(T[k]-m(i,k))*(T[k]-m(i,k))/(2.0*sigma*sigma))/sigma;
      b=ef-kappa;
      g=ef*v[i]*(T[k]-m(i,k))/(sigma*sigma);
      q=v[i]*smoothstep(n(i,k)/nmin);
      h2=q/(q+vbar); // Heritability
      // Assign calculated rates to vector of derivatives for output
      dvdt[i+k*S]=(n(i,k)*b+sumgr)*smoothstep(n(i,k)/nmin)+summig;
      dvdt[S*L+i+k*S]=h2*(g-bsumgr+bsummig);
      // Periodic boundary conditions
      if (k==0) {
        dvdt[i]+=d[i]*(mig(0,1)*n(i,1)-mig(1,0)*n(i,0));
        dvdt[S*L+i]+=d[i]*h2*mig(0,1)*n(i,1)*(m(i,1)-m(i,0))/(n(i,0)+1.0e-10);
      }
      if (k==(L-1)) {
        dvdt[i+k*S]+=d[i]*(mig(k,k-1)*n(i,k-1)-mig(k-1,k)*n(i,k));
        dvdt[S*L+i+k*S]+=d[i]*h2*mig(k,k-1)*n(i,k-1)*
          (m(i,k-1)-m(i,k))/(n(i,k)+1.0e-10);
      }
    }
  }
  return(List::create(dvdt));
}
