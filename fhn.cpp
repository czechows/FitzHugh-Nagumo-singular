#include <iostream>
#include "capd/capdlib.h"
#include "capd/dynsys/DiscreteDynSys.h"
#include "time.h"

using std::cout;

using namespace capd;
using namespace matrixAlgorithms;
using namespace dynsys;

// in all diagonalizations, unless otherwise stated, first variable is stable second unstable third (if present) neutral

const interval EPS = interval(1./1e15);  // small number greater than zero for coverings
const double accuracy = 1e-12;           // accuracy for nonrigorous numerics (i.e. approximation of the slow manifold)
const int order = 18;                    // order for all the Taylor integrators (high is fast)

IMap Fhn_vf("par:theta,eps;var:u,w,v;fun:w,(2/10)*(theta*w+u*(u-1)*(u-(1/10))+v),(eps/theta)*(u-v);"); 
// FitzHugh-Nagumo vector field is u'=w, w'=0.2*(theta*w +u*(u-1)*(u-0.1)+v, v'= eps/theta * (u-v)
IMap Fhn_vf_rev("par:theta,eps;var:u,w,v;fun:-w,(-2/10)*(theta*w+u*(u-1)*(u-(1/10))+v),(-eps/theta)*(u-v);"); 
// reversed field for backward integration
IMap Fhn_vf_withParams("var:u,w,v,theta,eps;fun:w,(2/10)*(theta*w+u*(u-1)*(u-(1/10))+v),(eps/theta)*(u-v),0,0;"); 
// the same vector field with parameters as variables of velocity 0
IMap Fhn_vf_withParams_rev("var:u,w,v,theta,eps;fun:-w,(-2/10)*(theta*w+u*(u-1)*(u-(1/10))+v),(-eps/theta)*(u-v),0,0;"); 
// again, the reversed vector field with parameters of velocity 0

#include "numerics.hpp"   // Warning! When changing the vector field, one needs to make manual changes in this header file (class FhnBifurcation)!
#include "poincare.hpp"
#include "segments.hpp"
#include "proof.hpp"


// ---------------------------------------------------------------------------------
// ----------------------------------- MAIN ----------------------------------------
// ---------------------------------------------------------------------------------




int main(){

  time_t start1,end1;

  time (&start1);
  cout.precision(15);

  interval theta = interval(61.)/100.;  
  interval eps = interval(0.,1.)/1e4;  
  bool verbose = 1; 
  bool with_params = 0; // allowing parameters to evolve as variables with velocity 0 does not improve significantly the results
  
  FhnVerifyExistenceOfPeriodicOrbit( theta, eps, verbose, with_params );

  time (&end1);
  double dif1 = difftime( end1, start1 );
  cout << "Elapsed time for the proof for parameter range eps = " << eps << " is " << dif1 << " seconds. \n";

  time_t start2,end2;
  time (&start2);
  eps = interval("1e-4","1.5e-4");  
  
  FhnVerifyExistenceOfPeriodicOrbit( theta, eps, verbose, with_params );
  
  time (&end2);
  double dif2 = difftime( end2, start2 );
  cout << "Elapsed time for the proof for parameter range eps = " << eps << " is " << dif2 << " seconds. \n";

 /* 
  eps = interval(1.5,2.)/1e4;  
  FhnVerifyExistenceOfPeriodicOrbit( theta, eps, verbose, with_params );
  // this already fails


  // other thetas:

  theta = interval(53.)/100;
  eps = interval(0.,1.)/2e4;

  FhnVerifyExistenceOfPeriodicOrbit( theta, eps, verbose, with_params );



  theta = interval(47.)/100;
  eps = interval(0.,1.)/2e4;

  FhnVerifyExistenceOfPeriodicOrbit( theta, eps, verbose, with_params );

  theta = interval(550., 554.)/1000;
  eps = interval(0.,1.)/2e4;

  FhnVerifyExistenceOfPeriodicOrbit( theta, eps, verbose, with_params );
  */


  return 0;
} 

