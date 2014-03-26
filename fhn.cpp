#include <iostream>
#include "capd/capdlib.h"
#include "capd/dynsys/DiscreteDynSys.h"

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

  cout.precision(9);

  interval theta = interval(61.)/100.; // theta = 0.53 also works
  interval eps = interval(0.,7.)/1e5; //interval(0.,1.)/2e6;
  bool verbose = 1; 
  bool with_params = 0;
  
  FhnVerifyExistenceOfPeriodicOrbit( theta, eps, verbose, with_params );
 // FhnVerifyExistenceOfPeriodicOrbit( theta, eps, verbose );

//  theta = interval(63.)/100;
 
//  FhnVerifyExistenceOfPeriodicOrbit( theta, eps, verbose );

  return 0;
} 

