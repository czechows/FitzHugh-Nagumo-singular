#include <iostream>
#include "capd/capdlib.h"
#include "capd/dynsys/DiscreteDynSys.h"
#include "time.h"

using std::cout;

using namespace capd;
using namespace matrixAlgorithms;
using namespace dynsys;

// in all diagonalizations, unless otherwise stated, first variable is stable/entry second unstable/exit third (if present) central
// in the paper we later changed the first variable to be unstable/exit second stable/entry third (if present) central 
// this was to match some conventions from previous papers on h-sets

IMap *Fhn_vf;
IMap *Fhn_vf_rev;
IMap *Fhn_vf_withParams;
IMap *Fhn_vf_withParams_rev;

const interval EPS = interval(1./1e15);  // small number greater than zero for coverings
const double accuracy = 1e-12;           // accuracy for nonrigorous numerics (i.e. approximation of the slow manifold)
const int order = 18;                    // order for all the Taylor integrators (high is fast)

#include "numerics.hpp"   // Warning! When changing the vector field, one needs to make manual changes in this header file (class FhnBifurcation)!
#include "auxiliaries.hpp"
#include "segments.hpp"
#include "poincare.hpp"
#include "block.hpp"  // Warning! parameter a=1/10 hardcoded there
#include "proof.hpp" // the proof for the periodic orbit
#include "homoclinic_proof.hpp" // the proof for the homoclinic orbit

// ---------------------------------------------------------------------------------
// ----------------------------------- MAIN ----------------------------------------
// ---------------------------------------------------------------------------------




int main(){

  time_t start1,end1;

  Fhn_vf = new IMap("par:theta,eps;var:u,v,w;fun:v,(2/10)*(theta*v+u*(u-1)*(u-(1/10))+w),(eps/theta)*(u-w);"); 
  // FitzHugh-Nagumo vector field is u'=v, v'=0.2*(theta*v +u*(u-1)*(u-0.1)+w, w'= eps/theta * (u-w)
  Fhn_vf_rev = new IMap("par:theta,eps;var:u,v,w;fun:-v,(-2/10)*(theta*v+u*(u-1)*(u-(1/10))+w),(-eps/theta)*(u-w);"); 
  // reversed field for backward integration
  Fhn_vf_withParams = new IMap("var:u,v,w,theta,eps;fun:v,(2/10)*(theta*v+u*(u-1)*(u-(1/10))+w),(eps/theta)*(u-w),0,0;"); 
  // the same vector field with parameters as variables of velocity 0
  Fhn_vf_withParams_rev = new IMap("var:u,v,w,theta,eps;fun:-v,(-2/10)*(theta*v+u*(u-1)*(u-(1/10))+w),(-eps/theta)*(u-w),0,0;"); 
  // again, the reversed vector field with parameters of velocity 0

  time (&start1);
  cout.precision(15);

  // A PROOF FOR THE HOMOCLINIC ORBIT

  interval thetaGuess = interval(111.)/100.;   // we only try to prove the fast wave, the slow one does not come from the singular perturbation [KSS]
  bool verbose = 1; 
  bool with_params = 0; // allowing parameters to evolve as variables with velocity 0 does not improve significantly the results -- OUT OF SUPPORT!
 
  interval eps = interval(0.,5.)/1e5;  
  FhnVerifyExistenceOfHomoclinicOrbit( thetaGuess, eps, verbose, with_params );
 
  time (&end1);
  double dif1 = difftime( end1, start1 );
  cout << "Elapsed time for the homoclinic orbit proof is " << dif1 << " seconds. \n";



  // THE PERIODIC ORBIT PROOF FROM THE ARXIV PAPER
   interval theta = interval(61.)/100.;  
   eps = interval(0.,1.)/1e4;  
  
  time_t start2,end2;
  time (&start2);

  FhnVerifyExistenceOfPeriodicOrbit( theta, eps, verbose, with_params );

  time (&end2);
  double dif2 = difftime( end2, start2 );
  cout << "Elapsed time for the periodic orbit proof for parameter range eps = " << eps << " is " << dif2 << " seconds. \n";


  time_t start3,end3;
  time (&start3);
  eps = interval("1e-4","1.5e-4");  
  
  FhnVerifyExistenceOfPeriodicOrbit( theta, eps, verbose, with_params );
  
  time (&end3);
  double dif3 = difftime( end3, start3 );
  cout << "Elapsed time for the periodic orbit proof for parameter range eps = " << eps << " is " << dif3 << " seconds. \n";


  // other thetas
/*
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
 /* 
  eps = interval(1.5,2.)/1e4;  
  FhnVerifyExistenceOfPeriodicOrbit( theta, eps, verbose, with_params );
  // this already fails


  // other thetas:
  */


  return 0;
} 

