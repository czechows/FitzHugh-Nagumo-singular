
/* -----------------------------------------------------------------------------------------
 * This is a header file to fhn.cpp providing numerical simulations necessary
 * for the computer assistance of proofs of periodic orbits in the FitzHugh-Nagumo system.
 * We find values of v,theta for which there exist heteroclinic connections between equilibria
 * in the fast subsystem via shooting from stable/unstable manifolds of equilibria
 * to the Poincare section given in between.
 * The validity of the proofs does not depend on validity of numerics given below.
 * If the full FitzHugh-Nagumo vector field is alternated in any other way than changing
 * parameters theta or eps, one needs to manually readjust the fast vector field
 * and the fast reversed vector field given in the class FhnBifurcation
 * ----------------------------------------------------------------------------------------*/


/* ----------------------------------------------------------------------------------------- */
/* ---------------------------- FAST SUBSYSTEM NUMERICS ------------------------------------ */
/* ----------------------------------------------------------------------------------------- */

class FhnBifurcation
{
public:
  DMap vectorField;
  DMap vectorFieldRev;    // reversed (minus) vector field for backward integration
  DTaylor solver;
  DTaylor solverRev;
  DVector EqU;            // "upper" equilibrium (guess)
  DVector EqD;            // "lower" equilibrium (guess)
  double DISP;             // displacement from the equilibra in stable/unstable direction
  DAffineSection section;
  DPoincareMap pm;
  DPoincareMap pmRev;     // reversed Poincare map for backward integration
  bool dir;               // direction - do we go from EqU to EqD or other way round
  bool homoclinic;        // 1 if we look for theta instead of v for the homoclinic connection
  double homTheta;           // used only for the homoclinic case

  FhnBifurcation (int order, double& _theta, const DVector& _EqU, const DVector& _EqD, double _DISP, bool _dir = 1, bool _homoclinic=0) 
    : vectorField("par:theta,w;var:u,v;fun:v,(2/10)*(theta*v+u*(u-1)*(u-(1/10))+w);"), // vector field is u'=w, w'=0.2*(theta*w +u*(u-1)*(u-0.1)+v, v is parameter
      vectorFieldRev("par:theta,w;var:u,v;fun:-v,(-2/10)*(theta*v+u*(u-1)*(u-(1/10))+w);"), //  minus vector field for reverse integration
      solver(vectorField,order),
      solverRev(vectorFieldRev,order),
      EqU(_EqU),
      EqD(_EqD),
      DISP(_DISP),
      section( DVector({0.2, 0. }), DVector({-1.,0.})), // an arbitrary choice of coordinate section here
      pm(solver,section),
      pmRev(solverRev,section),
      dir(_dir),
      homoclinic(_homoclinic),
      homTheta(_theta)
  {
    vectorField.setParameter("theta",_theta);
    vectorFieldRev.setParameter("theta",_theta);
  }

  DVector Eq_correct(DVector& guess, double w)              // corrects initial guesses of u so they are closer to real equilibria, w is always 0
  {
    vectorField.setParameter("w", w);
    double error = 1.;
    double result = guess[0];
    double oldresult = result;
    while(error > accuracy)
    {
      oldresult = result;
      result = result - ( vectorField(DVector({result, 0.}) )[1] / vectorField[ DVector({result, 0.}) ][1][0] );  // Newton algorithm to calculate zeroes of the vector field - v is always 0.,
                                                                                                                  // derivative is of the second equation wr. to first variable - u
      error = abs(oldresult - result);
    }
    DVector new_Eq(2);
    new_Eq[0] = result;
    new_Eq[1] = 0.;
    return new_Eq;
  }

  DMatrix J_correct(DVector Eq, double w)
  {
    vectorField.setParameter("w", w); 
    return( vectorField[Eq] );
  }

  double v_function(DVector guessEqU, DVector guessEqD, double w)  // returns distance in v variable on the Poincare section between integrated displacement in unstable direction from EqU
                                                                   // and integrated displacement in stable direction from EqD if dir = 1 or vice-versa elsewise
  {
    double return_time = 1.;
    vectorField.setParameter("w", w); 
    vectorFieldRev.setParameter("w", w);

    DMatrix EigenvectU_real(2,2),
            EigenvectU_im(2,2);
    DVector EigenvalU_real(2),
            EigenvalU_im(2);
  
    DMatrix EigenvectD_real(2,2),
            EigenvectD_im(2,2);
    DVector EigenvalD_real(2),
            EigenvalD_im(2);   
    
    computeEigenvaluesAndEigenvectors( J_correct(guessEqU,w), EigenvalU_real, EigenvalU_im, EigenvectU_real, EigenvectU_im );
    computeEigenvaluesAndEigenvectors( J_correct(guessEqD,w), EigenvalD_real, EigenvalD_im, EigenvectD_real, EigenvectD_im );

    if(EigenvalU_real[0]*EigenvalU_real[1] >= 0.)
      throw "EIGENVALUES OF FAST SUBSYSTEM AT STATIONARY POINTS NOT OF OPPOSITE SIGNS! \n";
    if(EigenvalD_real[0]*EigenvalD_real[1] >= 0.)
      throw "EIGENVALUES OF FAST SUBSYSTEM AT STATIONARY POINTS NOT OF OPPOSITE SIGNS! \n";
    if(EigenvalU_im[0] != 0.)
      throw "EIGENVALUES OF FAST SUBSYSTEM AT STATIONARY POINTS NOT REAL! \n";
    if(EigenvalD_im[0] != 0.)
      throw "EIGENVALUES OF FAST SUBSYSTEM AT STATIONARY POINTS NOT REAL! \n";

    // we want to have first eigenvector stable, second unstable. We do not swap eigenvalues since we don't use them later. Slows down everything by ~30s
    // commented because we already have the right ones
 /*   if( EigenvalU_real[0] > EigenvalU_real[1] )
      matrixAlgorithms::columnExchange( EigenvectU_real, 1, 2 );

    if( EigenvalD_real[0] > EigenvalD_real[1] )
      matrixAlgorithms::columnExchange( EigenvectD_real, 1, 2 );
*/

    // here one has to play manually with pluses and minuses so we are on the right side of stable/unstable manifolds and we catch the right eigenvectors
    if(dir == 1) 
        return pm(guessEqU + EigenvectU_real.column(1)*DISP, return_time)[1] - pmRev(guessEqD - EigenvectD_real.column(0)*DISP, return_time)[1]; 
    else
        return pmRev(guessEqU + EigenvectU_real.column(0)*DISP, return_time)[1] - pm(guessEqD - EigenvectD_real.column(1)*DISP, return_time)[1];   
  }

  double w_correct(double w) // secant method to correct w to the bifurcation point, as side effect corrects equilibria EqU and EqD to right positions 
  {
    double error = 1.;

    double w0 = w + 1e-4;
    double w1 = w;
    double w_temp = w;
 
    DVector EqU0 = Eq_correct(EqU, w0);
    DVector EqU1 = Eq_correct(EqU, w1);
    DVector EqD0 = Eq_correct(EqD, w0);
    DVector EqD1 = Eq_correct(EqD, w1);
 
   while(error > accuracy)
   {
    w_temp = w1;
    w1 = w1 - v_function(EqU1, EqD1, w1)*( (w1-w0) / ( v_function(EqU1, EqD1, w1) - v_function(EqU0, EqD0, w0) ) );
    EqU0 = Eq_correct( EqU0, w_temp );
    EqD0 = Eq_correct( EqD0, w_temp );
    EqU1 = Eq_correct( EqU1, w1 );
    EqD1 = Eq_correct( EqD1, w1 );
    w0 = w_temp;

    error = abs( v_function(EqU1, EqD1, w1) );
   }

   EqU = EqU1;
   EqD = EqD1;
    
   return w1;
  }

  double theta_correct() // corrects homTheta to one for which a heteroclinic connection between (0,0) and an equilibrium on the upper branch of the slow manifold exists
                         // as a side effect corrects EqD to (0,0) and EqU to a correct position on the equilibrium upper branch, returns homTheta
                         // all happens for w=0 (one needs manual readjustments if that happened for other w's)
  {
    if( !homoclinic )
      throw "Homoclinic option needs to be enabled";

    double error = 1.;

    double theta0 = homTheta + 1e-4;
    double theta1 = homTheta;
    double theta_temp = theta1;

    vectorField.setParameter("theta",theta0); 
    DVector EqU0 = Eq_correct(EqU, 0.);  
    DVector EqD0 = Eq_correct(EqD, 0.);  

    vectorField.setParameter("theta",theta1);
    DVector EqU1 = Eq_correct(EqU, 0.);  
    DVector EqD1 = Eq_correct(EqD, 0.);  


    while(error > accuracy)
    {
      theta_temp = theta1;

      vectorField.setParameter("theta",theta0);
      vectorFieldRev.setParameter("theta",theta0);
      EqU0 = Eq_correct(EqU, 0.);  
      double v0 = v_function( EqU0, EqD0, 0. );

      vectorField.setParameter("theta",theta1);
      vectorFieldRev.setParameter("theta",theta1);
      EqU1 = Eq_correct(EqU, 0.);  
      double v1 = v_function( EqU1, EqD1, 0. );

      theta1 = theta1 - v1*( ( theta1 - theta0 ) / ( v1 - v0 ) );
 
      vectorField.setParameter("theta",theta_temp);
      EqU0 = Eq_correct(EqU0, 0.);  
      EqD0 = Eq_correct(EqD0, 0.);  
 
      vectorField.setParameter("theta",theta1);
      EqU1 = Eq_correct(EqU1, 0.);  
      EqD1 = Eq_correct(EqD1, 0.);  

      theta0 = theta_temp;

      vectorField.setParameter("theta",theta1);
      error = abs( v_function(EqU1, EqD1, 0.) );
    }

    EqU = EqU1;
    EqD = EqD1;   // again, we know that this is (0,0) but we recompute numerically so code will work for other systems with non-explicit stationary points
    homTheta = theta1;

    return theta1;
  }
};


void GammaQuad_correct( const interval& _theta, IVector& _GammaUL, IVector& _GammaDL, IVector& _GammaUR, IVector& _GammaDR ) // corrects original guesses of Gammas for given theta
{
  double theta( _theta.leftBound() );
  double DISP(1e-12);

  double wR( _GammaUR[2].leftBound() );
  double wL( _GammaUL[2].leftBound() );

  DVector EqUL(2),
          EqUR(2),
          EqDR(2),
          EqDL(2);      // equilibria are first two coordinates of each Gamma, they will be corrected by Newton inside the program

  EqUL[0] = _GammaUL[0].leftBound();
  EqUL[1] = _GammaUL[1].leftBound();
  EqUR[0] = _GammaUR[0].leftBound();
  EqUR[1] = _GammaUR[1].leftBound();

  EqDR[0] = _GammaDR[0].leftBound();
  EqDR[1] = _GammaDR[1].leftBound();
  EqDL[0] = _GammaDL[0].leftBound();
  EqDL[1] = _GammaDL[1].leftBound();


  FhnBifurcation BifR(order, theta, EqUR, EqDR, DISP);
  FhnBifurcation BifL(order, theta, EqUL, EqDL, DISP, 0); 

  double wR_c( BifR.w_correct(wR) );        // corrected w values & equlibria coordinates
  double wL_c( BifL.w_correct(wL) );


  DVector EqUL_c( BifL.EqU ),
          EqUR_c( BifR.EqU ),
          EqDR_c( BifR.EqD ), 
          EqDL_c( BifL.EqD );

  _GammaUL[0] = EqUL_c[0];                  // assignment back to Gammas
  _GammaUL[1] = EqUL_c[1];
 
  _GammaUR[0] = EqUR_c[0];
  _GammaUR[1] = EqUR_c[1];
 
  _GammaDL[0] = EqDL_c[0];
  _GammaDL[1] = EqDL_c[1];
 
  _GammaDR[0] = EqDR_c[0];
  _GammaDR[1] = EqDR_c[1];

  _GammaUL[2] = _GammaDL[2] = wL_c;
  _GammaUR[2] = _GammaDR[2] = wR_c;
}

void GammaHom_correct( interval& _theta, IVector& _GammaUL, IVector& _GammaDL, IVector& _GammaUR, IVector& _GammaDR ) // corrects original guesses of Gammas for given theta, updates theta
{
  double theta( _theta.leftBound() );
  double DISP(1e-5);


  DVector EqUL(2),
          EqDL(2);      // equilibria are first two coordinates of each Gamma, they will be corrected by Newton inside the program

  EqUL[0] = _GammaUL[0].leftBound();
  EqUL[1] = _GammaUL[1].leftBound();
 
  EqDL[0] = _GammaDL[0].leftBound();
  EqDL[1] = _GammaDL[1].leftBound();

  FhnBifurcation BifHom(order, theta, EqUL, EqDL, DISP, 0, 1); // the homoclinic option on, dir = 0
  _theta = interval( BifHom.theta_correct() );


  DVector EqUL_c( BifHom.EqU ),
          EqDL_c( BifHom.EqD );
 
  _GammaUL[0] = EqUL_c[0];                  // assignment back to Gammas
  _GammaUL[1] = EqUL_c[1];
  _GammaUL[2] = 0.;

  _GammaDL[0] = EqDL_c[0];
  _GammaDL[1] = EqDL_c[1];
  _GammaDL[2] = 0.;

  GammaQuad_correct( _theta, _GammaUL, _GammaDL, _GammaUR, _GammaDR );  // only to fix GammaUR and GammaDR now
}


