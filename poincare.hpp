/* -----------------------------------------------------------------------------------------
 * This is a header file to fhn.cpp providing Poincare map class for the Poincare maps
 * between branches of the slow manifold and auxiliaries for rigorous verification of
 * 2-dim covering relations.
 * ----------------------------------------------------------------------------------------*/


/* ------------------------------------------------------------------------------------ */
/* ---------------------------- POINCARE MAPS ----------------------------------------- */
/* ------------------------------------------------------------------------------------ */

class FhnPoincareMap 
                       // this class was originally used to integrate from section near one corner point to section near another. However, this functionality was removed,
                    // and now it only serves as a base class for midPoincareMap which integrates a section midway between the corner points
{
public:
  int dim;          // phase space dimension
  IMap vectorField;
  ITaylor solver;

  IMatrix P1;
  IMatrix P2;
 
  IVector GammaLeft1;
  IVector GammaRight1;
  IVector GammaLeft2;
  IVector GammaRight2;
  IVector GammaCenter1;
  IVector GammaCenter2;

  IVector y1vector;
  IVector section1CenterVector;
  IVector y2vector;
  IVector section2CenterVector;

  int div;         // number of subdivisions in each dimension (supports two, does not support subdivisions in parameter space) for integration of h-sets
  IVector params;   // vector of parameters


  FhnPoincareMap( IMap _vectorField, const FhnIsolatingSegment _segment1, const FhnIsolatingSegment _segment2, interval dir = interval(1.), int _div=1 ) 
    : dim( 3 ),
      vectorField( _vectorField ),                               // a 3d vector field, needs to have a parameter named "theta"
      solver( vectorField, order ),

      P1( _segment1.P ),
      P2( _segment2.P ),

      GammaLeft1( _segment1.GammaLeft ),
      GammaRight1( _segment1.GammaRight ),
      GammaLeft2( _segment2.GammaLeft ),
      GammaRight2( _segment2.GammaRight ),

      GammaCenter1( GammaLeft1/2. + GammaRight1/2. ),            // these are only for 'nonrigorous' numerics now, i.e. finding the mid section
      GammaCenter2( GammaLeft2/2. + GammaRight2/2. ),            // we always choose the segments so the center vector is one of the corner points

      y1vector( 0., -dir*_segment1.rightFace[1].rightBound(), 0. ),           // here we choose which one of two faces we integrate, "dir is reversed"                   
      section1CenterVector( P1*y1vector + GammaCenter1 ),

      y2vector( dir*_segment2.rightFace[0].rightBound(), 0., 0. ),            // same here
      section2CenterVector( P2*y2vector + GammaCenter2 ),

      div( _div ),
      params( 1 ) 
  {
    if( !(_segment1.leftFace == _segment1.rightFace && _segment2.leftFace == _segment2.rightFace) )
      throw "Poincare maps only implemented for segments with leftFace == rightFace \n";
    
    if( !(dir == interval(1.)) && !(dir == interval(-1.)) )
      throw "dir must be plus minus 1";

    if( !(_segment1.leftFace == -_segment1.leftFace && _segment2.leftFace == -_segment2.leftFace) )
      throw "Only segments with faces symmetric with respect to 0 allowed \n";
  }


 // #include "todo/fhnPMwithparams.hpp"
};





/* this is a derived class, which allows to integrate forward from one branch of the slow manifold and backward from the other
 * to verify forward/backward covering of a set on a section halfway between them. This is more efficient as eliminates possible nontransversal
 * intersections with sections which occur close to fixed points / slow manifolds. The midsection and induced coordinate system
 * is created by integrating the equation from section 1 to a temporary section halfway (in u coord.) between the slow manifold branches.
 * The actual midSection is chosen as to be orthogonal to the vector field at its center.
 * Then we integrate the variational equation to induce a coordinate system good for checking coverings.
 * The placement of the midsection and coordinate system does not need to be rigorous and could have been done on doubles, but we did it on intervals
 * to avoid code repetition. */

class midPoincareMap : public FhnPoincareMap
{
public:
  IVector midCenterVector;
  IMatrix midP;
  IAffineSection midSection;
  IMap vectorFieldRev;
  interval thetaRange; // has to be same as theta for the vector field (there is no .getParameter!)
  interval epsRange; // same as thetaRange
  
  midPoincareMap( IMap _vectorField, IMap _vectorFieldRev, const FhnIsolatingSegment _segment1, const FhnIsolatingSegment _segment2, interval _thetaRange, interval _epsRange, 
      interval dir = interval(1.), int _div=1 ) 
  : FhnPoincareMap( _vectorField, _segment1, _segment2, dir, _div ), // if dir = -1 we integrate forward from right face, if dir=1 from left (why reversed? to improve readability maybe correct)
    midCenterVector( dim ),
    midP( dim, dim ),
    midSection( midCenterVector, midCenterVector ),
    vectorFieldRev( _vectorFieldRev ),
    thetaRange( _thetaRange ),
    epsRange( _epsRange )
  {
    ICoordinateSection tempSection( dim, 0, ( (80./100.)*GammaCenter1[0] + (20./100.)*GammaCenter2[0] ) );  

    if( _segment1.isABlock ) // a check on whether we do a homoclinic orbit proof
    {
      tempSection.setConstant( (20./100.)*midVector(GammaCenter1)[0] + (80./100.)*midVector(GammaCenter2)[0] );
    }

    interval midThetaRange( thetaRange.leftBound()/2. + thetaRange.rightBound()/2. );
    interval midEpsRange( epsRange.leftBound()/2. + epsRange.rightBound()/2. );

    vectorField.setParameter("theta", midThetaRange ); // in what is in this block we need only approximately good variables
    vectorField.setParameter("eps", midEpsRange );

    ITaylor tempSolver( vectorField, order );
    IPoincareMap tempPM( tempSolver, tempSection );
    
    interval returnTime;
    IMatrix tempMonodromyMatrix( dim, dim );
    C1Rect2Set tempCenterSet( section1CenterVector );

    midCenterVector = tempPM( tempCenterSet, tempMonodromyMatrix, returnTime );
    midSection.setOrigin( midVector(midCenterVector) );

    IVector normalVector( leftVector( midVector(vectorField( midCenterVector )) ) );
    
    IEuclNorm vectorNorm;

    // here some changes! (for old version see with params construction). they dont work that great. improving frontcovering worsens backcovering & vice versa
   // IVector initialNormal = Transpose(inverseMatrix(P2)).column(0);
  //  IVector normalVector( leftVector( midVector( Transpose( inverseMatrix( tempMonodromyMatrix ) )*initialNormal ) ) ); // THAT IS NEW! (SEE ROBERTS NOTES!)

    midSection.setNormalVector( normalVector );

    // some objects here (solver, returntime) could have been reused but for safety we create new ones
    interval returnTime2;    
    interval returnTime2Rev;

    IMatrix monodromyMatrix( dim, dim );
    IMatrix monodromyMatrixRev( dim, dim );

    ITaylor tempSolver2( vectorField, order );
    ITaylor tempSolver2Rev( vectorFieldRev, order );

    C1Rect2Set C1TempCenterSet( section1CenterVector );
    C1Rect2Set C1TempCenterSetRev( section2CenterVector );

    IPoincareMap tempPM2( tempSolver2, midSection );
    IPoincareMap tempPM2Rev( tempSolver2Rev, midSection );

    IVector tempVect = tempPM2( C1TempCenterSet, monodromyMatrix, returnTime2 );

    IVector tempVectRev = tempPM2Rev( C1TempCenterSetRev, monodromyMatrixRev, returnTime2Rev );

    IVector stableDir = ( tempPM2Rev.computeDP( tempVectRev, monodromyMatrixRev, returnTime2Rev ) )*P2.column(2);
    IVector unstableDir = ( tempPM2.computeDP( tempVect, monodromyMatrix, returnTime2 ) )*P1.column(2);
                                                                                     // new coordinates are variational equations eval. at identity matrix times original matrix P1
                                                                                     // we cannot eval variational equations at P1 because CAPD doesn't support that yet 
                                                                                     // variational equation is linear though so that is ok

    stableDir = leftVector( midVector( stableDir / vectorNorm( stableDir ) ) );
    unstableDir = leftVector( midVector ( unstableDir / vectorNorm( unstableDir ) ) );

    for( int i = 1; i <= dim; i++ )
    {
      midP(i,2) = normalVector[i-1];      // we insert the section normal vector as the second column of the coordinate change matrix
                                          // as this was the unstable row of P1 in direction of which we integrated
                                          // so before the insertion this column was ~0. This should make the matrix nonsingular.
      midP(i,1) = stableDir(i);
      midP(i,3) = unstableDir(i);         // todo: for homoclinic it would be better to compute dpm/dtheta as the unstable dir
    }

    orthogonalizeRelativeColumn(midP,1);
    vectorField.setParameter("theta", thetaRange );
    vectorField.setParameter("eps", epsRange );
  }

  // #include "todo/midPMwithparams.hpp" 

  IVector integrateToMidSection( const IVector& theSet, bool dir ) // 2-dim h-set is embedded into space and integrated forward from section 1 if dir = 0 and backward from section 2 elsewise
  {
    interval totalTime;
    ITaylor *midSolver;
    IPoincareMap *midPM;

    if( !dir )
      midSolver = new ITaylor( vectorField, order );
    else
      midSolver = new ITaylor( vectorFieldRev, order );

    midPM = new IPoincareMap( *midSolver, midSection );

    IVector resultArr(2);
 
    int div_i( div );
    int div_j( div );

    // check whether we integrate one of the unstable edges of an h-set the stable edges of an h-set
    interval edgeIntegrationSign(0.); // stores information about which side of a left/right stable edge we integrate: <0 left >0 right ==0 we integrate the whole h-set
    if( theSet[0].leftBound() == theSet[0].rightBound() )        
    {
      div_i=1;
      edgeIntegrationSign = theSet[0].leftBound();
    }
    if( theSet[1].leftBound() == theSet[1].rightBound() )        
    {
      div_i=1;
      edgeIntegrationSign = theSet[1].leftBound();
    }


    for(int i=1; i<=div_i; i++)
    {
      interval ti = interval(i-1, i)/div_i;
      IVector Gamma_div( dim );

      if( edgeIntegrationSign > 0. && !dir )
        Gamma_div = GammaRight1;
      else if( edgeIntegrationSign < 0. && !dir )
        Gamma_div = GammaLeft1;
      else if( edgeIntegrationSign == 0. && !dir )
        Gamma_div = ( GammaRight1 - GammaLeft1 )*ti + GammaLeft1;     // subdivision of w coordinate
      else if( edgeIntegrationSign > 0. && dir )
        Gamma_div = GammaRight2;
      else if( edgeIntegrationSign < 0. && dir )
        Gamma_div = GammaLeft2;
      else if( edgeIntegrationSign == 0. && dir )
        Gamma_div = ( GammaRight2 - GammaLeft2 )*ti + GammaLeft2;     // subdivision of w coordinate
      else 
        throw "SUBDIVISION ERROR";

      for(int j=1; j<=div_j; j++)
      {
        interval tj = interval(j-1, j)/div_j;

        IVector Set_ij( dim ); // the centered part of the set with expanded directions
        Set_ij.clear();        

        if( !dir )
        {
          Set_ij[0] = ( theSet[0].rightBound() - theSet[0].leftBound() )*tj + theSet[0].leftBound();    // subdivision of ys coordinate
          Set_ij[1] = y1vector[1];
        }
        else
        {
          Set_ij[0] = y2vector[0];
          Set_ij[1] = ( theSet[1].rightBound() - theSet[1].leftBound() )*tj + theSet[1].leftBound();    // subdivision of yu coordinate
        }

        C0Rect2Set *setAff;

        if( !dir )
          setAff = new C0Rect2Set( Gamma_div, P1, Set_ij ); // the set moved to default space, observe that parameters remain unchanged when withParams is on
        else
          setAff = new C0Rect2Set( Gamma_div, P2, Set_ij );
        

        interval returntime(0.);
        IVector result = (*midPM)( *setAff, midCenterVector, inverseMatrix(midP), returntime );     // result is moved back to local coordinates, yu should be close to 0 
                                                                                            // in other words midP^-1( PM(setAff) - midCenterVector ) is computed 
                                                                                            // WARNING! FOR WITH_PARAMS=1 THIS ZEROES PARAMETERS SO AS SUCH RESULT SHOULD NOT BE USED,
                                                                                            // ONLY FIRST 3 COORDINATES OF IT (RETURNED BY THIS FUNCTION) CAN BE USED

        delete setAff;
 
        if( i==1 && j==1)
        {
          resultArr[0] = result[0];   // midSection coordinates are given by midP - matrix P1 evolved by var. equation so similarly to P1 we project to ys, w coords, w "unstable"
          resultArr[1] = result[2];  
          totalTime = returntime;
        }
        else
        {
          resultArr[0] = intervalHull(resultArr[0], result[0]);
          resultArr[1] = intervalHull(resultArr[1], result[2]);
          totalTime = intervalHull( totalTime, returntime );
        }
      }
    }
    delete midPM;
    delete midSolver; 

     //  cout << "TOTAL TIME: " << totalTime << "\n";
    return resultArr;
  }

  bool checkCovering( const IVector& Set1, const IVector& Set2, bool _verbose )  // both Set1 and Set2 are 2-dim and have first variable stable second unstable ( Set1 : ys, w; Set2 : w, yu )
  {
    IVector PSet1( integrateToMidSection( Set1, 0 ) );
    IVector PSetUL1( integrateToMidSection( leftU(Set1), 0 ) );
    IVector PSetUR1( integrateToMidSection( rightU(Set1), 0 ) );

    IVector PSet2( integrateToMidSection( Set2 , 1 ) );
    IVector PSetSL2( integrateToMidSection( leftS(Set2), 1 ) );
    IVector PSetSR2( integrateToMidSection( rightS(Set2), 1 ) );


    IVector setToBackCover( shrinkAndExpand( PSet1, 1. + EPS ) );       // we define a set that is covered by midPM( Set1 ), shrinkAndExpand adjust stable direction
    setToBackCover[1] = interval( (PSetUL1[1] + EPS).rightBound(), (PSetUR1[1] - EPS).leftBound() );    // here we adjust the unstable direction
   
    if( _verbose )
    {
      cout << "Right bound of image of left stable edge for the backcovering set: " <<  PSetSL2[0].rightBound() 
        << "\nLeft bound of the stable direction for the set to be covered: " << setToBackCover[0].leftBound();
      cout << "\n --- \n";
      cout << "Left bound of image of right stable edge for the backcovering set: " <<  PSetSR2[0].leftBound() 
        << "\nRight bound of the stable direction for the set to be covered: " << setToBackCover[0].rightBound();
      cout << "\n --- \n";
      cout << "Bound of the unstable direction projection of the image of the backcovering set: " << PSet2[1] 
        << "\nUnstable direction of the set to be backcovered: " << setToBackCover[1] << "\n";
    }


    if( !( PSetUL1[1] + EPS < 0. && PSetUR1[1] - EPS > 0. && PSetSL2[0] + EPS < 0. && PSetSR2[0] - EPS > 0. ) )  // some reality checks for hyperbolicity
     throw "INTEGRATION TO MIDSECTION ERROR! \n";

    if( PSetSL2[0] < setToBackCover[0].leftBound() && PSetSR2[0] > setToBackCover[0].rightBound() && subsetInterior( PSet2[1], setToBackCover[1] ) )
      return 1;
    else 
      return 0;
  }
 
  bool shootWithTheta( C0Rect2Set uMan, C0Rect2Set uManLeft, C0Rect2Set uManRight, const IVector& Set2, bool _verbose=1 )  // both Set1 and Set2 are 2-dim and have first variable stable second unstable ( Set1 : ys, v; Set2 : v, yu )
  {
    // shooting with theta - theta is the unstable direction of a 1-dim h-set, rest is "error" in enclosure of the unstable manifold of (0,0,0)

    IPoincareMap midPM( solver, midSection );

    interval returnTime(0.);
    vectorField.setParameter("theta", thetaRange);
    IVector PSet1( midPM( uMan, midCenterVector, inverseMatrix(midP), returnTime ) );
      
    returnTime = 0.;
    vectorField.setParameter("theta", thetaRange.leftBound() );
    IVector PSetUL1( midPM( uManLeft, midCenterVector, inverseMatrix(midP), returnTime ) );
 
    returnTime = 0.;
    vectorField.setParameter("theta", thetaRange.rightBound() );
    IVector PSetUR1( midPM( uManRight, midCenterVector, inverseMatrix(midP), returnTime ) );

    // these two lines to come back to original parameters
    vectorField.setParameter("theta", thetaRange );
    vectorFieldRev.setParameter("theta", thetaRange );

    // here we do a backcovering
    IVector PSet2( integrateToMidSection( Set2 , 1 ) );
    IVector PSetSL2( integrateToMidSection( leftS(Set2), 1 ) );
    IVector PSetSR2( integrateToMidSection( rightS(Set2), 1 ) );

  //  cout << "PSet1 = " << PSet1 << "\n " << "PSetUL1 = " << PSetUL1 << "\n PSetUR1 = " << PSetUR1 << "\n \n";
  //  cout << "PSet2 = " << PSet2 << "\n " << "PSetSL2 = " << PSetSL2 << "\n PSetSR2 = " << PSetSR2 << "\n \n";

    IVector setToBackCover( shrinkAndExpand( PSet1, 1. + EPS ) );       // we define a set that is covered by midPM( Set1 ), shrinkAndExpand adjust stable direction
    setToBackCover[1] = interval( (PSetUL1[2] + EPS).rightBound(), (PSetUR1[2] - EPS).leftBound() );    // here we adjust the unstable direction -- we remember that the third coordinate in
                                                                                                        // midP variables is the "unstable" one
   
    if( _verbose )
    {
      cout << "Right bound of image of left stable edge for the backcovering set: " <<  PSetSL2[0].rightBound() 
        << "\nLeft bound of the stable direction for the set to be covered: " << setToBackCover[0].leftBound();
      cout << "\n --- \n";
      cout << "Left bound of image of right stable edge for the backcovering set: " <<  PSetSR2[0].leftBound() 
        << "\nRight bound of the stable direction for the set to be covered: " << setToBackCover[0].rightBound();
      cout << "\n --- \n";
      cout << "Bound of the unstable direction projection of the image of the backcovering set: " << PSet2[1] 
        << "\nUnstable direction of the set to be backcovered: " << setToBackCover[1] << "\n";
    }


    if( !( PSetUL1[2] + EPS < 0. && PSetUR1[2] - EPS > 0. && PSetSL2[0] + EPS < 0. && PSetSR2[0] - EPS > 0. ) )  // some reality checks for hyperbolicity
     throw "INTEGRATION TO MIDSECTION ERROR! \n";

    if( PSetSL2[0] < setToBackCover[0].leftBound() && PSetSR2[0] > setToBackCover[0].rightBound() && subsetInterior( PSet2[1], setToBackCover[1] ) )
      return 1;
    else 
      return 0;
  }

};



