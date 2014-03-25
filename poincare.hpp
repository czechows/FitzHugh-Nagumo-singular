
/* -----------------------------------------------------------------------------------------
 * This is a header file to fhn.cpp providing Poincare map classes for the Poincare maps
 * between branches of the slow manifold and auxiliaries for rigorous verification of
 * 2-dim covering relations.
 * ----------------------------------------------------------------------------------------*/




/* ----------------------------------------------------------------------------------------- */
/* ---------------------------- COVERING RELATIONS ----------------------------------------- */
/* ----------------------------------------------------------------------------------------- */

// left right stable unstable edges for covering - use with 2dim vectors with first variable stable second unstable, can have more dimensions of neutral variables

IVector leftU(const IVector &N)
{
  IVector _leftU( N.dimension() ); 
  _leftU = N;
  _leftU[1] =  N[1].leftBound();
  return _leftU;
}

IVector rightU(const IVector &N)
{
  IVector _rightU( N.dimension() );
  _rightU = N;
  _rightU[1] =  N[1].rightBound();
  return _rightU;
}

IVector leftS(const IVector &N)
{
  IVector _leftS( N.dimension() ); 
  _leftS = N;
  _leftS[0] =  N[0].leftBound();
  return _leftS;
}

IVector rightS(const IVector &N)
{
  IVector _rightS( N.dimension() );
  _rightS = N;
  _rightS[0] =  N[0].rightBound();
  return _rightS;
}

bool isCovering( IVector& setCovering, const IMatrix& setCoveringCoord, IVector& setToCover ) 
                                                      // verifies covering between image of setCovering by a matrix setCoveringCoord over setToCover 
                                                      // first variable stable second unstable
{
 bool leftCheck( (setCoveringCoord*( leftU(setCovering) ))[1] < setToCover[1].leftBound() );
 bool rightCheck( (setCoveringCoord*( rightU(setCovering) ))[1] > setToCover[1].rightBound() );

 if( leftCheck && rightCheck && subsetInterior( (setCoveringCoord*setCovering)[0], setToCover[0] ) ) 
  return 1;
 else
   return 0;
};

/*
bool isBackwardCovering( IVector& setCovering, const IMatrix& setCoveringCoord, IVector& setToCover ) 
                                                      // verifies backward covering between image of setCovering by a matrix setCoveringCoord over setToCover 
                                                      // first variable stable second unstable
{
 bool leftCheck( (setCoveringCoord*( leftS(setCovering) ))[0] < setToCover[0].leftBound() );
 bool rightCheck( (setCoveringCoord*( rightS(setCovering) ))[0] > setToCover[0].rightBound() );

 if( leftCheck && rightCheck && subsetInterior( (setCoveringCoord*setCovering)[1], setToCover[1] ) ) 
  return 1;
 else
   return 0;
};
*/

IVector shrinkAndExpand(IVector &N, interval factor)  // shrinks a rectangle in unstable direction and expands it in stable to get a covering (for example by original rectangle)
{
    IVector result(N);
    result[0] = N[0]*factor;
    result[1] = N[1]/factor;
    return result;
}




/* ------------------------------------------------------------------------------------ */
/* ---------------------------- POINCARE MAPS ----------------------------------------- */
/* ------------------------------------------------------------------------------------ */

class FhnPoincareMap
{
public:
  int dim;          // phase space dimension
  IMap vectorField;
  ITaylor solver;
  IVector y2vector;
  IVector section2CenterVector;
  IAffineSection section2;
  IVector y1vector;
  IVector section1CenterVector;
  IPoincareMap pm;
  IMatrix P1;
  IMatrix P2;
  int disc;         // number of subdivisions in each dimension (supports two, does not support subdivisions in parameter space) for integration of h-sets
  IVector params;   // vector of parameters
  IVector GammaU1;
  IVector GammaU2;

  FhnPoincareMap( IMap _vectorField, const IMatrix& _P1, const IMatrix& _P2, const IVector& _GammaU1, const IVector& _GammaU2, 
                                                                                        interval& _ru1, interval& _rs2, interval dir = interval(1.), int _disc=1 ) 
    : dim( 3 ),
      vectorField( _vectorField ),                               // a 3d vector field
      solver( vectorField, order ),
      y2vector( dir*_rs2 , 0., 0. ),                             // a section center given in linearization around slow manifold coordinates (normal is (-rs2,0,0))
      section2CenterVector( _P2*y2vector + _GammaU2 ),
      section2( section2CenterVector, Transpose(inverseMatrix(_P2))*y2vector ),      // a section moved to uwv space normal is _P2^(-T)*(-y2vector)
  
      y1vector( 0., -dir*_ru1, 0. ),                              // same for y1vector, the center of the set to integrate
      section1CenterVector( _P1*y1vector + _GammaU1 ),

      pm( solver, section2 ),
      P1( _P1 ),
      P2( _P2 ),
      disc( _disc ),
      params( 1 ),
      GammaU1( _GammaU1 ),
      GammaU2( _GammaU2 )
  {
  }

  FhnPoincareMap( IVector _params, IMap _vectorFieldWithParams, const IMatrix& _P1, const IMatrix& _P2, const IVector& _GammaU1, const IVector& _GammaU2, 
                                                                                        interval& _ru1, interval& _rs2, interval dir = interval(1.), int _disc=1 ) 
    : dim( 3 + _params.dimension() ),
      vectorField( _vectorFieldWithParams ),                               
      solver( vectorField, order ),
      y2vector( dim ),                             
      section2CenterVector( dim ),
      section2( IVector(dim), IVector(dim) ),      // uninitialized section
      y1vector( dim ),                              
      section1CenterVector( dim ),
      pm ( solver, section2 ),
      P1( dim, dim ),
      P2( dim, dim ),
      disc( _disc ),
      params( _params ),
      GammaU1( dim ),
      GammaU2( dim )
  {
    y1vector = IVector( dim );
    y1vector.clear();                 // ensures vector is all zeroes
    y1vector[1] = -dir*_ru1;

    y2vector = IVector( dim );
    y2vector.clear();
    y2vector[0] = dir*_rs2;

    for( int i=0; i < dim; i++ )
    {
      if( i < 3 )
      {
        GammaU1[i] = _GammaU1[i];
        GammaU2[i] = _GammaU2[i];
      }
      else
      {
        GammaU1[i] = params[i-3];     // we embed the parameters
        GammaU2[i] = params[i-3]; 
      }

      for( int j=0; j < dim; j++ )
      {
        if( j < 3 && i < 3 )
        {
          P1(i+1,j+1) = _P1(i+1,j+1);               // indexing of matrices is shifted by 1 (starts from 1, not 0)
          P2(i+1,j+1) = _P2(i+1,j+1);
        }
        else if( i == j )
          P1(i+1,j+1) = P2(i+1,j+1) = 1.;
        else
          P1(i+1,j+1) = P2(i+1,j+1) = 0.;
      }
    }

    section1CenterVector = IVector( P1*y1vector + GammaU1 );
    section2CenterVector = IVector( P2*y2vector + GammaU2 );

    section2.setOrigin( section2CenterVector );
    section2.setNormalVector( Transpose(inverseMatrix(P2)).column(1) );         // solver works by reference so gets autoupdated
                                                                               // normal vector is one of the rows of inverseMatrix(P2), that matches
                                                                               // with transforming the result of Poincare map by inverseMatrix(P2)
  }


  IVector operator()(const IVector& theSet) // we give a set in local variables on one section centered on 0 (ys & v_centered) return in variables on the other (v_centered & yu)
  {
    IVector resultArr(2);

    IMatrix P2Alt( P2 );
    for( int k = 1; k <= dim; k++ )
      P2Alt(k,1) = Transpose(inverseMatrix(P2))(k,1);  // we need to set the altered normal vector to the change of coordinates matrix DEPRECATED???

    for(int i=1; i<=disc; i++)
    {
      int disc1;
      if( theSet[0].leftBound() == theSet[0].rightBound() )   // check whether we integrate one of the unstable edges of an h-set
        disc1=1;
      else 
        disc1=disc;
 
      interval ti = interval(i-1, i)/disc;

      for(int j=1; j<=disc1; j++)
      {
        interval tj = interval(j-1, j)/disc1;

        IVector Set_ij( dim ); // the centered part of the set with expanded directions
        Set_ij.clear();        

        Set_ij[0] = ( theSet[0].rightBound() - theSet[0].leftBound() )*ti + theSet[0].leftBound();    // subdivision of ys coordinate
        Set_ij[2] = ( theSet[1].rightBound() - theSet[1].leftBound() )*tj + theSet[1].leftBound();    // subdivision of v coordinate

  /*      if( dim > 3 )  // checks whether we have parameters
          for( int k=3; k < dim; k++ )
            Set_ij[k] = 0.; // params[k-3];      // we embed parameters  */ // this part of code is probably deprecated since we embedded the parameters in the constructor

        C0HOTripletonSet setAff( section1CenterVector, P1, Set_ij ); // the set moved to default space, observe that parameters remain unchanged

        interval returntime(0.);
        IVector result = pm( setAff, GammaU2, inverseMatrix(P2), returntime ); // result is moved back to local coordinates, ys should be close to 0 
                                                                                            // in other words P2Alt^-1( PM(setAff) - section2center ) is computed 
                                                                                            // where section2center = _P2&y2vector + _GammaU2
                                                                                            //  IMPORTANT: I NEED TO CHANGE HERE TO GAMMA2?
        if( i==1 && j==1)
        {
          resultArr[0] = result[2];   // being on a section given by ys we only return v,yu coordinates, now v is the stable
          resultArr[1] = result[1];  
        }
        else
        {
          resultArr[0] = intervalHull(resultArr[0], result[2]);
          resultArr[1] = intervalHull(resultArr[1], result[1]);
        }
      }
    }
    return resultArr; 
  }
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
  
  midPoincareMap( IMap _vectorField, IMap _vectorFieldRev, const IMatrix& _P1, const IMatrix& _P2, const IVector& _GammaU1, const IVector& _GammaU2, 
                                                                                        interval& _ru1, interval& _rs2, interval dir = interval(1.), int _disc=1 ) 
  : FhnPoincareMap( _vectorField, _P1, _P2, _GammaU1, _GammaU2, _ru1, _rs2, dir, _disc ),
    midCenterVector( dim ),
    midP( dim, dim ),
    midSection( midCenterVector, midCenterVector ),
    vectorFieldRev( _vectorFieldRev )
  {
    ICoordinateSection tempSection( dim, 0, ( (93./100.)*_GammaU1[0] + (7./100.)*_GammaU2[0] ) ); // an auxiliary section u = ( GammaU1[0] + GammaU2[0] )/2
    ITaylor tempSolver( vectorField, order );
    IPoincareMap tempPM( tempSolver, tempSection );
    interval returnTime;
    C0HOTripletonSet C0TempCenterSet( section1CenterVector );

    midCenterVector = tempPM( C0TempCenterSet, returnTime );
    midSection.setOrigin( midVector(midCenterVector) );
    midSection.setNormalVector( midVector(vectorField( midCenterVector )) );

    interval returnTime2;                                       // some objects here (solver, returntime) could have been reused but for safety we create new ones
    IMatrix monodromyMatrix( dim, dim );
    ITaylor tempSolver2( vectorField, order );
    C1Rect2Set C1TempCenterSet( section1CenterVector );

    IPoincareMap tempPM2( tempSolver2, midSection );
    IVector tempVect = tempPM2( C1TempCenterSet, monodromyMatrix, returnTime2 );

    midP = ( tempPM2.computeDP( tempVect, monodromyMatrix, returnTime2 ) )*P1;        // new coordinates are variational equations eval. at identity matrix times original matrix P1
                                                                                     // we cannot eval variational equations at P1 because CAPD doesn't support that yet 
                                                                                     // variational equation is linear though so that is ok
 
    for( int i = 1; i <= dim; i++ )
      midP(i,2) = midVector( vectorField( midCenterVector ) )[i-1];  // we insert the section normal vector as the second column of the coordinate change matrix
                                                                     // as this was the unstable row of P1 in direction of which we integrated
                                                                     // so before the insertion this column was ~0. This should make the matrix nonsingular.
  }
  
  midPoincareMap( IVector _params, IMap _vectorField, IMap _vectorFieldRev, const IMatrix& _P1, const IMatrix& _P2, const IVector& _GammaU1, const IVector& _GammaU2, 
                                                                                        interval& _ru1, interval& _rs2, interval dir = interval(1.), int _disc=1 ) 
    // the same but with params treated as variables of velocity 0, same as with second constructor of FhnPoincareMap
  : FhnPoincareMap( _params, _vectorField, _P1, _P2, _GammaU1, _GammaU2, _ru1, _rs2, dir, _disc ),  
    midCenterVector( dim ),
    midP( dim, dim ),
    midSection( midCenterVector, midCenterVector ),
    vectorFieldRev( _vectorFieldRev )
  {
    ICoordinateSection tempSection( dim, 0, ( (945./1000.)*_GammaU1[0] + (55./1000.)*_GammaU2[0] ) ); // an auxiliary section u = ( GammaU1[0] + GammaU2[0] )/2
    ITaylor tempSolver( vectorField, order );
    IPoincareMap tempPM( tempSolver, tempSection );
    interval returnTime;
    C0HOTripletonSet C0TempCenterSet( section1CenterVector );

    midCenterVector = tempPM( C0TempCenterSet, returnTime );
    midSection.setOrigin( midVector(midCenterVector) );
    midSection.setNormalVector( midVector( vectorField( midCenterVector ) ) );

    interval returnTime2;
    IMatrix monodromyMatrix( dim, dim );
    ITaylor tempSolver2( vectorField, order );                  
    C1Rect2Set C1TempCenterSet( section1CenterVector );

    IPoincareMap tempPM2( tempSolver2, midSection );
    IVector tempVect = tempPM2( C1TempCenterSet, monodromyMatrix, returnTime2 );

    midP = ( tempPM2.computeDP( tempVect, monodromyMatrix, returnTime2 ) )*P1;  

    for( int i = 1; i <= dim; i++ )
    {
      if( i > 3 )
      {
        for( int j = 1; j <= dim; j++ )         // for the parameters we want the coordinate change matrix to be like identity
        {                                       // so we adjust it manually (observe that this agrees with that midSection normal vector is 0 on coordinates corresp. to parameters
          midP(i,j) = 0.;                       // We are allowed to perform such actions because the coordinate change does not need to be rigorous, just good enough to get coverings.
          midP(j,i) = 0.;                       // The only thing to keep in mind is that the section normal vector has to agree with one of the coordinate change columns.  
        }
        midP(i,i) = interval(1.);
      }
      midP(i,2) = midVector( vectorField( midCenterVector ) )[i-1];
    }
    
  }


  IVector integrateToMidSection( const IVector& theSet, bool dir ) // 2-dim h-set is embedded into space and integrated forward from section 1 if dir = 0 and backward from section 2 elsewise
  {
    ITaylor *midSolver;
    IPoincareMap *midPM;

    if( !dir )
      midSolver = new ITaylor( vectorField, order );
    else
      midSolver = new ITaylor( vectorFieldRev, order );

    midPM = new IPoincareMap( *midSolver, midSection );

    IVector resultArr(2);
 
    int disc_i;
    int disc_j;

    if( theSet[0].leftBound() == theSet[0].rightBound() )   // check whether we integrate one of the unstable edges of an h-set
       disc_j=1;
    else 
       disc_j=disc;
  
    if( theSet[1].leftBound() == theSet[1].rightBound() )   // check whether we integrate one of the stable edges of an h-set
       disc_i=1;
    else 
       disc_i=disc;

    for(int i=1; i<=disc_i; i++)
    {

      interval ti = interval(i-1, i)/disc_i;

      for(int j=1; j<=disc_j; j++)
      {
        interval tj = interval(j-1, j)/disc_j;

        IVector Set_ij( dim ); // the centered part of the set with expanded directions
        Set_ij.clear();        

        if( !dir )
        {
          Set_ij[0] = ( theSet[0].rightBound() - theSet[0].leftBound() )*ti + theSet[0].leftBound();    // subdivision of ys coordinate
          Set_ij[2] = ( theSet[1].rightBound() - theSet[1].leftBound() )*tj + theSet[1].leftBound();    // subdivision of v coordinate
        }
        else
        {
          Set_ij[2] = ( theSet[0].rightBound() - theSet[0].leftBound() )*ti + theSet[0].leftBound();    // subdivision of v coordinate
          Set_ij[1] = ( theSet[1].rightBound() - theSet[1].leftBound() )*tj + theSet[1].leftBound();    // subdivision of yu coordinate
        }
        

 /*       if( dim > 3 )  // checks whether we have parameters
          for( int k=3; k < dim; k++ )
            Set_ij[k] = params[k-3];      // we embed parameters  DEPRECATED
*/
        C0HOTripletonSet *setAff;

        if( !dir )
          setAff = new C0HOTripletonSet( section1CenterVector, P1, Set_ij ); // the set moved to default space, observe that parameters remain unchanged
        else
          setAff = new C0HOTripletonSet( section2CenterVector, P2, Set_ij );
        
  
        interval returntime(0.);
        IVector result = (*midPM)( *setAff, midCenterVector, inverseMatrix(midP), returntime );     // result is moved back to local coordinates, yu should be close to 0 
                                                                                            // in other words midP^-1( PM(setAff) - midCenterVector ) is computed 
                                                                                            // WARNING! THIS ZEROES PARAMETERS SO AS SUCH RESULT SHOULD NOT BE USED,
                                                                                            // ONLY FIRST 3 COORDINATES OF IT (RETURNED BY THIS FUNCTION) CAN BE USED
        

        delete setAff;
 
        if( i==1 && j==1)
        {
          resultArr[0] = result[0];   // midSection coordinates are given by midP - matrix P1 evolved by var. equation so similarly to P1 we project to ys, v coords, v "unstable"
          resultArr[1] = result[2];  
        }
        else
        {
          resultArr[0] = intervalHull(resultArr[0], result[0]);
          resultArr[1] = intervalHull(resultArr[1], result[2]);
        }
      }
    }
    delete midPM;
    delete midSolver; 
    return resultArr;
  }


  bool checkCovering( const IVector& Set1, const IVector& Set2 )  // both Set1 and Set2 are 2-dim and have first variable stable second unstable ( Set1 : ys, v; Set2 : v, yu )
  {
    IVector PSet1( integrateToMidSection( Set1, 0 ) );
    IVector PSetUL1( integrateToMidSection( leftU(Set1), 0 ) );
    IVector PSetUR1( integrateToMidSection( rightU(Set1), 0 ) );

    IVector PSet2( integrateToMidSection( Set2 , 1 ) );
    IVector PSetSL2( integrateToMidSection( leftS(Set2), 1 ) );
    IVector PSetSR2( integrateToMidSection( rightS(Set2), 1 ) );


    IVector setToBackCover( shrinkAndExpand( PSet1, 1. + EPS ) );       // we define a set that is covered by midPM( Set1 ), shrinkAndExpand adjust stable direction
    setToBackCover[1] = interval( (PSetUL1[1] + EPS).rightBound(), (PSetUR1[1] - EPS).leftBound() );    // here we adjust the unstable direction
   
    cout << "PSetSL2[0] <? setToBackCover[0].leftBound() " <<  PSetSL2[0].rightBound() << " " << setToBackCover[0].leftBound() << "\n";
    cout << "PSetSR2[0] >? setToBackCover[0].rightBound() " <<  PSetSR2[0].leftBound() << " " << setToBackCover[0].rightBound() << "\n";
    cout << "subsetInterior( PSet2[1], setToBackCover[1] )? " << PSet2[1] << " " << setToBackCover[1] << "\n";


    if( !( PSetUL1[1] + EPS < 0. && PSetUR1[1] - EPS > 0. && PSetSL2[0] + EPS < 0. && PSetSR2[0] - EPS > 0. ) )  // some reality checks for hyperbolicity
     throw "INTEGRATION TO MIDSECTION ERROR 1! \n";
   
 
    if( !( containsZero( PSet1 ) && containsZero( PSet2 ) && containsZero( setToBackCover ) ) )     // some of these "reality checks" may be redundant but better to have them
      throw "INTEGRATION TO MIDSECTION ERROR 2! \n";                                                  // to ensure all alignments are correct

    if( PSetSL2[0] < setToBackCover[0].leftBound() && PSetSR2[0] > setToBackCover[0].rightBound() && subsetInterior( PSet2[1], setToBackCover[1] ) )
      return 1;
    else 
      return 0;
  }

};



