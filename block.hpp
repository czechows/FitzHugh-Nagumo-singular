//  first coordinate stable second unstable third central

class FhnBlockWithCones 
{
  public:
    IMap vectorField;         // has to be the FHN vector field and has to have a parameter eps
    interval a;               // parameter for FHN
    interval deltaU;          // size in unstable dir
    interval deltaS;          // size in stable dir
    interval deltaMu;         // size in central (second stable) dir
    IMatrix CB;               // change of coordinates to diagonal basis
    IMatrix InvCB;            // inverse of the above
    IMatrix deltaM;            // matrix of deltas

    FhnBlockWithCones( IMap _vectorField, const interval& _deltaU, const interval& _deltaS, const interval& _deltaMu, const interval& _a = interval(1./10.) )
      : vectorField( _vectorField ),
        a( _a ),
        deltaU( _deltaU ),
        deltaS( _deltaS ),
        deltaMu( _deltaMu ),
        CB(3,3),
        InvCB(3,3),
        deltaM(3,3)
  {
    interval deltaMdat[] = { deltaS, 0., 0., 0., deltaU, 0., 0., 0., deltaMu };
    deltaM = IMatrix(3,3,deltaMdat);

    InvCB = coordChange( vectorField, IVector(0.,0.,0.) ); // equilibrium is at 0


    InvCB[0][2] = -1./a;
    InvCB[1][2] = 0.;
    InvCB[2][2] = 1;              // the last column is the tangent to the slow manifold, this should be automated

    InvCB = leftVector( midVector( InvCB ) );

    InvCB[2][0] = 0.;             // this is already done in coordChange, but just for clarity
    InvCB[2][1] = 0.; 

    InvCB = InvCB*deltaM;

    CB = krawczykInverse( InvCB );
    CB = leftVector( midVector( CB ) );

    InvCB = krawczykInverse( CB );
 
    InvCB[2][0] = 0.;             // this is already done in coordChange, but just for clarity
    InvCB[2][1] = 0.; 
    
    CB[2][0] = 0.;             
    CB[2][1] = 0.; 

    //if( vectalg::containsZero( IVector({det(CB)}) ) || vectalg::containsZero( IVector({det(InvCB)}) ) )
    //  throw "Change of coordinates not invertible!";
  }

  IVector evaluateVFinNewVariables( IVector x )
  {
    C0Rect2Set ev( IVector(0.,0.,0.), InvCB, x );
    DiscreteDynSys<IMap> vectorFieldEval( vectorField );
    ev.move( vectorFieldEval );
    IVector result( ev );

    result = CB*result; // the vector field in block variables
    return result; 
  }

  IVector evaluateLastRowDFcDivEps( IVector x ) // computes 1/eps (last row DF_c), that is the slow part of the vector field without epsilon -- only when the slow part is epsilon-independent!
  {
    IMap newVectorField( vectorField );
    newVectorField.setParameter("eps", interval(1.));      // from the form of C we have last row of Qeps * CB * DF * CB-1
                                                           // is equal to last row of Q1 * CB * DFdivEps * CB-1
                                                           // and only the last row of DFdivEps matters
    return ( CB*newVectorField[ InvCB*x ]*InvCB ).row(2);
  }

 
  IMatrix QepsDFc( IVector x ) // evaluates QepsDFc, Qeps = diag(1,-1,-1/eps)
  {
    interval Q1formdat[] = {-1.,0.,0.,0.,1.,0.,0.,0.,-1.};
    IMatrix Q1(3,3,Q1formdat);
    
    IMatrix result = Q1*CB*vectorField[ InvCB*x ]*InvCB;
    IVector resultLastRow = -evaluateLastRowDFcDivEps( x ); 

    result[2][0] = resultLastRow[0];
    result[2][1] = resultLastRow[1];
    result[2][2] = resultLastRow[2];

    return result; 
  }

  void coneConditionsVerification() // if it doesn't throw the cone conditions are verified
  {
    IVector Bc( interval(-1.,1.), interval(-1.,1.), interval(-1.,1.) );
    
    IMatrix SymDF( Transpose( QepsDFc( Bc ) ) + QepsDFc( Bc ) );

    if( !( SymDF(1,1) > 0. ) )
    {
      cout << "det(firstPrincipalMinor)=" << SymDF(1,1) << "\n";
      cout.flush(); 
      throw "Error in verification of cone conditions\n";
    }

    if( !( SymDF(1,1)*SymDF(2,2) - SymDF(2,1)*SymDF(1,2) > 0. ) )
    {
      cout << "det(secondPrincipalMinor)=" << SymDF(1,1)*SymDF(2,2) - SymDF(2,1)*SymDF(1,2) << "\n";
      cout.flush(); 
      throw "Error in verification of cone conditions\n";
    }


    if( !( det( SymDF ) > 0. ) )
    {
      cout << "det(SymDF)=" << det( SymDF ) << "\n";
      cout << "SymDF= " << SymDF << "\n";
      cout.flush(); 
      throw "Error in verification of cone conditions\n";
    }
  }

  IVector getEnclosureUnstableMan()
  {
    return InvCB*IVector(interval(-1.,1.), 1., interval(-1.,1.));
  }
 
  IMatrix getPstableMan()
  {
    IMatrix result( InvCB );
    result[2][0] = 0.;
    result[2][1] = 0.;
    result[0][2] = 0.;
    result[1][2] = 0.;
    result[2][2] = 1.;
    
    interval mnormDat[] = { 1./deltaS, 0., 0., 0., 1./deltaU, 0., 0., 0., 1. };
    IMatrix mnorm(3,3,mnormDat);

    result = result*mnorm; // we want to creat a matrix close to identity

    return result;
  }

  IVector getGammaRightStableMan()
  {
    return IVector( InvCB*IVector(0.,0.,1.) );
  }

  IVector getFaceStableMan()
  {
    IVector result(3);
    result[0] = deltaS*interval(-1.,1.);         // corrections due to normalization in getPstableMan
    result[1] = deltaU*interval(-1.,1.);
    result[2] = 0.;
    return result;
  }

  C0Rect2Set getUnstableManBound()
  {
    return C0Rect2Set( IVector(0.,0.,0.), InvCB, IVector(interval(-1.,1.), 1., interval(-1.,1.)) );
  }

  FhnIsolatingBlock createABlock( int divCount ) // creates a block for isolation verification
  {
    IVector Left( InvCB*IVector(0.,0.,-1.) );
    IVector Right( InvCB*IVector(0.,0.,1.) );
    IMatrix P = getPstableMan(); 

    return FhnIsolatingBlock( vectorField, Left, Right, P, getFaceStableMan(), getFaceStableMan(), divCount );
  }
};


class uManBlockWithCones : public FhnBlockWithCones  // a class that allows to propagate the unstable manifold, i.e. no need to verify the cone conditions on the whole set
{
public:
  FhnBlockWithCones shortBlock;

  uManBlockWithCones( IMap _vectorField, const interval& _deltaU, const interval& _deltaS, const interval& _deltaMu, const interval uProportion = interval(0.3), 
      const interval& _a = interval(1./10.) )
  : FhnBlockWithCones( _vectorField, _deltaU, _deltaS, _deltaMu, _a ), 
    shortBlock(  _vectorField, uProportion*_deltaU, _deltaS, _deltaMu, _a ) // it is a subblock of the uManBlock

  {
    shortBlock.coneConditionsVerification(); // we verify cone conditions on a shorter block

    IVector upLeftOverBlockSupport( interval(-1.,1.), interval( uProportion.leftBound(), 1. ), interval(-1.,1.) );
    if( !( evaluateVFinNewVariables( upLeftOverBlockSupport )[1] > 0. ) )
    {
      cout << "evaluateVFinNewVariables( upLeftOverBlockSupport )[1]= " << evaluateVFinNewVariables( upLeftOverBlockSupport )[1] << "\n";
      cout.flush();
      throw "Non-uniform flow in the unstable direction for the propagation of the unstable manifold \n";
    };
  }
};
