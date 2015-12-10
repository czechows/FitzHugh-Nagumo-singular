/* -----------------------------------------------------------------------------------------
 * This is a header file to fhn.cpp providing auxiliaries for rigorous verification of
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

bool isCovering( const IVector& setCovering, const IMatrix& setCoveringCoord, const IVector& setToCover ) 
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



IVector shrinkAndExpand(const IVector &N, interval factor)  // shrinks a rectangle in unstable direction and expands it in stable to get a covering (for example by original rectangle)
{
    IVector result(N);
    result[0] = N[0]*factor;
    result[1] = N[1]/factor;
    return result;
}

void orthogonalizeRelativeColumn( IMatrix& matrixToOrthogonalize, unsigned int columnNo )
{
  for( unsigned int i = 0; i <= matrixToOrthogonalize.numberOfColumns() - 1; i++ ) 
  { 
    IVector vectorInvariant( matrixToOrthogonalize.column( columnNo ) );
    if( i != columnNo )
    {
      IVector vectorToOrthogonalize( matrixToOrthogonalize.column(i) );
      vectorToOrthogonalize = leftVector( midVector( vectorToOrthogonalize ) );
      IVector projection = ( scalarProduct( vectorToOrthogonalize, vectorInvariant )/scalarProduct( vectorInvariant, vectorInvariant ) ) * vectorInvariant;

      for( unsigned int j = 1; j <= matrixToOrthogonalize.numberOfRows(); j++ )
      {
        matrixToOrthogonalize(j,i+1) = vectorToOrthogonalize(j) - projection(j);
      }
    }
  }
}


IVector Eq_correct(IMap vectorField, const IVector guess) 
// corrects initial guesses of u so they are closer to real equilibria of the slow flow; using Newton alg., w is always 0, computes on  intervals to avoid code repetition - but is not rigorous (! this is not an interval Newton operator !) -- ONLY FOR FITZHUGH-NAGUMO VECTOR FIELD
{
  interval error(1.);
  interval result(guess[0]);
  interval oldresult(result);
  while(error > accuracy)
  {
    oldresult = result;
    result = result - (vectorField( IVector(result, 0., guess[2]) )[1])/(vectorField[ IVector(result, 0., guess[2]) ][1][0]);  
                                                                        // Newton algorithm to calculate zeroes of the vector field - v is always 0.,
                                                                        // derivative is of the second equation with respect to first variable - u
    error = abs(oldresult - result);
  }
  IVector new_Eq(3);
  new_Eq[0] = result;
  new_Eq[1] = interval(0.);
  new_Eq[2] = guess[2];
  return new_Eq;
}

