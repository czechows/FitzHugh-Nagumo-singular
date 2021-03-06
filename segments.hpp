
/* -----------------------------------------------------------------------------------------
 * This is a header file to fhn.cpp providing classes for isolating segments: FhnIsolatingSegment,
 * which is used when a constructed segment is short and one can stick with one affine coordinate change
 * to check the isolation and longIsolatingSegment - a derived class
 * to construct a chain of (possibly rotating) isolating segments and check covering relations 
 * between their faces. We also provide a coordinate change function to straightened fast coordinates
 * along the slow manifold.
 * ----------------------------------------------------------------------------------------*/



/* ------------------------------------------------------------------------------------ */
/* ---------------------------- COORDINATE CHANGE ------------------------------------- */
/* ------------------------------------------------------------------------------------ */

IMatrix coordChange( IMap vectorField, const IVector& Gamma ) // matrix of coordinate change on slow manifold 
                                                              // from straightened stable/unstable (i.e. (1,0,0), (0,1,0)) coordinates to real ones; third "neutral" variable unchanged
                                                              // doesn't need to be rigorous (and isn't) 
{
  int vdim( Gamma.dimension() );
  DMatrix JacobianD( vdim, vdim );
 
  for(int i=0; i<vdim; i++)              // we have to convert to doubles to use computeEigenvaluesAndEigenvectors function
  {
    for(int j=0; j<vdim; j++)
      JacobianD[i][j] = ( ( vectorField[Gamma] )[i][j] ).leftBound();
  }

  // temporary vectors and matrices to hold eigenvalues & imaginary parts of eigenvectors
  DVector tempvect( vdim );   
  DMatrix tempmatrix( vdim, vdim );
  
  DMatrix P( vdim, vdim );

  computeEigenvaluesAndEigenvectors(JacobianD, tempvect, tempvect, P, tempmatrix);

  // next two lines depend on dimension and mean that we are only changing coordinates for the fast variables (2x2 matrix), slow remain unchanged (are treated as a parameter)
  // here we explicitly assume vdim = 3 and last variable is slow!
  P[0][2] = P[1][2] = P[2][0] = P[2][1] = 0.;
  P[2][2] = -1.;      // -1 because we add a minus in return

  return -IMatrix(P); // minus eigenvectors are also eigenvectors and such transformed matrix suits better our computations - one could also put minuses
                      // into displacements of sections and sets to integrate from slow manifolds
};



/* ----------------------------------------------------------------------------------------- */
/* ---------------------------- ISOLATING SEGMENTS ----------------------------------------- */
/* ----------------------------------------------------------------------------------------- */

class FhnIsolatingSegment                 // class for verification of existence of isolating segments
{
public:
  IMap vectorField;
  IMatrix P;                              // diagonalization matrix along given slow manifold branch
  IVector GammaLeft;                      // slow manifold left end point 
  IVector GammaRight;                     // slow manifold right end point 
  IVector leftFace;                       // left face of the box (ys x yu centered at 0)
  IVector rightFace;                      // right face of the box (ys x yu centered at 0)
  interval disc;                          // number of discretization points
  IMatrix InvP;                           // P^(-1)
  DiscreteDynSys<IMap> vectorFieldEval;   // this is only to evaluate the vector field on C0Rect2Set in most effective way - not a real dynamical system
  IVector segmentEnclosure;               // whether we are moving to the right or to the left on the slow variable     

  FhnIsolatingSegment( IMap _vectorField, const IVector& _GammaLeft, const IVector& _GammaRight, const IMatrix& _P, const IVector& _leftFace, const IVector& _rightFace, interval _disc )
    : vectorField(_vectorField), 
      P(_P),
      GammaLeft(_GammaLeft),
      GammaRight(_GammaRight),
      leftFace(_leftFace),
      rightFace(_rightFace),
      disc(_disc),
      InvP(inverseMatrix(P)),
      vectorFieldEval(vectorField),
      segmentEnclosure( intervalHull( GammaLeft + P*leftFace, GammaRight + P*rightFace ) ) // a rough enclosure for the isolating segment to check whether
                                                                                           // slow vector field is moving in one direction only
  {
    if( !intersectionIsEmpty( IVector( {segmentEnclosure[0]} ), IVector( {segmentEnclosure[2]} ) ) )   // check whether slow vector field goes in one direction, assumes nonlinearity
      throw "ZERO OF THE SLOW SUBSYSTEM DETECTED IN ONE OF THE SEGMENTS! \n";       // is const*(u-v), const>0
  }


  // ------------- entrance verification --------------------


  IVector entranceVerification() // all normals are outward pointing
  {
    IVector normalSL( -1., 0., -( (InvP*GammaRight)[0] + rightFace[0].leftBound() - ( (InvP*GammaLeft)[0] + leftFace[0].leftBound() ) )/( GammaRight[2] - GammaLeft[2] ) );
    IVector normalSR( 1., 0., -( (InvP*GammaRight)[0] + rightFace[0].rightBound() - ( (InvP*GammaLeft)[0] + leftFace[0].rightBound() ) )/( GammaRight[2] - GammaLeft[2] ) );   
        // outward normal to (t(b-a)+a, s, t(v2-v1)+v1) is (-1,0,-(b-a)/(v2-v1)), a < 0 (left)
        // here a = (P-1(gammaleft))[0] + leftface[0].leftbound, b = (P-1(gammaright))[0] + rightface[0].leftbound so later we need to transform whole segment by P
        // to obtain the normal stable "left" vector
        // for normal stable "right" vector we do the same 
    
    normalSL = Transpose(InvP)*normalSL; // normals under affine (linear = P) transformations are transformed under inverse transpose of the transformation
    normalSR = Transpose(InvP)*normalSR;
   
    interval NormalSLxVectorField;
    interval NormalSRxVectorField;

    for(int i=1; i <= disc; i++)
    {
      interval ti = interval(i-1, i)/disc;

      IVector Gamma_i( ( GammaRight - GammaLeft )*ti + GammaLeft );
  
      // stable left evaluation

      IVector faceSL_i(3);
      faceSL_i[0] = ( rightFace[0].leftBound() - leftFace[0].leftBound() )*ti + leftFace[0].leftBound(); 
      faceSL_i[1] = interval( ( ( rightFace[1].leftBound() - leftFace[1].leftBound() )*ti + leftFace[1].leftBound() ).leftBound(), // remove some leftBounds?
                                      ( ( rightFace[1].rightBound() - leftFace[1].rightBound() )*ti + leftFace[1].rightBound() ).rightBound() ); // remove some rightBounds?
      faceSL_i[2] = 0.; 

      // stable right evaluation

      IVector faceSR_i(3);
      faceSR_i[0] = ( rightFace[0].rightBound() - leftFace[0].rightBound() )*ti + leftFace[0].rightBound();
      faceSR_i[1] = interval( (  ( rightFace[1].leftBound() - leftFace[1].leftBound() )*ti + leftFace[1].leftBound() ).leftBound(), // remove some leftBounds?
                                       ( ( rightFace[1].rightBound() - leftFace[1].rightBound() )*ti + leftFace[1].rightBound() ).rightBound() ); // remove some rightBounds?
      faceSR_i[2] = 0.;
 
      for(int j=1; j <= disc; j++)
      {
        interval tj = interval(j-1, j)/disc;

        IVector faceSL_ij(faceSL_i);
        faceSL_ij[1] = ( faceSL_i[1].rightBound() - faceSL_i[1].leftBound() )*tj + faceSL_i[1].leftBound();

        IVector faceSR_ij(faceSR_i);
        faceSR_ij[1] = ( faceSR_i[1].rightBound() - faceSR_i[1].leftBound() )*tj + faceSR_i[1].leftBound();

        C0Rect2Set CfaceSL_ij(Gamma_i, P, faceSL_ij);
        C0Rect2Set CfaceSR_ij(Gamma_i, P, faceSR_ij);
        
        CfaceSL_ij.move(vectorFieldEval);
        CfaceSR_ij.move(vectorFieldEval);
         
        IVector vectorFieldSL_ij(CfaceSL_ij);
        IVector vectorFieldSR_ij(CfaceSR_ij);

        if( i==1 && j==1 )
        {
          NormalSLxVectorField = scalarProduct(vectorFieldSL_ij, normalSL);
          NormalSRxVectorField = scalarProduct(vectorFieldSR_ij, normalSR);
        }
        else
        {
          NormalSLxVectorField = intervalHull( NormalSLxVectorField, scalarProduct(vectorFieldSL_ij, normalSL) ); 
          NormalSRxVectorField = intervalHull( NormalSRxVectorField, scalarProduct(vectorFieldSR_ij, normalSR) );
        }

      }
    }

    return IVector({NormalSLxVectorField, NormalSRxVectorField});
  } 


  // ---------------------------- exit verification -------------------------------------------


  IVector exitVerification() // all normals are outward pointing
  {    
    IVector normalUL( 0., -1., -( (InvP*GammaRight)[1] + rightFace[1].leftBound() - ( (InvP*GammaLeft)[1] + leftFace[1].leftBound() ) )/( GammaRight[2] - GammaLeft[2] ) );
    IVector normalUR( 0., 1., -( (InvP*GammaRight)[1] + rightFace[1].rightBound() - ( (InvP*GammaLeft)[1] + leftFace[1].rightBound() ) )/( GammaRight[2] - GammaLeft[2] ) );
        // again, outward normal to (s, t(b-a)+a, t(v2-v1)+v1) is (0, -1,-(b-a)/(v2-v1)) for a < 0
        // here a = (P-1(gammaleft))[1] + leftface[1].leftbound, b = (P-1(gammaright))[1] + rightface[1].leftbound so later we need to transform whole segment by P
        // same for unstable right normal ( a > 0 )

    normalUL = Transpose(InvP)*normalUL; // normals under affine transformations are transformed under inverse transpose of the transformation
    normalUR = Transpose(InvP)*normalUR;
    
    interval NormalULxVectorField;
    interval NormalURxVectorField;

    for(int i=1; i <= disc; i++)
    {
      interval ti = interval(i-1, i)/disc;

      IVector Gamma_i( ( GammaRight - GammaLeft )*ti + GammaLeft );
  
      // unstable left evaluation
          
      IVector faceUL_i(3);
      faceUL_i[0] = interval( ( ( rightFace[0].leftBound() - leftFace[0].leftBound() )*ti + leftFace[0].leftBound() ).leftBound(), // remove some leftBounds?
                                       ( ( rightFace[0].rightBound() - leftFace[0].rightBound() )*ti + leftFace[0].rightBound() ).rightBound() ); // remove some rightBounds?
      faceUL_i[1] = ( rightFace[1].leftBound() - leftFace[1].leftBound() )*ti + leftFace[1].leftBound();
      faceUL_i[2] = 0.;
 
      // unstable right evaluation

      IVector faceUR_i(3);
      faceUR_i[0] = interval( ( ( rightFace[0].leftBound() - leftFace[0].leftBound() )*ti + leftFace[0].leftBound() ).leftBound(), // remove some leftBounds?
                                       ( ( rightFace[0].rightBound() - leftFace[0].rightBound() )*ti + leftFace[0].rightBound() ).rightBound() ); // remove some rightBounds?
      faceUR_i[1] = ( rightFace[1].rightBound() - leftFace[1].rightBound() )*ti + leftFace[1].rightBound();
      faceUR_i[2] = 0.;
      
      for(int j=1; j <= disc; j++)
      {
          interval tj = interval(j-1, j)/disc;

          IVector faceUL_ij(faceUL_i);
          faceUL_ij[0] = ( faceUL_i[0].rightBound() - faceUL_i[0].leftBound() )*tj + faceUL_i[0].leftBound();

          IVector faceUR_ij(faceUR_i);
          faceUR_ij[0] = ( faceUR_i[0].rightBound() - faceUR_i[0].leftBound() )*tj + faceUR_i[0].leftBound();

          C0Rect2Set CfaceUL_ij(Gamma_i, P, faceUL_ij);
          C0Rect2Set CfaceUR_ij(Gamma_i, P, faceUR_ij);

          CfaceUL_ij.move(vectorFieldEval);
          CfaceUR_ij.move(vectorFieldEval);
      
          IVector vectorFieldUR_ij(CfaceUR_ij);
          IVector vectorFieldUL_ij(CfaceUL_ij);
           
          if( i==1 && j==1 )
          {
           NormalULxVectorField = scalarProduct(vectorFieldUL_ij, normalUL);
           NormalURxVectorField = scalarProduct(vectorFieldUR_ij, normalUR);
          }
          else
          {
           NormalULxVectorField = intervalHull( NormalULxVectorField, scalarProduct(vectorFieldUL_ij, normalUL) );
           NormalURxVectorField = intervalHull( NormalURxVectorField, scalarProduct(vectorFieldUR_ij, normalUR) );
          }
      }
    }

    return IVector({ NormalULxVectorField, NormalURxVectorField });
  }
};


class longIsolatingSegment : public FhnIsolatingSegment
{
public:
  IMatrix endP;

  longIsolatingSegment( IMap _vectorField, const IVector& _GammaLeft, const IVector& _GammaRight, const IMatrix& _P, const IMatrix& _endP, 
                        const IVector& _leftFace, const IVector& _rightFace, interval _disc )
  : FhnIsolatingSegment( _vectorField, _GammaLeft, _GammaRight, _P, _leftFace, _rightFace, _disc ),
    endP(_endP)
  // here we store an end coordinate change to be able to verify the last covering
  {
  }
  
  IVector Eq_correct(IVector& guess)              // corrects initial guesses of u so they are closer to real equilibria using Newton alg., w is always 0, computes on 
                                                  // intervals to avoid code repetition - but is not rigorous (! this is not an interval Newton operator !)
  {
    interval error(1.);
    interval result(guess[0]);
    interval oldresult(result);
    while(error > accuracy)
    {
      oldresult = result;
      result = result - (vectorField( IVector(result, 0., guess[2]) )[1])/(vectorField[ IVector(result, 0., guess[2]) ][1][0]);  
                                                                          // Newton algorithm to calculate zeroes of the vector field - w is always 0.,
                                                                          // derivative is of the second equation with respect to first variable - u
      error = abs(oldresult - result);
    }
    IVector new_Eq(3);
    new_Eq[0] = result;
    new_Eq[1] = interval(0.);
    new_Eq[2] = guess[2];
    return new_Eq;
  }
  
  IVector entranceAndExitVerification(int N_Segments) // first two coordinates are hulls of normalSLxVectorField, normalSRxVectorField, then normalULxVectorField and normalURxVectorField
    // exit and entrance verification are done together here to speed up calculations, reduce amount of code and memory used, etc.
    // N_Segments is the number of subsegments of a long isolating segment; disc is then number of discretizations of each such subsegment
    // we do not "rotate" subsegments, we also do not need to widen and shorten them to get coverings - we treat them as a part of one long
    // partially smooth IS
  {
    IVector Gamma_i0( GammaLeft );
    IVector Gamma_i1(3);

    IVector Face_i0( leftFace );
    IVector Face_i1(3);
    IVector Face_i0_adj(3); // adjusted Face_i0 (widened in stable direction, shrinked in unstable) so that there are coverings between subsegments

    IMatrix P_i0( P );
    IMatrix P_i1(3,3);

    // hulls of all normals (unstable, stable, left, right) times vector fields of all subsegments
    interval NormalULxVectorFieldHull;
    interval NormalURxVectorFieldHull;   

    interval NormalSLxVectorFieldHull; 
    interval NormalSRxVectorFieldHull;


    for(int i=1; i<=N_Segments; i++)
    {
     if( i < N_Segments )
     {
      interval ti1( double(i)/double(N_Segments) );
      Gamma_i1 = ( GammaRight - GammaLeft )*ti1 + GammaLeft;
      Gamma_i1 = Eq_correct( Gamma_i1 ); // we correct linear approx. of a slow manifold point by Newtons method

      // we widen the faces by linearly extending/contracting width and length from leftFace to rightFace sizes
      Face_i1[0] = interval( ( ( rightFace[0].leftBound() - leftFace[0].leftBound() )*ti1 + leftFace[0].leftBound() ).leftBound(), // remove some leftBounds?
                                       ( ( rightFace[0].rightBound() - leftFace[0].rightBound() )*ti1 + leftFace[0].rightBound() ).rightBound() ); // remove some rightBounds?
      Face_i1[1] = interval( ( ( rightFace[1].leftBound() - leftFace[1].leftBound() )*ti1 + leftFace[1].leftBound() ).leftBound(), // remove some leftBounds?
                                       ( ( rightFace[1].rightBound() - leftFace[1].rightBound() )*ti1 + leftFace[1].rightBound() ).rightBound() ); // remove some rightBounds?      
      Face_i1[2] = 0.;

      P_i1 = coordChange( vectorField, Gamma_i1 ); // we rotate the subsegments

      Face_i0_adj = shrinkAndExpand( Face_i0, 1.1 ); // we shrink and expand the face by a fixed constant to get covering between subsegment faces

      if( !isCovering( Face_i0, inverseMatrix(P_i1)*P_i0, Face_i0_adj ) )            // checking whether Face_i0 covers Face_i0_adj by matrix P_i1^(-1)*P_i0 (so changing coordinates
                                                                                     // from P_i0 to P_i1)
                                                                                     // this will happen if our partition into subsegments is fine enough
         throw "NO COVERING BETWEEN SUBSEGMENTS! \n";
     }
     else
     {
      Gamma_i1 = GammaRight;  // we dont need to shrink and expand here, we will arrive at the exactly same face
      Face_i1 = rightFace;
      P_i1 = endP;
     }
     
     FhnIsolatingSegment Segment_i( vectorField, Gamma_i0, Gamma_i1, P_i1, Face_i0_adj, Face_i1, disc ); 
    
     if( i == 1 )
     {
      NormalSLxVectorFieldHull = Segment_i.entranceVerification()[0]; 
      NormalSRxVectorFieldHull = Segment_i.entranceVerification()[1];
      NormalULxVectorFieldHull = Segment_i.exitVerification()[0];
      NormalURxVectorFieldHull = Segment_i.exitVerification()[1];
     }
     else
     {
      NormalSLxVectorFieldHull = intervalHull( NormalSLxVectorFieldHull, Segment_i.entranceVerification()[0] ); 
      NormalSRxVectorFieldHull = intervalHull( NormalSRxVectorFieldHull, Segment_i.entranceVerification()[1] );
      NormalULxVectorFieldHull = intervalHull( NormalULxVectorFieldHull, Segment_i.exitVerification()[0] );
      NormalURxVectorFieldHull = intervalHull( NormalURxVectorFieldHull, Segment_i.exitVerification()[1] );
     }

     Gamma_i0 = Gamma_i1;  // we move to the next subsegment
     Face_i0 = Face_i1;
     P_i0 = P_i1;
    }

    return IVector({ NormalSLxVectorFieldHull, NormalSRxVectorFieldHull, NormalULxVectorFieldHull, NormalURxVectorFieldHull });
  }
  


};


