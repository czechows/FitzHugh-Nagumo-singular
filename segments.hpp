/* -----------------------------------------------------------------------------------------
 * This is a header file to fhn.cpp providing classes for isolating segments: FhnIsolatingSegment,
 * which is used when a constructed segment is short and one can stick with one affine coordinate change
 * to check the isolation and chainOfSegments - a derived class
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
  int vdim( 3 );   // should be used only in dimension 3!
  DMatrix JacobianD( vdim, vdim );

  // a patch to set eps to 0 for the vector field for computing coordinates around slow manifold
  IMap vectorFieldZeroEps( vectorField );          
  vectorFieldZeroEps.setParameter("eps", interval(0.) );
  // WARNING: specific to the (type of) the vector field. Takes care of the problem that proof does not go for subintervals of parameter epsilon away from 0.
 
  for(int i=0; i<vdim; i++)              // we have to convert to doubles to use computeEigenvaluesAndEigenvectors function
  {
    for(int j=0; j<vdim; j++)
      JacobianD[i][j] = ( ( vectorFieldZeroEps[Gamma] )[i][j] ).leftBound();
  }

  // temporary vectors and matrices to hold eigenvalues & imaginary parts of eigenvectors
  DVector tempvectRe( vdim );   
  DVector tempvectIm( vdim );

  DMatrix tempmatrix( vdim, vdim );
  
  DMatrix P( vdim, vdim );
  DMatrix P_result( vdim, vdim );

  computeEigenvaluesAndEigenvectors(JacobianD, tempvectRe, tempvectIm, P, tempmatrix);

  int i_min(42);  // i's initialized with any value
  int i_med(42);
  int i_max(42); 

  for( int i = 1; i <= vdim; i++ ) // sorting so we are sure we have stable coordinate first unstable second neutral third. ONLY FOR 3D vector fields!
  {
    if( tempvectRe(i) == min( tempvectRe(1), min(tempvectRe(2), tempvectRe(3)) ) )
      i_min = i;
    else if( tempvectRe(i) == max( tempvectRe(1), max(tempvectRe(2), tempvectRe(3)) ) )
      i_max = i;
    else 
      i_med = i;
  }


  for( int i = 1; i <=vdim; i++ )
  { 
    P_result(i,1) = P( i, i_min );
    P_result(i,2) = P( i, i_max );
    P_result(i,3) = P( i, i_med );
  }

  // next two lines depend on the dimension and mean that we are only changing coordinates for the fast variables (2x2 matrix), slow remain unchanged (are treated as a parameter)
  // here we explicitly assume vdim = 3 and last variable is slow!
  P_result[0][2] = P_result[1][2] = P_result[2][0] = P_result[2][1] = 0.;
  P_result[2][2] = -1.;      // -1 because we add a minus in return

  //cout << P_result << "\n";
  //cout << "eigval: " << tempvectRe << "\n";


  return -IMatrix(P_result); // minus eigenvectors are also eigenvectors and such transformed matrix suits better our computations - one could also put minuses
                            // into displacements of sections and sets to integrate from slow manifolds
};


/* ----------------------------------------------------------------------------------------- */
/* ---------------------------- ISOLATING SEGMENTS ----------------------------------------- */
/* ----------------------------------------------------------------------------------------- */

class FhnIsolatingSegment                 // class for verification of isolation in segments
{
public:
  IMap vectorField;
  IMatrix P;                              // diagonalization matrix along given slow manifold branch
  IVector GammaLeft;                      // slow manifold left end point - this variable name is Front in the paper for segments with u>v and Rear elsewise 
  IVector GammaRight;                     // slow manifold right end point - this variable name is Rear in the paper for segments with u>v and Front elsewise
  IVector leftFace;                       // left face of the box (ys x yu centered at 0) = [-b,b] x [-a,a] in the paper for segments with u>v and [-d,d] x [-c,c] elsewise
  IVector rightFace;                      // right face of the box (ys x yu centered at 0) = [-d,d] x [-c,c] in the paper for segments with u>v and [-b,b] x [-a,a] elsewise
  interval div;                          // number of subdivisions
  IMatrix InvP;                           // P^(-1)
  DiscreteDynSys<IMap> vectorFieldEval;   // this is only to evaluate the vector field on C0Rect2Set in most effective way - not a real dynamical system
  IVector segmentEnclosure;               // whether we are moving to the right or to the left on the slow variable     
  bool isABlock;

  FhnIsolatingSegment( IMap _vectorField, const IVector& _GammaLeft, const IVector& _GammaRight, const IMatrix& _P, const IVector& _leftFace, const IVector& _rightFace, interval _div,
      bool _isABlock = 0 )
    : vectorField(_vectorField), 
      P(_P),
      GammaLeft(_GammaLeft),
      GammaRight(_GammaRight),
      leftFace(_leftFace),
      rightFace(_rightFace),
      div(_div),
      InvP(inverseMatrix(P)),
      vectorFieldEval(vectorField),
      segmentEnclosure( intervalHull( GammaLeft + P*leftFace, GammaRight + P*rightFace ) ), // a rough enclosure for the isolating segment to check whether
                                                                                           // slow vector field is moving in one direction only
      isABlock( _isABlock )
  {
    // check whether the slow vector field goes in one direction, assumes the nonlinearity is const*(u-v), const>0
    if( !isABlock )
    {
      if( !intersectionIsEmpty( IVector( {segmentEnclosure[0]} ), IVector( {segmentEnclosure[2]} ) ) )
        throw "ZERO OF THE SLOW SUBSYSTEM DETECTED IN ONE OF THE SEGMENTS! \n";       
    }
  }


  // ------------- entrance verification --------------------


  IVector entranceVerification() // all normals are outward pointing
  {
    IVector normalSL( -1., 0., ( (InvP*GammaRight)[0] + rightFace[0].leftBound() - ( (InvP*GammaLeft)[0] + leftFace[0].leftBound() ) )/( GammaRight[2] - GammaLeft[2] ) );
    IVector normalSR( 1., 0., -( (InvP*GammaRight)[0] + rightFace[0].rightBound() - ( (InvP*GammaLeft)[0] + leftFace[0].rightBound() ) )/( GammaRight[2] - GammaLeft[2] ) );   
        // ex. outward normal to (t(b-a)+a, s, t(v2-v1)+v1) is (1,0,-(b-a)/(v2-v1)) (for normalSR, with - for normalSL)
        // here a = (P-1(gammaleft))[0] + leftface[0].leftbound, b = (P-1(gammaright))[0] + rightface[0].leftbound so later we need to transform whole segment by P
        // to obtain the normal stable "left" vector
        // for normal stable "right" vector we do the same 
    
    normalSL = Transpose(InvP)*normalSL; // normals under affine (linear = P) transformations are transformed under inverse transpose of the transformation
    normalSR = Transpose(InvP)*normalSR;
   
    interval NormalSLxVectorField;
    interval NormalSRxVectorField;

    for(int i=1; i <= div; i++)
    {
      interval ti = interval(i-1, i)/div;

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
 
      for(int j=1; j <= div; j++)
      {
        interval tj = interval(j-1, j)/div;

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
    IVector normalUL( 0., -1., ( (InvP*GammaRight)[1] + rightFace[1].leftBound() - ( (InvP*GammaLeft)[1] + leftFace[1].leftBound() ) )/( GammaRight[2] - GammaLeft[2] ) );
    IVector normalUR( 0., 1., -( (InvP*GammaRight)[1] + rightFace[1].rightBound() - ( (InvP*GammaLeft)[1] + leftFace[1].rightBound() ) )/( GammaRight[2] - GammaLeft[2] ) );
        // again, outward normal to (s, t(b-a)+a, t(v2-v1)+v1) is (0, 1,-(b-a)/(v2-v1)) for UR, minus that for UL
        // here a = (P-1(gammaleft))[1] + leftface[1].leftbound, b = (P-1(gammaright))[1] + rightface[1].leftbound so later we need to transform whole segment by P
        // same for unstable right normal ( a > 0 )

    normalUL = Transpose(InvP)*normalUL; // normals under affine transformations are transformed under inverse transpose of the transformation
    normalUR = Transpose(InvP)*normalUR;
    
    interval NormalULxVectorField;
    interval NormalURxVectorField;

    for(int i=1; i <= div; i++)
    {
      interval ti = interval(i-1, i)/div;

      IVector Gamma_i( ( GammaRight - GammaLeft )*ti + GammaLeft );
  
      // unstable left evaluation
          
      IVector faceUL_i(3);
      faceUL_i[0] = interval( ( ( rightFace[0].leftBound() - leftFace[0].leftBound() )*ti + leftFace[0].leftBound() ).leftBound(), 
                                       ( ( rightFace[0].rightBound() - leftFace[0].rightBound() )*ti + leftFace[0].rightBound() ).rightBound() ); 
      faceUL_i[1] = ( rightFace[1].leftBound() - leftFace[1].leftBound() )*ti + leftFace[1].leftBound();
      faceUL_i[2] = 0.;
 
      // unstable right evaluation

      IVector faceUR_i(3);
      faceUR_i[0] = interval( ( ( rightFace[0].leftBound() - leftFace[0].leftBound() )*ti + leftFace[0].leftBound() ).leftBound(), 
                                       ( ( rightFace[0].rightBound() - leftFace[0].rightBound() )*ti + leftFace[0].rightBound() ).rightBound() ); 
      faceUR_i[1] = ( rightFace[1].rightBound() - leftFace[1].rightBound() )*ti + leftFace[1].rightBound();
      faceUR_i[2] = 0.;
      
      for(int j=1; j <= div; j++)
      {
          interval tj = interval(j-1, j)/div;

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


// ------------------------ class chainOfSegments ------------------------


class chainOfSegments : public FhnIsolatingSegment
{
public:
  IMatrix endP;

  chainOfSegments( IMap _vectorField, const IVector& _GammaLeft, const IVector& _GammaRight, const IMatrix& _P, const IMatrix& _endP, 
                        const IVector& _leftFace, const IVector& _rightFace, interval _div )
  : FhnIsolatingSegment( _vectorField, _GammaLeft, _GammaRight, _P, _leftFace, _rightFace, _div ),
    endP(_endP)
  // here we store an end coordinate change to be able to verify the last covering
  {
  }
  

  
  IVector entranceAndExitVerification(int N_Segments) // first two coordinates are hulls of normalSLxVectorField, normalSRxVectorField, then normalULxVectorField and normalURxVectorField
    // exit and entrance verification are done together here to speed up calculations, reduce amount of code and memory used, etc.
    // N_Segments is the number of subsegments of a long isolating segment; div is then number of subdivisions of each such subsegment in each direction
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
      Gamma_i1 = Eq_correct( vectorField, Gamma_i1 ); // we correct linear approx. of a slow manifold point by Newtons method

      // we widen the faces by linearly extending/contracting width and length from leftFace to rightFace sizes
      Face_i1[0] = interval( ( ( rightFace[0].leftBound() - leftFace[0].leftBound() )*ti1 + leftFace[0].leftBound() ).leftBound(), 
                                       ( ( rightFace[0].rightBound() - leftFace[0].rightBound() )*ti1 + leftFace[0].rightBound() ).rightBound() ); 
      Face_i1[1] = interval( ( ( rightFace[1].leftBound() - leftFace[1].leftBound() )*ti1 + leftFace[1].leftBound() ).leftBound(), 
                                       ( ( rightFace[1].rightBound() - leftFace[1].rightBound() )*ti1 + leftFace[1].rightBound() ).rightBound() );       
      Face_i1[2] = 0.;

      P_i1 = coordChange( vectorField, Gamma_i1 ); // we rotate the subsegments

      Face_i0_adj = shrinkAndExpand( Face_i0, 1.05 ); // we shrink and expand the face by a fixed constant to get covering between subsegment faces

      if( !isCovering( Face_i0, inverseMatrix(P_i1)*P_i0, Face_i0_adj ) )            // checking whether Face_i0 covers Face_i0_adj by matrix P_i1^(-1)*P_i0 (so changing coordinates
                                                                                     // from P_i0 to P_i1)
                                                                                     // this will happen if our partition into subsegments is fine enough
      {
        cout << "ERROR AT i= " << i << "\n";
        cout << "Face_i0 = " << Face_i0 << "\n";
        cout << "inverseMatrix(P_i1)*P_i0: " << inverseMatrix(P_i1)*P_i0 << "\n";
        cout << "Face_i0_adj: " << Face_i0_adj << "\n"; 
         throw "NO COVERING BETWEEN SUBSEGMENTS ! \n";
      }
     }
     else
     {
      Gamma_i1 = GammaRight;  // we dont need to shrink and expand here, we will arrive at the exactly same face
      Face_i1 = rightFace;
      P_i1 = endP;
     }
     
     FhnIsolatingSegment Segment_i( vectorField, Gamma_i0, Gamma_i1, P_i1, Face_i0_adj, Face_i1, div ); 
    
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

     // UNCOMMENT FOR THROWING ISOLATION EXCEPTIONS ON THE RUN TO BREAK FROM PROGRAM FASTER - NOT NECESSARY FOR THE PROOF BUT SAVES TIME
     if( vectalg::containsZero( IVector( {intervalHull( NormalSLxVectorFieldHull, NormalSRxVectorFieldHull )} ) ) )
     {
       cout << "intervalHull( NormalSLxVectorFieldHull, NormalSRxVectorFieldHull ) = " << intervalHull( NormalSLxVectorFieldHull, NormalSRxVectorFieldHull ) << "\n";
       cout << "NO ISOLATION AT SEGMENT i= " << i << "\n";
       cout.flush();

       throw "ISOLATION ERROR FOR ONE OF THE REGULAR ISOLATING SEGMENTS! \n" ; 
     }

     if( vectalg::containsZero( IVector( {intervalHull( NormalULxVectorFieldHull, NormalURxVectorFieldHull )} ) ) )
     {
       cout << "intervalHull( NormalULxVectorFieldHull, NormalURxVectorFieldHull ) = " << intervalHull( NormalULxVectorFieldHull, NormalURxVectorFieldHull ) << "\n";
       cout << "NO ISOLATION AT SEGMENT i= " << i << "\n";
       cout.flush();

       throw "ISOLATION ERROR FOR ONE OF THE REGULAR ISOLATING SEGMENTS! \n" ; 
     }

     Gamma_i0 = Gamma_i1;  // we move to the next subsegment
     Face_i0 = Face_i1;
     P_i0 = P_i1;
    }

    return IVector({ NormalSLxVectorFieldHull, NormalSRxVectorFieldHull, NormalULxVectorFieldHull, NormalURxVectorFieldHull });
  }
  


};


class FhnIsolatingBlock : public FhnIsolatingSegment
{
public:
  FhnIsolatingBlock( IMap _vectorField, const IVector& _GammaLeft, const IVector& _GammaRight, const IMatrix& _P, const IVector& _leftFace, const IVector& _rightFace, interval _div )
  : FhnIsolatingSegment( _vectorField, _GammaLeft, _GammaRight, _P, _leftFace, _rightFace, _div, 1 ) // it is a block so we don't check whether VF is uniform in one direction
  {
    // to obtain isolation in the central (second entry) direction we check u>v on the left slow face and u<v on the right slow face, ONLY FOR THE FITZHUGH-NAGUMO VECTOR FIELD!
    if( !( ( GammaRight + P*rightFace )[0] < GammaRight[2] && ( GammaLeft + P*leftFace )[0] > GammaLeft[2] ) )  
    {
      throw "No isolation in the slow direction for the isolating block! \n";
    }
  }
};
