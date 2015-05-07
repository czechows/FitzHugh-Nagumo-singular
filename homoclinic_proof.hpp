/* -----------------------------------------------------------------------------------------
 * This is a header file to fhn.cpp providing a rigorous procedure 
 * of verification of existence of homoclinic orbits
 * in the FitzHugh-Nagumo system for given paramet eps, close to a given theta
 * (a codim 1 phenomenon, a guess of theta is updated to a range of theta containing the pulse). 
 * * ----------------------------------------------------------------------------------------*/

/* --------------------------------------------------------------------------------------- */
/* ----- VERIFICATION OF EXISTENCE OF A HOMOCLINIC ORBIT FOR GIVEN PARAMETER VALUES ------ */
/* --------------------------------------------------------------------------------------- */


void FhnVerifyExistenceOfHomoclinicOrbit( interval _theta, interval _eps, bool _verbose = 0, bool withParams = 0, interval _thetaVar=interval(-2.5e-3,2.5e-3), int _pMapDivCount = 20, 
     int _chainSubsegmentCountU = 200, int _chainSubsegmentCountD =400, int _chainSegmentDivCount = 110, int _cornerSegmentDivCount = 150 ) 
  // verbose on displays all the interval enclosures for Poincare maps / products of the vector fields with normals; other parameters control respectively: 
  // number of subdivisions of sets to integrate (in each dimension), number of subsegments along slow manifolds, number of subdivisions of regular/corner segments
  // for evaluation of the scalar product of the vector field with outward pointing normals
{
  try            // we check negations of all assumptions to throw exceptions, if no exception is thrown existence of the orbit is verified
  {
    IVector parameters({ _theta, _eps });

    IVector GammaUL(0.970345591417269, 0., 0.);                                                   // some guesses for the corner points which are equilibria
    IVector GammaDL(0., 0., 0.);                                                                  // of the fast subsystem for critical parameter v values (third variable)
                                                                                                  // where heteroclinics exist
    IVector GammaUR(1.0, 0., 0.12);                                                               // UR up right, DR down right, UL up left, DL down left
    IVector GammaDR(-0.3, 0., 0.12);

    GammaHom_correct( _theta, GammaUL, GammaDL, GammaUR, GammaDR );                              // we correct the initial guesses by nonrigorous shooting methods (see numerics.hpp)

    GammaDL = IVector({0.,0.,0.});       // this we know, i.e. we can compute the stationary point analytically
    GammaUL[2] = 0.;

    cout << "Initial guesses: \n";
    cout << "theta = " << _theta << "\nGammaUL = " << GammaUL << "\nGammaDL = " << GammaDL << "\nGammaUR = " << GammaUR << "\nGammaDR = " << GammaDR << "\n \n";
 
    if( !(GammaUL[0] > GammaDL[0] && GammaUR[0] > GammaDR[0] && GammaUR[2] > GammaUL[2] && GammaDR[2] > GammaDL[2] ) )
      throw "NEWTON CORRECTION METHOD FOR CORNER POINTS ERROR! \n";

    _theta = _theta + _thetaVar; // for shooting with theta

    (*Fhn_vf).setParameter("theta",_theta);
    (*Fhn_vf).setParameter("eps",_eps);
 
    (*Fhn_vf_rev).setParameter("theta",_theta);
    (*Fhn_vf_rev).setParameter("eps",_eps);

    IMatrix PUL( coordChange( *Fhn_vf, GammaUL ) ), 
            PUR( coordChange( *Fhn_vf, GammaUR ) ), 
            PDL( coordChange( *Fhn_vf, GammaDL ) ),  
            PDR( coordChange( *Fhn_vf, GammaDR ) ); 

    IVector setToIntegrateDL(2);
    IVector setToIntegrateUR(2);
    IVector setToBackIntegrateUL(2);
    IVector setToBackIntegrateDR(2);

    // below sizes of blocks and distances from appropriate sections in appropr. direction (stable for sections to integrate from, unstable for sections to integrate onto)

    // the stable manifold block sizes
    double sMan_ruDL(5.0e-4);          
    double sMan_rsDL(5.0e-4);
    interval sMan_vDL(6.0e-4);

    // the unstable manifold block sizes
    interval ruDL(8.0e-5);    // this is yu at the downleft corner      
    setToIntegrateDL[0] = 2.0e-5*interval(-1,1);  // this is ys at downleft corner  
    setToIntegrateDL[1] = 1.0e-5*interval(-1,1);  // this is v at downleft corner

    interval ruUR(5.0e-3);  // this is yu at the upright corner
    setToIntegrateUR[0] = 2.0e-3*interval(-1,1);  // this is ys at upright corner
    setToIntegrateUR[1] = 7.0e-4*interval(-1,1);  // this is v at upright corner

 /*   interval rsUL(1.5e-2);   // this is ys at the upleft corner THESE SIZES ARE DEPRECATED?
    setToBackIntegrateUL[0] = 5.0e-3*interval(-1,1);     // this is v at upleft corner
    setToBackIntegrateUL[1] = 1.0e-2*interval(-1,1);         // this is yu at upleft corner 1e-2 works for isolation */

    interval rsUL(1.5e-3);   // this is ys at the upleft corner
    setToBackIntegrateUL[0] = 7.0e-4*interval(-1,1);     // this is v at upleft corner
    setToBackIntegrateUL[1] = 5.0e-4*interval(-1,1);         // this is yu at upleft corner 

    interval rsDR(1.0e-2);   // this is ys at the downright corner
    setToBackIntegrateDR[0] = 2.0e-3*interval(-1,1);     // this is v at downright corner 
    setToBackIntegrateDR[1] = 2.0e-3*interval(-1,1);        // this is yu at downright corner


    IVector URface( setToIntegrateUR[0], ruUR*interval(-1,1), 0. );
    IVector DRface( rsDR*interval(-1,1), setToBackIntegrateDR[1], 0. ); 

    IVector GammaUR_left( Eq_correct( *Fhn_vf, GammaUR + IVector( 0., 0., setToIntegrateUR[1].leftBound() ) ) );
    GammaUR_left[2] = GammaUR[2] + setToIntegrateUR[1].leftBound();
    IVector GammaUR_right( Eq_correct( *Fhn_vf, GammaUR + IVector( 0., 0., setToIntegrateUR[1].rightBound() ) ) );
    GammaUR_right[2] = GammaUR[2] + setToIntegrateUR[1].rightBound();

    IVector GammaDR_left( Eq_correct( *Fhn_vf, GammaDR + IVector( 0., 0., setToBackIntegrateDR[0].leftBound() ) ) );
    GammaDR_left[2] = GammaDR[2] + setToBackIntegrateDR[0].leftBound();
    IVector GammaDR_right( Eq_correct( *Fhn_vf, GammaDR + IVector( 0., 0., setToBackIntegrateDR[0].rightBound() ) ) );
    GammaDR_right[2] = GammaDR[2] + setToBackIntegrateDR[0].rightBound();

    FhnIsolatingSegment URSegment( *Fhn_vf, GammaUR_left, GammaUR_right, PUR, URface, URface, _cornerSegmentDivCount );  
    FhnIsolatingSegment DRSegment( *Fhn_vf, GammaDR_left, GammaDR_right, PDR, DRface, DRface, _cornerSegmentDivCount );  

    midPoincareMap *rightMap; // this class implements the Poincare maps described in the paper as pmUR, pmDR onto rightSection

 /*   if( withParams )
    {
      rightMap = new midPoincareMap( parameters, *Fhn_vf_withParams, *Fhn_vf_withParams_rev, PUR, PDR, GammaUR, GammaDR, ruUR, rsDR, 1., _pMapDivCount );
    }
    else
    {*/
     rightMap = new midPoincareMap( *Fhn_vf, *Fhn_vf_rev, URSegment, DRSegment, _theta, _eps, 1., _pMapDivCount );
   // }

    // covering checks 

    if( _verbose )
      cout << "\n ------------------- RIGHT SIDE COVERING CHECKS: -------------------------- \n \n";
    if( !(*rightMap).checkCovering( setToIntegrateUR, setToBackIntegrateDR, _verbose )  )    
      throw "FAILURE TO CHECK COVERINGS IN THE FAST REGIME (RIGHT MAP)! \n";

    delete rightMap;
    // right isolating segments


    IVector URSegment_entranceVerification( URSegment.entranceVerification() );
    IVector URSegment_exitVerification( URSegment.exitVerification() );
    IVector DRSegment_entranceVerification( DRSegment.entranceVerification() );
    IVector DRSegment_exitVerification( DRSegment.exitVerification() );

    if( _verbose )
    {
      cout << "\n ------------------- UR, DR SEGMENTS ISOLATION: -------------------------- \n \n";

      cout << "Enclosures of scalar product of the vector field with entrance faces normals for UR segment: \n \n" << URSegment_entranceVerification << "\n";  
      cout << "\n --- \n";
      cout << "Enclosures of scalar product of the vector field with exit faces normals for UR segment: \n \n" << URSegment_exitVerification << "\n";

      cout << "\n --- \n";

      cout << "Enclosures of scalar product of the vector field with entrance faces normals for DR segment: \n \n" << DRSegment_entranceVerification << "\n";  
      cout << "\n --- \n";
      cout << "Enclosures of scalar product of the vector field with exit faces normals for DR segment: \n \n" << DRSegment_exitVerification << "\n";


      cout << "\n --- \n";
      cout << "\n --- \n";
    };

    if( !( URSegment_entranceVerification[0] < 0. && URSegment_entranceVerification[1] < 0. && URSegment_exitVerification[0] > 0. && URSegment_exitVerification[1] > 0. ) )
      throw "ISOLATION ERROR FOR UR CORNER SEGMENT! \n";
    if( !( DRSegment_entranceVerification[0] < 0. && DRSegment_entranceVerification[1] < 0. && DRSegment_exitVerification[0] > 0. && DRSegment_exitVerification[1] > 0. ) )
      throw "ISOLATION ERROR FOR DR CORNER SEGMENT! \n";


    // left isolating segments

    IVector ULface( rsUL*interval(-1,1), setToBackIntegrateUL[1], 0. ); 
    IVector uManDLface( setToIntegrateDL[0], ruDL*interval(-1,1), 0. ); 
    IVector sManDLface( interval(-sMan_rsDL, sMan_rsDL), interval(-sMan_ruDL, sMan_ruDL), 0. );
 
    IVector uManGammaDL_left( Eq_correct( *Fhn_vf, GammaDL + IVector( 0., 0., setToIntegrateDL[1].leftBound() ) ) );
    uManGammaDL_left[2] = GammaDL[2] + setToIntegrateDL[1].leftBound();
    IVector uManGammaDL_right( Eq_correct( *Fhn_vf, GammaDL + IVector( 0., 0., setToIntegrateDL[1].rightBound() ) ) );
    uManGammaDL_right[2] = GammaDL[2] + setToIntegrateDL[1].rightBound();

    IVector sManGammaDL_left( Eq_correct( *Fhn_vf, GammaDL + IVector( 0., 0., -sMan_vDL ) ) );
    sManGammaDL_left[2] = GammaDL[2] - sMan_vDL;
    IVector sManGammaDL_right( Eq_correct( *Fhn_vf, GammaDL + IVector( 0., 0., sMan_vDL ) ) );
    sManGammaDL_right[2] = GammaDL[2] + sMan_vDL;


    IVector GammaUL_left( Eq_correct( *Fhn_vf, GammaUL + IVector( 0., 0., setToBackIntegrateUL[0].leftBound() ) ) );
    GammaUL_left[2] = GammaUL[2] + setToBackIntegrateUL[0].leftBound();
    IVector GammaUL_right( Eq_correct( *Fhn_vf, GammaUL + IVector( 0., 0., setToBackIntegrateUL[0].rightBound() ) ) );
    GammaUL_right[2] = GammaUL[2] + setToBackIntegrateUL[0].rightBound();

   // cout << "uManGammaDL left = " << uManGammaDL_left << "\n uManGammaDL right " << uManGammaDL_right << "\n";

    FhnIsolatingSegment ULSegment( *Fhn_vf, GammaUL_left, GammaUL_right, PUL, ULface, ULface, _cornerSegmentDivCount ); 

    FhnIsolatingBlock uManDLBlock( *Fhn_vf, uManGammaDL_left, uManGammaDL_right, PDL, uManDLface, uManDLface, _cornerSegmentDivCount );  
    FhnIsolatingBlock sManDLBlock( *Fhn_vf, sManGammaDL_left, sManGammaDL_right, PDL, sManDLface, sManDLface, _cornerSegmentDivCount );  

    // left Poincare map - shooting with theta
    
    midPoincareMap *leftMap; // this class implements the Poincare maps described in the paper as pmUL, pmDL onto leftSection

    leftMap = new midPoincareMap( *Fhn_vf, *Fhn_vf_rev, uManDLBlock, ULSegment, _theta, _eps, -1., _pMapDivCount );

    // covering checks 

    if( _verbose )
      cout << "\n ------------------- LEFT SIDE COVERING CHECKS: -------------------------- \n \n";
    if( !(*leftMap).shootWithTheta( setToIntegrateDL, setToBackIntegrateUL )  )    
      throw "FAILURE TO CHECK COVERINGS IN THE FAST REGIME (LEFT MAP)! \n";

    delete leftMap;


    IVector ULSegment_entranceVerification( ULSegment.entranceVerification() );
    IVector ULSegment_exitVerification( ULSegment.exitVerification() );

    IVector uManDLBlock_entranceVerification( uManDLBlock.entranceVerification() );
    IVector uManDLBlock_exitVerification( uManDLBlock.exitVerification() );

    IVector sManDLBlock_entranceVerification( sManDLBlock.entranceVerification() );
    IVector sManDLBlock_exitVerification( sManDLBlock.exitVerification() );

    if( _verbose )
    {
      cout << "\n ------------------- UL, DL SEGMENTS/BLOCKS ISOLATION: -------------------------- \n \n";

      cout << "Enclosures of scalar product of the vector field with entrance faces normals for uManDL block: \n \n" << uManDLBlock_entranceVerification << "\n";  
      cout << "\n --- \n";
      cout << "Enclosures of scalar product of the vector field with exit faces normals for uManDL block: \n \n" << uManDLBlock_exitVerification << "\n";

      cout << "\n --- \n";
 
      cout << "Enclosures of scalar product of the vector field with entrance faces normals for sManDL block: \n \n" << sManDLBlock_entranceVerification << "\n";  
      cout << "\n --- \n";
      cout << "Enclosures of scalar product of the vector field with exit faces normals for sManDL block: \n \n" << sManDLBlock_exitVerification << "\n";

      cout << "\n --- \n";
      cout << "Enclosures of scalar product of the vector field with entrance faces normals for UL segment: \n \n" << ULSegment_entranceVerification << "\n";  
      cout << "\n --- \n";
      cout << "Enclosures of scalar product of the vector field with exit faces normals for UL segment: \n \n" << ULSegment_exitVerification << "\n";


      cout << "\n --- \n";
      cout << "\n --- \n";
    };

    if( !( ULSegment_entranceVerification[0] < 0. && ULSegment_entranceVerification[1] < 0. && ULSegment_exitVerification[0] > 0. && ULSegment_exitVerification[1] > 0. ) )
      throw "ISOLATION ERROR FOR UL CORNER SEGMENT! \n";
    if( !( uManDLBlock_entranceVerification[0] < 0. && uManDLBlock_entranceVerification[1] < 0. && uManDLBlock_exitVerification[0] > 0. && uManDLBlock_exitVerification[1] > 0. ) )
      throw "ISOLATION ERROR FOR uManDL CORNER BLOCK! \n";
    if( !( sManDLBlock_entranceVerification[0] < 0. && sManDLBlock_entranceVerification[1] < 0. && sManDLBlock_exitVerification[0] > 0. && sManDLBlock_exitVerification[1] > 0. ) )
      throw "ISOLATION ERROR FOR sManDL CORNER BLOCK! \n";


    // up down isolating segments

    chainOfSegments UpSegment( *Fhn_vf, ULSegment.GammaRight, URSegment.GammaLeft, PUL, PUR, ULface, URface, _chainSegmentDivCount );
    chainOfSegments DownSegment( *Fhn_vf, sManDLBlock.GammaRight, DRSegment.GammaLeft, PDL, PDR, sManDLface, DRface, _chainSegmentDivCount ); 

    if( !( ULSegment.segmentEnclosure[0] > ULSegment.segmentEnclosure[2] ) )
      throw "MISALIGNMENT OF ONE OF THE UPPER SEGMENTS! \n";
    if( !( DRSegment.segmentEnclosure[0] < DRSegment.segmentEnclosure[2] ) )
      throw "MISALIGNMENT OF ONE OF THE LOWER SEGMENTS! \n";      // checks on whether we are above/below u=v plane for upper/lower segments

    IVector UpSegment_entranceAndExitVerification( UpSegment.entranceAndExitVerification( _chainSubsegmentCountU ) );
    IVector DownSegment_entranceAndExitVerification( DownSegment.entranceAndExitVerification( _chainSubsegmentCountD ) );

    if( _verbose )
    {
      cout << "\n ---------------------------- UPPER, LOWER CHAINS OF SEGMENTS ISOLATION: ---------------------------- \n \n";

      cout << "Interval hull of enclosures of scalar products of the vector field with the upper chain of segments (not including corner ones, left/right entrance faces first, then exit faces): \n \n"
        << UpSegment_entranceAndExitVerification << "\n";  
      cout << "\n --- \n";
      cout << "Interval hull of enclosures of scalar products of the vector field with the lower chain segments (not including corner ones, left/right entrance faces first, then exit faces): \n \n"
        << DownSegment_entranceAndExitVerification << "\n \n";  
 
      cout << "\n --- \n";
      cout << "\n --- \n";   
    };

   //  ISOLATION CAN BE ALSO CHECKED NOW FOR EACH SUBSEGMENT IN SEGMENTS.HPP TO THROW AN EXCEPTION QUICKER, SEE COMMENTED LINES IN SEGMENTS.HPP
    if( !( UpSegment_entranceAndExitVerification[0] < 0. && UpSegment_entranceAndExitVerification[1] < 0. 
          && UpSegment_entranceAndExitVerification[2] > 0. && UpSegment_entranceAndExitVerification[3] > 0. ) )
      throw "ISOLATION ERROR FOR ONE OF THE UPPER REGULAR SEGMENTS! \n";
    if( !( DownSegment_entranceAndExitVerification[0] < 0. && DownSegment_entranceAndExitVerification[1] < 0. 
          && DownSegment_entranceAndExitVerification[2] > 0. && DownSegment_entranceAndExitVerification[3] > 0. ) )
      throw "ISOLATION ERROR FOR ONE OF THE LOWER REGULAR SEGMENTS! \n";

    cout << "Existence of a periodic orbit for the FitzHugh-Nagumo system with parameter values theta=" << _theta << " and eps=" << _eps << " verified! \n";



  }
  catch(const char* Message)
  {
    cout << Message << "EXISTENCE OF A HOMOCLINIC ORBIT FOR PARAMETER VALUES THETA=" << _theta << " AND EPS=" << _eps << " NOT VERIFIED! \n";
  }


};
