
/* -----------------------------------------------------------------------------------------
 * This is a header file to fhn.cpp providing a rigorous procedure 
 * of verification of existence of periodic orbits
 * in the FitzHugh-Nagumo system for given parameters theta, eps. 
 * * ----------------------------------------------------------------------------------------*/

/* ------------------------------------------------------------------------------------ */
/* ----- VERIFICATION OF EXISTENCE OF PERIODIC ORBITS FOR GIVEN PARAMETER VALUES ------ */
/* ------------------------------------------------------------------------------------ */


void FhnVerifyExistenceOfPeriodicOrbit( interval _theta, interval _eps, bool _verbose = 0, bool withParams = 0, int _pMapDivCount = 20, 
     int _chainSubsegmentCount = 80, int _chainSegmentDivCount = 110, int _cornerSegmentDivCount = 150 ) 
  // verbose on displays all the interval enclosures for Poincare maps / products of vector fields with normals; other parameters control respectively: 
  // number of subdivisions of sets to integrate (in each dimension), number of subsegments along slow manifolds, number of subdivisions of regular/corner segments
  // for evaluation of the scalar product of vector field with outward pointing normals
{
  try                   // we check negations of all assumptions to throw exceptions, if no exception is thrown existence of the orbit is verified
  {
    (*Fhn_vf).setParameter("theta",_theta);
    (*Fhn_vf).setParameter("eps",_eps);
 
    (*Fhn_vf_rev).setParameter("theta",_theta);
    (*Fhn_vf_rev).setParameter("eps",_eps);

    IVector parameters({ _theta, _eps });

    IVector GammaUL(0.970345591417269, 0., 0.0250442158334208);                                   // some guesses for the corner points which are equilibria
    IVector GammaDL(-0.108412947498862, 0., 0.0250442158334208);                                  // of the fast subsystem for critical parameter v values (third variable)
                                                                                                  // where heteroclinics exist
    IVector GammaUR(0.841746280832201, 0., 0.0988076360184288);                                   // UR up right, DR down right, UL up left, DL down left
    IVector GammaDR(-0.237012258083933, 0., 0.0988076360184288);

    GammaQuad_correct( _theta, GammaUL, GammaDL, GammaUR, GammaDR );                              // we correct the initial guesses by nonrigorous shooting methods (see numerics.hpp)

    if( !(GammaUL[0] > GammaDL[0] && GammaUR[0] > GammaDR[0] && GammaUR[2] > GammaUL[2] && GammaDR[2] > GammaDL[2] ) )
      throw "NEWTON CORRECTION METHOD FOR CORNER POINTS ERROR! \n";

    // below distances from appropriate sections in appropr. direction (stable for sections to integrate from, unstable for sections to integrate onto)
    interval ruDL(0.015);    // this is a,c in the paper      
    interval rsUL(0.015);   // this is b,d in the paper

    interval ruUR(0.029);  // this is a, c in the paper
    interval rsDR(0.03);   // this is b, d in the paper

    IMatrix PUL( coordChange( *Fhn_vf, GammaUL ) ), 
            PUR( coordChange( *Fhn_vf, GammaUR ) ), 
            PDL( coordChange( *Fhn_vf, GammaDL ) ),  
            PDR( coordChange( *Fhn_vf, GammaDR ) ); 
 
    IVector setToIntegrateDL(2);
    IVector setToIntegrateUR(2);

    setToIntegrateDL[0] = 1.2e-2*interval(-1,1);  // this is ys at downleft corner  
    setToIntegrateDL[1] = 5.0e-3*interval(-1,1);  // this is v at downleft corner

    setToIntegrateUR[0] = 1.9e-2*interval(-1,1);  // this is ys at upright corner
    setToIntegrateUR[1] = 5.0e-3*interval(-1,1);  // this is v at upright corner

    IVector setToBackIntegrateUL(2);
    IVector setToBackIntegrateDR(2);

    setToBackIntegrateUL[0] = 5.0e-3*interval(-1,1);     // this is v at upleft corner
    setToBackIntegrateUL[1] = 1.0e-2*interval(-1,1);         // this is yu at upleft corner

    setToBackIntegrateDR[0] = 5.0e-3*interval(-1,1);     // this is v at downright corner 
    setToBackIntegrateDR[1] = 0.7e-2*interval(-1,1);        // this is yu at downright corner


    // left isolating segments

    IVector ULface( rsUL*interval(-1,1), setToBackIntegrateUL[1], 0. ); 
    IVector DLface( setToIntegrateDL[0], ruDL*interval(-1,1), 0. ); 

    FhnIsolatingSegment ULSegment( *Fhn_vf, GammaUL + IVector( 0., 0., setToBackIntegrateUL[0].leftBound() ), 
        GammaUL + IVector( 0., 0., setToBackIntegrateUL[0].rightBound() ), PUL, ULface, ULface, _cornerSegmentDivCount ); 
    FhnIsolatingSegment DLSegment( *Fhn_vf, GammaDL + IVector( 0., 0., setToIntegrateDL[1].leftBound() ), 
      GammaDL + IVector( 0., 0., setToIntegrateDL[1].rightBound() ), PDL, DLface, DLface, _cornerSegmentDivCount );  


    // right isolating segments

    IVector URface( setToIntegrateUR[0], ruUR*interval(-1,1), 0. );
    IVector DRface( rsDR*interval(-1,1), setToBackIntegrateDR[1], 0. ); 
 
    FhnIsolatingSegment URSegment( *Fhn_vf, GammaUR + IVector( 0., 0., setToIntegrateUR[1].leftBound() ), 
        GammaUR + IVector( 0.,0.,setToIntegrateUR[1].rightBound() ), PUR, URface, URface, _cornerSegmentDivCount );  
    FhnIsolatingSegment DRSegment( *Fhn_vf, GammaDR + IVector( 0., 0., setToBackIntegrateDR[0].leftBound() ), 
        GammaDR + IVector( 0., 0., setToBackIntegrateDR[0].rightBound() ), PDR, DRface, DRface, _cornerSegmentDivCount );  
 

    midPoincareMap *leftMap; // this class implements the Poincare maps described in the paper as pmUL, pmDL onto leftSection
    midPoincareMap *rightMap; // this class implements the Poincare maps described in the paper as pmUR, pmDR onto rightSection

 /*   if( withParams )
    {
      leftMap = new midPoincareMap( parameters, *Fhn_vf_withParams, *Fhn_vf_withParams_rev, PDL, PUL, GammaDL, GammaUL, ruDL, rsUL, -1., _pMapDivCount );
      rightMap = new midPoincareMap( parameters, *Fhn_vf_withParams, *Fhn_vf_withParams_rev, PUR, PDR, GammaUR, GammaDR, ruUR, rsDR, 1., _pMapDivCount );
    }
    else
    {*/
      leftMap = new midPoincareMap( *Fhn_vf, *Fhn_vf_rev, DLSegment, ULSegment, _theta, _eps, -1., _pMapDivCount );
      rightMap = new midPoincareMap( *Fhn_vf, *Fhn_vf_rev, URSegment, DRSegment, _theta, _eps, 1., _pMapDivCount );
   // }

    // covering checks 


    if( _verbose )
      cout << "\n ------------------- LEFT SIDE COVERING CHECKS: -------------------------- \n \n";
    if( !(*leftMap).checkCovering( setToIntegrateDL, setToBackIntegrateUL, _verbose )  ) 
      throw "FAILURE TO CHECK COVERINGS IN THE FAST REGIME (LEFT MAP)! \n";

    if( _verbose )
      cout << "\n ------------------- RIGHT SIDE COVERING CHECKS: -------------------------- \n \n";
    if( !(*rightMap).checkCovering( setToIntegrateUR, setToBackIntegrateDR, _verbose )  )    
      throw "FAILURE TO CHECK COVERINGS IN THE FAST REGIME (RIGHT MAP)! \n";

    delete leftMap;
    delete rightMap;

    IVector ULSegment_entranceVerification( ULSegment.entranceVerification() );
    IVector ULSegment_exitVerification( ULSegment.exitVerification() );
    IVector DLSegment_entranceVerification( DLSegment.entranceVerification() );
    IVector DLSegment_exitVerification( DLSegment.exitVerification() );

    if( _verbose )
    {
      cout << "\n ------------------- UL, DL SEGMENTS ISOLATION: -------------------------- \n \n";

      cout << "Enclosures of scalar product of vector field with entrance faces normals for DL segment: \n \n" << DLSegment_entranceVerification << "\n";  
      cout << "\n --- \n";
      cout << "Enclosures of scalar product of vector field with exit faces normals for DL segment: \n \n" << DLSegment_exitVerification << "\n";

      cout << "\n --- \n";

      cout << "Enclosures of scalar product of vector field with entrance faces normals for UL segment: \n \n" << ULSegment_entranceVerification << "\n";  
      cout << "\n --- \n";
      cout << "Enclosures of scalar product of vector field with exit faces normals for UL segment: \n \n" << ULSegment_exitVerification << "\n";


      cout << "\n --- \n";
      cout << "\n --- \n";
    };

    if( !( ULSegment_entranceVerification[0] < 0. && ULSegment_entranceVerification[1] < 0. && ULSegment_exitVerification[0] > 0. && ULSegment_exitVerification[1] > 0. ) )
      throw "ISOLATION ERROR FOR UL CORNER SEGMENT! \n";
    if( !( DLSegment_entranceVerification[0] < 0. && DLSegment_entranceVerification[1] < 0. && DLSegment_exitVerification[0] > 0. && DLSegment_exitVerification[1] > 0. ) )
      throw "ISOLATION ERROR FOR DL CORNER SEGMENT! \n";


    IVector URSegment_entranceVerification( URSegment.entranceVerification() );
    IVector URSegment_exitVerification( URSegment.exitVerification() );
    IVector DRSegment_entranceVerification( DRSegment.entranceVerification() );
    IVector DRSegment_exitVerification( DRSegment.exitVerification() );

    if( _verbose )
    {
      cout << "\n ------------------- UR, DR SEGMENTS ISOLATION: -------------------------- \n \n";

      cout << "Enclosures of scalar product of vector field with entrance faces normals for UR segment: \n \n" << URSegment_entranceVerification << "\n";  
      cout << "\n --- \n";
      cout << "Enclosures of scalar product of vector field with exit faces normals for UR segment: \n \n" << URSegment_exitVerification << "\n";

      cout << "\n --- \n";

      cout << "Enclosures of scalar product of vector field with entrance faces normals for DR segment: \n \n" << DRSegment_entranceVerification << "\n";  
      cout << "\n --- \n";
      cout << "Enclosures of scalar product of vector field with exit faces normals for DR segment: \n \n" << DRSegment_exitVerification << "\n";


      cout << "\n --- \n";
      cout << "\n --- \n";
    };

    if( !( URSegment_entranceVerification[0] < 0. && URSegment_entranceVerification[1] < 0. && URSegment_exitVerification[0] > 0. && URSegment_exitVerification[1] > 0. ) )
      throw "ISOLATION ERROR FOR UR CORNER SEGMENT! \n";
    if( !( DRSegment_entranceVerification[0] < 0. && DRSegment_entranceVerification[1] < 0. && DRSegment_exitVerification[0] > 0. && DRSegment_exitVerification[1] > 0. ) )
      throw "ISOLATION ERROR FOR DR CORNER SEGMENT! \n";


    if( !( URSegment.segmentEnclosure[0] > DRSegment.segmentEnclosure[0] && ULSegment.segmentEnclosure[0] > DLSegment.segmentEnclosure[0] &&
          URSegment.segmentEnclosure[2] > ULSegment.segmentEnclosure[2] && DRSegment.segmentEnclosure[2] > DLSegment.segmentEnclosure[2]) )
      throw "CORNER SEGMENTS ALIGNMENT ERROR! \n";          // a check on whether corner segments are really up/down to the left/right of each other


    // up down isolating segments

    chainOfSegments UpSegment( *Fhn_vf, ULSegment.GammaRight, URSegment.GammaLeft, PUL, PUR, ULface, URface, _chainSegmentDivCount );
    chainOfSegments DownSegment( *Fhn_vf, DLSegment.GammaRight, DRSegment.GammaLeft, PDL, PDR, DLface, DRface, _chainSegmentDivCount ); 

    if( !( ULSegment.segmentEnclosure[0] > ULSegment.segmentEnclosure[2] ) )
      throw "MISALIGNMENT OF ONE OF THE UPPER SEGMENTS! \n";
    if( !( DLSegment.segmentEnclosure[0] < DLSegment.segmentEnclosure[2] ) )
      throw "MISALIGNMENT OF ONE OF THE LOWER SEGMENTS! \n";      // checks on whether we are above/below u=v plane for upper/lower segments

    IVector UpSegment_entranceAndExitVerification( UpSegment.entranceAndExitVerification( _chainSubsegmentCount ) );
    IVector DownSegment_entranceAndExitVerification( DownSegment.entranceAndExitVerification( _chainSubsegmentCount ) );

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

    // ISOLATION CAN BE ALSO CHECKED NOW FOR EACH SUBSEGMENT IN SEGMENTS.HPP TO THROW AN EXCEPTION QUICKER, SEE COMMENTED LINES IN SEGMENTS.HPP
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
    cout << Message << "EXISTENCE OF PERIODIC ORBIT FOR PARAMETER VALUES THETA=" << _theta << " AND EPS=" << _eps << " NOT VERIFIED! \n";
  }
};


