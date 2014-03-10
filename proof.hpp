
/* -----------------------------------------------------------------------------------------
 * This is a header file to fhn.cpp providing a non-rigorous coordinate change function 
 * to straightened variables and a rigorous proof of verification of periodic orbits
 * in the FitzHugh-Nagumo system for given parameters theta, eps. 
 * * ----------------------------------------------------------------------------------------*/

/* ------------------------------------------------------------------------------------ */
/* ----- VERIFICATION OF EXISTENCE OF PERIODIC ORBITS FOR GIVEN PARAMETER VALUES ------ */
/* ------------------------------------------------------------------------------------ */


void FhnVerifyExistenceOfPeriodicOrbit( interval _theta, interval _eps, bool _verbose = 0, bool withParams = 0, int _pMapDivCount = 20, 
     int _longSubsegmentCount = 100, int _longSegmentDivCount = 80, int _cornerSegmentDivCount = 200 ) 
  // verbose on displays all the interval enclosures for Poincare maps / products of vector fields with normals; other parameters control respectively: 
  // number of subdivisions of sets to integrate (in each dimension), number of subsegments along slow manifolds, number of subdivisions of regular/corner segments
  // for evaluation of the scalar product of vector field with outward pointing normals
{
  try                   // we check negations of all assumptions to throw exceptions, if no exception is thrown existence of the orbit is verified
  {
    Fhn_vf.setParameter("theta",_theta);
    Fhn_vf.setParameter("eps",_eps);

    IVector parameters({ _theta, _eps });

    IVector GammaUL(0.970345591417269, 0., 0.0250442158334208);                                   // some guesses for the corner points which are equilibria
    IVector GammaDL(-0.108412947498862, 0., 0.0250442158334208);                                  // of the fast subsystem for critical parameter v values (third variable)
                                                                                                  // where heteroclinics exist
    IVector GammaUR(0.841746280832201, 0., 0.0988076360184288);                                   // UR up right, DR down right, UL up left, DL down left
    IVector GammaDR(-0.237012258083933, 0., 0.0988076360184288);

    GammaQuad_correct( _theta, GammaUL, GammaDL, GammaUR, GammaDR );                              // we correct the initial guesses by nonrigorous Newtons methods (see numerics.hpp)

    if( !(GammaUL[0] > GammaDL[0] && GammaUR[0] > GammaDR[0] && GammaUR[2] > GammaUL[2] && GammaDR[2] > GammaDL[2] ) )
      throw "NEWTON CORRECTION METHOD FOR CORNER POINTS ERROR! \n";
        
    interval ruDL(0.011);           // distances from appropriate sections in appropr. direction (stable for sections to integrate from, unstable for sections to integrate onto)
    interval rsUL(0.01);         

    interval ruUR(0.0015);
    interval rsDR(0.028);

    IMatrix PUL( coordChange( Fhn_vf, GammaUL ) ), 
            PUR( coordChange( Fhn_vf, GammaUR ) ), 
            PDL( coordChange( Fhn_vf, GammaDL ) ),  
            PDR( coordChange( Fhn_vf, GammaDR ) ); 

    FhnPoincareMap *PMAPL;
    FhnPoincareMap *PMAPR;
    
    midPoincareMap testMap( parameters, Fhn_vf_withParams, Fhn_vf_withParams_rev, PDL, PUL, GammaDL, GammaUL, ruDL, rsUL, -1., 60 );
    


    if( withParams )
    {
      PMAPL = new FhnPoincareMap( Fhn_vf, PDL, PUL, GammaDL, GammaUL, ruDL, rsUL, -1., _pMapDivCount ); // -1 because the exit/entrance sections are aligned in a reversed order
      PMAPR = new FhnPoincareMap( Fhn_vf, PUR, PDR, GammaUR, GammaDR, ruUR, rsDR, 1., _pMapDivCount );
    }
    else
    {
      PMAPL = new FhnPoincareMap( parameters, Fhn_vf_withParams, PDL, PUL, GammaDL, GammaUL, ruDL, rsUL, -1., _pMapDivCount ); 
      PMAPR = new FhnPoincareMap( parameters, Fhn_vf_withParams, PUR, PDR, GammaUR, GammaDR, ruUR, rsDR, 1., _pMapDivCount );
    } 

    IVector setToIntegrateDL(2);
    IVector setToIntegrateUR(2);

    setToIntegrateDL[0] = 1.0e-3*interval(-1,1);  // this is ys at downleft corner  
    setToIntegrateDL[1] = 1.0e-3*interval(-1,1);  // this is v at downleft corner

    setToIntegrateUR[0] = 1.0e-3*interval(-1,1);  // this is ys at upright corner
    setToIntegrateUR[1] = 1.0e-4*interval(-1,1);  // this is v at upright corner

    // sets to integrate backwards - only with the parameter _midsection = 1 on

    IVector setToBackIntegrateUL(2);
    IVector setToBackIntegrateDR(2);

    setToBackIntegrateUL[0] = 0.4*1.0e-3*interval(-1,1);     // this is v at upleft corner
    setToBackIntegrateUL[1] = 1.0e-3*interval(-1,1);     // this is yu at upleft corner

    setToBackIntegrateDR[0] = 1.0e-3*interval(-1,1);      
    setToBackIntegrateDR[1] = 1.0e-4*interval(-1,1);     

    // left corner segments/coverings
  
    cout << testMap.checkCovering( setToIntegrateDL, setToBackIntegrateUL ) << "! \n";

    IVector PMAPL_leftU = (*PMAPL)( leftU(setToIntegrateDL) );
    IVector PMAPL_rightU = (*PMAPL)( rightU(setToIntegrateDL) );
    IVector PMAPL_all = (*PMAPL)( setToIntegrateDL );

    delete PMAPL;

    if( _verbose )
    {
      cout << "\n ----------------------------- LEFT POINCARE MAP: --------------------------------- \n \n";

      cout << "All enclosures in section coordinates! \n" << "Image enclosure of DL segment exit face through left Poincare map: \n \n" << PMAPL_all << "\n \n" 
        << "Image enclosure of left unstable edge of DL segment exit face through left Poincare map: \n \n" << PMAPL_leftU << "\n \n" << 
        "Image enclosure of right unstable edge of DL segment exit face through left Poincare map: \n \n" << PMAPL_rightU << "\n \n";
    }

   
    if( !( PMAPL_leftU[1] + EPS < 0. && PMAPL_rightU[1] - EPS > 0. && PMAPL_all[0].leftBound() < 0. && PMAPL_all[0].rightBound() > 0. ) )
      throw "LEFT POINCARE MAP COVERING ERROR! \n";

    // faces of two isolating segments around slow manifolds - determined by the stable/unstable distances from slow manifolds given above, used also in rigorous integration
    // unstable faces are shortened by a small number EPS - so that there is covering by image of Poincare map
    // zeroes at the third coordinate are just to make some algebra easier (adding to 3d vectors etc.)

    IVector ULface( rsUL*interval(-1,1), interval( (PMAPL_leftU[1] + EPS).rightBound(), (PMAPL_rightU[1] - EPS).leftBound() ), 0. ); 
    IVector DLface( setToIntegrateDL[0], ruDL*interval(-1,1), 0. ); 

    FhnIsolatingSegment ULSegment( Fhn_vf, GammaUL + IVector( 0., 0., PMAPL_all[0].leftBound()-EPS ), 
        GammaUL + IVector( 0., 0., PMAPL_all[0].rightBound()+EPS ), PUL, ULface, ULface, _cornerSegmentDivCount ); // v face is expanded by EPS to get stable face covering from Poincare map
    FhnIsolatingSegment DLSegment( Fhn_vf, GammaDL + IVector( 0., 0., setToIntegrateDL[1].leftBound() ), 
      GammaDL + IVector( 0., 0., setToIntegrateDL[1].rightBound() ), PDL, DLface, DLface, _cornerSegmentDivCount );  // TODO: add EPS?

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


    // right corner segments/coverings

    IVector PMAPR_leftU = (*PMAPR)( leftU(setToIntegrateUR) );
    IVector PMAPR_rightU = (*PMAPR)( rightU(setToIntegrateUR) );
    IVector PMAPR_all = (*PMAPR)( setToIntegrateUR );
    
    delete PMAPR;
  
    if( _verbose )
    {
      cout << "\n ----------------------------- RIGHT POINCARE MAP: --------------------------------- \n \n";

      cout << "All enclosures in section coordinates! \n" << "Image enclosure of UR segment exit face through right Poincare map: \n \n" << PMAPR_all << "\n \n" 
        << "Image enclosure of left unstable edge of UR segment exit face through right Poincare map: \n \n" << PMAPR_leftU << "\n \n" << 
        "Image enclosure of right unstable edge of UR segment exit face through right Poincare map: \n \n" << PMAPR_rightU << "\n \n";
    }
    
    if( !( PMAPR_leftU[1] + EPS < 0. && PMAPR_rightU[1] - EPS > 0. && PMAPR_all[0].leftBound() < 0. && PMAPR_all[0].rightBound() > 0.) )
      throw "RIGHT POINCARE MAP COVERING ERROR! \n";

    IVector URface( setToIntegrateUR[0], ruUR*interval(-1,1), 0. );
    IVector DRface( rsDR*interval(-1,1), interval( (PMAPR_leftU[1] + EPS).rightBound(), (PMAPR_rightU[1] - EPS).leftBound() ), 0. ); 
 
    FhnIsolatingSegment URSegment( Fhn_vf, GammaUR + IVector( 0., 0., setToIntegrateUR[1].leftBound() ), 
        GammaUR + IVector( 0.,0.,setToIntegrateUR[1].rightBound() ), PUR, URface, URface, _cornerSegmentDivCount );  // TODO: add EPS?
    FhnIsolatingSegment DRSegment( Fhn_vf, GammaDR + IVector( 0., 0., PMAPR_all[0].leftBound()-EPS ), 
        GammaDR + IVector( 0., 0., PMAPR_all[0].rightBound()+EPS ), PDR, DRface, DRface, _cornerSegmentDivCount );  // again, v face is expanded by EPS in both directions
 

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


    longIsolatingSegment UpSegment( Fhn_vf, ULSegment.GammaRight, URSegment.GammaLeft, PUL, PUR, ULface, URface, _longSegmentDivCount );
    longIsolatingSegment DownSegment( Fhn_vf, DLSegment.GammaRight, DRSegment.GammaLeft, PDL, PDR, DLface, DRface, _longSegmentDivCount ); 

    if( !( ULSegment.segmentEnclosure[0] > ULSegment.segmentEnclosure[2] && UpSegment.segmentEnclosure[0] > UpSegment.segmentEnclosure[2] && 
          URSegment.segmentEnclosure[0] > URSegment.segmentEnclosure[2] ) )
      throw "MISALIGNMENT OF ONE OF THE UPPER SEGMENTS! \n";
    if( !( DLSegment.segmentEnclosure[0] < DLSegment.segmentEnclosure[2] && DownSegment.segmentEnclosure[0] < DownSegment.segmentEnclosure[2] && 
          DRSegment.segmentEnclosure[0] < DRSegment.segmentEnclosure[2] ) )
      throw "MISALIGNMENT OF ONE OF THE LOWER SEGMENTS! \n";      // checks on whether we are above/below u=v plane for upper/lower segments

    IVector UpSegment_entranceAndExitVerification( UpSegment.entranceAndExitVerification( _longSubsegmentCount ) );
    IVector DownSegment_entranceAndExitVerification( DownSegment.entranceAndExitVerification( _longSubsegmentCount ) );

    if( _verbose )
    {
      cout << "\n ---------------------------- UP, DOWN SEGMENTS ISOLATION: ---------------------------- \n \n";

      cout << "Interval hull of enclosures of scalar products of the vector field with up segments (not including corner ones, left/right entrance faces first, then exit faces): \n \n"
        << UpSegment_entranceAndExitVerification << "\n";  
      cout << "\n --- \n";
      cout << "Interval hull of enclosures of scalar products of the vector field with down segments (not including corner ones, left/right entrance faces first, then exit faces): \n \n"
        << DownSegment_entranceAndExitVerification << "\n \n";  
 
      cout << "\n --- \n";
      cout << "\n --- \n";   
    };

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


