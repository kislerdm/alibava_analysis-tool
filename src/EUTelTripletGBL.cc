// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-

/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#if defined USE_GEAR

// EUTelescope includes:
#include "EUTelTripletGBL.h"

#include "EUTELESCOPE.h"
//#include "EUTelSparseDataImpl.h"
#include "EUTelBaseSparsePixel.h"
#include "EUTelGenericSparsePixel.h"
#include "EUTelTrackerDataInterfacerImpl.h"
#include "EUTelExceptions.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTelPStream.h" // process streams redi::ipstream

// AIDA histogram package (on top of ROOT):

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IAxis.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IProfile1D.h>
#include <AIDA/IProfile2D.h>
#include <AIDA/ITree.h>
#endif

// GBL:
#include "include/GblTrajectory.h"
#include "include/MilleBinary.h"


// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"
#include "marlin/ProcessorMgr.h"
#include "marlin/Global.h"

// gear includes <.h>
#include <gear/GearMgr.h>
#include <gear/SiPlanesParameters.h>

// LCIO includes <.h>
#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <Exceptions.h>

// system includes <>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <memory>
#include <string.h>
#include <map>
#include <cstdlib>
#include <limits>

// ROOT includes ".h"
#include <TMath.h>
#include <TVectorD.h>
#include <TF1.h>
#include <TMatrixD.h>
#include <TVector3.h>
#include <TRotation.h>
#include "TH1D.h"

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace eutelescope;


EUTelTripletGBL::EUTelTripletGBL() : Processor("EUTelTripletGBL"), _siPlanesParameters(), trk_counter(0),
    _siPlanesLayerLayout(), _inputCollectionTelescope(""), _isFirstEvent(0), _eBeam(0), _nEvt(0), _nTelPlanes(0),
    _dut_plane(3), _cutx(0.15), _cuty(0.1), _planeSort(), _planeID(), _planePosition(),
    _planeThickness(), _planeX0(), _planeResolution(), _res(3.44E-3),
    triplet_angle_cut(0.010), triplet_residual_cut(0.1), intersect_residual_cut(0.1) {

  // modify processor description
  _description = "Analysis for Alibava as DUT (PlaneID  = 7), CMSPixRef (PlaneID  = 8) in DATURA";

  // processor parameters
  registerInputCollection( LCIO::TRACKERHIT,
                           "InputHitCollection" ,
                           "Name of the input TrackerHit collection of the telescope",
                           _inputCollectionTelescope,
                           std::string("alignedhits") );

  registerProcessorParameter( "Ebeam",
                              "Beam energy [GeV]",
                  _eBeam, static_cast < double >( 5.0));


  registerProcessorParameter( "matching_cut_x",
                              "cut for matching in x coordinate in mm",
			      _cutx, static_cast < double >(0.15));
  registerProcessorParameter( "matching_cut_y",
                              "cut for matching in y coordinate in mm",
			      _cuty, static_cast < double >(0.10));

  registerProcessorParameter( "slope_cut_x",
                              "cut for track slopes in x coordinate in rad",
			      _slope_x, static_cast < double >(0.002));
  registerProcessorParameter( "slope_cut_y",
                              "cut for track slopes in y coordinate in rad",
			      _slope_y, static_cast < double >(0.002));
  registerProcessorParameter( "dut_plane",
			      "plane to be considered the DUT and excluded from the track fit",
			      _dut_plane, static_cast <int>(3));
  registerProcessorParameter( "telescope intrinsic resolution",
                  "intrinsic resoution of the telescope in mm (depends on the beam energy, plane spacing and the M26 threshold)",
                  _res, static_cast <double>(3.44E-3));
  registerProcessorParameter( "Triplet_matching_distance",
                  "Triplets matching paramter in mm: triplets belong to one track if the diff. between triples'' coordinates on DUT pos.",
  intersect_residual_cut, static_cast <double>(2));
  registerProcessorParameter( "Triplet_angle_cut",
                  "Angle matching paramter in radian: input collection hits are being assigned as a triplet if the angular scattering is not higher than this number",
                              triplet_angle_cut, static_cast <double>(0.010));
  registerProcessorParameter( "Triplet_residual_cut",
                  "Triplet residual matching paramter in mm: input collection hits are being assignes as a triplet if the residual between the middle plane's hit and the estimation from the line drawn through the surrounding planes' hits is not higher than this number",
                              triplet_residual_cut, static_cast <double>(0.1));

  
}


void EUTelTripletGBL::init() {

  // usually a good idea to
  printParameters();

  _nEvt = 0;
  trk_counter = 0;
  _isFirstEvent = true;

  streamlog_out(MESSAGE0) << "Beam energy " << _eBeam << " GeV" <<  std::endl;

  // Read geometry information from GEAR

  streamlog_out( MESSAGE0 ) << "Reading telescope geometry description from GEAR " << std::endl;

  _siPlanesParameters  = const_cast<gear::SiPlanesParameters* >( &(Global::GEAR->getSiPlanesParameters()));
  _siPlanesLayerLayout = const_cast<gear::SiPlanesLayerLayout*>(  &(_siPlanesParameters->getSiPlanesLayerLayout() ));


  // Take all layers defined in GEAR geometry
  _nTelPlanes = _siPlanesLayerLayout->getNLayers();

  _planeSort = new int[_nTelPlanes];
  _planePosition = new double[_nTelPlanes]; // z pos
  _planeID         = new int[_nTelPlanes];
  _planeThickness  = new double[_nTelPlanes];
  _planeX0         = new double[_nTelPlanes];
  _planeResolution = new double[_nTelPlanes];

  for(int i = 0; i < _nTelPlanes; i++) {
    
    _planeID[i]=_siPlanesLayerLayout->getID(i);
    _planePosition[i]=_siPlanesLayerLayout->getLayerPositionZ(i);
    _planeThickness[i]=_siPlanesLayerLayout->getLayerThickness(i);
    _planeX0[i]=_siPlanesLayerLayout->getLayerRadLength(i);
    _planeResolution[i] = _siPlanesLayerLayout->getSensitiveResolution(i);
//  streamlog_out(MESSAGE4) << "i = " << i << " PlaneID = " << _planeID[i]
//                             << "z = " << _planePosition[i] << endl;
  }

  streamlog_out( MESSAGE2 ) <<  "Telescope configuration with " << _nTelPlanes << " planes" << std::endl;

//  for( int ipl = 0; ipl < _nTelPlanes; ipl++) {
//    streamlog_out( MESSAGE4 ) << "IDz " << ipl << "  ID = " << _planeID[ipl]
//       << "  at Z [mm] = " << _planePosition[ipl]
//       << " dZ [um] = " << _planeThickness[ipl]*1000.
//       << "  Res [um] = " << _planeResolution[ipl]*1000. << endl;
//  }

  // Book histograms:
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
  bookHistos();
#endif
}//init

//------------------------------------------------------------------------------
void EUTelTripletGBL::processRunHeader( LCRunHeader* runHeader) {

  std::auto_ptr<EUTelRunHeaderImpl> eutelHeader( new EUTelRunHeaderImpl( runHeader ) );
  eutelHeader->addProcessor( type() );

  // Decode and print out Run Header information - just a check

  _nRun = runHeader->getRunNumber();

  streamlog_out( MESSAGE2 )  << "Processing run header"
                             << ", run nr " << runHeader->getRunNumber() << std::endl;

  const std::string detectorName = runHeader->getDetectorName();
  const std::string detectorDescription = runHeader->getDescription();

  streamlog_out( MESSAGE0 ) << detectorName << " : " << detectorDescription << std::endl;

} // processRunHeader

//----------------------------------------------------------------------------
void EUTelTripletGBL::processEvent( LCEvent * event ) {

  EUTelEventImpl * euEvent = static_cast<EUTelEventImpl*>( event );

  if( euEvent->getEventType() == kEORE ) {
    streamlog_out( DEBUG5 ) <<  "EORE found: nothing else to do." << std::endl;
    return;
  }

  int runNumber = event->getRunNumber();

  if( _isFirstEvent ) {

    // apply all GEAR/alignment offsets to get corrected X,Y,Z
    // position of the sensor center

    for( int iplane = 0; iplane < _siPlanesLayerLayout->getNLayers(); iplane++ ) {

      std::map< unsigned int , double > _planeCenter;
      std::map< unsigned int , double > _planeNormal;

      // start filling the map with Gear values:
      _planeCenter[ 0 ] =  _siPlanesLayerLayout->getLayerPositionX(iplane); // X
      _planeCenter[ 1 ] =  _siPlanesLayerLayout->getLayerPositionY(iplane); // Y
      _planeCenter[ 2 ] =  _siPlanesLayerLayout->getLayerPositionZ(iplane); // Z
      _planeNormal[ 0 ] =  0.; // X
      _planeNormal[ 1 ] =  0.; // Y
      _planeNormal[ 2 ] =  1.; // Z

      TVector3  _normalTVec( _planeNormal[0], _planeNormal[1], _planeNormal[2]); 

      // do initial rotation from GEAR
      try{
 	double gRotation[3] = { 0., 0., 0.}; // not rotated
                         
	gRotation[0] = _siPlanesLayerLayout->getLayerRotationXY(iplane); // Euler alpha
	gRotation[1] = _siPlanesLayerLayout->getLayerRotationZX(iplane); // Euler alpha
	gRotation[2] = _siPlanesLayerLayout->getLayerRotationZY(iplane); // Euler alpha
                          
	// input angles are in degree, translate into radian
	gRotation[0] =  gRotation[0]*3.1415926/180.;
	gRotation[1] =  gRotation[1]*3.1415926/180.;
	gRotation[2] =  gRotation[2]*3.1415926/180.;

	TRotation r;
	r.RotateX( gRotation[2] );
	r.RotateY( gRotation[1] );
	r.RotateZ( gRotation[0] );
	_normalTVec.Transform( r );
                                  
	_planeNormal[0] = _normalTVec[0];
	_planeNormal[1] = _normalTVec[1];
	_planeNormal[2] = _normalTVec[2];
      }
      catch(...) {
	streamlog_out(DEBUG5) << "no sensor rotation is given in the GEAR steering file, assume NONE" << std::endl;
      }
    }
    if( isFirstEvent() ) _isFirstEvent = false;

  } // FirstEvent
  
  // Increase event count
  _nEvt++;
     
  //----------------------------------------------------------------------------
  // check input collection (aligned hits):

  LCCollection *collection;
  try {
    collection = event->getCollection( _inputCollectionTelescope );
  }
  catch( lcio::DataNotAvailableException& e) {
    streamlog_out( DEBUG1 ) << "Not able to get collections "
			    << _inputCollectionTelescope << " "
			    << "\nfrom event " << event->getEventNumber()
			    << " in run " << runNumber  << std::endl;

    return;
  }


  //----------------------------------------------------------------------------
  // Copy hits to local table
  // Assign hits to sensor planes

  streamlog_out( DEBUG4 )  << "Total of " << collection->getNumberOfElements() << " tracker hits in input collection " << std::endl;
  nAllHitHisto->fill(collection->getNumberOfElements());

  //----------------------------------------------------------------------------

  std::vector<hit> * hits = new std::vector<hit>;

  // Extract all hits from the LCIO collection and calculate the
  // position uncertainty if necessary:
  for( int ihit = 0; ihit < collection->getNumberOfElements(); ihit++ ) {

    hit newhit;
    TrackerHit * meshit = dynamic_cast<TrackerHit*>( collection->getElementAt(ihit) );
    const double * pos = meshit->getPosition();
    const EVENT::FloatVec cov = meshit->getCovMatrix();

    // Write the position:
    newhit.x = pos[0];
    newhit.y = pos[1];
    newhit.z = pos[2];

    // Find Plane ID to which the hit belongs by minimizing its
    // distance in Z
    //FIXME to be fixed with more general geometry description!
    double distMin = 99; // [mm]

    bool foundplane = false;
    for( int ipl = 0; ipl < _nTelPlanes; ipl++ ) {

        if( abs(newhit.z - _planePosition[ipl]) < distMin ) {
//        if(_planeID[ipl] == 7) { newhit.plane = 6;} //alibava (DUT)
//        if(_planeID[ipl] == 8) { newhit.plane = 7;} //CMSPixRef\
//        if(_planeID[ipl] < 6) { newhit.plane = _planeID[ipl]; }
        newhit.plane = ipl;
        foundplane = true;
        distMin = abs(newhit.z - _planePosition[ipl]);
      }
        streamlog_out(DEBUG0) << "counter = " << ipl << " HitID = " << newhit.plane
                                    << " GearID = " << _planeID[ipl]
                                    << " Z = " << _planePosition[ipl] << " nexZ = " << newhit.z  << endl;
    }

    if(!foundplane) {
      streamlog_out(ERROR5) << "Could not associate hit with a telescope plane. Skipping event." << std::endl;
      return;
    }

    // ...and uncertainty:
    if( cov.at(0) > 0. ) newhit.ex = sqrt(cov.at(0));
    else newhit.ex = _planeResolution[newhit.plane];

    if( cov.at(2) > 0. ) newhit.ey = sqrt(cov.at(2));
    else newhit.ey = _planeResolution[newhit.plane];

    streamlog_out(DEBUG3) << "Hit " << ihit << ": " << newhit << std::endl;

    // Add it to the vector of telescope hits:
    hits->push_back(newhit);

  } // loop over hits

  streamlog_out(DEBUG4) << "Event " << event->getEventNumber()
			<< " contains " << hits->size() << " hits" << std::endl;


  // Fill the telescope plane correlation plots:
  TelescopeCorrelationPlots(hits);


  // -----------------------------------------------------------------
  // Downstream Telescope Triplets ("driplets")

  // Generate new triplet set for the Telescope Downstream Arm:
  std::vector<triplet> * downstream_triplets = new std::vector<triplet>;
  downstream_triplets = FindTriplets(hits, 3, 4, 5);

  // Iterate over all found downstream triplets to fill histograms and match them to the REF and DUT:
  for( std::vector<triplet>::iterator drip = downstream_triplets->begin(); drip != downstream_triplets->end(); drip++ ){

    // Fill some histograms for downstream triplets:
    dridxHisto->fill( (*drip).getdx(4)*1E3 ); //sig = 7 um at 5 GeV
    dridxvsx->fill( (*drip).base().x, (*drip).getdx(4)*1E3 ); // check for rot
    dridxvsy->fill( (*drip).base().y, (*drip).getdx(4)*1E3 );
    dridxvstx->fill( (*drip).slope().x*1E3, (*drip).getdx(4)*1E3 ); // check for z shift
    dridxvsty->fill( (*drip).slope().y*1E3, (*drip).getdx(4)*1E3 );

    dridyHisto->fill( (*drip).getdy(4)*1E3 );
    dridyvsx->fill( (*drip).base().x, (*drip).getdy(4)*1E3 );
    dridyvsy->fill( (*drip).base().y, (*drip).getdy(4)*1E3 );
    dridyvstx->fill( (*drip).slope().x*1E3, (*drip).getdy(4)*1E3 );
    dridyvsty->fill( (*drip).slope().y*1E3, (*drip).getdy(4)*1E3 );

    drixHisto->fill( -(*drip).gethit(4).x ); // -x = x_DP = out
    driyHisto->fill( -(*drip).gethit(4).y ); // -y = y_DP =  up
    drixyHisto->fill( -(*drip).gethit(4).x, -(*drip).gethit(4).y );
    dritxHisto->fill( (*drip).slope().x*1E3 );
    drityHisto->fill( (*drip).slope().y*1E3 );

  }//loop over driplets

  ndriHisto->fill( downstream_triplets->size() );


  //##########################  0-1-2 upstream triplets #######################
  // here again a telescope hit triplet is formed, from planes 0 and 2
  // and the correlation to plane 1. the found triplets are
  // extrapolated and matched to the DUT plane and exhaustive
  // histograms are written.
  // Furthermore the hits of a triplet matched to the DUT are
  // collected and fed into GBL for a track fit.

  // Generate new triplet set for the Telescope Downstream Arm:
  std::vector<triplet> * upstream_triplets = new std::vector<triplet>;
  upstream_triplets = FindTriplets(hits, 0, 1, 2);

  // Iterate over all found upstream triplets to fill histograms and match them to the REF and DUT:
  for( std::vector<triplet>::iterator trip = upstream_triplets->begin(); trip != upstream_triplets->end(); trip++ ) {

    // Fill some histograms for the upstream triplets:
    tridxHisto->fill( (*trip).getdx(1)*1E3 );
    tridyHisto->fill( (*trip).getdy(1)*1E3 );
    tridx1Histo->fill( (*trip).getdx(1)*1E0 );
    tridy1Histo->fill( (*trip).getdy(1)*1E0 );
    tridxvsx->fill( (*trip).base().x, (*trip).getdx(1)*1E3 ); // check for rot
    tridxvsy->fill( (*trip).base().y, (*trip).getdx(1)*1E3 );
    tridxvstx->fill( (*trip).slope().x*1E3, (*trip).getdx(1)*1E3 ); // check for z shift
    tridxvsty->fill( (*trip).slope().y*1E3, (*trip).getdx(1)*1E3 );
    tridyvsx->fill( (*trip).base().x, (*trip).getdy(1)*1E3 );
    tridyvsy->fill( (*trip).base().y, (*trip).getdy(1)*1E3 );
    tridyvstx->fill( (*trip).slope().x*1E3, (*trip).getdy(1)*1E3 );
    tridyvsty->fill( (*trip).slope().y*1E3, (*trip).getdy(1)*1E3 );
    trixHisto->fill( -(*trip).gethit(1).x );
    triyHisto->fill( -(*trip).gethit(1).y );
    trixyHisto->fill( -(*trip).gethit(1).x, -(*trip).gethit(1).y );
    tritxHisto->fill( (*trip).slope().x*1E3 );
    trityHisto->fill( (*trip).slope().y*1E3 );


    // Extrapolate Upstream triplet to Downstream planes 3,4,5: Resolution studies
    for( std::vector<hit>::iterator lhit = hits->begin(); lhit != hits->end(); lhit++ ){

      if( (*lhit).plane <= 2 ) continue; // want 3,4, or 5

      // Fill residuals of triplet and hit in the selected plane:
      if( (*lhit).plane == 3 ) {
	tridx3Histo->fill( (*trip).getdx((*lhit))*1E0 );
	tridy3Histo->fill( (*trip).getdy((*lhit))*1E0 ); // 65 um at 4.7 GeV with CMS
	tridx3bHisto->fill( (*trip).getdx((*lhit))*1E3 ); // finer binning
	tridy3bHisto->fill( (*trip).getdy((*lhit))*1E3 ); // 
      }
      else if( (*lhit).plane == 4 ) {
	tridx4Histo->fill( (*trip).getdx((*lhit))*1E0 );
	tridy4Histo->fill( (*trip).getdy((*lhit))*1E0 ); //174 um at 4.7 GeV
	tridx4bHisto->fill( (*trip).getdx((*lhit))*1E3 ); // finer binning
	tridy4bHisto->fill( (*trip).getdy((*lhit))*1E3 ); // 
      }
      else if( (*lhit).plane == 5 ) {
	tridx5Histo->fill( (*trip).getdx((*lhit))*1E0 );
	tridy5Histo->fill( (*trip).getdy((*lhit))*1E0 ); //273 um at 4.7 GeV
	tridx5bHisto->fill( (*trip).getdx((*lhit))*1E3 ); // finer binning
	tridy5bHisto->fill( (*trip).getdy((*lhit))*1E3 ); // 
      }
    }// Resolution studies

  }// iterate over upstream triplets

  ntriHisto->fill( upstream_triplets->size() );


  //----------------------------------------------------------------------------
  // six: triplets A and driplets B
  // matching and GBL fit
  // kinks: triplets A vs driplets B
  // scattering point = DUT:

  //FIXME: DUTz is the place where they should be matched!
  // now hardcoded to center of telescope...
  double DUTz = _planePosition[6]; //DUT position

  // Match the Telescope Upstream and Downstream Arm triplets to get tracks:
  std::vector<track> * telescope_tracks = new std::vector<track>;
  telescope_tracks = MatchTriplets(upstream_triplets,downstream_triplets,DUTz);

  int ngbl = 0;
  for( std::vector<track>::iterator tr = telescope_tracks->begin(); tr != telescope_tracks->end(); tr++ ){

    triplet trip = (*tr).get_upstream();
    triplet drip = (*tr).get_downstream();

    double probchi = 0;

    // Track kinks as difference in triplet slopes:
    //      double kx = (*drip).slope().x - (*trip).slope().x; //kink
    //      double ky = (*drip).slope().y - (*trip).slope().y;
    double kx = (*tr).kink_x();
    double ky = (*tr).kink_y();

    // Track impact position at DUT from Downstream:
    double xB = drip.getx_at(DUTz);
    double yB = drip.gety_at(DUTz);

    // Track impact position at DUT from Upstream:
    double xA = trip.getx_at(DUTz);
    double yA = trip.gety_at(DUTz);

    double dx = xB - xA; // driplet - triplet
    double dy = yB - yA;

    // GBL with triplet A as seed:

    std::vector<gbl::GblPoint> traj_points;

    // build up trajectory:
    std::vector<double> sPoint;

    // plane 0:
    double s = 0;

    TMatrixD proL2m(2,2);
    proL2m.UnitMatrix();

    TVectorD meas(2);

    TVectorD measPrec(2); // precision = 1/resolution^2
    measPrec[0] = 1.0 / _res / _res;
    measPrec[1] = 1.0 / _res / _res;

    // scatter:
    TVectorD scat(2);
    scat.Zero(); //mean is zero

    double p = _eBeam; // beam momentum
    double X0Si = 65e-3 / 94; // Si + Kapton
    double tetSi = 0.0136 * sqrt(X0Si) / p * ( 1 + 0.038*std::log(X0Si) );

    TVectorD wscatSi(2);
    wscatSi[0] = 1.0 / ( tetSi * tetSi ); //weight
    wscatSi[1] = 1.0 / ( tetSi * tetSi );

    TMatrixD alDer( 2, 3 ); // alignment derivatives
    alDer[0][0] = 1.0; // dx/dx
    alDer[1][0] = 0.0; // dy/dx
    alDer[0][1] = 0.0; // dx/dy
    alDer[1][1] = 1.0; // dy/dy

    std::vector<unsigned int> ilab; // 0-5 = telescope

    // plane 0-5:

    double rx[6];
    double ry[6];
    double zprev = _planePosition[0];

    int DUT_label;
    for( int ipl = 0; ipl < 6; ++ipl ){

      // Get the corresponding hit from the track:
      hit trackhit = (*tr).gethit(ipl);

      double step = _planePosition[ipl] - zprev;
      zprev = _planePosition[ipl];

      gbl::GblPoint * point = new gbl::GblPoint( JacobianPointToPoint( step ) );
      s += step;

      double dz = trackhit.z - trip.base().z;
      double xs = trip.base().x + trip.slope().x * dz; // Ax at plane
      double ys = trip.base().y + trip.slope().y * dz; // Ay at plane

      rx[ipl] = trackhit.x - xs;
      ry[ipl] = trackhit.y - ys;

      if( ipl == 0 ) {
	sixx0Histo->fill( -trackhit.x );
	sixy0Histo->fill( -trackhit.y );
      }
      if( ipl == 1 ) {
	sixx1Histo->fill( -trackhit.x );
	sixy1Histo->fill( -trackhit.y );
      }
      if( ipl == 2 ) {
	sixx2Histo->fill( -trackhit.x );
	sixy2Histo->fill( -trackhit.y );
      }
      if( ipl == 3 ) {
	sixx3Histo->fill( -trackhit.x );
	sixy3Histo->fill( -trackhit.y );
      }
      if( ipl == 4 ) {
	sixx4Histo->fill( -trackhit.x );
	sixy4Histo->fill( -trackhit.y );
      }
      if( ipl == 5 ) {
	sixx5Histo->fill( -trackhit.x );
	sixy5Histo->fill( -trackhit.y );
      }

      meas[0] = rx[ipl];
      meas[1] = ry[ipl];

      if(ipl != _dut_plane) { point->addMeasurement( proL2m, meas, measPrec ); }

      point->addScatterer( scat, wscatSi );

      std::vector<int> globalLabels(3);
      globalLabels[0] = 10 + ipl; // dx
      globalLabels[1] = 20 + ipl; // dy
      globalLabels[2] = 40 + ipl; // rot
      alDer[0][2] = -ys; // dx/rot
      alDer[1][2] =  xs; // dy/rot

      traj_points.push_back(*point);
      sPoint.push_back( s );
      DUT_label = sPoint.size();
      delete point;

    } // loop over planes

    // monitor what we put into GBL:
    selxHisto->fill( -xA ); // triplet at DUT
    selyHisto->fill( -yA );
    selaxHisto->fill( trip.slope().x*1E3 );
    selayHisto->fill( trip.slope().y*1E3 );
    seldxHisto->fill( dx*1E3 ); // triplet-driplet match
    seldyHisto->fill( dy*1E3 );
    selkxHisto->fill( kx*1E3 ); // triplet-driplet kink
    selkyHisto->fill( ky*1E3 );

    seldx1Histo->fill( rx[1]*1E3 ); // triplet interpol
    seldy1Histo->fill( ry[1]*1E3 );
    seldx3Histo->fill( rx[3]*1E3 ); // triplet extrapol
    seldy3Histo->fill( ry[3]*1E3 );
    seldx4Histo->fill( rx[4]*1E3 );
    seldy4Histo->fill( ry[4]*1E3 );
    seldx5Histo->fill( rx[5]*1E3 );
    seldy5Histo->fill( ry[5]*1E3 );

    double Chi2;
    int Ndf;
    double lostWeight;

    gbl::GblTrajectory traj(traj_points, false ); // curvature = false
    traj.fit( Chi2, Ndf, lostWeight );
    traj.getLabels(ilab);

    ngbl++;

    gblndfHisto->fill( Ndf );
    if( Ndf == 8 ) 
      gblchi2aHisto->fill( Chi2 );
    else
      gblchi2bHisto->fill( Chi2 );
    probchi = TMath::Prob( Chi2, Ndf );
    gblprbHisto->fill( probchi );

    // bad fits:

    if( probchi < 0.01 ) {

      badxHisto->fill( -xA ); // triplet at DUT
      badyHisto->fill( -yA );
      badaxHisto->fill( trip.slope().x*1E3 );
      badayHisto->fill( trip.slope().y*1E3 );
      baddxHisto->fill( dx*1E3 ); // triplet-driplet match
      baddyHisto->fill( dy*1E3 );
      badkxHisto->fill( kx*1E3 ); // triplet-driplet kink
      badkyHisto->fill( ky*1E3 );

      baddx1Histo->fill( rx[1]*1E3 ); // triplet interpol
      baddy1Histo->fill( ry[1]*1E3 );
      baddx3Histo->fill( rx[3]*1E3 ); // triplet extrapol
      baddy3Histo->fill( ry[3]*1E3 );
      baddx4Histo->fill( rx[4]*1E3 );
      baddy4Histo->fill( ry[4]*1E3 );
      baddx5Histo->fill( rx[5]*1E3 );
      baddy5Histo->fill( ry[5]*1E3 );
    }// bad fit
    else {

      goodx1Histo->fill( rx[1]*1E3 ); // triplet interpol
      goody1Histo->fill( ry[1]*1E3 );
    } // OK fit

    // look at fit:

    TVectorD aCorrection(5);
    TMatrixDSym aCovariance(5);

    double ax[8];
    // double ay[8];
    unsigned int k = 0;


    unsigned int ndim = 2;
    TVectorD aResiduals(ndim);
    TVectorD aMeasErrors(ndim);
    TVectorD aResErrors(ndim);
    TVectorD aDownWeights(ndim);


    TVectorD aKinks(ndim);
    TVectorD aKinkErrors(ndim);
    TVectorD kResErrors(ndim);
    TVectorD kDownWeights(ndim);

    //track = q/p, x', y', x, y
    //        0,   1,  2,  3, 4
    // at plane DUT:
    unsigned int ipos = ilab[1];
    traj.getResults( ipos, aCorrection, aCovariance );
    traj.getMeasResults( ipos, ndim, aResiduals, aMeasErrors, aResErrors, aDownWeights );
    traj.getScatResults( ipos, ndim, aKinks, aKinkErrors, kResErrors, kDownWeights );
    gblax0Histo->fill( aCorrection[1]*1E3 ); // angle x [mrad]
    gbldx0Histo->fill( aCorrection[3]*1E3 ); // shift x [um]
    gblrx0Histo->fill( ( rx[0] - aCorrection[3] ) * 1E3 ); // residual x [um]
    gblpx0Histo->fill( aResiduals[0] / aResErrors[0] ); // pull
    gblqx0Histo->fill( aKinks[0]*1E3 ); // kink
    ax[k] = aCorrection[1]; // angle correction at plane, for kinks
    // ay[k] = aCorrection[2]; // angle correction at plane, for kinks
    k++;

    ipos = ilab[1];
    traj.getResults( ipos, aCorrection, aCovariance );
    traj.getMeasResults(static_cast<unsigned int>(ipos), ndim, aResiduals, aMeasErrors, aResErrors, aDownWeights );
    traj.getScatResults(static_cast<unsigned int>(ipos), ndim, aKinks, aKinkErrors, kResErrors, kDownWeights );
    gblax1Histo->fill( aCorrection[1]*1E3 ); // angle x [mrad]
    gbldx1Histo->fill( aCorrection[3]*1E3 ); // shift x [um]
    gblrx1Histo->fill( ( rx[1] - aCorrection[3] ) * 1E3 ); // residual x [um]
    gblpx1Histo->fill( aResiduals[0] / aResErrors[0] ); // pull
    gblqx1Histo->fill( aKinks[0]*1E3 ); // kink
    gblsx1Histo->fill( aKinks[0]/aKinkErrors[0] ); // x kink pull
    gbltx1Histo->fill( aKinks[0]/kResErrors[0] ); // x kink pull
    ax[k] = aCorrection[1]; // angle correction at plane, for kinks
    //ay[k] = aCorrection[2]; // angle correction at plane, for kinks
    k++;

    ipos = ilab[2];
    traj.getResults( ipos, aCorrection, aCovariance );
    traj.getMeasResults(static_cast<unsigned int>(ipos), ndim, aResiduals, aMeasErrors, aResErrors, aDownWeights );
    traj.getScatResults(static_cast<unsigned int>(ipos), ndim, aKinks, aKinkErrors, kResErrors, kDownWeights );
    gblax2Histo->fill( aCorrection[1]*1E3 ); // angle x [mrad]
    gbldx2Histo->fill( aCorrection[3]*1E3 ); // shift x [um]
    gblrx2Histo->fill( ( rx[2] - aCorrection[3] ) * 1E3 ); // residual x [um]
    gblpx2Histo->fill( aResiduals[0] / aResErrors[0] ); // pull
    gblqx2Histo->fill( aKinks[0]*1E3 ); // kink
    ax[k] = aCorrection[1]; // angle correction at plane, for kinks
    //ay[k] = aCorrection[2]; // angle correction at plane, for kinks
    k++;

    ipos = ilab[3];
    traj.getResults( ipos, aCorrection, aCovariance );
    traj.getMeasResults(static_cast<unsigned int>(ipos), ndim, aResiduals, aMeasErrors, aResErrors, aDownWeights );
    traj.getScatResults(static_cast<unsigned int>(ipos), ndim, aKinks, aKinkErrors, kResErrors, kDownWeights );
    gblax3Histo->fill( aCorrection[1]*1E3 ); // angle x [mrad]
    gbldx3Histo->fill( aCorrection[3]*1E3 ); // shift x [um]
    gblrx3Histo->fill( ( rx[3] - aCorrection[3] ) * 1E3 ); // residual x [um]
    gblpx3Histo->fill( aResiduals[0] / aResErrors[0] ); // pull
    gblqx3Histo->fill( aKinks[0]*1E3 ); // kink
    ax[k] = aCorrection[1]; // angle correction at plane, for kinks
    //ay[k] = aCorrection[2]; // angle correction at plane, for kinks
    k++;

    ipos = ilab[4];
    traj.getResults( ipos, aCorrection, aCovariance );
    traj.getMeasResults(static_cast<unsigned int>(ipos), ndim, aResiduals, aMeasErrors, aResErrors, aDownWeights );
    traj.getScatResults(static_cast<unsigned int>(ipos), ndim, aKinks, aKinkErrors, kResErrors, kDownWeights );
    gblax4Histo->fill( aCorrection[1]*1E3 ); // angle x [mrad]
    gbldx4Histo->fill( aCorrection[3]*1E3 ); // shift x [um]
    gblrx4Histo->fill( ( rx[4] - aCorrection[3] ) * 1E3 ); // residual x [um]
    gblpx4Histo->fill( aResiduals[0] / aResErrors[0] ); // pull
    gblqx4Histo->fill( aKinks[0]*1E3 ); // kink
    ax[k] = aCorrection[1]; // angle correction at plane, for kinks
    //ay[k] = aCorrection[2]; // angle correction at plane, for kinks
    k++;

    ipos = ilab[5];
    traj.getResults( ipos, aCorrection, aCovariance );
    traj.getMeasResults(static_cast<unsigned int>(ipos), ndim, aResiduals, aMeasErrors, aResErrors, aDownWeights );
    traj.getScatResults(static_cast<unsigned int>(ipos), ndim, aKinks, aKinkErrors, kResErrors, kDownWeights );
    gblax5Histo->fill( aCorrection[1]*1E3 ); // angle x [mrad]
    gbldx5Histo->fill( aCorrection[3]*1E3 ); // shift x [um]
    gblrx5Histo->fill( ( rx[5] - aCorrection[3] ) * 1E3 ); // residual x [um]
    gblpx5Histo->fill( aResiduals[0] / aResErrors[0] ); // pull
    gblqx5Histo->fill( aKinks[0]*1E3 ); // kink
    ax[k] = aCorrection[1]; // angle correction at plane, for kinks
    //ay[k] = aCorrection[2]; // angle correction at plane, for kinks
    k++;

    // kinks: 1,2 = tele, 3 = DUT, 4,5 = tele

    gblkx1Histo->fill( (ax[1] - ax[0])*1E3 ); // kink at 1 [mrad]
    gblkx2Histo->fill( (ax[2] - ax[1])*1E3 ); // kink at 2 [mrad]
    gblkx3Histo->fill( (ax[3] - ax[2])*1E3 ); // kink at 3 [mrad]
    gblkx4Histo->fill( (ax[4] - ax[3])*1E3 ); // kink at 4 [mrad]
    gblkx5Histo->fill( (ax[5] - ax[4])*1E3 ); // kink at 5 [mrad]
    gblkx6Histo->fill( (ax[6] - ax[5])*1E3 ); // kink at 6 [mrad]


    //------------------------------------------------------------------------
    // intersect point in z:
    hit intersect = (*tr).intersect();
    double zx = intersect.x;
    double zy = intersect.y;

    if( abs(dy) < 0.1 ) { // no cut on dx
      if( abs( kx ) > 0.003 ) { // 
	sixzx3Histo->fill( zx - _planePosition[2] );
      }
      if( abs( kx ) > 0.002 ) { // 
	sixzx2Histo->fill( zx - _planePosition[2] );
      }
    }

    if( abs(dx) < 0.1 ) { // no cut on dy
      if( abs( ky ) > 0.003 ) { // 
	sixzy3Histo->fill( zy - _planePosition[2] );
      }
      if( abs( ky ) > 0.002 ) { // 
	sixzy2Histo->fill( zy - _planePosition[2] );
      }
    }

    //------------------------------------------------------------------------
    // z intersect:

    if( abs(dx) < 0.2 && abs(dy) < 0.2 ) { // looser cut allows more z range

      // measure scattering angle x after cuts in y:
      // cut on ky creates bias in kx

      if( abs( kx ) > 0.001 ) {
	sixzx1Histo->fill( zx - _planePosition[2] );
	if( abs( zx - DUTz ) < 30 ) {
	  sixkyzxHisto->fill( ky*1E3 );
	  sixkxzxHisto->fill( kx*1E3 ); // plot with gap, fittp0g.C("sixkxzx")
	}
      }

      if( abs( ky ) > 0.001 ) {
	sixzy1Histo->fill( zy - _planePosition[2] );
	if( abs( zy - DUTz ) < 30 ) {
	  sixkxzyHisto->fill( kx*1E3 );
	  sixkyzyHisto->fill( ky*1E3 ); // plot with gap
	}
      }

    }//match

     //------------------------------------------------------------------------

  } // Loop over found tracks

  nsixHisto->fill( telescope_tracks->size() );

  //prevdutrefddt = dutrefddt;
  // Clear memory
  delete telescope_tracks;
  delete downstream_triplets;
  delete upstream_triplets;
  delete hits;
}


//------------------------------------------------------------------------------
void EUTelTripletGBL::check( LCEvent * /* evt */  ) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


//------------------------------------------------------------------------------
void EUTelTripletGBL::end(){

  // Print the summary:
  streamlog_out(MESSAGE5) << "Processed events: " << _nEvt << " | Tracks found: " << trk_counter << endl
                          << "Triplets matching distance = " << intersect_residual_cut << " [mm]" << endl
                             << "Hits matching distance (triplet building) = " << triplet_residual_cut << " [mm]" << endl
                                << "Hits matching angle (triplet building) = " << triplet_angle_cut << " [rad.]" << endl
                          << endl;

} // end end


//------------------------------------------------------------------------------
void EUTelTripletGBL::bookHistos()
{

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

  // telescope hits per plane:
  AIDAProcessor::tree(this)->mkdir("Telescope");

  nAllHitHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/nallhit", 201, -0.5, 200.5 );
  nAllHitHisto->setTitle( "Telescope hits/event;telescope hits;events" );

  // telescope dx:

  dx01Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dx01", 100, -1, 1 );
  dx01Histo->setTitle( "x1-x0;x_{1}-x_{0} [mm];hit pairs" );

  dy01Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dy01", 100, -1, 1 );
  dy01Histo->setTitle( "y1-y0;y_{1}-y_{0} [mm];hit pairs" );

  du01Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/du01", 100, -1, 1 );
  du01Histo->setTitle( "x1-x0, |dy| < 1;x_{1}-x_{0} [mm];hit pairs" );

  dx02Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dx02", 100, -1, 1 );
  dx02Histo->setTitle( "x2-x0;x_{2}-x_{0} [mm];hit pairs" );

  dx03Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dx03", 100, -1, 1 );
  dx03Histo->setTitle( "x3-x0;x_{3}-x_{0} [mm];hit pairs" );

  dx04Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dx04", 100, -1, 1 );
  dx04Histo->setTitle( "x4-x0;x_{4}-x_{0} [mm];hit pairs" );

  dx05Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dx05", 100, -1, 1 );
  dx05Histo->setTitle( "x5-x0;x_{5}-x_{0} [mm];hit pairs" );

  dx12Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dx12", 100, -1, 1 );
  dx12Histo->setTitle( "x2-x1;x_{2}-x_{1} [mm];hit pairs" );

  dy12Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dy12", 100, -1, 1 );
  dy12Histo->setTitle( "y2-y1;y_{2}-y_{1} [mm];hit pairs" );

  du12Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/du12", 100, -1, 1 );
  du12Histo->setTitle( "x2-x1, |dy| < 1;x_{2}-x_{1} [mm];hit pairs" );

  dx23Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dx23", 100, -1, 1 );
  dx23Histo->setTitle( "x3-x2;x_{3}-x_{2} [mm];hit pairs" );

  dy23Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dy23", 100, -1, 1 );
  dy23Histo->setTitle( "y3-y2;y_{3}-y_{2} [mm];hit pairs" );

  du23Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/du23", 100, -1, 1 );
  du23Histo->setTitle( "x3-x2, |dy| < 1;x_{3}-x_{2} [mm];hit pairs" );

  dx34Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dx34", 100, -1, 1 );
  dx34Histo->setTitle( "x4-x3;x_{4}-x_{3} [mm];hit pairs" );

  dy34Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dy34", 100, -1, 1 );
  dy34Histo->setTitle( "y4-y3;y_{4}-y_{3} [mm];hit pairs" );

  du34Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/du34", 100, -1, 1 );
  du34Histo->setTitle( "x4-x3, |dy| < 1;x_{4}-x_{3} [mm];hit pairs" );

  dx45Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dx45", 100, -1, 1 );
  dx45Histo->setTitle( "x5-x4;x_{5}-x_{4} [mm];hit pairs" );

  dy45Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dy45", 100, -1, 1 );
  dy45Histo->setTitle( "y5-y4;y_{5}-y_{4} [mm];hit pairs" );

  du45Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/du45", 100, -1, 1 );
  du45Histo->setTitle( "x5-x4, |dy| < 1;x_{5}-x_{4} [mm];hit pairs" );

  // triplets:
  AIDAProcessor::tree(this)->mkdir("Upstream");

  dzcvsxy = AIDAProcessor::histogramFactory(this)->
    createProfile2D( "Upstream/dzcvsxy", 120, -12, 12, 60, -6, 6, -999, 999 );
  dzcvsxy->setTitle( "DUT plane;telescope track x_{DUT} [mm];telescope track y_{DUT} [mm];z_{DUT} [mm]" );

  z3vsxy = AIDAProcessor::histogramFactory(this)->
    createProfile2D( "Upstream/z3vsxy", 120, -12, 12, 60, -6, 6, -999, 999 );
  z3vsxy->setTitle( "DUT plane;telescope track x_{DUT} [mm];telescope track y_{DUT} [mm];z_{DUT} [mm]" );

  tridxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridx", 100, -100, 100 );
  tridxHisto->setTitle( "triplet dx;x_{1}-x_{m} [#mum];telescope triplets" );

  tridyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridy", 100, -100, 100 );
  tridyHisto->setTitle( "triplet dy;y_{1}-y_{m} [#mum];telescope triplets" );

  tridx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridx1", 100, -1, 1 );
  tridx1Histo->setTitle( "triplet dx;x_{1}-x_{t} [mm];telescope triplets" );

  tridy1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridy1", 100, -1, 1 );
  tridy1Histo->setTitle( "triplet dy;y_{1}-y_{t} [mm];telescope triplets" );

  tridxvsx = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Upstream/tridxvsx", 100, -10, 10, -100, 100 );
  tridxvsx->setTitle( "triplet x resid vs x;x [mm];triplet <#Deltax> [#mum]" );

  tridxvsy = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Upstream/tridxvsy", 50, -5, 5, -100, 100 );
  tridxvsy->setTitle( "triplet x resid vs y;y [mm];triplet <#Deltax> [#mum]" );

  tridxvstx = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Upstream/tridxvstx", 80, -2, 2, -100, 100 );
  tridxvstx->setTitle( "triplet x resid vs tx;t_{x} [mrad];triplet <#Deltax> [#mum]" );

  tridxvsty = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Upstream/tridxvsty", 80, -2, 2, -100, 100 );
  tridxvsty->setTitle( "triplet x resid vs ty;t_{y} [mrad];triplet <#Deltax> [#mum]" );

  tridyvsx = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Upstream/tridyvsx", 100, -10, 10, -100, 100 );
  tridyvsx->setTitle( "triplet y resid vs x;x [mm];triplet <#Deltay> [#mum]" );

  tridyvsy = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Upstream/tridyvsy", 50, -5, 5, -100, 100 );
  tridyvsy->setTitle( "triplet y resid vs y;y [mm];triplet <#Deltay> [#mum]" );

  tridyvstx = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Upstream/tridyvstx", 80, -2, 2, -100, 100 );
  tridyvstx->setTitle( "triplet y resid vs tx;t_{x} [mrad];triplet <#Deltay> [#mum]" );

  tridyvsty = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Upstream/tridyvsty", 80, -2, 2, -100, 100 );
  tridyvsty->setTitle( "triplet y resid vs ty;t_{y} [mrad];triplet <#Deltay> [#mum]" );

  tridx3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridx3", 100, -1, 1 );
  tridx3Histo->setTitle( "triplet dx;x_{3}-x_{t} [mm];telescope triplets" );

  tridy3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridy3", 100, -1, 1 );
  tridy3Histo->setTitle( "triplet dy;y_{3}-y_{t} [mm];telescope triplets" );

  tridx3bHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridx3b", 100, -250, 250 );
  tridx3bHisto->setTitle( "triplet dx;x_{3}-x_{t} [um];telescope triplets" );

  tridy3bHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridy3b", 100, -250, 250 );
  tridy3bHisto->setTitle( "triplet dy;y_{3}-y_{t} [um];telescope triplets" );

  tridx4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridx4", 100, -2, 2 );
  tridx4Histo->setTitle( "triplet dx;x_{4}-x_{t} [mm];telescope triplets" );

  tridy4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridy4", 100, -2, 2 );
  tridy4Histo->setTitle( "triplet dy;y_{4}-y_{t} [mm];telescope triplets" );

  tridx4bHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridx4b", 100, -400, 400 );
  tridx4bHisto->setTitle( "triplet dx;x_{4}-x_{t} [um];telescope triplets" );

  tridy4bHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridy4b", 100, -400, 400 );
  tridy4bHisto->setTitle( "triplet dy;y_{4}-y_{t} [um];telescope triplets" );

  tridx5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridx5", 100, -3, 3 );
  tridx5Histo->setTitle( "triplet dx;x_{5}-x_{t} [mm];telescope triplets" );

  tridy5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridy5", 100, -3, 3 );
  tridy5Histo->setTitle( "triplet dy;y_{5}-y_{t} [mm];telescope triplets" );

  tridx5bHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridx5b", 100, -1000, 1000 );
  tridx5bHisto->setTitle( "triplet dx;x_{5}-x_{t} [um];telescope triplets" );

  tridy5bHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridy5b", 100, -1000, 1000 );
  tridy5bHisto->setTitle( "triplet dy;y_{5}-y_{t} [um];telescope triplets" );

  trixHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/trix", 240, -12, 12 );
  trixHisto->setTitle( "triplet x1;x1_{out} [mm];telescope triplets" );

  triyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/triy", 120, -6, 6 );
  triyHisto->setTitle( "triplet y1;y1_{up} [mm];telescope triplets" );

  trixyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "Upstream/trixy", 240, -12, 12, 120, -6, 6 );
  trixyHisto->setTitle( "triplet y1 vs x1;x1_{out} [mm];y1_{up} [mm];telescope triplets" );

  tritxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tritx", 100, -10, 10 );
  tritxHisto->setTitle( "triplet slope x;#theta_{x} [mrad];telescope triplets" );

  trityHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/trity", 100, -10, 10 );
  trityHisto->setTitle( "triplet slope y;#theta_{y} [mrad];telescope triplets" );

  trixdutHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/trixdut", 240, -12, 12 );
  trixdutHisto->setTitle( "triplet at DUT;triplet x_{out} at DUT [mm];telescope triplets" );

  triydutHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/triydut", 120, -6, 6 );
  triydutHisto->setTitle( "triplet at DUT;triplet y_{up} at DUT [mm];telescope triplets" );

  trixydutHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "Upstream/trixydut", 240, -12, 12, 120, -6, 6 );
  trixydutHisto->setTitle( "triplet at DUT;triplet x_{out} at DUT [mm];triplet y_{up} at DUT [mm];telescope triplets" );

  triddaMindutHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "triddaMindut", 120, 0, 10 );
  triddaMindutHisto->setTitle( "minimal triplet distance at DUT;triplet distance at DUT [mm];telescope triplets" );


  ntriHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "ntri", 31, -0.5, 30.5 );
  ntriHisto->setTitle( "telescope triplets;0-1-2 triplets;events" );

  // driplets:
  AIDAProcessor::tree(this)->mkdir("Downstream");

  dridxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Downstream/dridx", 100, -100, 100 );
  dridxHisto->setTitle( "driplet dx;x_{4}-x_{m} [#mum];driplets" );

  dridyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Downstream/dridy", 100, -100, 100 );
  dridyHisto->setTitle( "driplet dy;y_{4}-y_{m} [#mum];driplets" );

  dridxvsx = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Downstream/dridxvsx", 100, -10, 10, -100, 100 );
  dridxvsx->setTitle( "driplet x resid vs x;x [mm];driplet <#Deltax> [#mum]" );

  dridxvsy = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Downstream/dridxvsy", 50, -5, 5, -100, 100 );
  dridxvsy->setTitle( "driplet x resid vs y;y [mm];driplet <#Deltax> [#mum]" );

  dridxvstx = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Downstream/dridxvstx", 80, -2, 2, -100, 100 );
  dridxvstx->setTitle( "driplet x resid vs tx;t_{x} [mrad];driplet <#Deltax> [#mum]" );

  dridxvsty = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Downstream/dridxvsty", 80, -2, 2, -100, 100 );
  dridxvsty->setTitle( "driplet x resid vs ty;t_{y} [mrad];driplet <#Deltax> [#mum]" );

  dridyvsx = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Downstream/dridyvsx", 100, -10, 10, -100, 100 );
  dridyvsx->setTitle( "driplet y resid vs x;x [mm];driplet <#Deltay> [#mum]" );

  dridyvsy = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Downstream/dridyvsy", 50, -5, 5, -100, 100 );
  dridyvsy->setTitle( "driplet y resid vs y;y [mm];driplet <#Deltay> [#mum]" );

  dridyvstx = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Downstream/dridyvstx", 80, -2, 2, -100, 100 );
  dridyvstx->setTitle( "driplet y resid vs tx;t_{x} [mrad];driplet <#Deltay> [#mum]" );

  dridyvsty = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Downstream/dridyvsty", 80, -2, 2, -100, 100 );
  dridyvsty->setTitle( "driplet y resid vs ty;t_{y} [mrad];driplet <#Deltay> [#mum]" );

  drixHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Downstream/drix", 240, -12, 12 );
  drixHisto->setTitle( "driplet x4;x4_{out} [mm];driplets" );

  driyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Downstream/driy", 120, -6, 6 );
  driyHisto->setTitle( "driplet y4;y4_{up} [mm];driplets" );

  drixyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "Downstream/drixy", 240, -12, 12, 120, -6, 6 );
  drixyHisto->setTitle( "driplet y4 vs x4;x4_{out} [mm];y4_{up} [mm];driplets" );

  dritxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Downstream/dritx", 100, -10, 10 );
  dritxHisto->setTitle( "driplet slope x;#theta_{x} [mrad];driplets" );

  drityHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Downstream/drity", 100, -10, 10 );
  drityHisto->setTitle( "driplet slope y;#theta_{y} [mrad];driplets" );



  bacsxaHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "bacsxa", 100, -10, 10 );
  bacsxaHisto->setTitle( "DUT + driplet x;DUT cluster + driplet #Sigmax [mm];DUT clusters" );

  bacdyaHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "bacdya", 100, -10, 10 );
  bacdyaHisto->setTitle( "DUT - driplet y;DUT cluster - driplet #Deltay [mm];DUT clusters" );

  bacsxcHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "bacsxc", 200, -500, 500 );
  bacsxcHisto->setTitle( "DUT + driplet x;DUT cluster + driplet #Sigmax [mm];DUT clusters" );

  bacdycHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "bacdyc", 200, -500, 500 );
  bacdycHisto->setTitle( "DUT - driplet y;DUT cluster - driplet #Deltay [mm];DUT clusters" );

  bacsxcqHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "bacsxcq", 200, -500, 500 );
  bacsxcqHisto->setTitle( "DUT + driplet x;DUT cluster + driplet #Sigmax [mm];DUT clusters" );

  bacdycqHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "bacdycq", 100, -500, 500 );
  bacdycqHisto->setTitle( "DUT - driplet y;DUT cluster - driplet #Deltay [mm];DUT clusters" );


  ndriHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Downstream/ndri", 31, -0.5, 30.5 );
  ndriHisto->setTitle( "telescope driplets;3-4-5 driplets;events" );

  //driplets-triplets
  // Tracks
  AIDAProcessor::tree(this)->mkdir("Tracks");

  nsixHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/nsix", 21, -0.5, 20.5 );
  nsixHisto->setTitle( "telescope six-plane-tracks;six-plane-tracks;events" );

  sixkxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixkx", 100, -10, 10 );
  sixkxHisto->setTitle( "kink x;kink x [mrad];triplet-driplet pairs" );

  sixkyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixky", 100, -10, 10 );
  sixkyHisto->setTitle( "kink y;kink y [mrad];triplet-driplet pairs" );

  sixdxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixdx", 100, -1, 1 );
  sixdxHisto->setTitle( "six match x;match x [mm];triplet-driplet pairs" );

  sixdyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixdy", 100, -1, 1 );
  sixdyHisto->setTitle( "six match y;match y [mm];triplet-driplet pairs" );

  sixdxcHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixdxc", 100, -250, 250 );
  sixdxcHisto->setTitle( "six match x;track #Deltax[#mum];triplet-driplet pairs" );

  sixdycHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixdyc", 100, -250, 250 );
  sixdycHisto->setTitle( "six match y;track #Deltay[#mum];triplet-driplet pairs" );

  sixkxcHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixkxc", 100, -10, 10 );
  sixkxcHisto->setTitle( "kink x, x-y matched;kink x [mrad];triplet-driplet pairs" );

  sixkycHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixkyc", 100, -10, 10 );
  sixkycHisto->setTitle( "kink y, x-y matched;kink y [mrad];triplet-driplet pairs" );

  sixxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixx", 240, -12, 12 );
  sixxHisto->setTitle( "six x at DUT;six x_{out} at DUT [mm];six-plane tracks" );

  sixyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixy", 120, -6, 6 );
  sixyHisto->setTitle( "six y at DUT;six y_{up} at DUT [mm];six-plane tracks" );

  sixxyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "Tracks/sixxy", 240, -12, 12, 120, -6, 6 );
  sixxyHisto->setTitle( "six at DUT;six x_{out} at DUT [mm];six y_{up} at DUT [mm];six-plane tracks" );

  sixxycHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "Tracks/sixxyc", 240, -12, 12, 120, -6, 6 );
  sixxycHisto->setTitle( "six large kink;six x_{out} at DUT [mm];six y_{up} at DUT [mm];large kink tracks" );

  sixxylkHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "Tracks/sixxylk", 240, -12, 12, 120, -6, 6 );
  sixxylkHisto->setTitle( "six with REF link at DUT;six x_{out} at DUT [mm];six y_{up} at DUT [mm];six-plane tracks with REF link" );

  kinkvsxy = AIDAProcessor::histogramFactory(this)->
    createProfile2D( "Tracks/kinkvsxy", 120, -12, 12, 60, -6, 6, 0, 100 );
  kinkvsxy->setTitle( "kink;six x_{out} at DUT [mm];six y_{up} at DUT [mm];<kink^{2}> [mrad^{2}]" );

  sixx0Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixx0", 240, -12, 12 );
  sixx0Histo->setTitle( "six x at 0;six x_{out} at 0 [mm];six-plane tracks" );

  sixy0Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixy0", 120, -6, 6 );
  sixy0Histo->setTitle( "six y at 0;six y_{up} at 0 [mm];six-plane tracks" );

  sixx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixx1", 240, -12, 12 );
  sixx1Histo->setTitle( "six x at 1;six x_{out} at 1 [mm];six-plane tracks" );

  sixy1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixy1", 120, -6, 6 );
  sixy1Histo->setTitle( "six y at 1;six y_{up} at 1 [mm];six-plane tracks" );

  sixx2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixx2", 240, -12, 12 );
  sixx2Histo->setTitle( "six x at 2;six x_{out} at 2 [mm];six-plane tracks" );

  sixy2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixy2", 120, -6, 6 );
  sixy2Histo->setTitle( "six y at 2;six y_{up} at 2 [mm];six-plane tracks" );

  sixx3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixx3", 240, -12, 12 );
  sixx3Histo->setTitle( "six x at 3;six x_{out} at 3 [mm];six-plane tracks" );

  sixy3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixy3", 120, -6, 6 );
  sixy3Histo->setTitle( "six y at 3;six y_{up} at 3 [mm];six-plane tracks" );

  sixx4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixx4", 240, -12, 12 );
  sixx4Histo->setTitle( "six x at 4;six x_{out} at 4 [mm];six-plane tracks" );

  sixy4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixy4", 120, -6, 6 );
  sixy4Histo->setTitle( "six y at 4;six y_{up} at 4 [mm];six-plane tracks" );

  sixx5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixx5", 240, -12, 12 );
  sixx5Histo->setTitle( "six x at 5;six x_{out} at 5 [mm];six-plane tracks" );

  sixy5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixy5", 120, -6, 6 );
  sixy5Histo->setTitle( "six y at 5;six y_{up} at 5 [mm];six-plane tracks" );

  // GBL:
  AIDAProcessor::tree(this)->mkdir("GBL");

  derxtiltHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/derxtilt", 100, -0.1, 0.1 );
  derxtiltHisto->setTitle( "ddx/dtilt;ddx/dtilt [mm/rad];align hits" );

  derytiltHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/derytilt", 100, -10, 10 );
  derytiltHisto->setTitle( "ddy/dtilt;ddy/dtilt [mm/rad];align hits" );

  derxturnHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/derxturn", 100, -10, 10 );
  derxturnHisto->setTitle( "ddx/dturn;ddx/dturn [mm/rad];align hits" );

  deryturnHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/deryturn", 100, -1, 1 );
  deryturnHisto->setTitle( "ddy/dturn;ddy/dturn [mm/rad];align hits" );

  selxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/selx", 240, -12, 12 );
  selxHisto->setTitle( "x at DUT, sel GBL;six x_{out} at DUT [mm];selected tracks" );

  selyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/sely", 120, -6, 6 );
  selyHisto->setTitle( "y at DUT, sel GBL;six y_{up} at DUT [mm];selected tracks" );

  selaxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/selax", 100, -5, 5 );
  selaxHisto->setTitle( "track angle x, sel GBL;x angle [mrad];tracks" );

  selayHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/selay", 100, -5, 5 );
  selayHisto->setTitle( "track angle y, sel GBL;y angle [mrad];tracks" );

  seldxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/seldx", 100, -150, 150 );
  seldxHisto->setTitle( "track match x, sel GBL;#Deltax [#mum];tracks" );

  seldyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/seldy", 100, -150, 150 );
  seldyHisto->setTitle( "track match y, sel GBL;#Deltay [#mum];tracks" );

  selkxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/selkx", 100, -10, 10 );
  selkxHisto->setTitle( "kink x, sel GBL;kink x [mrad];tracks" );

  selkyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/selky", 100, -10, 10 );
  selkyHisto->setTitle( "kink y, sel GBL;kink y [mrad];tracks" );

  seldx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/seldx1", 100, -100, 100 );
  seldx1Histo->setTitle( "triplet resid x at 1, sel GBL;#Deltax [#mum];tracks" );

  seldy1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/seldy1", 100, -100, 100 );
  seldy1Histo->setTitle( "triplet resid y at 1, sel GBL;#Deltay [#mum];tracks" );

  seldx3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/seldx3", 100, -1000, 1000 );
  seldx3Histo->setTitle( "triplet resid x at 3, sel GBL;#Deltax [#mum];tracks" );

  seldy3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/seldy3", 100, -1000, 1000 );
  seldy3Histo->setTitle( "triplet resid y at 3, sel GBL;#Deltay [#mum];tracks" );

  seldx4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/seldx4", 100, -1500, 1500 );
  seldx4Histo->setTitle( "triplet resid x at 4, sel GBL;#Deltax [#mum];tracks" );

  seldy4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/seldy4", 100, -1500, 1500 );
  seldy4Histo->setTitle( "triplet resid y at 4, sel GBL;#Deltay [#mum];tracks" );

  seldx5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/seldx5", 100, -3000, 3000 );
  seldx5Histo->setTitle( "triplet resid x at 5, sel GBL;#Deltax [#mum];tracks" );

  seldy5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/seldy5", 100, -3000, 3000 );
  seldy5Histo->setTitle( "triplet resid y at 5, sel GBL;#Deltay [#mum];tracks" );

  seldx6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/seldx6", 100, -500, 500 );
  seldx6Histo->setTitle( "triplet resid x at DUT, sel GBL;#Deltax [#mum];tracks" );

  seldy6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/seldy6", 100, -500, 500 );
  seldy6Histo->setTitle( "triplet resid y at DUT, sel GBL;#Deltay [#mum];tracks" );

  gblndfHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblndf", 16, -0.5, 15.5 );
  gblndfHisto->setTitle( "GBL fit NDF;GBL NDF;tracks" );

  gblchi2aHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblchi2a", 100, 0, 100 );
  gblchi2aHisto->setTitle( "GBL fit chi2, DoF 8;GBL chi2;tracks" );

  gblchi2bHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblchi2b", 100, 0, 100 );
  gblchi2bHisto->setTitle( "GBL fit chi2;GBL chi2;tracks" );

  gblprbHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblprb", 100, 0, 1 );
  gblprbHisto->setTitle( "GBL fit probability;GBL fit probability;tracks" );

  // bad fits:

  badxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/badx", 240, -12, 12 );
  badxHisto->setTitle( "x at DUT, bad GBL;six x_{out} at DUT [mm];bad tracks" );

  badyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/bady", 120, -6, 6 );
  badyHisto->setTitle( "y at DUT, bad GBL;six y_{up} at DUT [mm];bad tracks" );

  badaxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/badax", 100, -5, 5 );
  badaxHisto->setTitle( "track angle x, bad GBL;x angle [mrad];tracks" );

  badayHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baday", 100, -5, 5 );
  badayHisto->setTitle( "track angle y, bad GBL;y angle [mrad];tracks" );

  baddxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baddx", 100, -150, 150 );
  baddxHisto->setTitle( "track match x, bad GBL;#Deltax [#mum];tracks" );

  baddyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baddy", 100, -150, 150 );
  baddyHisto->setTitle( "track match y, bad GBL;#Deltay [#mum];tracks" );

  badkxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/badkx", 100, -10, 10 );
  badkxHisto->setTitle( "kink x, bad GBL;kink x [mrad];tracks" );

  badkyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/badky", 100, -10, 10 );
  badkyHisto->setTitle( "kink y, bad GBL;kink y [mrad];tracks" );

  baddx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baddx1", 100, -100, 100 );
  baddx1Histo->setTitle( "triplet resid x at 1, bad GBL;#Deltax [#mum];tracks" );

  baddy1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baddy1", 100, -100, 100 );
  baddy1Histo->setTitle( "triplet resid y at 1, bad GBL;#Deltay [#mum];tracks" );

  baddx3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baddx3", 100, -1000, 1000 );
  baddx3Histo->setTitle( "triplet resid x at 3, bad GBL;#Deltax [#mum];tracks" );

  baddy3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baddy3", 100, -1000, 1000 );
  baddy3Histo->setTitle( "triplet resid y at 3, bad GBL;#Deltay [#mum];tracks" );

  baddx4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baddx4", 100, -1500, 1500 );
  baddx4Histo->setTitle( "triplet resid x at 4, bad GBL;#Deltax [#mum];tracks" );

  baddy4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baddy4", 100, -1500, 1500 );
  baddy4Histo->setTitle( "triplet resid y at 4, bad GBL;#Deltay [#mum];tracks" );

  baddx5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baddx5", 100, -3000, 3000 );
  baddx5Histo->setTitle( "triplet resid x at 5, bad GBL;#Deltax [#mum];tracks" );

  baddy5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baddy5", 100, -3000, 3000 );
  baddy5Histo->setTitle( "triplet resid y at 5, bad GBL;#Deltay [#mum];tracks" );

  baddx6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baddx6", 100, -250, 250 );
  baddx6Histo->setTitle( "triplet resid x at DUT, bad GBL;#Deltax [#mum];tracks" );

  baddy6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baddy6", 100, -250, 250 );
  baddy6Histo->setTitle( "triplet resid y at DUT, bad GBL;#Deltay [#mum];tracks" );

  // good fits:

  goodx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/goodx1", 100, -100, 100 );
  goodx1Histo->setTitle( "triplet resid x at 1, good GBL;#Deltax [#mum];tracks" );

  goody1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/goody1", 100, -100, 100 );
  goody1Histo->setTitle( "triplet resid y at 1, good GBL;#Deltay [#mum];tracks" );

  goodx6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/goodx6", 100, -250, 250 );
  goodx6Histo->setTitle( "triplet resid x at 6, good GBL;#Deltax [#mum];tracks" );

  goody6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/goody6", 100, -250, 250 );
  goody6Histo->setTitle( "triplet resid y at 6, good GBL;#Deltay [#mum];tracks" );

  // look at fit:

  gblax0Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblax0", 100, -1, 1 );
  gblax0Histo->setTitle( "GBL x angle at plane 0;x angle at plane 0 [mrad];tracks" );

  gbldx0Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gbldx0", 100, -10, 10 );
  gbldx0Histo->setTitle( "GBL x shift at plane 0;x shift at plane 0 [#mum];tracks" );

  gblrx0Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblrx0", 100, -25, 25 );
  gblrx0Histo->setTitle( "GBL x resid at plane 0;x resid at plane 0 [#mum];tracks" );

  gblpx0Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpx0", 100, -10, 10 );
  gblpx0Histo->setTitle( "GBL x pull at plane 0;x pull at plane 0 [#sigma];tracks" );

  gblqx0Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblqx0", 100, -1, 1 );
  gblqx0Histo->setTitle( "GBL x kink at plane 0;x kink at plane 0 [mrad];tracks" );


  gblax1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblax1", 100, -1, 1 );
  gblax1Histo->setTitle( "GBL x angle at plane 1;x angle at plane 1 [mrad];tracks" );

  gbldx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gbldx1", 100, -100, 100 );
  gbldx1Histo->setTitle( "GBL x shift at plane 1;x shift at plane 1 [#mum];tracks" );

  gblrx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblrx1", 100, -25, 25 );
  gblrx1Histo->setTitle( "GBL x resid at plane 1;x resid at plane 1 [#mum];tracks" );

  gblpx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpx1", 100, -10, 10 );
  gblpx1Histo->setTitle( "GBL x pull at plane 1;x pull at plane 1 [#sigma];tracks" );

  gblqx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblqx1", 100, -1, 1 );
  gblqx1Histo->setTitle( "GBL x kink at plane 1;x kink at plane 1 [mrad];tracks" );

  gblsx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblsx1", 100, -10, 10 );
  gblsx1Histo->setTitle( "GBL x kink at plane 1/error;x kink at plane 1/error;tracks" );

  gbltx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gbltx1", 100, -10, 10 );
  gbltx1Histo->setTitle( "GBL x kink pull at plane 1;x kink pull at plane 1;tracks" );


  gblax2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblax2", 100, -1, 1 );
  gblax2Histo->setTitle( "GBL x angle at plane 2;x angle at plane 2 [mrad];tracks" );

  gbldx2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gbldx2", 100, -20, 20 );
  gbldx2Histo->setTitle( "GBL x shift at plane 2;x shift at plane 2 [#mum];tracks" );

  gblrx2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblrx2", 100, -25, 25 );
  gblrx2Histo->setTitle( "GBL x resid at plane 2;x resid at plane 2 [#mum];tracks" );

  gblpx2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpx2", 100, -10, 10 );
  gblpx2Histo->setTitle( "GBL x pull at plane 2;x pull at plane 2 [#sigma];tracks" );

  gblqx2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblqx2", 100, -1, 1 );
  gblqx2Histo->setTitle( "GBL x kink at plane 2;x kink at plane 2 [mrad];tracks" );


  gblax3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblax3", 100, -1, 1 );
  gblax3Histo->setTitle( "GBL x angle at plane 3;x angle at plane 3 [mrad];tracks" );

  gbldx3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gbldx3", 100, -250, 250 );
  gbldx3Histo->setTitle( "GBL x shift at plane 3;x shift at plane 3 [#mum];tracks" );

  gblrx3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblrx3", 100, -25, 25 );
  gblrx3Histo->setTitle( "GBL x resid at plane 3;x resid at plane 3 [#mum];tracks" );

  gblpx3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpx3", 100, -10, 10 );
  gblpx3Histo->setTitle( "GBL x pull at plane 3;x pull at plane 3 [#sigma];tracks" );

  gblqx3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblqx3", 100, -1, 1 );
  gblqx3Histo->setTitle( "GBL x kink at plane 3;x kink at plane 3 [mrad];tracks" );


  gblax4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblax4", 100, -1, 1 );
  gblax4Histo->setTitle( "GBL x angle at plane 4;x angle at plane 4 [mrad];tracks" );

  gbldx4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gbldx4", 100, -500, 500 );
  gbldx4Histo->setTitle( "GBL x shift at plane 4;x shift at plane 4 [#mum];tracks" );

  gblrx4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblrx4", 100, -25, 25 );
  gblrx4Histo->setTitle( "GBL x resid at plane 4;x resid at plane 4 [#mum];tracks" );

  gblpx4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpx4", 100, -10, 10 );
  gblpx4Histo->setTitle( "GBL x pull at plane 4;x pull at plane 4 [#sigma];tracks" );

  gblqx4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblqx4", 100, -1, 1 );
  gblqx4Histo->setTitle( "GBL x kink at plane 4;x kink at plane 4 [mrad];tracks" );


  gblax5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblax5", 100, -1, 1 );
  gblax5Histo->setTitle( "GBL x angle at plane 5;x angle at plane 5 [mrad];tracks" );

  gbldx5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gbldx5", 100, -1000, 1000 );
  gbldx5Histo->setTitle( "GBL x shift at plane 5;x shift at plane 5 [#mum];tracks" );

  gblrx5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblrx5", 100, -25, 25 );
  gblrx5Histo->setTitle( "GBL x resid at plane 5;x resid at plane 5 [#mum];tracks" );

  gblpx5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpx5", 100, -10, 10 );
  gblpx5Histo->setTitle( "GBL x pull at plane 5;x pull at plane 5 [#sigma];tracks" );

  gblqx5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblqx5", 100, -1, 1 );
  gblqx5Histo->setTitle( "GBL x kink at plane 5;x kink at plane 5 [mrad];tracks" );


  gblax6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblax6", 100, -1, 1 );
  gblax6Histo->setTitle( "GBL x angle at DUT;x angle at DUT [mrad];tracks" );

  gbldx6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gbldx6", 100, -1000, 1000 );
  gbldx6Histo->setTitle( "GBL x shift at DUT;x shift at DUT [#mum];tracks" );

  gbldy6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gbldy6", 100, -1000, 1000 );
  gbldy6Histo->setTitle( "GBL y shift at DUT;y shift at DUT [#mum];tracks" );

  gblrx6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblrx6", 100, -250, 250 );
  gblrx6Histo->setTitle( "GBL x resid at DUT;x resid at DUT [#mum];tracks" );

  gblry6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblry6", 100, -100, 100 );
  gblry6Histo->setTitle( "GBL y resid at DUT;y resid at DUT [#mum];tracks" );

  gblpx6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpx6", 100, -10, 10 );
  gblpx6Histo->setTitle( "GBL x pull at DUT;x pull at DUT [#sigma];tracks" );

  gblpy6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpy6", 100, -10, 10 );
  gblpy6Histo->setTitle( "GBL y pull at DUT;y pull at DUT [#sigma];tracks" );

  gblqx6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblqx6", 100, -10, 10 );
  gblqx6Histo->setTitle( "GBL x kink at DUT;x kink at DUT [mrad];tracks" );

  gblsx6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblsx6", 100, -10, 10 );
  gblsx6Histo->setTitle( "GBL x kink at DUT/error;x kink at DUT/error;tracks" );

  gbltx6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gbltx6", 100, -10, 10 );
  gbltx6Histo->setTitle( "GBL x kink pull at DUT;x kink pull at DUT;tracks" );


  gblkx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblkx1", 100, -1, 1 );
  gblkx1Histo->setTitle( "GBL kink angle at plane 1;plane 1 kink [mrad];tracks" );

  gblkx2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblkx2", 100, -1, 1 );
  gblkx2Histo->setTitle( "GBL kink angle at plane 2;plane 2 kink [mrad];tracks" );

  gblkx3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblkx3", 100, -10, 10 );
  gblkx3Histo->setTitle( "GBL kink angle at plane 3;plane 3 kink [mrad];tracks" );

  gblkx4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblkx4", 100, -1, 1 );
  gblkx4Histo->setTitle( "GBL kink angle at plane 4;plane 4 kink [mrad];tracks" );

  gblkx5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblkx5", 100, -1, 1 );
  gblkx5Histo->setTitle( "GBL kink angle at plane 5;plane 5 kink [mrad];tracks" );

  gblkx6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblkx6", 100, -1, 1 );
  gblkx6Histo->setTitle( "GBL kink angle at plane 6;plane 6 kink [mrad];tracks" );

  // z intersect:

  sixzx3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixzx3", 100, -50, 250 );
  sixzx3Histo->setTitle( "intersect z-x, kink > 3 mrad;intersect z(x) - z_{2} [mm];tracks" );

  sixzy3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixzy3", 100, -50, 250 );
  sixzy3Histo->setTitle( "intersect z-y, kink > 3 mrad;intersect z(y) - z_{2} [mm];tracks" );

  sixzx2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixzx2", 100, -50, 250 );
  sixzx2Histo->setTitle( "intersect z-x, kink > 2 mrad;intersect z(x) - z_{2} [mm];tracks" );

  sixzy2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixzy2", 100, -50, 250 );
  sixzy2Histo->setTitle( "intersect z-y, kink > 2 mrad;intersect z(y) - z_{2} [mm];tracks" );

  sixzx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixzx1", 100, -50, 250 );
  sixzx1Histo->setTitle( "intersect z-x, kink > 1 mrad;intersect z(x) - z_{2} [mm];tracks" );

  sixzy1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixzy1", 100, -50, 250 );
  sixzy1Histo->setTitle( "intersect z-y, kink > 1 mrad;intersect z(y) - z_{2} [mm];tracks" );

  sixkxzyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixkxzy", 200, -20, 20 );
  sixkxzyHisto->setTitle( "kink x at DUT;kink x [mrad];tracks" );

  sixkyzxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixkyzx", 200, -20, 20 );
  sixkyzxHisto->setTitle( "kink y at DUT;kink y [mrad];tracks" );

  sixkxzxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixkxzx", 200, -20, 20 );
  sixkxzxHisto->setTitle( "kink x at DUT;kink x [mrad];tracks" );

  sixkyzyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixkyzy", 200, -20, 20 );
  sixkyzyHisto->setTitle( "kink y at DUT;kink y [mrad];tracks" );

  streamlog_out( DEBUG2 ) << "Histogram booking completed \n\n" << std::endl;

#else

  streamlog_out( MESSAGE4 ) << "No histogram produced because Marlin doesn't use AIDA" << std::endl;

#endif

  return;
}

TMatrixD EUTelTripletGBL::JacobianPointToPoint( double ds ) {
  /* for GBL:
     Jacobian for straight line track
     track = q/p, x', y', x, y
     0,   1,  2,  3, 4
  */
  TMatrixD jac( 5, 5 );
  jac.UnitMatrix();
  jac[3][1] = ds; // x = x0 + xp * ds
  jac[4][2] = ds; // y = y0 + yp * ds
  return jac;
}


void EUTelTripletGBL::TelescopeCorrelationPlots(std::vector<hit> * telescopehits) {
  for( std::vector<hit>::iterator ihit = telescopehits->begin(); ihit != telescopehits->end(); ihit++ ){

    int ipl = (*ihit).plane;

    for( std::vector<hit>::iterator jhit = telescopehits->begin(); jhit != telescopehits->end(); jhit++ ){

      int jpl = (*jhit).plane;
      double dx = (*jhit).x - (*ihit).x;
      double dy = (*jhit).y - (*ihit).y;

      if( ipl == 0 ) {
	if( jpl == 1 ) {
	  dx01Histo->fill( dx );
	  dy01Histo->fill( dy );
	  if( abs(dy) < 1 ) du01Histo->fill( dx );
	}
	if( jpl == 2 ) {
	  dx02Histo->fill( dx );
	}
	if( jpl == 3 ) {
	  dx03Histo->fill( dx );
	}
	if( jpl == 4 ) {
	  dx04Histo->fill( dx );
	}
	if( jpl == 5 ) {
	  dx05Histo->fill( dx );
	}

      }//ipl==0

      if( ipl == 1 ) {
	if( jpl == 2 ) {
	  dx12Histo->fill( dx );
	  dy12Histo->fill( dy );
	  if( abs(dy) < 1 ) du12Histo->fill( dx );
	}
      }

      if( ipl == 2 ) {
	if( jpl == 3 ) {
	  dx23Histo->fill( dx );
	  dy23Histo->fill( dy );
	  if( abs(dy) < 1 ) du23Histo->fill( dx );
	}
      }

      if( ipl == 3 ) {
	if( jpl == 4 ) {
	  dx34Histo->fill( dx );
	  dy34Histo->fill( dy );
	  if( abs(dy) < 1 ) du34Histo->fill( dx );
	}
      }

      if( ipl == 4 ) {
	if( jpl == 5 ) {
	  dx45Histo->fill( dx );
	  dy45Histo->fill( dy );
	  if( abs(dy) < 1 ) du45Histo->fill( dx );
	}
      }
    }
  }
}


std::vector<EUTelTripletGBL::triplet> * EUTelTripletGBL::FindTriplets(std::vector<hit> * hits, unsigned int plane0, unsigned int plane1, unsigned int plane2) {

  std::vector<triplet> * triplets = new std::vector<triplet>;

  for( std::vector<hit>::iterator ihit = hits->begin(); ihit != hits->end(); ihit++ ){
    if( (*ihit).plane != plane0 ) continue; // First plane

    for( std::vector<hit>::iterator jhit = hits->begin(); jhit != hits->end(); jhit++ ){
      if( (*jhit).plane != plane2 ) continue; // Last plane

      for( std::vector<hit>::iterator khit = hits->begin(); khit != hits->end(); khit++ ){
	if( (*khit).plane != plane1 ) continue; // Middle plane

	// Create new preliminary triplet from the three hits:
	triplet new_triplet((*ihit),(*khit),(*jhit));

	// Setting cuts on the triplet track angle:
	if( abs(new_triplet.getdy()) > triplet_angle_cut * new_triplet.getdz()) continue;
	if( abs(new_triplet.getdx()) > triplet_angle_cut * new_triplet.getdz()) continue;

	// Setting cuts on the triplet residual on the middle plane
	if( abs(new_triplet.getdx(plane1)) > triplet_residual_cut) continue;
	if( abs(new_triplet.getdy(plane1)) > triplet_residual_cut) continue;

	// The triplet is accepted, push it back:
	triplets->push_back(new_triplet);
	streamlog_out(DEBUG2) << new_triplet;

      }//loop over hits
    }//loop over hits
  }// loop over hits

  return triplets;
}


std::vector<EUTelTripletGBL::track> * EUTelTripletGBL::MatchTriplets(std::vector<triplet> * up, std::vector<triplet> * down, double z_match) {
  
  std::vector<track> * tracks = new std::vector<track>;

  for( std::vector<triplet>::iterator trip = up->begin(); trip != up->end(); trip++ ){

    for( std::vector<triplet>::iterator drip = down->begin(); drip != down->end(); drip++ ){

      // Build a track candidate from one upstream and one downstream triplet:
      track newtrack((*trip),(*drip));

      // Track kinks as difference in triplet slopes:
      double kx = newtrack.kink_x();
      double ky = newtrack.kink_y();

      // Track impact position at Matching Point from Downstream:
      double xB = (*drip).getx_at(z_match);
      double yB = (*drip).gety_at(z_match);

      // Track impact position at Matching Point from Upstream:
      double xA = (*trip).getx_at(z_match);
      double yA = (*trip).gety_at(z_match);
      
      double dx = xB - xA; // driplet - triplet
      double dy = yB - yA;

      sixkxHisto->fill( kx*1E3 );
      sixkyHisto->fill( ky*1E3 );
      sixdxHisto->fill( dx );
      sixdyHisto->fill( dy );
      
      
      if( abs(dy) < 0.4 ) sixdxcHisto->fill( dx*1E3 ); // sig = 17 um at 5 GeV
      if( abs(dx) < 0.4 ) sixdycHisto->fill( dy*1E3 );

      // match driplet and triplet:
      if( abs(dx) > intersect_residual_cut) continue;
      if( abs(dy) > intersect_residual_cut) continue;
	
      sixkxcHisto->fill( kx*1E3 );
      sixkycHisto->fill( ky*1E3 );
      sixxHisto->fill( -xA ); // -xA = x_DP = out
      sixyHisto->fill( -yA ); // -yA = y_DP = up
      sixxyHisto->fill( -xA, -yA ); // DP: x_out, y_up

      // Fill kink map histogram:
      if( abs( kx ) > 0.002 || abs( ky ) > 0.002 ) sixxycHisto->fill( -xA, -yA );

      kinkvsxy->fill( -xA, -yA, (kx*kx + ky*ky)*1E6 ); //<kink^2> [mrad^2]

      // Add the track to the vector, the triplets are matched:
      tracks->push_back(newtrack);

    } // Downstream
  } // Upstream

  streamlog_out(MESSAGE4) << "Found " << tracks->size() << " tracks from matched triplets." << std::endl;
  trk_counter += tracks->size();
  return tracks;
}



EUTelTripletGBL::track::track(triplet up, triplet down) : upstream(up), downstream(down) {}

double EUTelTripletGBL::track::kink_x() {
  return (downstream.slope().x - upstream.slope().x);
}

double EUTelTripletGBL::track::kink_y() {
  return (downstream.slope().y - upstream.slope().y);
}

EUTelTripletGBL::hit EUTelTripletGBL::track::intersect() {
  hit inter;
  // Re-check what this actually is...
  // and simplifie using triplet class members...
  inter.x = ( upstream.base().x - upstream.slope().x * upstream.base().z - downstream.base().x + downstream.slope().x * downstream.base().z ) / kink_x();
  inter.y = ( upstream.base().y - upstream.slope().y * upstream.base().z - downstream.base().y + downstream.slope().y * downstream.base().z ) / kink_y();
  return inter;
}

EUTelTripletGBL::triplet EUTelTripletGBL::track::get_upstream() {
  return upstream;
}

EUTelTripletGBL::triplet EUTelTripletGBL::track::get_downstream() {
  return downstream;
}

EUTelTripletGBL::hit EUTelTripletGBL::track::gethit(int plane) {
  if(plane < 3) return upstream.gethit(plane);
  else return downstream.gethit(plane);
}



EUTelTripletGBL::triplet::triplet() : linked_dut(false), hits() {
  // Empty default constructor
}

EUTelTripletGBL::triplet::triplet(hit hit0, hit hit1, hit hit2) : linked_dut(false), hits() {
  triplet();
  filltriplet(hit0, hit1, hit2);
}

EUTelTripletGBL::hit EUTelTripletGBL::triplet::getpoint_at(double z) {
  hit impact;
  impact.z = z - base().z;
  impact.x = base().x + slope().x * impact.z;
  impact.y = base().y + slope().y * impact.z;
  return impact;
}

double EUTelTripletGBL::triplet::getx_at(double z) {
  return this->base().x + this->slope().x * (z - this->base().z);
}

double EUTelTripletGBL::triplet::getdx() {
  return this->hits.rbegin()->second.x - this->hits.begin()->second.x;
}

double EUTelTripletGBL::triplet::getdx(int ipl) {
  return this->hits[ipl].x - this->base().x - this->slope().x * (this->hits[ipl].z - this->base().z);
}

double EUTelTripletGBL::triplet::getdx(hit point) {
  return point.x - this->base().x - this->slope().x * (point.z - this->base().z);
}

double EUTelTripletGBL::triplet::gety_at(double z) {
  return this->base().y + this->slope().y * (z - this->base().z);
}

double EUTelTripletGBL::triplet::getdy() {
  return this->hits.rbegin()->second.y - this->hits.begin()->second.y;
}

double EUTelTripletGBL::triplet::getdy(int ipl) {
  return this->hits[ipl].y - this->base().y - this->slope().y * (this->hits[ipl].z - this->base().z);
}

double EUTelTripletGBL::triplet::getdy(hit point) {
  return point.y - this->base().y - this->slope().y * (point.z - this->base().z);
}

double EUTelTripletGBL::triplet::getdz() {
  return this->hits.rbegin()->second.z - this->hits.begin()->second.z;
}

EUTelTripletGBL::hit EUTelTripletGBL::triplet::gethit(int plane) {
  return this->hits[plane];
}

EUTelTripletGBL::hit EUTelTripletGBL::triplet::base() {
  hit center;
  center.x = 0.5*( this->hits.begin()->second.x + this->hits.rbegin()->second.x );
  center.y = 0.5*( this->hits.begin()->second.y + this->hits.rbegin()->second.y );
  center.z = 0.5*( this->hits.begin()->second.z + this->hits.rbegin()->second.z );
  return center;
}

EUTelTripletGBL::hit EUTelTripletGBL::triplet::slope() {
  hit sl;
  double dz = (this->hits.rbegin()->second.z - this->hits.begin()->second.z);
  sl.x = (this->hits.rbegin()->second.x - this->hits.begin()->second.x) / dz;
  sl.y = (this->hits.rbegin()->second.y - this->hits.begin()->second.y) / dz;
  return sl;
}

#endif
