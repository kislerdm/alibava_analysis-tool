// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
// Authors:
// Philipp Roloff, DESY <mailto:philipp.roloff@desy.de>
// Joerg Behr, Hamburg Uni/DESY <joerg.behr@desy.de>
// Slava Libov, DESY <mailto:vladyslav.libov@desy.de>
// Igor Rubinskiy, DESY <mailto:igorrubinsky@gmail.com>
// Daniel Pitzl, DESY <mailto:daniel.pitzl@desy.de>
// Dmitry Kisler, DESY <mailto:dmitry.kisler@desy.de>
//
// Version: $Id: EUTelMilleGBL.cc,v 1.48 2009-08-01 10:49:46 bulgheroni Exp $
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
// built only if GEAR and MARLINUTIL are used
//debug
#if defined(USE_GEAR) && defined(USE_MARLINUTIL)

// eutelescope includes ".h"
#include "EUTelMilleGBL.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTELESCOPE.h"
#include "EUTelVirtualCluster.h"
#include "EUTelFFClusterImpl.h"
#include "EUTelDFFClusterImpl.h"
#include "EUTelBrickedClusterImpl.h"
#include "EUTelSparseClusterImpl.h"
#include "EUTelExceptions.h"
#include "EUTelPStream.h" // process streams redi::ipstream
#include "EUTelAlignmentConstant.h"

// GBL:

#include "include/GblTrajectory.h"
#include "include/MilleBinary.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Global.h"
#include "marlin/Exceptions.h"
#include "marlin/AIDAProcessor.h"

// marlin util includes
#include "mille/Mille.h"

// gear includes <.h>
#include <gear/GearMgr.h>
#include <gear/SiPlanesParameters.h>

// aida includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/ITree.h>
#endif

// lcio includes <.h>
#include <IO/LCWriter.h>
#include <UTIL/LCTime.h>
#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <Exceptions.h>
#include <UTIL/CellIDEncoder.h>
#include <IMPL/TrackerRawDataImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/TrackerPulseImpl.h>
#include <IMPL/LCCollectionVec.h>

// ROOT includes
#if defined(USE_ROOT) || defined(MARLIN_USE_ROOT)
#include <TRandom.h>
//#include <TMinuit.h>
#include <TSystem.h>
#include <TMath.h>
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TVector3.h>
#include <TRotation.h>
#endif

// system includes <>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <memory>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace eutelescope;

#if defined(USE_ROOT) || defined(MARLIN_USE_ROOT)
//global:
//EUTelMilleGBL::hit *hitsarray_gbl;
#endif

// definition of static members mainly used to name histograms
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

AIDA::IHistogram1D * dx01Hist;
AIDA::IHistogram1D * dy01Hist;
AIDA::IHistogram1D * dx02Hist;
AIDA::IHistogram1D * dy02Hist;
AIDA::IHistogram1D * dx03Hist;
AIDA::IHistogram1D * dy03Hist;
AIDA::IHistogram1D * dx04Hist;
AIDA::IHistogram1D * dy04Hist;
AIDA::IHistogram1D * dx05Hist;
AIDA::IHistogram1D * dy05Hist;

AIDA::IHistogram1D * tridxHist;
AIDA::IHistogram1D * tridyHist;
AIDA::IHistogram1D * tridx_noCut_Hist;
AIDA::IHistogram1D * tridy_noCut_Hist;
AIDA::IHistogram1D * tridr_noCut_Hist;
AIDA::IHistogram1D * ntriHist;

AIDA::IHistogram1D * dridxHist;
AIDA::IHistogram1D * dridyHist;
AIDA::IHistogram1D * dridx_noCut_Hist;
AIDA::IHistogram1D * dridy_noCut_Hist;
AIDA::IHistogram1D * dridr_noCut_Hist;
AIDA::IHistogram1D * ndriHist;

AIDA::IHistogram1D * sixdxHist;
AIDA::IHistogram1D * sixdyHist;
AIDA::IHistogram1D * sixdx_noCut_Hist;
AIDA::IHistogram1D * sixdy_noCut_Hist;
AIDA::IHistogram1D * sixdr_noCut_Hist;
AIDA::IHistogram1D * sixkxHist;
AIDA::IHistogram1D * sixkyHist;

AIDA::IHistogram1D * selxHist;
AIDA::IHistogram1D * selyHist;
AIDA::IHistogram1D * selaxHist;
AIDA::IHistogram1D * selayHist;
AIDA::IHistogram1D * seldxHist;
AIDA::IHistogram1D * seldyHist;
AIDA::IHistogram1D * selkxHist;
AIDA::IHistogram1D * selkyHist;

AIDA::IHistogram1D * seldx1Hist;
AIDA::IHistogram1D * seldy1Hist;
AIDA::IHistogram1D * seldx3Hist;
AIDA::IHistogram1D * seldy3Hist;
AIDA::IHistogram1D * seldx4Hist;
AIDA::IHistogram1D * seldy4Hist;
AIDA::IHistogram1D * seldx5Hist;
AIDA::IHistogram1D * seldy5Hist;
AIDA::IHistogram1D * seldx6Hist;
AIDA::IHistogram1D * seldy6Hist;

AIDA::IHistogram1D * gblndfHist;
AIDA::IHistogram1D * gblchi2Hist;
AIDA::IHistogram1D * gblprbHist;

AIDA::IHistogram1D * badxHist;
AIDA::IHistogram1D * badyHist;
AIDA::IHistogram1D * badaxHist;
AIDA::IHistogram1D * badayHist;
AIDA::IHistogram1D * baddxHist;
AIDA::IHistogram1D * baddyHist;
AIDA::IHistogram1D * badkxHist;
AIDA::IHistogram1D * badkyHist;

AIDA::IHistogram1D * baddx1Hist;
AIDA::IHistogram1D * baddy1Hist;
AIDA::IHistogram1D * baddx3Hist;
AIDA::IHistogram1D * baddy3Hist;
AIDA::IHistogram1D * baddx4Hist;
AIDA::IHistogram1D * baddy4Hist;
AIDA::IHistogram1D * baddx5Hist;
AIDA::IHistogram1D * baddy5Hist;
AIDA::IHistogram1D * baddx6Hist;
AIDA::IHistogram1D * baddy6Hist;

AIDA::IHistogram1D * goodx1Hist;
AIDA::IHistogram1D * goody1Hist;
AIDA::IHistogram1D * goodx6Hist;
AIDA::IHistogram1D * goody6Hist;

AIDA::IHistogram1D * gblax0Hist;
AIDA::IHistogram1D * gbldx0Hist;
AIDA::IHistogram1D * gblrx0Hist;

AIDA::IHistogram1D * gblax1Hist;
AIDA::IHistogram1D * gbldx1Hist;
AIDA::IHistogram1D * gblrx1Hist;

AIDA::IHistogram1D * gblax2Hist;
AIDA::IHistogram1D * gbldx2Hist;
AIDA::IHistogram1D * gblrx2Hist;

AIDA::IHistogram1D * gblax3Hist;
AIDA::IHistogram1D * gbldx3Hist;
AIDA::IHistogram1D * gblrx3Hist;

AIDA::IHistogram1D * gblax4Hist;
AIDA::IHistogram1D * gbldx4Hist;
AIDA::IHistogram1D * gblrx4Hist;

AIDA::IHistogram1D * gblax5Hist;
AIDA::IHistogram1D * gbldx5Hist;
AIDA::IHistogram1D * gblrx5Hist;

AIDA::IHistogram1D * gblax6Hist;
AIDA::IHistogram1D * gbldx6Hist;
AIDA::IHistogram1D * gbldy6Hist;
AIDA::IHistogram1D * gblrx6Hist;
AIDA::IHistogram1D * gblry6Hist;

AIDA::IHistogram1D * gblkx1Hist;
AIDA::IHistogram1D * gblkx2Hist;
AIDA::IHistogram1D * gblkx3Hist;
AIDA::IHistogram1D * gblkx4Hist;
AIDA::IHistogram1D * gblkx5Hist;
AIDA::IHistogram1D * gblkx6Hist;

AIDA::IHistogram1D * nmHist;

//How many tracks (matched triplests) were found per event when the hits on DUT are presented
AIDA::IHistogram1D * hTripletsMatched_DUT;
AIDA::IHistogram2D * hTripletsMatched_DUT_2D;
//How many tracks (matched triplests) were found per event when the hits on REF plane are presented
AIDA::IHistogram1D * hTripletsMatched_REF;
AIDA::IHistogram2D * hTripletsMatched_REF_2D;
//How many tracks (matched triplests) were found per event when the hits on DUT and REF plane are presented
AIDA::IHistogram1D * hTripletsMatched_DUT_REF;
AIDA::IHistogram2D * hTripletsMatched_DUT_REF_2D;

//track matching to REF
AIDA::IHistogram1D * hTrk_REF;
AIDA::IHistogram1D * hTrk_REF_dx_noCut;
AIDA::IHistogram1D * hTrk_REF_dy_noCut;
AIDA::IHistogram1D * hTrk_REF_dr_noCut;
AIDA::IHistogram1D * hTrk_REF_dx;
AIDA::IHistogram1D * hTrk_REF_dy;
AIDA::IHistogram1D * hTrk_REF_dr;
//track matching to DUT
AIDA::IHistogram1D * hTrk_DUT;

#endif

gbl::MilleBinary * milleGBL; // for producing MillePede-II binary file

//------------------------------------------------------------------------------
EUTelMilleGBL::EUTelMilleGBL(): Processor("EUTelMilleGBL") {

    //some default values
    FloatVec MinimalResidualsX;
    FloatVec MinimalResidualsY;
    FloatVec MaximalResidualsX;
    FloatVec MaximalResidualsY;

    FloatVec PedeUserStartValuesX;
    FloatVec PedeUserStartValuesY;

    FloatVec PedeUserStartValuesGamma;

    FloatVec SensorZPositions;

    FloatVec SensorXShifts;
    FloatVec SensorYShifts;

    FloatVec SensorGamma;

    FloatVec SensorAlpha;

    FloatVec SensorBeta;

    //maybe one has to chose a larger value than 6?

    //  for( int i = 0; i < 8; i++ ) {

    //    MinimalResidualsX.push_back(0.0);
    //    MinimalResidualsY.push_back(0.0);
    //    MaximalResidualsX.push_back(0.0);
    //    MaximalResidualsY.push_back(0.0);

    //    PedeUserStartValuesX.push_back(0.0);
    //    PedeUserStartValuesY.push_back(0.0);

    //    PedeUserStartValuesGamma.push_back(0.0);

    //    float zpos = 20000.0 +  20000.0 * (float)i; // [um]
    //    SensorZPositions.push_back(zpos);

    //    SensorXShifts.push_back(0.0);
    //    SensorYShifts.push_back(0.0);

    //    SensorGamma.push_back(0.0);
    //    SensorAlpha.push_back(0.0);
    //    SensorBeta.push_back(0.0);
    //  }

    // modify processor description
    _description =
            "EUTelMilleGBL uses the MILLE program to write data files for MILLEPEDE II.";

    // choose input mode
    registerOptionalParameter("InputMode","Selects the source of input hits."
                              "\n0 - hits read from hitfile with simple trackfinding. "
                              "\n1 - hits read from output of tracking processor. "
                              "\n2 - Test mode. Simple internal simulation and simple trackfinding. "
                              "\n3 - Mixture of a track collection from the telescope and hit collections for the DUT (only one DUT layer can be used unfortunately)",
                              _inputMode, static_cast <int> (0));

    // input collections
    std::vector<std::string > HitCollectionNameVecExample;
    HitCollectionNameVecExample.push_back("corrhits");

    registerInputCollections(LCIO::TRACKERHIT,"HitCollectionName",
                             "Hit collections name",
                             _hitCollectionName,HitCollectionNameVecExample);

    registerInputCollection(LCIO::TRACK,"TrackCollectionName",
                            "Track collection name",
                            _trackCollectionName,std::string("fittracks"));

    // parameters


    registerProcessorParameter( "EventMaxNumb",
                                "Max number of events to be processed",
                                _nmax, int( 500000));

    registerProcessorParameter( "Ebeam",
                                "Beam energy [GeV]",
                                _eBeam, static_cast < double >( 4.0));

    registerOptionalParameter("DistanceMax","Maximal allowed distance between hits entering the fit per 10 cm space between the planes.",
                              _distanceMax, static_cast <float> (2000.0));

    registerOptionalParameter("DistanceMaxVec","Maximal allowed distance between hits entering the fit per 10 cm space between the planes. One value for each neighbor planes. DistanceMax will be used for each pair if this vector is empty.",
                              _distanceMaxVec, std::vector<float> ());

    registerOptionalParameter("ExcludePlanes","Exclude planes from fit according to their sensor ids.",_excludePlanes_sensorIDs ,std::vector<int>());

    registerOptionalParameter("FixedPlanes","Fix sensor planes in the fit according to their sensor ids.",_FixedPlanes_sensorIDs ,std::vector<int>());


    registerOptionalParameter("MaxTrackCandidatesTotal","Maximal number of track candidates (Total).",_maxTrackCandidatesTotal, static_cast <int> (1000000));
    registerOptionalParameter("MaxTrackCandidates","Maximal number of track candidates.",_maxTrackCandidates, static_cast <int> (4000));

    registerOptionalParameter("BinaryFilename","Name of the Millepede binary file.",_binaryFilename, string ("mille.bin"));

    registerOptionalParameter("TelescopeResolution","Resolution of the telescope for Millepede (sigma_x=sigma_y).",_telescopeResolution, static_cast <float> (3.0));

    registerOptionalParameter("DUTResolution","Resolution of the DUT for Millepede (sigma_y/sigma_x, if y/x is the activeaxis).",_DUTResolution, static_cast <float> (10.0));

    registerOptionalParameter("REFResolution","Resolution of the REF plane for Millepede (sigma_x=sigma_y).",_REFResolution, static_cast <float> (10.0));

    registerOptionalParameter("OnlySingleHitEvents","Use only events with one hit in every plane.",_onlySingleHitEvents, static_cast <int> (0));

    registerOptionalParameter("OnlySingleTrackEvents","Use only events with one track candidate.",_onlySingleTrackEvents, static_cast <int> (0));

    registerOptionalParameter("AlignMode","Number of alignment constants used. Available mode are: "
                              "\n2 - shifts in X and Y"
                              "\n3 - shifts in X and Y and rotation around the Z axis,"
                              "\n4 - shifts in X,Y and Z and rotation around the Z axis",
                              _alignMode, static_cast <int> (3));

    registerOptionalParameter("UseResidualCuts","Use cuts on the residuals to reduce the combinatorial background. 0 for off (default), 1 for on",_useResidualCuts,
                              static_cast <int> (0));

    registerOptionalParameter("AlignmentConstantLCIOFile","This is the name of the LCIO file name with the output alignment"
                              "constants (add .slcio)",_alignmentConstantLCIOFile, static_cast< string > ( "alignment.slcio" ) );

    registerOptionalParameter("AlignmentConstantCollectionName", "This is the name of the alignment collection to be saved into the slcio file",
                              _alignmentConstantCollectionName, static_cast< string > ( "alignment" ));

    registerOptionalParameter("ResidualsXMin","Minimal values of the hit residuals in the X direction for a track. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.",_residualsXMin,MinimalResidualsX);

    registerOptionalParameter("ResidualsYMin","Minimal values of the hit residuals in the Y direction for a track. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.",_residualsYMin,MinimalResidualsY);

    registerOptionalParameter("ResidualsXMax","Maximal values of the hit residuals in the X direction for a track. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.",_residualsXMax,MaximalResidualsX);

    registerOptionalParameter("ResidualsYMax","Maximal values of the hit residuals in the Y direction for a track. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.",_residualsYMax,MaximalResidualsY);

    registerOptionalParameter("ResolutionX","X resolution parameter for each plane. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.",_resolutionX,  std::vector<float> (static_cast <int> (6), 10.));

    registerOptionalParameter("ResolutionY","Y resolution parameter for each plane. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.",_resolutionY,std::vector<float> (static_cast <int> (6), 10.));

    registerOptionalParameter("ResolutionZ","Z resolution parameter for each plane. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.",_resolutionZ,std::vector<float> (static_cast <int> (6), 10.));

    registerOptionalParameter("FixParameter","Fixes the given alignment parameters in the fit if alignMode==6 is used. For each sensor an integer must be specified (If no value is given, then all parameters will be free). bit 0 = x shift, bit 1 = y shift, bit 2 = z shift, bit 3 = alpha, bit 4 = beta, bit 5 = gamma. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.",_FixParameter, std::vector<int> (static_cast <int> (6), 24));

    registerOptionalParameter("GeneratePedeSteerfile","Generate a steering file for the pede program.",_generatePedeSteerfile, static_cast <int> (1));

    registerOptionalParameter("PedeSteerfileName","Name of the steering file for the pede program.",_pedeSteerfileName, string("steer_mille.txt"));

    registerOptionalParameter("RunPede","Execute the pede program using the generated steering file.",_runPede, static_cast <int> (0));

    registerOptionalParameter("UsePedeUserStartValues","Give start values for pede by hand (0 - automatic calculation of start values, 1 - start values defined by user).", _usePedeUserStartValues, static_cast <int> (0));

    registerOptionalParameter("PedeUserStartValuesX","Start values for the alignment for shifts in the X direction.",_pedeUserStartValuesX,PedeUserStartValuesX);

    registerOptionalParameter("PedeUserStartValuesY","Start values for the alignment for shifts in the Y direction.",_pedeUserStartValuesY,PedeUserStartValuesY);

    registerOptionalParameter("PedeUserStartValuesZ","Start values for the alignment for shifts in the Z direction.",_pedeUserStartValuesZ,std::vector<float> (static_cast <int> (6), 0.0));

    registerOptionalParameter("PedeUserStartValuesAlpha","Start values for the alignment for the angle alpha.",_pedeUserStartValuesAlpha,std::vector<float> (static_cast <int> (6), 0.0));

    registerOptionalParameter("PedeUserStartValuesBeta","Start values for the alignment for the angle beta.",_pedeUserStartValuesBeta,std::vector<float> (static_cast <int> (6), 0.0));

    registerOptionalParameter("PedeUserStartValuesGamma","Start values for the alignment for the angle gamma.",_pedeUserStartValuesGamma,PedeUserStartValuesGamma);

    registerOptionalParameter("TestModeSensorResolution","Resolution assumed for the sensors in test mode.",_testModeSensorResolution, static_cast <float> (3.0));

    registerOptionalParameter("TestModeXTrackSlope","Width of the track slope distribution in the x direction",_testModeXTrackSlope, static_cast <float> (0.0005));

    registerOptionalParameter("TestModeYTrackSlope","Width of the track slope distribution in the y direction",_testModeYTrackSlope, static_cast <float> (0.0005));

    registerOptionalParameter("TestModeSensorZPositions","Z positions of the sensors in test mode.",_testModeSensorZPositions,SensorZPositions);

    registerOptionalParameter("TestModeSensorXShifts","X shifts of the sensors in test mode (to be determined by the alignment).",
                              _testModeSensorXShifts,SensorXShifts);

    registerOptionalParameter("TestModeSensorYShifts","Y shifts of the sensors in test mode (to be determined by the alignment).",
                              _testModeSensorYShifts,SensorYShifts);

    registerOptionalParameter("TestModeSensorGamma","Rotation around the z axis of the sensors in test mode (to be determined by the alignment).",
                              _testModeSensorGamma,SensorGamma);

    registerOptionalParameter("TestModeSensorAlpha","Rotation around the x axis of the sensors in test mode (to be determined by the alignment).",
                              _testModeSensorAlpha,SensorAlpha);

    registerOptionalParameter("TestModeSensorBeta","Rotation around the y axis of the sensors in test mode (to be determined by the alignment).",
                              _testModeSensorBeta,SensorBeta);

    std::vector<int> initRect;
    registerOptionalParameter("UseSensorRectangular","Do not use all pixels for alignment, only these in the rectangular (A|B) e.g. (0,0) and (C|D) e.g. (100|100) of sensor S. Type in the way S1 A1 B1 C1 D1 S2 A2 B2 C2 D2 ...",
                              _useSensorRectangular,initRect);

    registerOptionalParameter("DUTSensetiveAxis", "What is the DUT sensetive axis", _sensetiveAxis, string("y"));

    registerOptionalParameter("DUTscatterer", "Shall the DUT be takken to account as a scatterer? (1-yes, 0-no)", DUT_scatterer, int(0));

    registerOptionalParameter("DUTX0", "Radiation length of the DUT material", DUTX0, double(0.08));

    registerOptionalParameter("Chi2ToNdfCut", "The Chi2/Ndf cut", _chi2tondfcut, double(0.01));

    //CUTS
    registerOptionalParameter( "triCut", "Upstream triplet residual cut [um]", _triCut, double (300.0) );
    registerOptionalParameter( "driCut", "Downstream triplet residual cut [um]", _driCut, double (400.0) );
    registerOptionalParameter( "sixCut", "Upstream-Downstream Track matching cut [um]", _sixCut, double (600.0) );
    registerOptionalParameter( "TriDutCut", "Upstream Track - DUT matching cut [um]", _tridutCut, double (50.0) );
    registerOptionalParameter( "DriRefCut", "Downstream Track - REF matching cut [um]", _drirefCut, double (50.0) );

    //Tel of DUT alignment? (assume that tel planes are fixed if DUT and REF are being alligned) 0 (default) - tel align, 1 - REF align, 2 - DUT align
    registerOptionalParameter( "AlignType", "Type of alignment procedure", _aligntype, int (0) );

}


//------------------------------------------------------------------------------
void EUTelMilleGBL::init() {
    // check if Marlin was built with GEAR support or not
#ifndef USE_GEAR

    streamlog_out( ERROR2 ) << "Marlin was not built with GEAR support." << endl;
    streamlog_out( ERROR2 ) << "You need to install GEAR and recompile Marlin with -DUSE_GEAR before continue." << endl;

    exit(-1);

#else

    // check if the GEAR manager pointer is not null!
    if( Global::GEAR == 0x0 ) {
        streamlog_out( ERROR2) << "The GearMgr is not available, for an unknown reason." << endl;
        exit(-1);
    }

    _siPlanesParameters  = const_cast<gear::SiPlanesParameters* > (&(Global::GEAR->getSiPlanesParameters()));
    _siPlanesLayerLayout = const_cast<gear::SiPlanesLayerLayout*> ( &(_siPlanesParameters->getSiPlanesLayerLayout() ));

    _histogramSwitch = true;

    streamlog_out(MESSAGE4) << "assumed beam energy " << _eBeam << " GeV" <<  endl;
    //lets guess the number of planes
    if( _inputMode == 0 || _inputMode == 2 ) { // we use inputMode 0

        // the number of planes is got from the GEAR description and is
        // the sum of the telescope reference planes and the DUT (if
        // any)
        _nPlanes = _siPlanesLayerLayout->getNLayers();
        streamlog_out (MESSAGE4) << "n planes = " << _nPlanes <<endl;

        if( _useSensorRectangular.size() == 0 ) {
            streamlog_out(MESSAGE4) << "No rectangular limits on pixels of sensorplanes applied" << endl;
        }
        else {
            if( _useSensorRectangular.size() % 5 != 0) {
                streamlog_out(WARNING2) << "Wrong number of arguments in RectangularLimits! Ignoring this cut!" << endl;
            }
            else {
                streamlog_out(MESSAGE4) << "Reading in SensorRectangularCuts: " << endl;
                int sensorcuts = _useSensorRectangular.size()/5;
                for( int i = 0; i < sensorcuts; ++i) {
                    int sensor = _useSensorRectangular.at(5*i+0);
                    int A = _useSensorRectangular.at(5*i+1);
                    int B = _useSensorRectangular.at(5*i+2);
                    int C = _useSensorRectangular.at(5*i+3);
                    int D = _useSensorRectangular.at(5*i+4);
                    SensorRectangular r(sensor,A,B,C,D);
                    r.print();
                    _rect.addRectangular(r);
                }
            }
        }

    }
    else if( _inputMode == 1 ) {
        _nPlanes = _siPlanesParameters->getSiPlanesNumber();
    }
    else if(_inputMode == 3) {

        // the number of planes is got from the GEAR description and is
        // the sum of the telescope reference planes and the DUT (if
        // any)
        _nPlanes = _siPlanesParameters->getSiPlanesNumber();
        if( _siPlanesParameters->getSiPlanesType() == _siPlanesParameters->TelescopeWithDUT ) {
            ++_nPlanes;
        }
    }
    else {
        streamlog_out (MESSAGE4) << "unknown input mode " << _inputMode << endl;
        exit(-1);
    }

    // an associative map for getting also the sensorID ordered
    map< double, int > sensorIDMap;

    //lets create an array with the z positions of each layer
    for(  int iPlane = 0 ; iPlane < _nPlanes; iPlane++ )
    {
        _siPlaneZPosition.push_back(_siPlanesLayerLayout->getLayerPositionZ(iPlane));
        sensorIDMap.insert( make_pair( _siPlanesLayerLayout->getLayerPositionZ(iPlane), _siPlanesLayerLayout->getID(iPlane) ) );
    }

    //lets sort the array with increasing z
    //sort(_siPlaneZPosition.begin(), _siPlaneZPosition.end());

    streamlog_out (DEBUG1) << "planes sorted along z, tel planes come first, then DUT, REF is the last :" << endl;
    for( size_t i = 0; i < _siPlaneZPosition.size(); i++ ) {
        _orderedSensorID_wo_excluded.push_back(_siPlanesLayerLayout->getID(i));
        streamlog_out (DEBUG1) << i << "geomID = " << _orderedSensorID_wo_excluded.at(i) << "  Z [mm] = " << _siPlaneZPosition[i] << endl;
    }

    //during the DUT and REF alignment procedure, only the DUT and REF planes are being aligned, the rest stays fixed
    if (_aligntype > 0)
    {
        for(int iExcl = 0; iExcl < 6; iExcl++) {_FixedPlanes_sensorIDs[iExcl] = iExcl;}
    }

    streamlog_out (DEBUG1) << "FixedPlanes " << _FixedPlanes_sensorIDs.size() << endl;

    for( size_t i = 0; i < _FixedPlanes_sensorIDs.size(); i++ ) {
        streamlog_out (DEBUG1) << "plane " << _FixedPlanes_sensorIDs[i] << " fixed\n";
        map< double, int >::iterator iter = sensorIDMap.begin();
        int counter = 0;
        while ( iter != sensorIDMap.end() ) {
            if( iter->second == _FixedPlanes_sensorIDs[i] ) {
                _FixedPlanes.push_back(counter);
                break;
            }
            ++iter;
            ++counter;
        }
    }
    //    int nexcl_pl = _excludePlanes_sensorIDs.size();
    //    streamlog_out (MESSAGE4) << "ExcludedPlanes " << nexcl_pl << endl;
    //    for( size_t i = 0; i < _excludePlanes_sensorIDs.size(); i++ ) {
    //        streamlog_out (MESSAGE4) << "plane " << _excludePlanes_sensorIDs[i] << " excluded\n";
    //        map< double, int >::iterator iter = sensorIDMap.begin();
    //        int counter = 0;
    //        while ( iter != sensorIDMap.end() ) {
    //            if( iter->second == _excludePlanes_sensorIDs[i]) {
    //                _excludePlanes.push_back(counter);
    //                break;
    //            }
    //            ++iter;
    //            ++counter;
    //        }
    //    }

    //    // strip from the map the sensor id already sorted.
    //    map< double, int >::iterator iter = sensorIDMap.begin();
    //    int counter = 0;
    //    while ( iter != sensorIDMap.end() ) {

    //        bool excluded = false;
    //        for( size_t i = 0; i < _excludePlanes.size(); i++) {
    //            if(_excludePlanes[i] == counter) {
    //                // printf("excludePlanes %2d of %2d (%2d) \n", i, _excludePlanes_sensorIDs.size(), counter );
    //                excluded = true;
    //                break;
    //            }
    //        }
    //        if(!excluded) {_orderedSensorID_wo_excluded.push_back( iter->second );}

    //        _orderedSensorID.push_back( iter->second );

    //        ++iter;
    //        ++counter;
    //    }
    //check
    //    for(int i = 0; i < _nPlanes; i++)
    //    {
    //        streamlog_out(DEBUG0) << "Check the _orderedSensorID_wo_excluded GeomID order"<<endl
    //                                 <<"counter = " << i << "GearID = " << _orderedSensorID_wo_excluded.at(i) << endl;
    //        streamlog_out(DEBUG0) << "Check the _orderedSensorID order"<<endl
    //                                 <<"counter = " << i << "GearID = " << _orderedSensorID.at(i) << endl;
    //    }

    //consistency
    //  if( ( (int)_siPlaneZPosition.size() - nexcl_pl ) !=  _nPlanes   ) {
    //    streamlog_out( ERROR2 ) << "the number of detected planes (after exclusion) is " << _nPlanes << " but only " << _siPlaneZPosition.size() << " layer z positions were found!"  << endl;
    //    exit(-1);
    //  }

#endif

    // this method is called only once even when the rewind is active
    // usually a good idea to
    printParameters ();

    // set to zero the run and event counters
    _iRun = 0;
    _iEvt = 0;

    // Initialize number of excluded planes
    _nExcludePlanes = _excludePlanes.size();

    streamlog_out( MESSAGE2 ) << "Number of planes excluded from the alignment fit: " << _nExcludePlanes << endl;

    // Initialise Mille statistics:
    _nMilleTracksperEvent = 0;
    _nMilleTracksTotal = 0;
    ntrk_REF_Total = 0;
    ntrk_REF_Event = 0;
    ntrk_DUT_Total = 0;
    ntrk_DUT_Event = 0;


    _waferResidX = new double[_nPlanes];
    _waferResidY = new double[_nPlanes];
    _waferResidZ = new double[_nPlanes];


    _xFitPos = new double[_nPlanes];
    _yFitPos = new double[_nPlanes];

    _telescopeResolX = new double[_nPlanes];
    _telescopeResolY = new double[_nPlanes];
    _telescopeResolZ = new double[_nPlanes];

    //  printf("print resolution  X: %2d, Y: %2d, Z: %2d \n", _resolutionX.size(), _resolutionY.size(), _resolutionZ.size() );

    // booking histograms
    bookHistos();

    streamlog_out( MESSAGE2 ) << "Initialising Mille..." << endl;

    unsigned int reserveSize = 800000;
    milleGBL = new gbl::MilleBinary( _binaryFilename, reserveSize );

    streamlog_out( MESSAGE2 ) << "The filename for the binary file is: " << _binaryFilename.c_str() << endl;

    for(int i = 0; i < _maxTrackCandidates; i++) {
        _xPos.push_back(std::vector<double>(_nPlanes,0.0));
        _yPos.push_back(std::vector<double>(_nPlanes,0.0));
        _zPos.push_back(std::vector<double>(_nPlanes,0.0));
    }

    if( _distanceMaxVec.size() > 0 ) {
        if(_distanceMaxVec.size() !=  static_cast<unsigned int>(_nPlanes )-1 ) {
            streamlog_out( WARNING2 ) << "Consistency check of the DistanceMaxVec array failed. Its size is different compared to the number of planes! Will now use _distanceMax for each pair of planes." << endl;
            _distanceMaxVec.clear();
            for(int i = 0; i < _nPlanes-1; i++) {
                _distanceMaxVec.push_back(_distanceMax);
            }
        }
    }
    else {
        _distanceMaxVec.clear();
        for(int i = 0; i < _nPlanes-1; i++) {
            _distanceMaxVec.push_back(_distanceMax);
        }
    }
    streamlog_out( MESSAGE2 ) << "end of init" << endl;
}


//------------------------------------------------------------------------------
void EUTelMilleGBL::processRunHeader( LCRunHeader * rdr ) {

    auto_ptr<EUTelRunHeaderImpl> header ( new EUTelRunHeaderImpl (rdr) );
    header->addProcessor( type() ) ;

    // this is the right place also to check the geometry ID. This is a
    // unique number identifying each different geometry used at the
    // beam test. The same number should be saved in the run header and
    // in the xml file. If the numbers are different, instead of barely
    // quitting ask the user what to do.

    if( header->getGeoID() != _siPlanesParameters->getSiPlanesID() ) {
        streamlog_out( ERROR2 ) << "Error during the geometry consistency check: " << endl;
        streamlog_out( ERROR2 ) << "The run header says the GeoID is " << header->getGeoID() << endl;
        streamlog_out( ERROR2 ) << "The GEAR description says is     " << _siPlanesParameters->getSiPlanesNumber() << endl;

#ifdef EUTEL_INTERACTIVE
        string answer;
        while (true) {
            streamlog_out( ERROR2 ) << "Type Q to quit now or C to continue using the actual GEAR description anyway [Q/C]" << endl;
            cin >> answer;
            // put the answer in lower case before making the comparison.
            transform( answer.begin(), answer.end(), answer.begin(), ::tolower );
            if( answer == "q" ) {
                exit(-1);
            } else if( answer == "c" ) {
                break;
            }
        }
#endif

    }

    // increment the run counter
    ++_iRun;
}

//------------------------------------------------------------------------------
TMatrixD Jac55( double ds ) {
    /* for GBL:
                                        Jacobian for straight line track
                                        track = q/p, x', y', x, y
                                                0,   1,  2,  3, 4
                                     */
    TMatrixD jac(5, 5);
    jac.UnitMatrix();
    jac[3][1] = ds; // x = x0 + xp * ds
    jac[4][2] = ds; // y = y0 + yp * ds
    return jac;
}


//------------------------------------------------------------------------------
void EUTelMilleGBL::processEvent( LCEvent * event ) {

    if( _iEvt % 1000 == 0 ) {
        streamlog_out( MESSAGE4 ) << "Processing event "
                                  << setw(6) << setiosflags(ios::right)
                                  << event->getEventNumber() << " in run "
                                  << setw(6) << setiosflags(ios::right)
                                  << event->getRunNumber()
                                  << ", currently having "
                                  << _nMilleTracksTotal << " tracks "
                                  << endl;
    }

    if(_iEvt > _nmax)
    {
        streamlog_out(MESSAGE4) << "Max set event number was reached" << endl;
        return;
    }

    if( _nMilleTracksTotal > _maxTrackCandidatesTotal ) {
        throw StopProcessingException(this);
    }

    // fill resolution arrays
    for( int help = 0; help < _nPlanes; help++ ) {
        _telescopeResolX[help] = _telescopeResolution;
        _telescopeResolY[help] = _telescopeResolution;
        if(_siPlanesLayerLayout->getID(help) == 7) { //dut resolution
            if(_sensetiveAxis == "y")
            {
                _telescopeResolY[help] = _DUTResolution;
                _telescopeResolX[help] = _telescopeResolution;
            }
            if(_sensetiveAxis == "x")
            {
                _telescopeResolX[help] = _DUTResolution;
                _telescopeResolY[help] = _telescopeResolution;
            }
        }
        if(_siPlanesLayerLayout->getID(help) == 8) { //ref resolution
            _telescopeResolX[help] = _REFResolution;
            _telescopeResolY[help] = _REFResolution;}
    }

    EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event) ;

    if( evt->getEventType() == kEORE ) {
        streamlog_out( DEBUG2 ) << "EORE found: nothing else to do." << endl;
        return;
    }

    std::vector<std::vector<EUTelMilleGBL::HitsInPlane> > _hitsArray(_nPlanes - _nExcludePlanes, std::vector<EUTelMilleGBL::HitsInPlane>() );

    std::vector<int> indexconverter (_nPlanes,-1);
    streamlog_out (DEBUG0) << "checking index converter : " << endl;
    for( int i = 0; i < _nPlanes; i++ ) {

        if(_siPlanesLayerLayout->getID(i) < 6){ indexconverter[i] = i; }
        else if(_siPlanesLayerLayout->getID(i) == 7) { indexconverter[i] = _nPlanes-1; } //dut
        else if(_siPlanesLayerLayout->getID(i) == 8) { indexconverter[i] = _nPlanes-2; } //ref

       streamlog_out (DEBUG0)  << "indexconverter(i=" <<i<< ") : " << indexconverter[i] <<
                                  " geomID = " << _siPlanesLayerLayout->getID(i) <<endl;
    }



    //  printf("887: excludePlanes %2d  \n", _nExcludePlanes );
    //  if( _nExcludePlanes > 0 )
    //  {
//    int icounter = 0;
//    for( int i = 0; i < _nPlanes; i++ ) {
//        int excluded = 0; //0 - not excluded, 1 - excluded
//        if( _nExcludePlanes > 0 ) {
//            for( int helphelp = 0; helphelp < _nExcludePlanes; helphelp++ ) {
//                if( i == _excludePlanes[helphelp] ) {
//                    excluded = 1;
//                    break;//leave the for loop
//                }
//            }
//        }
//        if( excluded == 1 )
//            indexconverter[i] = -1;
//        else {
//            indexconverter[i] = icounter;
//            icounter++;
//        }
//    }
    //  }

    if( _inputMode != 1 && _inputMode != 3 ) // we use _inputMode 0

        for( size_t i = 0; i < _hitCollectionName.size(); i++ ) {

            LCCollection* collection;
            try {
                collection = event->getCollection(_hitCollectionName[i]);
            } catch (DataNotAvailableException& e) {
                //	streamlog_out( WARNING2 ) << "No input collection "
                //			  << _hitCollectionName[i] << " found for event "
                //			  << event->getEventNumber()
                //			  << " in run " << event->getRunNumber() << endl;
                throw SkipEventException(this);
            }
            int layerIndex = -1;
            //array of hits on plane
            HitsInPlane hitsInPlane;

            // check if running in input mode 0 or 2
            if( _inputMode == 0 ) {

                // loop over all hits in collection:

                for( int iHit = 0; iHit < collection->getNumberOfElements(); iHit++ ) {

                    //taking a hit from the hit collection(s)
                    TrackerHitImpl * hit = static_cast<TrackerHitImpl*> ( collection->getElementAt(iHit) );

                    LCObjectVec clusterVector = hit->getRawHits();

                    double minDistance =  numeric_limits< double >::max();
                    double * hitPosition = const_cast<double * > (hit->getPosition());

                    //forming the vectors of hits corresponding to each plane
                    for(  int i = 0; i < (int)_siPlaneZPosition.size(); i++ ) {
                        double distance = std::abs( hitPosition[2] - _siPlaneZPosition[i] );
                        if( distance < minDistance ) {
                            minDistance = distance;
                            layerIndex = i;
                        }
                    }
                    if( minDistance > 30 /* mm */ ) {
                        // advise the user that the guessing wasn't successful
                        streamlog_out( WARNING3 ) << "A hit was found " << minDistance << " mm far from the nearest plane\n"
                                                     "Please check the consistency of the data with the GEAR file" << endl;
                    }

                    // Getting positions of the hits.
                    // ------------------------------
                    hitsInPlane.measuredX = 1000 * hit->getPosition()[0]; // microns
                    hitsInPlane.measuredY = 1000 * hit->getPosition()[1]; // microns
                    hitsInPlane.measuredZ = 1000 * hit->getPosition()[2]; // microns

                    if( indexconverter[layerIndex] != -1 )
                        _hitsArray[indexconverter[layerIndex]].push_back( hitsInPlane );
                } // end loop over all hits in collection

            }//mode 0

            else if( _inputMode == 2) {

#if defined( USE_ROOT ) || defined(MARLIN_USE_ROOT)

                const float resolX = _testModeSensorResolution;
                const float resolY = _testModeSensorResolution;

                const float xhitpos = gRandom->Uniform(-3500.0,3500.0);
                const float yhitpos = gRandom->Uniform(-3500.0,3500.0);

                const float xslope = gRandom->Gaus(0.0,_testModeXTrackSlope);
                const float yslope = gRandom->Gaus(0.0,_testModeYTrackSlope);

                // loop over all planes
                for( int help = 0; help < _nPlanes; help++ ) {

                    // The x and y positions are given by the sums of the measured
                    // hit positions, the detector resolution, the shifts of the
                    // planes and the effect due to the track slopes.
                    hitsInPlane.measuredX = xhitpos + gRandom->Gaus(0.0,resolX) + _testModeSensorXShifts[help] + _testModeSensorZPositions[help] * tan(xslope) - _testModeSensorGamma[help] * yhitpos - _testModeSensorBeta[help] * _testModeSensorZPositions[0];
                    hitsInPlane.measuredY = yhitpos + gRandom->Gaus(0.0,resolY) + _testModeSensorYShifts[help] + _testModeSensorZPositions[help] * tan(yslope) + _testModeSensorGamma[help] * xhitpos - _testModeSensorAlpha[help] * _testModeSensorZPositions[help];
                    hitsInPlane.measuredZ = _testModeSensorZPositions[help];
                    if( indexconverter[help] != -1 )
                        _hitsArray[indexconverter[help]].push_back( hitsInPlane );

                    // ADDITIONAL PRINTOUTS
                    // printf("plane:%3d hit:%3d %8.3f %8.3f %8.3f \n", help, indexconverter[help], hitsInPlane.measuredX,hitsInPlane.measuredY,hitsInPlane.measuredZ);
                    _hitsArray[help].push_back(hitsInPlane);
                    _telescopeResolX[help] = resolX;
                    _telescopeResolY[help] = resolY;
                } // end loop over all planes

#else // USE_ROOT

                throw MissingLibraryException( this, "ROOT" );

#endif

            } // end if check running in input mode 0 or 2

        }//loop over hits

    //debugging
    streamlog_out( DEBUG0 ) << "hits per plane: " << endl;
    for(size_t j = 0; j < _nPlanes; j++)
    {
        streamlog_out( DEBUG0 ) << "Plane ID = " << j << " : " << "NumbOfHits : " << _hitsArray[j].size() << "  " << endl;
        for( size_t i = 0; i < _hitsArray[j].size(); i++ )
        {
            streamlog_out( DEBUG0 ) << "x ; y ; z : " << _hitsArray[j][i].measuredX << " ; "
                                    << _hitsArray[j][i].measuredY << " ; " << _hitsArray[j][i].measuredZ << endl;
        }
    }

    // mode 1 or 3: tracks are read in
    std::vector<int> fitplane(_nPlanes, 0);

    for( int help = 0; help < _nPlanes; help++ ) {
        fitplane[help] = 1;
    }

    // DP: correlate telescope hits from different planes
    int iA = indexconverter[0]; // plane 0

    if( iA >= 0 ) { // not excluded
        for( size_t jA = 0; jA < _hitsArray[iA].size(); jA++ ) { // hits in plane

            int iB = indexconverter[1]; // plane 1
            if( iB >= 0 ) { // not excluded
                for( size_t jB = 0; jB < _hitsArray[iB].size(); jB++ ) {
                    double dx = _hitsArray[iB][jB].measuredX - _hitsArray[iA][jA].measuredX;
                    double dy = _hitsArray[iB][jB].measuredY - _hitsArray[iA][jA].measuredY;
                    dx01Hist->fill( dx );
                    dy01Hist->fill( dy );
                }//loop hits jB
            }//iB valid

            iB = indexconverter[2]; // plane 2
            if( iB >= 0 ) { // not excluded
                for( size_t jB = 0; jB < _hitsArray[iB].size(); jB++ ) {
                    double dx = _hitsArray[iB][jB].measuredX - _hitsArray[iA][jA].measuredX;
                    double dy = _hitsArray[iB][jB].measuredY - _hitsArray[iA][jA].measuredY;
                    dx02Hist->fill( dx );
                    dy02Hist->fill( dy );
                }//loop hits jB
            }//iB valid

            iB = indexconverter[3]; // plane 3
            if( iB >= 0 ) { // not excluded
                for( size_t jB = 0; jB < _hitsArray[iB].size(); jB++ ) {
                    double dx = _hitsArray[iB][jB].measuredX - _hitsArray[iA][jA].measuredX;
                    double dy = _hitsArray[iB][jB].measuredY - _hitsArray[iA][jA].measuredY;
                    dx03Hist->fill( dx );
                    dy03Hist->fill( dy );
                }//loop hits jB
            }//iB valid

            iB = indexconverter[4]; // plane 4
            if( iB >= 0 ) { // not excluded
                for( size_t jB = 0; jB < _hitsArray[iB].size(); jB++ ) {
                    double dx = _hitsArray[iB][jB].measuredX - _hitsArray[iA][jA].measuredX;
                    double dy = _hitsArray[iB][jB].measuredY - _hitsArray[iA][jA].measuredY;
                    dx04Hist->fill( dx );
                    dy04Hist->fill( dy );
                }//loop hits jB
            }//iB valid

            iB = indexconverter[5]; // plane 5
            if( iB >= 0 ) { // not excluded
                for( size_t jB = 0; jB < _hitsArray[iB].size(); jB++ ) {
                    double dx = _hitsArray[iB][jB].measuredX - _hitsArray[iA][jA].measuredX;
                    double dy = _hitsArray[iB][jB].measuredY - _hitsArray[iA][jA].measuredY;
                    dx05Hist->fill( dx );
                    dy05Hist->fill( dy );
                }//loop hits jB
            }//iB valid

        }//loop hits jA
    }// iA valid

    // triplets 0-1-2:
    // driplets 3-4-5:

    // i is plane index
    // j is hit index
    // k is triplet index

    int i0 = indexconverter[0]; // Tel plane 0
    int i1 = indexconverter[1]; // Tel plane 1
    int i2 = indexconverter[2]; // Tel plane 2

    int i3 = indexconverter[3]; // Tel plane 3
    int i4 = indexconverter[4]; // Tel plane 4
    int i5 = indexconverter[5]; // Tel plane 5

    int i6 = indexconverter[7]; //DUT, Alibava
    int i7 = indexconverter[6]; //REF plane

    // How many tri-dri pairs were matched per event when hits on DUT presented
    int Counter_trk_DUT = 0;
    // How many tri-dri pairs were matched per event when hits on REF presented
    int Counter_trk_REF = 0;
    // How many tri-dri pairs were matched per event when hits on DUT and REF presented
    int Counter_trk_DUT_REF = 0;


    streamlog_out (DEBUG0) << "Triplets matching plane indeces = " << i0 << " " << i1 << " " << i2 << " " << i3 << " " << i4 << " " << i5 <<endl;
    streamlog_out (DEBUG0) << "DUT index = " << i6 << " ; REF index = " << i7 << endl;

    //check the align type
    int align;

    if(_aligntype == 0) {align = i0*i1*i2*i3*i4*i5;} //tel alignment
    if(_aligntype == 1) {align = i0*i1*i2*i3*i4*i5*i7;} //REF alignment
    if(_aligntype == 2) {align = i0*i1*i2*i3*i4*i5*i6;} //DUT alignment

    if( align >= 0 ) { // tel planes not excluded

        int ntri = 0;
        double xmA[_maxTrackCandidates];
        double ymA[_maxTrackCandidates];
        double zmA[_maxTrackCandidates];
        double sxA[_maxTrackCandidates];
        double syA[_maxTrackCandidates];
        int hts[8][_maxTrackCandidates]; // 8 planes

        //triplet:
        for( size_t j0 = 0; j0 < _hitsArray[i0].size(); j0++ ) { // hits on plane 0

            for( size_t j2 = 0; j2 < _hitsArray[i2].size(); j2++ ) { // hits on plane 2

                double dx02 = _hitsArray[i2][j2].measuredX - _hitsArray[i0][j0].measuredX;
                double dy02 = _hitsArray[i2][j2].measuredY - _hitsArray[i0][j0].measuredY;
                double dz02 = _hitsArray[i2][j2].measuredZ - _hitsArray[i0][j0].measuredZ; //difference between z pos

                double tx = dx02 / dz02; // angle theta x
                double ty = dy02 / dz02; // angle theta y

                double avx = 0.5 * ( _hitsArray[i0][j0].measuredX + _hitsArray[i2][j2].measuredX ); // average x why do we need to look at av numb here?
                double avy = 0.5 * ( _hitsArray[i0][j0].measuredY + _hitsArray[i2][j2].measuredY ); // average y
                double avz = 0.5 * ( _hitsArray[i0][j0].measuredZ + _hitsArray[i2][j2].measuredZ ); // average z

                // middle plane for the triplet:

                for( size_t j1 = 0; j1 < _hitsArray[i1].size(); j1++ ) {
                    // triplet residual:
                    double zs = _hitsArray[i1][j1].measuredZ - avz;
                    double xs = avx + tx * zs; // track at 1
                    double ys = avy + ty * zs;

                    double dx = _hitsArray[i1][j1].measuredX - xs;
                    double dy = _hitsArray[i1][j1].measuredY - ys;
                    double dr = sqrt (_hitsArray[i1][j1].measuredX*_hitsArray[i1][j1].measuredX + _hitsArray[i1][j1].measuredY*_hitsArray[i1][j1].measuredY ) -
                            sqrt ( xs*xs + ys*ys );
                    tridx_noCut_Hist->fill(dx);
                    tridy_noCut_Hist->fill(dy);
                    tridr_noCut_Hist->fill(dr);

                    if( abs(dy) < _triCut ) tridxHist->fill( dx );
                    if( abs(dx) < _triCut ) tridyHist->fill( dy );

                    if( abs(dx) < _triCut  && abs(dy) < _triCut ) {
                        //form an estimated hit on the middle plane
                        if( ntri < _maxTrackCandidates ) {
                            xmA[ntri] = avx;
                            ymA[ntri] = avy;
                            zmA[ntri] = avz;
                            sxA[ntri] = tx;
                            syA[ntri] = ty;
                            hts[i0][ntri] = j0;
                            hts[i1][ntri] = j1;
                            hts[i2][ntri] = j2;
                        }
                        ntri++;
                    }//valid triplet

                }//loop hits j1
            }//loop hits j2
        }//loop hits j0

        ntriHist->fill( ntri );

        if( ntri >= _maxTrackCandidates ) {
            streamlog_out( WARNING2 ) << "Maximum number of track candidates reached. Maybe further tracks were skipped" << endl;
        }

        // driplets 3-4-5:

        int ndri = 0;
        double xmB[_maxTrackCandidates];
        double ymB[_maxTrackCandidates];
        double zmB[_maxTrackCandidates];
        double sxB[_maxTrackCandidates];
        double syB[_maxTrackCandidates];

        for( size_t j3 = 0; j3 < _hitsArray[i3].size(); j3++ ) { // hits in plane

            for( size_t j5 = 0; j5 < _hitsArray[i5].size(); j5++ ) {

                double dx35 = _hitsArray[i5][j5].measuredX - _hitsArray[i3][j3].measuredX;
                double dy35 = _hitsArray[i5][j5].measuredY - _hitsArray[i3][j3].measuredY;
                double dz35 = _hitsArray[i5][j5].measuredZ - _hitsArray[i3][j3].measuredZ;

                double tx = dx35 / dz35; // angle theta x
                double ty = dy35 / dz35; // angle theta x

                double avx = 0.5 * ( _hitsArray[i3][j3].measuredX + _hitsArray[i5][j5].measuredX ); // average
                double avy = 0.5 * ( _hitsArray[i3][j3].measuredY + _hitsArray[i5][j5].measuredY ); // average
                double avz = 0.5 * ( _hitsArray[i3][j3].measuredZ + _hitsArray[i5][j5].measuredZ ); // average

                // middle plane for the driplet:

                for( size_t j4 = 0; j4 < _hitsArray[i4].size(); j4++ ) {

                    // triplet residual:

                    double zs = _hitsArray[i4][j4].measuredZ - avz;
                    double xs = avx + tx * zs; // track at 4
                    double ys = avy + ty * zs;

                    double dx = _hitsArray[i4][j4].measuredX - xs;
                    double dy = _hitsArray[i4][j4].measuredY - ys;
                    double dr = sqrt (_hitsArray[i4][j4].measuredX*_hitsArray[i4][j4].measuredX + _hitsArray[i4][j4].measuredY*_hitsArray[i4][j4].measuredY ) -
                            sqrt ( xs*xs + ys*ys );

                    dridx_noCut_Hist->fill(dx);
                    dridy_noCut_Hist->fill(dy);
                    dridr_noCut_Hist->fill(dr);

                    if( abs(dx) < _driCut ) dridxHist->fill( dy );
                    if( abs(dy) < _driCut ) dridyHist->fill( dx );

                    if( abs(dx) < _driCut  && abs(dy) < _driCut ) {

                        if( ndri < _maxTrackCandidates ) {
                            xmB[ndri] = avx;
                            ymB[ndri] = avy;
                            zmB[ndri] = avz;
                            sxB[ndri] = tx;
                            syB[ndri] = ty;
                            hts[i3][ndri] = j3;
                            hts[i4][ndri] = j4;
                            hts[i5][ndri] = j5;
                        }
                        ndri++;

                    }//valid driplet

                }//loop hits j4
            }//loop hits j5
        }//loop hits j3

        ndriHist->fill( ndri );

        if( ndri >= _maxTrackCandidates ) {
            streamlog_out( WARNING2 ) << "Maximum number of track candidates reached. Maybe further tracks were skipped" << endl;
        }
        streamlog_out (DEBUG1) << "numb of triples = " << ntri << " | number of driplets = " << ndri << endl;

        //matching driplets with REF


        // match triplets A to driplets B:
        //how many tracks do we have, how many triplets matched
        _nMilleTracksperEvent = 0;
        //counting the number of trk matching REF/DUT
        ntrk_REF_Event = 0;
        ntrk_DUT_Event = 0;

        // tracks finder
        //trk coordinates
        double xmA_trk[_maxTrackCandidates];
        double ymA_trk[_maxTrackCandidates];
        double zmA_trk[_maxTrackCandidates];
        double sxA_trk[_maxTrackCandidates];
        double syA_trk[_maxTrackCandidates];
        double xmB_trk[_maxTrackCandidates];
        double ymB_trk[_maxTrackCandidates];
        double zmB_trk[_maxTrackCandidates];
        double sxB_trk[_maxTrackCandidates];
        double syB_trk[_maxTrackCandidates];
        double dx_trk[_maxTrackCandidates];
        double dy_trk[_maxTrackCandidates];
        double xA_trk[_maxTrackCandidates];
        double yA_trk[_maxTrackCandidates];
        double xB_trk[_maxTrackCandidates];
        double yB_trk[_maxTrackCandidates];
        int trk[8][_maxTrackCandidates]; //which hit on the plane i# belongs to a track
        //residual arrays for each plane, each track
        double rx[8][_maxTrackCandidates];
        double ry[8][_maxTrackCandidates];
        //hits on how many planes do we have for further analysis
        int nPl = 6;

        for( int kA = 0; (kA < ntri) && (kA < _maxTrackCandidates); ++kA ) { // i = A

            int j2 = hts[i2][kA];
            for( int kB = 0; (kB < ndri) && (kB < _maxTrackCandidates); ++kB ) { // j = B = REF

                int j3 = hts[i3][kB];
                double zmid = 0.5 * ( _hitsArray[i2][j2].measuredZ + _hitsArray[i3][j3].measuredZ );
                //                double zmid = 1000*_siPlaneZPosition[3];  //DUT pos in um

                double xA = xmA[kA] + sxA[kA] * ( zmid - zmA[kA] ); // A at zmid
                double yA = ymA[kA] + syA[kA] * ( zmid - zmA[kA] );

                double xB = xmB[kB] + sxB[kB] * ( zmid - zmB[kB] ); // B at zmid
                double yB = ymB[kB] + syB[kB] * ( zmid - zmB[kB] );

                double dx = xB - xA; // matching residual
                double dy = yB - yA;

                double rA = sqrt(xA*xA + yA*yA);
                double rB = sqrt(xB*xB + yB*yB);
                double dr = rB - rA; //distance difference

                sixdx_noCut_Hist->fill(dx);
                sixdy_noCut_Hist->fill(dy);
                sixdr_noCut_Hist->fill(dr);

                if( abs(dy) < _sixCut ) sixdxHist->fill( dx );
                if( abs(dx) < _sixCut ) sixdyHist->fill( dy );

                if( abs(dx) < _sixCut  && abs(dy) < _sixCut ) { // triplet-driplet match

                    if(_hitsArray[i6].size() > 0) {Counter_trk_DUT++;}
                    if(_hitsArray[i7].size() > 0) {Counter_trk_REF++;}
                    if(_hitsArray[i6].size()*_hitsArray[i7].size() > 0) {Counter_trk_DUT_REF++;}

                    double kx = sxB[kB] - sxA[kA]; //kink x
                    double ky = syB[kB] - syA[kA]; //kink y

                    sixkxHist->fill( kx*1E3 );
                    sixkyHist->fill( ky*1E3 );
                    selxHist->fill( xA*1E-3 ); // triplet at mid
                    selyHist->fill( yA*1E-3 );
                    selaxHist->fill( sxA[kA]*1E3 );//track slope
                    selayHist->fill( syA[kA]*1E3 );

                    //residuals calculation
                    // extrapolate triplet to each plane:
                    //temp res for tel planes
                    double res_x[8];
                    double res_y[8];

                    for( int ipl = 0; ipl < nPl; ipl ++)
                    {
                        int whichhit;
                        if(ipl < 3) { whichhit = hts[ipl][kA];}
                        else if(ipl > 2) { whichhit = hts[ipl][kB];}

                        double dz = _hitsArray[ipl][whichhit].measuredZ - zmA[kA];
                        double xs = xmA[kA] + sxA[kA] * dz; // Ax at plane
                        double ys = ymA[kA] + syA[kA] * dz; // Ay at plane
                        res_x[ipl] = _hitsArray[ipl][whichhit].measuredX - xs; // resid hit-track
                        res_y[ipl] = _hitsArray[ipl][whichhit].measuredY - ys; // resid
                    }

                    seldx1Hist->fill( res_x[i1] ); // triplet interpol
                    seldy1Hist->fill( res_y[i1] );
                    seldx3Hist->fill( res_x[i3] ); // triplet extrapol
                    seldy3Hist->fill( res_y[i3] );
                    seldx4Hist->fill( res_x[i4] );
                    seldy4Hist->fill( res_y[i4] );
                    seldx5Hist->fill( res_x[i5] );
                    seldy5Hist->fill( res_y[i5] );

                    //track sicription:
                    if(_aligntype == 0)
                    {
                        //indeces showing which hit on which plane belongs to the track
                        trk[0][_nMilleTracksperEvent] = hts[i0][kA];
                        trk[1][_nMilleTracksperEvent] = hts[i1][kA];
                        trk[2][_nMilleTracksperEvent] = hts[i2][kA];
                        trk[3][_nMilleTracksperEvent] = hts[i3][kB];
                        trk[4][_nMilleTracksperEvent] = hts[i4][kB];
                        trk[5][_nMilleTracksperEvent] = hts[i5][kB];
                        //writing the arrays of
                        xmA_trk[_nMilleTracksperEvent] = xmA[kA];
                        ymA_trk[_nMilleTracksperEvent] = ymA[kA];
                        zmA_trk[_nMilleTracksperEvent] = zmA[kA];
                        sxA_trk[_nMilleTracksperEvent] = sxA[kA];
                        syA_trk[_nMilleTracksperEvent] = syA[kA];
                        xmB_trk[_nMilleTracksperEvent] = xmB[kB];
                        ymB_trk[_nMilleTracksperEvent] = ymB[kB];
                        zmB_trk[_nMilleTracksperEvent] = zmB[kB];
                        sxB_trk[_nMilleTracksperEvent] = sxB[kB];
                        syB_trk[_nMilleTracksperEvent] = syB[kB];
                        dx_trk[_nMilleTracksperEvent] = dx;
                        dy_trk[_nMilleTracksperEvent] = dy;
                        xA_trk[_nMilleTracksperEvent] = xA;
                        yA_trk[_nMilleTracksperEvent] = yA;
                        xB_trk[_nMilleTracksperEvent] = xB;
                        yB_trk[_nMilleTracksperEvent] = yB;

                        //residuals:
                        for( int ipl = 0; ipl < nPl; ipl ++)
                        {
                            rx[ipl][_nMilleTracksperEvent] = res_x[ipl];
                            ry[ipl][_nMilleTracksperEvent] = res_y[ipl];
                        }

                    }

                    //match the REF hits to the track
                    if(_aligntype == 1)
                    {
                        for(int iREF = 0; iREF < _hitsArray[i7].size(); iREF++)
                        {
                            double dz = _hitsArray[i7][iREF].measuredZ - zmB[kB];
                            double xDrip = xmB[kB] + sxB[kB] * dz; // estimation of driplet to the REF plane
                            double yDrip = ymB[kB] + syB[kB] * dz;

                            double dxREF = _hitsArray[i7][iREF].measuredX - xDrip ;
                            double dyREF = _hitsArray[i7][iREF].measuredY - yDrip ;
                            double drREF = sqrt(_hitsArray[i7][iREF].measuredX*_hitsArray[i7][iREF].measuredX + _hitsArray[i7][iREF].measuredY*_hitsArray[i7][iREF].measuredY) -
                                    sqrt(xDrip*xDrip + yDrip*yDrip);

                            hTrk_REF_dx_noCut->fill(dxREF);
                            hTrk_REF_dy_noCut->fill(dyREF);
                            hTrk_REF_dr_noCut->fill(drREF);

                            if(dxREF < _drirefCut) {hTrk_REF_dy->fill(dyREF);}
                            if(dyREF < _drirefCut) {hTrk_REF_dx->fill(dxREF);}

                            if( (dxREF < _drirefCut) && (dyREF < _drirefCut) )
                            {
                                //indeces showing which hit on which plane belongs to the track
                                trk[0][ntrk_REF_Event] = hts[i0][kA];
                                trk[1][ntrk_REF_Event] = hts[i1][kA];
                                trk[2][ntrk_REF_Event] = hts[i2][kA];
                                trk[3][ntrk_REF_Event] = hts[i3][kB];
                                trk[4][ntrk_REF_Event] = hts[i4][kB];
                                trk[5][ntrk_REF_Event] = hts[i5][kB];
                                trk[6][ntrk_REF_Event] = iREF;
                                //writing the arrays of
                                xmA_trk[ntrk_REF_Event] = xmA[kA];
                                ymA_trk[ntrk_REF_Event] = ymA[kA];
                                zmA_trk[ntrk_REF_Event] = zmA[kA];
                                sxA_trk[ntrk_REF_Event] = sxA[kA];
                                syA_trk[ntrk_REF_Event] = syA[kA];
                                xmB_trk[ntrk_REF_Event] = xmB[kB];
                                ymB_trk[ntrk_REF_Event] = ymB[kB];
                                zmB_trk[ntrk_REF_Event] = zmB[kB];
                                sxB_trk[ntrk_REF_Event] = sxB[kB];
                                syB_trk[ntrk_REF_Event] = syB[kB];
                                dx_trk[ntrk_REF_Event] = dx;
                                dy_trk[ntrk_REF_Event] = dy;
                                xA_trk[ntrk_REF_Event] = xA;
                                yA_trk[ntrk_REF_Event] = yA;
                                xB_trk[ntrk_REF_Event] = xB;
                                yB_trk[ntrk_REF_Event] = yB;

                                nPl = 7;

                                for( int ipl = 0; ipl < nPl; ipl ++)
                                {
                                    if(ipl == 6) {
                                        rx[ipl][ntrk_REF_Event] = dxREF;
                                        ry[ipl][ntrk_REF_Event] = dyREF;}
                                    else if(ipl < 6){
                                        rx[ipl][ntrk_REF_Event] = res_x[ipl];
                                        ry[ipl][ntrk_REF_Event] = res_y[ipl];
                                    }
                                }

                                ntrk_REF_Event++;
                            }
                        }
                    }

                    _nMilleTracksperEvent++;
                    if(_nMilleTracksperEvent > _maxTrackCandidates)
                    {
                        streamlog_out (MESSAGE4) << " Current number of track candidates exceeds the set limit value of " << _maxTrackCandidates << ". Please check the steering file." << endl;
                        break;
                    }

                } //triplet-driplet match
            } //driplet hits
        } //triplet hits

        _nMilleTracksTotal += _nMilleTracksperEvent;
        ntrk_REF_Total += ntrk_REF_Event;

        streamlog_out( DEBUG1 ) << "Triplets matched (tracks) found: " << " event = " << _iEvt << " : " <<  _nMilleTracksperEvent << endl;
        streamlog_out( DEBUG1 ) << "Tracks including REF found: " << " event = " << _iEvt << " : " <<  ntrk_REF_Event << endl;

        nmHist->fill( _nMilleTracksperEvent );
        hTrk_REF->fill(ntrk_REF_Event);

        hTripletsMatched_DUT->fill(Counter_trk_DUT);
        hTripletsMatched_DUT_2D->fill(_iEvt, Counter_trk_DUT);
        hTripletsMatched_REF->fill(Counter_trk_REF);
        hTripletsMatched_REF_2D->fill(_iEvt, Counter_trk_DUT);
        hTripletsMatched_DUT_REF->fill(Counter_trk_DUT_REF);
        hTripletsMatched_DUT_REF_2D->fill(_iEvt, Counter_trk_DUT);

        int tracks;
        if (_aligntype == 0) {tracks = _nMilleTracksperEvent;}
        if (_aligntype == 1) {tracks = ntrk_REF_Event;}


        //forming GBL trajectories from the found above tracks

        for(int iTrk = 0; iTrk < tracks; iTrk++) { // form GBL traj from tracks
            // GBL point vector for the trajectory
            std::vector<gbl::GblPoint> traj_points;

            // build up trajectory:
            std::vector<unsigned int> ilab; // 0-5 = telescope, 6 = DUT, 7 = REF
            vector<double> sPoint;

            double s = 0;

            TMatrixD jacPointToPoint( 5, 5 );

            TMatrixD proL2m(2,2);
            proL2m.UnitMatrix();

            TVectorD meas(2);

            // scatter:
            TVectorD scat(2);
            scat.Zero(); //mean is zero

            double p = _eBeam; // beam momentum
            double X0Si = 65e-3 / 94; // Si + Kapton
            double tetSi = 0.0136 * sqrt(X0Si) / p * ( 1 + 0.038*std::log(X0Si) );

            TVectorD wscatSi(2);
            wscatSi[0] = 1.0 / ( tetSi * tetSi ); //weight
            wscatSi[1] = 1.0 / ( tetSi * tetSi );

            TMatrixD alDer2( 2, 2 ); // alignment derivatives
            alDer2[0][0] = 1.0; // dx/dx GBL sign convetion
            alDer2[1][0] = 0.0; // dy/dx
            alDer2[0][1] = 0.0; // dx/dy
            alDer2[1][1] = 1.0; // dy/dy

            TMatrixD alDer3( 2, 3 ); // alignment derivatives
            alDer3[0][0] = 1.0; // dx/dx
            alDer3[1][0] = 0.0; // dy/dx
            alDer3[0][1] = 0.0; // dx/dy
            alDer3[1][1] = 1.0; // dy/dy

            TMatrixD alDer4( 2, 4 ); // alignment derivatives
            alDer4[0][0] = 1.0; // dx/dx
            alDer4[1][0] = 0.0; // dy/dx
            alDer4[0][1] = 0.0; // dx/dy
            alDer4[1][1] = 1.0; // dy/dy
            alDer4[0][3] = sxA_trk[iTrk]; // dx/dz
            alDer4[1][3] = syA_trk[iTrk]; // dx/dz

            // telescope planes 0-5:

            double zprev = _hitsArray[0][iTrk].measuredZ; // id0 tel. plane, including any pre-alignment

            //precision for adding a track to trajectory
            TVectorD measPrec(2); // precision = 1/resolution^2
            //plane resolutions
            double resx;
            double resy;

            for( int ipl = 0; ipl < nPl; ipl++ ) {
                //tel alignment
                if(_aligntype == 0){

                    resx = _telescopeResolution; // [um] telescope initial resolution
                    resy = _telescopeResolution; // [um] telescope initial resolution
                }
                //REF alignment
                if(_aligntype == 1){

                    if(ipl < 6){
                        resx = _telescopeResolution; // [um] telescope initial resolution
                        resy = _telescopeResolution; // [um] telescope initial resolution
                    }
                    if(ipl == 6){
                        resx = _REFResolution; // [um] telescope initial resolution
                        resy = _REFResolution; // [um] telescope initial resolution
                    }

                }

                measPrec[0] = 1.0 / resx / resx;
                measPrec[1] = 1.0 / resy / resy;

                double zz = _hitsArray[ipl][iTrk].measuredZ; // [um]
                double step;

                step = zz - zprev;

                jacPointToPoint = Jac55( step );
                gbl::GblPoint *point = new gbl::GblPoint( jacPointToPoint );
                s += step;
                zprev = zz;

                meas[0] =  rx[ipl][iTrk];
                meas[1] =  ry[ipl][iTrk];

                streamlog_out (DEBUG0) << "GBL points formation :"<< endl
                                          << "plane = " << ipl << " | geomID = " << _siPlanesLayerLayout->getID(ipl) << " : X/Y/Z = "
                                             << _hitsArray[ipl][iTrk].measuredX <<"/" << _hitsArray[ipl][iTrk].measuredY <<"/" <<
                                                _hitsArray[ipl][iTrk].measuredZ << endl
                                                << "meas[0] = " << meas[0] << " ;  meas[1] = " << meas[1] << " | step = " << step << " s = " << s << endl;

                point->addMeasurement( proL2m, meas, measPrec );

                point->addScatterer( scat, wscatSi );

                if( _alignMode == 2 ) { // only x and y shifts
                    // global labels for MP:
                    std::vector<int> globalLabels(2);
                    globalLabels[0] = 1 + 2*ipl;
                    globalLabels[1] = 2 + 2*ipl;
                    point->addGlobals( globalLabels, alDer2 ); // for MillePede alignment
                }
                else if( _alignMode == 3 ) { // with rot
                    std::vector<int> globalLabels(3);
                    globalLabels[0] = 1 + 3*ipl; // x
                    globalLabels[1] = 2 + 3*ipl; // y
                    globalLabels[2] = 3 + 3*ipl; // rot
                    //alDer3[0][2] = -_hitsArray[ipl][jhit].measuredY; // dx/dphi
                    //alDer3[1][2] =  _hitsArray[ipl][jhit].measuredX; // dy/dphi
                    double dz = _hitsArray[ipl][trk[ipl][iTrk]].measuredZ - zmA_trk[iTrk];
                    double xs = xmA_trk[iTrk] + sxA_trk[iTrk] * dz; // Ax at plane
                    double ys = ymA_trk[iTrk] + syA_trk[iTrk] * dz; // Ay at plane
                    alDer3[0][2] = -ys; // dx/dphi
                    alDer3[1][2] =  xs; // dy/dphi
                    point->addGlobals( globalLabels, alDer3 ); // for MillePede alignment
                }
                else if( _alignMode == 4 ) { // with rot and z shift
                    std::vector<int> globalLabels(4);
                    globalLabels[0] = 1 + 4*ipl;
                    globalLabels[1] = 2 + 4*ipl;
                    globalLabels[2] = 3 + 4*ipl;
                    globalLabels[3] = 4 + 4*ipl; // z
                    //alDer4[0][2] = -_hitsArray[ipl][jhit].measuredY; // dx/dphi
                    //alDer4[1][2] =  _hitsArray[ipl][jhit].measuredX; // dy/dphi
                    double dz = _hitsArray[ipl][trk[ipl][iTrk]].measuredZ - zmA_trk[iTrk];
                    double xs = xmA_trk[iTrk] + sxA_trk[iTrk] * dz; // Ax at plane
                    double ys = ymA_trk[iTrk] + syA_trk[iTrk] * dz; // Ay at plane
                    alDer4[0][2] = -ys; // dx/dphi
                    alDer4[1][2] =  xs; // dy/dphi
                    point->addGlobals( globalLabels, alDer4 ); // for MillePede alignment
                }

                //unsigned int iLabel = traj.addPoint(*point);
                traj_points.push_back(*point);

                //ilab[ipl] = iLabel;
                sPoint.push_back( s );

                delete point;

                if( (DUT_scatterer == 1) && (ipl == 2) ) { // insert DUT

                    zz =  _siPlaneZPosition[i6]*1000; //DUT pos in um
                    step = zz - _hitsArray[i2][0].measuredZ;


                    jacPointToPoint = Jac55( step );
                    point = new gbl::GblPoint( jacPointToPoint );
                    s += step;
                    sPoint.push_back( s );
                    zprev = zz;

                    streamlog_out (DEBUG0) << "Add DUT as scatterer :"<< endl
                                             << " | geomID = " << _siPlanesLayerLayout->getID(i6) << " : zz = " << zz << " step = " <<step << " s = " << s << endl;

                    double tetDUT = 0.0136 * sqrt(DUTX0) / p * ( 1 + 0.038*std::log(DUTX0) );

                    TVectorD wscatDUT(2);
                    wscatDUT[0] = 1.0 / ( tetDUT * tetDUT ); //weight
                    wscatDUT[1] = 1.0 / ( tetDUT * tetDUT );

                    point->addScatterer( scat, wscatDUT );

                    //iLabel = traj.addPoint(*point);
                    traj_points.push_back(*point);

                    // ilab[6] = iLabel;

                    delete point;

                }//DUT present

            } // loop over planes


            double Chi2;
            int Ndf;
            double lostWeight;

            gbl::GblTrajectory traj(traj_points, false); // curvature = false
            traj.fit( Chi2, Ndf, lostWeight );
            traj.getLabels(ilab);

            // debug:
            //                        streamlog_out (DEBUG0) << "check 3: the traj points:" << endl;
            //                        for(int i =0; i < 7; i++)
            //                        {
            //                            streamlog_out (DEBUG0) << ilab[i] << endl;
            //                        }

            /*
                                                 *
                                       cout << "traj with " << traj.getNumPoints() << " points:" << endl;
                                       for( int ipl = 0; ipl < 6; ++ipl ){
                                         cout << "  plane " << ipl << ", lab " << ilab[ipl];
                                         cout << "  z " << sPoint[ipl]*1E-3;
                                         cout << " dx " << rx[ipl];
                                         cout << " dy " << ry[ipl];
                                         cout << endl;
                                       }
                                       */

            //cout << " chi2 " << Chi2 << ", ndf " << Ndf << endl;

            gblndfHist->fill( Ndf );
            gblchi2Hist->fill( Chi2 );
            double probchi = TMath::Prob( Chi2, Ndf );
            gblprbHist->fill( probchi );

            // bad fits:

            if( probchi < _chi2tondfcut ) {

                badxHist->fill( xA_trk[iTrk]*1E-3 ); // triplet at DUT
                badyHist->fill( yA_trk[iTrk]*1E-3 );
                badaxHist->fill( sxA_trk[iTrk]*1E3 );
                badayHist->fill( syA_trk[iTrk]*1E3 );
                baddxHist->fill( dx_trk[iTrk] ); // triplet-driplet match
                baddyHist->fill( dy_trk[iTrk] );
                badkxHist->fill( (sxB_trk[iTrk] - sxA_trk[iTrk])*1E3 ); // triplet-driplet kink
                badkyHist->fill( (syB_trk[iTrk] - syA_trk[iTrk])*1E3 );

                baddx1Hist->fill( rx[i1][iTrk] ); // triplet interpol
                baddy1Hist->fill( ry[i1][iTrk] );
                baddx3Hist->fill( rx[i3][iTrk] ); // triplet extrapol
                baddy3Hist->fill( ry[i3][iTrk] );
                baddx4Hist->fill( rx[i4][iTrk] );
                baddy4Hist->fill( ry[i4][iTrk] );
                baddx5Hist->fill( rx[i5][iTrk] );
                baddy5Hist->fill( ry[i5][iTrk] );

            }// bad fit

            else {

                goodx1Hist->fill( rx[i1][iTrk] ); // triplet interpol
                goody1Hist->fill( ry[i1][iTrk] );

            } // OK fit

            // look at fit:

            TVectorD aCorrection(5);
            TMatrixDSym aCovariance(5);

            double ax[8];
            double ay[8];
            unsigned int k = 0;

            // at plane 0:

            int ipos = ilab[0];
            traj.getResults( ipos, aCorrection, aCovariance );

            //track = q/p, x', y', x, y
            //        0,   1,  2,  3, 4

            gblax0Hist->fill( aCorrection[1]*1E3 ); // angle x [mrad]
            gbldx0Hist->fill( aCorrection[3] ); // shift x [um]
            gblrx0Hist->fill( rx[i0][iTrk] - aCorrection[3] ); // residual x [um]
            ax[k] = aCorrection[1]; // angle correction at plane, for kinks
            ay[k] = aCorrection[2]; // angle correction at plane, for kinks
            k++;

            ipos = ilab[1];
            traj.getResults( ipos, aCorrection, aCovariance );
            gblax1Hist->fill( aCorrection[1]*1E3 ); // angle x [mrad]
            gbldx1Hist->fill( aCorrection[3] ); // shift x [um]
            gblrx1Hist->fill( rx[i1][iTrk] - aCorrection[3] ); // residual x [um]
            ax[k] = aCorrection[1]; // angle correction at plane, for kinks
            ay[k] = aCorrection[2]; // angle correction at plane, for kinks
            k++;

            ipos = ilab[2];
            traj.getResults( ipos, aCorrection, aCovariance );
            gblax2Hist->fill( aCorrection[1]*1E3 ); // angle x [mrad]
            gbldx2Hist->fill( aCorrection[3] ); // shift x [um]
            gblrx2Hist->fill( rx[i2][iTrk] - aCorrection[3] ); // residual x [um]
            ax[k] = aCorrection[1]; // angle correction at plane, for kinks
            ay[k] = aCorrection[2]; // angle correction at plane, for kinks
            k++;

            //                        if( DUT_scatterer == 1) {
            //                            ipos = ilab[6]; // 6 = DUT
            //                            traj.getResults( ipos, aCorrection, aCovariance );
            //                            gblax6Hist->fill( aCorrection[1]*1E3 ); // angle x [mrad]
            //                            gbldx6Hist->fill( aCorrection[3] ); // shift x [um]
            //                            gbldy6Hist->fill( aCorrection[4] ); // shift x [um]
            //                            ax[k] = aCorrection[1]; // angle correction at plane, for kinks
            //                            ay[k] = aCorrection[2]; // angle correction at plane, for kinks
            //                            k++;
            //                        }

            ipos = ilab[3];
            traj.getResults( ipos, aCorrection, aCovariance );
            gblax3Hist->fill( aCorrection[1]*1E3 ); // angle x [mrad]
            gbldx3Hist->fill( aCorrection[3] ); // shift x [um]
            gblrx3Hist->fill( rx[i3][iTrk] - aCorrection[3] ); // residual x [um]
            ax[k] = aCorrection[1]; // angle correction at plane, for kinks
            ay[k] = aCorrection[2]; // angle correction at plane, for kinks
            k++;

            ipos = ilab[4];
            traj.getResults( ipos, aCorrection, aCovariance );
            gblax4Hist->fill( aCorrection[1]*1E3 ); // angle x [mrad]
            gbldx4Hist->fill( aCorrection[3] ); // shift x [um]
            gblrx4Hist->fill( rx[i4][iTrk] - aCorrection[3] ); // residual x [um]
            ax[k] = aCorrection[1]; // angle correction at plane, for kinks
            ay[k] = aCorrection[2]; // angle correction at plane, for kinks
            k++;

            ipos = ilab[5];
            traj.getResults( ipos, aCorrection, aCovariance );
            gblax5Hist->fill( aCorrection[1]*1E3 ); // angle x [mrad]
            gbldx5Hist->fill( aCorrection[3] ); // shift x [um]
            gblrx5Hist->fill( rx[i5][iTrk] - aCorrection[3] ); // residual x [um]
            ax[k] = aCorrection[1]; // angle correction at plane, for kinks
            ay[k] = aCorrection[2]; // angle correction at plane, for kinks
            k++;

            // kinks: 1,2 = tele, 3,4 = tele,

            gblkx1Hist->fill( (ax[1] - ax[0])*1E3 ); // kink at 1 [mrad]
            gblkx2Hist->fill( (ax[2] - ax[1])*1E3 ); // kink at 2 [mrad]
            gblkx3Hist->fill( (ax[3] - ax[2])*1E3 ); // kink at 3 [mrad]
            gblkx4Hist->fill( (ax[4] - ax[3])*1E3 ); // kink at 4 [mrad]
            gblkx5Hist->fill( (ax[5] - ax[4])*1E3 ); // kink at 5 [mrad]
            //                        if( DUT_scatterer == 1 ) gblkx6Hist->fill( (ax[6] - ax[5])*1E3 ); // kink at 6 [mrad]

            traj.milleOut( *milleGBL );
        }
    }


    streamlog_out (DEBUG0) << "--------- new event ----" << endl;
    // count events:

    ++_iEvt;
    if( isFirstEvent() ) _isFirstEvent = false;

}

//------------------------------------------------------------------------------
void EUTelMilleGBL::end() {

    delete [] _telescopeResolY;
    delete [] _telescopeResolX;
    delete [] _telescopeResolZ;
    delete [] _yFitPos;
    delete [] _xFitPos;
    delete [] _waferResidY;
    delete [] _waferResidX;
    delete [] _waferResidZ;

    // close the output file:
    delete milleGBL;

    // if write the pede steering file
    if( _generatePedeSteerfile ) {

        streamlog_out( MESSAGE4 ) << endl << "Generating the steering file for the pede program..." << endl;

        ofstream steerFile;
        steerFile.open(_pedeSteerfileName.c_str());

        if( steerFile.is_open()) {

            // find first and last excluded plane
            int firstnotexcl = _nPlanes;
            int lastnotexcl = 0;

            // loop over all planes:

            for( int ipl = 0; ipl < _nPlanes; ipl++) {

                int excluded = 0;

                // loop over excluded planes:

                for( int jpl = 0; jpl < _nExcludePlanes; jpl++ ) {
                    if( ipl == _excludePlanes[jpl] ) excluded = 1;
                }

                if( excluded == 0 && firstnotexcl > ipl ) firstnotexcl = ipl;

                if( excluded == 0 && lastnotexcl < ipl ) lastnotexcl = ipl;

            } // end loop over all planes

            steerFile << "Cfiles" << endl;
            steerFile << _binaryFilename << endl;
            steerFile << endl;

            steerFile << "Parameter" << endl;

            int counter = 0;
            int nfix = 0;

            // loop over all planes:

            for( int ipl = 0; ipl < _nPlanes; ipl++) {

                //                    streamlog_out(DEBUG0) << "mille writing parameters :" <<endl
                //                                          << "Counter = " << ipl << " ; gearPlaneID = " << _orderedSensorID_wo_excluded.at(ipl) << endl;

                int excluded = 0; // flag for excluded planes

                //                    streamlog_out (DEBUG0) << "Plane " << ipl << " exclude = " << excluded << endl;

                if( excluded == 0 ) {

                    bool fixed = false;
                    for( size_t i = 0;i < _FixedPlanes.size(); i++ ) {
                        if( _FixedPlanes[i] == ipl )
                        {
                            fixed = true;
                            break;
                        }
                    }

                    // if fixed planes
                    //if non of the planes are fixed, the first id0 and the last planes are being fixed
                    if( fixed || (_FixedPlanes.size() == 0 && (ipl == firstnotexcl || ipl == lastnotexcl) ) ) {
                        nfix++;
                        if( _alignMode == 2 ) {
                            steerFile << (counter * 2 + 1) << "  0.0 -1.0" << endl;
                            steerFile << (counter * 2 + 2) << "  0.0 -1.0" << endl;
                        }
                        if( _alignMode == 3 ) {
                            steerFile << (counter * 3 + 1) << "  0.0 -1.0" << endl; // fix x
                            steerFile << (counter * 3 + 2) << "  0.0 -1.0" << endl; // fix y
                            steerFile << (counter * 3 + 3) << "  0.0 -1.0" << endl; // fix rot
                        }
                        if( _alignMode == 4 ) {
                            steerFile << (counter * 4 + 1) << "  0.0 -1.0" << endl;
                            steerFile << (counter * 4 + 2) << "  0.0 -1.0" << endl;
                            steerFile << (counter * 4 + 3) << "  0.0 -1.0" << endl;
                        }
                    }

                    else {

                        if( _alignMode == 2 ) {
                            steerFile << (counter * 2 + 1) << "  0.0  0.0" << endl;
                            steerFile << (counter * 2 + 2) << "  0.0  0.0" << endl;
                        }

                        if( _alignMode == 3 ) {
                            steerFile << (counter * 3 + 1) << "  0.0  0.0" << endl;
                            steerFile << (counter * 3 + 2) << "  0.0  0.0" << endl;
                            steerFile << (counter * 3 + 3) << "  0.0  0.0" << endl;
                        }

                        if( _alignMode == 4 ) {
                            steerFile << (counter * 4 + 1) << "  0.0  0.0" << endl;
                            steerFile << (counter * 4 + 2) << "  0.0  0.0" << endl;
                            steerFile << (counter * 4 + 3) << "  0.0  0.0" << endl;
                        }

                    }// not fixed

                    // special for z shift:

                    if( _alignMode == 4 ) {
                        if( ipl == 1 )
                            steerFile << (counter * 4 + 4) << "  0.0 -1.0" << endl;
                        else if( ipl == 4 )
                            steerFile << (counter * 4 + 4) << "  0.0 -1.0" << endl;
                        else
                            steerFile << (counter * 4 + 4) << "  0.0  0.0" << endl;
                    }

                    ++counter;

                } // end if plane not excluded

            } // end loop over all planes

            if( nfix < 2 ) {

                if( _alignMode == 2 ) {

                    steerFile << "Constraint 0 ! sum dx = 0" << endl;
                    steerFile << " 1  1.0" << endl;
                    steerFile << " 3  1.0" << endl;
                    steerFile << " 5  1.0" << endl;
                    steerFile << " 7  1.0" << endl;
                    steerFile << " 9  1.0" << endl;
                    steerFile << "11  1.0" << endl;

                    steerFile << "Constraint 0 ! sum dy = 0" << endl;
                    steerFile << " 2  1.0" << endl;
                    steerFile << " 4  1.0" << endl;
                    steerFile << " 6  1.0" << endl;
                    steerFile << " 8  1.0" << endl;
                    steerFile << "10  1.0" << endl;
                    steerFile << "12  1.0" << endl;
                }

                if( _alignMode == 3 ) {

                    steerFile << "Constraint 0 ! sum dx = 0" << endl;
                    steerFile << " 1  1.0" << endl;
                    steerFile << " 4  1.0" << endl;
                    steerFile << " 7  1.0" << endl;
                    steerFile << "10  1.0" << endl;
                    steerFile << "13  1.0" << endl;
                    steerFile << "16  1.0" << endl;

                    steerFile << "Constraint 0 ! sum dy = 0" << endl;
                    steerFile << " 2  1.0" << endl;
                    steerFile << " 5  1.0" << endl;
                    steerFile << " 8  1.0" << endl;
                    steerFile << "11  1.0" << endl;
                    steerFile << "14  1.0" << endl;
                    steerFile << "17  1.0" << endl;

                    steerFile << "Constraint 0 ! sum dphi = 0" << endl;
                    steerFile << " 3  1.0" << endl;
                    steerFile << " 6  1.0" << endl;
                    steerFile << " 9  1.0" << endl;
                    steerFile << "12  1.0" << endl;
                    steerFile << "15  1.0" << endl;
                    steerFile << "18  1.0" << endl;
                }

                if( _alignMode == 4 ) {

                    steerFile << "Constraint 0 ! sum dx = 0" << endl;
                    steerFile << " 1  1.0" << endl;
                    steerFile << " 5  1.0" << endl;
                    steerFile << " 9  1.0" << endl;
                    steerFile << "13  1.0" << endl;
                    steerFile << "17  1.0" << endl;
                    steerFile << "21  1.0" << endl;

                    steerFile << "Constraint 0 ! sum dy = 0" << endl;
                    steerFile << " 2  1.0" << endl;
                    steerFile << " 6  1.0" << endl;
                    steerFile << "10  1.0" << endl;
                    steerFile << "14  1.0" << endl;
                    steerFile << "18  1.0" << endl;
                    steerFile << "22  1.0" << endl;

                    steerFile << "Constraint 0 ! sum dphi = 0" << endl;
                    steerFile << " 3  1.0" << endl;
                    steerFile << " 7  1.0" << endl;
                    steerFile << "11  1.0" << endl;
                    steerFile << "15  1.0" << endl;
                    steerFile << "19  1.0" << endl;
                    steerFile << "23  1.0" << endl;

                    steerFile << "Constraint 0 ! sum dz = 0" << endl;
                    steerFile << " 4  1.0" << endl;
                    steerFile << " 8  1.0" << endl;
                    steerFile << "12  1.0" << endl;
                    steerFile << "16  1.0" << endl;
                    steerFile << "20  1.0" << endl;
                    steerFile << "24  1.0" << endl;
                }

            }//nfix < 2

            steerFile << endl;
            steerFile << "! chiscut 5.0 2.5" << endl;
            steerFile << "outlierdownweighting 4" << endl;
            steerFile << "dwfractioncut 0.1" << endl;
            steerFile << endl;
            steerFile << "method inversion 10  0.1" << endl;
            // Use 10 OpenMP threads to process the data:
            steerFile << "threads 10 1" << endl;
            steerFile << endl;
            steerFile << "histprint" << endl;
            steerFile << endl;
            steerFile << "end" << endl;

            steerFile.close();

            streamlog_out( MESSAGE2 ) << "File " << _pedeSteerfileName << " written." << endl;

        }
        else {
            streamlog_out( ERROR2 ) << "Could not open steering file." << endl;
        }

    } // end if write the pede steering file

    streamlog_out( MESSAGE2 ) << endl;
    if(_aligntype == 0){  streamlog_out( MESSAGE4 ) << "Number of tracks for millipede: " << _nMilleTracksTotal << endl; }
    if(_aligntype == 1){  streamlog_out( MESSAGE4 ) << "Number of tracks for millipede: " << ntrk_REF_Total << endl; }

    // if running pede using the generated steering file

    if( _runPede == 1 ) {

        // check if steering file exists

        if( _generatePedeSteerfile == 1 ) {

            std::string command = "pede " + _pedeSteerfileName;

            // before starting pede, let's check if it is in the path
            bool isPedeInPath = true;

            // create a new process
            redi::ipstream which("which pede");

            // wait for the process to finish
            which.close();

            // get the status
            // if it 255 then the program wasn't found in the path
            isPedeInPath = !( which.rdbuf()->status() == 255 );

            if( !isPedeInPath ) {
                streamlog_out( ERROR ) << "Pede cannot be executed because not found in the path" << endl;
            }
            else {

                streamlog_out( MESSAGE2 ) << endl;
                streamlog_out( MESSAGE2 ) << "Starting pede..." << endl;
                streamlog_out( MESSAGE2 ) << command.c_str() << endl;

                redi::ipstream pede( command.c_str() );
                string output;
                while ( getline( pede, output ) ) {
                    streamlog_out( MESSAGE2 ) << output << endl;
                }

                // wait for the pede execution to finish
                pede.close();

                // check the exit value of pede
                if( pede.rdbuf()->status() == 0 ) {
                    streamlog_out( MESSAGE2 ) << "Pede successfully finished" << endl;
                }

                // reading back the millepede.res file:

                string millepedeResFileName = "millepede.res";

                streamlog_out( MESSAGE2 ) << "Reading back the " << millepedeResFileName << endl
                                          << "Saving the alignment constant into " << _alignmentConstantLCIOFile << endl;

                // open the millepede ASCII output file
                ifstream millepede( millepedeResFileName.c_str() );

                // reopen the LCIO file this time in append mode
                LCWriter * lcWriter = LCFactory::getInstance()->createLCWriter();

                try {
                    lcWriter->open( _alignmentConstantLCIOFile, LCIO::WRITE_NEW );
                }
                catch ( IOException& e ) {
                    streamlog_out( ERROR4 ) << e.what() << endl
                                            << "Sorry for quitting. " << endl;
                    exit(-1);
                }

                // write an almost empty run header
                LCRunHeaderImpl * lcHeader  = new LCRunHeaderImpl;
                lcHeader->setRunNumber( 0 );

                lcWriter->writeRunHeader(lcHeader);

                delete lcHeader;

                LCEventImpl * event = new LCEventImpl;
                event->setRunNumber( 0 );
                event->setEventNumber( 0 );

                LCTime * now = new LCTime;
                event->setTimeStamp( now->timeStamp() );
                delete now;

                LCCollectionVec * constantsCollection = new LCCollectionVec( LCIO::LCGENERICOBJECT );

                if( millepede.bad() || !millepede.is_open() ) {
                    streamlog_out( ERROR4 ) << "Error opening the " << millepedeResFileName << endl
                                            << "The alignment slcio file cannot be saved" << endl;
                }
                else {
                    vector<double > tokens;
                    stringstream tokenizer;
                    string line;
                    double buffer;

                    // get the first line and throw it away since it is a comment!

                    getline( millepede, line );
                    std::cout << "line:" <<  line  << std::endl;

                    int counter = 0;

                    while( ! millepede.eof() ) {

                        EUTelAlignmentConstant * constant = new EUTelAlignmentConstant;

                        bool goodLine = true;

                        unsigned int numpars = 2; // align pars per plane in Pede
                        if(      _alignMode == 3 )
                            numpars = 3;
                        else if( _alignMode == 4 )
                            numpars = 4;

                        for( unsigned int iParam = 0; iParam < numpars; ++iParam ) {

                            getline( millepede, line );

                            if( line.empty() ) {
                                goodLine = false;
                                continue;
                            }

                            tokens.clear();
                            tokenizer.clear();
                            tokenizer.str( line );

                            while( tokenizer >> buffer ) {
                                tokens.push_back( buffer );
                            }

                            if( ( tokens.size() == 3 ) || ( tokens.size() == 6 ) || (tokens.size() == 5) ) {
                                goodLine = true;
                            }
                            else goodLine = false;

                            bool isFixed = ( tokens.size() == 3 );
                            if( isFixed ) {
                                streamlog_out( DEBUG1 )
                                        << "Parameter " << tokens[0]
                                        << " is at " << ( tokens[1] / 1000 )
                                        << " (fixed)"  << endl;
                            }
                            else {
                                streamlog_out( DEBUG1 )
                                        << "Parameter " << tokens[0]
                                        << " is at " << (tokens[1] / 1000 )
                                        << " +/- " << ( tokens[4] / 1000 )  << endl;
                            }

                            if( iParam == 0 ) {
                                constant->setXOffset( tokens[1] / 1000 );
                                if( ! isFixed ) constant->setXOffsetError( tokens[4] / 1000 );
                            }
                            if( iParam == 1 ) {
                                constant->setYOffset( tokens[1] / 1000 );
                                if( ! isFixed ) constant->setYOffsetError( tokens[4] / 1000 );
                            }
                            if( iParam == 2 ) {
                                constant->setGamma( tokens[1]  );
                                if( ! isFixed ) constant->setGammaError( tokens[4] );
                            }
                            if( iParam == 3 ) {
                                constant->setZOffset( tokens[1] / 1000 );
                                if( ! isFixed ) constant->setZOffsetError( tokens[4] / 1000 );
                            }

                        }//loop param

                        // right place to add the constant to the collection:

                        if( goodLine ) {
                            constant->setSensorID( _orderedSensorID_wo_excluded.at( counter ) );
                            streamlog_out (DEBUG0) << "check the IDs after mille :" << endl
                                                   << "counter = " << counter << " ; PlaneID = " << _orderedSensorID_wo_excluded.at( counter ) <<endl;
                            ++counter;
                            constantsCollection->push_back( constant );
                            streamlog_out( MESSAGE4 ) << (*constant) << endl;

                        }
                        else delete constant;

                    }//millepede eof

                }//millepede OK

                event->addCollection( constantsCollection, _alignmentConstantCollectionName );
                lcWriter->writeEvent( event );
                delete event;

                lcWriter->close();

                millepede.close();

            }//PedeInPath

        }// Pede steering exist

        else {
            streamlog_out( ERROR2 ) << "Unable to run pede. No steering file has been generated." << endl;
        }

    } // Pede wanted

    streamlog_out( MESSAGE2 ) << endl;
    streamlog_out( MESSAGE2 ) << "Successfully finished" << endl;

}//end

//------------------------------------------------------------------------------
void EUTelMilleGBL::bookHistos() {

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

    try {
        streamlog_out( MESSAGE2 ) << "Booking histograms..." << endl;

        int lim = 1000; //the lim (from left and right) of res in um
        int bin = 10 ; //bin size in um
        int Nbin = 2*lim/bin+1;

        //DP:

        dx01Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "dx01", Nbin, -lim, lim );
        dx01Hist->setTitle( "dx01;x_{1}-x_{0} [#mum];hit pairs" );

        dy01Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "dy01", Nbin, -lim, lim );
        dy01Hist->setTitle( "dy01;y_{1}-y_{0} [#mum];hit pairs" );

        dx02Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "dx02", Nbin, -lim, lim );
        dx02Hist->setTitle( "dx02;x_{2}-x_{0} [#mum];hit pairs" );

        dy02Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "dy02", Nbin, -lim, lim );
        dy02Hist->setTitle( "dy02;y_{2}-y_{0} [#mum];hit pairs" );

        dx03Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "dx03", Nbin, -lim, lim );
        dx03Hist->setTitle( "dx03;x_{3}-x_{0} [#mum];hit pairs" );

        dy03Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "dy03", Nbin, -lim, lim );
        dy03Hist->setTitle( "dy03;y_{3}-y_{0} [#mum];hit pairs" );

        dx04Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "dx04", Nbin, -lim, lim );
        dx04Hist->setTitle( "dx04;x_{4}-x_{0} [#mum];hit pairs" );

        dy04Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "dy04", Nbin, -lim, lim );
        dy04Hist->setTitle( "dy04;y_{4}-y_{0} [#mum];hit pairs" );

        dx05Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "dx05", Nbin, -lim, lim );
        dx05Hist->setTitle( "dx05;x_{5}-x_{0} [#mum];hit pairs" );

        dy05Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "dy05", Nbin, -lim, lim );
        dy05Hist->setTitle( "dy05;y_{5}-y_{0} [#mum];hit pairs" );

        //triplet-driplet histos

        //triplet
        tridxHist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "tridx", 4*_triCut+1, -2*_triCut, 2*_triCut );
        tridxHist->setTitle( "tridx;x_{1}-x_{02} [#mum];triplet" );

        tridyHist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "tridy", 4*_triCut+1, -2*_triCut, 2*_triCut );
        tridyHist->setTitle( "tridy;y_{1}-y_{02} [#mum];triplet" );

        tridx_noCut_Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "tridx_noCut", 10*_triCut+1, -5*_triCut, 5*_triCut );
        tridx_noCut_Hist->setTitle( "tridx no triplet cut applied;x_{1}-x_{02} [#mum];triplet" );

        tridy_noCut_Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "tridy_noCut", 10*_triCut+1, -5*_triCut, 5*_triCut );
        tridy_noCut_Hist->setTitle( "tridy no triplet cut applied;y_{1}-y_{02} [#mum];triplet" );

        tridr_noCut_Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "tridr_noCut", 10*_triCut+1, -5*_triCut, 5*_triCut );
        tridr_noCut_Hist->setTitle( "tridr, no driplet cut applied; r_{1}-r_{02} [#mum]; triplet" );

        ntriHist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "ntri", 21, -0.5, 20.5 );
        ntriHist->setTitle( "Number of triplets distr.; Number of triplets per event; Numb. of Events" );

        //driplet
        dridxHist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "dridx", 4*_driCut+1, -2*_driCut, 2*_driCut );
        dridxHist->setTitle( "dridx;x_{4}-x_{35} [#mum];triplet" );

        dridyHist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "dridy", 4*_driCut+1, -2*_driCut, 2*_driCut );
        dridyHist->setTitle( "dridy;y_{4}-y_{35} [#mum];triplet" );

        dridx_noCut_Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "dridx_noCut", 10*_driCut+1, -5*_driCut, 5*_driCut );
        dridx_noCut_Hist->setTitle( "dridx, no driplet cut applied;x_{4}-x_{35} [#mum];triplet" );

        dridy_noCut_Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "dridy_noCut", 10*_driCut+1, -5*_driCut, 5*_driCut );
        dridy_noCut_Hist->setTitle( "dridy, no driplet cut applied;y_{4}-y_{35} [#mum];triplet" );

        dridr_noCut_Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "dridr_noCut", 10*_driCut+1, -5*_driCut, 5*_driCut );
        dridr_noCut_Hist->setTitle( "dridr, no driplet cut applied; r_{4}-r_{35} [#mum]; triplet" );

        ndriHist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "ndri", 21, -0.5, 20.5 );
        ndriHist->setTitle( "Number of driplets distr.; Number of driplets per event; Numb. of Events" );

        //tracks (triplet-driplet matching)
        sixdxHist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "sixdx", 4*_sixCut+1, -2*_sixCut, 2*_sixCut );
        sixdxHist->setTitle( "Triplet-driplet X displacement at DUT pos. ; x_{upstream}-x_{downstream} [#mum]; Number of tracks" );

        sixdyHist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "sixdy", 4*_sixCut+1, -2*_sixCut, 2*_sixCut );
        sixdyHist->setTitle( "Triplet-driplet Y displacement at the DUT pos. ; y_{upstream}-y_{downstream} [#mum]; Number of tracks" );

        sixdx_noCut_Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "sixdx_noCut", 10*_sixCut+1, -5*_sixCut, 5*_sixCut );
        sixdx_noCut_Hist->setTitle( "Triplet-driplet X displacement at the DUT pos. (no trk Cuts applied) ; x_{upstream}-x_{downstream} [#mum]; Number of tracks" );

        sixdy_noCut_Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "sixdy_noCut", 10*_sixCut+1, -5*_sixCut, 5*_sixCut );
        sixdy_noCut_Hist->setTitle( "Triplet-driplet Y displacement at the DUT pos. (no trk Cuts applied) ; y_{upstream}-y_{downstream} [#mum]; Number of tracks" );

        sixdr_noCut_Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "sixdr_noCut", 10*_sixCut+1, -5*_sixCut, 5*_sixCut );
        sixdr_noCut_Hist->setTitle( "Triplet-driplet displacement at the DUT pos. (no trk Cuts applied) ; r_{upstream}-r_{downstream} [#mum]; Number of tracks" );

        sixkxHist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "sixkx", 100, -25, 25 );
        sixkxHist->setTitle( "Triplet-driplet X kink at the DUT pos.;x kink angle [mrad];Number of tracks" );

        sixkyHist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "sixky", 100, -25, 25 );
        sixkyHist->setTitle( "Triplet-driplet Y kink at the DUT pos.;y kink angle [mrad];Number of tracks" );

        // GBL:

        selxHist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "selx", 150, -15, 15 );
        selxHist->setTitle( "x at DUT, sel GBL;x [mm];tracks" );

        selyHist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "sely", 100, -10, 10 );
        selyHist->setTitle( "y at DUT, sel GBL;y [mm];tracks" );

        selaxHist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "selax", 100, -25, 25 );
        selaxHist->setTitle( "track angle x, sel GBL;x angle [mrad];tracks" );

        selayHist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "selay", 100, -25, 25 );
        selayHist->setTitle( "track angle y, sel GBL;y angle [mrad];tracks" );

        seldxHist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "seldx", 100, -5000, 5000 );
        seldxHist->setTitle( "track match x, sel GBL;#Deltax [#mum];tracks" );

        seldyHist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "seldy", 100, -2000, 2000 );
        seldyHist->setTitle( "track match y, sel GBL;#Deltay [#mum];tracks" );

        selkxHist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "selkx", 100, -25, 25 );
        selkxHist->setTitle( "kink x, sel GBL;kink x [mrad];tracks" );

        selkyHist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "selky", 100, -25, 25 );
        selkyHist->setTitle( "kink y, sel GBL;kink y [mrad];tracks" );

        seldx1Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "seldx1", 100, -2500, 2500 );
        seldx1Hist->setTitle( "triplet resid x at 1, sel GBL;#Deltax [#mum];tracks" );

        seldy1Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "seldy1", 100, -2500, 2500 );
        seldy1Hist->setTitle( "triplet resid y at 1, sel GBL;#Deltay [#mum];tracks" );

        seldx3Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "seldx3", 100, -2500, 2500 );
        seldx3Hist->setTitle( "triplet resid x at 3, sel GBL;#Deltax [#mum];tracks" );

        seldy3Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "seldy3", 100, -2500, 2500 );
        seldy3Hist->setTitle( "triplet resid y at 3, sel GBL;#Deltay [#mum];tracks" );

        seldx4Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "seldx4", 100, -2500, 2500 );
        seldx4Hist->setTitle( "triplet resid x at 4, sel GBL;#Deltax [#mum];tracks" );

        seldy4Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "seldy4", 100, -2500, 2500 );
        seldy4Hist->setTitle( "triplet resid y at 4, sel GBL;#Deltay [#mum];tracks" );

        seldx5Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "seldx5", 100, -5000, 5000 );
        seldx5Hist->setTitle( "triplet resid x at 5, sel GBL;#Deltax [#mum];tracks" );

        seldy5Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "seldy5", 100, -5000, 5000 );
        seldy5Hist->setTitle( "triplet resid y at 5, sel GBL;#Deltay [#mum];tracks" );

        seldx6Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "seldx6", 100, -5000, 5000 );
        seldx6Hist->setTitle( "triplet resid x at DUT, sel GBL;#Deltax [#mum];tracks" );

        seldy6Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "seldy6", 100, -5000, 5000 );
        seldy6Hist->setTitle( "triplet resid y at DUT, sel GBL;#Deltay [#mum];tracks" );

        gblndfHist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "gblndf", 16, -0.5, 15.5 );
        gblndfHist->setTitle( "GBL fit NDF;GBL NDF;tracks" );

        gblchi2Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "gblchi2", 101, 0, 100 );
        gblchi2Hist->setTitle( "GBL fit chi2;GBL chi2;tracks" );

        gblprbHist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "gblprb", 101, 0, 1 );
        gblprbHist->setTitle( "GBL fit probability;GBL fit probability;tracks" );

        // bad fits:

        badxHist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "badx", 150, -15, 15 );
        badxHist->setTitle( "x at DUT, bad GBL;x [mm];tracks" );

        badyHist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "bady", 100, -10, 10 );
        badyHist->setTitle( "y at DUT, bad GBL;y [mm];tracks" );

        badaxHist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "badax", 100, -25, 25 );
        badaxHist->setTitle( "track angle x, bad GBL;x angle [mrad];tracks" );

        badayHist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "baday", 100, -25, 25 );
        badayHist->setTitle( "track angle y, bad GBL;y angle [mrad];tracks" );

        baddxHist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "baddx", 100, -5000, 5000 );
        baddxHist->setTitle( "track match x, bad GBL;#Deltax [#mum];tracks" );

        baddyHist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "baddy", 100, -5000, 5000 );
        baddyHist->setTitle( "track match y, bad GBL;#Deltay [#mum];tracks" );

        badkxHist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "badkx", 100, -25, 25 );
        badkxHist->setTitle( "kink x, bad GBL;kink x [mrad];tracks" );

        badkyHist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "badky", 100, -25, 25 );
        badkyHist->setTitle( "kink y, bad GBL;kink y [mrad];tracks" );

        baddx1Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "baddx1", 100, -2500, 2500 );
        baddx1Hist->setTitle( "triplet resid x at 1, bad GBL;#Deltax [#mum];tracks" );

        baddy1Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "baddy1", 100, -2500, 2500 );
        baddy1Hist->setTitle( "triplet resid y at 1, bad GBL;#Deltay [#mum];tracks" );

        baddx3Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "baddx3", 100, -2500, 2500 );
        baddx3Hist->setTitle( "triplet resid x at 3, bad GBL;#Deltax [#mum];tracks" );

        baddy3Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "baddy3", 100, -2500, 2500 );
        baddy3Hist->setTitle( "triplet resid y at 3, bad GBL;#Deltay [#mum];tracks" );

        baddx4Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "baddx4", 100, -1500, 1500 );
        baddx4Hist->setTitle( "triplet resid x at 4, bad GBL;#Deltax [#mum];tracks" );

        baddy4Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "baddy4", 100, -1500, 1500 );
        baddy4Hist->setTitle( "triplet resid y at 4, bad GBL;#Deltay [#mum];tracks" );

        baddx5Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "baddx5", 100, -3000, 3000 );
        baddx5Hist->setTitle( "triplet resid x at 5, bad GBL;#Deltax [#mum];tracks" );

        baddy5Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "baddy5", 100, -3000, 3000 );
        baddy5Hist->setTitle( "triplet resid y at 5, bad GBL;#Deltay [#mum];tracks" );

        baddx6Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "baddx6", 100, -250, 250 );
        baddx6Hist->setTitle( "triplet resid x at DUT, bad GBL;#Deltax [#mum];tracks" );

        baddy6Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "baddy6", 100, -250, 250 );
        baddy6Hist->setTitle( "triplet resid y at DUT, bad GBL;#Deltay [#mum];tracks" );

        // good fits:

        goodx1Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "goodx1", 2*_sixCut+1, -_sixCut, _sixCut );
        goodx1Hist->setTitle( "Triplet residual x at 1, good GBL;#Deltax (residual) on plane 1 [#mum];Number of tracks" );

        goody1Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "goody1", 2*_sixCut+1, -_sixCut, _sixCut  );
        goody1Hist->setTitle( "Triplet residual y at 1, good GBL;#Deltay (residual) on plane 1 [#mum];Number of tracks" );

        goodx6Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "goodx6", 2*_sixCut+1, -_sixCut, _sixCut  );
        goodx6Hist->setTitle( "Triplet residual x at 6, good GBL;#Deltax (residual) on plane 1 [#mum];Number of tracks" );

        goody6Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "goody6", 2*_sixCut+1, -_sixCut, _sixCut  );
        goody6Hist->setTitle( "Triplet residual y at 6, good GBL;#Deltay (residual) on plane 1 [#mum];Number of tracks" );

        // look at fit:

        gblax0Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "gblax0", 100, -25, 25 );
        gblax0Hist->setTitle( "GBL angle at plane 0;x angle at plane 0 [mrad];tracks" );

        gbldx0Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "gbldx0", 100, -2500, 2500 );
        gbldx0Hist->setTitle( "GBL shift at plane 0;x shift at plane 0 [#mum];tracks" );

        gblrx0Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "gblrx0", 100, -250, 250 );
        gblrx0Hist->setTitle( "GBL resid at plane 0;x resid at plane 0 [#mum];tracks" );


        gblax1Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "gblax1", 2*_sixCut+1, -_sixCut, _sixCut  );
        gblax1Hist->setTitle( "GBL angle at plane 1;x angle at plane 1 [mrad];tracks" );

        gbldx1Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "gbldx1", 2*_sixCut+1, -_sixCut, _sixCut  );
        gbldx1Hist->setTitle( "GBL shift at plane 1;x shift at plane 1 [#mum];tracks" );

        gblrx1Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "gblrx1", 2*_sixCut+1, -_sixCut, _sixCut  );
        gblrx1Hist->setTitle( "GBL resid at plane 1;x resid at plane 1 [#mum];tracks" );


        gblax2Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "gblax2", 1000, -25, 25 );
        gblax2Hist->setTitle( "GBL angle at plane 2;x angle at plane 2 [mrad];tracks" );

        gbldx2Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "gbldx2", 2*_sixCut+1, -_sixCut, _sixCut  );
        gbldx2Hist->setTitle( "GBL shift at plane 2;x shift at plane 2 [#mum];tracks" );

        gblrx2Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "gblrx2", 2*_sixCut+1, -_sixCut, _sixCut );
        gblrx2Hist->setTitle( "GBL resid at plane 2;x resid at plane 2 [#mum];tracks" );


        gblax3Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "gblax3", 1000, -25, 25 );
        gblax3Hist->setTitle( "GBL angle at plane 3;x angle at plane 3 [mrad];tracks" );

        gbldx3Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "gbldx3", 2*_sixCut+1, -_sixCut, _sixCut  );
        gbldx3Hist->setTitle( "GBL shift at plane 3;x shift at plane 3 [#mum];tracks" );

        gblrx3Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "gblrx3", 2*_sixCut+1, -_sixCut, _sixCut  );
        gblrx3Hist->setTitle( "GBL resid at plane 3;x resid at plane 3 [#mum];tracks" );


        gblax4Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "gblax4", 1000, -25, 25 );
        gblax4Hist->setTitle( "GBL angle at plane 4;x angle at plane 4 [mrad];tracks" );

        gbldx4Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "gbldx4", 2*_sixCut+1, -_sixCut, _sixCut  );
        gbldx4Hist->setTitle( "GBL shift at plane 4;x shift at plane 4 [#mum];tracks" );

        gblrx4Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "gblrx4", 2*_sixCut+1, -_sixCut, _sixCut  );
        gblrx4Hist->setTitle( "GBL resid at plane 4;x resid at plane 4 [#mum];tracks" );


        gblax5Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "gblax5", 1000, -25, 25 );
        gblax5Hist->setTitle( "GBL angle at plane 5;x angle at plane 5 [mrad];tracks" );

        gbldx5Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "gbldx5", 2*_sixCut+1, -_sixCut, _sixCut  );
        gbldx5Hist->setTitle( "GBL shift at plane 5;x shift at plane 5 [#mum];tracks" );

        gblrx5Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "gblrx5", 2*_sixCut+1, -_sixCut, _sixCut  );
        gblrx5Hist->setTitle( "GBL resid at plane 5;x resid at plane 5 [#mum];tracks" );


        gblax6Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "gblax6", 1000, -25, 25 );
        gblax6Hist->setTitle( "GBL angle at DUT;x angle at DUT [mrad];tracks" );

        gbldx6Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "gbldx6", 2*_sixCut+1, -_sixCut, _sixCut );
        gbldx6Hist->setTitle( "GBL shift at DUT;x shift at DUT [#mum];tracks" );

        gbldy6Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "gbldy6", 2*_sixCut+1, -_sixCut, _sixCut );
        gbldy6Hist->setTitle( "GBL shift at DUT;y shift at DUT [#mum];tracks" );

        gblrx6Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "gblrx6", 2*_sixCut+1, -_sixCut, _sixCut );
        gblrx6Hist->setTitle( "GBL resid at DUT;x resid at DUT [#mum];tracks" );

        gblry6Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "gblry6", 2*_sixCut+1, -_sixCut, _sixCut );
        gblry6Hist->setTitle( "GBL resid at DUT;y resid at DUT [#mum];tracks" );


        gblkx1Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "gblkx1", 1000, -5, 5 );
        gblkx1Hist->setTitle( "GBL kink angle at plane 1;plane 1 kink [mrad];tracks" );

        gblkx2Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "gblkx2", 1000, -5, 5 );
        gblkx2Hist->setTitle( "GBL kink angle at plane 2;plane 2 kink [mrad];tracks" );

        gblkx3Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "gblkx3", 1000, -5, 5 );
        gblkx3Hist->setTitle( "GBL kink angle at plane 3;plane 3 kink [mrad];tracks" );

        gblkx4Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "gblkx4", 1000, -5, 5 );
        gblkx4Hist->setTitle( "GBL kink angle at plane 4;plane 4 kink [mrad];tracks" );

        gblkx5Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "gblkx5", 1000, -5, 5 );
        gblkx5Hist->setTitle( "GBL kink angle at plane 5;plane 5 kink [mrad];tracks" );

        gblkx6Hist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "gblkx6", 1000, -5, 5 );
        gblkx6Hist->setTitle( "GBL kink angle at plane 6;plane 6 kink [mrad];tracks" );

        nmHist = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "nm", 21, -0.5, 20.5 );
        nmHist->setTitle( "track matches;track matches;events" );

        hTripletsMatched_DUT = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "hTripletsMatched_DUT", 21, -0.5, 20.5 );
        hTripletsMatched_DUT->setTitle( "Tracks matched per event, DUT hits are presented;Number of matched tracks;Numb. of Events" );

        hTripletsMatched_DUT_2D = AIDAProcessor::histogramFactory(this)->
                createHistogram2D( "hTripletsMatched_DUTvsEv", _nmax, 0, _nmax - 1, 21, -0.5, 20.5 );
        hTripletsMatched_DUT_2D->setTitle( "Tracks matched per event vs. event, DUT hits are presented; Event number ; Number of matched tracks" );

        hTripletsMatched_REF = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "hTripletsMatched_REF", 21, -0.5, 20.5 );
        hTripletsMatched_REF->setTitle( "Tracks matched per event, REF hits are presented;Number of matched tracks;Numb. of Events" );

        hTripletsMatched_REF_2D = AIDAProcessor::histogramFactory(this)->
                createHistogram2D( "hTripletsMatched_REFvsEv", _nmax, 0, _nmax - 1, 21, -0.5, 20.5 );
        hTripletsMatched_REF_2D->setTitle( "Tracks matched per event vs. event, REF hits are presented; Event number ; Number of matched tracks" );

        hTripletsMatched_DUT_REF = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "hTripletsMatched_DUT_REF", 21, -0.5, 20.5 );
        hTripletsMatched_DUT_REF->setTitle( "Tracks matched per event, DUT and REF hits are presented;Number of matched tracks;Numb. of Events" );

        hTripletsMatched_DUT_REF_2D = AIDAProcessor::histogramFactory(this)->
                createHistogram2D( "hTripletsMatched_DUT_REFvsEv", _nmax, 0, _nmax - 1, 21, -0.5, 20.5 );
        hTripletsMatched_DUT_REF_2D->setTitle( "Tracks matched per event vs. event, DUT and REF hits are presented; Event number ; Number of matched tracks" );

        //tracks matched to REF
        hTrk_REF = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "hTrk_REF", 21, -0.5, 20.5 );
        hTrk_REF->setTitle( "Tracks matched to REF hits;Number of matched tracks; ; Numb. of Events" );

        hTrk_REF_dx_noCut = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "hTrk_REF_dx_noCut", 10*_drirefCut+1, -5*_drirefCut, 5*_drirefCut );
        hTrk_REF_dx_noCut->setTitle( "Track-REF X displacement (no DriRefCut applied) ; x_{REF}-x_{driplet} [#mum]; Number of tracks" );

        hTrk_REF_dy_noCut = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "hTrk_REF_dy_noCut", 10*_drirefCut+1, -5*_drirefCut, 5*_drirefCut );
        hTrk_REF_dy_noCut->setTitle( "Track-REF Y displacement (no DriRefCut applied) ; y_{REF}-y_{driplet} [#mum]; Number of tracks" );

        hTrk_REF_dr_noCut = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "hTrk_REF_dr_noCut", 10*_drirefCut+1, -5*_drirefCut, 5*_drirefCut );
        hTrk_REF_dr_noCut->setTitle( "Track-REF displacement (no DriRefCut applied) ; r_{REF}-r_{driplet} [#mum]; Number of tracks" );

        hTrk_REF_dx = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "hTrk_REF_dx", 10*_drirefCut+1, -5*_drirefCut, 5*_drirefCut );
        hTrk_REF_dx->setTitle( "Track-REF X displacement ; x_{REF}-x_{driplet} [#mum]; Number of tracks" );

        hTrk_REF_dy = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "hTrk_REF_dy", 10*_drirefCut+1, -5*_drirefCut, 5*_drirefCut );
        hTrk_REF_dy->setTitle( "Track-REF Y displacement ; y_{REF}-y_{driplet} [#mum]; Number of tracks" );

        //tracks matched to DUT
        hTrk_DUT = AIDAProcessor::histogramFactory(this)->
                createHistogram1D( "hTrk_DUT", 21, -0.5, 20.5 );
        hTrk_DUT->setTitle( "Tracks matched to DUT hits;Number of matched tracks; ; Numb. of Events" );


    }//try
    catch( lcio::Exception& e ) {

#ifdef EUTEL_INTERACTIVE
        streamlog_out( ERROR2 ) << "No AIDAProcessor initialized. Type q to exit or c to continue without histogramming" << endl;
        string answer;
        while ( true ) {
            streamlog_out( ERROR2 ) << "[q]/[c]" << endl;
            cin >> answer;
            transform( answer.begin(), answer.end(), answer.begin(), ::tolower );
            if( answer == "q" ) {
                exit(-1);
            }
            else if( answer == "c" )
                _histogramSwitch = false;
            break;
        }
#else
        streamlog_out( WARNING2 ) << "No AIDAProcessor initialized. Continue without histogramming" << endl;

#endif

    }
#endif

}

#endif // USE_GEAR
