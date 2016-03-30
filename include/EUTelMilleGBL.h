// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTELMULTILINEFIT_H
#define EUTELMULTILINEFIT_H

// built only if GEAR is available
#ifdef USE_GEAR
// eutelescope includes ".h"

// marlin includes ".h"
#include "marlin/Processor.h"

// marlin util includes
#include "mille/Mille.h"
#include "include/MilleBinary.h"

// gear includes <.h>
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>

// lcio includes <.h>
#include <EVENT/LCRunHeader.h>
#include <EVENT/LCEvent.h>
#include <IMPL/TrackerHitImpl.h>


// AIDA includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IBaseHistogram.h>
#endif

// system includes <>
#include <string>
#include <vector>
#include <map>

#if defined(USE_ROOT) || defined(MARLIN_USE_ROOT)
#include "TMinuit.h"
#include <TSystem.h>
#include <TMath.h>
class TMinuit;
#endif


namespace eutelescope {

 
  //Specify a Rectangular in a sensor
  class SensorRectangular {
	protected:
	// SensorID
	int sensor;
	// lowest pixel in X direction
	int A; 
	// lowest pixel in Y direction
	int B;
	// highest pixel in X direction
	int C;
	// highest pixel in Y direction
	int D;
	public:
	SensorRectangular(int s, int a, int b, int c, int d) : sensor(s), A(a), B(b), C(c), D(d) {};
	SensorRectangular() : sensor(0), A(0), B(0), C(0), D(0) {};
	int getSensor() const {return sensor; }
	//look if x and y are inside the foreseen rectangular
	bool isInside(int x, int y) const {  return (x >=A && x <=C && y >=B && y <=D); }
	void print() { streamlog_out(MESSAGE4) << "Sensor: " << sensor << ": (" << A << "|" << B << ") to (" << C << "|" << D << ")" << std::endl; }
	
  };
  
  class RectangularArray {
  protected:
    std::map<int,SensorRectangular > _rect;
	
  public:
    void addRectangular(SensorRectangular &s) { _rect[s.getSensor() ] = s;}
	
    bool isInside(int s, int x, int y) {
      std::map<int,SensorRectangular >::iterator it = _rect.find(s);
      if (it == _rect.end()) { // not in the map means no limit on this sensor -> always true
	return true;
      }
      SensorRectangular cSensor = _rect[s];
      return cSensor.isInside(x,y);
    }
		
  };


  //! Straight line fit processor
  /*!
   *
   */

  class EUTelMilleGBL : public marlin::Processor {

  public:
#if defined(USE_ROOT) || defined(MARLIN_USE_ROOT)
    class hit
    {
    public:
      hit()
      {
      }
      hit(double tx, double ty, double tz, double rx, double ry, double rz,int i)
      {
        x = tx;
        y = ty;
        z = tz;
        resolution_x = rx;
        resolution_y = ry;
        resolution_z = rz;
    
        planenumber = i;
      }
      double x;
      double y;
      double z;
      double resolution_x;
      double resolution_y;
      double resolution_z;
  
      int planenumber;
    };

#endif

    //! Variables for hit parameters
    class HitsInPlane {
    public:
      HitsInPlane(){
        measuredX = 0.0;
        measuredY = 0.0;
        measuredZ = 0.0;
      }
      HitsInPlane(double x, double y, double z)
      {
        measuredX = x;
        measuredY = y;
        measuredZ = z;
      }
      bool operator<(const HitsInPlane& b) const
      {
        return (measuredZ < b.measuredZ);
      }
      double measuredX;
      double measuredY;
      double measuredZ;
    };

    /*DP
    virtual void FitTrack(
                          int nPlanesFitter,
                          double xPosFitter[],
                          double yPosFitter[],
                          double zPosFitter[],
                          double xResFit[],
                          double yResFit[],
                          double chi2Fit[2],
                          double residXFit[],
                          double residYFit[],
                          double angleFit[2]
                          );
    */


    //recursive method which searches for track candidates

    /*
    virtual void findtracks(
                            std::vector<std::vector<int> > &indexarray, //resulting vector of hit indizes
                            std::vector<int> vec, //for internal use
                            std::vector<std::vector<EUTelMilleGBL::HitsInPlane> > &_hitsArray, //contains all hits for each plane
                            int i, //plane number
                            int y //hit index number
                            );
    */

    //! Returns a new instance of EUTelMilleGBL
    /*! This method returns a new instance of this processor.  It is
     *  called by Marlin execution framework and it shouldn't be
     *  called/used by the final user.
     *
     *  @return a new EUTelMilleGBL.
     */
    virtual Processor * newProcessor() {
      return new EUTelMilleGBL;
    }

    //DP virtual bool hitContainsHotPixels( TrackerHitImpl   * hit) ;


    //! Default constructor
    EUTelMilleGBL ();

    //! Called at the job beginning.
    /*! This is executed only once in the whole execution. It prints
     *  out the processor parameters and check that the GEAR
     *  environment is properly set up and accessible from Marlin.
     */
    virtual void init ();

    //! Called for every run.
    /*! It is called for every run, and consequently the run counter
     *  is incremented. The geometry ID of the file is compared with
     *  the one provided by the GEAR geometry description. In case the
     *  two are different, the user is asked to decide to quit or to
     *  continue with a description that might be wrong.
     *
     *  @param run the LCRunHeader of the this current run
     */
    virtual void processRunHeader (LCRunHeader * run);

    //! Called for first event per run
    /*! Reads hotpixel information from hotPixelCollection into hotPixelMap
     * to be used in the sensor exclusion area logic 
     */
    //DP virtual void  FillHotPixelMap(LCEvent *event);

    //! Called every event
    /*! This is called for each event in the file. Each element of the
     *  pulse collection is scanned and the center of the cluster is
     *  translated into the external frame of reference thanks to the
     *  GEAR geometry description.
     *
     *  The cluster center might be calculate using a standard linear
     *  charge center of gravity algortihm or applying a more
     *  sophisticated non linear eta function. This behaviour is
     *  regulated by the user from the steering file.
     *
     *  @throw UnknownDataTypeException if the cluster type is unknown
     *
     *  @param evt the current LCEvent event as passed by the
     *  ProcessMgr
     */
    virtual void processEvent (LCEvent * evt);


    //! Called after data processing.
    /*! This method is called when the loop on events is
     *  finished.
     */
    virtual void end();


    //! Histogram booking
    /*! Some control histograms are filled during this procedure in
     *  order to be able to perform easy check on the quality of the
     *  output hits and also to understand if the frame of reference
     *  conversion has been properly done. Of course this method is
     *  effectively doing something only in the case MARLIN_USE_AIDA.
     */
    void bookHistos();

    //!rootObject where to store histos for putting histos to sub-folders
    std::map<std::string , TObject * > _rootObjectMap;


  protected:


    //! Ordered sensor ID
    /*! Within the processor all the loops are done up to _nPlanes and
     *  according to their position along the Z axis (beam axis).
     *
     *  This vector is containing the sensorID sorted according to the
     *  same rule.
     */
    std::vector< int > _orderedSensorID;
    std::vector< int > _orderedSensorID_wo_excluded;



    //! TrackerHit collection name
    /*! Input collection with hits.
     */
    std::vector<std::string > _hitCollectionName;

    //! TRACK collection name
    /*! Output collection with fitted tracks.
     */
    std::string _trackCollectionName;

    //! Hot pixel collection name.
    /*! 
     * this collection is saved in a db file to be used at the clustering level
     */
    //DP std::string _hotPixelCollectionName;

    //! Vector of map arrays, keeps record of hit pixels 
    /*! The vector elements are sorted by Detector ID
     *  For each Detector unique ID element a map of pixels is created. 
     *  first level key   sensor unique 
     *              value sensor map
     *  sensor map key    unique row number
     *             value  vector of column numbers.
     */
    
    //DP std::map<std::string, bool > _hotPixelMap;


    // parameters
    //shall mille be initialized or not
    int _milleInit;

    //DUT parameters:
    std::string _sensetiveAxis;

    //analysis parameters:
    double _chi2tondfcut;
    int DUT_scatterer;
    double DUTX0;
    int _aligntype;

    float _distanceMax;
    std::vector<float> _distanceMaxVec;
    std::vector<int > _excludePlanes; //only for internal usage
    std::vector<int > _excludePlanes_sensorIDs; //this is going to be
                                                //set by the user.
    std::vector<int > _FixedPlanes; //only for internal usage
    std::vector<int > _FixedPlanes_sensorIDs; //this is going to be
    //GEAR planes ID
    int _GEARidDUT;
    int _GEARidREF;
    //set by the user.

    double _eBeam; // DP
    double _triCut;
    double _driCut;
    double _sixCut;
    double _tridutCut;
    double _dridutCut;
    double _drirefCut;

    int _maxTrackCandidates;
    int _maxTrackCandidatesTotal;

    std::string _binaryFilename;

    float _telescopeResolution;
    float _DUTResolution;
    float _REFResolution;
    int _onlySingleHitEvents;
    int _onlySingleTrackEvents;
    int _alignMode;
    int _useResidualCuts;

    std::vector<float > _residualsXMin;
    std::vector<float > _residualsYMin;
    std::vector<float > _residualsXMax;
    std::vector<float > _residualsYMax;

    std::vector<float> _resolutionX;
    std::vector<float> _resolutionY;
    std::vector<float> _resolutionZ;

    std::vector<int> _FixParameter;

    int _generatePedeSteerfile;
    std::string _pedeSteerfileName;
    int _runPede;
    int _usePedeUserStartValues;
    std::vector<float > _pedeUserStartValuesX;
    std::vector<float > _pedeUserStartValuesY;
    std::vector<float > _pedeUserStartValuesZ;
    
    std::vector<float > _pedeUserStartValuesAlpha;
    std::vector<float > _pedeUserStartValuesBeta;
    std::vector<float > _pedeUserStartValuesGamma;

    int _inputMode;
    float _testModeSensorResolution;
    float _testModeXTrackSlope;
    float _testModeYTrackSlope;

    std::vector<float > _testModeSensorZPositions;

    std::vector<float > _testModeSensorXShifts;
    std::vector<float > _testModeSensorYShifts;
    std::vector<float > _testModeSensorGamma;
    std::vector<float > _testModeSensorAlpha;
    std::vector<float > _testModeSensorBeta;

    std::string _alignmentConstantLCIOFile;
    std::string _alignmentConstantCollectionName;

    std::vector<int> _useSensorRectangular;

  private:

    //! Run number
    int _iRun;

    //! Event number
    int _iEvt;

    //! Max number of events to be processed
    int _nmax;

    // Excluded planes
    int _nExcludePlanes;

    // Statistics
    int _nMilleTracksTotal;
    int _nMilleTracksperEvent;
    int ntrk_REF_Event;
    int ntrk_DUT_Event;
    int ntrk_REF_Total;
    int ntrk_DUT_Total;

    // Mille
    //DP Mille * _mille;

    //! Conversion ID map.
    /*! In the data file, each cluster is tagged with a detector ID
     *  identify the sensor it belongs to. In the geometry
     *  description, there are along with the sensors also "passive"
     *  layers and other stuff. Those are identify by a layerindex. So
     *  we need a conversion table to go from the detectorID to the
     *  layerindex.
     */
    std::map< int, int > _conversionIdMap;

    //! Silicon planes parameters as described in GEAR
    /*! This structure actually contains the following:
     *  @li A reference to the telescope geoemtry and layout
     *  @li An integer number saying if the telescope is w/ or w/o DUT
     *  @li An integer number saying the number of planes in the
     *  telescope.
     *
     *  This object is provided by GEAR during the init() phase and
     *  stored here for local use.
     */
    gear::SiPlanesParameters * _siPlanesParameters;

    //! Silicon plane layer layout
    /*! This is the real geoemetry description. For each layer
     *  composing the telescope the relevant information are
     *  available.
     *
     *  This object is taken from the _siPlanesParameters during the
     *  init() phase and stored for local use
     */
    gear::SiPlanesLayerLayout * _siPlanesLayerLayout;

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
    //! AIDA histogram map
    /*! Instead of putting several pointers to AIDA histograms as
     *  class members, histograms are booked in the init() method and
     *  their pointers are inserted into this map keyed by their
     *  names.
     *  The histogram filling can proceed recalling an object through
     *  its name
     */
    std::map<std::string, AIDA::IBaseHistogram * > _aidaHistoMap;

    static std::string _numberTracksLocalname;

    static std::string _chi2XLocalname;
    static std::string _chi2YLocalname;

    static std::string _residualXLocalname;
    static std::string _residualYLocalname;
    static std::string _residualZLocalname;
#endif
    //numb of all planes taken from the gear file
    int _nPlanes;
    //numb. of planes which hits are being used for analysis
    int nPl;

    std::vector<std::vector<double> > _xPos;
    std::vector<std::vector<double> > _yPos;
    std::vector<std::vector<double> > _zPos;

    double * _xPosHere;
    double * _yPosHere;
    double * _zPosHere;
    double * _waferResidX;
    double * _waferResidY;
    double * _waferResidZ;
    double * _telescopeResolX;
    double * _telescopeResolY;
    double * _telescopeResolZ;
    double * _xFitPos;
    double * _yFitPos;

    std::vector<double> _siPlaneZPosition;

    //! Fill histogram switch
    /*! Only for debug reason
     */
    bool _histogramSwitch;

    //! Limits the pixels on each sensor-plane to a sub-rectangular  
    RectangularArray _rect;

  };

  //! A global instance of the processor
  EUTelMilleGBL gEUTelMilleGBL;

}
#endif
#endif
