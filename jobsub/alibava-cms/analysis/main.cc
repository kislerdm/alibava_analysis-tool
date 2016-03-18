//////////////////////
//
//	Analysis of Testbeam data
//
//	This loops over some ntuples created with EUTelescope and then does magic
//
//////////////////////

// start with:		root +x -q -l -b 'main.cc("runselection")'
// compile with:	g++ -I `root-config --incdir` -o asdf main.cc `root-config --libs`


// differences: g0 lower limit to 1 in langaus as it sometimes goes to fail low values...

/*
 *
 * Rescue Pig to the rescue!

                                 _
    _._ _..._ .-',     _.._(`))
   '-. `     '  /-._.-'    ',/
      )         \            '.
     / _    _    |             \
    |  a    a    /              |
    \   .-.                     ;  
     '-('' ).-'       ,'       ;
        '-;           |      .'
           \           \    /
           | 7  .__  _.-\   \
           | |  |  ``/  /`  /
          /,_|  |   /,_/   /
             /,_/      '`-'


*/




//C++ headers
#include <iostream>
#include <algorithm>
#include <fstream>
#include <cstring>
#include <string>
#include <vector>
#include <list>
#include <ctime>
#include <cmath>
#include <iomanip>
#include <stdio.h>
#include <utility>
#include <map>
#include <sstream>

//Root headers
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TF1.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TAxis.h"
#include "TMath.h"
#include "TLine.h"
#include "TLegend.h"
#include "TTree.h"
#include "TObject.h"
#include "TH1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TGaxis.h"
#include "TFile.h"
#include "TSpectrum.h"
#include "TStopwatch.h"
#include "TLatex.h"

#define PI 3.14159265

using namespace std;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	CUTS
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// the number of runs we read from the runlist -> 230 max for all runs

// 91 for all 0deg
// 11 for good n
// 33 for 0deg p
// 25 for 0deg y
// 33 for 0deg n (all)

// 54 for 0,25 p
// 41 for 0,25 y
# define _runcount 231

// the maximum number of tracks we allow in an event
int _cut_maxtracksperevent = 20;

// only singletrack events?
bool _cut_onlysingletrackevents = false;

// the divergence cut in x
float _cut_divergencex = 1.0;

// do we require the track signal to be highest in the middle 3 strips?
bool _cut_highchannelmiddle = true;

// the +- distance for a track to be considered "on" a strip
// 0.25 * 0.08 pitch -> 20 um
float _cut_onstrip = 0.25;

// the debug level
int _debug = 1;

// the maximum number of telescope triplet tracks in a run
int _maxtotaltracks = 99999999;
//int _maxtotaltracks = 20;

// the telescope resolution in um - this is used to plot a line in the residual vs voltage plot, this is not really a cut
float _cut_telescoperesolution = 5.0;

// the noise range to plot
float _cut_minnoise = 0.0;
float _cut_maxnoise = 50.0;

// the residual range to plot
float _cut_minXresidual = 0.0;
float _cut_maxXresidual = 0.5;
float _cut_minYresidual = 0.02;
float _cut_maxYresidual = 0.04;

// the signal range to plot
float _cut_minsignal = 0.0;
float _cut_maxsignal = 100.0;

// do we want an extra comparison between the eta integral of matched and unmatched etas?
bool _cut_drawetaintegralcomparison = false;

// the eta value above (and 1- below) we consider as chargesharing
float _cut_eta_chargeshare = 0.2;

// the maximum channel noise value to calculate rghs
// if lots of channels in a run are over this noise then rgh calculation makes no sense, since the sensor is too noisy
float _cut_maxrghnoise = 75.0;

// the maximum rgh percent value to plot
float _cut_maxrghplot = 10.0;

// the maximum cluster count per track we want to plot
float _cut_maxclustercount = 1.0;

// the cut times noise to define charge sharing
float _cut_chargesharing = 1.5;

// apply a temperature correction?
bool _cut_applytempcorrection = true;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	Global settings
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// the sensor features
// pitch in mm
float _pitchx = 0.018402778;
float _pitchy = 0.08;

// absolute error on the voltage in V:
float _voltageerror = 5.0;

// absolute error on the clustercount per run:
float _clustercounterror = 0.0;

// the runlist if nothing is selected
string _runfile = "../runlist.csv";

// the suffix of our root files
string _filesuffix = "-alibava-tracking-3.root";

// the suffix of the root cluster files
string _clusterfilesuffix = "-alibava-clustering-2.root";

// the vectors to store runlist information in
// _total gives the number of distinct entries in the vectors for plotting etc.
std::vector<int> _counter;
std::vector<int> _runnumber;
std::vector<int> _pedestal;
std::vector<int> _thickness;
std::vector<int> _thickness_list;
int _thickness_total = 0;
std::vector<int> _dutrotation;
std::vector<int> _dutrotation_list;
int _dutrotation_total = 0;
std::vector<int> _biasvoltage;
std::vector<int> _biasvoltage_list;
int _biasvoltage_total = 0;
std::vector<int> _temperature;
std::vector<int> _temperature_list;
int _temperature_total = 0;
std::vector<double> _sensorcurrent;
std::vector<double> _sensorcurrent_list;
int _sensorcurrent_total = 0;
std::vector<int> _polarity;
std::vector<int> _polarity_list;
int _polarity_total = 0;
std::vector<string> _sensorname;
std::vector<string> _sensorname_list;
int _sensorname_total = 0;
std::vector<char> _type;
std::vector<char> _type_list;
int _type_total = 0;
std::vector<float> _irradfluence;
std::vector<float> _irradfluence_list;
int _irradfluence_total = 0;
int _channels[_runcount][128];
float _noise[_runcount][128];
std::vector<int> _goodtracks;
std::vector<float> _sensorminx;
std::vector<float> _sensormaxx;
std::vector<float> _sensormintdc;
std::vector<float> _sensormaxtdc;
std::vector<float> _sensoralignx;
std::vector<float> _sensoraligny;
std::vector<float> _sensoralignz;
std::vector<float> _sensoraligna;
std::vector<float> _sensoralignb;
std::vector<float> _sensoralignc;

// the value which is written to the tuple by EUTelescope if there is no entry
int _missingvalue = -999;

// map for the root objects
std::map < string , TObject * > _rootObjectMap;

// the global run time
TStopwatch _totaltime;

// how many we have read from file
int _actual_runcount = 0;

string _thefilename = "test.root";

TFile* _outputFile;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// declarations

void readrunlist();

void bookhistos();

void openfile();

TF1* lanGausFit(TH1* inHist, double negSigmaFit, double posSigmaFit);

TF1* gausLanGausFitFixGausNoise(TH1* inHist, double negSigmaFit, double posSigmaFit, double mean, double sigma);



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// let's go!

//int main(string typetorun){
  
  
  
int main(int argc, char** argv){

 string  typetorun = "fail";
 stringstream astream;
  if (argc>1)
  {
    astream << argv[1];
    typetorun = astream.str();
    cout << " you selected " << typetorun << endl;
  } else {
    cout << "Please select a type to run!" << endl;
  }


  // root pre stuff
  //TApplication* test = new TApplication();

  gROOT->SetBatch();
  gStyle->SetLabelSize(0.035, "x");
  gStyle->SetLabelSize(0.035, "y");
  gStyle->SetTitleSize(0.05, "x");
  gStyle->SetTitleSize(0.05, "y");
  gStyle->SetTitleOffset(0.95, "x");
  gStyle->SetTitleOffset(0.95, "y");
  gStyle->SetOptFit(1111);

  _totaltime.Start();

  cout << "##############################################" << endl;
  cout << " " << endl;
  cout << "Hello!" << endl;
  cout << " " << endl;

  // get the bulkselection from user:
  if (typetorun == "n")
  {

    cout << "You selected all n-types!" << endl;
    _runfile = "../runlists/n.csv";
    _thefilename = "epi_n.root";

  } else if (typetorun == "p") {

    cout << "You selected all p-stop-types!" << endl;
    _runfile = "../runlists/p.csv";
    _thefilename = "epi_p.root";

  } else if (typetorun == "y") {

    cout << "You selected all p-spray-types!" << endl;
    _runfile = "../runlists/y.csv";
    _thefilename = "epi_y.root";

  } else if (typetorun == "n0") {

    cout << "You selected 0° n!" << endl;
    _runfile = "../runlists/0_n.csv";
    _thefilename = "epi_n0.root";

  } else if (typetorun == "p0") {

    cout << "You selected 0° p!" << endl;
    _runfile = "../runlists/0_p.csv";
    _thefilename = "epi_p0.root";

  } else if (typetorun == "y0") {

    cout << "You selected 0° y!" << endl;
    _runfile = "../runlists/0_y.csv";
    _thefilename = "epi_y0.root";

  } else if (typetorun == "n25") {

    cout << "You selected 25° n!" << endl;
    _runfile = "../runlists/25_n.csv";
    _thefilename = "epi_n25.root";

  } else if (typetorun == "p25") {

    cout << "You selected 25° p!" << endl;
    _runfile = "../runlists/25_p.csv";
    _thefilename = "epi_p25.root";

  } else if (typetorun == "y25") {

    cout << "You selected 25° y!" << endl;
    _runfile = "../runlists/25_y.csv";
    _thefilename = "epi_y25.root";

  } else if (typetorun == "n51") {

    cout << "You selected 51° n!" << endl;
    _runfile = "../runlists/51_n.csv";
    _thefilename = "epi_n51.root";

  } else if (typetorun == "p51") {

    cout << "You selected 51° p!" << endl;
    _runfile = "../runlists/51_p.csv";
    _thefilename = "epi_p51.root";

  } else if (typetorun == "y51") {

    cout << "You selected 51° y!" << endl;
    _runfile = "../runlists/51_y.csv";
    _thefilename = "epi_y51.root";

  } else if (typetorun == "npy025") {

    cout << "You selected all 0° and 25°" << endl;
    _runfile = "../runlists/025_npy.csv";
    _thefilename = "epi_npy025.root";
    
  } else if (typetorun == "n025") {

    cout << "You selected all n025 runs!" << endl;
    _runfile = "../runlists/n025.csv";
    _thefilename = "epi_n025.root";

  } else if (typetorun == "p025") {

    cout << "You selected all p025 runs!" << endl;
    _runfile = "../runlists/p025.csv";
    _thefilename = "epi_p025.root";
    
  } else if (typetorun == "y025") {

    cout << "You selected all y025 runs!" << endl;
    _runfile = "../runlists/y025.csv";
    _thefilename = "epi_y025.root";


  } else if (typetorun == "a") {

    cout << "You selected all runs!" << endl;
    _runfile = "../runlist.csv";
    _thefilename = "epi_all.root";


  } else if (typetorun == "t") {

    cout << "You selected test!" << endl;
    _runfile = "../runlists/test.csv";
    _thefilename = "test.root";

  } else {
    cout << "You selected nothing!" << endl;
  }

  // the output
  _outputFile = new TFile(_thefilename.c_str(), "RECREATE");

  // read the runlist into the vectors
  readrunlist();

  // book ALL histograms
  bookhistos();

  // open the root files and fill the histograms
  openfile();

  // done
  _totaltime.Stop();

  cout << " " << endl;
  cout << "##############################################" << endl;
  cout << " " << endl;
  cout << "Time elapsed: " << _totaltime.RealTime() << " s!" << endl;
  cout << " " << endl;
  cout << "Goodbye!" << endl;
  return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// this reads the runlist to get filenames of tuples
void readrunlist()
{

  cout << "##############################################" << endl;
  cout << " " << endl;
  cout << "Reading " << _runfile << " !" << endl;
  cout << " " << endl;

  // the strings we parse into
  string fail;
  string sensor;
  string fluence;
  string rotation;
  string voltage;
  string temp;
  string current;
  string minx;
  string maxx;
  string tdcmin;
  string tdcmax;
  string bonds;
  string pede;

  // the ints
  int ped;
  int run;
  int pol;
  int failint;

  // char
  char failchar[200];

  // the floats
  float flu;
  float minxval;
  float maxxval;
  float mintdcval;
  float maxtdcval;

  // the line to read
  string line;

  // the stream
  ifstream fileRead;

  // open
  fileRead.open(_runfile.c_str());

  // count on the found sensors
  int i = 0;

  // flag for finding information
  bool foundinfo = false;
  bool foundpede = false;

  // loop over the lines in the file
  while (std::getline(fileRead, line))
  {
    string startpart = line.substr(0,2);
    if (startpart == "#!")
    {

      // we found something
      foundinfo = true;

      // the input string of this line
      string input;
      std::istringstream iss(line);
      input = iss.str();

      // the delimiter of the comment (space)
      string delimiter = " ";

      // the position we found of the delimiter
      size_t pos = 0;

      // the iterator for counting elements
      std::vector<int>::iterator it;

      // first the #!
      pos = input.find(delimiter);
      fail = input.substr(0, pos);
      input.erase(0, pos + delimiter.length());

      // now the sensor
      // push back into _polarity, sensor name, _thickness and _type
      pos = input.find(delimiter);
      sensor = input.substr(0, pos);
      input.erase(0, pos + delimiter.length());
      int thick = 0;
      thick = atoi(sensor.c_str());
      _thickness.push_back(thick);
      if (sensor.at(3) == 'N')
      {
	pol = 1;
      } else {
	pol = -1;
      }

      it = find (_polarity.begin(), _polarity.end(), pol);
      if (it != _polarity.end())
      {
	if (_debug <=1)
	{
	  cout << "Polarity " << pol << " already in list!" << endl;
	}
      } else {
	_polarity_total++;
	_polarity_list.push_back(pol);
	if (_debug <=2)
	{
	  cout << "Polarity " << pol << " not in list! Total polarities now " << _polarity_total << endl;
	}
      }

      // string iterator
      std::vector<string>::iterator its;

      its = find (_sensorname.begin(), _sensorname.end(), sensor);
      if (its != _sensorname.end())
      {
	if (_debug <=1)
	{
	  cout << "Sensor " << sensor << " already in list!" << endl;
	}
      } else {
	_sensorname_total++;
	_sensorname_list.push_back(sensor);
	if (_debug <=2)
	{
	  cout << "Sensor " << sensor << " not in list! Total sensors now " << _sensorname_total << endl;
	}
      }

      // char iterator
      std::vector<char>::iterator itc;

      itc = find (_type.begin(), _type.end(), sensor.at(3));
      if (itc != _type.end())
      {
	if (_debug <=1)
	{
	  cout << "Type " << sensor.at(3) << " already in list!" << endl;
	}
      } else {
	_type_total++;
	_type_list.push_back(sensor.at(3));
	if (_debug <=2)
	{
	  cout << "Type " << sensor.at(3) << " not in list! Total types now " << _type_total << endl;
	}
      }

      _polarity.push_back(pol);
      _sensorname.push_back(sensor);
      _type.push_back(sensor.at(3));

      // now the fluence
      pos = input.find(delimiter);
      fluence = input.substr(0, pos);
      input.erase(0, pos + delimiter.length());
      if (fluence == "unirr,")
      {
	flu = 0.0;
      } else {
	flu = atof(fluence.c_str());
      }

      // float iterator
      std::vector<float>::iterator itf;

      itf = find (_irradfluence.begin(), _irradfluence.end(), flu);
      if (itf != _irradfluence.end())
      {
	if (_debug <=1)
	{
	  cout << "Irradiation " << flu << " already in list!" << endl;
	}
      } else {
	_irradfluence_total++;
	_irradfluence_list.push_back(flu);
	if (_debug <=2)
	{
	  cout << "Irradiation " << flu << " not in list! Total fluences now " << _irradfluence_total << endl;
	}
      }

      _irradfluence.push_back(flu);

      // now the rotation
      pos = input.find(delimiter);
      rotation = input.substr(0, pos);
      input.erase(0, pos + delimiter.length());
      int rot = 0;
      rot = atoi(rotation.c_str());

      it = find (_dutrotation.begin(), _dutrotation.end(), rot);
      if (it != _dutrotation.end())
      {
	if (_debug <=1)
	{
	  cout << "Rotation " << rot << "deg already in list!" << endl;
	}
      } else {
	_dutrotation_total++;
	_dutrotation_list.push_back(rot);
	if (_debug <=2)
	{
	  cout << "Rotation " << rot << "deg not in list! Total rotations now " << _dutrotation_total << endl;
	}
      }

      _dutrotation.push_back(rot);

      // now the voltage
      pos = input.find(delimiter);
      voltage = input.substr(0, pos);
      input.erase(0, pos + delimiter.length());
      int vol = 0;
      vol = atoi(voltage.c_str());


      it = find (_biasvoltage.begin(), _biasvoltage.end(), vol);
      if (it != _biasvoltage.end())
      {
	if (_debug <=1)
	{
	  cout << "Voltage " << vol << "V already in list!" << endl;
	}
      } else {
	_biasvoltage_total++;
	_biasvoltage_list.push_back(vol);
	if (_debug <=2)
	{
	  cout << "Voltage " << vol << "V not in list! Total voltages now " << _biasvoltage_total << endl;
	}
      }

      _biasvoltage.push_back(vol);


      // now the temp
      pos = input.find(delimiter);
      temp = input.substr(0, pos);
      input.erase(0, pos + delimiter.length());
      int tem = 0;
      tem = atoi(temp.c_str());

      it = find (_temperature.begin(), _temperature.end(), tem);
      if (it != _temperature.end())
      {
	if (_debug <=1)
	{
	  cout << "Temperature " << tem << "°C already in list!" << endl;
	}
      } else {
	_temperature_total++;
	_temperature_list.push_back(tem);
	if (_debug <=2)
	{
	  cout << "Temperature " << tem << "°C not in list! Total temperatures now " << _temperature_total << endl;
	}
      }

      _temperature.push_back(tem);


      // now the current
      pos = input.find(delimiter);
      current = input.substr(0, pos);
      input.erase(0, pos + delimiter.length());
      float cur = 0;
      cur = atof(current.c_str());

      // double iterator
      std::vector<double>::iterator itd;

      itd = find (_sensorcurrent.begin(), _sensorcurrent.end(), cur);
      if (itd != _sensorcurrent.end())
      {
	if (_debug <=1)
	{
	  cout << "Sensor Current " << cur << "mA already in list!" << endl;
	}
      } else {
	_sensorcurrent_total++;
	_sensorcurrent_list.push_back(cur);
	if (_debug <=2)
	{
	  cout << "Sensor Current " << cur << "mA not in list! Total sensor currents now " << _sensorcurrent_total << endl;
	}
      }

      _sensorcurrent.push_back(cur);


      // now the min/max x

      pos = input.find(delimiter);
      minx = input.substr(0, pos);
      input.erase(0, pos + delimiter.length());
      minxval = 0;
      minxval = atof(minx.c_str());
      _sensorminx.push_back(minxval);

      pos = input.find(delimiter);
      maxx = input.substr(0, pos);
      input.erase(0, pos + delimiter.length());
      maxxval = 0;
      maxxval = atof(maxx.c_str());
      _sensormaxx.push_back(maxxval);


      // the min/max tdc
      pos = input.find(delimiter);
      tdcmin = input.substr(0, pos);
      input.erase(0, pos + delimiter.length());
      mintdcval = 0;
      mintdcval = atof(tdcmin.c_str());
      _sensormintdc.push_back(mintdcval);

      pos = input.find(delimiter);
      tdcmax = input.substr(0, pos);
      input.erase(0, pos + delimiter.length());
      maxtdcval = 0;
      maxtdcval = atof(tdcmax.c_str());
      _sensormaxtdc.push_back(maxtdcval);


      // now the aligments
      string tempstring;
      float tempfloat;
      pos = input.find(delimiter);
      tempstring = input.substr(0, pos);
      input.erase(0, pos + delimiter.length());
      tempfloat = 0;
      tempfloat = atof(tempstring.c_str());
      _sensoralignx.push_back(tempfloat);

      pos = input.find(delimiter);
      tempstring = input.substr(0, pos);
      input.erase(0, pos + delimiter.length());
      tempfloat = 0;
      tempfloat = atof(tempstring.c_str());
      _sensoraligny.push_back(tempfloat);

      pos = input.find(delimiter);
      tempstring = input.substr(0, pos);
      input.erase(0, pos + delimiter.length());
      tempfloat = 0;
      tempfloat = atof(tempstring.c_str());
      _sensoralignz.push_back(tempfloat);

      pos = input.find(delimiter);
      tempstring = input.substr(0, pos);
      input.erase(0, pos + delimiter.length());
      tempfloat = 0;
      tempfloat = atof(tempstring.c_str());
      _sensoraligna.push_back(tempfloat);

      pos = input.find(delimiter);
      tempstring = input.substr(0, pos);
      input.erase(0, pos + delimiter.length());
      tempfloat = 0;
      tempfloat = atof(tempstring.c_str());
      _sensoralignb.push_back(tempfloat);

      pos = input.find(delimiter);
      tempstring = input.substr(0, pos);
      input.erase(0, pos + delimiter.length());
      tempfloat = 0;
      tempfloat = atof(tempstring.c_str());
      _sensoralignc.push_back(tempfloat);

      // output
      if (_debug <= 2)
      {
	cout << "Sensor " << i << ": " << sensor << " " << fluence << " " << rotation << " " << voltage << " " << temp << " " << current << endl;
	cout << " Range: X: " << minxval << " " << maxxval << " T: " << mintdcval << " " << maxtdcval << endl;
      }

    } else {

      // if the line is not a special comment, we look if the bools are true.
      // first run, then ped, so we don't find things twice!

      if (foundpede == true)
      {

	// read run numbers
	std::istringstream iss(line);
	while (iss >> run)
	{
	  if (_debug <= 3)
	  {
	    cout << "Run " << i << ": " << run << endl;
	  }
	}

	// set back to false
	foundpede = false;

	_counter.push_back(i);
	_runnumber.push_back(run);

	i++;
      }

      if (foundinfo == true)
      {
	string input;
	std::istringstream iss(line);
	input = iss.str();
	string delimiter = ",";
	size_t pos = input.find(delimiter);
	pede = input.substr(0, pos);
	input.erase(0, pos + delimiter.length());
	ped = atoi(pede.c_str());
	_pedestal.push_back(ped);

	if (_debug <= 1)
	{
	  cout << "Ped " << i << " " << ped << endl;
	}

	delimiter = "'";
	pos = input.find(delimiter);
	fail = input.substr(0, pos);
	input.erase(0, pos + delimiter.length());

	pos = input.find(delimiter);
	bonds = input.substr(0, pos);
	input.erase(0, pos + delimiter.length());

	// set all channels to off
	for (int ii=0;ii<128;ii++)
	{
	  _channels[i][ii] = 0;
	}

	int start, end;
	pos=0;
	while (pos<bonds.size()-1)
	{

	  delimiter = ":";
	  pos = bonds.find(delimiter);
	  fail = bonds.substr(0, pos);
	  bonds.erase(0, pos + delimiter.length());
	  start = atoi(bonds.c_str());

	  delimiter = "-";
	  pos = bonds.find(delimiter);
	  fail = bonds.substr(0, pos);
	  bonds.erase(0, pos + delimiter.length());
	  end = atoi(bonds.c_str());

	  for (int ii=start;ii<=end;ii++)
	  {
	    _channels[i][ii] = 1;
	  }

	}

	if (_debug <= 2)
	{
	  cout << "Good Channels:" << endl;
	  cout << "  0 - 31: ";
	  for (int ii=0;ii<32;ii++)
	  {
	    cout << _channels[i][ii];
	  }
	  cout << endl << " 32 - 63: ";
	  for (int ii=32;ii<64;ii++)
	  {
	    cout << _channels[i][ii];
	  }
	  cout << endl << " 64 - 95: ";
	  for (int ii=64;ii<96;ii++)
	  {
	    cout << _channels[i][ii];
	  }
	  cout << endl << "96 - 127: ";
	  for (int ii=96;ii<128;ii++)
	  {
	    cout << _channels[i][ii];
	  }
	  cout << endl;
	}

	// set back to false
	foundinfo = false;
	foundpede = true;

	_pedestal.push_back(ped);
      }

    }

    // limit the reading to save time...
    if (i == _runcount)
      break;
  }

  cout << "Done reading after " << i << " found sensors!" << endl;
  _actual_runcount = i;
  fileRead.close();

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// this reads the nth file in the main _runnumber vector
void openfile()
{

  cout << "##############################################" << endl;
  cout << " " << endl;
  cout << "Opening root files and filling histograms!" << endl;
  cout << " " << endl;

  // the names and titles for the histograms
  char name[100];
  char title[300];

  // the string for the name
  string histoname;
  stringstream sstream;

  // the fractions used later for mod pitch
  double fractpartx, fractparty, intpartx, intparty;

  // this will count the entries in a graph
  int tempcount = 0;

  // the vector of all graphs if we do vs voltage
  std::vector<int> voltagegraphs_total;
  // the individual ones
  std::vector<int> voltagegraphs;
  // the point count in each graph
  int voltagepoint[500];
  for (int i = 0; i< 500; i++)
  {
    voltagepoint[i] = 0;
  }

  // these global histos will get the alignments of each run
  histoname = "xshift";
  sstream << histoname;
  TH1D* xshift = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
  sstream.str(string());

  histoname = "yshift";
  sstream << histoname;
  TH1D* yshift = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
  sstream.str(string());

  histoname = "zshift";
  sstream << histoname;
  TH1D* zshift = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
  sstream.str(string());

  histoname = "ashift";
  sstream << histoname;
  TH1D* ashift = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
  sstream.str(string());

  histoname = "bshift";
  sstream << histoname;
  TH1D* bshift = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
  sstream.str(string());

  histoname = "cshift";
  sstream << histoname;
  TH1D* cshift = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
  sstream.str(string());


  // now open each file!
  for ( int irun = 0 ; irun < _actual_runcount ; irun++ )
  {

    // set the noise map in this run to zero
    for (int ii=0;ii<128;ii++)
    {
      _noise[irun][ii] = 0.0;
    }

    // get the actual filenumber from the vector
    int filenumbertoopen = _runnumber.at(irun);

    // output
    if (_debug <= 3)
    {
      cout << " " << endl;
      cout << "**********************************************" << endl;
      cout << "Reading run " << irun << " of " << _actual_runcount - 1 << " , nr. " << _runnumber.at(irun) << endl;
    }

    // base string
    string filename("../output/histograms/000");
    string cfilename("../output/histograms/000");

    // leading zero
    if (filenumbertoopen <= 99)
    {
      filename+="0";
      cfilename+="0";
    }

    // move the number to a string
    stringstream ss;
    ss << filenumbertoopen;
    string numbertoopen = ss.str();

    // construct the filename
    filename+=numbertoopen;
    filename+=_filesuffix;
    cfilename+=numbertoopen;
    cfilename+=_clusterfilesuffix;

    if (_debug <= 2)
    {
      cout << "Opening file: " << filename << endl;
    }

    // enter the file...
    TFile *f1 = TFile::Open(filename.c_str());
    f1->cd();

    TFile *f2 = TFile::Open(cfilename.c_str());
    f2->cd();

    if (_debug <= 2)
    {
      cout << "File open!" << endl;
    }

    // clone 3d hitmap
    TH3D* DUThitmap3D = (TH3D*) (f1->Get("Ntuple/DUTHitmap"))->Clone();

    // fail:
    //get the final alignment from this
    float min3dx = 0.0;
    float max3dx = 0.0;
    float min3dy = 0.0;
    float max3dy = 0.0;
    float min3dz = 0.0;
    float max3dz = 0.0;
/*
    // the hitmap is filled xzy!
    min3dx = DUThitmap3D->FindFirstBinAbove(0,1);
    max3dx = DUThitmap3D->FindLastBinAbove(0,1);
    min3dy = DUThitmap3D->FindFirstBinAbove(0,3);
    max3dy = DUThitmap3D->FindLastBinAbove(0,3);
    min3dz = DUThitmap3D->FindFirstBinAbove(0,2);
    max3dz = DUThitmap3D->FindLastBinAbove(0,2);
    
    cout << " " << endl;
    cout << "x: " << max3dx << " to " << min3dx << endl;
    cout << "yshift: " << max3dy << " to " << min3dy << endl;
    cout << "zshift: " << max3dz << " to " << min3dz << endl;

    cout << " " << endl;
    */

    // since alpha is turned, subtract the nominal position to get the "alignment"
    xshift->Fill(_sensoralignx.at(irun));
    yshift->Fill(_sensoraligny.at(irun));
    zshift->Fill(_sensoralignz.at(irun));
    ashift->Fill(_sensoraligna.at(irun) - _dutrotation.at(irun));
    bshift->Fill(_sensoralignb.at(irun));
    cshift->Fill(_sensoralignc.at(irun));
    

    // get the cluster info
    TH1D* clustersizehisto = (TH1D*) (f2->Get("MyAlibavaClustering/ClusterSize"))->Clone();
    TH1D* clusteretadistribution = (TH1D*) (f2->Get("MyAlibavaClustering/EtaDistribution"))->Clone();
    TH1D* clustersignalhisto = (TH1D*) (f2->Get("MyAlibavaClustering/SignalfromClusters"))->Clone();

    // set the run vars we write out later on...

    double clustersize = 0.0;
    double clustersizeerror = 0.0;
    double etamean = 0.0;
    double etameanerror = 0.0;
    double clustersignal = 0.0;
    double clustersignalerror = 0.0;
    double clustercount = 0.0;
    double clustercounterror = 0.0;
    double signaltonoise = 0.0;
    double signaltonoiseerror = 0.0;

    clustersize = clustersizehisto->GetMean();
    clustersizeerror = clustersizehisto->GetRMS();
    etamean = clusteretadistribution->GetMean();
    etameanerror = clusteretadistribution->GetRMS();
    clustersignal = clustersignalhisto->GetMean();
    clustersignalerror = clustersignalhisto->GetRMS();

    // the temperature correction
    float tempcorscale= 1.0;
    if (_cut_applytempcorrection == true)
    {
      if (_temperature.at(irun) == 20)
      {
	tempcorscale = 1.19110;
	cout << "Applying temperature correction!" << endl;
      }
    }



    // declare all the vars
    // general
    int Event, RunNr, EvtNr, Ndf;
    float Chi2;

    // measurements and fits
    double measX_0, measY_0, measZ_0, measQ_0, fitX_0, fitY_0;
    double measX_1, measY_1, measZ_1, measQ_1, fitX_1, fitY_1;
    double measX_2, measY_2, measZ_2, measQ_2, fitX_2, fitY_2;
    double measX_3, measY_3, measZ_3, measQ_3, fitX_3, fitY_3;
    double measX_4, measY_4, measZ_4, measQ_4, fitX_4, fitY_4;
    double measX_5, measY_5, measZ_5, measQ_5, fitX_5, fitY_5;
    double measX_6, measY_6, measZ_6, measQ_6, fitX_6, fitY_6;

    // the track point on the dut
    double dutTrackX_global, dutTrackY_global;
    double dutTrackX_local, dutTrackY_local;
    double dutTrackX_pixel, dutTrackY_pixel;

    // the hit on the dut if it was matched
    double dutHitX_global, dutHitY_global;
    double dutHitX_local, dutHitY_local;
    double dutHitX_pixel, dutHitY_pixel;

    // misc hit info
    double dutHitR, dutHitQ;

    // alibava header stuff
    float alibava_tdc, alibava_temp;

    // and the alibava reco data
    double alibava_reco_ch_[128];

    // go to position...
    TTree *ttel = (TTree*)f1->Get("Ntuple/EUFit");

    // ... and read!
    ttel->SetBranchAddress("Event", &Event);
    ttel->SetBranchAddress("RunNr", &RunNr);
    ttel->SetBranchAddress("EvtNr", &EvtNr);
    ttel->SetBranchAddress("Ndf", &Ndf);
    ttel->SetBranchAddress("Chi2", &Chi2);
    ttel->SetBranchAddress("measX_0", &measX_0);
    ttel->SetBranchAddress("measY_0", &measY_0);
    ttel->SetBranchAddress("measZ_0", &measZ_0);
    ttel->SetBranchAddress("measQ_0", &measQ_0);
    ttel->SetBranchAddress("fitX_0", &fitX_0);
    ttel->SetBranchAddress("fitY_0", &fitY_0);
    ttel->SetBranchAddress("measX_1", &measX_1);
    ttel->SetBranchAddress("measY_1", &measY_1);
    ttel->SetBranchAddress("measZ_1", &measZ_1);
    ttel->SetBranchAddress("measQ_1", &measQ_1);
    ttel->SetBranchAddress("fitX_1", &fitX_1);
    ttel->SetBranchAddress("fitY_1", &fitY_1);
    ttel->SetBranchAddress("measX_2", &measX_2);
    ttel->SetBranchAddress("measY_2", &measY_2);
    ttel->SetBranchAddress("measZ_2", &measZ_2);
    ttel->SetBranchAddress("measQ_2", &measQ_2);
    ttel->SetBranchAddress("fitX_2", &fitX_2);
    ttel->SetBranchAddress("fitY_2", &fitY_2);
    ttel->SetBranchAddress("measX_3", &measX_3);
    ttel->SetBranchAddress("measY_3", &measY_3);
    ttel->SetBranchAddress("measZ_3", &measZ_3);
    ttel->SetBranchAddress("measQ_3", &measQ_3);
    ttel->SetBranchAddress("fitX_3", &fitX_3);
    ttel->SetBranchAddress("fitY_3", &fitY_3);
    ttel->SetBranchAddress("measX_4", &measX_4);
    ttel->SetBranchAddress("measY_4", &measY_4);
    ttel->SetBranchAddress("measZ_4", &measZ_4);
    ttel->SetBranchAddress("measQ_4", &measQ_4);
    ttel->SetBranchAddress("fitX_4", &fitX_4);
    ttel->SetBranchAddress("fitY_4", &fitY_4);
    ttel->SetBranchAddress("measX_5", &measX_5);
    ttel->SetBranchAddress("measY_5", &measY_5);
    ttel->SetBranchAddress("measZ_5", &measZ_5);
    ttel->SetBranchAddress("measQ_5", &measQ_5);
    ttel->SetBranchAddress("fitX_5", &fitX_5);
    ttel->SetBranchAddress("fitY_5", &fitY_5);
    ttel->SetBranchAddress("measX_6", &measX_6);
    ttel->SetBranchAddress("measY_6", &measY_6);
    ttel->SetBranchAddress("measZ_6", &measZ_6);
    ttel->SetBranchAddress("measQ_6", &measQ_6);
    ttel->SetBranchAddress("fitX_6", &fitX_6);
    ttel->SetBranchAddress("fitY_6", &fitY_6);
    ttel->SetBranchAddress("dutTrackX_global", &dutTrackX_global);
    ttel->SetBranchAddress("dutTrackY_global", &dutTrackY_global);
    ttel->SetBranchAddress("dutTrackX_local",  &dutTrackX_local);
    ttel->SetBranchAddress("dutTrackY_local",  &dutTrackY_local);
    ttel->SetBranchAddress("dutTrackX_pixel",  &dutTrackX_pixel);
    ttel->SetBranchAddress("dutTrackY_pixel",  &dutTrackY_pixel);
    ttel->SetBranchAddress("dutHitX_global",   &dutHitX_global);
    ttel->SetBranchAddress("dutHitY_global",   &dutHitY_global);
    ttel->SetBranchAddress("dutHitX_local",    &dutHitX_local);
    ttel->SetBranchAddress("dutHitY_local",    &dutHitY_local);
    ttel->SetBranchAddress("dutHitX_pixel",    &dutHitX_pixel);
    ttel->SetBranchAddress("dutHitY_pixel",    &dutHitY_pixel);
    ttel->SetBranchAddress("dutHitR",          &dutHitR);
    ttel->SetBranchAddress("dutHitQ",          &dutHitQ);
    ttel->SetBranchAddress("alibava_tdc",      &alibava_tdc);
    ttel->SetBranchAddress("alibava_temp",     &alibava_temp);

    for(int i = 0; i < 128; ++i)
    {
      sprintf(name, "alibava_reco_ch_%i", i);
      ttel->SetBranchAddress(name, &alibava_reco_ch_[i]);
    }

    // set all to inactive, only the ones we need to 1:
    ttel->SetBranchStatus("*", 0);
    ttel->SetBranchStatus("EvtNr", 1);
    ttel->SetBranchStatus("alibava*", 1);
    ttel->SetBranchStatus("dutTrackX_*", 1);
    ttel->SetBranchStatus("dutTrackY_*", 1);
    ttel->SetBranchStatus("dutHitX_global", 1);
    ttel->SetBranchStatus("dutHitY_global", 1);
    ttel->SetBranchStatus("dutHitY_pixel", 1);
    ttel->SetBranchStatus("dutHitQ", 1);
    ttel->SetBranchStatus("measX_*", 1);
    ttel->SetBranchStatus("measX_*", 1);
    ttel->SetBranchStatus("fitX_*", 1);
    ttel->SetBranchStatus("fitX_*", 1);


    // how many entries in this tuple? -> tracks!
    long int tupleentrycount = ttel->GetEntries();

    // define the loop limit
    int tuplelimit = 0;
    if (tupleentrycount >= _maxtotaltracks)
    {
      tuplelimit = _maxtotaltracks; 
    }
    if (_maxtotaltracks > tupleentrycount)
    {
      tuplelimit = tupleentrycount;
    }

    // the run we are in
    int j = _runnumber.at(irun);

    // the histos to be filled

    TH1D* noise[128];
    for (int ii=0;ii<128;ii++)
    {
      histoname = "noise_";
      sstream << histoname << j << "_chan_" << ii;
      noise[ii] = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
      sstream.str(string());
    }

    histoname = "allnoise_";
    sstream << histoname << j ;
    TH1D* allnoise = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "adcnoise_";
    sstream << histoname << j ;
    TH1D* adcnoise = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "fivechannoise_";
    sstream << histoname << j ;
    TH1D* fivechannoise = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());


    histoname = "residualsX_";
    sstream << histoname << j ;
    TH1D* residualsX = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "residualsY_";
    sstream << histoname << j ;
    TH1D* residualsY = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "residualsXY_";
    sstream << histoname << j ;
    TH2D* residualsXY = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "residualsY_Q_";
    sstream << histoname << j ;
    TH2D* residualsY_Q = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "residualsY_clu1_";
    sstream << histoname << j ;
    TH1D* residualsY_clu1 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "residualsY_clu2_";
    sstream << histoname << j ;
    TH1D* residualsY_clu2 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "residualsY_clu3_";
    sstream << histoname << j ;
    TH1D* residualsY_clu3 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "residualsY_clu4_";
    sstream << histoname << j ;
    TH1D* residualsY_clu4 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());


    histoname = "residualsDXvsX_";
    sstream << histoname << j ;
    TH2D* residualsDXvsX = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "residualsDXvsY_";
    sstream << histoname << j ;
    TH2D* residualsDXvsY = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "residualsDYvsX_";
    sstream << histoname << j ;
    TH2D* residualsDYvsX = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "residualsDYvsY_";
    sstream << histoname << j ;
    TH2D* residualsDYvsY = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());


    histoname = "residualProfileDXvsX_";
    sstream << histoname << j ;
    TProfile* residualProfileDXvsX = dynamic_cast<TProfile*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "residualProfileDXvsY_";
    sstream << histoname << j ;
    TProfile* residualProfileDXvsY = dynamic_cast<TProfile*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "residualProfileDYvsX_";
    sstream << histoname << j ;
    TProfile* residualProfileDYvsX = dynamic_cast<TProfile*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "residualProfileDYvsY_";
    sstream << histoname << j ;
    TProfile* residualProfileDYvsY = dynamic_cast<TProfile*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());


    histoname = "residualsXvsEvt_";
    sstream << histoname << j ;
    TH2D* residualsXvsEvt = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "residualsYvsEvt_";
    sstream << histoname << j ;
    TH2D* residualsYvsEvt = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());


    histoname = "residualMapTemp_";
    sstream << histoname << j ;
    TH3D* residualMapTemp = dynamic_cast<TH3D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "residualMapY_";
    sstream << histoname << j ;
    TH2D* residualMapY = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());


    histoname = "hitmapDUT_";
    sstream << histoname << j ;
    TH2D* hitmapDUT = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "hitmapDUTmodpitch_";
    sstream << histoname << j ;
    TH2D* hitmapDUTmodpitch = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "hitmapTracks_";
    sstream << histoname << j ;
    TH2D* hitmapTracks = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "hitmapTracksmodpitch_";
    sstream << histoname << j ;
    TH2D* hitmapTracksmodpitch = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "hitmapMatch_";
    sstream << histoname << j ;
    TH2D* hitmapMatch = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "hitmapMatchmodpitch_";
    sstream << histoname << j ;
    TH2D* hitmapMatchmodpitch = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "hitmapMatchTracks_";
    sstream << histoname << j ;
    TH2D* hitmapMatchTracks = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "hitmapChargeShared0_";
    sstream << histoname << j ;
    TH2D* hitmapChargeShared0 = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "hitmapChargeShared1_";
    sstream << histoname << j ;
    TH2D* hitmapChargeShared1 = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "hitmapChargeShared2_";
    sstream << histoname << j ;
    TH2D* hitmapChargeShared2 = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "hitmapChargeShared3_";
    sstream << histoname << j ;
    TH2D* hitmapChargeShared3 = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "hitmapChargeShared4_";
    sstream << histoname << j ;
    TH2D* hitmapChargeShared4 = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "hitmapClusterSizeA_";
    sstream << histoname << j ;
    TH1D* hitmapClusterSizeA = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "hitmapClusterSize1_";
    sstream << histoname << j ;
    TH1D* hitmapClusterSize1 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "hitmapClusterSize2_";
    sstream << histoname << j ;
    TH1D* hitmapClusterSize2 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "hitmapClusterSize3_";
    sstream << histoname << j ;
    TH1D* hitmapClusterSize3 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "hitmapClusterSize4_";
    sstream << histoname << j ;
    TH1D* hitmapClusterSize4 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "signalClusterSize1_";
    sstream << histoname << j ;
    TH1D* signalClusterSize1 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "signalClusterSize2_";
    sstream << histoname << j ;
    TH1D* signalClusterSize2 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "signalClusterSize3_";
    sstream << histoname << j ;
    TH1D* signalClusterSize3 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "signalClusterSize4_";
    sstream << histoname << j ;
    TH1D* signalClusterSize4 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "divergencehistoX_";
    sstream << histoname << j ;
    TH1D* divergencehistoX = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "divergencehistoY_";
    sstream << histoname << j ;
    TH1D* divergencehistoY = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());


    histoname = "matchedEta_";
    sstream << histoname << j ;
    TH1D* matchedEta = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "matchedEta2D_";
    sstream << histoname << j ;
    TH2D* matchedEta2D = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "matchedEtaIntegral_";
    sstream << histoname << j ;
    TH1D* matchedEtaIntegral = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "unmatchedEta_";
    sstream << histoname << j ;
    TH1D* unmatchedEta = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "unmatchedEta2D_";
    sstream << histoname << j ;
    TH2D* unmatchedEta2D = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "unmatchedEtaIntegral_";
    sstream << histoname << j ;
    TH1D* unmatchedEtaIntegral = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "striphitEta_";
    sstream << histoname << j ;
    TH1D* striphitEta = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "striphitEta2D_";
    sstream << histoname << j ;
    TH2D* striphitEta2D = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "striphitEtaIntegral_";
    sstream << histoname << j ;
    TH1D* striphitEtaIntegral = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "noiseEta_";
    sstream << histoname << j ;
    TH1D* noiseEta = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "noisesubtractedEta_";
    sstream << histoname << j ;
    TH1D* noisesubtractedEta = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "EtaR10_";
    sstream << histoname << j ;
    TH1D* EtaR10 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "EtaR20_";
    sstream << histoname << j ;
    TH1D* EtaR20 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "EtaR30_";
    sstream << histoname << j ;
    TH1D* EtaR30 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());
    
    histoname = "EtaR40_";
    sstream << histoname << j ;
    TH1D* EtaR40 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "EtaR50_";
    sstream << histoname << j ;
    TH1D* EtaR50 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "EtaL10_";
    sstream << histoname << j ;
    TH1D* EtaL10 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "EtaL20_";
    sstream << histoname << j ;
    TH1D* EtaL20 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());


    histoname = "trackSignal_";
    sstream << histoname << j ;
    TH1D* trackSignal = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "trackSignalMap_";
    sstream << histoname << j ;
    TH2D* trackSignalMap = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "trackSignalMapModPitch_";
    sstream << histoname << j ;
    TH2D* trackSignalMapModPitch = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "trackSignalTDC_";
    sstream << histoname << j ;
    TH2D* trackSignalTDC = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "tracksperevent_";
    sstream << histoname << j ;
    TH1D* tracksperevent = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "tracksvsevents_";
    sstream << histoname << j ;
    TH2D* tracksvsevents = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "highchanneldistribution_";
    sstream << histoname << j ;
    TH1D* highchanneldistribution = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());


    histoname = "selectedtrackSignal_";
    sstream << histoname << j ;
    TH1D* selectedtrackSignal = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "selectedtrackSignalMap_";
    sstream << histoname << j ;
    TH2D* selectedtrackSignalMap = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "selectedtrackSignalTDC_";
    sstream << histoname << j ;
    TH2D* selectedtrackSignalTDC = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());


    histoname = "signalLeft2_";
    sstream << histoname << j ;
    TH1D* signalLeft2 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "signalLeft1_";
    sstream << histoname << j ;
    TH1D* signalLeft1 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "signalCenter_";
    sstream << histoname << j ;
    TH1D* signalCenter = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "signalRight1_";
    sstream << histoname << j ;
    TH1D* signalRight1 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "signalRight2_";
    sstream << histoname << j ;
    TH1D* signalRight2 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "signalGoodEvents_";
    sstream << histoname << j ;
    TH1D* signalGoodEvents = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "signalmapA_";
    sstream << histoname << j ;
    TH1D* signalmapA = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "signalmapB_";
    sstream << histoname << j ;
    TH1D* signalmapB = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "signalmapC_";
    sstream << histoname << j ;
    TH1D* signalmapC = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "signalmapD_";
    sstream << histoname << j ;
    TH1D* signalmapD = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "signalareaplot_";
    sstream << histoname << j;
    TGraphErrors* signalareaplot = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());


    histoname = "fiducial_discard_";
    sstream << histoname << j;
    TH1D* fiducial_discard = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "goodchannel_discard_";
    sstream << histoname << j;
    TH1D* goodchannel_discard = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "goodevent_discard_";
    sstream << histoname << j;
    TH1D* goodevent_discard = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "trackselection_discard_";
    sstream << histoname << j;
    TH1D* trackselection_discard = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "timecut_discard_";
    sstream << histoname << j;
    TH1D* timecut_discard = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "highchannel_discard_";
    sstream << histoname << j;
    TH1D* highchannel_discard = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());


    histoname = "fiducial_allow_";
    sstream << histoname << j;
    TH1D* fiducial_allow = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "goodchannel_allow_";
    sstream << histoname << j;
    TH1D* goodchannel_allow = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "goodevent_allow_";
    sstream << histoname << j;
    TH1D* goodevent_allow = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "trackselection_allow_";
    sstream << histoname << j;
    TH1D* trackselection_allow = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "timecut_allow_";
    sstream << histoname << j;
    TH1D* timecut_allow = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "highchannel_allow_";
    sstream << histoname << j;
    TH1D* highchannel_allow = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());


    // group some runs together!
    // find the right histo!
    std::vector<string>::iterator its;
    its = find(_sensorname_list.begin(), _sensorname_list.end(), _sensorname.at(irun));
    std::string::size_type pos1 = its - _sensorname_list.begin();

    std::vector<float>::iterator itf;
    itf = find(_irradfluence_list.begin(), _irradfluence_list.end(), _irradfluence.at(irun));
    size_t pos2 = itf - _irradfluence_list.begin();

    std::vector<int>::iterator it;
    it = find(_dutrotation_list.begin(), _dutrotation_list.end(), _dutrotation.at(irun));
    size_t pos3 = it - _dutrotation_list.begin();


    // if plotting quantity vs voltage the encoding is: sensor - fluence - rotation
    int vsvoltage_histonr = pos1*_irradfluence_total*_dutrotation_total + pos2*_dutrotation_total + pos3;

    // posN must be inc'ed for sensible output, as it starts at 0
    pos1++;
    pos2++;
    pos3++;

    if (_debug<=0)
    {
      cout << " " << endl;
      cout << "Sensor:   " << pos1 << " of " << _sensorname_total << endl;
      cout << "Fluence:  " << pos2 << " of " << _irradfluence_total << endl;
      cout << "Rotation: " << pos3 << " of " << _dutrotation_total << endl;
    }

    voltagegraphs_total.push_back(vsvoltage_histonr);

    it = find (voltagegraphs.begin(), voltagegraphs.end(), vsvoltage_histonr);
    if (it != voltagegraphs.end())
    {
      if (_debug <=0)
      {
	cout << "Voltage graph " << vsvoltage_histonr << " already in list!" << endl;
      }
    } else {
      voltagegraphs.push_back(vsvoltage_histonr);
      if (_debug <=0)
      {
	cout << "Voltage graph " << vsvoltage_histonr << " not in list!" << endl;
      }
    }

    histoname = "residualsXvoltage_";
    sstream << histoname << vsvoltage_histonr;
    TGraphErrors* residualsXvoltage = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "residualsYvoltage_";
    sstream << histoname << vsvoltage_histonr;
    TGraphErrors* residualsYvoltage = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());


    histoname = "signalvoltage_";
    sstream << histoname << vsvoltage_histonr;
    TGraphErrors* signalvoltage = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "snvoltage_";
    sstream << histoname << vsvoltage_histonr;
    TGraphErrors* snvoltage = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());


    histoname = "noisevoltage_";
    sstream << histoname << vsvoltage_histonr;
    TGraphErrors* noisevoltage = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "rghvoltage_";
    sstream << histoname << vsvoltage_histonr;
    TGraphErrors* rghvoltage = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());


    histoname = "clustersizevoltage_";
    sstream << histoname << vsvoltage_histonr;
    TGraphErrors* clustersizevoltage = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "clustercountvoltage_";
    sstream << histoname << vsvoltage_histonr;
    TGraphErrors* clustercountvoltage = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "channelcountvoltage_";
    sstream << histoname << vsvoltage_histonr;
    TGraphErrors* channelcountvoltage = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());


    histoname = "currentvoltage_";
    sstream << histoname << vsvoltage_histonr;
    TGraphErrors* currentvoltage = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "volumecurrent_";
    sstream << histoname << vsvoltage_histonr;
    TGraphErrors* volumecurrent = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());


    histoname = "chargesharingvoltage_";
    sstream << histoname << vsvoltage_histonr;
    TGraphErrors* chargesharingvoltage= dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "etachargesharingvoltage_";
    sstream << histoname << vsvoltage_histonr;
    TGraphErrors* etachargesharingvoltage= dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());


    histoname = "signaldistancevoltage_";
    sstream << histoname << vsvoltage_histonr;
    TGraphErrors* signaldistancevoltage= dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());


    histoname = "signalareamap_";
    sstream << histoname << vsvoltage_histonr;
    TH2D* signalareamap = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    histoname = "etamap_";
    sstream << histoname << vsvoltage_histonr;
    TH1D* etamap = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());



    // var for the event nr so we can cut
    int currentevent = -1;
    int tracksinthisevent = 1;
    bool isthisfirstevent = false;
    bool goodevent = false;
//    int goodtracksinthisevent = -1;

    TStopwatch runtime;

    if (_debug <= 2)
    {
      cout << " " << endl;
      cout << "Starting tuple loop..." << endl;
    }
    runtime.Start();

    float goodtrackcount = 0.0;
    float matchedhitcount = 0.0;

    // loop over the tuple contents and fill
    for(int i = 0; i< tuplelimit;i++)
    {

      ttel->GetEntry(i);

      // start filling histos
      
      float eventtemp = 0.0;
      eventtemp = alibava_temp;
      
      if (_temperature.at(irun) < 0)
      {
	float targetgain = 0.004700854;
	tempcorscale = targetgain / (0.00429259 + 7.9121e-5 * eventtemp - 3.82946e-6 * eventtemp*eventtemp);
	
      }


      // first find out where we are
      if (EvtNr == currentevent)
      {
	// this is not the first track in the event
	tracksinthisevent++;
	
	// this is now false
	isthisfirstevent = false;

      } else {
	// this IS the first track in the event
	// tracksinthisevent still has the track count from the previous event, so we fill the histo now!
	tracksperevent->Fill(tracksinthisevent);
	tracksvsevents->Fill(EvtNr-1,tracksinthisevent);


	//goodtracksinthisevent = 0;
	
	// did all tracks in the previous event fullfill the cuts? then we have a good event!
	if (goodevent == true)
	{
	  for (int j=0;j<tracksinthisevent;j++)
	  {
	    _goodtracks.push_back(i-1-j);
	  }

	  goodevent_allow->Fill(1.0);

	} else {
	  goodevent_discard->Fill(1.0);
	}

	// first event: reset
	goodevent = true;

	// this is now true
	isthisfirstevent = true;

	// reset
	tracksinthisevent = 1;
	currentevent = EvtNr;

      }
      // adding 0.5 makes the rounding correct!
      // this will hold the integer channel nearest to the track

      // we also have to correct for rotations
      // define the pixel as the center of the sensor, not the impact
      float angularoffset = 0.0;//_thickness.at(irun) / 2.0 * tan( _dutrotation.at(irun));

      float tempfloat = dutTrackY_pixel + 0.5 + angularoffset;
      int pixelpos = static_cast <int>(tempfloat);

      // pixelbin: hit is between .0 and .999 -> this is used for eta left right
      // this will hold the integer channel before the decimal point
      int pixelbin = static_cast <int>(dutTrackY_pixel + angularoffset);

      // this will be the signal in this track
      float tempsignal = 0.0;

      // the number of channels we have summed for 5chan noise_
      int channelssummed = 0;
      float sumnoise = 0.0;

      // noise calculation: only take the first track per event, exclude area under track + 1
      // noise eta: additionally require the channel and its neighbour to be good
 
      // noise does not get temp corrected!
      if (isthisfirstevent == true)
      {
	for (int ii=0;ii<128;ii++)
	{
	  if (ii != (pixelpos-3) && ii != (pixelpos-2) && ii != (pixelpos-1) && ii != (pixelpos) && ii != (pixelpos+1) && ii != (pixelpos+2) + ii != (pixelpos+3))
	  {
	    noise[ii]->Fill(alibava_reco_ch_[ii]);

	    if (ii > 1 && ii < 127 && _channels[irun][ii] == 1 && _channels[irun][ii+1] == 1)
	    {
	      adcnoise->Fill(alibava_reco_ch_[ii]);
	      noiseEta->Fill(alibava_reco_ch_[ii+1]/(alibava_reco_ch_[ii]+alibava_reco_ch_[ii+1]));
	    }
	    if (channelssummed <5)
	    {
	      if (_channels[irun][ii] == 1)
	      {
		sumnoise+=alibava_reco_ch_[ii];
		channelssummed++;
	      }

	    }
	    if (channelssummed == 5)
	    {
	      fivechannoise->Fill(sumnoise);
	      sumnoise = 0.0;
	      channelssummed = 0;
	    }

	  }
	}



      }

      double divergenceX = _missingvalue;
      double divergenceY = _missingvalue;

      // calculate the divergence of the triplet track...
      if (dutHitX_global != _missingvalue)
      {
	divergenceX = dutHitX_global - dutTrackX_global;
	divergenceY = dutHitY_global - dutTrackY_global;
	divergencehistoX->Fill(divergenceX);
	divergencehistoY->Fill(divergenceY);
      }

      // cut: tracks per event
      if (tracksinthisevent <= _cut_maxtracksperevent)
      {

	// fiducial x on the track
	if (dutTrackX_global >= _sensorminx.at(irun) && dutTrackX_global <= _sensormaxx.at(irun))
	{

	  // cut: existing and good channel pointed to:
	  if (pixelpos >= 2 && pixelpos <= 125 && _channels[irun][pixelpos] == 1 && _channels[irun][pixelpos+1] == 1 && _channels[irun][pixelpos-1] == 1 && _channels[irun][pixelpos+2] == 1 && _channels[irun][pixelpos-2] == 1)
	  {

	    // dutHits are already preselected in time, if this is wide, then it will not affect the hits...
	    // cut: in time
//	    if (alibava_tdc > _sensormintdc.at(irun) && alibava_tdc < _sensormaxtdc.at(irun))
//	    {

	      // cut: matched hit
	      // hit is from EUTel definition!!!
	      // cut: not missingvalue!
	      // cut on x divergence
	      //  && (dutHitX_global) > _cut_minx && (dutHitX_global) < _cut_maxx
	      if (dutHitX_global != _missingvalue && dutHitY_global != _missingvalue && dutTrackX_global != _missingvalue && dutTrackY_global != _missingvalue && fabs(divergenceX) < _cut_divergencex)
	      {
		residualsX->Fill(dutHitX_global-dutTrackX_global);
		residualsY->Fill(dutHitY_global-dutTrackY_global);
		residualsXY->Fill(dutHitX_global-dutTrackX_global,dutHitY_global-dutTrackY_global);

		residualsDXvsX->Fill(dutHitX_global-dutTrackX_global,dutTrackX_global);
		residualsDXvsY->Fill(dutHitX_global-dutTrackX_global,dutTrackY_global);
		residualsDYvsX->Fill(dutHitY_global-dutTrackY_global,dutTrackX_global);
		residualsDYvsY->Fill(dutHitY_global-dutTrackY_global,dutTrackY_global);
		residualsXvsEvt->Fill(EvtNr,dutHitX_global-dutTrackX_global);
		residualsYvsEvt->Fill(EvtNr,dutHitY_global-dutTrackY_global);

		residualProfileDXvsX->Fill(dutTrackX_global,dutHitX_global-dutTrackX_global,1);
		residualProfileDXvsY->Fill(dutTrackY_global,dutHitX_global-dutTrackX_global,1);
		residualProfileDYvsX->Fill(dutTrackX_global,dutHitY_global-dutTrackY_global,1);
		residualProfileDYvsY->Fill(dutTrackY_global,dutHitY_global-dutTrackY_global,1);

		residualMapTemp->Fill(dutTrackX_global,dutTrackY_global,dutHitY_global-dutTrackY_global);

		hitmapMatch->Fill(dutHitX_global,dutHitY_global);
		hitmapMatchTracks->Fill(dutTrackX_global,dutTrackY_global);

		fractpartx = modf (dutHitX_global/_pitchx , &intpartx);
		fractparty = modf (dutHitY_global/_pitchy , &intparty);
		hitmapMatchmodpitch->Fill(fabs(fractpartx),fabs(fractparty));

		matchedEta->Fill(tempcorscale*alibava_reco_ch_[pixelbin+1]/(tempcorscale*alibava_reco_ch_[pixelbin]+tempcorscale*alibava_reco_ch_[pixelbin+1]));
		fractparty = modf ((dutTrackY_pixel + angularoffset), &intparty);
		matchedEta2D->Fill(tempcorscale*alibava_reco_ch_[pixelbin+1]/(tempcorscale*alibava_reco_ch_[pixelbin]+tempcorscale*alibava_reco_ch_[pixelbin+1]),fractparty);

		matchedhitcount++;

	      } // done EUTel hit... from now on we don't care about them...

	      // cut: track passes through a strip
	      fractparty = modf ((dutTrackY_pixel + angularoffset), &intparty);
	      if (fractparty <= _cut_onstrip)
	      {
		striphitEta->Fill(tempcorscale*alibava_reco_ch_[pixelbin+1]/(tempcorscale*alibava_reco_ch_[pixelbin-1]+tempcorscale*alibava_reco_ch_[pixelbin+1]));
		striphitEta2D->Fill(tempcorscale*alibava_reco_ch_[pixelbin+1]/(tempcorscale*alibava_reco_ch_[pixelbin-1]+tempcorscale*alibava_reco_ch_[pixelbin+1]),fractparty);
	      }
	      if (fractparty >= _cut_onstrip+0.5)
	      {
		striphitEta->Fill(tempcorscale*alibava_reco_ch_[pixelbin+2]/(tempcorscale*alibava_reco_ch_[pixelbin]+tempcorscale*alibava_reco_ch_[pixelbin+2]));
		striphitEta2D->Fill(tempcorscale*alibava_reco_ch_[pixelbin+2]/(tempcorscale*alibava_reco_ch_[pixelbin]+tempcorscale*alibava_reco_ch_[pixelbin+2]),fractparty);
	      }

	      
	      goodtrackcount++;
	      
//	    } // done time cut, this does not kill a event!

	      goodchannel_allow->Fill(1.0);

	  } else { // done existing good channel

	    goodevent = false;

	    goodchannel_discard->Fill(1.0);

	  }

	  fiducial_allow->Fill(1.0);

	} else {// done fiducial x cut

	  goodevent = false;

	  fiducial_discard->Fill(1.0);

	}

	// all tracks: no cuts at all:
	hitmapTracks->Fill(dutTrackX_global,dutTrackY_global);
	fractpartx = modf (dutTrackX_global/_pitchx , &intpartx);
	fractparty = modf (dutTrackY_global/_pitchy , &intparty);
	hitmapTracksmodpitch->Fill(fabs(fractpartx),fabs(fractparty));


	// cut on non-missingvalue: all dut hits, if matched or not
	if (dutHitX_global != _missingvalue && dutHitY_global != _missingvalue)
	{

	  hitmapDUT->Fill(dutHitX_global,dutHitY_global);
	  fractpartx = modf (dutHitX_global/_pitchx , &intpartx);
	  fractparty = modf (dutHitY_global/_pitchy , &intparty);
	  hitmapDUTmodpitch->Fill(fabs(fractpartx),fabs(fractparty));

	} // done missingvalue cut

      } else { // tracks per event cut

	goodevent = false;

      }

      //////////
      //
      //	All cuts should be completed!
      //
      //////////


    } // done tuple contents loop


    runtime.Stop();
    if (_debug <= 2)
    {
      cout << "Done tuple loop after " << runtime.RealTime() << " s!" << endl;
    }



    // define cluster count as matched hits per good track
    //clustercount = clustersizehisto->GetEntries()/goodtrackcount;
    clustercount = matchedhitcount / goodtrackcount;
    clustercounterror = _clustercounterror;
    if (_debug <=0)
    {
      cout << "Matched cluster count / good track is: " << clustercount << " !" << endl;
    }



    // set the current point and scale it
    double scaledcurrent = 0.0;
    double scaledcurrenterror = 0.0;
    double kboltz = 0.00008617;
    double reftemp = 253.0;
    double kelvintemp = 273.15 + _temperature.at(irun);

    // go from milli amp to micro amp:
    scaledcurrent = 1000*_sensorcurrent.at(irun)*_polarity.at(irun) * reftemp * reftemp / kelvintemp / kelvintemp * exp (-1.12 / 2 / kboltz * ((1/reftemp)-(1/kelvintemp)));

    // 10% as error:
    scaledcurrenterror = scaledcurrent * 0.1;

    tempcount = currentvoltage->GetN();

    currentvoltage->SetPoint(tempcount,_biasvoltage.at(irun)*_polarity.at(irun),scaledcurrent);
    currentvoltage->SetPointError(tempcount,_voltageerror,scaledcurrenterror);

    tempcount = 0;
    tempcount = currentvoltage->GetN();

    double volcur = 0.0;
    double volcurerror = 0.0;
    
    //volcur = scaledcurrent / (_thickness.at(irun)/1e4 * 64*0.008);
    volcur = scaledcurrent / (2.5*0.512*_thickness.at(irun)/1e4);
    cout << "Scaled current: " << volcur << endl;
    volcurerror = scaledcurrenterror;

    tempcount = 0;
    tempcount = volumecurrent->GetN();

    volumecurrent->SetPoint(tempcount,_biasvoltage.at(irun)*_polarity.at(irun),volcur);
    volumecurrent->SetPointError(tempcount,_voltageerror,volcurerror);

    tempcount = 0;
    tempcount = clustersizevoltage->GetN();

    // set the cluster points
    clustersizevoltage->SetPoint(tempcount,_biasvoltage.at(irun)*_polarity.at(irun),clustersize);
    clustersizevoltage->SetPointError(tempcount,_voltageerror,clustersizeerror);

    tempcount = 0;
    tempcount = clustercountvoltage->GetN();

    if (clustercount < _cut_maxclustercount)
    {
      clustercountvoltage->SetPoint(tempcount,_biasvoltage.at(irun)*_polarity.at(irun),clustercount);
      clustercountvoltage->SetPointError(tempcount,_voltageerror,clustercounterror);
    }


    // calculate noise and rgh of this run

    int channelsfornoise = 0;
    double noisecount = 0.0;
    double noiseerrorcount = 0.0;
    double sensornoise = 0.0;
    double sensornoisecenter = 0.0;
    double sensornoiseerror = 0.0;
    double globalrgh = 0.0;
    double globalnonrgh = 0.0;
    double rghratio = 0.0;
    bool rghfail = false;
    TF1 * noisefit[128];
    for (int ii=0;ii<128;ii++)
    {
      noisefit[ii] = new TF1("noisefit","gaus");
      noise[ii]->Fit(noisefit[ii],"Q");
      double tempnoise = 0.0;
      double tempnoiseerror = 0.0;
      tempnoise = noisefit[ii]->GetParameter(2);
      tempnoiseerror = noisefit[ii]->GetParError(2);
      _noise[irun][ii] = tempnoise;
      allnoise->SetBinContent(ii,tempnoise);
      if (tempnoise > 1)
      {
	channelsfornoise++;
	noisecount+=tempnoise;
	noiseerrorcount+=tempnoiseerror;
      }

      // calculate RGH ratio for the good channels
      // loop over this channel's noise and count the entries outside of 5 * sigma
      // 1000 bins
      double rgh = 0.0;
      double nonrgh = 0.0;

      // if the noise is far too high, something already went wrong in EUTel reconstruction
      // assume a high noise to estimate the (probably high) rgh ratio
      if (tempnoise > _cut_maxrghnoise)
      {
	cout << "Fail in RGH! Too much noise in channel " << ii << " ! Breaking noise loop!" << endl;
	rghfail = true;
	break;
      }

      // only good channels
      // the maxnoise cut is redundant, but stays incase the above fixed high noise definition changes
      if (_channels[irun][ii] == 1 && tempnoise < _cut_maxrghnoise)
      {
	// loop the 1k bins
	for (int jj=1; jj<1000;jj++)
	{
	  // count the bin contents below center bin - 5sigma and above
	  if (jj < (500.0-5.0*tempnoise) || jj > (500.0+5.0*tempnoise))
	  {
	    rgh += noise[ii]->GetBinContent(jj);
	  } else {
	    nonrgh += noise[ii]->GetBinContent(jj);
	  }
	}

	globalrgh += rgh;
	globalnonrgh += nonrgh;

      }
    }

    TF1 * adcfit = new TF1("adcnoisefit","gaus");
    adcnoise->Fit(adcfit,"Q");

    if (_debug <= 0)
    {
      cout << "RGH events: " << globalrgh << " , non-RGH events: " << globalnonrgh << " !" << endl;
    }

    // catch div by zero
    if ((globalnonrgh + globalrgh) >0 && rghfail == false)
    {
      rghratio = globalrgh / (globalnonrgh + globalrgh);
    }
    if (globalnonrgh + globalrgh == 0)
    {
      cout << " " << endl;
      cout << "FAIL in RGH estimation! No hits at all!" << endl;
      cout << " " << endl;
      rghratio = 1.0;
    }
    if (rghfail == true)
    {
      cout << " " << endl;
      cout << "FAIL in RGH estimation! A channel had too much noise!" << endl;
      cout << " " << endl;
      rghratio = 1.0;
    }

    float rghpercent = rghratio * 100;

    tempcount = 0;
    tempcount = rghvoltage->GetN();

    // set the point for the graph
    if (rghpercent <= _cut_maxrghplot)
    {
      rghvoltage->SetPoint(tempcount,_biasvoltage.at(irun)*_polarity.at(irun),rghpercent);
      rghvoltage->SetPointError(tempcount,_voltageerror,rghpercent*0.1);
    }

    // the noise of the sensor, divided by the used channels
    sensornoise = noisecount / channelsfornoise;
    sensornoiseerror = noiseerrorcount / channelsfornoise;
    sensornoisecenter = adcfit->GetParameter(1);

    tempcount = 0;
    tempcount = noisevoltage->GetN();

    // cut for the histogram: only plot sensible values, drop others
    if (sensornoise > _cut_minnoise && sensornoise < _cut_maxnoise)
    {
      noisevoltage->SetPoint(tempcount,_biasvoltage.at(irun)*_polarity.at(irun),sensornoise);
      noisevoltage->SetPointError(tempcount,_voltageerror,sensornoiseerror);
    } else {
      cout << "Failed cut in max sensor noise! The noise is " << sensornoise << " ADCs, but the cut is " << _cut_maxnoise << " ADCs!" << endl;
    }

    // done noise and rgh of this run


    if (_debug <= 2)
    {
      cout << "Starting clustersize vs residual!" << endl;
    }
    int noseedcount = 0;
    int distanceseed = 0;

    for(int i = 0; i< tuplelimit;i++)
    {

      ttel->GetEntry(i);
      
      if (dutHitY_pixel != _missingvalue)
      {

      int tempclustersize = 0;
      bool goodsearch = true;
      int foundseed = -1;
      float tempsignal = 0;

      // this only finds the last seed
      for (int ii=0; ii<128; ii++)
      {
	if (alibava_reco_ch_[ii]*_polarity.at(irun) > _noise[irun][ii]*5 && _channels[irun][ii] == 1)
	{
	  foundseed = ii;
	}
      }

      // drop seed if too far from the hit, 
      if (fabs(foundseed-dutHitY_pixel) > 5)
      {
	foundseed = -1; 
	distanceseed++;
      }

      if (foundseed > 0)
      {
	tempclustersize = 1;
	tempsignal += alibava_reco_ch_[foundseed]*_polarity.at(irun);

	for (int ii=1; ii<3; ii++)
	{

	  if (alibava_reco_ch_[foundseed+ii]*_polarity.at(irun) > _noise[irun][foundseed+ii]*2.5 && goodsearch == true && _channels[irun][foundseed+ii] == 1)
	  {
	    tempclustersize++;
	    tempsignal += alibava_reco_ch_[foundseed+ii]*_polarity.at(irun);
	  } else {
	    goodsearch = false;
	  }
	}

	goodsearch = true;

	for (int ii=1; ii<3; ii++)
	{
	  if (alibava_reco_ch_[foundseed-ii]*_polarity.at(irun) > _noise[irun][foundseed-ii]*2.5 && goodsearch == true && _channels[irun][foundseed-ii] == 1)
	  {
	    tempclustersize++;
	    tempsignal += alibava_reco_ch_[foundseed-ii]*_polarity.at(irun);
	  } else {
	    goodsearch = false;
	  }
	}

      } else {
	noseedcount++;

      }
		
	fractparty = modf (dutTrackY_pixel, &intparty);
	if (tempclustersize>0)
	{
	  hitmapClusterSizeA->Fill(fractparty);
	}

		if (tempclustersize == 1)
		{
		  residualsY_clu1->Fill(dutHitY_global-dutTrackY_global);
		  hitmapClusterSize1->Fill(fractparty);
		  signalClusterSize1->Fill(tempsignal);
		}
		if (tempclustersize == 2)
		{
		  residualsY_clu2->Fill(dutHitY_global-dutTrackY_global);
		  hitmapClusterSize2->Fill(fractparty);
		  signalClusterSize2->Fill(tempsignal);
		}
		if (tempclustersize == 3)
		{
		  residualsY_clu3->Fill(dutHitY_global-dutTrackY_global);
		  hitmapClusterSize3->Fill(fractparty);
		  signalClusterSize3->Fill(tempsignal);
		}
		if (tempclustersize > 3)
		{
		  residualsY_clu4->Fill(dutHitY_global-dutTrackY_global);
		  hitmapClusterSize4->Fill(fractparty);
		  signalClusterSize4->Fill(tempsignal);
		}
		
		
		residualsY_Q->Fill(dutHitQ/1000.0,dutHitY_global-dutTrackY_global);
		
		
		
      }
    }
    
    TF1 *fit_clsz1 = lanGausFit(signalClusterSize1,3.0,7.0);
    TF1 *fit_clsz2 = lanGausFit(signalClusterSize2,3.0,7.0);
    

      TF1 * residualsYfit_1 = new TF1("residualsYfit_1","gaus + [3]", -0.1, 0.1);
      residualsYfit_1->SetParNames("Constant", "Mean", "Sigma", "Offset");
      residualsYfit_1->SetLineColor(kRed);
      residualsYfit_1->SetParameter(0, residualsY_clu1->GetMaximum());
      residualsYfit_1->SetParameter(1, residualsY_clu1->GetMean());
      residualsYfit_1->SetParameter(2, residualsY_clu1->GetRMS() * 0.5);
      residualsYfit_1->SetParameter(3, 0);
      residualsYfit_1->SetParLimits(0, 0, residualsY_clu1->GetEntries());
      residualsYfit_1->SetParLimits(1, residualsY_clu1->GetMean() - residualsY_clu1->GetRMS(), residualsY_clu1->GetMean() + residualsY_clu1->GetRMS());
      residualsYfit_1->SetParLimits(2, 0, 2 * residualsY_clu1->GetRMS());
      residualsYfit_1->SetParLimits(3, 0, residualsY_clu1->GetEntries());
      residualsY_clu1->Fit(residualsYfit_1,"QR");
      
      TF1 * residualsYfit_2 = new TF1("residualsYfit_2","gaus + [3]", -0.1, 0.1);
      residualsYfit_2->SetParNames("Constant", "Mean", "Sigma", "Offset");
      residualsYfit_2->SetLineColor(kRed);
      residualsYfit_2->SetParameter(0, residualsY_clu2->GetMaximum());
      residualsYfit_2->SetParameter(1, residualsY_clu2->GetMean());
      residualsYfit_2->SetParameter(2, residualsY_clu2->GetRMS() * 0.5);
      residualsYfit_2->SetParameter(3, 0);
      residualsYfit_2->SetParLimits(0, 0, residualsY_clu2->GetEntries());
      residualsYfit_2->SetParLimits(1, residualsY_clu2->GetMean() - residualsY_clu2->GetRMS(), residualsY_clu2->GetMean() + residualsY_clu2->GetRMS());
      residualsYfit_2->SetParLimits(2, 0, 2 * residualsY_clu2->GetRMS());
      residualsYfit_2->SetParLimits(3, 0, residualsY_clu2->GetEntries());
      residualsY_clu2->Fit(residualsYfit_2,"QR");
      
      TF1 * residualsYfit_3 = new TF1("residualsYfit_3","gaus + [3]", -0.1, 0.1);
      residualsYfit_3->SetParNames("Constant", "Mean", "Sigma", "Offset");
      residualsYfit_3->SetLineColor(kRed);
      residualsYfit_3->SetParameter(0, residualsY_clu3->GetMaximum());
      residualsYfit_3->SetParameter(1, residualsY_clu3->GetMean());
      residualsYfit_3->SetParameter(2, residualsY_clu3->GetRMS() * 0.5);
      residualsYfit_3->SetParameter(3, 0);
      residualsYfit_3->SetParLimits(0, 0, residualsY_clu3->GetEntries());
      residualsYfit_3->SetParLimits(1, residualsY_clu3->GetMean() - residualsY_clu3->GetRMS(), residualsY_clu3->GetMean() + residualsY_clu3->GetRMS());
      residualsYfit_3->SetParLimits(2, 0, 2 * residualsY_clu3->GetRMS());
      residualsYfit_3->SetParLimits(3, 0, residualsY_clu3->GetEntries());
      residualsY_clu3->Fit(residualsYfit_3,"QR");
    
    cout << "no seeds: " << noseedcount << endl;
    cout << "dist fail:" << distanceseed << endl;
    
    
    
    
    
    
    
    

    if (_debug <= 2)
    {
      cout << "Starting track selection loop!" << endl;
    }

    // limit to one track per event, find the one which gives more charge
    // can help with the 3GeV runs and kills some 0adc counts

    std::vector<int> selectedtracks;
    int tempeventnr = -1;
    float highsignal = -999.0;
    int trackselection = 0;
    bool singletrackevent = true;

    // loop on the good tracks
    for (int j=0;j<_goodtracks.size();j++)
    {

      float tempsignal = 0.0;
      ttel->GetEntry(_goodtracks.at(j));
      float angularoffset = 0.0;//_thickness.at(irun) / 2.0 * tan( _dutrotation.at(irun));
      float tempfloat = dutTrackY_pixel + 0.5 + angularoffset;
      int pixelpos = static_cast <int>(tempfloat);
      float eventtemp = alibava_temp;
      if (_temperature.at(irun) < 0)
      {
	float targetgain = 0.004700854;
	tempcorscale = targetgain / (0.00429259 + 7.9121e-5 *eventtemp - 3.82946e-6 *eventtemp*eventtemp);
      }

      tempsignal += tempcorscale*alibava_reco_ch_[pixelpos-2];
      tempsignal += tempcorscale*alibava_reco_ch_[pixelpos-1];
      tempsignal += tempcorscale*alibava_reco_ch_[pixelpos];
      tempsignal += tempcorscale*alibava_reco_ch_[pixelpos+1];
      tempsignal += tempcorscale*alibava_reco_ch_[pixelpos+2];

      // only works for j!= 0
      if (j>0)
      {

	// if this track is from the same (aka previous) event
	if (tempeventnr == EvtNr)
	{

	  singletrackevent = false;

	  // trackselection is a shift to select the track with highest signal
	  if (tempsignal*_polarity.at(irun) > highsignal*_polarity.at(irun))
	  {
	    trackselection = 0;
	  } else {
	    trackselection++;
	  }

	// this is a new event
	} else {

	  if (_cut_onlysingletrackevents == false)
	  {
	    // push the previous selected track into the vector
	    selectedtracks.push_back(_goodtracks.at(j-1-trackselection));
	    trackselection = 0;
	  }
	  
	  if (_cut_onlysingletrackevents == true && singletrackevent == true)
	  {
	    // push the previous selected track into the vector
	    selectedtracks.push_back(_goodtracks.at(j-1-trackselection));
	    trackselection = 0;
	  }

	  singletrackevent = true;
	}

      }

      // set these to compare the next iteration against
      tempeventnr = EvtNr;
      highsignal = tempsignal;

    } // done good track loop


    if (_debug <= 2)
    {
      cout << "Starting good track loop!" << endl;
    }

    int goodtracksintimecut = 0;

    float signalL2 = 0.0;
    float signalL1 = 0.0;
    float signalC = 0.0;
    float signalR1 = 0.0;
    float signalR2 = 0.0;

    // loop over selected good tracks
    for (int j=0;j<selectedtracks.size();j++)
    {

      ttel->GetEntry(selectedtracks.at(j));
      float angularoffset =0.0;// _thickness.at(irun) / 2.0 * tan( _dutrotation.at(irun));

      float tempfloat = dutTrackY_pixel + 0.5 + angularoffset;
      int pixelpos = static_cast <int>(tempfloat);
      float tempsignal = 0.0;
      int pixelbin = static_cast <int>(dutTrackY_pixel + angularoffset);

      float eventtemp = 0.0;
      eventtemp = alibava_temp;
      
      if (_temperature.at(irun) < 0)
      {
	float targetgain = 0.004700854;
	tempcorscale = targetgain / (0.00429259 + 7.9121e-5 *eventtemp - 3.82946e-6 *eventtemp*eventtemp);
      }

      // add the signal of the hit channels
      if (pixelpos <=125 && pixelpos >= 2)
      {
	
	if (alibava_tdc >= _sensormintdc.at(irun) && alibava_tdc <= _sensormaxtdc.at(irun))
	{

	  // cut: the strip with the highest pulse height is one of the center 3:
	  int highchannel = -999;
	  float highpulseheight = 0.0;
	  for (int jj = -9;jj<=9;jj++)
	  {
	    // only co		 cout << "no seed found in evt " << EvtNr << endl; nsider good channels for this:
	    if ( _channels[irun][pixelpos+jj] == 1 )
	    {
	      if ( ( tempcorscale*alibava_reco_ch_[pixelpos+jj] * _polarity.at(irun) ) > highpulseheight )
	      {
		highpulseheight = tempcorscale*alibava_reco_ch_[pixelpos+jj] * _polarity.at(irun);
		highchannel = jj;
	      }
	    }
	  }
	  highchanneldistribution->Fill(highchannel);

	  if (_cut_highchannelmiddle == false)
	  {
	    highchannel = 0;
	  }

	  if (highchannel == -1 || highchannel == 0 || highchannel == 1)
	  {

	    tempsignal = 0.0;
	    // add the signal of the hit channels
	    tempsignal += tempcorscale*alibava_reco_ch_[pixelpos-2];
	    tempsignal += tempcorscale*alibava_reco_ch_[pixelpos-1];
	    tempsignal += tempcorscale*alibava_reco_ch_[pixelpos];
	    tempsignal += tempcorscale*alibava_reco_ch_[pixelpos+1];
	    tempsignal += tempcorscale*alibava_reco_ch_[pixelpos+2];

	    // positivise
	    tempsignal = tempsignal*_polarity.at(irun);

	    // now we can fill the positivised signal
	    trackSignal->Fill(tempsignal);
	    trackSignalMap->Fill(dutTrackY_global,tempsignal);

	    //cout << " " << selectedtracks.at(j) << endl;

	    goodtracksintimecut++;

	    signalGoodEvents->Fill(tempsignal);

	    signalL2 += tempcorscale*alibava_reco_ch_[pixelpos-2]*_polarity.at(irun);
	    signalLeft2->Fill(tempcorscale*alibava_reco_ch_[pixelpos-2]*_polarity.at(irun));

	    signalL1 += tempcorscale*alibava_reco_ch_[pixelpos-1]*_polarity.at(irun);
	    signalLeft1->Fill(tempcorscale*alibava_reco_ch_[pixelpos-1]*_polarity.at(irun));

	    signalC += tempcorscale*alibava_reco_ch_[pixelpos]*_polarity.at(irun);
	    signalCenter->Fill(tempcorscale*alibava_reco_ch_[pixelpos]*_polarity.at(irun));

	    signalR1 += tempcorscale*alibava_reco_ch_[pixelpos+1]*_polarity.at(irun);
	    signalRight1->Fill(tempcorscale*alibava_reco_ch_[pixelpos+1]*_polarity.at(irun));

	    signalR2 += tempcorscale*alibava_reco_ch_[pixelpos+2]*_polarity.at(irun);
	    signalRight2->Fill(tempcorscale*alibava_reco_ch_[pixelpos+2]*_polarity.at(irun));

	    fractparty = modf ((dutTrackY_pixel + angularoffset), &intparty);
	    unmatchedEta->Fill(tempcorscale*alibava_reco_ch_[pixelbin+1]/(tempcorscale*alibava_reco_ch_[pixelbin]+tempcorscale*alibava_reco_ch_[pixelbin+1]));
	    unmatchedEta2D->Fill(tempcorscale*alibava_reco_ch_[pixelbin+1]/(tempcorscale*alibava_reco_ch_[pixelbin]+tempcorscale*alibava_reco_ch_[pixelbin+1]),fractparty);

	    int temppixelbin = static_cast <int>(dutTrackY_pixel + 0.5);
	    EtaR50->Fill(tempcorscale*alibava_reco_ch_[temppixelbin+1]/(tempcorscale*alibava_reco_ch_[temppixelbin]+tempcorscale*alibava_reco_ch_[temppixelbin+1]));
	    temppixelbin = static_cast <int>(dutTrackY_pixel + 0.4);
	    EtaR40->Fill(tempcorscale*alibava_reco_ch_[temppixelbin+1]/(tempcorscale*alibava_reco_ch_[temppixelbin]+tempcorscale*alibava_reco_ch_[temppixelbin+1]));
	    temppixelbin = static_cast <int>(dutTrackY_pixel + 0.3);
	    EtaR30->Fill(tempcorscale*alibava_reco_ch_[temppixelbin+1]/(tempcorscale*alibava_reco_ch_[temppixelbin]+tempcorscale*alibava_reco_ch_[temppixelbin+1]));
	    temppixelbin = static_cast <int>(dutTrackY_pixel + 0.2);
	    EtaR20->Fill(tempcorscale*alibava_reco_ch_[temppixelbin+1]/(tempcorscale*alibava_reco_ch_[temppixelbin]+tempcorscale*alibava_reco_ch_[temppixelbin+1]));
	    temppixelbin = static_cast <int>(dutTrackY_pixel + 0.1);
	    EtaR10->Fill(tempcorscale*alibava_reco_ch_[temppixelbin+1]/(tempcorscale*alibava_reco_ch_[temppixelbin]+tempcorscale*alibava_reco_ch_[temppixelbin+1]));
	    temppixelbin = static_cast <int>(dutTrackY_pixel - 0.1);
	    EtaL10->Fill(tempcorscale*alibava_reco_ch_[temppixelbin+1]/(tempcorscale*alibava_reco_ch_[temppixelbin]+tempcorscale*alibava_reco_ch_[temppixelbin+1]));
	    temppixelbin = static_cast <int>(dutTrackY_pixel - 0.2);
	    EtaL20->Fill(tempcorscale*alibava_reco_ch_[temppixelbin+1]/(tempcorscale*alibava_reco_ch_[temppixelbin]+tempcorscale*alibava_reco_ch_[temppixelbin+1]));


	    // divide the interstrip region into 8 parts and fill the signal accordingly
	    if(fractparty>=0.0 && fractparty<0.25)
	    {
	      signalmapA->Fill(tempsignal);
	      trackSignalMapModPitch->Fill(0.0,tempsignal);
	    }
	    if(fractparty>=0.25 && fractparty<0.5)
	    {
	      signalmapB->Fill(tempsignal);
	      trackSignalMapModPitch->Fill(0.25,tempsignal);
	    }
	    if(fractparty>=0.5 && fractparty<0.75)
	    {
	      signalmapC->Fill(tempsignal);
	      trackSignalMapModPitch->Fill(0.5,tempsignal);
	    }
	    if(fractparty>=0.75 && fractparty<1.0)
	    {
	      signalmapD->Fill(tempsignal);
	      trackSignalMapModPitch->Fill(0.75,tempsignal);
	    }

	    // charge sharing map:
	    // identify the shared channels:
	    bool c = false;
	    bool l1 = false;
	    bool l2 = false;
	    bool l3 = false;
	    bool r1 = false;
	    bool r2 = false;
	    bool r3 = false;

	    // if we have charge in the pointed-to chans:
	    if (tempcorscale*alibava_reco_ch_[pixelpos]*_polarity.at(irun) > _noise[irun][pixelpos]*_cut_chargesharing)
	    {
	      c = true;
	    }
	    if (tempcorscale*alibava_reco_ch_[pixelpos-1]*_polarity.at(irun) > _noise[irun][pixelpos-1]*_cut_chargesharing)
	    {
	      l1 = true;
	    }
	    if (tempcorscale*alibava_reco_ch_[pixelpos-2]*_polarity.at(irun) > _noise[irun][pixelpos-2]*_cut_chargesharing)
	    {
	      l2 = true;
	    }
	    if (tempcorscale*alibava_reco_ch_[pixelpos-3]*_polarity.at(irun) > _noise[irun][pixelpos-3]*_cut_chargesharing)
	    {
	      l3 = true;
	    }
	    if (tempcorscale*alibava_reco_ch_[pixelpos+1]*_polarity.at(irun) > _noise[irun][pixelpos+1]*_cut_chargesharing)
	    {
	      r1 = true;
	    }
	    if (tempcorscale*alibava_reco_ch_[pixelpos+2]*_polarity.at(irun) > _noise[irun][pixelpos+2]*_cut_chargesharing)
	    {
	      r2 = true;
	    }
	    if (tempcorscale*alibava_reco_ch_[pixelpos+3]*_polarity.at(irun) > _noise[irun][pixelpos+3]*_cut_chargesharing)
	    {
	      r3 = true;
	    }

	    // fill maps
	    if (c==true)
	    {
	      // 1:
	      if (l1==false && r1==false)
	      {
		hitmapChargeShared1->Fill(dutTrackX_global,dutTrackY_global);
	      }

	      // 2:
	      if (l1==false && r1==true && r2==false)
	      {
		hitmapChargeShared2->Fill(dutTrackX_global,dutTrackY_global);
	      }
	      if (l2=false && l1==true && r1==false)
	      {
		hitmapChargeShared2->Fill(dutTrackX_global,dutTrackY_global);
	      }

	      // 3:
	      if (l2==false && l1==true && r1==true && r2 == false)
	      {
		hitmapChargeShared3->Fill(dutTrackX_global,dutTrackY_global);
	      }
	      if (l3==false && l2==true && l1==true && r1 == false)
	      {
		hitmapChargeShared3->Fill(dutTrackX_global,dutTrackY_global);
	      }
	      if (l1==false && r1==true && r2==true && r3 == false)
	      {
		hitmapChargeShared3->Fill(dutTrackX_global,dutTrackY_global);
	      }

	      // 4:
	      if (l3==true && l2==true && l1==true && r1==false)
	      {
		hitmapChargeShared4->Fill(dutTrackX_global,dutTrackY_global);
	      }
	      if (l3==false && l2==true && l1==true && r1 == true && r2==false)
	      {
		hitmapChargeShared4->Fill(dutTrackX_global,dutTrackY_global);
	      }
	      if (l2==false && l1==true && r1==true && r2 == true && r3==false)
	      {
		hitmapChargeShared4->Fill(dutTrackX_global,dutTrackY_global);
	      }
	      if (l1==false && r1==true && r2==true && r3 == true)
	      {
		hitmapChargeShared4->Fill(dutTrackX_global,dutTrackY_global);
	      }


	    } else {
	      hitmapChargeShared0->Fill(dutTrackX_global,dutTrackY_global);
	    } // done charge in pointed to channel

	    highchannel_allow->Fill(1.0);

	  } else {

	    highchannel_discard->Fill(1.0);

	  }

	  timecut_allow->Fill(1.0);

	} else { // done timecut

	  timecut_discard->Fill(1.0);

	}

	// this kills runtime!

	if (alibava_tdc>0.1)
	{

	  tempsignal = 0.0;
	  tempsignal += tempcorscale*alibava_reco_ch_[pixelpos-2];
	  tempsignal += tempcorscale*alibava_reco_ch_[pixelpos-1];
	  tempsignal += tempcorscale*alibava_reco_ch_[pixelpos];
	  tempsignal += tempcorscale*alibava_reco_ch_[pixelpos+1];
	  tempsignal += tempcorscale*alibava_reco_ch_[pixelpos+2];

	  // positivise
	  tempsignal = tempsignal*_polarity.at(irun);
	  trackSignalTDC->Fill(alibava_tdc,tempsignal);

	} // 0.1 drops some entries to save time...

      } // pixelpos safety check

    }

    // the amount of charge sharing
    float chargeshared = 0.0;
    float tracksshared = 0.0;
    tracksshared += hitmapChargeShared4->GetEntries();
    tracksshared += hitmapChargeShared3->GetEntries();
    tracksshared += hitmapChargeShared2->GetEntries();

    if (_debug <1)
    {
      cout << "Tracks sharing charge: " << tracksshared << endl;
      cout << "Total tracks         :  " << goodtracksintimecut << endl;
    }

    chargeshared = tracksshared / goodtracksintimecut * 100;
    tempcount = 0;
    tempcount = chargesharingvoltage->GetN();
    chargesharingvoltage->SetPoint(tempcount,_biasvoltage.at(irun)*_polarity.at(irun),chargeshared);
    chargesharingvoltage->SetPointError(tempcount,_voltageerror,0.1*chargeshared);

    if (_debug <1)
    {
      cout << "shared charge: " << chargeshared << endl;
      float temptemp = 0.0;
      temptemp = (signalL1 + signalR1) / (signalL1 + signalC + signalR1) * 100.0;
      cout << "alternative: " << temptemp << endl;

      //chargesharingvoltage->SetPoint(tempcount,_biasvoltage.at(irun)*_polarity.at(irun),temptemp);
      //chargesharingvoltage->SetPointError(tempcount,_voltageerror,0.1*temptemp);

    }

    _goodtracks.clear();
    selectedtracks.clear();


    // integrate the eta distributions
    double counts = 0.0;
    double integral = 0.0;
    for (int ii = 1; ii<matchedEta->GetNbinsX(); ii++)
    {
      counts = matchedEta->GetBinContent(ii);
      integral += counts ;
      matchedEtaIntegral->SetBinContent(ii,integral);
    }

    counts = 0.0;
    integral = 0.0;
    for (int ii = 1; ii<striphitEta->GetNbinsX(); ii++)
    {
      counts = striphitEta->GetBinContent(ii);
      integral += counts ;
      striphitEtaIntegral->SetBinContent(ii,integral);
    }

    counts = 0.0;
    integral = 0.0;
    for (int ii = 1; ii<unmatchedEta->GetNbinsX(); ii++)
    {
      counts = unmatchedEta->GetBinContent(ii);
      integral += counts ;
      unmatchedEtaIntegral->SetBinContent(ii,integral);
    }

    // since we have calculated the eta distribution of noise, we can subtract it from unmatchedeta...

    // scale to -0.5 and +1.5, average and subtract
    float leftnoise = noiseEta->GetBinContent(20);
    float rightnoise = noiseEta->GetBinContent(100);
    float avgetanoise = (leftnoise+rightnoise)/2.0;
    float etaleftsig = unmatchedEta->GetBinContent(20);
    float etarightsig = unmatchedEta->GetBinContent(100);
    float avgetasig = (etaleftsig+etarightsig)/2.0;
    float etascale = avgetasig/avgetanoise;

    if (_debug<1)
    {
      cout << "eta scale is " << etascale << endl;
    }

    for (int ii=1;ii<unmatchedEta->GetNbinsX();ii++)
    {
      float tempeta = unmatchedEta->GetBinContent(ii);
      float tempeta2 = noiseEta->GetBinContent(ii);
      if ((tempeta - etascale*tempeta2) > 0)
      {
	noisesubtractedEta->SetBinContent(ii,(tempeta - etascale*tempeta2));
      } else {
	noisesubtractedEta->SetBinContent(ii,0);
      }
    }


    // define charge sharing out of eta:
    float totaleta = unmatchedEta->GetEntries();
    // noise are the counts <0 and >1
    float etanoise = 0.0;
    float etachargeshared = 0.0;

    // which eta distribution to use? unmatched or noisesubtractedEta ?


    for (int ii = 1; ii<unmatchedEta->GetNbinsX(); ii++)
    {
      if (unmatchedEta->GetBinCenter(ii) < 0 || unmatchedEta->GetBinCenter(ii) > 1.0)
      {
	etanoise+=unmatchedEta->GetBinContent(ii);
      }
      if (unmatchedEta->GetBinCenter(ii) > _cut_eta_chargeshare && unmatchedEta->GetBinCenter(ii) < (1.0 - _cut_eta_chargeshare) )
      {
	etachargeshared+=unmatchedEta->GetBinContent(ii);
      }
    }

    /*

    for (int ii = 1; ii<noisesubtractedEta->GetNbinsX(); ii++)
    {
      if (noisesubtractedEta->GetBinCenter(ii) < 0 || noisesubtractedEta->GetBinCenter(ii) > 1.0)
      {
	etanoise+=noisesubtractedEta->GetBinContent(ii);
      }
      if (noisesubtractedEta->GetBinCenter(ii) > _cut_eta_chargeshare && noisesubtractedEta->GetBinCenter(ii) < (1.0 - _cut_eta_chargeshare) )
      {
	etachargeshared+=noisesubtractedEta->GetBinContent(ii);
      }
    }

    */


    //cout << "total eta counts: " << totaleta << endl;
    //cout << "shared counts: " << etachargeshared << endl;
    if ((totaleta) > 0)
    {
      //etanoise = etanoise / totaleta * 100.0;
      etachargeshared = etachargeshared / (totaleta) * 100.0;
    } else {
      etanoise = 0.0;
      etachargeshared = 0.0;
    }

    tempcount = 0;
    tempcount = etachargesharingvoltage->GetN();
    etachargesharingvoltage->SetPoint(tempcount,_biasvoltage.at(irun)*_polarity.at(irun),etachargeshared);
    etachargesharingvoltage->SetPointError(tempcount,_voltageerror,0.1*etachargeshared);

    if (_debug < 1)
    {
      cout << "Done eta charge sharing!" << endl;
    }

    // residuals of this run

    // cut for the histograms: only plot sensible values, drop others
    // assume that a correct residual needs 500 entries...
    if (residualsY->GetEntries() > 500)
    {

      if (_debug < 1)
      {
	cout << "Calculating residual!" << endl;
      }

      TF1 * residualsXfit = new TF1("residualsXfit","gaus + [3]", -0.1, 0.1);
      TF1 * residualsYfit = new TF1("residualsYfit","gaus + [3]", -0.1, 0.1);
      residualsXfit->SetParNames("Constant", "Mean", "Sigma", "Offset");
      residualsXfit->SetLineColor(kRed);
      residualsXfit->SetParameter(0, residualsX->GetMaximum());
      residualsXfit->SetParameter(1, residualsX->GetMean());
      residualsXfit->SetParameter(2, residualsX->GetRMS() * 0.5);
      residualsXfit->SetParameter(3, 0);
      residualsXfit->SetParLimits(0, 0, residualsX->GetEntries());
      residualsXfit->SetParLimits(1, residualsX->GetMean() - residualsX->GetRMS(), residualsX->GetMean() + residualsX->GetRMS());
      residualsXfit->SetParLimits(2, 0, 2 * residualsX->GetRMS());
      residualsXfit->SetParLimits(3, 0, residualsX->GetEntries());

      residualsYfit->SetParNames("Constant", "Mean", "Sigma", "Offset");
      residualsYfit->SetLineColor(kRed);
      residualsYfit->SetParameter(0, residualsY->GetMaximum());
      residualsYfit->SetParameter(1, residualsY->GetMean());
      residualsYfit->SetParameter(2, residualsY->GetRMS() * 0.5);
      residualsYfit->SetParameter(3, 0);
      residualsYfit->SetParLimits(0, 0, residualsY->GetEntries());
      residualsYfit->SetParLimits(1, residualsY->GetMean() - residualsY->GetRMS(), residualsY->GetMean() + residualsY->GetRMS());
      residualsYfit->SetParLimits(2, 0, 2 * residualsY->GetRMS());
      residualsYfit->SetParLimits(3, 0, residualsY->GetEntries());

      residualsX->Fit(residualsXfit,"QR");
      residualsY->Fit(residualsYfit,"QR");

      // set the graph points for resi vs voltage
      double fitparx1 = 0.0;
      double fitpary1 = 0.0;
      double fitparx2 = 0.0;
      double fitpary2 = 0.0;
      fitparx1 = residualsXfit->GetParameter(2);
      fitparx2 = residualsXfit->GetParError(2);
      fitpary1 = residualsYfit->GetParameter(2);
      fitpary2 = residualsYfit->GetParError(2);

      tempcount = 0;
      tempcount = residualsXvoltage->GetN();
      //cout << " " << endl;
      //cout << tempcount << " points in Residuals Y graph, adding one..." << endl;
      if (fitparx1 < _cut_maxXresidual)
      {
	residualsXvoltage->SetPoint(tempcount,_biasvoltage.at(irun)*_polarity.at(irun),fitparx1);
	residualsXvoltage->SetPointError(tempcount,_voltageerror,fitparx2);
      } else {
	if (_debug <= 2)
	{
	  cout << "Dropping residual because it is too high!" << endl;
	}
      }

      tempcount = 0;
      tempcount = residualsYvoltage->GetN();

      if (fitpary1 < _cut_maxYresidual)
      {
	residualsYvoltage->SetPoint(tempcount,_biasvoltage.at(irun)*_polarity.at(irun),fitpary1);
	residualsYvoltage->SetPointError(tempcount,_voltageerror,fitpary2);
      } else {
	if (_debug <= 2)
	{
	  cout << "Dropping residual because it is too high!" << endl;
	}
      }
    } else {
      if (_debug <= 2)
      {
	cout << "Dropping residual because of entry count!" << endl;
      }
    }

    // done this runs residuals


    // signal  and fit of this run:

    // the clone-histos can only be written later on if they exist -> the (if tracksignal ) has to be true
    bool dosubtrackhistosexist = false;
    bool skipsignals = false;

    TH1D *subtractedtrackSignal;

    // limit this to good track signals -> at least 100 entires
    if (trackSignal->GetEntries() > 100)
    {

      if (_debug < 1)
      {
	cout << "Calculating signal!" << endl;
      }

      dosubtrackhistosexist = true;
      // a gaus to get the first noise peak
      TF1 * gausnoise = new TF1("gausnoise","gaus",-150,150.0);

      // constrain the parameters:
      // first set all to "sensible" values: the constant to the bin entry at 0 ADCs (500th bin)
      // the mean gets set to the sensor noise mean and fixed
      // the sigma is assumed to be 1.5 times the channel sigma
      gausnoise->SetParameter(0,trackSignal->GetBinContent(500));
      gausnoise->SetParameter(1,sensornoisecenter);
      gausnoise->SetParameter(2,sensornoise*1.5);

      // limits on the parameters:
      // if we have five channels with only noise contributions, their width will be:
      float fivechannoisetest = sqrt(5*sensornoise*sensornoise);
      gausnoise->SetParLimits(0,10,2000);
      gausnoise->FixParameter(1,sensornoisecenter);
      gausnoise->SetParLimits(2,sensornoise,fivechannoisetest);

      // fit to 5adcs:
      trackSignal->Fit(gausnoise,"RBQ","",-50.0,5.0);
      gausnoise->Draw();

      // copy this histo, extend the fit drawing range
      TH1D *selectedtrackSignal = (TH1D*)trackSignal->Clone("selectedtrackSignal");
      if (selectedtrackSignal->GetListOfFunctions()->Contains("gausnoise") )
      {
	TF1 *totalgausnoise = selectedtrackSignal->GetFunction("gausnoise");
	totalgausnoise->SetRange(-150,150);
	totalgausnoise->Draw();

	// now copy the copy and subtract the fit
	subtractedtrackSignal = (TH1D*)selectedtrackSignal->Clone("subtractedtrackSignal");

	// now we can subtract the noise peak from the signal...
	subtractedtrackSignal->Add(totalgausnoise,-1);

	if (_debug <= 2)
	{
	  cout << " " << endl;
	  cout << "Call to Landau-Gaussian fit on noise subtracted data!" << endl;
	}

	// the call to the landau fit on the noise-cleaned signal:
	TF1 *fitsnr2 = lanGausFit(subtractedtrackSignal,1.0,7.0);

	// so that the histo is accessible outside of f1
	subtractedtrackSignal->SetDirectory(0);

      } else {
	cout << "Fail in gaus noise fit!" << endl;
	dosubtrackhistosexist = false;
      }

      // fit this guy with a landau-gaus
      //TF1 *fit_goodevents = lanGausFit(signalGoodEvents,3.0,7.0);
      TF1 *fit_goodevents = gausLanGausFitFixGausNoise(signalGoodEvents,1.0,7.0,fivechannoise->GetMean()*_polarity.at(irun),fivechannoise->GetRMS());


      // set the graph points for signal vs voltage
      double fitparx1 = 0.0;
      double fitparx2 = 0.0;
      fitparx1 = fit_goodevents->GetParameter(4);
      fitparx2 = fit_goodevents->GetParError(4);

      // cut for the histogram: only plot sensible values, drop others
      if (fitparx1 > _cut_minsignal && fitparx1 < _cut_maxsignal)
      {

	tempcount = 0;
	tempcount = signalvoltage->GetN();

	signalvoltage->SetPoint(tempcount,_biasvoltage.at(irun)*_polarity.at(irun),fitparx1);
	signalvoltage->SetPointError(tempcount,_voltageerror,fitparx2);

	tempcount = 0;
	tempcount = signaldistancevoltage->GetN();

	signaldistancevoltage->SetPoint(tempcount,_biasvoltage.at(irun)*_polarity.at(irun),fitparx1/(_thickness.at(irun)/cos(_dutrotation.at(irun)*PI/180.0)));
	signaldistancevoltage->SetPointError(tempcount,_voltageerror,fitparx2/(_thickness.at(irun)/cos(_dutrotation.at(irun)*PI/180.0)));
	
	//cout << "point is " << fitparx1/(_thickness.at(irun)/cos(_dutrotation.at(irun)*PI/180.0)) << endl;
	//cout << "distance is " << cos(_dutrotation.at(irun)*PI/180.0) << endl;

	// signal to noise:
	if (sensornoise != 0)
	{
	  signaltonoise = fitparx1 / sensornoise;
	  signaltonoiseerror = sqrt(fitparx2*fitparx2 + sensornoiseerror*sensornoiseerror);

	  tempcount = 0;
	  tempcount = snvoltage->GetN();

	  snvoltage->SetPoint(tempcount,_biasvoltage.at(irun)*_polarity.at(irun),signaltonoise);
	  snvoltage->SetPointError(tempcount,_voltageerror,signaltonoiseerror);
	} else {
	  cout << "Dropping S/N point: Noise is 0!" << endl;
	}
      } else {

	cout << "Dropping signal, as it is out of range! ( " << fitparx1 << " )" << endl;

      }

      // estimate how many channels (of the 5) are giving signal and how many are noise:
      // get the fit's noise sigma and compare it to the expected width of 5 noisy channels
      float effectivechan = 0.0;
      float effectivechanerror = 0.0;
      if (sensornoise != 0)
      {
	//effectivechan = 5.0 - (fit_goodevents->GetParameter(3))*(fit_goodevents->GetParameter(3))/sensornoise/sensornoise;
	effectivechan = (fit_goodevents->GetParameter(6))*(fit_goodevents->GetParameter(6)) - 5*sensornoise*sensornoise;
	if (_debug < 1)
	{
	  cout << "Effective channel count is " << effectivechan << endl;
	}
      } else {
	cout << "Dropping channel point! Noise is 0!" << endl;
      }
      if (sensornoiseerror != 0)
      {
	effectivechanerror = sqrt(fit_goodevents->GetParError(6)*fit_goodevents->GetParError(6) +sensornoiseerror*sensornoiseerror);
      } else {
	cout << "Dropping channel point error! Noiseerror is 0!" << endl;
      }

      tempcount = 0;
      tempcount = channelcountvoltage->GetN();

      channelcountvoltage->SetPoint(tempcount,_biasvoltage.at(irun)*_polarity.at(irun),effectivechan);
      channelcountvoltage->SetPointError(tempcount,_voltageerror,0.0);

    } else {
      if (_debug <= 2)
      {
	cout << "Dropping signal, not enough entries!" << endl;
	skipsignals = true;
      }
    }

    if (skipsignals == false)
    {

    TF1 *fit_signalmapA = gausLanGausFitFixGausNoise(signalmapA,1.0,7.0,fivechannoise->GetMean()*_polarity.at(irun),fivechannoise->GetRMS());
    TF1 *fit_signalmapB = gausLanGausFitFixGausNoise(signalmapB,1.0,7.0,fivechannoise->GetMean()*_polarity.at(irun),fivechannoise->GetRMS());
    TF1 *fit_signalmapC = gausLanGausFitFixGausNoise(signalmapC,1.0,7.0,fivechannoise->GetMean()*_polarity.at(irun),fivechannoise->GetRMS());
    TF1 *fit_signalmapD = gausLanGausFitFixGausNoise(signalmapD,1.0,7.0,fivechannoise->GetMean()*_polarity.at(irun),fivechannoise->GetRMS());

    double fitpara = 0.0;
    double fitparb = 0.0;
    double fitparc = 0.0;
    double fitpard = 0.0;
    double fitparera = 0.0;
    double fitparerb = 0.0;
    double fitparerc = 0.0;
    double fitparerd = 0.0;

    fitpara = fit_signalmapA->GetParameter(4);
    fitparb = fit_signalmapB->GetParameter(4);
    fitparc = fit_signalmapC->GetParameter(4);
    fitpard = fit_signalmapD->GetParameter(4);
    fitparera = fit_signalmapA->GetParError(4);
    fitparerb = fit_signalmapB->GetParError(4);
    fitparerc = fit_signalmapC->GetParError(4);
    fitparerd = fit_signalmapD->GetParError(4);

    int thisvoltage = _biasvoltage.at(irun)*_polarity.at(irun) / 100;

    // signalareamap is the correct histo, so we can fill it now:
    signalareamap->SetBinContent(thisvoltage,1,fitpara);
    signalareamap->SetBinContent(thisvoltage,2,fitparb);
    signalareamap->SetBinContent(thisvoltage,3,fitparc);
    signalareamap->SetBinContent(thisvoltage,4,fitpard);


    int tempcount = signalareaplot->GetN();
    signalareaplot->SetPoint(tempcount,0.125,fitpara);
    signalareaplot->SetPointError(tempcount,0.05,fitparera);
    tempcount = signalareaplot->GetN();
    signalareaplot->SetPoint(tempcount,0.375,fitparb);
    signalareaplot->SetPointError(tempcount,0.05,fitparerb);
    tempcount = signalareaplot->GetN();
    signalareaplot->SetPoint(tempcount,0.625,fitparc);
    signalareaplot->SetPointError(tempcount,0.05,fitparerc);
    tempcount = signalareaplot->GetN();
    signalareaplot->SetPoint(tempcount,0.875,fitpard);
    signalareaplot->SetPointError(tempcount,0.05,fitparerd);
    }
    if (_debug <= 2)
    {
      cout << " " << endl;
      cout << "Done signal processing, writing this run to file!" << endl;
    }

    // done with the signal processing












    // starting with the output of this run

    // done with the loaded ntuple file
    f1->Close();

    // doesn't work? FIXME
    // may need a root object/application loaded at the start... 
    gStyle->SetOptFit(1111);

    // write all the individual histos to file
    _outputFile->cd();

    // each run gets its own directory
    int therun = _runnumber.at(irun);
    sprintf(name,"run%i",therun);
    TDirectory* rundirectory = _outputFile->mkdir(name);

    rundirectory->cd();
    TDirectory* noisedirectory = rundirectory->mkdir("noise");
    noisedirectory->cd();
    TDirectory* noisechandirectory = noisedirectory->mkdir("channels");
    noisechandirectory->cd();
    for(int ii=0;ii<128;ii++)
    {
      noise[ii]->Write();
    }
    noisedirectory->cd();
    allnoise->Write();
    adcnoise->Write();
    fivechannoise->Write();

    rundirectory->cd();
    TDirectory* resdirectory = rundirectory->mkdir("residuals");
    resdirectory->cd();

    residualsX->Write();
    residualsY->Write();
    residualsXY->Write();
    residualsY_Q->Write();

    TDirectory* resdirectorysub = resdirectory->mkdir("Clustersizes");
    resdirectorysub->cd();

    residualsY_clu1->Write();
    residualsY_clu2->Write();
    residualsY_clu3->Write();
    residualsY_clu4->Write();

    resdirectory->cd();

    residualsDXvsX->Write();
    residualsDXvsY->Write();
    residualsDYvsX->Write();
    residualsDYvsY->Write();

    residualsXvsEvt->Write();
    residualsYvsEvt->Write();

    residualMapY->Write();

    residualProfileDXvsX->Write();
    residualProfileDXvsY->Write();
    residualProfileDYvsX->Write();
    residualProfileDYvsY->Write();

    rundirectory->cd();
    TDirectory* hitmapdirectory = rundirectory->mkdir("hitmaps");
    hitmapdirectory->cd();

    hitmapDUT->Write();
    hitmapDUTmodpitch->Write();
    DUThitmap3D->Write();

    hitmapTracks->Write();
    hitmapTracksmodpitch->Write();

    hitmapMatch->Write();
    hitmapMatchmodpitch->Write();
    hitmapMatchTracks->Write();

    divergencehistoX->Write();
    divergencehistoY->Write();

    hitmapChargeShared0->Write();
    hitmapChargeShared1->Write();
    hitmapChargeShared2->Write();
    hitmapChargeShared3->Write();
    hitmapChargeShared4->Write();

    hitmapClusterSizeA->Write();
    hitmapClusterSize1->Write();
    hitmapClusterSize2->Write();
    hitmapClusterSize3->Write();
    hitmapClusterSize4->Write();

    signalClusterSize1->Write();
    signalClusterSize2->Write();
    signalClusterSize3->Write();
    signalClusterSize4->Write();

    rundirectory->cd();
    TDirectory* etadirectory = rundirectory->mkdir("eta");
    etadirectory->cd();

    // a canvas for the eta integrals
    if (_cut_drawetaintegralcomparison == true)
    {

      TCanvas *canv_integral = new TCanvas("Eta Integrals","transparent pad",600,400);
      gStyle->SetOptStat(kFALSE);
      unmatchedEtaIntegral->SetStats(kFALSE);
      unmatchedEtaIntegral->Draw();
      canv_integral->Update();

      Float_t rightmax = 1.05*matchedEtaIntegral->GetMaximum();
      Float_t scale = gPad->GetUymax()/rightmax;
      matchedEtaIntegral->SetLineColor(kRed);
      matchedEtaIntegral->Scale(scale);
      matchedEtaIntegral->SetStats(kFALSE);
      matchedEtaIntegral->Draw("same");

      //draw an axis on the right side
      TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
      gPad->GetUxmax(), gPad->GetUymax(),0,rightmax,510,"+L");
      axis->SetLineColor(kRed);
      axis->SetLabelColor(kRed);
      axis->SetTitle("Events");
      axis->Draw();

      canv_integral->Write();

    }

    matchedEta->Write();
    matchedEta2D->Write();
    matchedEtaIntegral->Write();
    unmatchedEta->Write();
    unmatchedEta2D->Write();
    unmatchedEtaIntegral->Write();
    striphitEta->Write();
    striphitEta2D->Write();
    striphitEtaIntegral->Write();
    noiseEta->Write();
    noisesubtractedEta->Write();
    TDirectory* etashiftdirectory = etadirectory->mkdir("shifts");
    etashiftdirectory->cd();
    EtaL20->Write();
    EtaL10->Write();
    EtaR10->Write();
    EtaR20->Write();
    EtaR30->Write();
    EtaR40->Write();
    EtaR50->Write();

    rundirectory->cd();
    TDirectory* trackdirectory = rundirectory->mkdir("tracks");
    trackdirectory->cd();

    trackSignal->Write();
    trackSignalMap->Write();
    trackSignalMapModPitch->Write();
    trackSignalTDC->Write();
    tracksperevent->Write();
    tracksvsevents->Write();
    highchanneldistribution->Write();

    signalLeft2->Write();
    signalLeft1->Write();
    signalCenter->Write();
    signalRight1->Write();
    signalRight2->Write();
    signalGoodEvents->Write();
    signalmapA->Write();
    signalmapB->Write();
    signalmapC->Write();
    signalmapD->Write();

    if (dosubtrackhistosexist == true)
    {
      selectedtrackSignal->Write();
      subtractedtrackSignal->Write();
    }

    // FIXME
    //selectedtrackSignalMap->Write();
    //selectedtrackSignalTDC->Write();

    TDirectory* trackcutdirectory = trackdirectory->mkdir("cuts");
    trackcutdirectory->cd();

    fiducial_discard->Write();
    fiducial_allow->Write();
    goodchannel_discard->Write();
    goodchannel_allow->Write();
    goodevent_discard->Write();
    goodevent_allow->Write();
    trackselection_discard->Write();
    trackselection_allow->Write();
    timecut_discard->Write();
    timecut_allow->Write();
    highchannel_discard->Write();
    highchannel_allow->Write();


    rundirectory->cd();
    TDirectory* clusterdirectory = rundirectory->mkdir("clusters");
    clusterdirectory->cd();

    clustersizehisto->Write();
    clusteretadistribution->Write();
    clustersignalhisto->Write();

    // done with cluster file
    f2->Close();


    // inc the counter for runs to be plotted vs voltage
    // note: residual, noise and signal have cuts, so this number may be too high!
    voltagepoint[vsvoltage_histonr]++;

  } // done run loop

  if (_debug <= 3)
  {
    cout << " " << endl;
    cout << "**********************************************" << endl;
    cout << " " << endl;
    cout << "Done run loop!" << endl;
  }


  // The multi-run things

  if (_debug <= 1)
  {
    cout << " " << endl;
    cout << "Starting Signal Area Mapping..." << endl;
  }

  // signal area map
  for (int i=0;i<voltagegraphs.size();i++)
  {

    sstream.str(string());
    histoname = "signalareamap_";
    sstream << histoname << voltagegraphs.at(i);
    int temp = voltagegraphs.at(i) % _dutrotation_total;
    int n = _dutrotation_list.at(temp);
    int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

    double l = _irradfluence_list.at(temp2);
    int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
    string k = _sensorname_list.at(temp3).c_str();

    // jump yellow:
    if (temp2==4)
    {
     temp2++; 
    }

    TH2D* signalareamap = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    char o [100];
    char title[300];
    sprintf(o,"Epi");
    sstream << o << k;

    if (k.at(0) == '2')
    {
      if (k.at(3) == 'F')
      {
	sstream.str(string());
	sprintf(o,"Fth200Y");
	sstream << o;
      } else {
	sstream.str(string());
	sprintf(o,"Mcz");
	sstream << o << k;
      }
    }
    string sensname = sstream.str();
    sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);

    signalareamap->SetXTitle("Bias Voltage [|V|]");
    signalareamap->SetYTitle("Track mod Pitch [Pitch]");
    signalareamap->SetTitle(title);

    _outputFile->cd();
    signalareamap->Write();

  }


  if (_debug <= 1)
  {
    cout << " " << endl;
    cout << "Starting Signal Area Plot Comparison..." << endl;
  }

  // draw all the signal area plots
  // do this for each sensor and each fluence
  // get the inner loop count (aka voltage points) from clustersize, as this is uncut and comes directly from the ntuple

  // the position counter in the runlist, since we will loop over all from runcount
  int pointcounter = 0;

  // book runcount^2 memory
  TCanvas* canv_signalareaplot[500];
  TLegend* lsignalareaplot[500];
  TH2F* hsignalareaplot[500];

  // the count on sensors, fluences and rotations, so the number of histograms
  for (int jj=0;jj<voltagegraphs.size();jj++)
  {

    sstream.str(string());
    histoname = "clustersizevoltage_";
    sstream << histoname << voltagegraphs.at(jj);
    TGraphErrors* clustersizevoltage = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    int looplimit = clustersizevoltage->GetN();

    canv_signalareaplot[jj] = new TCanvas("Signal Area Plots","transparent pad",600,400);
    char thename [100];
    sprintf(thename, "hsignalareaplot_%i", jj);
    hsignalareaplot[jj] = new TH2F(thename,"Signal vs. Track mod Pitch",20,0,1,20,0,100);
    hsignalareaplot[jj]->SetXTitle("Track mod Pitch [Pitch]");
    hsignalareaplot[jj]->SetYTitle("Landau MPV [ADCs]");
    hsignalareaplot[jj]->SetStats(0000);
    hsignalareaplot[jj]->Draw();

    lsignalareaplot[jj] = new TLegend(0.59,0.55,0.90,0.85);
    lsignalareaplot[jj]->SetBorderSize(1);
    lsignalareaplot[jj]->SetFillColor(0);
    lsignalareaplot[jj]->SetFillStyle(0);

    // for counting subruns
    int tempint = 0;

    for (int i=pointcounter;i<pointcounter+looplimit;i++)
    {

      sstream.str(string());
      histoname = "signalareaplot_";
      sstream << histoname << _runnumber.at(i);
      TGraphErrors* signalareaplot = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);

      gStyle->SetOptStat(kFALSE);

      signalareaplot->SetLineColor(tempint+1);
      signalareaplot->SetLineWidth(3);

      // jump yellow
      if ( tempint >= 4 )
      {
	signalareaplot->SetLineColor(tempint+2);
      }


      char o [100];
      char title[300];
      char title2[300];
      sprintf(o,"Epi");
      sstream.str(string());
      string k = _sensorname.at(i);
      sstream << o << k;

      if (k.at(0) == '2')
      {
	if (k.at(3) == 'F')
	{
	  sstream.str(string());
	  sprintf(o,"Fth200Y");
	  sstream << o;
	} else {
	  sstream.str(string());
	  sprintf(o,"Mcz");
	  sstream << o << k;
	}
      }
      string sensname = sstream.str();
      double fl = _irradfluence.at(i);
      int m = _biasvoltage.at(i);
      int n = _dutrotation.at(i);

      sprintf(title2, "%s, F=%.1e, %irot", sensname.c_str(),fl,n);
      if (tempint==0)
      {
	signalareaplot->Draw("P");
	signalareaplot->SetTitle(title2);

      } else {
	signalareaplot->Draw("P");
      }

      canv_signalareaplot[jj]->Update();

      sprintf(title, "%s, F=%.1e, %iV, %irot", sensname.c_str(),fl,m,n);
      lsignalareaplot[jj]->AddEntry(signalareaplot,title,"lep");

      tempint++;


    } // done subrun loop

    lsignalareaplot[jj]->Draw();
    canv_signalareaplot[jj]->Write();
    canv_signalareaplot[jj]->Close();

    pointcounter=pointcounter+looplimit;

  } // done voltagegraph loop 

  
  
  
  
  




  if (_debug <= 1)
  {
    cout << " " << endl;
    cout << "Starting Eta Comparison..." << endl;
  }

  // draw all the eta distributions
  // do this for each sensor and each fluence
  // get the inner loop count (aka voltage points) from clustersize, as this is uncut and comes directly from the ntuple

  // the position counter in the runlist, since we will loop over all from runcount
  pointcounter = 0;

  // book runcount^2 memory
  TCanvas* canv_eta[500];
  TLegend* leta[500];

  // the count on sensors, fluences and rotations, so the number of histograms
  for (int jj=0;jj<voltagegraphs.size();jj++)
  {

    sstream.str(string());
    histoname = "clustersizevoltage_";
    sstream << histoname << voltagegraphs.at(jj);
    TGraphErrors* clustersizevoltage = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());

    int looplimit = clustersizevoltage->GetN();

    canv_eta[jj] = new TCanvas("Eta Distributions","transparent pad",600,400);

    // get the highest integral of this canvas for scale
    // get the ymax of the normal distri for axis range
    Double_t highestentrycount = 0;
    int highrun = 0;
    Double_t yscale = 0;

    for (int i=pointcounter;i<pointcounter+looplimit;i++)
    {
      sstream.str(string());
      histoname = "unmatchedEtaIntegral_";
      sstream << histoname << _runnumber.at(i);
      TH1D* thiseta = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
      Double_t temp = thiseta->GetMaximum();
      if (temp>highestentrycount)
      {
	highestentrycount = temp;
	highrun = i;
      }

      sstream.str(string());
      histoname = "unmatchedEta_";
      sstream << histoname << _runnumber.at(i);
      TH1D* thiseta2 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
      Double_t temp2 = thiseta2->GetMaximum();
      if (temp2>yscale)
      {
	yscale = temp2;
      }

    }


    leta[jj] = new TLegend(0.59,0.55,0.90,0.85);
    leta[jj]->SetBorderSize(1);
    leta[jj]->SetFillColor(0);
    leta[jj]->SetFillStyle(0);

    // for counting subruns
    int tempint = 0;

    for (int i=pointcounter;i<pointcounter+looplimit;i++)
    {

      sstream.str(string());
      histoname = "unmatchedEta_";
      sstream << histoname << _runnumber.at(i);
      TH1D* thiseta = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);

      sstream.str(string());
      histoname = "unmatchedEtaIntegral_";
      sstream << histoname << _runnumber.at(i);
      TH1D* thiseta2 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);

      Float_t scale = highestentrycount / thiseta2->GetMaximum();
      //cout << "the scale here is : " << scale << endl;

      thiseta->Scale(scale);

      gStyle->SetOptStat(kFALSE);
      thiseta->SetStats(kFALSE);
      thiseta->SetLineColor(tempint+1);

      // jump yellow
      if ( tempint >= 4 )
      {
	thiseta->SetLineColor(tempint+2);
      }


      char o [100];
      char title[300];
      char title2[300];
      sprintf(o,"Epi");
      sstream.str(string());
      string k = _sensorname.at(i);
      sstream << o << k;

      if (k.at(0) == '2')
      {
	if (k.at(3) == 'F')
	{
	  sstream.str(string());
	  sprintf(o,"Fth200Y");
	  sstream << o;
	} else {
	  sstream.str(string());
	  sprintf(o,"Mcz");
	  sstream << o << k;
	}
      }
      string sensname = sstream.str();
      double fl = _irradfluence.at(i);
      int m = _biasvoltage.at(i);
      int n = _dutrotation.at(i);

      sprintf(title2, "%s, F=%.1e, %irot", sensname.c_str(),fl,n);
      if (tempint==0)
      {
	thiseta->Draw("E");
	thiseta->SetTitle(title2);

      } else {
	thiseta->Draw("E same");
      }
      thiseta->SetMaximum(yscale);
      canv_eta[jj]->Update();

      sprintf(title, "%s, F=%.1e, %iV, %irot", sensname.c_str(),fl,m,n);
      leta[jj]->AddEntry(thiseta,title,"lep");

      tempint++;

      double counts = 0.0;
      double integral = 0.0;
      for (int kk = 1; kk<thiseta->GetNbinsX(); kk++)
      {
	counts = thiseta->GetBinContent(kk);
	integral += counts ;
      }

    } // done subrun loop

    leta[jj]->Draw();
    canv_eta[jj]->Write();

    pointcounter=pointcounter+looplimit;

  } // done voltagegraph loop 


  if (_debug <= 1)
  {
    cout << " " << endl;
    cout << "Starting Residual X vs. Voltage Plot..." << endl;
  }


  // residual vs voltage

  // first x
  TCanvas* canvas_residualsvsvoltage_X = new TCanvas("canvas_residualsvsvoltage_X","Residuals in X vs. Voltage",200,10,700,500);
  canvas_residualsvsvoltage_X->cd();
  TH2F *hresixvolt = new TH2F("hresixvolt","Residuals in X vs. Voltage",20,0,1000,20,_cut_minXresidual,_cut_maxXresidual);
  hresixvolt->SetXTitle("Bias Voltage [|V|]");
  hresixvolt->SetYTitle("Residual in X [mm]");
  hresixvolt->SetStats(0000);
  hresixvolt->Draw();

  TLegend *lresixvolt = new TLegend(0.59,0.55,0.90,0.85);
  lresixvolt->SetBorderSize(1);
  lresixvolt->SetFillColor(0);
  lresixvolt->SetFillStyle(0);

  for (int i=0;i<voltagegraphs.size();i++)
  {

    sstream.str(string());
    histoname = "residualsXvoltage_";
    sstream << histoname << voltagegraphs.at(i);
    int temp = voltagegraphs.at(i) % _dutrotation_total;
    int n = _dutrotation_list.at(temp);
    int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

    double l = _irradfluence_list.at(temp2);
    int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
    string k = _sensorname_list.at(temp3).c_str();

    // jump yellow:
    if (temp2==4)
    {
     temp2++; 
    }

    TGraphErrors* residualsXvoltage = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());
    residualsXvoltage->SetMarkerStyle(temp3+20);
    residualsXvoltage->SetMarkerColor(temp2+1);
    residualsXvoltage->SetMarkerSize(2);
    residualsXvoltage->SetLineColor(temp2+1);
    residualsXvoltage->SetLineWidth(2);
    residualsXvoltage->SetLineStyle(temp+1);

    char o [100];
    char title[300];
    sprintf(o,"Epi");
    sstream << o << k;

    if (k.at(0) == '2')
    {
      if (k.at(3) == 'F')
      {
	sstream.str(string());
	sprintf(o,"Fth200Y");
	sstream << o;
      } else {
	sstream.str(string());
	sprintf(o,"Mcz");
	sstream << o << k;
      }
    }
    string sensname = sstream.str();
    sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
    residualsXvoltage->SetTitle(title);
    residualsXvoltage->Draw("LP");
    lresixvolt->AddEntry(residualsXvoltage,title,"lep");

  }

  lresixvolt->Draw();
  canvas_residualsvsvoltage_X->Update();
  _outputFile->cd();
  canvas_residualsvsvoltage_X->Write();
  canvas_residualsvsvoltage_X->Close();


  if (_debug <= 1)
  {
    cout << " " << endl;
    cout << "Starting Residual Y vs. Voltage Plot..." << endl;
  }


  // then y
  TCanvas* canvas_residualsvsvoltage_Y = new TCanvas("canvas_residualsvsvoltage_Y","Residuals in Y vs. Voltage",200,10,700,500);
  canvas_residualsvsvoltage_Y->cd();
  TH2F *hresiyvolt = new TH2F("hresiyvolt","Residuals in Y vs. Voltage",20,0,1000,20,_cut_minYresidual,_cut_maxYresidual);
  hresiyvolt->SetXTitle("Bias Voltage [|V|]");
  hresiyvolt->SetYTitle("Residual in Y [mm]");
  hresiyvolt->SetStats(0000);
  hresiyvolt->Draw();

  // for the y plot we also add a simple line to guide the eye where the binary resolution should be.
  // binary is : pitch / sqrt(12) plus telescope resolution
  float binaryresolution = (sqrt(80.0*80.0/12 + _cut_telescoperesolution*_cut_telescoperesolution))/1000.0;
  TLine *binaryline = new TLine(1,binaryresolution,999,binaryresolution);
  binaryline->SetLineColor(kRed);
  binaryline->SetLineStyle(3);
  binaryline->SetLineWidth(3);
  binaryline->Draw();

  TLegend *lresiyvolt = new TLegend(0.59,0.55,0.90,0.85);
  lresiyvolt->SetBorderSize(1);
  lresiyvolt->SetFillColor(0);
  lresiyvolt->SetFillStyle(0);

  for (int i=0;i<voltagegraphs.size();i++)
  {

    sstream.str(string());
    histoname = "residualsYvoltage_";
    sstream << histoname << voltagegraphs.at(i);
    int temp = voltagegraphs.at(i) % _dutrotation_total;
    int n = _dutrotation_list.at(temp);
    int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

    double l = _irradfluence_list.at(temp2);
    int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
    string k = _sensorname_list.at(temp3).c_str();

    // jump yellow:
    if (temp2==4)
    {
     temp2++; 
    }

    TGraphErrors* residualsYvoltage = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());
    residualsYvoltage->SetMarkerStyle(temp3+20);
    residualsYvoltage->SetMarkerColor(temp2+1);
    residualsYvoltage->SetMarkerSize(2);
    residualsYvoltage->SetLineColor(temp2+1);
    residualsYvoltage->SetLineWidth(2);
    residualsYvoltage->SetLineStyle(temp+1);

    char o [100];
    char title[300];
    sprintf(o,"Epi");
    sstream << o << k;

    if (k.at(0) == '2')
    {
      if (k.at(3) == 'F')
      {
	sstream.str(string());
	sprintf(o,"Fth200Y");
	sstream << o;
      } else {
	sstream.str(string());
	sprintf(o,"Mcz");
	sstream << o << k;
      }
    }
    string sensname = sstream.str();
    sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
    residualsYvoltage->SetTitle(title);
    residualsYvoltage->Draw("LP");
    lresiyvolt->AddEntry(residualsYvoltage,title,"lep");

  }


  lresiyvolt->Draw();
  canvas_residualsvsvoltage_Y->Update();
  _outputFile->cd();
  canvas_residualsvsvoltage_Y->Write();
  canvas_residualsvsvoltage_Y->Close();

  // done residual vs voltage


  if (_debug <= 1)
  {
    cout << " " << endl;
    cout << "Starting Signal vs. Voltage Plot..." << endl;
  }


  // signal vs voltage

  TCanvas* canvas_signalvoltage = new TCanvas("canvas_signalvoltage","Signal vs. Voltage",200,10,700,500);
  canvas_signalvoltage->cd();
  TH2F *hsigvolt = new TH2F("hsigvolt","Signal vs. Voltage",20,0,1000,20,0,100);
  hsigvolt->SetXTitle("Bias Voltage [|V|]");
  hsigvolt->SetYTitle("Landau MPV [ADCs]");
  hsigvolt->SetStats(0000);
  hsigvolt->Draw();

  TLegend *lsigvolt = new TLegend(0.59,0.55,0.90,0.85);
  lsigvolt->SetBorderSize(1);
  lsigvolt->SetFillColor(0);
  lsigvolt->SetFillStyle(0);

  for (int i=0;i<voltagegraphs.size();i++)
  {

    sstream.str(string());
    histoname = "signalvoltage_";
    sstream << histoname << voltagegraphs.at(i);
    int temp = voltagegraphs.at(i) % _dutrotation_total;
    int n = _dutrotation_list.at(temp);
    int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

    double l = _irradfluence_list.at(temp2);
    int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
    string k = _sensorname_list.at(temp3).c_str();

    // jump yellow:
    if (temp2==4)
    {
     temp2++; 
    }

    TGraphErrors* signalvoltage = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());
    signalvoltage->SetMarkerStyle(temp3+20);
    signalvoltage->SetMarkerColor(temp2+1);
    signalvoltage->SetMarkerSize(2);
    signalvoltage->SetLineColor(temp2+1);
    signalvoltage->SetLineWidth(2);
    signalvoltage->SetLineStyle(temp+1);

    char o [100];
    char title[300];
    sprintf(o,"Epi");
    sstream << o << k;

    if (k.at(0) == '2')
    {
      if (k.at(3) == 'F')
      {
	sstream.str(string());
	sprintf(o,"Fth200Y");
	sstream << o;
      } else {
	sstream.str(string());
	sprintf(o,"Mcz");
	sstream << o << k;
      }
    }
    string sensname = sstream.str();
    sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
    signalvoltage->SetTitle(title);
    signalvoltage->Draw("LP");
    lsigvolt->AddEntry(signalvoltage,title,"lep");

  }

  lsigvolt->Draw();
  canvas_signalvoltage->Update();
  _outputFile->cd();
  canvas_signalvoltage->Write();
  canvas_signalvoltage->Close();

  // done signal vs voltage


  if (_debug <= 1)
  {
    cout << " " << endl;
    cout << "Starting Signal/Noise vs. Voltage Plot..." << endl;
  }


  // signal/noise vs voltage

  TCanvas* canvas_snvoltage = new TCanvas("canvas_snvoltage","Signal/Noise vs. Voltage",200,10,700,500);
  canvas_snvoltage->cd();
  TH2F *hsnvolt = new TH2F("hsnvolt","Signal/noise vs. Voltage",20,0,1000,20,0,25);
  hsnvolt->SetXTitle("Bias Voltage [|V|]");
  hsnvolt->SetYTitle("Signal to Noise [1]");
  hsnvolt->SetStats(0000);
  hsnvolt->Draw();

  TLegend *lsnvolt = new TLegend(0.59,0.55,0.90,0.85);
  lsnvolt->SetBorderSize(1);
  lsnvolt->SetFillColor(0);
  lsnvolt->SetFillStyle(0);

  for (int i=0;i<voltagegraphs.size();i++)
  {

    sstream.str(string());
    histoname = "snvoltage_";
    sstream << histoname << voltagegraphs.at(i);
    int temp = voltagegraphs.at(i) % _dutrotation_total;
    int n = _dutrotation_list.at(temp);
    int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

    double l = _irradfluence_list.at(temp2);
    int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
    string k = _sensorname_list.at(temp3).c_str();

    // jump yellow:
    if (temp2==4)
    {
     temp2++; 
    }

    TGraphErrors* snvoltage = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());
    snvoltage->SetMarkerStyle(temp3+20);
    snvoltage->SetMarkerColor(temp2+1);
    snvoltage->SetMarkerSize(2);
    snvoltage->SetLineColor(temp2+1);
    snvoltage->SetLineWidth(2);
    snvoltage->SetLineStyle(temp+1);

    char o [100];
    char title[300];
    sprintf(o,"Epi");
    sstream << o << k;

    if (k.at(0) == '2')
    {
      if (k.at(3) == 'F')
      {
	sstream.str(string());
	sprintf(o,"Fth200Y");
	sstream << o;
      } else {
	sstream.str(string());
	sprintf(o,"Mcz");
	sstream << o << k;
      }
    }
    string sensname = sstream.str();
    sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
    snvoltage->SetTitle(title);
    snvoltage->Draw("LP");
    lsnvolt->AddEntry(snvoltage,title,"lep");

  }

  lsnvolt->Draw();
  canvas_snvoltage->Update();
  _outputFile->cd();
  canvas_snvoltage->Write();
  canvas_snvoltage->Close();

  // done signal vs voltage


  if (_debug <= 1)
  {
    cout << " " << endl;
    cout << "Starting Signal/Effective Thickness vs. Voltage Plot..." << endl;
  }


  // signal/noise vs voltage

  TCanvas* canvas_signaleffthickness = new TCanvas("canvas_signaleffthickness","Signal/Effective Thickness vs. Voltage",200,10,700,500);
  canvas_signaleffthickness->cd();
  TH2F *hseffthickvolt = new TH2F("hseffthickvolt","Signal/Effective Thickness vs. Voltage",20,0,1000,20,0,1);
  hseffthickvolt->SetXTitle("Bias Voltage [|V|]");
  hseffthickvolt->SetYTitle("Signal / Effective Thickness [ADC / #mum]");
  hseffthickvolt->SetStats(0000);
  hseffthickvolt->Draw();

  TLegend *leffthickvolt = new TLegend(0.59,0.55,0.90,0.85);
  leffthickvolt->SetBorderSize(1);
  leffthickvolt->SetFillColor(0);
  leffthickvolt->SetFillStyle(0);

  for (int i=0;i<voltagegraphs.size();i++)
  {

    sstream.str(string());
    histoname = "signaldistancevoltage_";
    sstream << histoname << voltagegraphs.at(i);
    int temp = voltagegraphs.at(i) % _dutrotation_total;
    int n = _dutrotation_list.at(temp);
    int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

    double l = _irradfluence_list.at(temp2);
    int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
    string k = _sensorname_list.at(temp3).c_str();

    // jump yellow:
    if (temp2==4)
    {
     temp2++; 
    }

    TGraphErrors* signaldistancevoltage = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());
    signaldistancevoltage->SetMarkerStyle(temp3+20);
    signaldistancevoltage->SetMarkerColor(temp2+1);
    signaldistancevoltage->SetMarkerSize(2);
    signaldistancevoltage->SetLineColor(temp2+1);
    signaldistancevoltage->SetLineWidth(2);
    signaldistancevoltage->SetLineStyle(temp+1);

    char o [100];
    char title[300];
    sprintf(o,"Epi");
    sstream << o << k;

    if (k.at(0) == '2')
    {
      if (k.at(3) == 'F')
      {
	sstream.str(string());
	sprintf(o,"Fth200Y");
	sstream << o;
      } else {
	sstream.str(string());
	sprintf(o,"Mcz");
	sstream << o << k;
      }
    }
    string sensname = sstream.str();
    sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
    signaldistancevoltage->SetTitle(title);
    signaldistancevoltage->Draw("LP");
    leffthickvolt->AddEntry(signaldistancevoltage,title,"lep");

  }

  leffthickvolt->Draw();
  canvas_signaleffthickness->Update();
  _outputFile->cd();
  canvas_signaleffthickness->Write();
  canvas_signaleffthickness->Close();

  // done signal vs voltage


  if (_debug <= 1)
  {
    cout << " " << endl;
    cout << "Starting Noise vs. Voltage Plot..." << endl;
  }


  // noise vs voltage
  TCanvas* canvas_noisevoltage = new TCanvas("canvas_noisevoltage","Noise vs. Voltage",200,10,700,500);
  canvas_noisevoltage->cd();
  TH2F *hnoisevolt = new TH2F("hnoisevolt","Noise vs. Voltage",20,0,1000,20,0,50);
  hnoisevolt->SetXTitle("Bias Voltage [|V|]");
  hnoisevolt->SetYTitle("Average Good Channel Noise [ADCs]");
  hnoisevolt->SetStats(0000);
  hnoisevolt->Draw();

  TLegend *lnoisevolt = new TLegend(0.59,0.55,0.90,0.85);
  lnoisevolt->SetBorderSize(1);
  lnoisevolt->SetFillColor(0);
  lnoisevolt->SetFillStyle(0);

  for (int i=0;i<voltagegraphs.size();i++)
  {

    sstream.str(string());
    histoname = "noisevoltage_";
    sstream << histoname << voltagegraphs.at(i);
    int temp = voltagegraphs.at(i) % _dutrotation_total;
    int n = _dutrotation_list.at(temp);
    int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

    double l = _irradfluence_list.at(temp2);
    int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
    string k = _sensorname_list.at(temp3).c_str();

    // jump yellow:
    if (temp2==4)
    {
     temp2++; 
    }

    TGraphErrors* noisevoltage = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());
    noisevoltage->SetMarkerStyle(temp3+20);
    noisevoltage->SetMarkerColor(temp2+1);
    noisevoltage->SetMarkerSize(2);
    noisevoltage->SetLineColor(temp2+1);
    noisevoltage->SetLineWidth(2);
    noisevoltage->SetLineStyle(temp+1);

    char o [100];
    char title[300];
    sprintf(o,"Epi");
    sstream << o << k;

    if (k.at(0) == '2')
    {
      if (k.at(3) == 'F')
      {
	sstream.str(string());
	sprintf(o,"Fth200Y");
	sstream << o;
      } else {
	sstream.str(string());
	sprintf(o,"Mcz");
	sstream << o << k;
      }
    }
    string sensname = sstream.str();
    sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
    noisevoltage->SetTitle(title);
    noisevoltage->Draw("LP");
    lnoisevolt->AddEntry(noisevoltage,title,"lep");

  }

  lnoisevolt->Draw();
  canvas_noisevoltage->Update();
  _outputFile->cd();
  canvas_noisevoltage->Write();
  canvas_noisevoltage->Close();

  // done noise vs voltage



    if (_debug <= 1)
  {
    cout << " " << endl;
    cout << "Starting RGH vs. Voltage Plot..." << endl;
  }


  // rgh vs voltage
  TCanvas* canvas_rghvoltage = new TCanvas("canvas_rghvoltage","Noise vs. Voltage",200,10,700,500);
  canvas_rghvoltage->cd();
  TH2F *hrghvolt = new TH2F("hrghvolt","RGH Percentage vs. Voltage",20,0,1000,20,0,_cut_maxrghplot);
  hrghvolt->SetXTitle("Bias Voltage [|V|]");
  hrghvolt->SetYTitle("Percentage of Random Ghost Hits [%]");
  hrghvolt->SetStats(0000);
  hrghvolt->Draw();

  TLegend *lrghvolt = new TLegend(0.59,0.55,0.90,0.85);
  lrghvolt->SetBorderSize(1);
  lrghvolt->SetFillColor(0);
  lrghvolt->SetFillStyle(0);

  for (int i=0;i<voltagegraphs.size();i++)
  {

    sstream.str(string());
    histoname = "rghvoltage_";
    sstream << histoname << voltagegraphs.at(i);
    int temp = voltagegraphs.at(i) % _dutrotation_total;
    int n = _dutrotation_list.at(temp);
    int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

    double l = _irradfluence_list.at(temp2);
    int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
    string k = _sensorname_list.at(temp3).c_str();

    // jump yellow:
    if (temp2==4)
    {
     temp2++; 
    }

    TGraphErrors* rghvoltage = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());
    rghvoltage->SetMarkerStyle(temp3+20);
    rghvoltage->SetMarkerColor(temp2+1);
    rghvoltage->SetMarkerSize(2);
    rghvoltage->SetLineColor(temp2+1);
    rghvoltage->SetLineWidth(2);
    rghvoltage->SetLineStyle(temp+1);

    char o [100];
    char title[300];
    sprintf(o,"Epi");
    sstream << o << k;

    if (k.at(0) == '2')
    {
      if (k.at(3) == 'F')
      {
	sstream.str(string());
	sprintf(o,"Fth200Y");
	sstream << o;
      } else {
	sstream.str(string());
	sprintf(o,"Mcz");
	sstream << o << k;
      }
    }
    string sensname = sstream.str();
    sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
    rghvoltage->SetTitle(title);
    rghvoltage->Draw("LP");
    lrghvolt->AddEntry(rghvoltage,title,"lep");

  }

  lrghvolt->Draw();
  canvas_rghvoltage->Update();
  _outputFile->cd();
  canvas_rghvoltage->Write();
  canvas_rghvoltage->Close();

  // done rgh vs voltage



  if (_debug <= 1)
  {
    cout << " " << endl;
    cout << "Starting Clustersize vs. Voltage Plot..." << endl;
  }


  // clustersize vs voltage
  TCanvas* canvas_clustersizevoltage = new TCanvas("canvas_clustersizevoltage","Cluster Size vs. Voltage",200,10,700,500);
  canvas_clustersizevoltage->cd();
  TH2F *hclustervolt = new TH2F("hclustervolt","Cluster Size vs. Voltage",20,0,1000,20,0,5);
  hclustervolt->SetXTitle("Bias Voltage [|V|]");
  hclustervolt->SetYTitle("Cluster Size [1]");
  hclustervolt->SetStats(0000);
  hclustervolt->Draw();

  TLegend *lclustervolt = new TLegend(0.59,0.55,0.90,0.85);
  lclustervolt->SetBorderSize(1);
  lclustervolt->SetFillColor(0);
  lclustervolt->SetFillStyle(0);

  for (int i=0;i<voltagegraphs.size();i++)
  {

    sstream.str(string());
    histoname = "clustersizevoltage_";
    sstream << histoname << voltagegraphs.at(i);
    int temp = voltagegraphs.at(i) % _dutrotation_total;
    int n = _dutrotation_list.at(temp);
    int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

    double l = _irradfluence_list.at(temp2);
    int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
    string k = _sensorname_list.at(temp3).c_str();

    // jump yellow:
    if (temp2==4)
    {
     temp2++; 
    }

    TGraphErrors* clustersizevoltage = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());
    clustersizevoltage->SetMarkerStyle(temp3+20);
    clustersizevoltage->SetMarkerColor(temp2+1);
    clustersizevoltage->SetMarkerSize(2);
    clustersizevoltage->SetLineColor(temp2+1);
    clustersizevoltage->SetLineWidth(2);
    clustersizevoltage->SetLineStyle(temp+1);

    char o [100];
    char title[300];
    sprintf(o,"Epi");
    sstream << o << k;

    if (k.at(0) == '2')
    {
      if (k.at(3) == 'F')
      {
	sstream.str(string());
	sprintf(o,"Fth200Y");
	sstream << o;
      } else {
	sstream.str(string());
	sprintf(o,"Mcz");
	sstream << o << k;
      }
    }
    string sensname = sstream.str();
    sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
    clustersizevoltage->SetTitle(title);
    clustersizevoltage->Draw("LP");
    lclustervolt->AddEntry(clustersizevoltage,title,"lep");

  }

  lclustervolt->Draw();
  canvas_clustersizevoltage->Update();
  _outputFile->cd();
  canvas_clustersizevoltage->Write();
  canvas_clustersizevoltage->Close();

  // done clustersize vs voltage



  if (_debug <= 1)
  {
    cout << " " << endl;
    cout << "Starting Matched Clusters vs. Voltage Plot..." << endl;
  }


  // clustercount vs voltage
  TCanvas* canvas_clustercountvoltage = new TCanvas("canvas_clustercountvoltage","Cluster Count vs. Voltage",200,10,700,500);
  canvas_clustercountvoltage->cd();
  TH2F *hclustercountvolt = new TH2F("hclustercountvolt","Matched Cluster Count vs. Voltage",20,0,1000,20,0,_cut_maxclustercount);
  hclustercountvolt->SetXTitle("Bias Voltage [|V|]");
  hclustercountvolt->SetYTitle("Matched Clusters per good Track [1]");
  hclustercountvolt->SetStats(0000);
  hclustercountvolt->Draw();

  TLegend *lclustercountvolt = new TLegend(0.59,0.55,0.90,0.85);
  lclustercountvolt->SetBorderSize(1);
  lclustercountvolt->SetFillColor(0);
  lclustercountvolt->SetFillStyle(0);

  for (int i=0;i<voltagegraphs.size();i++)
  {

    sstream.str(string());
    histoname = "clustercountvoltage_";
    sstream << histoname << voltagegraphs.at(i);
    int temp = voltagegraphs.at(i) % _dutrotation_total;
    int n = _dutrotation_list.at(temp);
    int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

    double l = _irradfluence_list.at(temp2);
    int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
    string k = _sensorname_list.at(temp3).c_str();

    // jump yellow:
    if (temp2==4)
    {
     temp2++; 
    }

    TGraphErrors* clustercountvoltage = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());
    clustercountvoltage->SetMarkerStyle(temp3+20);
    clustercountvoltage->SetMarkerColor(temp2+1);
    clustercountvoltage->SetMarkerSize(2);
    clustercountvoltage->SetLineColor(temp2+1);
    clustercountvoltage->SetLineWidth(2);
    clustercountvoltage->SetLineStyle(temp+1);

    char o [100];
    char title[300];
    sprintf(o,"Epi");
    sstream << o << k;

    if (k.at(0) == '2')
    {
      if (k.at(3) == 'F')
      {
	sstream.str(string());
	sprintf(o,"Fth200Y");
	sstream << o;
      } else {
	sstream.str(string());
	sprintf(o,"Mcz");
	sstream << o << k;
      }
    }
    string sensname = sstream.str();
    sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
    clustercountvoltage->SetTitle(title);
    clustercountvoltage->Draw("LP");
    lclustercountvolt->AddEntry(clustercountvoltage,title,"lep");

  }

  lclustercountvolt->Draw();
  canvas_clustercountvoltage->Update();
  _outputFile->cd();
  canvas_clustercountvoltage->Write();
  canvas_clustercountvoltage->Close();

  // done clustercount vs voltage



  if (_debug <= 1)
  {
    cout << " " << endl;
    cout << "Starting Effective Channels vs. Voltage Plot..." << endl;
  }


  // clustercount vs voltage
  TCanvas* canvas_channelcountvoltage = new TCanvas("canvas_channelcountvoltage","Channel Count vs. Voltage",200,10,700,500);
  canvas_channelcountvoltage->cd();
  TH2F *hchannelcountvolt = new TH2F("hchannelcountvolt","Effective Channels in Signal vs. Voltage",20,0,1000,20,0,500);
  hchannelcountvolt->SetXTitle("Bias Voltage [|V|]");
  hchannelcountvolt->SetYTitle("Effective Channels in Signal [1]");
  hchannelcountvolt->SetStats(0000);
  hchannelcountvolt->Draw();

  TLegend *lchannelcountvolt = new TLegend(0.59,0.55,0.90,0.85);
  lchannelcountvolt->SetBorderSize(1);
  lchannelcountvolt->SetFillColor(0);
  lchannelcountvolt->SetFillStyle(0);

  for (int i=0;i<voltagegraphs.size();i++)
  {

    sstream.str(string());
    histoname = "channelcountvoltage_";
    sstream << histoname << voltagegraphs.at(i);
    int temp = voltagegraphs.at(i) % _dutrotation_total;
    int n = _dutrotation_list.at(temp);
    int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

    double l = _irradfluence_list.at(temp2);
    int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
    string k = _sensorname_list.at(temp3).c_str();

    // jump yellow:
    if (temp2==4)
    {
     temp2++; 
    }

    TGraphErrors* channelcountvoltage = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());
    channelcountvoltage->SetMarkerStyle(temp3+20);
    channelcountvoltage->SetMarkerColor(temp2+1);
    channelcountvoltage->SetMarkerSize(2);
    channelcountvoltage->SetLineColor(temp2+1);
    channelcountvoltage->SetLineWidth(2);
    channelcountvoltage->SetLineStyle(temp+1);

    char o [100];
    char title[300];
    sprintf(o,"Epi");
    sstream << o << k;

    if (k.at(0) == '2')
    {
      if (k.at(3) == 'F')
      {
	sstream.str(string());
	sprintf(o,"Fth200Y");
	sstream << o;
      } else {
	sstream.str(string());
	sprintf(o,"Mcz");
	sstream << o << k;
      }
    }
    string sensname = sstream.str();
    sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
    channelcountvoltage->SetTitle(title);
    channelcountvoltage->Draw("LP");
    lchannelcountvolt->AddEntry(channelcountvoltage,title,"lep");

  }

  lchannelcountvolt->Draw();
  canvas_channelcountvoltage->Update();
  _outputFile->cd();
  canvas_channelcountvoltage->Write();
  canvas_channelcountvoltage->Close();

  // done channelcount vs voltage



  if (_debug <= 1)
  {
    cout << " " << endl;
    cout << "Starting Current vs. Voltage Plot..." << endl;
  }


  // current vs voltage
  TCanvas* canvas_currentvoltage = new TCanvas("canvas_currentvoltage","Sensor Current vs. Voltage",200,10,700,500);
  canvas_currentvoltage->cd();
  TH2F *hcurrentvolt = new TH2F("hcurrentvolt","Sensor Current vs. Voltage",20,0,1000,20,0,1000);
  hcurrentvolt->SetXTitle("Bias Voltage [|V|]");
  hcurrentvolt->SetYTitle("Sensor Current (at 253K) [| #muA|]");
  hcurrentvolt->SetStats(0000);
  hcurrentvolt->Draw();

  TLegend *lcurrentcountvolt = new TLegend(0.59,0.55,0.90,0.85);
  lcurrentcountvolt->SetBorderSize(1);
  lcurrentcountvolt->SetFillColor(0);
  lcurrentcountvolt->SetFillStyle(0);

  for (int i=0;i<voltagegraphs.size();i++)
  {

    sstream.str(string());
    histoname = "currentvoltage_";
    sstream << histoname << voltagegraphs.at(i);
    int temp = voltagegraphs.at(i) % _dutrotation_total;
    int n = _dutrotation_list.at(temp);
    int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

    double l = _irradfluence_list.at(temp2);
    int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
    string k = _sensorname_list.at(temp3).c_str();

    // jump yellow:
    if (temp2==4)
    {
     temp2++; 
    }

    TGraphErrors* currentvoltage = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());
    currentvoltage->SetMarkerStyle(temp3+20);
    currentvoltage->SetMarkerColor(temp2+1);
    currentvoltage->SetMarkerSize(2);
    currentvoltage->SetLineColor(temp2+1);
    currentvoltage->SetLineWidth(2);
    currentvoltage->SetLineStyle(temp+1);

    char o [100];
    char title[300];
    sprintf(o,"Epi");
    sstream << o << k;

    if (k.at(0) == '2')
    {
      if (k.at(3) == 'F')
      {
	sstream.str(string());
	sprintf(o,"Fth200Y");
	sstream << o;
      } else {
	sstream.str(string());
	sprintf(o,"Mcz");
	sstream << o << k;
      }
    }
    string sensname = sstream.str();
    sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);

    // skip rotated sensors, as they have the same current
    if (n < 10)
    {
      currentvoltage->SetTitle(title);
      currentvoltage->Draw("LP");
      lcurrentcountvolt->AddEntry(currentvoltage,title,"lep");
    }

  }

  lcurrentcountvolt->Draw();
  canvas_currentvoltage->Update();
  _outputFile->cd();
  canvas_currentvoltage->Write();
  canvas_currentvoltage->Close();

  // done current vs voltage


  // current vs voltage
  TCanvas* canvas_volumecurrent = new TCanvas("canvas_volumecurrent","Sensor Current / Thickness vs. Voltage",200,10,700,500);
  canvas_volumecurrent->cd();
  TH2F *hvolumecurrent = new TH2F("hvolumecurrent","Sensor Current vs. Voltage",20,0,1000,20,0,100000);
  hvolumecurrent->SetXTitle("Bias Voltage [|V|]");
  hvolumecurrent->SetYTitle("Sensor Current (at 253K) [| #muA| / cm^{-3}]");
  hvolumecurrent->SetStats(0000);
  hvolumecurrent->Draw();

  TLegend *lvolumecurrent = new TLegend(0.59,0.55,0.90,0.85);
  lvolumecurrent->SetBorderSize(1);
  lvolumecurrent->SetFillColor(0);
  lvolumecurrent->SetFillStyle(0);

  for (int i=0;i<voltagegraphs.size();i++)
  {

    sstream.str(string());
    histoname = "volumecurrent_";
    sstream << histoname << voltagegraphs.at(i);
    int temp = voltagegraphs.at(i) % _dutrotation_total;
    int n = _dutrotation_list.at(temp);
    int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

    double l = _irradfluence_list.at(temp2);
    int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
    string k = _sensorname_list.at(temp3).c_str();

    // jump yellow:
    if (temp2==4)
    {
     temp2++; 
    }

    TGraphErrors* volumecurrent = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());
    volumecurrent->SetMarkerStyle(temp3+20);
    volumecurrent->SetMarkerColor(temp2+1);
    volumecurrent->SetMarkerSize(2);
    volumecurrent->SetLineColor(temp2+1);
    volumecurrent->SetLineWidth(2);
    volumecurrent->SetLineStyle(temp+1);

    char o [100];
    char title[300];
    sprintf(o,"Epi");
    sstream << o << k;

    if (k.at(0) == '2')
    {
      if (k.at(3) == 'F')
      {
	sstream.str(string());
	sprintf(o,"Fth200Y");
	sstream << o;
      } else {
	sstream.str(string());
	sprintf(o,"Mcz");
	sstream << o << k;
      }
    }
    string sensname = sstream.str();
    sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);

    // skip rotated sensors, as they have the same current
    if (n < 10)
    {
      volumecurrent->SetTitle(title);
      volumecurrent->Draw("LP");
      lvolumecurrent->AddEntry(volumecurrent,title,"lep");
    }

  }

  lvolumecurrent->Draw();
  canvas_volumecurrent->Update();
  _outputFile->cd();
  canvas_volumecurrent->Write();
  canvas_volumecurrent->Close();

  // done current vs voltage





  if (_debug <= 1)
  {
    cout << " " << endl;
    cout << "Starting Charge Sharing vs. Voltage Plot..." << endl;
  }


  // chargesharing vs voltage
  TCanvas* canvas_chargesharingvoltage = new TCanvas("canvas_chargesharingvoltage","Shared Charge vs. Voltage",200,10,700,500);
  canvas_chargesharingvoltage->cd();
  TH2F *hchargesharevolt = new TH2F("hchargesharevolt","Shared Charge vs. Voltage",20,0,1000,20,0,100);
  hchargesharevolt->SetXTitle("Bias Voltage [|V|]");
  hchargesharevolt->SetYTitle("Shared Charge [%]");
  hchargesharevolt->SetStats(0000);
  hchargesharevolt->Draw();

  TLegend *lchargesharevolt = new TLegend(0.59,0.55,0.90,0.85);
  lchargesharevolt->SetBorderSize(1);
  lchargesharevolt->SetFillColor(0);
  lchargesharevolt->SetFillStyle(0);

  for (int i=0;i<voltagegraphs.size();i++)
  {

    sstream.str(string());
    histoname = "chargesharingvoltage_";
    sstream << histoname << voltagegraphs.at(i);
    int temp = voltagegraphs.at(i) % _dutrotation_total;
    int n = _dutrotation_list.at(temp);
    int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

    double l = _irradfluence_list.at(temp2);
    int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
    string k = _sensorname_list.at(temp3).c_str();

    // jump yellow:
    if (temp2==4)
    {
     temp2++; 
    }

    TGraphErrors* chargesharingvoltage = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());
    chargesharingvoltage->SetMarkerStyle(temp3+20);
    chargesharingvoltage->SetMarkerColor(temp2+1);
    chargesharingvoltage->SetMarkerSize(2);
    chargesharingvoltage->SetLineColor(temp2+1);
    chargesharingvoltage->SetLineWidth(2);
    chargesharingvoltage->SetLineStyle(temp+1);

    char o [100];
    char title[300];
    sprintf(o,"Epi");
    sstream << o << k;

    if (k.at(0) == '2')
    {
      if (k.at(3) == 'F')
      {
	sstream.str(string());
	sprintf(o,"Fth200Y");
	sstream << o;
      } else {
	sstream.str(string());
	sprintf(o,"Mcz");
	sstream << o << k;
      }
    }
    string sensname = sstream.str();
    sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
    chargesharingvoltage->SetTitle(title);
    chargesharingvoltage->Draw("LP");
    lchargesharevolt->AddEntry(chargesharingvoltage,title,"lep");

  }

  lchargesharevolt->Draw();
  canvas_chargesharingvoltage->Update();
  _outputFile->cd();
  canvas_chargesharingvoltage->Write();
  canvas_chargesharingvoltage->Close();

  // done chargesharing vs voltage


  // eta chargesharing vs voltage
  TCanvas* canvas_etachargesharingvoltage = new TCanvas("canvas_etachargesharingvoltage","Shared Charge vs. Voltage",200,10,700,500);
  canvas_etachargesharingvoltage->cd();
  TH2F *hetachargesharevolt = new TH2F("hetachargesharevolt","Shared Charge vs. Voltage",20,0,1000,20,0,100);
  hetachargesharevolt->SetXTitle("Bias Voltage [|V|]");
  hetachargesharevolt->SetYTitle("Shared Charge [%]");
  hetachargesharevolt->SetStats(0000);
  hetachargesharevolt->Draw();

  TLegend *letachargesharevolt = new TLegend(0.59,0.55,0.90,0.85);
  letachargesharevolt->SetBorderSize(1);
  letachargesharevolt->SetFillColor(0);
  letachargesharevolt->SetFillStyle(0);

  for (int i=0;i<voltagegraphs.size();i++)
  {

    sstream.str(string());
    histoname = "etachargesharingvoltage_";
    sstream << histoname << voltagegraphs.at(i);
    int temp = voltagegraphs.at(i) % _dutrotation_total;
    int n = _dutrotation_list.at(temp);
    int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

    double l = _irradfluence_list.at(temp2);
    int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
    string k = _sensorname_list.at(temp3).c_str();

    // jump yellow:
    if (temp2==4)
    {
     temp2++; 
    }

    TGraphErrors* etachargesharingvoltage = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
    sstream.str(string());
    etachargesharingvoltage->SetMarkerStyle(temp3+20);
    etachargesharingvoltage->SetMarkerColor(temp2+1);
    etachargesharingvoltage->SetMarkerSize(2);
    etachargesharingvoltage->SetLineColor(temp2+1);
    etachargesharingvoltage->SetLineWidth(2);
    etachargesharingvoltage->SetLineStyle(temp+1);

    char o [100];
    char title[300];
    sprintf(o,"Epi");
    sstream << o << k;

    if (k.at(0) == '2')
    {
      if (k.at(3) == 'F')
      {
	sstream.str(string());
	sprintf(o,"Fth200Y");
	sstream << o;
      } else {
	sstream.str(string());
	sprintf(o,"Mcz");
	sstream << o << k;
      }
    }
    string sensname = sstream.str();
    sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
    etachargesharingvoltage->SetTitle(title);
    etachargesharingvoltage->Draw("LP");
    letachargesharevolt->AddEntry(etachargesharingvoltage,title,"lep");

  }

  letachargesharevolt->Draw();
  canvas_etachargesharingvoltage->Update();
  _outputFile->cd();
  canvas_etachargesharingvoltage->Write();
  canvas_etachargesharingvoltage->Close();


  // write out the alignment histos
  TDirectory* positiondirectory = _outputFile->mkdir("Sensor Positions");
  positiondirectory->cd();
  xshift->Write();
  yshift->Write();
  zshift->Write();
  ashift->Write();
  bshift->Write();
  cshift->Write();

  // done all multi-run things

} // done read and fill function


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// this books all the histograms needed
void bookhistos()
{

  cout << "##############################################" << endl;
  cout << " " << endl;
  cout << "Booking all the histograms!" << endl;
  cout << " " << endl;

  // the names and titles for the histograms
  char name[100];
  char title[300];
  stringstream sstream;

  // we have _runcount many runs
  // details further down on which histo does what
  TH1D* noise[_runcount][128];
  TH1D* allnoise[_runcount];
  TH1D* adcnoise[_runcount];
  TH1D* fivechannoise[_runcount];

  TH1D* residualsX[_runcount];
  TH1D* residualsY[_runcount];
  TH2D* residualsXY[_runcount];
  TH2D* residualsY_Q[_runcount];

  TH1D* residualsY_clu1[_runcount];
  TH1D* residualsY_clu2[_runcount];
  TH1D* residualsY_clu3[_runcount];
  TH1D* residualsY_clu4[_runcount];

  TH2D* residualsDXvsX[_runcount];
  TH2D* residualsDXvsY[_runcount];
  TH2D* residualsDYvsX[_runcount];
  TH2D* residualsDYvsY[_runcount];

  TH2D* residualsXvsEvt[_runcount];
  TH2D* residualsYvsEvt[_runcount];

  TH3D* residualMapTemp[_runcount];
  TH2D* residualMapY[_runcount];

  TProfile* residualProfileDXvsX[_runcount];
  TProfile* residualProfileDXvsY[_runcount];
  TProfile* residualProfileDYvsX[_runcount];
  TProfile* residualProfileDYvsY[_runcount];

  TH2D* hitmapDUT[_runcount];
  TH2D* hitmapDUTmodpitch[_runcount];

  TH2D* hitmapTracks[_runcount];
  TH2D* hitmapTracksmodpitch[_runcount];

  TH2D* hitmapMatch[_runcount];
  TH2D* hitmapMatchmodpitch[_runcount];
  TH2D* hitmapMatchTracks[_runcount];

  TH1D* divergencehistoX[_runcount];
  TH1D* divergencehistoY[_runcount];

  TH2D* hitmapChargeShared0[_runcount];
  TH2D* hitmapChargeShared1[_runcount];
  TH2D* hitmapChargeShared2[_runcount];
  TH2D* hitmapChargeShared3[_runcount];
  TH2D* hitmapChargeShared4[_runcount];

  TH1D* hitmapClusterSizeA[_runcount];
  TH1D* hitmapClusterSize1[_runcount];
  TH1D* hitmapClusterSize2[_runcount];
  TH1D* hitmapClusterSize3[_runcount];
  TH1D* hitmapClusterSize4[_runcount];

  TH1D* signalClusterSize1[_runcount];
  TH1D* signalClusterSize2[_runcount];
  TH1D* signalClusterSize3[_runcount];
  TH1D* signalClusterSize4[_runcount];

  TH1D* matchedEta[_runcount];
  TH2D* matchedEta2D[_runcount];
  TH1D* unmatchedEta[_runcount];
  TH2D* unmatchedEta2D[_runcount];
  TH1D* matchedEtaIntegral[_runcount];
  TH1D* unmatchedEtaIntegral[_runcount];
  TH1D* striphitEta[_runcount];
  TH2D* striphitEta2D[_runcount];
  TH1D* striphitEtaIntegral[_runcount];
  TH1D* noiseEta[_runcount];
  TH1D* noisesubtractedEta[_runcount];
  TH1D* EtaL20[_runcount];
  TH1D* EtaL10[_runcount];
  TH1D* EtaR10[_runcount];
  TH1D* EtaR20[_runcount];
  TH1D* EtaR30[_runcount];
  TH1D* EtaR40[_runcount];
  TH1D* EtaR50[_runcount];

  TH1D* trackSignal[_runcount];
  TH2D* trackSignalMap[_runcount];
  TH2D* trackSignalMapModPitch[_runcount];
  TH2D* trackSignalTDC[_runcount];
  TH1D* tracksperevent[_runcount];
  TH2D* tracksvsevents[_runcount];
  TH1D* highchanneldistribution[_runcount];

  TH1D* selectedtrackSignal[_runcount];
  TH2D* selectedtrackSignalMap[_runcount];
  TH2D* selectedtrackSignalTDC[_runcount];

  TH1D* signalLeft2[_runcount];
  TH1D* signalLeft1[_runcount];
  TH1D* signalCenter[_runcount];
  TH1D* signalRight1[_runcount];
  TH1D* signalRight2[_runcount];
  TH1D* signalGoodEvents[_runcount];
  TH1D* signalmapA[_runcount];
  TH1D* signalmapB[_runcount];
  TH1D* signalmapC[_runcount];
  TH1D* signalmapD[_runcount];

  TGraphErrors* signalareaplot[_runcount];

  TH1D* fiducial_discard[_runcount];
  TH1D* fiducial_allow[_runcount];
  TH1D* goodchannel_discard[_runcount];
  TH1D* goodchannel_allow[_runcount];
  TH1D* goodevent_discard[_runcount];
  TH1D* goodevent_allow[_runcount];
  TH1D* trackselection_discard[_runcount];
  TH1D* trackselection_allow[_runcount];
  TH1D* timecut_discard[_runcount];
  TH1D* timecut_allow[_runcount];
  TH1D* highchannel_discard[_runcount];
  TH1D* highchannel_allow[_runcount];


  // the ones for more than one run: assume runcount^2 as maximum amount...

  // sort by voltage:
  TGraphErrors* residualsXvoltage[500];
  TGraphErrors* residualsYvoltage[500];

  TGraphErrors* signalvoltage[500];
  TGraphErrors* snvoltage[500];

  TGraphErrors* noisevoltage[500];

  TGraphErrors* rghvoltage[500];

  TGraphErrors* clustersizevoltage[500];

  TGraphErrors* clustercountvoltage[500];

  TGraphErrors* channelcountvoltage[500];

  TGraphErrors* currentvoltage[500];

  TGraphErrors* volumecurrent[500];

  TGraphErrors* chargesharingvoltage[500];

  TGraphErrors* etachargesharingvoltage[500];

  TGraphErrors* signaldistancevoltage[500];

  TH2D* signalareamap[500];

  TH1D* etamap[500];

  // only once:
  TH1D* xshift = new TH1D("xshift","Sensor Offset in X;Position [mm];Entries", 100, -5.0, 5.0);
  _rootObjectMap["xshift"] = xshift;
  TH1D* yshift = new TH1D("yshift","Sensor Offset in Y;Position [mm];Entries", 100, -5.0, 5.0);
  _rootObjectMap["yshift"] = yshift;
  TH1D* zshift = new TH1D("zshift","Sensor Offset in Z;Position [mm];Entries", 100, -10.0, 10.0);
  _rootObjectMap["zshift"] = zshift;
  TH1D* ashift = new TH1D("ashift","Sensor Offset in A;Rotation [#circ];Entries", 100, -5.0, 5.0);
  _rootObjectMap["ashift"] = ashift;
  TH1D* bshift = new TH1D("bshift","Sensor Offset in B;Rotation [#circ];Entries", 100, -5.0, 5.0);
  _rootObjectMap["bshift"] = bshift;
  TH1D* cshift = new TH1D("cshift","Sensor Offset in C;Rotation [#circ];Entries", 100, -5.0, 5.0);
  _rootObjectMap["cshift"] = cshift;



  // j is the loop var, i is the actual run number
  for ( int j = 0 ; j < _actual_runcount; j++ )
  {

    // the sensor information goes into the histogram title
    int i = _runnumber.at(j);
    double l = _irradfluence.at(j);
    int m = _biasvoltage.at(j);
    int n = _dutrotation.at(j);
    char o[100] = "Fail";
    int z = _thickness.at(j);
    string k = _sensorname.at(j).c_str();

    // slight hack: material is deducted from thickness... FIXME
    if (z <=150)
    {
      sprintf(o,"Epi");
    }
   /* if (z > 150)
    {
      // fth sensors from runnumber
      if (i>=638 && i<=644)
      {
	sprintf(o,"Fth");
      } else {
	sprintf(o,"MCz");
      }
    }*/
    stringstream sstream;
    sstream << o << k;

    if (k.at(0) == '2')
    {
      if (k.at(3) == 'F')
      {
	sstream.str(string());
	sprintf(o,"Fth200Y");
	sstream << o;
      } else {
	sstream.str(string());
	sprintf(o,"Mcz");
	sstream << o << k;
      }
    }

    // the sensor name for the title

    string sensname = sstream.str();


    // noise * 128 channels
    for (int ii=0;ii<128;ii++)
    {
      sprintf(name, "noise_%i_chan_%i", i,ii);
      sprintf(title, "DUT Off-Beam Noise, Channel %i, %s, F=%.1e, %iV, %irot, run %i;Signal [ADCs];Events", ii, sensname.c_str(),l,m,n,i);
      noise[j][ii] = new TH1D(name, title, 1000, -500, 500);
      _rootObjectMap[name] = noise[j][ii] ;
    }

    // noise of a sensor
    sprintf(name, "allnoise_%i", i);
    sprintf(title, "DUT Off-Beam Noise, %s, F=%.1e, %iV, %irot, run %i;Channel;Noise [ADCs]", sensname.c_str(),l,m,n,i);
    allnoise[j] = new TH1D(name, title, 128, 0, 127);
    _rootObjectMap[name] = allnoise[j] ;

    // noise of a sensor, all channels
    sprintf(name, "adcnoise_%i", i);
    sprintf(title, "DUT Off-Beam Noise, %s, F=%.1e, %iV, %irot, run %i;Noise [ADCs];Events", sensname.c_str(),l,m,n,i);
    adcnoise[j] = new TH1D(name, title, 1000, -500, 500);
    _rootObjectMap[name] = adcnoise[j] ;

    // noise of a sensor, 5 channels summed
    sprintf(name, "fivechannoise_%i", i);
    sprintf(title, "DUT Off-Beam Noise in 5 Channels, %s, F=%.1e, %iV, %irot, run %i;Noise [ADCs];Events", sensname.c_str(),l,m,n,i);
    fivechannoise[j] = new TH1D(name, title, 1000, -500, 500);
    _rootObjectMap[name] = fivechannoise[j] ;


    // residuals
    sprintf(name, "residualsX_%i", i);
    sprintf(title, "DUT Residual in X, %s, F=%.1e, %iV, %irot, run %i;X_{hit - track} [mm];Events", sensname.c_str(),l,m,n,i);
    residualsX[j] = new TH1D(name, title, 100, -0.5, 0.5);
    _rootObjectMap[name] = residualsX[j] ;

    sprintf(name, "residualsY_%i", i);
    sprintf(title, "DUT Residual in Y, %s, F=%.1e, %iV, %irot, run %i;Y_{hit - track} [mm];Events", sensname.c_str(),l,m,n,i);
    residualsY[j] = new TH1D(name, title, 100, -0.5, 0.5);
    _rootObjectMap[name] = residualsY[j] ;

    sprintf(name, "residualsXY_%i", i);
    sprintf(title, "DUT Residual in XY, %s, F=%.1e, %iV, %irot, run %i;X_{hit - track} [mm];Y_{hit - track} [mm];Events", sensname.c_str(),l,m,n,i);
    residualsXY[j] = new TH2D(name, title, 50, -0.5, 0.5, 50, -0.5, 0.5);
    _rootObjectMap[name] = residualsXY[j] ;

    sprintf(name, "residualsY_Q_%i", i);
    sprintf(title, "DUT Residual in Y vs. Hit Charge, %s, F=%.1e, %iV, %irot, run %i;Q_{hit} [ADCs];Y_{hit - track} [mm];Events", sensname.c_str(),l,m,n,i);
    residualsY_Q[j] = new TH2D(name, title, 200, -100, 100, 50, -0.5, 0.5);
    _rootObjectMap[name] = residualsY_Q[j] ;

    sprintf(name, "residualsY_clu1_%i", i);
    sprintf(title, "DUT Residual in Y from Clustersize 1, %s, F=%.1e, %iV, %irot, run %i;Y_{hit - track} [mm];Events", sensname.c_str(),l,m,n,i);
    residualsY_clu1[j] = new TH1D(name, title, 100, -0.5, 0.5);
    _rootObjectMap[name] = residualsY_clu1[j] ;

    sprintf(name, "residualsY_clu2_%i", i);
    sprintf(title, "DUT Residual in Y from Clustersize 2, %s, F=%.1e, %iV, %irot, run %i;Y_{hit - track} [mm];Events", sensname.c_str(),l,m,n,i);
    residualsY_clu2[j] = new TH1D(name, title, 100, -0.5, 0.5);
    _rootObjectMap[name] = residualsY_clu2[j] ;

    sprintf(name, "residualsY_clu3_%i", i);
    sprintf(title, "DUT Residual in Y from Clustersize 3, %s, F=%.1e, %iV, %irot, run %i;Y_{hit - track} [mm];Events", sensname.c_str(),l,m,n,i);
    residualsY_clu3[j] = new TH1D(name, title, 100, -0.5, 0.5);
    _rootObjectMap[name] = residualsY_clu3[j] ;

    sprintf(name, "residualsY_clu4_%i", i);
    sprintf(title, "DUT Residual in Y from Clustersize 4, %s, F=%.1e, %iV, %irot, run %i;Y_{hit - track} [mm];Events", sensname.c_str(),l,m,n,i);
    residualsY_clu4[j] = new TH1D(name, title, 100, -0.5, 0.5);
    _rootObjectMap[name] = residualsY_clu4[j] ;

    // residuals vs xy
    sprintf(name, "residualsDXvsX_%i" ,i);
    sprintf(title, "DUT Residual in X vs. Track X, %s, F=%.1e, %iV, %irot, run %i;X_{hit - track} [mm];X_{track} [mm];Events", sensname.c_str(),l,m,n,i);
    residualsDXvsX[j] = new TH2D(name, title, 100, -0.5, 0.5, 100, -15.0, 15.0);
    _rootObjectMap[name] = residualsDXvsX[j] ;

    sprintf(name, "residualsDXvsY_%i" ,i);
    sprintf(title, "DUT Residual in X vs. Track Y, %s, F=%.1e, %iV, %irot, run %i;X_{hit - track} [mm];Y_{track} [mm];Events", sensname.c_str(),l,m,n,i);
    residualsDXvsY[j] = new TH2D(name, title, 100, -0.5, 0.5, 100, -10.0, 10.0);
    _rootObjectMap[name] = residualsDXvsY[j] ;

    sprintf(name, "residualsDYvsX_%i" ,i);
    sprintf(title, "DUT Residual in Y vs. Track X, %s, F=%.1e, %iV, %irot, run %i;Y_{hit - track} [mm];X_{track} [mm];Events", sensname.c_str(),l,m,n,i);
    residualsDYvsX[j] = new TH2D(name, title, 100, -0.5, 0.5, 100, -15.0, 15.0);
    _rootObjectMap[name] = residualsDYvsX[j] ;

    sprintf(name, "residualsDYvsY_%i" ,i);
    sprintf(title, "DUT Residual in Y vs. Track Y, %s, F=%.1e, %iV, %irot, run %i;Y_{hit - track} [mm];Y_{track} [mm];Events", sensname.c_str(),l,m,n,i);
    residualsDYvsY[j] = new TH2D(name, title, 100, -0.5, 0.5, 100, -10.0, 10.0);
    _rootObjectMap[name] = residualsDYvsY[j] ;

    // residual profiles
    sprintf(name, "residualProfileDXvsX_%i", i);
    sprintf(title, "DUT Residual in X vs. Track X - Profile, %s, F=%.1e, %iV, %irot, run %i;X_{track} [mm];X_{hit - track} [mm]", sensname.c_str(),l,m,n,i);
    residualProfileDXvsX[j] =  new TProfile (name, title,150,-15,15,-0.2,0.2,"s");
    _rootObjectMap[name] = residualProfileDXvsX[j] ;

    sprintf(name, "residualProfileDXvsY_%i", i);
    sprintf(title, "DUT Residual in X vs. Track Y - Profile, %s, F=%.1e, %iV, %irot, run %i;Y_{track} [mm];X_{hit - track} [mm]", sensname.c_str(),l,m,n,i);
    residualProfileDXvsY[j] =  new TProfile (name, title,150,-15,15,-0.2,0.2,"s");
    _rootObjectMap[name] = residualProfileDXvsY[j] ;

    sprintf(name, "residualProfileDYvsX_%i", i);
    sprintf(title, "DUT Residual in Y vs. Track X - Profile, %s, F=%.1e, %iV, %irot, run %i;X_{track} [mm];Y_{hit - track} [mm]", sensname.c_str(),l,m,n,i);
    residualProfileDYvsX[j] =  new TProfile (name, title,150,-15,15,-0.2,0.2,"s");
    _rootObjectMap[name] = residualProfileDYvsX[j] ;

    sprintf(name, "residualProfileDYvsY_%i", i);
    sprintf(title, "DUT Residual in Y vs. Track Y - Profile, %s, F=%.1e, %iV, %irot, run %i;Y_{track} [mm];Y_{hit - track} [mm]", sensname.c_str(),l,m,n,i);
    residualProfileDYvsY[j] =  new TProfile (name, title,150,-15,15,-0.2,0.2,"s");
    _rootObjectMap[name] = residualProfileDYvsY[j] ;

    // residuals vs event
    sprintf(name, "residualsXvsEvt_%i" ,i);
    sprintf(title, "DUT Residual in X vs. Event Nr., %s, F=%.1e, %iV, %irot, run %i;Event Nr.;X_{hit - track} [mm];Events", sensname.c_str(),l,m,n,i);
    residualsXvsEvt[j] = new TH2D(name, title, 500, 0, 5e5, 100, -0.5, 0.5);
    _rootObjectMap[name] = residualsXvsEvt[j] ;

    sprintf(name, "residualsYvsEvt_%i" ,i);
    sprintf(title, "DUT Residual in Y vs. Event Nr., %s, F=%.1e, %iV, %irot, run %i;Event Nr.;X_{hit - track} [mm];Events", sensname.c_str(),l,m,n,i);
    residualsYvsEvt[j] = new TH2D(name, title, 500, 0, 5e5, 100, -0.5, 0.5);
    _rootObjectMap[name] = residualsYvsEvt[j] ;

    // residualmaps
    sprintf(name, "residualMapTemp_%i" ,i);
    sprintf(title, "DUT Residual in X vs. Track Position, %s, F=%.1e, %iV, %irot, run %i;Event Nr.;X_{track} [mm];Y_{track} [mm];X_{hit - track} [mm];Events", sensname.c_str(),l,m,n,i);
    residualMapTemp[j] = new TH3D(name, title, 100, -10, 10, 100, -10, 10, 100, -0.1, 0.1);
    _rootObjectMap[name] = residualMapTemp[j] ;

    sprintf(name, "residualMapY_%i" ,i);
    sprintf(title, "DUT Residual in Y vs. Track Position, %s, F=%.1e, %iV, %irot, run %i;Event Nr.;X_{track} [mm];Y_{track} [mm];Residual", sensname.c_str(),l,m,n,i);
    residualMapY[j] = new TH2D(name, title, 100, -10, 10, 100, -10, 10);
    _rootObjectMap[name] = residualMapY[j] ;


    // DUT hitmap
    sprintf(name, "hitmapDUT_%i" ,i);
    sprintf(title, "Measured Hits on the DUT, %s, F=%.1e, %iV, %irot, run %i;X_{hit} [mm];Y_{hit} [mm];Events", sensname.c_str(),l,m,n,i);
    hitmapDUT[j] = new TH2D(name, title, 1000, -20, 20, 1000, -20, 20);
    _rootObjectMap[name] = hitmapDUT[j] ;

    sprintf(name, "hitmapDUTmodpitch_%i" ,i);
    sprintf(title, "Measured Hits on the DUT mod Pitch, %s, F=%.1e, %iV, %irot, run %i;X_{hit} [Pitch];Y_{hit} [Pitch];Events", sensname.c_str(),l,m,n,i);
    hitmapDUTmodpitch[j] = new TH2D(name, title, 100, 0, 1, 100, 0, 1);
    _rootObjectMap[name] = hitmapDUTmodpitch[j] ;

    // track hitmap
    sprintf(name, "hitmapTracks_%i" ,i);
    sprintf(title, "Tracks on the DUT, %s, F=%.1e, %iV, %irot, run %i;X_{track} [mm];Y_{track} [mm];Events", sensname.c_str(),l,m,n,i);
    hitmapTracks[j] = new TH2D(name, title, 1000, -20, 20, 1000, -20, 20);
    _rootObjectMap[name] = hitmapTracks[j] ;

    sprintf(name, "hitmapTracksmodpitch_%i" ,i);
    sprintf(title, "Tracks on the DUT mod Pitch, %s, F=%.1e, %iV, %irot, run %i;X_{track} [mm];Y_{track} [mm];Events", sensname.c_str(),l,m,n,i);
    hitmapTracksmodpitch[j] = new TH2D(name, title, 100, 0, 1, 100, 0, 1);
    _rootObjectMap[name] = hitmapTracksmodpitch[j] ;

    // matched hitmap
    sprintf(name, "hitmapMatch_%i" ,i);
    sprintf(title, "Matched Hits on the DUT, %s, F=%.1e, %iV, %irot, run %i;X_{hit} [mm];Y_{hit} [mm];Events", sensname.c_str(),l,m,n,i);
    hitmapMatch[j] = new TH2D(name, title, 1000, -20, 20, 1000, -20, 20);
    _rootObjectMap[name] = hitmapMatch[j] ;

    sprintf(name, "hitmapMatchmodpitch_%i" ,i);
    sprintf(title, "Matched Hits on the DUT mod Pitch, %s, F=%.1e, %iV, %irot, run %i;X_{hit} [Pitch];Y_{hit} [Pitch];Events", sensname.c_str(),l,m,n,i);
    hitmapMatchmodpitch[j] = new TH2D(name, title, 100, 0, 1, 100, 0, 1);
    _rootObjectMap[name] = hitmapMatchmodpitch[j];

    sprintf(name, "hitmapMatchTracks_%i" ,i);
    sprintf(title, "Matched Hits on the DUT, Track Positions, %s, F=%.1e, %iV, %irot, run %i;X_{track} [mm];Y_{track} [mm];Events", sensname.c_str(),l,m,n,i);
    hitmapMatchTracks[j] = new TH2D(name, title, 1000, -20, 20, 1000, -20, 20);
    _rootObjectMap[name] = hitmapMatchTracks[j] ;

    // divergence of the beam
    sprintf(name, "divergencehistoX_%i" ,i);
    sprintf(title, "Beam Divergence, Upstream Triplet X, %s, F=%.1e, %iV, %irot, run %i;X_{hit 0} - X_{hit 2} [mm];Events", sensname.c_str(),l,m,n,i);
    divergencehistoX[j] = new TH1D(name, title, 100, -5, 5);
    _rootObjectMap[name] = divergencehistoX[j] ;

    sprintf(name, "divergencehistoY_%i" ,i);
    sprintf(title, "Beam Divergence, Upstream Triplet Y, %s, F=%.1e, %iV, %irot, run %i;Y_{hit 0} - Y_{hit 2} [mm];Events", sensname.c_str(),l,m,n,i);
    divergencehistoY[j] = new TH1D(name, title, 100, -5, 5);
    _rootObjectMap[name] = divergencehistoY[j] ;

    // charge sharing to 0, 1, 2, 3, 4 neighbours
    sprintf(name, "hitmapChargeShared0_%i" ,i);
    sprintf(title, "Tracks with Charge Sharing to No Neighbour, %s, F=%.1e, %iV, %irot, run %i;X_{track} [mm];Y_{track} [mm];Events", sensname.c_str(),l,m,n,i);
    hitmapChargeShared0[j] = new TH2D(name, title, 1000, -20, 20, 1000, -20, 20);
    _rootObjectMap[name] = hitmapChargeShared0[j] ;

    sprintf(name, "hitmapChargeShared1_%i" ,i);
    sprintf(title, "Tracks with Charge Sharing to 1 Neighbour, %s, F=%.1e, %iV, %irot, run %i;X_{track} [mm];Y_{track} [mm];Events", sensname.c_str(),l,m,n,i);
    hitmapChargeShared1[j] = new TH2D(name, title, 1000, -20, 20, 1000, -20, 20);
    _rootObjectMap[name] = hitmapChargeShared1[j] ;

    sprintf(name, "hitmapChargeShared2_%i" ,i);
    sprintf(title, "Tracks with Charge Sharing to 2 Neighbours, %s, F=%.1e, %iV, %irot, run %i;X_{track} [mm];Y_{track} [mm];Events", sensname.c_str(),l,m,n,i);
    hitmapChargeShared2[j] = new TH2D(name, title, 1000, -20, 20, 1000, -20, 20);
    _rootObjectMap[name] = hitmapChargeShared2[j] ;

    sprintf(name, "hitmapChargeShared3_%i" ,i);
    sprintf(title, "Tracks with Charge Sharing to 3 Neighbours, %s, F=%.1e, %iV, %irot, run %i;X_{track} [mm];Y_{track} [mm];Events", sensname.c_str(),l,m,n,i);
    hitmapChargeShared3[j] = new TH2D(name, title, 1000, -20, 20, 1000, -20, 20);
    _rootObjectMap[name] = hitmapChargeShared3[j] ;

    sprintf(name, "hitmapChargeShared4_%i" ,i);
    sprintf(title, "Tracks with Charge Sharing to 4 Neighbours, %s, F=%.1e, %iV, %irot, run %i;X_{track} [mm];Y_{track} [mm];Events", sensname.c_str(),l,m,n,i);
    hitmapChargeShared4[j] = new TH2D(name, title, 1000, -20, 20, 1000, -20, 20);
    _rootObjectMap[name] = hitmapChargeShared4[j] ;

    sprintf(name, "hitmapClusterSizeA_%i" ,i);
    sprintf(title, "Any cluster size on the DUT mod Pitch, %s, F=%.1e, %iV, %irot, run %i;Y_{track} [Pitch];Events", sensname.c_str(),l,m,n,i);
    hitmapClusterSizeA[j] = new TH1D(name, title, 100, 0, 1);
    _rootObjectMap[name] = hitmapClusterSizeA[j] ;

    sprintf(name, "hitmapClusterSize1_%i" ,i);
    sprintf(title, "CS1 on the DUT mod Pitch, %s, F=%.1e, %iV, %irot, run %i;Y_{track} [Pitch];Events", sensname.c_str(),l,m,n,i);
    hitmapClusterSize1[j] = new TH1D(name, title, 100, 0, 1);
    _rootObjectMap[name] = hitmapClusterSize1[j] ;

    sprintf(name, "hitmapClusterSize2_%i" ,i);
    sprintf(title, "CS2 on the DUT mod Pitch, %s, F=%.1e, %iV, %irot, run %i;Y_{track} [Pitch];Events", sensname.c_str(),l,m,n,i);
    hitmapClusterSize2[j] = new TH1D(name, title, 100, 0, 1);
    _rootObjectMap[name] = hitmapClusterSize2[j] ;

    sprintf(name, "hitmapClusterSize3_%i" ,i);
    sprintf(title, "CS3 on the DUT mod Pitch, %s, F=%.1e, %iV, %irot, run %i;Y_{track} [Pitch];Events", sensname.c_str(),l,m,n,i);
    hitmapClusterSize3[j] = new TH1D(name, title, 100, 0, 1);
    _rootObjectMap[name] = hitmapClusterSize3[j] ;

    sprintf(name, "hitmapClusterSize4_%i" ,i);
    sprintf(title, "CS4 on the DUT mod Pitch, %s, F=%.1e, %iV, %irot, run %i;Y_{track} [Pitch];Events", sensname.c_str(),l,m,n,i);
    hitmapClusterSize4[j] = new TH1D(name, title, 100, 0, 1);
    _rootObjectMap[name] = hitmapClusterSize4[j] ;

    sprintf(name, "signalClusterSize1_%i" ,i);
    sprintf(title, "CS1 Signal, %s, F=%.1e, %iV, %irot, run %i;Signal [ADCs];Events", sensname.c_str(),l,m,n,i);
    signalClusterSize1[j] = new TH1D(name, title, 1000, -500, 500);
    _rootObjectMap[name] = signalClusterSize1[j] ;

    sprintf(name, "signalClusterSize2_%i" ,i);
    sprintf(title, "CS2 Signal, %s, F=%.1e, %iV, %irot, run %i;Signal [ADCs];Events", sensname.c_str(),l,m,n,i);
    signalClusterSize2[j] = new TH1D(name, title, 1000, -500, 500);
    _rootObjectMap[name] = signalClusterSize2[j] ;

    sprintf(name, "signalClusterSize3_%i" ,i);
    sprintf(title, "CS3 Signal, %s, F=%.1e, %iV, %irot, run %i;Signal [ADCs];Events", sensname.c_str(),l,m,n,i);
    signalClusterSize3[j] = new TH1D(name, title, 1000, -500, 500);
    _rootObjectMap[name] = signalClusterSize3[j] ;

    sprintf(name, "signalClusterSize4_%i" ,i);
    sprintf(title, "CS4 Signal, %s, F=%.1e, %iV, %irot, run %i;Signal [ADCs];Events", sensname.c_str(),l,m,n,i);
    signalClusterSize4[j] = new TH1D(name, title, 1000, -500, 500);
    _rootObjectMap[name] = signalClusterSize4[j] ;


    // eta distribution if a hit is matched to a track
    sprintf(name, "matchedEta_%i", i);
    sprintf(title, "Eta Distribution of Tracks, if DUT Hit Matched, %s, F=%.1e, %iV, %irot, run %i;Eta [1];Events", sensname.c_str(),l,m,n,i);
    matchedEta[j] = new TH1D(name, title, 120, -1, 2);
    _rootObjectMap[name] = matchedEta[j] ;

    sprintf(name, "matchedEta2D_%i" ,i);
    sprintf(title, "Eta Distribution of Tracks, if DUT Hit Matched vs. Track mod Pitch , %s, F=%.1e, %iV, %irot, run %i;Eta [1];Y_{track} [Pitch];Events", sensname.c_str(),l,m,n,i);
    matchedEta2D[j] = new TH2D(name, title, 120, -1, 2, 100, 0, 1);
    _rootObjectMap[name] = matchedEta2D[j];

    sprintf(name, "matchedEtaIntegral_%i", i);
    sprintf(title, "Eta Distribution of Tracks, if DUT Hit Matched, Integral, %s, F=%.1e, %iV, %irot, run %i;Eta [1];Events", sensname.c_str(),l,m,n,i);
    matchedEtaIntegral[j] = new TH1D(name, title, 120, -1, 2);
    _rootObjectMap[name] = matchedEtaIntegral[j] ;

    // eta distribution if a hit is NOT matched to a track
    sprintf(name, "unmatchedEta_%i", i);
    sprintf(title, "Eta Distribution of Tracks, %s, F=%.1e, %iV, %irot, run %i;Eta [1];Events", sensname.c_str(),l,m,n,i);
    unmatchedEta[j] = new TH1D(name, title, 120, -1, 2);
    _rootObjectMap[name] = unmatchedEta[j] ;

    sprintf(name, "unmatchedEta2D_%i" ,i);
    sprintf(title, "Eta Distribution of Tracks vs. Track mod Pitch , %s, F=%.1e, %iV, %irot, run %i;Eta [1];Y_{track} [Pitch];Events", sensname.c_str(),l,m,n,i);
    unmatchedEta2D[j] = new TH2D(name, title, 120, -1, 2, 100, 0, 1);
    _rootObjectMap[name] = unmatchedEta2D[j];

    sprintf(name, "unmatchedEtaIntegral_%i", i);
    sprintf(title, "Eta Distribution of Tracks, Integral, %s, F=%.1e, %iV, %irot, run %i;Eta [1];Events", sensname.c_str(),l,m,n,i);
    unmatchedEtaIntegral[j] = new TH1D(name, title, 120, -1, 2);
    _rootObjectMap[name] = unmatchedEtaIntegral[j] ;

    // eta distribution if a track hits a strip
    sprintf(name, "striphitEta_%i", i);
    sprintf(title, "Eta Distribution of Tracks Hitting a Strip, %s, F=%.1e, %iV, %irot, run %i;Eta [1];Events", sensname.c_str(),l,m,n,i);
    striphitEta[j] = new TH1D(name, title, 120, -1, 2);
    _rootObjectMap[name] = striphitEta[j] ;

    sprintf(name, "striphitEta2D_%i" ,i);
    sprintf(title, "Eta Distribution of Tracks Hitting a Strip vs. Track mod Pitch , %s, F=%.1e, %iV, %irot, run %i;Eta [1];Y_{track} [Pitch];Events", sensname.c_str(),l,m,n,i);
    striphitEta2D[j] = new TH2D(name, title, 120, -1, 2, 100, 0, 1);
    _rootObjectMap[name] = striphitEta2D[j];

    sprintf(name, "striphitEtaIntegral_%i", i);
    sprintf(title, "Eta Distribution of Tracks, Integral, %s, F=%.1e, %iV, %irot, run %i;Eta [1];Events", sensname.c_str(),l,m,n,i);
    striphitEtaIntegral[j] = new TH1D(name, title, 120, -1, 2);
    _rootObjectMap[name] = striphitEtaIntegral[j] ;

    // eta based on noise
    sprintf(name, "noiseEta_%i", i);
    sprintf(title, "Eta Distribution of Noise, Integral, %s, F=%.1e, %iV, %irot, run %i;Eta [1];Events", sensname.c_str(),l,m,n,i);
    noiseEta[j] = new TH1D(name, title, 120, -1, 2);
    _rootObjectMap[name] = noiseEta[j] ;

    sprintf(name, "noisesubtractedEta_%i", i);
    sprintf(title, "Eta Distribution of Tracks with noise subtracted, %s, F=%.1e, %iV, %irot, run %i;Eta [1];Events", sensname.c_str(),l,m,n,i);
    noisesubtractedEta[j] = new TH1D(name, title, 120, -1, 2);
    _rootObjectMap[name] = noisesubtractedEta[j] ;

    sprintf(name, "EtaL20_%i", i);
    sprintf(title, "Eta Distribution of Tracks, Track impact shifted left 0.2*pitch, %s, F=%.1e, %iV, %irot, run %i;Eta [1];Events", sensname.c_str(),l,m,n,i);
    EtaL20[j] = new TH1D(name, title, 120, -1, 2);
    _rootObjectMap[name] = EtaL20[j] ;

    sprintf(name, "EtaL10_%i", i);
    sprintf(title, "Eta Distribution of Tracks, Track impact shifted left 0.1*pitch, %s, F=%.1e, %iV, %irot, run %i;Eta [1];Events", sensname.c_str(),l,m,n,i);
    EtaL10[j] = new TH1D(name, title, 120, -1, 2);
    _rootObjectMap[name] = EtaL10[j] ;

    sprintf(name, "EtaR10_%i", i);
    sprintf(title, "Eta Distribution of Tracks, Track impact shifted right 0.1*pitch, %s, F=%.1e, %iV, %irot, run %i;Eta [1];Events", sensname.c_str(),l,m,n,i);
    EtaR10[j] = new TH1D(name, title, 120, -1, 2);
    _rootObjectMap[name] = EtaR10[j] ;

    sprintf(name, "EtaR20_%i", i);
    sprintf(title, "Eta Distribution of Tracks, Track impact shifted right 0.2*pitch, %s, F=%.1e, %iV, %irot, run %i;Eta [1];Events", sensname.c_str(),l,m,n,i);
    EtaR20[j] = new TH1D(name, title, 120, -1, 2);
    _rootObjectMap[name] = EtaR20[j] ;

    sprintf(name, "EtaR30_%i", i);
    sprintf(title, "Eta Distribution of Tracks, Track impact shifted right 0.3*pitch, %s, F=%.1e, %iV, %irot, run %i;Eta [1];Events", sensname.c_str(),l,m,n,i);
    EtaR30[j] = new TH1D(name, title, 120, -1, 2);
    _rootObjectMap[name] = EtaR30[j] ;

    sprintf(name, "EtaR40_%i", i);
    sprintf(title, "Eta Distribution of Tracks, Track impact shifted right 0.4*pitch, %s, F=%.1e, %iV, %irot, run %i;Eta [1];Events", sensname.c_str(),l,m,n,i);
    EtaR40[j] = new TH1D(name, title, 120, -1, 2);
    _rootObjectMap[name] = EtaR40[j] ;

    sprintf(name, "EtaR50_%i", i);
    sprintf(title, "Eta Distribution of Tracks, Track impact shifted right 0.5*pitch, %s, F=%.1e, %iV, %irot, run %i;Eta [1];Events", sensname.c_str(),l,m,n,i);
    EtaR50[j] = new TH1D(name, title, 120, -1, 2);
    _rootObjectMap[name] = EtaR50[j] ;



    // the adc count under a track
    sprintf(name, "trackSignal_%i", i);
    sprintf(title, "Signal Under a Track, %s, F=%.1e, %iV, %irot, run %i;(-1) * Signal [ADCs];Events", sensname.c_str(),l,m,n,i);
    trackSignal[j] = new TH1D(name, title, 1000, -500, 500);
    _rootObjectMap[name] = trackSignal[j] ;

    // the adc count under a track vs position
    sprintf(name, "trackSignalMap_%i" ,i);
    sprintf(title, "ADCs Under a Track vs. Track Position, %s, F=%.1e, %iV, %irot, run %i;Y_{track} [mm];(-1) * Signal [ADCs];Events", sensname.c_str(),l,m,n,i);
    trackSignalMap[j] = new TH2D(name, title, 1000, -20, 20, 1000, -500, 500);
    _rootObjectMap[name] = trackSignalMap[j];

    // the adc count under a track vs position mod pitch
    sprintf(name, "trackSignalMapModPitch_%i" ,i);
    sprintf(title, "ADCs Under a Track vs. Track Position Mod Pitch, %s, F=%.1e, %iV, %irot, run %i;Y_{track} [Pitch];(-1) * Signal [ADCs];Events", sensname.c_str(),l,m,n,i);
    trackSignalMapModPitch[j] = new TH2D(name, title, 4, 0, 1, 35, -50, 300);
    _rootObjectMap[name] = trackSignalMapModPitch[j];

    // the adc count vs tdc of the event
    sprintf(name, "trackSignalTDC_%i" ,i);
    sprintf(title, "ADCs Under a Track vs. Event TDC, %s, F=%.1e, %iV, %irot, run %i;TDC [ns];(-1) * Signal [ADCs];Events", sensname.c_str(),l,m,n,i);
    trackSignalTDC[j] = new TH2D(name, title, 100, 0, 100, 1000, -500, 500);
    _rootObjectMap[name] = trackSignalTDC[j];

    // the track count per event
    sprintf(name, "tracksperevent_%i", i);
    sprintf(title, "Tracks per Event, %s, F=%.1e, %iV, %irot, run %i;Tracks;Events", sensname.c_str(),l,m,n,i);
    tracksperevent[j] = new TH1D(name, title, 20, 0, 20);
    _rootObjectMap[name] = tracksperevent[j] ;

    sprintf(name, "tracksvsevents_%i" ,i);
    sprintf(title, "Tracks per Event vs. Event Nr., %s, F=%.1e, %iV, %irot, run %i;Event Nr.;Tracks per Event;Events", sensname.c_str(),l,m,n,i);
    tracksvsevents[j] = new TH2D(name, title, 500, 0, 5e5, 15, 0, 14);
    _rootObjectMap[name] = tracksvsevents[j];

    // the adc count under a track, selected
    sprintf(name, "selectedtrackSignal_%i", i);
    sprintf(title, "Signal Under Selected Track, %s, F=%.1e, %iV, %irot, run %i;(-1) * Signal [ADCs];Events", sensname.c_str(),l,m,n,i);
    selectedtrackSignal[j] = new TH1D(name, title, 1000, -500, 500);
    _rootObjectMap[name] = selectedtrackSignal[j] ;

    // the adc count under a track vs position, selected
    sprintf(name, "selecttrackSignalMap_%i" ,i);
    sprintf(title, "ADCs Under Selected Track vs. Track Position, %s, F=%.1e, %iV, %irot, run %i;Y_{track} [mm];(-1) * Signal [ADCs];Events", sensname.c_str(),l,m,n,i);
    selectedtrackSignalMap[j] = new TH2D(name, title, 1000, -20, 20, 1000, -500, 500);
    _rootObjectMap[name] = selectedtrackSignalMap[j];

    // the adc count vs tdc of the event, selected
    sprintf(name, "selectedtrackSignalTDC_%i" ,i);
    sprintf(title, "ADCs Under Selected Track vs. Event TDC, %s, F=%.1e, %iV, %irot, run %i;TDC [ns];(-1) * Signal [ADCs];Events", sensname.c_str(),l,m,n,i);
    selectedtrackSignalTDC[j] = new TH2D(name, title, 100, 0, 100, 1000, -500, 500);
    _rootObjectMap[name] = selectedtrackSignalTDC[j];

    // the distribution of the highest channel in a track
    sprintf(name, "highchanneldistribution_%i", i);
    sprintf(title, "Channel in a Track with Highest Signal, %s, F=%.1e, %iV, %irot, run %i;Channel Distance from Seed;Events", sensname.c_str(),l,m,n,i);
    highchanneldistribution[j] = new TH1D(name, title, 21, -10, 10);
    _rootObjectMap[name] = highchanneldistribution[j] ;


    // the adc count under a track, interstrip sections
    sprintf(name, "signalLeft2_%i", i);
    sprintf(title, "Signal Under Selected Track, Left2 Strip, %s, F=%.1e, %iV, %irot, run %i;(-1) * Signal [ADCs];Events", sensname.c_str(),l,m,n,i);
    signalLeft2[j] = new TH1D(name, title, 500, -500, 500);
    _rootObjectMap[name] = signalLeft2[j] ;

    sprintf(name, "signalLeft1_%i", i);
    sprintf(title, "Signal Under Selected Track, Left1 Strip, %s, F=%.1e, %iV, %irot, run %i;(-1) * Signal [ADCs];Events", sensname.c_str(),l,m,n,i);
    signalLeft1[j] = new TH1D(name, title, 500, -500, 500);
    _rootObjectMap[name] = signalLeft1[j] ;

    sprintf(name, "signalCenter_%i", i);
    sprintf(title, "Signal Under Selected Track, Center Strip, %s, F=%.1e, %iV, %irot, run %i;(-1) * Signal [ADCs];Events", sensname.c_str(),l,m,n,i);
    signalCenter[j] = new TH1D(name, title, 500, -500, 500);
    _rootObjectMap[name] = signalCenter[j] ;

    sprintf(name, "signalRight1_%i", i);
    sprintf(title, "Signal Under Selected Track, Right1 Strip, %s, F=%.1e, %iV, %irot, run %i;(-1) * Signal [ADCs];Events", sensname.c_str(),l,m,n,i);
    signalRight1[j] = new TH1D(name, title, 500, -500, 500);
    _rootObjectMap[name] = signalRight1[j] ;

    sprintf(name, "signalRight2_%i", i);
    sprintf(title, "Signal Under Selected Track, Right2 Strip, %s, F=%.1e, %iV, %irot, run %i;(-1) * Signal [ADCs];Events", sensname.c_str(),l,m,n,i);
    signalRight2[j] = new TH1D(name, title, 500, -500, 500);
    _rootObjectMap[name] = signalRight2[j] ;

    sprintf(name, "signalGoodEvents_%i", i);
    sprintf(title, "Signal Under Selected Track, %s, F=%.1e, %iV, %irot, run %i;(-1) * Signal [ADCs];Events", sensname.c_str(),l,m,n,i);
    signalGoodEvents[j] = new TH1D(name, title, 500, -500, 500);
    _rootObjectMap[name] = signalGoodEvents[j] ;

    sprintf(name, "signalmapA_%i", i);
    sprintf(title, "Signal Under Selected Track, Area A, %s, F=%.1e, %iV, %irot, run %i;(-1) * Signal [ADCs];Events", sensname.c_str(),l,m,n,i);
    signalmapA[j] = new TH1D(name, title, 500, -500, 500);
    _rootObjectMap[name] = signalmapA[j] ;

    sprintf(name, "signalmapB_%i", i);
    sprintf(title, "Signal Under Selected Track, Area B, %s, F=%.1e, %iV, %irot, run %i;(-1) * Signal [ADCs];Events", sensname.c_str(),l,m,n,i);
    signalmapB[j] = new TH1D(name, title, 500, -500, 500);
    _rootObjectMap[name] = signalmapB[j] ;

    sprintf(name, "signalmapC_%i", i);
    sprintf(title, "Signal Under Selected Track, Area C, %s, F=%.1e, %iV, %irot, run %i;(-1) * Signal [ADCs];Events", sensname.c_str(),l,m,n,i);
    signalmapC[j] = new TH1D(name, title, 500, -500, 500);
    _rootObjectMap[name] = signalmapC[j] ;

    sprintf(name, "signalmapD_%i", i);
    sprintf(title, "Signal Under Selected Track, Area D, %s, F=%.1e, %iV, %irot, run %i;(-1) * Signal [ADCs];Events", sensname.c_str(),l,m,n,i);
    signalmapD[j] = new TH1D(name, title, 500, -500, 500);
    _rootObjectMap[name] = signalmapD[j] ;

    sprintf(name, "signalareaplot_%i",i);
    signalareaplot[j] = new TGraphErrors;
    _rootObjectMap[name] = signalareaplot[j];


    // histos to see what we have cut:
    sprintf(name, "fiducial_discard_%i", i);
    sprintf(title, "Tracks discarded by fiducial cut, %s, F=%.1e, %iV, %irot, run %i;1;Events", sensname.c_str(),l,m,n,i);
    fiducial_discard[j] = new TH1D(name, title, 10, 0, 10);
    _rootObjectMap[name] = fiducial_discard[j] ;

    sprintf(name, "goodchannel_discard_%i", i);
    sprintf(title, "Tracks discarded by good channel cut, %s, F=%.1e, %iV, %irot, run %i;1;Events", sensname.c_str(),l,m,n,i);
    goodchannel_discard[j] = new TH1D(name, title, 10, 0, 10);
    _rootObjectMap[name] = goodchannel_discard[j] ;

    sprintf(name, "goodevent_discard_%i", i);
    sprintf(title, "Tracks discarded by good event cut, %s, F=%.1e, %iV, %irot, run %i;1;Events", sensname.c_str(),l,m,n,i);
    goodevent_discard[j] = new TH1D(name, title, 10, 0, 10);
    _rootObjectMap[name] = goodevent_discard[j] ;

    sprintf(name, "trackselection_discard_%i", i);
    sprintf(title, "Tracks discarded by track selection cut, %s, F=%.1e, %iV, %irot, run %i;1;Events", sensname.c_str(),l,m,n,i);
    trackselection_discard[j] = new TH1D(name, title, 10, 0, 10);
    _rootObjectMap[name] = trackselection_discard[j] ;

    sprintf(name, "timecut_discard_%i", i);
    sprintf(title, "Tracks discarded by time cut, %s, F=%.1e, %iV, %irot, run %i;1;Events", sensname.c_str(),l,m,n,i);
    timecut_discard[j] = new TH1D(name, title, 10, 0, 10);
    _rootObjectMap[name] = timecut_discard[j] ;

    sprintf(name, "highchannel_discard_%i", i);
    sprintf(title, "Tracks discarded by central high channel cut, %s, F=%.1e, %iV, %irot, run %i;1;Events", sensname.c_str(),l,m,n,i);
    highchannel_discard[j] = new TH1D(name, title, 10, 0, 10);
    _rootObjectMap[name] = highchannel_discard[j] ;

    // and not cut:
    sprintf(name, "fiducial_allow_%i", i);
    sprintf(title, "Tracks allowed by fiducial cut, %s, F=%.1e, %iV, %irot, run %i;1;Events", sensname.c_str(),l,m,n,i);
    fiducial_allow[j] = new TH1D(name, title, 10, 0, 10);
    _rootObjectMap[name] = fiducial_allow[j] ;

    sprintf(name, "goodchannel_allow_%i", i);
    sprintf(title, "Tracks allowed by good channel cut, %s, F=%.1e, %iV, %irot, run %i;1;Events", sensname.c_str(),l,m,n,i);
    goodchannel_allow[j] = new TH1D(name, title, 10, 0, 10);
    _rootObjectMap[name] = goodchannel_allow[j] ;

    sprintf(name, "goodevent_allow_%i", i);
    sprintf(title, "Tracks allowed by good event cut, %s, F=%.1e, %iV, %irot, run %i;1;Events", sensname.c_str(),l,m,n,i);
    goodevent_allow[j] = new TH1D(name, title, 10, 0, 10);
    _rootObjectMap[name] = goodevent_allow[j] ;

    sprintf(name, "trackselection_allow_%i", i);
    sprintf(title, "Tracks allowed by track selection cut, %s, F=%.1e, %iV, %irot, run %i;1;Events", sensname.c_str(),l,m,n,i);
    trackselection_allow[j] = new TH1D(name, title, 10, 0, 10);
    _rootObjectMap[name] = trackselection_allow[j] ;

    sprintf(name, "timecut_allow_%i", i);
    sprintf(title, "Tracks allowed by time cut, %s, F=%.1e, %iV, %irot, run %i;1;Events", sensname.c_str(),l,m,n,i);
    timecut_allow[j] = new TH1D(name, title, 10, 0, 10);
    _rootObjectMap[name] = timecut_allow[j] ;

    sprintf(name, "highchannel_allow_%i", i);
    sprintf(title, "Tracks allowed by central high channel cut, %s, F=%.1e, %iV, %irot, run %i;1;Events", sensname.c_str(),l,m,n,i);
    highchannel_allow[j] = new TH1D(name, title, 10, 0, 10);
    _rootObjectMap[name] = highchannel_allow[j] ;


  } // done loop over all runs

  // only once -> write in main dir
  _outputFile->cd();

  for (int i = 0;i<500;i++)
  {
    // resolution vs voltage
    sprintf(name, "residualsXvoltage_%i",i);
    residualsXvoltage[i] = new TGraphErrors;
    _rootObjectMap[name] = residualsXvoltage[i];

    sprintf(name, "residualsYvoltage_%i",i);
    residualsYvoltage[i] = new TGraphErrors;
    _rootObjectMap[name] = residualsYvoltage[i];


    // signal vs. voltage
    sprintf(name, "signalvoltage_%i",i);
    signalvoltage[i] = new TGraphErrors;
    _rootObjectMap[name] = signalvoltage[i];

    // signal/noise vs. voltage
    sprintf(name, "snvoltage_%i",i);
    snvoltage[i] = new TGraphErrors;
    _rootObjectMap[name] = snvoltage[i];


    // signal/distance vs. voltage
    sprintf(name, "signaldistancevoltage_%i",i);
    signaldistancevoltage[i] = new TGraphErrors;
    _rootObjectMap[name] = signaldistancevoltage[i];


    // noise vs. voltage
    sprintf(name, "noisevoltage_%i",i);
    noisevoltage[i] = new TGraphErrors;
    _rootObjectMap[name] = noisevoltage[i];

    // rgh vs. voltage
    sprintf(name, "rghvoltage_%i",i);
    rghvoltage[i] = new TGraphErrors;
    _rootObjectMap[name] = rghvoltage[i];


    // clustersize vs. voltage
    sprintf(name, "clustersizevoltage_%i",i);
    clustersizevoltage[i] = new TGraphErrors;
    _rootObjectMap[name] = clustersizevoltage[i];


    // clustercount vs. voltage
    sprintf(name, "clustercountvoltage_%i",i);
    clustercountvoltage[i] = new TGraphErrors;
    _rootObjectMap[name] = clustercountvoltage[i];


    // channelcount vs. voltage
    sprintf(name, "channelcountvoltage_%i",i);
    channelcountvoltage[i] = new TGraphErrors;
    _rootObjectMap[name] = channelcountvoltage[i];


    // IV
    sprintf(name, "currentvoltage_%i",i);
    currentvoltage[i] = new TGraphErrors;
    _rootObjectMap[name] = currentvoltage[i];

    sprintf(name, "volumecurrent_%i",i);
    volumecurrent[i] = new TGraphErrors;
    _rootObjectMap[name] = volumecurrent[i];


    // charge sharing
    sprintf(name, "chargesharingvoltage_%i",i);
    chargesharingvoltage[i] = new TGraphErrors;
    _rootObjectMap[name] = chargesharingvoltage[i];

    // eta charge sharing
    sprintf(name, "etachargesharingvoltage_%i",i);
    etachargesharingvoltage[i] = new TGraphErrors;
    _rootObjectMap[name] = etachargesharingvoltage[i];


    // signal area map
    sprintf(name, "signalareamap_%i",i);
    signalareamap[i] = new TH2D(name, "a title", 10,0,1000,4,0,1);
    _rootObjectMap[name] = signalareamap[i];


    // eta map
    sprintf(name, "etamap_%i",i);
    etamap[i] = new TH1D(name, "a title", 120,-1,2);
    _rootObjectMap[name] = etamap[i];


  }

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


Double_t langaufun(Double_t* x, Double_t* par) // the landau gaussian convolution
{

  //Fit parameters:
  //par[0]=Width (scale) parameter of Landau density
  //par[1]=Most Probable (MP, location) parameter of Landau density
  //par[2]=Total area (integral -inf to inf, normalization constant)
  //par[3]=Width (sigma) of convoluted Gaussian function
  //
  //In the Landau distribution (represented by the CERNLIB approximation), 
  //the maximum is located at x=-0.22278298 with the location parameter=0.
  //This shift is corrected within this function, so that the actual
  //maximum is identical to the MP parameter.

  // Numeric constants
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t mpshift  = -0.22278298;       // Landau maximum location

  // Control constants
  Double_t np = 100.0;      // number of convolution steps
  Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;


  // MP shift correction
  mpc = par[1] - mpshift * par[0]; 

  // Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];

  step = (xupp-xlow) / np;

  // Convolution integral of Landau and Gaussian by sum
  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);

    xx = xupp - (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }

  return (par[2] * step * sum * invsq2pi / par[3]);
}

// the function performs a gaussian fit around the highest bin content of the histo and uses the gaus sigma and center to set the 
// fit range and starting parameters of the landau gaus fit
// the parameters are: histo to be fitted, number of sigma in the negative direction for the fit range and the same in positive direction
TF1* lanGausFit(TH1* inHist, double negSigmaFit, double posSigmaFit) // function to be used
{
  int histMaxBin = inHist->GetMaximumBin();
  double histMax = inHist->GetXaxis()->GetBinCenter(histMaxBin);

  const double minMaxPosition = 0.5; // minimum position accepted for the max
  if(histMax < minMaxPosition) // if the max is too low it probably it is in the noise
    {
      int startBin = 1 + (minMaxPosition - inHist->GetXaxis()->GetXmin()) / inHist->GetXaxis()->GetBinWidth(5);
      int endBin = inHist->GetXaxis()->GetNbins() - 2;
      double yMax = 0;
      for(int iBin = startBin; iBin < endBin; ++iBin)
	if(inHist->GetBinContent(iBin) >= yMax)
	  {
	    yMax = inHist->GetBinContent(iBin);
	    histMax = inHist->GetXaxis()->GetBinCenter(iBin);
	  }
    }

  double halfRange = 20; // guess for the range of the gaus fit

  TF1* gausFit = new TF1("gausFit", "gaus", histMax - halfRange, histMax + halfRange);
  gausFit->SetLineColor(kBlue);

  inHist->Fit(gausFit, "RQL");

  double* Gpar = gausFit->GetParameters();

  double fitR1 = Gpar[1] - Gpar[2] * negSigmaFit; // fit range from the gaus mean and sigma
  double fitR2 = Gpar[1] + Gpar[2] * posSigmaFit; //...
  double gausSig = Gpar[2]; // initial guess of the gaus sigma

  // Setting fit range and start values
  Double_t fr[2];
  Double_t sv[4], pllo[4], plhi[4];
  fr[0]=fitR1;
  fr[1]=fitR2;

  // find mpv and integral start value
  int binMin; // max and min bin number (corresponding to the range)
  int binMax;
  double intStart = 0; // start value of the integral
  double mpvStart = Gpar[1]; // start value of the mpv

  double binW = inHist->GetXaxis()->GetBinWidth(5); // bin width from a random bin
  double xMin = inHist->GetXaxis()->GetXmin();
  binMin = 1 + (fitR1 - xMin) / binW;
  binMax = 1 + (fitR2 - xMin) / binW;

  double binCont;
  // double yMax = 0; // variable used to look for the maximum (mpv start)
  for(int iBn = binMin; iBn < binMax; ++iBn) // integral in the fit range to get the start value
    {
      binCont = inHist->GetBinContent(iBn);
      intStart += binCont;
      // if(binCont > yMax) 
      // 	 {
      // 	   yMax = binCont;
      // 	   mpvStart = inHist->GetXaxis()->GetBinCenter(iBn);
      // 	 }
    }

  // starting parameters
  sv[0] = 5;//landau width
  sv[1] = mpvStart; // mpv landau
  sv[2] = intStart; // integral
  if(gausSig > sv[0])
    sv[3] = sqrt(gausSig * gausSig - sv[0] * sv[0]); // gaussian width
  else
    sv[3] = gausSig;
    
  // std::cout << "Fitting histogram " << inHist->GetName() << std::endl;
  // std::cout << "Starting parameters" << std::endl;
  // std::cout << "Landau width " << sv[0] << std::endl;
  // std::cout << "MPV          " << sv[1] << std::endl;
  // std::cout << "Area         " << sv[2] << std::endl;
  // std::cout << "Gaus sigma   " << sv[3] << std::endl;

  // parameter limits
  pllo[0]=0.01; pllo[1]=-15.0; pllo[2]=1.0; pllo[3]=gausSig * 0.1;
  plhi[0]=20.0; plhi[1]=200.0; plhi[2]=10000000.0; plhi[3]=gausSig;

  TF1* ffit = new TF1("lanGausFit", langaufun, fr[0], fr[1], 4);
  //  ffit->SetNpx(1e4);
  ffit->SetParameters(sv);
  ffit->SetParNames("Width","MPV","Area","GSigma");
  ffit->SetLineColor(kRed);
   
  for (int i = 0; i < 4; i++)
    {
      ffit->SetParLimits(i, pllo[i], plhi[i]);
    }

    //ffit->SetRange(15,150);
  inHist->Fit(ffit,"RQL"); // fit within specified range

  return ffit;
}

Double_t gausLangaufun(Double_t* x, Double_t* par) // a peak at 0 and a landau gaussian convolution
{
  Double_t gausPart = par[0] * TMath::Gaus(*x, par[1], par[5]);
  Double_t langauPart = langaufun(x, &par[2]); // par 0 to 2 belong to the gauss part

  return langauPart + gausPart;
}

// the sigma of the convolution and the one of the noise are constrained to be the same
TF1* gausLanGausFit(TH1* inHist, double negSigmaFit, double posSigmaFit)
{
  TF1* gausFunc = new TF1("gausFunc", "gaus", -30, 5); // referred as g0 in the next comments
  inHist->Fit(gausFunc, "RQL");

  TF1* langauFunc = lanGausFit(inHist, negSigmaFit, posSigmaFit);

  const int nPars = 6;
  double par[nPars] = {0};

  for(int i = 0; i < 2; ++i) par[i] = gausFunc->GetParameter(i); // get the gaus fit parameters
  for(int i = 2; i < nPars; ++i) par[i] = langauFunc->GetParameter(i - 2); // get parameters form langaus fit

  double parLimHi[nPars] = {0};
  double parLimLo[nPars] = {0};

  parLimLo[0] = 1; // g0 const
  parLimHi[0] = 1000000;
  parLimLo[1] = par[1] - 0.5 * gausFunc->GetParameter(2); // g0 mean
  parLimHi[1] = par[1] + 0.5 * gausFunc->GetParameter(2);

  for(int i = 2; i < nPars; ++i) // allow a 50% variation on the already fitted parameters
    {
      parLimLo[i] = par[i] - 0.5 * fabs(par[i]);
      parLimHi[i] = par[i] + 0.5 * fabs(par[i]);
    }

  const char* parNames[nPars] = {"ConstG0", "MeanG0", "Width", "MPV", "Area", "GSigma"};
  std::cout << "Start parameters and limits\n";
  for(int i = 0; i < nPars; ++i)
    std::cout << parNames[i] << "\t\t" << par[i] << "    " << parLimLo[i] << "   " << parLimHi[i] << " \n";
  std::cout << std::endl;

  double fitR1 = inHist->GetXaxis()->GetXmin();
  double fitR2 = inHist->GetXaxis()->GetXmax();

  TF1* gausLang = new TF1("gausLang", gausLangaufun, fitR1, fitR2, nPars);
  gausLang->SetNpx(1e4);
  gausLang->SetParameters(par);
  gausLang->SetParNames("ConstG0", "MeanG0", "Width", "MPV", "Area", "GSigma");
  for(int i = 0; i < nPars; ++i)
    gausLang->SetParLimits(i, parLimLo[i], parLimHi[i]);

  inHist->Fit(gausLang, "RL");

  return gausLang;
}

// the sigma of the convolution and the one of the noise are constrained to be the same, mean and sigma of the noise are fixed
TF1* gausLanGausFitFixGaus(TH1* inHist, double negSigmaFit, double posSigmaFit, double mean, double sigma) // gauss parameters (mean and sigma) from another histo
{
  TF1* langauFunc = lanGausFit(inHist, negSigmaFit, posSigmaFit);

  const int nPars = 6;
  double par[nPars] = {0};

  // set the gaus fit parameters
  par[0] = inHist->GetBinContent(inHist->FindBin(0)); // constant gets the value of the bin at 0
  par[1] = mean; // these 2 remain fixed
  par[5] = sigma;
  for(int i = 2; i < nPars - 1; ++i) par[i] = langauFunc->GetParameter(i - 2); // get parameters form langaus fit, except gaus sigma

  double parLimHi[nPars] = {0};
  double parLimLo[nPars] = {0};

  parLimLo[0] = 1; // g0 const
  parLimHi[0] = 1000000;
  parLimLo[1] = mean; // g0 mean
  parLimHi[1] = mean;
  parLimLo[5] = sigma; // g0 sigma
  parLimHi[5] = sigma;

  for(int i = 2; i < nPars - 1; ++i) // allow a 50% variation on the already fitted parameters, fix gaus sigma
    {
      parLimLo[i] = par[i] - 0.5 * fabs(par[i]);
      parLimHi[i] = par[i] + 0.5 * fabs(par[i]);
    }
/*
  const char* parNames[nPars] = {"ConstG0", "MeanG0", "Width", "MPV", "Area", "GSigma"};
  std::cout << "Start parameters and limits\n";
  for(int i = 0; i < nPars; ++i)
    std::cout << parNames[i] << "\t\t" << par[i] << "    " << parLimLo[i] << "   " << parLimHi[i] << " \n";
  std::cout << std::endl;
*/
  double fitR1 = inHist->GetXaxis()->GetXmin();
  double fitR2 = inHist->GetXaxis()->GetXmax();

  TF1* gausLang = new TF1("gausLang", gausLangaufun, fitR1, fitR2, nPars);
  gausLang->SetNpx(1e4);
  gausLang->SetParameters(par);
  gausLang->SetParNames("ConstG0", "MeanG0", "Width", "MPV", "Area", "GSigma");
  for(int i = 0; i < nPars; ++i)
    gausLang->SetParLimits(i, parLimLo[i], parLimHi[i]);

  inHist->Fit(gausLang, "RQL");

  return gausLang;
}

Double_t gausNoiseLangaufun(Double_t* x, Double_t* par) // a peak at 0 and a landau gaussian convolution
{
  Double_t gausPart = par[0] * TMath::Gaus(*x, par[1], par[2]);
  Double_t langauPart = langaufun(x, &par[3]); // par 0 to 2 belong to the gauss part

  return langauPart + gausPart;
}

// just fix the gaus parameters of the gaussian close to 0, the rest is free
TF1* gausLanGausFitFixGausNoise(TH1* inHist, double negSigmaFit, double posSigmaFit, double mean, double sigma) // gauss parameters (mean and sigma) from another histo
{
  TF1* langauFunc = lanGausFit(inHist, negSigmaFit, posSigmaFit);

  const int nPars = 7;
  double par[nPars] = {0};

  // set the gaus fit parameters
  par[0] = inHist->GetBinContent(inHist->FindBin(0)); // constant gets the value of the bin at 0
  par[1] = mean; // these 2 remain fixed
  par[2] = sigma;
  for(int i = 3; i < nPars; ++i) par[i] = langauFunc->GetParameter(i - 3); // get parameters form langaus fit, except gaus sigma

  double parLimHi[nPars] = {0};
  double parLimLo[nPars] = {0};

  parLimLo[0] = 1; // g0 const
  parLimHi[0] = 1000000;
  parLimLo[1] = mean; // g0 mean
  parLimHi[1] = mean;
  parLimLo[2] = sigma; // g0 sigma
  parLimHi[2] = sigma;

  for(int i = 3; i < nPars; ++i) // allow a 50% variation on the already fitted parameters, fix gaus sigma
    {
      parLimLo[i] = par[i] - 0.5 * fabs(par[i]);
      parLimHi[i] = par[i] + 0.5 * fabs(par[i]);
    }
/*
  const char* parNames[nPars] = {"ConstG0", "MeanG0", "SigmaG0", "Width", "MPV", "Area", "GSigma"};
  std::cout << "Start parameters and limits\n";
  for(int i = 0; i < nPars; ++i)
    std::cout << parNames[i] << "\t\t" << par[i] << "    " << parLimLo[i] << "   " << parLimHi[i] << " \n";
  std::cout << std::endl;
*/
  double fitR1 = inHist->GetXaxis()->GetXmin();
  double fitR2 = inHist->GetXaxis()->GetXmax();

  TF1* gausLang = new TF1("gausLang", gausNoiseLangaufun, fitR1, fitR2, nPars);
  gausLang->SetNpx(1e4);
  gausLang->SetParameters(par);
  gausLang->SetParNames("ConstG0", "MeanG0", "SigmaG0", "Width", "MPV", "Area", "GSigma");
  for(int i = 0; i < nPars; ++i)
    gausLang->SetParLimits(i, parLimLo[i], parLimHi[i]);

  inHist->Fit(gausLang, "RQL");

  return gausLang;
}


void gausmagic(TH1* thehisto)
{

      // now we can subtract the noise peak from the signal...
    TF1 * gausnoise = new TF1("gausnoise","gaus",-50,5.0);
    thehisto->Fit(gausnoise,"QR");
    gausnoise->SetRange(-100.0,100.0);
    thehisto->Add(gausnoise,-1);

    // edit the histo to set all bins <0 to 0... 500 bins 
    for (int ii=0;ii<500;ii++)
    {
      thehisto->SetBinContent(ii,0.0);
    }

    // also catch adc>0 negative bins
    for (int ii=500;ii<515;ii++)
    {
      double tempdouble = thehisto->GetBinContent(ii);
      if (tempdouble < 0)
      {
	thehisto->SetBinContent(ii,0.0);
      }
    }

    // the call to the fit:
    TF1 *fitsnr2 = lanGausFit(thehisto,3.0,7.0);

}
