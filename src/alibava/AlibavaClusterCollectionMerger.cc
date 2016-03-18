/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *  email:eda.yildirim@cern.ch
 * Edited by Dmitry Kisler
 *  (2016 DESY)
 *  email:dmitry.kisler@desy.de
 */
// personal includes ".h"
#include "ALIBAVA.h"
#include "AlibavaClusterCollectionMerger.h"
#include "AlibavaRunHeaderImpl.h"
#include "AlibavaEventImpl.h"
#include "AlibavaBaseProcessor.h"
// eutelescope includes ".h"
#include "EUTELESCOPE.h"
#include "EUTelEventImpl.h"
// marlin includes
#include "marlin/Global.h"
#include "marlin/Exceptions.h"
#include "marlin/Processor.h"
#include "marlin/DataSourceProcessor.h"
#include "marlin/ProcessorMgr.h"
// aida includes <.h>
#include <marlin/AIDAProcessor.h>
#include <AIDA/ITree.h>
// lcio includes
#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/TrackerPulseImpl.h>
#include <IMPL/LCGenericObjectImpl.h>
#include <IMPL/LCRunHeaderImpl.h>
#include <IMPL/LCEventImpl.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h> //? did find in v01-17-05/lcio/v02-04-03/src/cpp/src/UTIL ?
// system includes
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cassert>
#include <memory>
#include <stdlib.h>
#include <algorithm>
#include <vector>
//ROOT includes
#include "TH2D.h"
#include "TObject.h"

using namespace std;
using namespace marlin;
using namespace lcio;
using namespace alibava;
using namespace eutelescope;
//AlibavaClusterCollectionMerger processor initiation:
AlibavaClusterCollectionMerger::AlibavaClusterCollectionMerger ():
    //rootObject, where we store histos
    _rootObjectMap(),
    //processor's name as it's required in the steering-template
    DataSourceProcessor("AlibavaClusterCollectionMerger"),
    //telescope data file :
    //_telescope_lcReader(),
    //which telescope data file to be opened
    _telescopeFileName(ALIBAVA::NOTSET),
    //which collection of puslses has to be selected from the telescope data file
    _telescopePulseCollectionName(ALIBAVA::NOTSET),
    //which sparse collection has to be selected from the telescope data file
    _telescopeSparseCollectionName(ALIBAVA::NOTSET),
    // alibava data file :
    //_alibava_lcReader(),
    _alibavaFileName(ALIBAVA::NOTSET),
    _alibavaPulseCollectionName(ALIBAVA::NOTSET),
    _alibavaSparseCollectionName(ALIBAVA::NOTSET),
    // output
    _outputPulseCollectionName(ALIBAVA::NOTSET),
    _outputSparseCollectionName(ALIBAVA::NOTSET)
{
    //variable initialization
    _description = "Merges alibava and telescope cluster collections";

    ////////////////
    //   INPUT    //
    ////////////////

    registerProcessorParameter ("MaxEventNumber",
                                "How many events do we want to process?",
                                _MaxEventNum , int(10000));
    ///////////////
    // Telescope //
    ///////////////

    // file name - telescope *.slcio file to be used as the input data for this, merging, step.
    //it was produced on the prev. steps: telescope-converter+telescope-clustering
    registerProcessorParameter("InputTelescopeFileName", "This is the input file name that telescope cluster collections stored",
                               _telescopeFileName, string("runXXXXXX.slcio") );

    // pulse collection - which telescope pulse collection to be taken for merging
    registerInputCollection (LCIO::TRACKERPULSE, "TelescopeClusterPulseCollectionName",
                             "Name of the cluster pulse collection of telescope data",
                             _telescopePulseCollectionName, string("telescope_cluster_pulse") );

    // sparse collection - which telescope sparse collection to be taken for merging
    registerInputCollection (LCIO::TRACKERDATA, "TelescopeSparseClusterCollectionName",
                             "Name of the sparse cluster collection of telescope data",
                             _telescopeSparseCollectionName, string("telescope_sparse_cluster") );

    /////////////
    // Alibava //
    /////////////

    // file name - alibava *.slcio file to be used as the input data for this, merging, step.
    //it was produced on the prev. steps: converter+reco+clustering-1+clustering-2
    registerProcessorParameter("InputAlibavaFileName", "This is the input file name that alibava cluster collections stored",
                               _alibavaFileName, string("runXXXXXX.slcio") );

    // pulse collection - which alibava pulse collection to be taken for merging
    registerInputCollection (LCIO::TRACKERPULSE, "AlibavaClusterPulseCollectionName",
                             "Name of the cluster pulse collection of alibava data",
                             _alibavaPulseCollectionName, string("alibava_cluster_pulse") );

    // sparse collection - which alibava sparse collection to be taken for merging
    registerInputCollection (LCIO::TRACKERDATA, "AlibavaSparseClusterCollectionName",
                             "Name of the sparse cluster collection of alibava data",
                             _alibavaSparseCollectionName, string("alibava_sparse_cluster") );

    ////////////
    // Output //
    ////////////
    // pulse collection
    registerOutputCollection (LCIO::TRACKERPULSE, "OutputClusterPulseCollectionName",
                              "Name of the merged/output cluster pulse collection",
                              _outputPulseCollectionName, string("merged_cluster_pulse") );

    // sparse collection
    registerOutputCollection (LCIO::TRACKERDATA, "OutputSparseClusterCollectionName",
                              "Name of the merged/output sparse cluster collection. DO NOT Change this. This is hard coded in other  ",
                              _outputSparseCollectionName, string("original_zsdata") );
    // differents in event counting between the telescope and alibava data
    registerProcessorParameter ("EventIDDifference",
                                "AlibavaEventNumber - TelescopeEventNumber",
                                _eventIDDiff , int(0));

    //Fake Alibava missing coordinate
    registerProcessorParameter ("MaxFakeDistance",
                                "Max difference of the mean clusters pix coordinates between the tel. plane 2 and plane 3 to fake the missing DUT coordinate",
                                _maxFakeDist , double(20));
    registerProcessorParameter ("DUTSensetiveCoordinate",
                                "Which DUT coordinate is position sensetive (which axis the DUT strips are perpedicular to?)",
                                _sensetiveAxis , string("y"));

}
//the processor is being initiated here
AlibavaClusterCollectionMerger * AlibavaClusterCollectionMerger::newProcessor ()
{
    return new AlibavaClusterCollectionMerger;
}
//funtion to print initial parameters which were set in the steering-template, config file and runlist
void AlibavaClusterCollectionMerger::init ()
{
    printParameters ();
}
//funtion to read the data corresponding to a certain event number.
//this event numb. is being sent as the input to this funciton:
//event sum
int eventCounter = 0;
//numb of events which contain clusters on each, telescope+CMSPixRef+DUT, plane
int event_cluster_at_each_plane = 0;
//numb of events which contain clusters on telescope+CMSPixRef planes
int event_cluster_at_tel_CMSRef = 0;
//numb of events which contain clusters on telescope+DUT planes
int event_cluster_at_tel_DUT = 0;
//numb of events which contain clusters on telescope(only telescope, 6 planes) planes
int event_cluster_at_tel = 0;
//array of indicators showing if cluster were/n't found on each of 8 planes: 1 - found, 0 - not found
int plane[8];

////CORRELATION PART defininitions
//clusters distribution on each telescope plane: y and x directions arrays [plane_ID][y] and [plane_ID][x]
int tel_y[6][576];
int tel_x[6][1152];
//clusters distr on DUT: y pos. coordinates (if x is unsensetive axis) is the real distr., estimation from the tel planes 2 and 3 along x coordinate
//Attention: ali_x and ali_y sizes should be changed to 128 and 576 correspondingly if DUT would be oriented along the global Y axis
int ali_y[128];
int ali_x[1152];
//clusters distr on CMSPixRef
int ref_y[80];
int ref_x[56];
//arrays to store mean X and Y coordinates of the clusters' pixels of the tel. planes 2 and 3 to fake the missing coordinate on DUT(Alibava)
//Attention: I limit the number of tracks to 1000 (D. Kisler)
double x_id2[1000];
double x_id3[1000];
double y_id2[1000];
double y_id3[1000];
double x_ali[1000];
int trk_counter;
//let's estimate number of tracks passing through the DUT
int trk_id2;
int trk_id3;

void AlibavaClusterCollectionMerger::readDataSource(int /* numEvents */)
{
    //booking hists
    bookHistos();
    // open telescope file
    LCReader* telescope_lcReader = LCFactory::getInstance()->createLCReader();
    try
    {
        telescope_lcReader->open(_telescopeFileName);
    }catch( IOException& e )
    {
        streamlog_out ( ERROR1 ) << "Can't open the telescope file: " << e.what() << endl ;
    }
    // open alibava file
    LCReader* alibava_lcReader = LCFactory::getInstance()->createLCReader();
    try
    {
        alibava_lcReader->open( _alibavaFileName );
    }catch( IOException& e )
    {
        streamlog_out ( ERROR1 ) << "Can't open the alibava file: " << e.what() << endl ;
    }
    // we will copy alibava run header to as the header of output file.
    try {
        LCRunHeader* alibava_runHeader = alibava_lcReader->readNextRunHeader();
        ProcessorMgr::instance()->processRunHeader( alibava_runHeader ) ;
    } catch( IOException& e ){
        streamlog_out ( ERROR1 ) << "Can't access run header of the alibava file: " << e.what() << endl ;
    }
    EUTelEventImpl*  telescopeEvent;
    AlibavaEventImpl* alibavaEvent;
    //matching the event numbers between the alibava and telescope
    if (_eventIDDiff<0 ){
        for (int i=0; i<abs(_eventIDDiff); i++)
            telescopeEvent = static_cast<EUTelEventImpl*> ( telescope_lcReader->readNextEvent() );
    }
    else if (_eventIDDiff>0){
        for (int i=0; i<_eventIDDiff; i++)
            alibavaEvent = static_cast<AlibavaEventImpl*> ( alibava_lcReader->readNextEvent() );
    }
    //boolean to check if the colleciton is presented in the *.slcio file
    bool noCollectionFound = false;
    //reading a telescope event. CHECK IF THE COLLECTIONS PRESENTED
    while( ((telescopeEvent = static_cast<EUTelEventImpl*> ( telescope_lcReader->readNextEvent())) != 0 )
           && ((alibavaEvent = static_cast<AlibavaEventImpl*> (alibava_lcReader->readNextEvent())) != 0 )
           && eventCounter < _MaxEventNum)
    {
        noCollectionFound = false;
        //stop reading if the end of run's event was reached
        if (telescopeEvent->getEventType() == kEORE)
        {
            streamlog_out ( MESSAGE5 ) << "Reached EORE of telescope data"<< endl;
            break;
        }
        if ( alibavaEvent->getEventNumber() % 1000 == 0 ){streamlog_out ( MESSAGE4 ) << "Looping events "<<alibavaEvent->getEventNumber() << endl;}
        //initiation of the vectors to store the input collections' data
        LCCollectionVec * alibavaPulseColVec = 0;
        LCCollectionVec * alibavaSparseColVec = 0;
        LCCollectionVec * telescopePulseColVec = 0;
        LCCollectionVec * telescopeSparseColVec = 0;
        // alibava collections are being treated in this try-catch
        try
        {
            //why we take sparse and pulse collections? do we need them both? (look to the manual)
            // get alibava collections
            alibavaPulseColVec = dynamic_cast< LCCollectionVec * > ( alibavaEvent->getCollection( _alibavaPulseCollectionName ) ) ;
            alibavaSparseColVec = dynamic_cast< LCCollectionVec * > ( alibavaEvent->getCollection( _alibavaSparseCollectionName ) ) ;

        }catch ( DataNotAvailableException& e )
        {
            noCollectionFound = true;
            streamlog_out( WARNING5 ) <<"No input collection " << _alibavaPulseCollectionName << " or "<< _alibavaSparseCollectionName << " found on alibava event " << alibavaEvent->getEventNumber() << endl;
        }
        // telescope collections are being treated in this try-catch
        try
        {
            // get telescope collections
            telescopePulseColVec = dynamic_cast< LCCollectionVec * > ( telescopeEvent->getCollection( _telescopePulseCollectionName ) ) ;
            telescopeSparseColVec = dynamic_cast< LCCollectionVec * > ( telescopeEvent->getCollection( _telescopeSparseCollectionName ) ) ;
        } catch ( DataNotAvailableException& e ) {
            noCollectionFound = true;
            streamlog_out( WARNING5 ) <<"No input collection " << _telescopePulseCollectionName << " or "<< _telescopeSparseCollectionName << " found on telescope event " << telescopeEvent->getEventNumber() << endl;
        }
        // create the vectors to store the data to be written as output afterwards
        LCCollectionVec * outputPulseColVec = new LCCollectionVec(LCIO::TRACKERPULSE);
        LCCollectionVec * outputSparseColVec = new LCCollectionVec(LCIO::TRACKERDATA);

        if ( !noCollectionFound )
        {
            //array to check is the event contains clusters at each plane:
            for(int i = 0; i < 8; i++){plane[i] = 0;}
            //clear the arreys for correlation
            memset(tel_x,0,sizeof(tel_x));
            memset(tel_y,0,sizeof(tel_y));
            memset(ali_x,0,sizeof(ali_x));
            memset(ali_y,0,sizeof(ali_y));
            memset(ref_x,0,sizeof(ref_x));
            memset(ref_y,0,sizeof(ref_y));
            memset(x_ali,0,sizeof(x_ali));
            //av. x coordinate on the tel planes 2 and 3 to estimate the x coordinate on DUT plane
            memset(x_id2, 0, sizeof(x_id2));
            memset(x_id3, 0, sizeof(x_id3));
            memset(y_id2, 0, sizeof(y_id2));
            memset(y_id3, 0, sizeof(y_id3));
            //let's estimate number of tracks passing through the DUT
            trk_id2 = 0;
            trk_id3 = 0;
            trk_counter = 0;

            // copy the telescope collections, Pusle and Sparcse, to the output collection
            copyClustersInCollection(outputPulseColVec, outputSparseColVec, telescopePulseColVec, telescopeSparseColVec);
            //check if we have clusters only on tel planes
            if((plane[0]+plane[1]+plane[2]+plane[3]+plane[4]+plane[5]) == 6)
            {
                event_cluster_at_tel++;
            }
            //check if we have clusters on tel+CMSPixRef planes
            if((plane[0]+plane[1]+plane[2]+plane[3]+plane[4]+plane[5]+plane[7]) == 7)
            {
                event_cluster_at_tel_CMSRef++;
            }
            // copy the alibava collections, Pusle and Sparcse, to the output collection
            copyClustersInCollection(outputPulseColVec, outputSparseColVec, alibavaPulseColVec, alibavaSparseColVec);
            //check if we have clusters on tel+DUT planes
            if((plane[0]+plane[1]+plane[2]+plane[3]+plane[4]+plane[5]+plane[6]) == 7)
            {
                event_cluster_at_tel_DUT++;
            }
            //check if we have clusters on each, tel+DUT+CMSPixRef, plane
            if((plane[0]+plane[1]+plane[2]+plane[3]+plane[4]+plane[5]+plane[6]+plane[7]) == 8)
            {
                event_cluster_at_each_plane++;
            }
            //fill the correlation histos
            fillingHistos();

            //the name and alibava event paraters, e.g. size, time, are being written to the output collections here
            try
            {
                // the pointer-vector to store the data from this events which will be written to the output collection it defined here
                AlibavaEventImpl* outputEvent = new AlibavaEventImpl();
                // the vector defined in the prev. line if being filled below
                outputEvent->setRunNumber( alibavaEvent->getRunNumber() );
                outputEvent->setEventNumber(eventCounter);
                outputEvent->setEventType( alibavaEvent->getEventType() );
                outputEvent->setEventSize( alibavaEvent->getEventSize() );
                outputEvent->setEventValue( alibavaEvent->getEventValue() );
                outputEvent->setEventTime( alibavaEvent->getEventTime() );
                outputEvent->setEventTemp( alibavaEvent->getEventTemp() );
                outputEvent->setCalCharge( alibavaEvent->getCalCharge() );
                outputEvent->setCalDelay( alibavaEvent->getCalDelay() );

                if (alibavaEvent->isEventMasked()) // the condition is being checked in AlibavaEventImpl.h
                    //masking event, see AlibavaEventImpl.cc/.h files for details
                    outputEvent->maskEvent();
                else
                    //unmasking event, see AlibavaEventImpl.cc/.h files for details
                    outputEvent->unmaskEvent();

                if ( !noCollectionFound ) {
                    outputEvent->addCollection(outputPulseColVec, _outputPulseCollectionName);
                    outputEvent->addCollection(outputSparseColVec, _outputSparseCollectionName);
                }
                ProcessorMgr::instance()->processEvent( static_cast<LCEventImpl*> ( outputEvent ) ) ;
                // delete outputEvent;
                //streamlog_out ( MESSAGE1 ) << "Successfully copied Alibava collections to output event" << endl ;
            } catch ( IOException& e)
            {
                streamlog_out( ERROR5 ) << e.what() << endl;
            }
            //add 1 to the eventCounter to jumt to next event
            eventCounter++;
        }// end of the loop over events

    }
}
// functoin to copy the input of both, telescope and alibava, collections to the output collection
void AlibavaClusterCollectionMerger::copyClustersInCollection(LCCollectionVec * outputPulseColVec, LCCollectionVec * outputSparseColVec, LCCollectionVec * inputPulseColVec, LCCollectionVec * inputSparseColVec)
{
    // Here is the Cell ID Encodes for pulseFrame and sparseFrame
    // CellID Encodes are introduced in eutelescope::EUTELESCOPE
    // for sparseFrame (usually called cluster collection)
    CellIDEncoder<TrackerDataImpl> outputSparseColEncoder ( eutelescope::EUTELESCOPE::ZSCLUSTERDEFAULTENCODING, outputSparseColVec );
    // for pulseFrame
    CellIDEncoder<TrackerPulseImpl> outputPulseColEncoder ( eutelescope::EUTELESCOPE::PULSEDEFAULTENCODING, outputPulseColVec );
    // for input
    CellIDDecoder<TrackerPulseImpl> inputPulseColDecoder(inputPulseColVec);
    CellIDDecoder<TrackerDataImpl> inputSparseColDecoder(inputSparseColVec);

    // how many clusters does telescope event have, declare the variable to store that
    unsigned int noOfClusters;
    // go through input clusters and copy them to output cluster collection
    noOfClusters = inputPulseColVec->getNumberOfElements();
    int ID;
    //calculating the
    for ( size_t i = 0; i < noOfClusters; ++i )
    {
        TrackerPulseImpl * outputPulseFrame = new TrackerPulseImpl();
        TrackerDataImpl * outputSparseFrame = new TrackerDataImpl();

        TrackerPulseImpl* inputPulseFrame = dynamic_cast<TrackerPulseImpl*>(inputPulseColVec->getElementAt(i));
        TrackerDataImpl* inputSparseFrame = dynamic_cast<TrackerDataImpl*>(inputPulseFrame->getTrackerData());
        ID = static_cast<int> (inputSparseColDecoder(inputSparseFrame)["sensorID"]) ;

        for(int counter = 0; counter < (int)(inputSparseFrame->getChargeValues().size()/4); counter ++ )
        {
            if(ID < 6) //(only) tel planes
            {
                plane[ID] = 1;
                tel_x[ID][(int)(inputSparseFrame->getChargeValues().at(0 + (int)(counter*4)))] ++;
                tel_y[ID][(int)(inputSparseFrame->getChargeValues().at(1 + (int)(counter*4)))] ++;
            }
            //mean X and Y on the tel. planes 2 and 3 for faking the missing DUT coordinate
            if(ID == 2)
            {
                x_id2[trk_id2] += (inputSparseFrame->getChargeValues().at(0 + (int)(counter*4)))/((int)(inputSparseFrame->getChargeValues().size()/4));
                y_id2[trk_id2] += (inputSparseFrame->getChargeValues().at(1 + (int)(counter*4)))/((int)(inputSparseFrame->getChargeValues().size()/4));
            }
            if(ID == 3)
            {
                x_id3[trk_id3] += (inputSparseFrame->getChargeValues().at(0 + (int)(counter*4)))/((int)(inputSparseFrame->getChargeValues().size()/4));
                y_id3[trk_id3] += (inputSparseFrame->getChargeValues().at(1 + (int)(counter*4)))/((int)(inputSparseFrame->getChargeValues().size()/4));
            }
            if(ID == 8) //cmspix
            {
                plane[7] = 1;
                ref_x[(int)(inputSparseFrame->getChargeValues().at(0 + (int)(counter*4)))] ++;
                ref_y[(int)(inputSparseFrame->getChargeValues().at(1 + (int)(counter*4)))] ++;
            }
            if(ID == 7) //DUT
            {
                plane[6] = 1;
                if(_sensetiveAxis == "y")
                {
                    ali_y[(int)(inputSparseFrame->getChargeValues().at(1 + (int)(counter*4)))] ++;
                }
                if(_sensetiveAxis == "x")
                {
                    ali_x[(int)(inputSparseFrame->getChargeValues().at(0 + (int)(counter*4)))] ++;
                }
            }
        }

        if(ID == 2) {trk_id2 ++;}
        if(ID == 3) {trk_id3 ++;}

        //estimating the DUT missing coordinate
        if(ID == 7)
        {
            streamlog_out (DEBUG0) << "numbOfClus = " << noOfClusters <<endl;

            if (i == 0)
            {
                double dist = 0;
                if(_sensetiveAxis == "y")
                {
                    for(int i = 0; i < trk_id2; i++)
                    {
                        for(int j = 0; j < trk_id3; j++)
                        {
                            dist = abs(sqrt(x_id2[i]*x_id2[i] + y_id2[i]*y_id2[i]) - sqrt(x_id3[j]*x_id3[j] + y_id3[j]*y_id3[j]));
                            //if the distance between the hits on id2 and id3 is less than the set limit, estimate the hit coordinate on DUT
                            if( dist < _maxFakeDist)
                            {
                                streamlog_out (DEBUG0) << "faking check. in the loop" << endl;
                                ali_x[(int)((x_id2[i] + x_id3[j])/2)] ++;
                                x_ali[trk_counter] = (int)((x_id2[i] + x_id3[j])/2);
                                streamlog_out (DEBUG0) << "checking the ali missing coordinate estimation" << endl
                                                       << "event = " << eventCounter << endl<<
                                                          "tel. ID2: "<< " trk_id2 = " <<trk_id2 << " : trk# = " << i << " av. x_id2 = " << x_id2[i] <<endl
                                                       << "tel. ID3: "<< " trk_id3 = " << trk_id3 << " : trk# = " << j << " av. x_id3 = " << x_id3[j] << endl
                                                       << "estimationDist = " << dist << " < " << _maxFakeDist << "(estimation limitDist)"
                                                       <<"; tot_trk# = " << trk_counter << " | x = " << (int)((x_id2[i] + x_id3[j])/2) <<
                                                         " | x (written) = " << x_ali[trk_counter] <<endl;
                                trk_counter++ ;
                            }
                        }
                    }
                }

                if(_sensetiveAxis == "x")
                {
                    for(int i = 0; i < trk_id2; i++)
                    {
                        for(int j = 0; j < trk_id3; j++)
                        {
                            dist = abs(sqrt(x_id2[i]*x_id2[i] + y_id2[i]*y_id2[i]) - sqrt(x_id3[j]*x_id3[j] + y_id3[j]*y_id3[j]));
                            if(dist < _maxFakeDist)
                            {
                                ali_y[(int)((y_id2[i] + y_id3[j])/2)] ++;
                                x_ali[trk_counter] = (int)((y_id2[i] + y_id3[j])/2);
                                streamlog_out (DEBUG0) << "checking the ali missing coordinate estimation" << endl
                                                       << "event = " << eventCounter << endl<<
                                                          "tel. ID2: "<< " trk_id2 = " <<trk_id2 << " : trk# = " << i << " av. y_id2 = " << y_id2[i] << endl
                                                       << "tel. ID3: "<< " trk_id3 = " << trk_id3 << " : trk# = " << j << " av. y_id3 = " << y_id3[i] <<endl
                                                       << "estimationDist = " << dist << " < " << _maxFakeDist << "(estimation limitDist)"
                                                       <<"; tot_trk# = " << trk_counter << " | x = " << (int)((y_id2[i] + y_id3[j])/2) <<
                                                         " | x (written) = " << x_ali[trk_counter] <<endl;
                                trk_counter++ ;
                            }
                        }
                    }
                }
                streamlog_out (DEBUG0) << "tot numb. of trk = " <<trk_counter << endl;
            }
        }

        //writing data to the output collection
        //sparse collection
        //if not DUT:
        //hardcoded: temporarly switched off separate writing of the DUT collections
        if(ID != 100)
        {
            // set Cell ID for sparse collection
            outputSparseColEncoder["sensorID"] = static_cast<int>(inputSparseColDecoder(inputSparseFrame) ["sensorID"]);
            outputSparseColEncoder["sparsePixelType"] =static_cast<int>(inputSparseColDecoder(inputSparseFrame)["sparsePixelType"]);
            outputSparseColEncoder["quality"] = static_cast<int>(inputSparseColDecoder(inputSparseFrame)["quality"]);
            outputSparseColEncoder.setCellID( outputSparseFrame );
            // copy tracker data
            outputSparseFrame->setChargeValues(inputSparseFrame->getChargeValues());
            // add it to the cluster collection
            outputSparseColVec->push_back( outputSparseFrame );
        }
        //if DUT:
//        if (ID == 7)
//        {
//            //final vector for setChargeValues method while forming the output sparse collection
//            const int size = inputSparseFrame->getChargeValues().size();
//            float aliVec[size];
//            //            EVENT::FloatVec Vec = new EVENT::FloatVec;

//            for(int itrk = 0; itrk < trk_counter; itrk++)
//            {
//                // set Cell ID for sparse collection
//                outputSparseColEncoder["sensorID"] = static_cast<int>(inputSparseColDecoder(inputSparseFrame) ["sensorID"]);
//                outputSparseColEncoder["sparsePixelType"] =static_cast<int>(inputSparseColDecoder(inputSparseFrame)["sparsePixelType"]);
//                outputSparseColEncoder["quality"] = static_cast<int>(inputSparseColDecoder(inputSparseFrame)["quality"]);
//                outputSparseColEncoder.setCellID( outputSparseFrame );
//                if(_sensetiveAxis == "y")
//                {
//                    for(int counter = 0; counter < (int)(size/4); counter ++ )
//                    {
//                        aliVec[0 + (int)(counter*4)] = x_ali[itrk];
//                        aliVec[1 + (int)(counter*4)] = inputSparseFrame->getChargeValues().at(1 + (int)(counter*4));
//                        aliVec[2 + (int)(counter*4)] = inputSparseFrame->getChargeValues().at(2 + (int)(counter*4));
//                        aliVec[3 + (int)(counter*4)] = inputSparseFrame->getChargeValues().at(3 + (int)(counter*4));
//                    }
//                }
//                if(_sensetiveAxis == "x")
//                {
//                    for(int counter = 0; counter < (int)(size/4); counter ++ )
//                    {
//                        aliVec[0 + (int)(counter*4)] = inputSparseFrame->getChargeValues().at(0 + (int)(counter*4));
//                        aliVec[1 + (int)(counter*4)] = x_ali[itrk];
//                        aliVec[2 + (int)(counter*4)] = inputSparseFrame->getChargeValues().at(2 + (int)(counter*4));
//                        aliVec[3 + (int)(counter*4)] = inputSparseFrame->getChargeValues().at(3 + (int)(counter*4));
//                    }
//                }
//                //                Vec-> push_back(aliVec);
//                //                outputSparseFrame->setChargeValues(Vec); //have to convert aliVec to the type required by setChargeValues function; look in TrackerDataImpl.cc
//                outputSparseFrame->setChargeValues(inputSparseFrame->getChargeValues());
//                outputSparseColVec->push_back( outputSparseFrame );
//            }
//        }

        //pulse collection
        // prepare a pulse for this cluster
        outputPulseColEncoder["sensorID"] = static_cast<int> (inputPulseColDecoder(inputPulseFrame) ["sensorID"]);
        outputPulseColEncoder["type"] = static_cast<int>(inputPulseColDecoder(inputPulseFrame) ["type"]);
        outputPulseColEncoder.setCellID( outputPulseFrame );

        outputPulseFrame->setCharge( inputPulseFrame->getCharge() );
        outputPulseFrame->setTrackerData( outputSparseFrame);
        outputPulseColVec->push_back( outputPulseFrame );

    } // end of loop over input clusters
}
//fillin in the histos with the correlations, DUT vs. tel. planes and CMSPixRef
void AlibavaClusterCollectionMerger::fillingHistos()
{
    if(_sensetiveAxis == "y")
    {
        //along Y axis
        for(int y = 0; y < 128; y++)
        {
            if(ali_y[y] > 0)
            {
                //correlation DUT vs. each tel. plane
                for(int id = 0; id < 6; id++)
                {
                    for(int yy = 0; yy < 576; yy++)
                    {
                        if(tel_y[id][yy] > 0)
                        {
                            if(id == 0)
                            {TH2D * hCorr_ali_0 = dynamic_cast<TH2D*> (_rootObjectMap[getHistoName("hCorrY_ali_vs_", id)]) ;
                                // hCorr_ali_0->Fill(y, yy, tel_y[id][yy]);
                                hCorr_ali_0->Fill(y, yy);
                            }
                            if(id == 1)
                            {TH2D * hCorr_ali_1 = dynamic_cast<TH2D*> (_rootObjectMap[getHistoName("hCorrY_ali_vs_", id)]) ;
                                //hCorr_ali_1->Fill(y, yy, tel_y[id][yy]);
                                hCorr_ali_1->Fill(y, yy);
                            }
                            if(id == 2)
                            {TH2D * hCorr_ali_2 = dynamic_cast<TH2D*> (_rootObjectMap[getHistoName("hCorrY_ali_vs_", id)]) ;
                                //hCorr_ali_2->Fill(y, yy, tel_y[id][yy]);
                                hCorr_ali_2->Fill(y, yy);
                            }
                            if(id == 3)
                            {TH2D * hCorr_ali_3 = dynamic_cast<TH2D*> (_rootObjectMap[getHistoName("hCorrY_ali_vs_", id)]) ;
                                //hCorr_ali_3->Fill(y, yy, tel_y[id][yy]);
                                hCorr_ali_3->Fill(y, yy);
                            }
                            if(id == 4)
                            {TH2D * hCorr_ali_4 = dynamic_cast<TH2D*> (_rootObjectMap[getHistoName("hCorrY_ali_vs_", id)]) ;
                                // hCorr_ali_4->Fill(y, yy, tel_y[id][yy]);
                                hCorr_ali_4->Fill(y, yy);
                            }
                            if(id == 5)
                            {TH2D * hCorr_ali_5 = dynamic_cast<TH2D*> (_rootObjectMap[getHistoName("hCorrY_ali_vs_", id)]) ;
                                //hCorr_ali_5->Fill(y, yy, tel_y[id][yy]);}
                                hCorr_ali_5->Fill(y, yy);}
                        }
                    }
                }
                //correlation DUT vs. CMSPixRef
                for(int yy = 0; yy < 80; yy++)
                {
                    if(ref_y[yy] > 0)
                    {
                        TH2D * hCorr_ali_6 = dynamic_cast<TH2D*> (_rootObjectMap[getHistoName("hCorrY_ali_vs_", 6)]) ;
                        hCorr_ali_6->Fill(y, yy); //, ref_y[yy]);
                        TH2D * hCorr_Y_ref_6 = dynamic_cast<TH2D*> (_rootObjectMap[getHistoName("hCorrY_CMSPixRef_vs_", 6)]) ;
                        hCorr_Y_ref_6->Fill(yy, y, ref_y[yy]);
                    }
                }
            }
        }

        //along X axis
        for(int x = 0; x < 1152; x++)
        {
            if(ali_x[x] > 0)
            {
                //correlation DUT vs. each tel. plane
                for(int id = 0; id < 6; id++)
                {
                    for(int xx = 0; xx < 1152; xx++)
                    {
                        if(tel_x[id][xx] > 0)
                        {
                            if(id == 0)
                            {TH2D * hCorr_X_ali_0 = dynamic_cast<TH2D*> (_rootObjectMap[getHistoName("hCorrX_ali_vs_", id)]) ;
                                hCorr_X_ali_0->Fill(x, xx);
                            }
                            if(id == 1)
                            {TH2D * hCorr_X_ali_1 = dynamic_cast<TH2D*> (_rootObjectMap[getHistoName("hCorrX_ali_vs_", id)]) ;
                                hCorr_X_ali_1->Fill(x, xx);
                            }
                            if(id == 2)
                            {TH2D * hCorr_X_ali_2 = dynamic_cast<TH2D*> (_rootObjectMap[getHistoName("hCorrX_ali_vs_", id)]) ;
                                hCorr_X_ali_2->Fill(x, xx);
                            }
                            if(id == 3)
                            {TH2D * hCorr_X_ali_3 = dynamic_cast<TH2D*> (_rootObjectMap[getHistoName("hCorrX_ali_vs_", id)]) ;
                                hCorr_X_ali_3->Fill(x, xx);
                            }
                            if(id == 4)
                            {TH2D * hCorr_X_ali_4 = dynamic_cast<TH2D*> (_rootObjectMap[getHistoName("hCorrX_ali_vs_", id)]) ;
                                hCorr_X_ali_4->Fill(x, xx);
                            }
                            if(id == 5)
                            {TH2D * hCorr_X_ali_5 = dynamic_cast<TH2D*> (_rootObjectMap[getHistoName("hCorrX_ali_vs_", id)]) ;
                                hCorr_X_ali_5->Fill(x, xx);
                            }
                        }
                    }
                }
                //correlation DUT vs. CMSPixRef
                for(int xx = 0; xx < 52; xx++)
                {
                    if(ref_x[xx] > 0)
                    {
                        TH2D * hCorr_X_ali_6 = dynamic_cast<TH2D*> (_rootObjectMap[getHistoName("hCorrX_ali_vs_", 6)]) ;
                        hCorr_X_ali_6->Fill(x, xx); //, ref_y[yy]);
                        TH2D * hCorr_X_ref_6 = dynamic_cast<TH2D*> (_rootObjectMap[getHistoName("hCorrX_CMSPixRef_vs_", 6)]) ;
                        hCorr_X_ref_6->Fill(xx, x, ref_x[xx]);
                    }
                }
            }
        }

    }
    //correlations between CMSPixRef and the tel. planes
    //along Y axis
    for(int y = 0; y < 80; y++)
    {
        if(ref_y[y] > 0)
        {
            //correlation CMSPixRef vs. each tel. plane
            for(int id = 0; id < 6; id++)
            {
                for(int yy = 0; yy < 576; yy++)
                {
                    if(tel_y[id][yy] > 0)
                    {
                        if(id == 0)
                        {TH2D * hCorr_Y_ref_0 = dynamic_cast<TH2D*> (_rootObjectMap[getHistoName("hCorrY_CMSPixRef_vs_", id)]) ;
                            hCorr_Y_ref_0->Fill(y, yy, tel_y[id][yy]);}
                        if(id == 1)
                        {TH2D * hCorr_Y_ref_1 = dynamic_cast<TH2D*> (_rootObjectMap[getHistoName("hCorrY_CMSPixRef_vs_", id)]) ;
                            hCorr_Y_ref_1->Fill(y, yy, tel_y[id][yy]);}
                        if(id == 2)
                        {TH2D * hCorr_Y_ref_2 = dynamic_cast<TH2D*> (_rootObjectMap[getHistoName("hCorrY_CMSPixRef_vs_", id)]) ;
                            hCorr_Y_ref_2->Fill(y, yy, tel_y[id][yy]);}
                        if(id == 3)
                        {TH2D * hCorr_Y_ref_3 = dynamic_cast<TH2D*> (_rootObjectMap[getHistoName("hCorrY_CMSPixRef_vs_", id)]) ;
                            hCorr_Y_ref_3->Fill(y, yy, tel_y[id][yy]);}
                        if(id == 4)
                        {TH2D * hCorr_Y_ref_4 = dynamic_cast<TH2D*> (_rootObjectMap[getHistoName("hCorrY_CMSPixRef_vs_", id)]) ;
                            hCorr_Y_ref_4->Fill(y, yy, tel_y[id][yy]);}
                        if(id == 5)
                        {TH2D * hCorr_Y_ref_5 = dynamic_cast<TH2D*> (_rootObjectMap[getHistoName("hCorrY_CMSPixRef_vs_", id)]) ;
                            hCorr_Y_ref_5->Fill(y, yy, tel_y[id][yy]);}
                    }
                }
            }
        }
    }
    //along X axis
    for(int x = 0; x < 52; x++)
    {
        if(ref_x[x] > 0)
        {
            //correlation CMSPixRef vs. each tel. plane
            for(int id = 0; id < 6; id++)
            {
                for(int xx = 0; xx < 1152; xx++)
                {
                    if(tel_x[id][xx] > 0)
                    {
                        if(id == 0)
                        {TH2D * hCorr_X_ref_0 = dynamic_cast<TH2D*> (_rootObjectMap[getHistoName("hCorrX_CMSPixRef_vs_", id)]) ;
                            hCorr_X_ref_0->Fill(x, xx, tel_x[id][xx]);}
                        if(id == 1)
                        {TH2D * hCorr_X_ref_1 = dynamic_cast<TH2D*> (_rootObjectMap[getHistoName("hCorrX_CMSPixRef_vs_", id)]) ;
                            hCorr_X_ref_1->Fill(x, xx, tel_x[id][xx]);}
                        if(id == 2)
                        {TH2D * hCorr_X_ref_2 = dynamic_cast<TH2D*> (_rootObjectMap[getHistoName("hCorrX_CMSPixRef_vs_", id)]) ;
                            hCorr_X_ref_2->Fill(x, xx, tel_x[id][xx]);}
                        if(id == 3)
                        {TH2D * hCorr_X_ref_3 = dynamic_cast<TH2D*> (_rootObjectMap[getHistoName("hCorrX_CMSPixRef_vs_", id)]) ;
                            hCorr_X_ref_3->Fill(x, xx, tel_x[id][xx]);}
                        if(id == 4)
                        {TH2D * hCorr_X_ref_4 = dynamic_cast<TH2D*> (_rootObjectMap[getHistoName("hCorrX_CMSPixRef_vs_", id)]) ;
                            hCorr_X_ref_4->Fill(x, xx, tel_x[id][xx]);}
                        if(id == 5)
                        {TH2D * hCorr_X_ref_5 = dynamic_cast<TH2D*> (_rootObjectMap[getHistoName("hCorrX_CMSPixRef_vs_", id)]) ;
                            hCorr_X_ref_5->Fill(x, xx, tel_x[id][xx]);}
                    }
                }
            }
        }
    }
}

//creating the hists
void AlibavaClusterCollectionMerger::bookHistos()
{
    //temp string variables for hists names and titles
    string histoName;
    string title;

    ////correlation of the DUT clusters with other planes' clusters along Y coordinate
    AIDAProcessor::tree(this)->cd(this->name());
    string basePath = "CorrY";
    AIDAProcessor::tree(this)->mkdir(basePath);
    basePath.append("/");
    basePath = "CorrY/Alibava";
    AIDAProcessor::tree(this)->mkdir(basePath);
    basePath.append("/");
    AIDAProcessor::tree(this)->cd(basePath);
    //id0
    histoName = getHistoName("hCorrY_ali_vs_",0);
    TH2D * hCorr_ali_0 = new TH2D (histoName.c_str(),"", 128, 0, 127, 576, 0, 575);
    title = "Clusters Y correlation: tel. plane 0 vs. DUT;Alibava channel number; Tel. plane Y coordinate";
    hCorr_ali_0->SetTitle(title.c_str());
    _rootObjectMap.insert(make_pair(histoName, hCorr_ali_0));
    //id1
    histoName = getHistoName("hCorrY_ali_vs_",1);
    TH2D * hCorr_ali_1 = new TH2D (histoName.c_str(),"", 128, 0, 127, 576, 0, 575);
    title = "Clusters Y correlation: tel. plane 1 vs. DUT;Alibava channel number; Tel. plane Y coordinate";
    hCorr_ali_1->SetTitle(title.c_str());
    _rootObjectMap.insert(make_pair(histoName, hCorr_ali_1));
    //id2
    histoName = getHistoName("hCorrY_ali_vs_",2);
    TH2D * hCorr_ali_2 = new TH2D (histoName.c_str(),"", 128, 0, 127, 576, 0, 575);
    title = "Clusters Y correlation: tel. plane 2 vs. DUT;Alibava channel number; Tel. plane Y coordinate";
    hCorr_ali_2->SetTitle(title.c_str());
    _rootObjectMap.insert(make_pair(histoName, hCorr_ali_2));
    //id3
    histoName = getHistoName("hCorrY_ali_vs_",3);
    TH2D * hCorr_ali_3 = new TH2D (histoName.c_str(),"", 128, 0, 127, 576, 0, 575);
    title = "Clusters Y correlation: tel. plane 3 vs. DUT;Alibava channel number; Tel. plane Y coordinate";
    hCorr_ali_3->SetTitle(title.c_str());
    _rootObjectMap.insert(make_pair(histoName, hCorr_ali_3));
    //id4
    histoName = getHistoName("hCorrY_ali_vs_",4);
    TH2D * hCorr_ali_4 = new TH2D (histoName.c_str(),"", 128, 0, 127, 576, 0, 575);
    title = "Clusters Y correlation: tel. plane 4 vs. DUT;Alibava channel number; Tel. plane Y coordinate";
    hCorr_ali_4->SetTitle(title.c_str());
    _rootObjectMap.insert(make_pair(histoName, hCorr_ali_4));
    //id5
    histoName = getHistoName("hCorrY_ali_vs_",5);
    TH2D * hCorr_ali_5 = new TH2D (histoName.c_str(),"", 128, 0, 127, 576, 0, 575);
    title = "Clusters Y correlation: tel. plane 5 vs. DUT;Alibava channel number; Tel. plane Y coordinate";
    hCorr_ali_5->SetTitle(title.c_str());
    _rootObjectMap.insert(make_pair(histoName, hCorr_ali_5));
    //CMSPixRef
    histoName = getHistoName("hCorrY_ali_vs_",6);
    TH2D * hCorr_ali_6 = new TH2D (histoName.c_str(),"", 128, 0, 127, 80, 0, 79);
    title = "Clusters Y correlation: CMSPixRef vs. DUT;Alibava channel number; CMSPixRef plane Y coordinate";
    hCorr_ali_6->SetTitle(title.c_str());
    _rootObjectMap.insert(make_pair(histoName, hCorr_ali_6));
    //along Y axis
    AIDAProcessor::tree(this)->cd(this->name());
    basePath = "CorrX";
    AIDAProcessor::tree(this)->mkdir(basePath);
    basePath.append("/");
    basePath = "CorrX/Alibava";
    AIDAProcessor::tree(this)->mkdir(basePath);
    basePath.append("/");
    AIDAProcessor::tree(this)->cd(basePath);
    //id0
    histoName = getHistoName("hCorrX_ali_vs_",0);
    TH2D * hCorr_X_ali_0 = new TH2D (histoName.c_str(),"", 1152, 0, 1151, 1152, 0, 1151);
    title = "Clusters X correlation: tel. plane 0 vs. DUT; Alibava 'reconstructed' coordinate; Tel. plane X coordinate";
    hCorr_X_ali_0->SetTitle(title.c_str());
    _rootObjectMap.insert(make_pair(histoName, hCorr_X_ali_0));
    //id1
    histoName = getHistoName("hCorrX_ali_vs_",1);
    TH2D * hCorr_X_ali_1 = new TH2D (histoName.c_str(),"", 1152, 0, 1151, 1152, 0, 1151);
    title = "Clusters X correlation: tel. plane 1 vs. DUT;Alibava 'reconstructed' coordinate; Tel. plane X coordinate";
    hCorr_X_ali_1->SetTitle(title.c_str());
    _rootObjectMap.insert(make_pair(histoName, hCorr_X_ali_1));
    //id2
    histoName = getHistoName("hCorrX_ali_vs_",2);
    TH2D * hCorr_X_ali_2 = new TH2D (histoName.c_str(),"", 1152, 0, 1151, 1152, 0, 1151);
    title = "Clusters X correlation: tel. plane 2 vs. DUT;Alibava 'reconstructed' coordinate; Tel. plane X coordinate";
    hCorr_X_ali_2->SetTitle(title.c_str());
    _rootObjectMap.insert(make_pair(histoName, hCorr_X_ali_2));
    //id3
    histoName = getHistoName("hCorrX_ali_vs_",3);
    TH2D * hCorr_X_ali_3 = new TH2D (histoName.c_str(),"", 1152, 0, 1151, 1152, 0, 1151);
    title = "Clusters X correlation: tel. plane 3 vs. DUT;Alibava 'reconstructed' coordinate; Tel. plane X coordinate";
    hCorr_X_ali_3->SetTitle(title.c_str());
    _rootObjectMap.insert(make_pair(histoName, hCorr_X_ali_3));
    //id4
    histoName = getHistoName("hCorrX_ali_vs_",4);
    TH2D * hCorr_X_ali_4 = new TH2D (histoName.c_str(),"", 1152, 0, 1151, 1152, 0, 1151);
    title = "Clusters X correlation: tel. plane 4 vs. DUT;Alibava 'reconstructed' coordinate; Tel. plane X coordinate";
    hCorr_X_ali_4->SetTitle(title.c_str());
    _rootObjectMap.insert(make_pair(histoName, hCorr_X_ali_4));
    //id5
    histoName = getHistoName("hCorrX_ali_vs_",5);
    TH2D * hCorr_X_ali_5 = new TH2D (histoName.c_str(),"", 1152, 0, 1151, 1152, 0, 1151);
    title = "Clusters X correlation: tel. plane 5 vs. DUT;Alibava 'reconstructed' coordinate; Tel. plane X coordinate";
    hCorr_X_ali_5->SetTitle(title.c_str());
    _rootObjectMap.insert(make_pair(histoName, hCorr_X_ali_5));
    //CMSPixRef
    histoName = getHistoName("hCorrX_ali_vs_",6);
    TH2D * hCorr_X_ali_6 = new TH2D (histoName.c_str(),"", 1152, 0, 1151, 52, 0, 51);
    title = "Clusters X correlation: CMSPixRef vs. DUT;Alibava 'reconstructed' coordinate; CMSPixRef plane X coordinate";
    hCorr_X_ali_6->SetTitle(title.c_str());
    _rootObjectMap.insert(make_pair(histoName, hCorr_X_ali_6));

    //Correlation of the clusters on CMSPixRef with the other planes' cluster coordinates
    //along X axis
    AIDAProcessor::tree(this)->cd(this->name());
    basePath = "CorrX";
    AIDAProcessor::tree(this)->mkdir(basePath);
    basePath.append("/");
    basePath = "CorrX/CMSPixRef";
    AIDAProcessor::tree(this)->mkdir(basePath);
    basePath.append("/");
    AIDAProcessor::tree(this)->cd(basePath);
    //id0
    histoName = getHistoName("hCorrX_CMSPixRef_vs_",0);
    TH2D * hCorr_X_ref_0 = new TH2D (histoName.c_str(),"", 52, 0, 51, 1152, 0, 1151);
    title = "Clusters X correlation: tel. plane 0 vs. CMSPixRef;CMSPixRef X coordinate; Tel. plane X coordinate";
    hCorr_X_ref_0->SetTitle(title.c_str());
    _rootObjectMap.insert(make_pair(histoName, hCorr_X_ref_0));
    //id1
    histoName = getHistoName("hCorrX_CMSPixRef_vs_",1);
    TH2D * hCorr_X_ref_1 = new TH2D (histoName.c_str(),"", 52, 0, 51, 1152, 0, 1151);
    title = "Clusters X correlation: tel. plane 1 vs. CMSPixRef;CMSPixRef X coordinate; Tel. plane X coordinate";
    hCorr_X_ref_1->SetTitle(title.c_str());
    _rootObjectMap.insert(make_pair(histoName, hCorr_X_ref_1));
    //id2
    histoName = getHistoName("hCorrX_CMSPixRef_vs_",2);
    TH2D * hCorr_X_ref_2 = new TH2D (histoName.c_str(),"", 52, 0, 51, 1152, 0, 1151);
    title = "Clusters X correlation: tel. plane 2 vs. CMSPixRef;CMSPixRef X coordinate; Tel. plane X coordinate";
    hCorr_X_ref_2->SetTitle(title.c_str());
    _rootObjectMap.insert(make_pair(histoName, hCorr_X_ref_2));
    //id3
    histoName = getHistoName("hCorrX_CMSPixRef_vs_",3);
    TH2D * hCorr_X_ref_3 = new TH2D (histoName.c_str(),"", 52, 0, 51, 1152, 0, 1151);
    title = "Clusters X correlation: tel. plane 3 vs. CMSPixRef;CMSPixRef X coordinate; Tel. plane X coordinate";
    hCorr_X_ref_3->SetTitle(title.c_str());
    _rootObjectMap.insert(make_pair(histoName, hCorr_X_ref_3));
    //id4
    histoName = getHistoName("hCorrX_CMSPixRef_vs_",4);
    TH2D * hCorr_X_ref_4 = new TH2D (histoName.c_str(),"", 52, 0, 51, 1152, 0, 1151);
    title = "Clusters X correlation: tel. plane 4 vs. CMSPixRef;CMSPixRef X coordinate; Tel. plane X coordinate";
    hCorr_X_ref_4->SetTitle(title.c_str());
    _rootObjectMap.insert(make_pair(histoName, hCorr_X_ref_4));
    //id5
    histoName = getHistoName("hCorrX_CMSPixRef_vs_",5);
    TH2D * hCorr_X_ref_5 = new TH2D (histoName.c_str(),"", 52, 0, 51, 1152, 0, 1151);
    title = "Clusters X correlation: tel. plane 5 vs. CMSPixRef;CMSPixRef X coordinate; Tel. plane X coordinate";
    hCorr_X_ref_5->SetTitle(title.c_str());
    _rootObjectMap.insert(make_pair(histoName, hCorr_X_ref_5));
    //DUT
    histoName = getHistoName("hCorrX_CMSPixRef_vs_",6);
    TH2D * hCorr_X_ref_6 = new TH2D (histoName.c_str(),"", 52, 0, 51, 1152, 0, 1151);
    title = "Clusters X correlation: DUT vs. CMSPixRef;CMSPixRef X coordinate; Alibava 'reconstructed' coordinate";
    hCorr_X_ref_6->SetTitle(title.c_str());
    _rootObjectMap.insert(make_pair(histoName, hCorr_X_ref_6));
    //along Y axis
    AIDAProcessor::tree(this)->cd(this->name());
    basePath = "CorrY";
    AIDAProcessor::tree(this)->mkdir(basePath);
    basePath.append("/");
    basePath = "CorrY/CMSPixRef";
    AIDAProcessor::tree(this)->mkdir(basePath);
    basePath.append("/");
    AIDAProcessor::tree(this)->cd(basePath);
    //id0
    histoName = getHistoName("hCorrY_CMSPixRef_vs_",0);
    TH2D * hCorr_Y_ref_0 = new TH2D (histoName.c_str(),"", 80, 0, 79, 576, 0, 575);
    title = "Clusters Y correlation: tel. plane 0 vs. CMSPixRef;CMSPixRef Y coordinate; Tel. plane Y coordinate";
    hCorr_Y_ref_0->SetTitle(title.c_str());
    _rootObjectMap.insert(make_pair(histoName, hCorr_Y_ref_0));
    //id1
    histoName = getHistoName("hCorrY_CMSPixRef_vs_",1);
    TH2D * hCorr_Y_ref_1 = new TH2D (histoName.c_str(),"", 80, 0, 79, 576, 0, 575);
    title = "Clusters Y correlation: tel. plane 1 vs. CMSPixRef;CMSPixRef Y coordinate; Tel. plane Y coordinate";
    hCorr_Y_ref_1->SetTitle(title.c_str());
    _rootObjectMap.insert(make_pair(histoName, hCorr_Y_ref_1));
    //id2
    histoName = getHistoName("hCorrY_CMSPixRef_vs_",2);
    TH2D * hCorr_Y_ref_2 = new TH2D (histoName.c_str(),"", 80, 0, 79, 576, 0, 575);
    title = "Clusters Y correlation: tel. plane 2 vs. CMSPixRef;CMSPixRef Y coordinate; Tel. plane Y coordinate";
    hCorr_Y_ref_2->SetTitle(title.c_str());
    _rootObjectMap.insert(make_pair(histoName, hCorr_Y_ref_2));
    //id3
    histoName = getHistoName("hCorrY_CMSPixRef_vs_",3);
    TH2D * hCorr_Y_ref_3 = new TH2D (histoName.c_str(),"", 80, 0, 79, 576, 0, 575);
    title = "Clusters Y correlation: tel. plane 3 vs. CMSPixRef;CMSPixRef Y coordinate; Tel. plane Y coordinate";
    hCorr_Y_ref_3->SetTitle(title.c_str());
    _rootObjectMap.insert(make_pair(histoName, hCorr_Y_ref_3));
    //id4
    histoName = getHistoName("hCorrY_CMSPixRef_vs_",4);
    TH2D * hCorr_Y_ref_4 = new TH2D (histoName.c_str(),"", 80, 0, 79, 576, 0, 575);
    title = "Clusters Y correlation: tel. plane 4 vs. CMSPixRef;CMSPixRef Y coordinate; Tel. plane Y coordinate";
    hCorr_Y_ref_4->SetTitle(title.c_str());
    _rootObjectMap.insert(make_pair(histoName, hCorr_Y_ref_4));
    //id5
    histoName = getHistoName("hCorrY_CMSPixRef_vs_",5);
    TH2D * hCorr_Y_ref_5 = new TH2D (histoName.c_str(),"", 80, 0, 79, 576, 0, 575);
    title = "Clusters Y correlation: tel. plane 5 vs. CMSPixRef;CMSPixRef Y coordinate; Tel. plane Y coordinate";
    hCorr_Y_ref_5->SetTitle(title.c_str());
    _rootObjectMap.insert(make_pair(histoName, hCorr_Y_ref_5));
    //DUT
    histoName = getHistoName("hCorrY_CMSPixRef_vs_",6);
    TH2D * hCorr_Y_ref_6 = new TH2D (histoName.c_str(),"", 80, 0, 79, 128, 0, 127);
    title = "Clusters Y correlation: DUT vs. CMSPixRef;CMSPixRef Y coordinate; Alibava channel number";
    hCorr_Y_ref_6->SetTitle(title.c_str());
    _rootObjectMap.insert(make_pair(histoName, hCorr_Y_ref_6));
}
//end of calculations, outout with the fractions of the clusters per event
void AlibavaClusterCollectionMerger::end ()
{
    streamlog_out ( MESSAGE5 )  << "AlibavaClusterCollectionMerger Successfully finished" << endl
                                << "Fraction tel+DUT+CMSPixRef = " << event_cluster_at_each_plane << endl
                                << "Fraction tel+CMSPixRef = " << event_cluster_at_tel_CMSRef << endl
                                << "Fraction tel+DUT = " << event_cluster_at_tel_DUT << endl
                                << "Fraction only tel = " << event_cluster_at_tel << endl;
}
//naming the histograms
string AlibavaClusterCollectionMerger::getHistoName(string histoName, int ID)
{
    stringstream s;
    if(ID < 6){s << histoName << "tel_"<< ID;}
    else if((ID == 6)&&((histoName == "hCorrY_ali_vs_") || (histoName == "hCorrX_ali_vs_"))){s << histoName << "CMSPixRef";}
    else if((ID == 6)&&((histoName == "hCorrX_CMSPixRef_vs_") || (histoName == "hCorrY_CMSPixRef_vs_"))){s << histoName << "ali";}
    return s.str();
}
