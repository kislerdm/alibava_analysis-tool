/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *  email:eda.yildirim@cern.ch
 * Edited by Dmitry Kisler
 *  (2016 DESY)
 *  email: dkisler@cern.ch
*/

// alibava includes ".h"
#include "AlibavaSeedClustering.h"
#include "AlibavaRunHeaderImpl.h"
#include "AlibavaEventImpl.h"
#include "ALIBAVA.h"
#include "AlibavaPedNoiCalIOManager.h"
#include "AlibavaCluster.h"
// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"
#include "marlin/Global.h"
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
// aida includes <.h>
#include <marlin/AIDAProcessor.h>
#include <AIDA/ITree.h>
#endif
// lcio includes <.h>
#include <lcio.h>
#include <UTIL/LCTOOLS.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/LCTime.h>
#include <IMPL/LCEventImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/LCEventImpl.h>
// ROOT includes ".h"
#include "TH1D.h"
#include "TH2D.h"
// system includes <>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <memory>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace alibava;

AlibavaSeedClustering::AlibavaSeedClustering () :
    AlibavaBaseProcessor("AlibavaSeedClustering"),
    _seedCut(3),
    _neighCut(2),
    _sensitiveAxisX(1),
    _signalPolarity(-1),
    _etaHistoName("hEta"),
    _clusterSizeHistoName("hClusterSize"),
    _isSensitiveAxisX(true),
    _crosstalk_c1(0),
    _crosstalk_c2(0),
    _event_num(500000),
    _run_num("000000")
{

    // modify processor description
    _description =
            "AlibavaSeedClustering finds clusters using seed and neighbour cuts ";

    // first of register the input /output collection
    registerInputCollection (LCIO::TRACKERDATA, "InputCollectionName",
                             "Input collection name, it should be pedestal subtracted",
                             _inputCollectionName, string("recodata") );

    registerOutputCollection (LCIO::TRACKERDATA, "OutputCollectionName",
                              "Output data collection name",
                              _outputCollectionName, string("alibava_clusters") );

    // if needed one can change these to optional parameters

    registerProcessorParameter ("NoiseInputFile",
                                "The filename where the pedestal and noise values stored",
                                _pedestalFile , string("pedestal.slcio"));

    registerProcessorParameter ("NoiseCollectionName",
                                "Noise collection name, better not to change",
                                _noiseCollectionName, string ("noise"));

    registerProcessorParameter ("SeedSNRCut",
                                "The signal/noise ratio that channels have to pass to be considered as seed channel",
                                _seedCut, float (3));

    registerProcessorParameter ("NeighbourSNRCut",
                                "The signal/noise ratio that neigbour channels have to pass to be added to the cluster",
                                _neighCut, float (2));

    registerProcessorParameter ("IsSensitiveAxisX",
                                "The default sensitive axis of the strip sensor(s) according to telescope is X. If sensitive axis is Y then set this parameter to zero (0). Any other value will be disregarded and sensitive axis will assumed to be \"X\" ",
                                _sensitiveAxisX, int(1));

    registerProcessorParameter ("SignalPolarity",
                                "Polarity of the signal. Set this parameter to -1 for negative signals, any other value will be disregarded and the signal will be assumed to be positive ",
                                _signalPolarity, int (-1));
    //cross-talk noise correction parameters:
    registerProcessorParameter ("Crosstalk_c1",
                                "Cross-talk noise correction coefficient c1 (see details in the Thomas Eichhorn's PhD thesis, p. 5.3.5.2",
                                _crosstalk_c1, float (0));
    registerProcessorParameter ("Crosstalk_c2",
                                "Cross-talk noise correction coefficient c2 (see details in the Thomas Eichhorn's PhD thesis, p. 5.3.5.2",
                                _crosstalk_c2, float (0));
    registerProcessorParameter ("NumEventsHist",
                                "Number of events we want ot process. This number is being used for scaling the hists.",
                                _event_num, int(500000));
    registerProcessorParameter ("RunNum",
                                "The number of the run we analyze",
                                _run_num, string ("000000"));

}
void AlibavaSeedClustering::init ()
{
    streamlog_out ( MESSAGE4 ) << "Running init" << endl;

    // this method is called only once even when the rewind is active

    //     To set of channels to be used
    //     ex.The format should be like $ChipNumber:StartChannel-EndChannel$
    //     ex. $0:5-20$ $0:30-100$ $1:50-70$
    //     means from chip 0 channels between 5-20 and 30-100, from chip 1 channels between 50-70 will be used (all numbers included). the rest will be masked and not used
    //     Note that the numbers should be in ascending order and there should be no space between two $ character

    if (Global::parameters->isParameterSet(ALIBAVA::CHANNELSTOBEUSED))
        Global::parameters->getStringVals(ALIBAVA::CHANNELSTOBEUSED,_channelsToBeUsed);
    else {
        streamlog_out ( MESSAGE4 ) << "The Global Parameter "<< ALIBAVA::CHANNELSTOBEUSED <<" is not set!" << endl;
    }
    /* To choose if processor should skip masked events
     ex. Set the value to 0 for false, to 1 for true
     */
    if (Global::parameters->isParameterSet(ALIBAVA::SKIPMASKEDEVENTS))
        _skipMaskedEvents = bool ( Global::parameters->getIntVal(ALIBAVA::SKIPMASKEDEVENTS) );
    else {
        streamlog_out ( MESSAGE4 ) << "The Global Parameter "<< ALIBAVA::SKIPMASKEDEVENTS <<" is not set! Masked events will be used!" << endl;
    }

    // check signal Polarity (-1 for p-/Y- type sensors)
    if (_signalPolarity != -1)
        _signalPolarity = 1;

    // check sensitive axis (0: strips are oriented horizontally, 1: strips are oriented vertically)
    if (_sensitiveAxisX == 0)
        _isSensitiveAxisX = false;
    else
        _isSensitiveAxisX = true;

    // usually a good idea to
    printParameters ();

}

void AlibavaSeedClustering::processRunHeader (LCRunHeader * rdr)
{
    streamlog_out ( MESSAGE4 ) << "Running processRunHeader" << endl;

    // Add processor name to the runheader
    auto_ptr<AlibavaRunHeaderImpl> arunHeader ( new AlibavaRunHeaderImpl(rdr)) ;
    arunHeader->addProcessor(type());

    // get and set selected chips
    setChipSelection(arunHeader->getChipSelection());

    // set channels to be used (if it is defined)
    setChannelsToBeUsed();

    // set pedestal and noise values
    setPedestals();

    // if you want
    bookHistos();

    // set number of skipped events to zero (defined in AlibavaBaseProcessor)
    _numberOfSkippedEvents = 0;
}

void AlibavaSeedClustering::processEvent (LCEvent * anEvent)
{

    AlibavaEventImpl * alibavaEvent = static_cast<AlibavaEventImpl*> (anEvent);

    if (_skipMaskedEvents && (alibavaEvent->isEventMasked()) ) {
        _numberOfSkippedEvents++;
        return;
    }

    // As an input collection we get pedestal subtracted alibava data
    // This collection contains AlibavaEvents
    LCCollectionVec * inputColVec;

    // The Alibava Cluster collection
    LCCollectionVec * clusterColVec = new LCCollectionVec(LCIO::TRACKERDATA);
    // cell id encode for AlibavaCluster
    CellIDEncoder<TrackerDataImpl> clusterIDEncoder(ALIBAVA::ALIBAVACLUSTER_ENCODE,clusterColVec);

    unsigned int noOfChip;
    try
    {
        //how collection is being accessed here? what is getInputCollectionName()?
        inputColVec = dynamic_cast< LCCollectionVec * > ( alibavaEvent->getCollection( getInputCollectionName() ) ) ;
        noOfChip = inputColVec->getNumberOfElements();

        for ( size_t i = 0; i < noOfChip; ++i ){
            // get your data from the collection and do what ever you want

            TrackerDataImpl * trkdata = dynamic_cast< TrackerDataImpl * > ( inputColVec->getElementAt( i ) ) ;
            vector<AlibavaCluster> clusters = findClusters(trkdata, alibavaEvent);

            // loop over clusters
            for (unsigned int icluster=0; icluster<clusters.size(); icluster++) {
                AlibavaCluster acluster = clusters[icluster];
                // create a TrackerDataImpl for each cluster
                TrackerDataImpl * alibavaCluster = new TrackerDataImpl();
                acluster.createTrackerData(alibavaCluster);

                // now store chip number, seed channel, cluster ID, cluster size, sensitive axis and signal polarity in CellIDEncode
                // set cluster ID
                clusterIDEncoder[ALIBAVA::ALIBAVACLUSTER_ENCODE_CLUSTERID] = acluster.getClusterID();
                // set sensitive axis
                if (acluster.getIsSensitiveAxisX())
                    clusterIDEncoder[ALIBAVA::ALIBAVACLUSTER_ENCODE_ISSENSITIVEAXISX] = 1;
                else
                    clusterIDEncoder[ALIBAVA::ALIBAVACLUSTER_ENCODE_ISSENSITIVEAXISX] = 0;
                // set signal polarity
                if (acluster.getSignalPolarity() == -1)
                    clusterIDEncoder[ALIBAVA::ALIBAVACLUSTER_ENCODE_ISSIGNALNEGATIVE] = 1;
                else
                    clusterIDEncoder[ALIBAVA::ALIBAVACLUSTER_ENCODE_ISSIGNALNEGATIVE] = 0;
                // set chip number
                clusterIDEncoder[ALIBAVA::ALIBAVACLUSTER_ENCODE_CHIPNUM] = acluster.getChipNum();
                // set cluster size
                clusterIDEncoder[ALIBAVA::ALIBAVACLUSTER_ENCODE_CLUSTERSIZE] = acluster.getClusterSize();
                // set seed channel number
                clusterIDEncoder[ALIBAVA::ALIBAVACLUSTER_ENCODE_SEED] = acluster.getSeedChanNum();
                clusterIDEncoder.setCellID(alibavaCluster);

                clusterColVec->push_back(alibavaCluster);

            } // end of loop over clusters

        } // end of loop ever detectors

        alibavaEvent->addCollection(clusterColVec, getOutputCollectionName());

    } catch ( lcio::DataNotAvailableException ) {
        // do nothing again
        streamlog_out( ERROR5 ) << "Collection ("<<getInputCollectionName()<<") not found! " << endl;
    }
}
vector<AlibavaCluster> AlibavaSeedClustering::findClusters(TrackerDataImpl * trkdata, AlibavaEventImpl * alibavaEvent)
{
    // first get chip number
    int chipnum = getChipNum(trkdata);
    // then get the data vector before cross-talk noise correction
    FloatVec dataVec; //_noncorr;
    dataVec = trkdata->getChargeValues();
    //    dataVec_noncorr = trkdata->getChargeValues();
    //create a vector for storing the data after the correction of cross-talk noise
    //    FloatVec dataVec;
    //    dataVec.clear();
    // we will need noise vector too
    FloatVec noiseVec = getNoiseOfChip(chipnum);
    //data after the cross-talk correction
    //    FloatVec dataVec_corr;
    //    dataVec_corr.clear();
    
    // then check which channels we can add to a cluster
    // obviously not the ones masked
    vector<bool> channel_can_be_used;
    //how many strips do we take to account
    int useful_ch = 128;
    for (int ichan=0; ichan<int(dataVec.size()); ichan++)
    {
        channel_can_be_used.push_back(!isMasked(chipnum,ichan));
        if(isMasked(chipnum,ichan)) useful_ch--;
    }
    // Now mask channels that cannot pass NeighbourSNRCut
    // And add the channel numbers as seed candidate if it pass SeedSNRCut
    vector<int> seedCandidates;
    seedCandidates.clear();
    //how many signals do we omit
    int out_of_clustering_counter = 0;
    for (int ichan=0; ichan<int(channel_can_be_used.size()); ichan++)
    {
        // if this channel is already masked, jump to next channel ichan
        if (!channel_can_be_used[ichan]) continue;
        //cross-talk correction:
        //        dataVec_corr[ichan] = dataVec[ichan]*_crosstalk_c1-_crosstalk_c2);
        //        dataVec_corr[ichan + 1 ] =;

        // now calculate snr = signal/noise
        float snr = (_signalPolarity * dataVec[ichan])/noiseVec[ichan];
        //debugging
        //streamlog_out (MESSAGE4) << "signal = " << dataVec[ichan] <<" | noise = " << noiseVec[ichan] << " | snr =" << snr <<endl;
        //filling the SNR distr. hist
        TH2D * hSNRvsCh = dynamic_cast<TH2D*> (_rootObjectMap[getHistoNameForChip("hSNRvsCh", chipnum)]) ;
        hSNRvsCh->Fill(ichan, snr);
        //filling the signal dist hist
        TH1D * hSigvsCh = dynamic_cast<TH1D*> (_rootObjectMap[getHistoNameForChip("hSigvsCh", chipnum)]) ;
        hSigvsCh->Fill(ichan, _signalPolarity*dataVec[ichan]);
        //filling the 2d signal dist hist
        TH2D * hSigvsCh2d = dynamic_cast<TH2D*> (_rootObjectMap[getHistoNameForChip("hSigvsCh2d", chipnum)]) ;
        hSigvsCh2d->Fill(ichan, _signalPolarity*dataVec[ichan]);
        //filling the noise dist hist
        TH1D * hNoisevsCh = dynamic_cast<TH1D*> (_rootObjectMap[getHistoNameForChip("hNoisevsCh", chipnum)]) ;
        hNoisevsCh->Fill(ichan, noiseVec[ichan]);
        //filling the noise2d dist hist
        TH2D * hNoisevsCh2d = dynamic_cast<TH2D*> (_rootObjectMap[getHistoNameForChip("hNoisevsCh2d", chipnum)]) ;
        hNoisevsCh2d->Fill(ichan, noiseVec[ichan]);
        // mask channels that cannot pass NeighbourSNRCut
        if (snr < _neighCut)
        {
            channel_can_be_used[ichan] = false;
            out_of_clustering_counter++;
        }
        else if (snr > _seedCut) { seedCandidates.push_back(ichan); }
    }

    //filling the strip exclustion histogram
    //which shows the fraction of strips we excluded above from being taken to account for clustering

    if(useful_ch > 0)
    {
    TH2D * hStripExclusion = dynamic_cast<TH2D*> (_rootObjectMap[getHistoNameForChip("hStripExclusion", chipnum)]) ;
    hStripExclusion->Fill(alibavaEvent->getEventNumber(), ((float)out_of_clustering_counter/(float)useful_ch*100));
    }
    // sort seed channels according to their SNR, highest comes first!
    for (unsigned int i=0; i+1<seedCandidates.size(); i++)
    {
        int ichan = seedCandidates[i];
        float isnr = (_signalPolarity*dataVec[ichan])/noiseVec[ichan];
        int next_chan = seedCandidates[i+1];
        float next_snr = (_signalPolarity*dataVec[next_chan])/noiseVec[next_chan];
        if(isnr < next_snr){ // swap
            seedCandidates[i] = next_chan;
            seedCandidates[i+1] = ichan;
            //start from beginning
            i=0;
        }
    }
    // now form clusters starting from the seed channel that has highest SNR
    int clusterID = 0;
    //how many clusters did we find
    int number_of_clusters = 0;
    // form clusters and store them in a vector
    vector<AlibavaCluster> clusterVector;
    for (unsigned int iseed=0; iseed<seedCandidates.size(); iseed++)
    {
        // if this seed channel used in another cluster, skip it
        int seedChan = seedCandidates[iseed];
        if (!channel_can_be_used[seedChan]) continue;

        //define the cluster's borders
        int cluSizemin = 0;
        int cluSizemax = 0;

        AlibavaCluster acluster;
        acluster.setChipNum(chipnum);
        acluster.setSeedChanNum(seedChan);
        acluster.setEta( calculateEta(trkdata,seedChan) );
        acluster.setIsSensitiveAxisX(_isSensitiveAxisX);
        acluster.setSignalPolarity(_signalPolarity);
        acluster.setClusterID(clusterID);
        clusterID++;

        // add seed channel to the cluster
        acluster.add(seedChan, dataVec[seedChan]);
        //counting the numb of clusters per event
        number_of_clusters++;
        // now mask seed channel so no other cluster can use it!
        channel_can_be_used[seedChan]=false;
        //filling the Seed SNR distr hist
        float SeedSNR = _signalPolarity*dataVec[seedChan]/noiseVec[seedChan];
        TH1D * hSeedSNR = dynamic_cast<TH1D*> (_rootObjectMap[getHistoNameForChip("hSeedSNR", chipnum)]) ;
        hSeedSNR -> Fill(SeedSNR);
        //filling the Seed Hit Map
        TH1D * hSeedHitMap = dynamic_cast<TH1D*> (_rootObjectMap[getHistoNameForChip("hSeedHitMap", chipnum)]) ;
        hSeedHitMap -> Fill(seedChan);
        // We will check if one of the neighbours not bonded. If it is not this is not a good cluster we will not use this
        bool thereIsNonBondedChan = false;
        // add channels on the left
        int ichan = seedChan-1;
        //initally cluster consists of only the seed strip
        cluSizemin = seedChan;
        bool continueSearch = true;
        while (continueSearch)
        {
            // first if ichan < 0
            if (ichan < 0)
            {
                // stop search
                continueSearch = false;
                // this also means that left neighbour is not bonded
                thereIsNonBondedChan = true;
                // then exit while loop
                break;
            }
            // if chan masked
            if (isMasked(chipnum,ichan)==true)
            {
                // this means that it is not bonded.
                thereIsNonBondedChan = true;
                // then we are sure that there is no other channel on the left
                break;
                // We will still form the cluster but we will not use this cluster.
                // the members on the right that should be in this cluster, should not be used by some other cluster.
                // so to mark them not to be used we will continue search
            }
            // now check if this channel can be added to the cluster
            if ( channel_can_be_used[ichan] )
            {
                // if yes, add it
                acluster.add(ichan, dataVec[ichan]);
                //cluster "left" border is being shifter to the left by one strip
                cluSizemin--;
                // and mask it so that it will not be added to any other cluster
                channel_can_be_used[ichan]=false;
            }
            else{
                // if you cannot use this channel there is no other channel on the left
                break;
            }
            // go to the next channel on left
            ichan--;
        }

        continueSearch = true;
        ichan = seedChan+1;
        //initally cluster consists of only the seed strip
        cluSizemax = seedChan;
        // add channels on right
        while (continueSearch)
        {
            // first if ichan < 0
            if (ichan >= int(channel_can_be_used.size()) )
            {
                // stop search
                continueSearch = false;
                // this also means that right neighbour is not bonded
                thereIsNonBondedChan = true;
                // then exit while loop
                break;
            }
            // if chan masked
            if (isMasked(chipnum,ichan)==true)
            {
                // this means that it is not bonded.
                thereIsNonBondedChan = true;
                // then we are sure that there is no other channel on the left
                break;
            }
            // now check if this channel can be added to the cluster
            if ( channel_can_be_used[ichan] )
            {
                // if yes, add it
                acluster.add(ichan, dataVec[ichan]);
                //cluster "right" border is being shifter to the rigth by one strip
                cluSizemax++;
                // and mask it so that it will not be added to any other cluster
                channel_can_be_used[ichan]=false;
            }
            else
            {
                // if you cannot use this channel there is no other channel on the right
                break;
            }
            // go to the next channel on right
            ichan++;
        }
        // now if there is no neighbour not bonded
        if(thereIsNonBondedChan == false)
        {
            // fill the histograms and add them to the cluster vector
            fillHistos(acluster);
            clusterVector.push_back(acluster);
            //filling the hist hClusterSNR
            float clusterCharge = 0;
            float clusterNoise = 0;
            float clusterSNR = 0;
            for(int i = cluSizemin; i <=cluSizemax; i++)
            {
                clusterCharge += _signalPolarity*dataVec[i];
                clusterNoise += noiseVec[i];
            }
            clusterSNR = clusterCharge/clusterNoise;
            //filling the cluster SNR dist hist
            TH1D * hClusterSNR = dynamic_cast<TH1D*> (_rootObjectMap[getHistoNameForChip("hClusterSNR", chipnum)]) ;
            hClusterSNR -> Fill(clusterSNR);
            //filling the clusterCharge vs TDC distr hist
            float TDC = alibavaEvent->getEventTime();
            TH2D * hClusterChargevsTDC = dynamic_cast<TH2D*> (_rootObjectMap[getHistoNameForChip("hClusterChargevsTDC", chipnum)]) ;
            hClusterChargevsTDC -> Fill(TDC, clusterCharge);
        }
    }

    //filling the number of clusters per event vs. event hist
    TH2D * hClustPerEvvsEv = dynamic_cast<TH2D*> (_rootObjectMap[getHistoNameForChip("hClustPerEvvsEv", chipnum)]) ;
    hClustPerEvvsEv->Fill(alibavaEvent->getEventNumber(), number_of_clusters);

    return clusterVector;
}

float AlibavaSeedClustering::calculateEta(TrackerDataImpl *trkdata, int seedChan)
{
    // first get chip number
    int chipnum = getChipNum(trkdata);
    // then get the data vector
    FloatVec dataVec;
    dataVec = trkdata->getChargeValues();

    // we will multiply all signal values by _signalPolarity to work on positive signal always
    float seedSignal = _signalPolarity * dataVec.at(seedChan);

    // now we make an unrealistic signal
    float unrealisticSignal = -10000; // you cannot get this signal from alibava

    int leftChan = seedChan - 1;
    float leftSignal = unrealisticSignal;
    // check if the channel on the left is masked
    if ( leftChan >= 0 && isMasked(chipnum,leftChan)==false ) {
        leftSignal = _signalPolarity * dataVec.at(leftChan);
    }

    int rightChan = seedChan+1;
    float rightSignal = unrealisticSignal;
    // check if the channel on the right is masked
    if ( rightChan < int(dataVec.size()) && (isMasked(chipnum, rightChan) == false ))
    {
        rightSignal = _signalPolarity * dataVec.at(rightChan);
    }
    // now compare right and left channel to see which one is higher

    // if both right anf left channel is masked. Simply return -1
    // this case should not be saved by clustering algorithm anyways
    if (rightSignal == unrealisticSignal && leftSignal == unrealisticSignal )
    {
        streamlog_out (DEBUG1) << "Both neighbours are masked!"<<endl;
        return -1;
    }
    float eta = -1;
    // compare left and right signals
    // here both signal has to be positive (if not noise)
    // if one of the channel is masked, since it will be set to unrealisticSignal other signal will always be higher then the masked one.

    if ( leftSignal > rightSignal) {
        // then seed channel is on the right
        eta = leftSignal / ( leftSignal + seedSignal );
    }
    else {
        // seed channel is on the left
        eta = seedSignal / (seedSignal + rightSignal);
    }
    // Eta calculation: chargeOnLeftChannel / (chargeOnLeftChannel + chargeOnRightChannel)
    return eta;
}

void AlibavaSeedClustering::check (LCEvent * /* evt */ ) {
    // nothing to check here - could be used to fill check plots in reconstruction processor
}

void AlibavaSeedClustering::end() {

    if (_numberOfSkippedEvents > 0)
        streamlog_out ( MESSAGE5 ) << _numberOfSkippedEvents<<" events skipped since they are masked" << endl;
    streamlog_out ( MESSAGE4 ) << "Successfully finished" << endl;

}

void AlibavaSeedClustering::fillHistos(AlibavaCluster anAlibavaCluster){

    int ichip = anAlibavaCluster.getChipNum();

    float eta = anAlibavaCluster.getEta();
    TH1D * hEta = dynamic_cast<TH1D*> (_rootObjectMap[getHistoNameForChip(_etaHistoName,ichip)]);
    hEta->Fill(eta);

    int clusterSize = anAlibavaCluster.getClusterSize();
    TH1D * hClusterSize = dynamic_cast<TH1D*> (_rootObjectMap[getHistoNameForChip(_clusterSizeHistoName,ichip)]);
    hClusterSize->Fill( clusterSize );
}

void AlibavaSeedClustering::bookHistos(){

    AIDAProcessor::tree(this)->cd(this->name());

    EVENT::IntVec chipSelection = getChipSelection();
    string histoName;
    string title;
    for ( unsigned int i=0; i<chipSelection.size(); i++)
    {
        int ichip=chipSelection[i];

        // Clustersize histogram
        histoName = getHistoNameForChip(_clusterSizeHistoName,ichip);
        TH1D * hClusterSize = new TH1D (histoName.c_str(),"",10, 0, 10);
        title = string("Cluster Size (chip ") +to_string(ichip)+string(").Run:") + to_string(_run_num)
                +string(";Cluster Size;Number of Entries");
        hClusterSize->SetTitle(title.c_str());
        _rootObjectMap.insert(make_pair(histoName, hClusterSize));
        //Seed Hit Map
        histoName = getHistoNameForChip("hSeedHitMap",ichip);
        TH1D * hSeedHitMap = new TH1D (histoName.c_str(),"", 128, 0, 127);
        title = string("How many times cluster's seed was assignel to a certain ch. (chip ")
                +to_string(ichip)+string(").Run:")+ to_string(_run_num)+string(";Alibava channel number; Number of Entries");
        _rootObjectMap.insert(make_pair(histoName, hSeedHitMap));
        hSeedHitMap->SetTitle(title.c_str());
        // Eta histogram ClusterSize > 1
        histoName = getHistoNameForChip(_etaHistoName,ichip);
        TH1D * hEta = new TH1D (histoName.c_str(),"",100, -0.1, 1.1);
        title = string("Eta distribution ClusterSize > 1 (chip ")
                +to_string(ichip)+string(").Run:")+to_string(_run_num)+string(";Eta; Number of Entries");
        hEta->SetTitle(title.c_str());
        _rootObjectMap.insert(make_pair(histoName, hEta));
        //SNRvsChnum distribution histogram
        histoName = getHistoNameForChip("hSNRvsCh",ichip);
        int binnum = 100;
        TH2D * hSNRvsCh = new TH2D (histoName.c_str(),"", 128, 0, 127, 2*binnum, - binnum + 1, binnum);
        title = string("Signal-to-Noise ratio distribution among all events (chip ")
                +to_string(ichip)+string(").Run:")+to_string(_run_num)+string("; Alibava channel number; Signal-to-Noise ratio");
        _rootObjectMap.insert(make_pair(histoName, hSNRvsCh));
        hSNRvsCh->SetTitle(title.c_str());
        //SignalvsChnum distribution histogram
        histoName = getHistoNameForChip("hSigvsCh",ichip);
        TH1D * hSigvsCh = new TH1D (histoName.c_str(),"", 128, 0, 127);
        title = string("Signal distribution among all events (chip ")
                +to_string(ichip)+string(").Run:")+to_string(_run_num)+string(");Alibava channel number; Signal [ADC]");
        _rootObjectMap.insert(make_pair(histoName, hSigvsCh));
        hSigvsCh->SetTitle(title.c_str());
        //SignalvsChnum2d distribution histogram
        histoName = getHistoNameForChip("hSigvsCh2d",ichip);
        TH2D * hSigvsCh2d = new TH2D (histoName.c_str(),"", 128, 0, 127, 2000, -99, 100);
        title = string("Signal vs. Ch numb. distribution among all events (chip ")
                +to_string(ichip)+string(").Run:")+to_string(_run_num)+string(";Alibava channel number; Signal [ADC]");
        _rootObjectMap.insert(make_pair(histoName, hSigvsCh2d));
        hSigvsCh2d->SetTitle(title.c_str());
        //NoisevsChnum distribution histogram
        histoName = getHistoNameForChip("hNoisevsCh",ichip);
        TH1D * hNoisevsCh = new TH1D (histoName.c_str(),"", 128, 0, 127);
        title = string("Noise distribution among all events (chip ")
                +to_string(ichip)+string(").Run:")+to_string(_run_num)+string(";Alibava channel number; Noise [ADC]");
        _rootObjectMap.insert(make_pair(histoName, hNoisevsCh));
        hNoisevsCh->SetTitle(title.c_str());
        //NoisevsChnum2d distribution histogram
        histoName = getHistoNameForChip("hNoisevsCh2d",ichip);
        TH2D * hNoisevsCh2d = new TH2D (histoName.c_str(),"", 128, 0, 127, 2000, -99, 100);
        title = string("Noise vs. Ch numb. distribution among all events (chip ")
                +to_string(ichip)+string(").Run:")+to_string(_run_num)+string(";Alibava channel number; Noise [ADC]");
        _rootObjectMap.insert(make_pair(histoName, hNoisevsCh2d));
        hNoisevsCh2d->SetTitle(title.c_str());
        //ClusterSNR distribution histogram
        histoName = getHistoNameForChip("hClusterSNR",ichip);
        TH1D * hClusterSNR = new TH1D (histoName.c_str(),"", 500, 0, 49);
        title = string("Cluster SNR distribution (chip ")
                +to_string(ichip)+string(").Run:")+to_string(_run_num)+string(";Cluster SNR; Number of Entries");
        _rootObjectMap.insert(make_pair(histoName, hClusterSNR));
        hClusterSNR->SetTitle(title.c_str());
        //SeedSNR distribution histogram
        histoName = getHistoNameForChip("hSeedSNR",ichip);
        TH1D * hSeedSNR = new TH1D (histoName.c_str(),"", 500, 0, 49);
        title = string("Cluster seed SNR distribution (chip ")
                +to_string(ichip)+string(").Run:")+to_string(_run_num)+string(";Cluster seed SNR; Number of Entries");
        _rootObjectMap.insert(make_pair(histoName, hSeedSNR));
        hSeedSNR->SetTitle(title.c_str());
        //ClusterChargevsTDC distribution histogram
        histoName = getHistoNameForChip("hClusterChargevsTDC",ichip);
        TH2D * hClusterChargevsTDC = new TH2D (histoName.c_str(),"", 100, 0, 99, 100, 0, 99);
        title = string("Cluster charge vs. TDC distribution (chip ")
                +to_string(ichip)+string(").Run:")+to_string(_run_num)+string(";TDC [ns]; Cluster Charge [ADC]");
        _rootObjectMap.insert(make_pair(histoName, hClusterChargevsTDC));
        hClusterChargevsTDC->SetTitle(title.c_str());
        //how many clusters per event do we have
        histoName = getHistoNameForChip("hClustPerEvvsEv",ichip);
        TH2D * hClustPerEvvsEv = new TH2D (histoName.c_str(),"", _event_num, 0, (_event_num - 1), 10, 0, 9);
        title = string("Number of clusters per event (chip ")
                +to_string(ichip)+string(").Run:")+to_string(_run_num)+string("; Event Number; Number of Clusters per Event");
        _rootObjectMap.insert(make_pair(histoName, hClustPerEvvsEv));
        hClustPerEvvsEv->SetTitle(title.c_str());
        //how many strips (fraction) we exclude per event from clustering because of clustering condition
        histoName = getHistoNameForChip("hStripExclusion",ichip);
        TH2D * hStripExclusion = new TH2D (histoName.c_str(),"", _event_num, 0, (_event_num - 1), 1010, 0, 100);
        title = string("Fraction of strips per event we excluded by cuts (chip ")
                +to_string(ichip)+string(").Run:")+to_string(_run_num)+string("; Event Number; Fraction of excluded strips [%]");
        _rootObjectMap.insert(make_pair(histoName, hStripExclusion));
        hStripExclusion->SetTitle(title.c_str());

    } // end of loop over selected chips
    streamlog_out ( MESSAGE1 )  << "End of Booking histograms. " << endl;
}

string AlibavaSeedClustering::getHistoNameForChip(string histoName, int ichip)
{
    stringstream s;
    s<< histoName <<"_chip_" << ichip;
    return s.str();
}
