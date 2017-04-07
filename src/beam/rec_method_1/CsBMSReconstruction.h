#ifndef CsBMSReconstruction_h
#define CsBMSReconstruction_h

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <stdio.h>

#include "CsTypes.h"
#include "CsEvent.h"
#include "CsGeom.h"
#include "CsErrLog.h"
#include "CsOpt.h"
#include "DaqDataDecoding/DaqEvent.h"

#include "CsBTReconstruction.h"
#include "CsBeamReconstructionSupport.h"
#include "CsBeamMomCalculation.h"
#include "CsBeamBackPropagation.h"

class CsBeamMomCalculation;
class CsBeamBackPropagation;
class CsBTReconstruction;
class CsBMSslice;
class CsBMScluster;
class CsBeamCombination;
class CsBMStrack;


class CsBMSReconstruction{

   private:

   //pointer
   static CsBMSReconstruction* instance_;

   CsBeamMomCalculation *BeamMomCalculation;
   CsBeamBackPropagation *BeamBackPropagation;
   CsBTReconstruction *BT;

   //constants from file
   bool _useBM05;                               //use these planes
   bool _useBM06;

   std::string _planeName[6];                   //name of each plane

   double _clusterMaxDistMult;                  //max distance difference between two hits
   double _clusterMaxTimeMult;                  //max time difference between two hits (factor * time resolution)

   double _timeResolutionBMS[6];                //time resolution
   double _spaceResolutionBMS[6];               //space resolution

   double _hitTimeWindow;                       //max track time wrt BT track

   double _BTvsBMSTimeMult;
   double _BTvsBMSResolution[6];
   double _BTvsBMSCorrection[6];

   int _LHCalculationMode;
   double _reasonable_LH;

   double _slice_mom;                           //momentum range in which momntum is sliced 
   int _slice_n;                                //number of momentum slices 
   double _slice_overlap;                       //slice overlap

   double _s_LH_cut;
   double _t_LH_cut;

   bool _save_cs_clusters;                      //save CsCluster objects in case when BMS clusters are stored in mDST

   //constants set by functions
   int _planeSize[6];                           //size of each plane
   double* _planeDecodingMap[6];                //decoding map of each plane
   double _planeZcorrection[64][4];             //corrections to z position

   CsDetector* _planeReference[6];              //reference to each BMS plane

   //variables
   CsBMScluster* _clusterEmpty;

   std::vector<CsBMScluster*> _clustersBM01;
   std::vector<CsBMScluster*> _clustersBM02;
   std::vector<CsBMScluster*> _clustersBM03;
   std::vector<CsBMScluster*> _clustersBM04;
   std::vector<CsBMScluster*> _clustersBM05;
   std::vector<CsBMScluster*> _clustersBM06;

   std::vector<CsBeamCombination*> _BeamCombinations;  //BT <-> BMS combinations

   //set size of decetcors and decoding map
   void set_BMS_decoding_map();

   //check hit channel range
   bool check_BMS_channel_address(int ipl, int channel);

   //make clasterization of old planes (DAQ digits -> clusters)
   void make_BMS_old_clusterization();

   //get clusters of new planes
   void get_BMS_new_clusters();

   //histogram residuals
   void histogram_residuals(CsBMStrack* BMStrack);

   //search BMS tracks
   void search_BMS_tracks();

   //histograms
   int _his_mode;
   void set_histograms();

   CsHist2D* his_cluster_01[6];                 //channel vs. time
   CsHist2D* his_cluster_02[6];                 //hit position vs. time
   CsHist2D* his_cluster_03[6];                 //cluster position vs. mean time
   CsHist2D* his_cluster_04;                    //number of hits/event vs. plane
   CsHist2D* his_cluster_05;                    //number of clusters/event vs. plane
   CsHist2D* his_cluster_06;                    //number of hits/cluster vs. plane
   CsHist2D* his_hit_time_vs_trigger[6];        //hit time vs. trigger

   CsHist2D* his_BMS_vs_BT0;                     //BMS clusters vs. BT track clusters
   CsHist2D* his_BMS_vs_BT1;                     //BMS clusters vs. BT track clusters, tLH, sLH > 0.01

   CsHist1D* his_track_01[2];                   //number of fired planes/track
   CsHist2D* his_track_02[2];                   //fired planes vs. track class
   CsHist2D* his_track_03[2];                   //time vs. track class
   CsHist2D* his_track_04[2];                   //momentum vs. track class
   CsHist2D* his_track_05[2];                   //dt vs. track class
   CsHist2D* his_track_06[2];                   //tchi2 vs. track class
   CsHist2D* his_track_07[2];                   //schi2 vs. track class
   CsHist2D* his_track_08[2];                   //chi2 vs. track class

   CsHist1D* his_residual_01[6];                //residual
   CsHist2D* his_residual_02[6];                //y vs. residual
   CsHist2D* his_residual_03[6];                //dy vs. residual

   public:

   //constructor and desctructor 
   CsBMSReconstruction();
   ~CsBMSReconstruction();

   //clear class for next event
   void clear();

   //get planes references
   CsDetector** get_planeReference(){  return _planeReference; }

   //main reconstruction function
   void make_BMS_reconstruction(std::vector<const CsTrack*> BTtracks);

   //get BMS tracks
   const std::vector<CsBeamCombination*> get_beam_combinations(){ return _BeamCombinations; }

   //instance
   static CsBMSReconstruction* Instance();

};

#endif
