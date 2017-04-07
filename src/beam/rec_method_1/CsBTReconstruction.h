#ifndef CsBTReconstruction_h
#define CsBTReconstruction_h

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

#include "CsBeamReconstructionSupport.h"

class CsBTReconstruction{

   private:

   //pointer
   static CsBTReconstruction* instance_;

   //constants from file
   double _trackTimeWindow;

   //variables
   std::vector<const CsTrack*> _BTtracks;                  //reconstructed BT tracks 

   //search BT tracks
   void search_BT_tracks();

   //histograms
   int _his_mode;
   void set_histograms();

   //planes
   std::vector<const CsDetector*> _planeReference;

   CsHist2D* his_track_time;                                //track time 
   CsHist2D** his_track_diff;                               //track time - hit time for various planes
   CsHist2D** his_hit_time;                                 //hit time for various planes 

   public:

   //constructor and descructor 
   CsBTReconstruction();
   ~CsBTReconstruction();

   //clear class for next event
   void clear();

   //main reconstruction function
   void make_BT_reconstruction();

   //get BT tracks
   const std::vector<const CsTrack*> get_BT_tracks(){ return _BTtracks; }

   //instance
   static CsBTReconstruction* Instance();

};

#endif
