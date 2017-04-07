#ifndef CsBeamReconstruction_h
#define CsBeamReconstruction_h

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

#include "TTree.h"

#include "CsBMSReconstruction.h" 
#include "CsBTReconstruction.h"
#include "CsBeamBackPropagation.h"
#include "CsBeam.h"

class CsBMSReconstruction;
class CsBTReconstruction;
class CsBeamBackPropagation;
class CsBMStrack;

class CsBeamReconstruction{

   private:

   //pointer
   static CsBeamReconstruction* instance_;

   CsBMSReconstruction *BMS;                   //pointers to BMS and BT classes
   CsBTReconstruction *BT;

   CsBeamBackPropagation *BeamBackPropagation;

   //constants from file
   double _nominalMomentum;            //nominal momentum
   int _beamCharge;                    //beam charge
   bool _do_alignment;                 //prepare alignment tree switch

   double _momResolution;              //momentum resolution

   //variables
   std::vector<CsBeam*> _beamtracks;        //reconstructed beam tracks 
   CsZone *BMS_zone;                   //BMS zone

   //set hit pattern
   bool set_hit_pattern(CsBeam* beam, CsBMStrack* bms_track);

   //tree for alignment procedure 
   TTree* alignment_tree;

   Float_t alignment_y[6];
   Float_t alignment_z[6];
   Float_t alignment_LH;
   Int_t alignment_event[3];

   //histograms
   int _his_mode;
   void set_histograms();

   public:

   //constructor and destructor
   CsBeamReconstruction();
   ~CsBeamReconstruction();

   //clear class for next event
   void clear();

   //main reconstruction function
   void make_beam_reconstruction();

   //get beam tracks (function is compatibile with rest of the coral code)
   std::list<CsBeam*> getBeam();

   //instance
   static CsBeamReconstruction* Instance();

};

#endif
