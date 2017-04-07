#include "CsBeamReconstruction.h"

using namespace std;

//=== CONSTRUCTOR ====================================================================
CsBeamReconstruction :: CsBeamReconstruction(){

   //veriables
   int i;
   
   //read configuration file 
   read_configuration_key("beam", "histogram_level", _his_mode);

   read_configuration_key("beam", "nominal_mom", _nominalMomentum);
   read_configuration_key("beam", "nominal_res", _momResolution);
   read_configuration_key("beam", "beam_charge", _beamCharge);

   read_configuration_key("beam", "prepare_alignment", _do_alignment);

   if(_nominalMomentum != 160.){
      cout << "CsBeamReconstruction.cc: Warning: Beam reconstruction not optimized for given nominal momentum value" << endl;
      cout << "CsBeamReconstruction.cc: Reconstructed momentum will be rescaled by factor " << _nominalMomentum << "/160 (it should be ok)" << endl;
   }

   if(_beamCharge != 1){
      cout << "CsBeamReconstruction.cc: Warning: Beam reconstruction not optimized for given sign" << endl;
      cout << "CsBeamReconstruction.cc: Coefficients for positive sign are taken (it should be ok)" << endl;
   }

   //set pointers
   BMS = CsBMSReconstruction :: Instance();
   BT = CsBTReconstruction :: Instance();
   BeamBackPropagation = CsBeamBackPropagation :: Instance();

   //create BMS zone
   list<CsDetector*> BMS_zone_detectors;

   for( i = 0; i < 6; i++ ){
      if( BMS->get_planeReference()[i] != NULL ) BMS_zone_detectors.push_back( BMS->get_planeReference()[i] );
   }

   BMS_zone = new CsZone(-13719, -6128, "BMSzone", BMS_zone_detectors);

   //set histograms
   set_histograms();

   //clear variables
   _beamtracks.clear();
}

//=== DESTRUCTOR=======================================================================
CsBeamReconstruction :: ~CsBeamReconstruction(){
   
   delete BMS_zone;
   clear();
}

//=== CLEAR VARIABLES ================================================================
void CsBeamReconstruction :: clear(){

   vector<CsBeam*>::iterator beam;

   for(beam = _beamtracks.begin(); beam != _beamtracks.end(); beam++){ delete *beam;  *beam = NULL;  }  
   _beamtracks.clear();

   BMS->clear();
   BT->clear();

}

//=== MAKE BEAM RECONSTRUCTION =======================================================
void CsBeamReconstruction :: make_beam_reconstruction(){

   int ipl;

   const CsTrack* BT_track;
   CsBMStrack* BMS_track; 
   double LH;

   CsHelix* tmp_helix;
   vector<CsHelix> tmp_helices;   

   clear();
   BT->make_BT_reconstruction();
   BMS->make_BMS_reconstruction(BT->get_BT_tracks());

   //save tracks
   vector<CsBeamCombination*>::iterator BeamCombination;
   vector<CsBeamCombination*> BeamCombinations = BMS->get_beam_combinations();

   for(BeamCombination = BeamCombinations.begin(); BeamCombination != BeamCombinations.end(); BeamCombination++){
      
      //set null
      BT_track = NULL;
      BMS_track = NULL;

      //check if BMS and BT tracks corelate and get the best combination
      BMS_track = (*BeamCombination)->get_best_BMStrack(&LH);

      if(BMS_track == NULL) continue;
      BT_track = (*BeamCombination)->get_BTtrack(); 

      //create object
      CsBeam* new_beam = new CsBeam(*BT_track);

      //additional informations
      new_beam->setBmdata(BMS_track->get_t(), BMS_track->get_nPlanes(), BMS_track->get_nPlanes(), BMS_track->get_t(), BMS_track->get_tchi2(), BMS_track->get_schi2(), 99, 9999, BT_track->getMeanTime(), 99);

      //add hit pattern
      set_hit_pattern(new_beam, BMS_track);

      //add zone
      new_beam->addZone(*BMS_zone);

      //set momentum modify helices
      tmp_helices.clear();

      vector<CsHelix>helices = BT_track->getHelices();
      vector<CsHelix>::iterator helix;
      
      for(helix = helices.begin(); helix != helices.end(); helix++){

         tmp_helix = new CsHelix((*helix).getX(), (*helix).getY(), (*helix).getZ(), 
                     (*helix).getDXDZ(), (*helix).getDYDZ(), 
                     (_beamCharge*160./_nominalMomentum)/BMS_track->get_mom(),   //include sign and rescale momentum if nominal momentum diffrent then 160GeV
                     (*helix).getCov());

         tmp_helix->getCov()[14] = pow(_momResolution*_beamCharge/_nominalMomentum, 2);

         tmp_helices.push_back(*tmp_helix);
         delete tmp_helix;
      }

      new_beam->setHelices(tmp_helices);

      //add LH
      if(LH < 0.0001) LH = 0.;
      new_beam->setBackPropLH(LH);

      //residuals
      BeamBackPropagation->histogram_residuals((*BeamCombination)->get_lastBThelix(), BMS_track->get_mom(), BMS_track);

      //prepare alignment data
      if(_do_alignment){
         for(ipl = 0; ipl < 6; ipl++){
            alignment_y[ipl] = (((BMS_track->get_clusters())[ipl])->isEmpty())?(1.E6):(((BMS_track->get_clusters())[ipl])->get_y());
            alignment_z[ipl] = (((BMS_track->get_clusters())[ipl])->isEmpty())?(1.E6):(((BMS_track->get_clusters())[ipl])->get_z());
         }
         alignment_LH = LH;
         alignment_event[0] = CsEvent::Instance()->getRunNumber();
         alignment_event[1] = CsEvent::Instance()->getBurstNumber();
         alignment_event[2] = CsEvent::Instance()->getEventNumberInBurst();
      }

      //add CsClusters
      for(ipl = 0; ipl < 6; ipl++){
         if( ((BMS_track->get_clusters())[ipl])->get_cs_cluster() != NULL ) new_beam->addCluster(*(((BMS_track->get_clusters())[ipl])->get_cs_cluster()));
      }

      //push back
      _beamtracks.push_back(new_beam);

      new_beam = NULL;
   }

   //save alignment data if only one reconstructed beam track
      if(_do_alignment && _beamtracks.size() == 1) alignment_tree->Fill();
}

//=== SET HIT PATTERN =========================================================
bool CsBeamReconstruction :: set_hit_pattern(CsBeam* beam, CsBMStrack* bms_track){

   //few variables
   int ipl;
   int word, bit;

   unsigned int det_position;
   unsigned int* expected;
   unsigned int* hited;

   //get data
   map<CsDetector*,int> det2bit = CsEvent::Instance()->getDetectorMap();

   expected = const_cast<unsigned int*>(beam->getExpectedDetsBitmap());
   hited = const_cast<unsigned int*>(beam->getFiredDetsBitmap());

   //loop over planes
   for(ipl = 0; ipl < 6; ipl++){
      
      //check if plane is used
      if(BMS->get_planeReference()[ipl] == NULL){

         if( ipl < 4 ){
            cout << "CsBeamReconstruction::setHitPattern: Error: Undefined reference to BM0" << ipl << endl;
            return false;
         }

         continue;
      }
      
      //position in the detector map
      det_position = (*(det2bit.find( BMS->get_planeReference()[ipl] ))).second;

      word = int(det_position/32);
      bit =      det_position%32;

      //check range
      if( word >= CSTRACK_MAPSIZE ){
         cout << "CsBeamReconstruction::setHitPattern: Error: array _expectedDets/_firedDets is to small to store hit pattern" << endl;
         return false;
      }

      //set 0
      expected[word] = (expected[word] | 1<<bit) ^ 1<<bit;
      hited[word] = (hited[word] | 1<<bit) ^ 1<<bit;

      //set expected
      expected[word] = expected[word] | 1<<bit;

      //check if plane fired
      if( (bms_track->get_clusters()[ipl])->isEmpty() ) continue;

      //set hited
      hited[word] = hited[word] | 1<<bit;
   }

   return true;
}

//=== SET HISTOGRAMS ============================================================
void CsBeamReconstruction :: set_histograms(){

   //"first" marker
   static bool first = true;

   if(first){

      alignment_tree = NULL;

      string histograms_path =  "/BeamReconstruction/BeamReconstruction";
      CsHistograms::SetCurrentPath(histograms_path);

      if(_do_alignment){

         alignment_tree = new TTree("beam_alignment", "Data for beam alignment procedure");

         alignment_tree->Branch("run",       &alignment_event[0], "run/I");
         alignment_tree->Branch("spill",     &alignment_event[1], "spill/I");
         alignment_tree->Branch("inspill",   &alignment_event[2], "inspill/I");

         alignment_tree->Branch("LH", &alignment_LH, "LH/F");

         alignment_tree->Branch("y", alignment_y, "y[6]/F");
         alignment_tree->Branch("z", alignment_z, "z[6]/F");
      }

      first = false;
   }

}



//=== GET BEAM TRACKS ================================================================
list<CsBeam*> CsBeamReconstruction :: getBeam(){

   vector<CsBeam*>::iterator beamtrack;
   list<CsBeam*> beamtracks;

   for(beamtrack = _beamtracks.begin(), beamtracks.clear(); beamtrack != _beamtracks.end(); beamtrack++){
      beamtracks.push_back(*beamtrack);
   }

   return beamtracks;
}

//=== INSTANCE =======================================================================
CsBeamReconstruction* CsBeamReconstruction :: instance_ = NULL;

CsBeamReconstruction* CsBeamReconstruction :: Instance() {
   if( instance_ == NULL ){
      instance_ = new CsBeamReconstruction();
   }
       
   return( instance_ );
}

