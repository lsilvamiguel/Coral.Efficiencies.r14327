#include "CsBTReconstruction.h"

using namespace std;

//=== CONSTRUCTOR ====================================================================
CsBTReconstruction :: CsBTReconstruction(){

   //read configuration file 
   read_configuration_key("beam", "histogram_level", _his_mode);
   read_configuration_key("beam_bt", "time_window", _trackTimeWindow);

   //set pointers to BT detectors
   list<CsDetector*>  det = CsGeom::Instance()->getDetectors();
   list<CsDetector*>::iterator idet;

   for(idet = det.begin(); idet != det.end(); idet++){

      //check position w.r.t. target
      if( (*idet)->getZcm() > CsGeom::Instance()->getTargetCenter() ) continue;

      //no BMS
      if( !( ((*idet)->GetTBName()).find("BM") == std::string::npos ) ) continue;

      //no veto
      if( ((*idet)->GetTBName()).find("V") == 0 ) continue;

      //add
      _planeReference.push_back(*idet);

   }

   //set histograms
   set_histograms();

   //clear variables 
   _BTtracks.clear();
}

//=== DESTRUCTOR =====================================================================
CsBTReconstruction :: ~CsBTReconstruction(){
   clear();
}

//=== CLAER VARIABLES ================================================================
void CsBTReconstruction :: clear(){

   _BTtracks.clear();

}

//=== MAKE RECONSTRUCTION ============================================================
void CsBTReconstruction :: make_BT_reconstruction(){

   //search BT tracks
   search_BT_tracks();
}

//=== GET BMS TRACKS =================================================================
void CsBTReconstruction :: search_BT_tracks(){

   //few variables
   bool isBT;

   //get list of all reconstructed by traffic tracks
   list<CsTrack*> tracks = CsEvent::Instance()->getTracks();
   list<CsTrack*>::iterator track;

   //loop over tracks
   for(track = tracks.begin(); track != tracks.end(); track++){

      //get list of helices
      vector<CsHelix> helices = (*track)->getHelices();
      vector<CsHelix>::iterator helix;

      //check if BT
      for(helix = helices.begin(), isBT = false; helix != helices.end(); helix++){
         if( helix->getZ() < CsGeom::Instance()->getTargetCenter() ) isBT = true;
      }

      //if not
      if( !isBT ) continue;

      //histogram
      for(int trig = 0; trig < 12; trig++){

         if( ((CsEvent::Instance()->getTriggerMask()) & (1<<trig)) == 0 ) continue;

         if( his_track_time != NULL) his_track_time->Fill((*track)->getMeanTime(), trig);

         //loop over clusters
         list<CsCluster*> clusters = (*track)->getClusters();
         list<CsCluster*>::iterator cluster;

         for(cluster = clusters.begin(); cluster != clusters.end(); cluster++){

            //loop over detectors
            list<CsDetector*> det = (*cluster)->getDetsList();
            list<CsDetector*>::iterator idet;

            for(idet = det.begin(); idet != det.end(); idet++) {

               //loop over BT planes
               for(unsigned int ipl = 0; ipl < _planeReference.size(); ipl++){

                  if( (*idet) != _planeReference[ipl]) continue;

                  double time;
                  if( !(*cluster)->getTime(time) ) continue;
                  if( his_track_diff[ipl] != NULL ) his_track_diff[ipl]->Fill((*track)->getMeanTime() - time, trig);
                  if( his_hit_time[ipl] != NULL ) his_hit_time[ipl]->Fill(time, trig);
               }
            }
         }
      }

      //check time
      if( (!(*track)->hasMeanTime()) || fabs((*track)->getMeanTime()) > _trackTimeWindow ) continue;

      //add
      _BTtracks.push_back(*track);
   }
}

//=== SET HISTOGRAMS ============================================================
void CsBTReconstruction :: set_histograms(){
   
   //"first" marker
   static bool first = true;
   char title[100], name[100];

   string histograms_path =  "/BeamReconstruction/BTReconstruction";
   CsHistograms::SetCurrentPath(histograms_path);

   if(first){

      his_track_time = NULL;
      
      his_track_diff = new CsHist2D*[_planeReference.size()];
      his_hit_time = new CsHist2D*[_planeReference.size()];

      for(unsigned int ipl = 0; ipl < _planeReference.size(); ipl++){
         his_track_diff[ipl] = NULL;
         his_hit_time[ipl] = NULL;
      }

      //track time vs trigger
      if( _his_mode > 0 ) his_track_time =  new CsHist2D("BTtrack_time",   "BT track time vs. trigger",  600, -30, 30, 12, 0, 12);

      for(unsigned int ipl = 0; ipl < _planeReference.size(); ipl++){

         //track time - cluster time for various planes vs. trigger
         strcpy(name, "BTtrack_diff_");
         strcat(name, (_planeReference[ipl]->GetTBName()).c_str());

         strcpy(title, "BT track time - cluster time vs. trigger for ");
         strcat(title, (_planeReference[ipl]->GetTBName()).c_str());
         if( _his_mode > 0 ) his_track_diff[ipl] = new  CsHist2D(name, title, 600, -30, 30, 12, 0, 12);

         //hit time vs. trigger 
         strcpy(name, "BThit_time_");
         strcat(name, (_planeReference[ipl]->GetTBName()).c_str());

         strcpy(title, "BT hit time vs. trigger for ");
         strcat(title, (_planeReference[ipl]->GetTBName()).c_str());
         if( _his_mode > 1 ) his_hit_time[ipl] = new  CsHist2D(name, title, 600, -30, 30, 12, 0, 12);
      }

      first = false;
   }
}

//=== INSTANCE ====================================================================
CsBTReconstruction* CsBTReconstruction :: instance_ = NULL;

CsBTReconstruction* CsBTReconstruction :: Instance() {
   if( instance_ == NULL ){
      instance_ = new CsBTReconstruction();
   }
       
   return( instance_ );
}

