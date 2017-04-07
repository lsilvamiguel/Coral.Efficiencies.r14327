#include "CsBMSReconstruction.h"

using namespace std;

//=== CONSTRUCTOR ====================================================================
CsBMSReconstruction :: CsBMSReconstruction(){

   //read configuration file 
   read_configuration_key("beam", "histogram_level", _his_mode);

   read_configuration_key("beam_bms", "use_BM05", _useBM05);
   read_configuration_key("beam_bms", "use_BM06", _useBM06);

   read_configuration_key("beam_bms", "BM01_name", _planeName[0]);
   read_configuration_key("beam_bms", "BM02_name", _planeName[1]);
   read_configuration_key("beam_bms", "BM03_name", _planeName[2]);
   read_configuration_key("beam_bms", "BM04_name", _planeName[3]);
   read_configuration_key("beam_bms", "BM05_name", _planeName[4]);
   read_configuration_key("beam_bms", "BM06_name", _planeName[5]);
   
   read_configuration_key("beam_bms", "time_resolution_BM01", _timeResolutionBMS[0]);
   read_configuration_key("beam_bms", "time_resolution_BM02", _timeResolutionBMS[1]);
   read_configuration_key("beam_bms", "time_resolution_BM03", _timeResolutionBMS[2]);
   read_configuration_key("beam_bms", "time_resolution_BM04", _timeResolutionBMS[3]);
   read_configuration_key("beam_bms", "time_resolution_BM05", _timeResolutionBMS[4]);
   read_configuration_key("beam_bms", "time_resolution_BM06", _timeResolutionBMS[5]);

   read_configuration_key("beam_bms", "space_resolution_BM01", _spaceResolutionBMS[0]);
   read_configuration_key("beam_bms", "space_resolution_BM02", _spaceResolutionBMS[1]);
   read_configuration_key("beam_bms", "space_resolution_BM03", _spaceResolutionBMS[2]);
   read_configuration_key("beam_bms", "space_resolution_BM04", _spaceResolutionBMS[3]);
   read_configuration_key("beam_bms", "space_resolution_BM05", _spaceResolutionBMS[4]);
   read_configuration_key("beam_bms", "space_resolution_BM06", _spaceResolutionBMS[5]);

   read_configuration_key("beam_bms", "time_window", _hitTimeWindow);
   read_configuration_key("beam_bms", "cluster_time_mult", _clusterMaxTimeMult);
   read_configuration_key("beam_bms", "cluster_space_mult", _clusterMaxDistMult);

   read_configuration_key("beam_bms", "btvsbms_time_resolution_BM01", _BTvsBMSResolution[0]);
   read_configuration_key("beam_bms", "btvsbms_time_resolution_BM02", _BTvsBMSResolution[1]);
   read_configuration_key("beam_bms", "btvsbms_time_resolution_BM03", _BTvsBMSResolution[2]);
   read_configuration_key("beam_bms", "btvsbms_time_resolution_BM04", _BTvsBMSResolution[3]);
   read_configuration_key("beam_bms", "btvsbms_time_resolution_BM05", _BTvsBMSResolution[4]);
   read_configuration_key("beam_bms", "btvsbms_time_resolution_BM06", _BTvsBMSResolution[5]);

   read_configuration_key("beam_bms", "btvsbms_time_correction_BM01", _BTvsBMSCorrection[0]);
   read_configuration_key("beam_bms", "btvsbms_time_correction_BM02", _BTvsBMSCorrection[1]);
   read_configuration_key("beam_bms", "btvsbms_time_correction_BM03", _BTvsBMSCorrection[2]);
   read_configuration_key("beam_bms", "btvsbms_time_correction_BM04", _BTvsBMSCorrection[3]);
   read_configuration_key("beam_bms", "btvsbms_time_correction_BM05", _BTvsBMSCorrection[4]);
   read_configuration_key("beam_bms", "btvsbms_time_correction_BM06", _BTvsBMSCorrection[5]);

   read_configuration_key("beam_bms", "btvsbms_time_mult", _BTvsBMSTimeMult);

   read_configuration_key("beam_bms", "slice_mom", _slice_mom);
   read_configuration_key("beam_bms", "slice_n", _slice_n);
   read_configuration_key("beam_bms", "slice_overlap", _slice_overlap);

   read_configuration_key("beam_bms", "lh_calculation_mode", _LHCalculationMode);
   read_configuration_key("beam_bms", "reasonable_lh", _reasonable_LH);

   read_configuration_key("beam_bms", "s_lh_cut", _s_LH_cut);
   read_configuration_key("beam_bms", "t_lh_cut", _t_LH_cut);
   
   //save CsClusters? 
   _save_cs_clusters = false;

   list<string> det_in_mDST;
   list<string>::iterator idet_in_mDST;

   CsOpt* opt = CsOpt::Instance();

   opt->getOpt("mDST", "hits", det_in_mDST);
   for(idet_in_mDST = det_in_mDST.begin(); idet_in_mDST != det_in_mDST.end(); idet_in_mDST++){
      if( (*idet_in_mDST) == "BM" ) _save_cs_clusters = true;
   }

   //few variables
   int i;

   //Set pointers to detectors
   list<CsDetector*>  det = CsGeom::Instance()->getDetectors();
   list<CsDetector*>::iterator idet;

   for(i = 0; i < 6; i++) {

      _planeReference[i] = NULL;
      
      for(idet = det.begin(); idet != det.end(); idet++) {
         if(*idet == NULL) continue;
         if((*idet)->GetTBName() == _planeName[i]){
            _planeReference[i] = *idet;
            break;
         }
      }

      if(_planeReference[i] == NULL){
         cout << "CsBMSReconstruction.cc: Error: Reference to BMS plane is not existed" << endl;
         cout << "CsBMSReconstruction.cc: Plane name = " << _planeName[i] << endl;
         CsErrLog::Instance()->mes(elFatal, "Exiting ...");
      }
   }

   //set size of decetcors and decoding table 
   set_BMS_decoding_map();

   //set beam calculation
   BeamMomCalculation = CsBeamMomCalculation :: Instance();

   //set back propagation
   BeamBackPropagation = CsBeamBackPropagation :: Instance();

   //set BT reconstruction
   BT = CsBTReconstruction :: Instance(); 

   //set histograms
   set_histograms();

   //clar variables 
   _clusterEmpty = NULL;
   _clustersBM01.clear();
   _clustersBM02.clear();
   _clustersBM03.clear();
   _clustersBM04.clear();
   _clustersBM05.clear();
   _clustersBM06.clear();
   
}

//=== DESTRUCTOR =====================================================================
CsBMSReconstruction :: ~CsBMSReconstruction(){
   clear();
}

//=== CLEAR VARIABLES ================================================================
void CsBMSReconstruction :: clear(){

   vector<CsBMScluster*>::iterator cluster;
   vector<CsBeamCombination*>::iterator combination;

   delete _clusterEmpty;
   _clusterEmpty = NULL;
   
   for(cluster = _clustersBM01.begin(); cluster != _clustersBM01.end(); cluster++){ delete *cluster;  *cluster = NULL;  }  _clustersBM01.clear();
   for(cluster = _clustersBM02.begin(); cluster != _clustersBM02.end(); cluster++){ delete *cluster;  *cluster = NULL;  }  _clustersBM02.clear();
   for(cluster = _clustersBM03.begin(); cluster != _clustersBM03.end(); cluster++){ delete *cluster;  *cluster = NULL;  }  _clustersBM03.clear();
   for(cluster = _clustersBM04.begin(); cluster != _clustersBM04.end(); cluster++){ delete *cluster;  *cluster = NULL;  }  _clustersBM04.clear();
   for(cluster = _clustersBM05.begin(); cluster != _clustersBM05.end(); cluster++){ delete *cluster;  *cluster = NULL;  }  _clustersBM05.clear();
   for(cluster = _clustersBM06.begin(); cluster != _clustersBM06.end(); cluster++){ delete *cluster;  *cluster = NULL;  }  _clustersBM06.clear();

   for(combination = _BeamCombinations.begin(); combination != _BeamCombinations.end(); combination++){ delete *combination;  *combination = NULL;  }  
   _BeamCombinations.clear();
}

//=== MAKE RECONSTRUCTION ============================================================
void CsBMSReconstruction :: make_BMS_reconstruction(vector<const CsTrack*> BTtracks){

   //create empty cluster
   CsBMScluster *_clusterEmpty = new CsBMScluster();

   //set combinations
   vector<const CsTrack*>::iterator BTtrack; 
   CsBeamCombination *new_beam_combination;

   for(BTtrack = BTtracks.begin(); BTtrack != BTtracks.end(); BTtrack++){

      new_beam_combination = new CsBeamCombination(*BTtrack, _hitTimeWindow, _slice_mom, _slice_n, _slice_overlap, _clusterEmpty);
      _BeamCombinations.push_back(new_beam_combination);

      new_beam_combination = NULL;
   }

   //make clusterization of old planes
   make_BMS_old_clusterization();

   //get clusters of new planes
   get_BMS_new_clusters();

   //search BMS tracks
   search_BMS_tracks();
}


//=== CHECK HIT RANGE ================================================================
bool CsBMSReconstruction :: check_BMS_channel_address(int ipl, int channel){

   if( channel < 0 || channel > _planeSize[ipl] ){
      cout << "CsBMSReconstruction.cc: Warning: Hit channel is out of range" << endl;
      cout << "CsBMSReconstruction.cc: Plane name = " << _planeName[ipl] << " Channel:Limit = " << channel << ":" << _planeSize[ipl] << endl;

      return false;
   } 

   return true;
}


//=== MAKE CLUSTERIZATION OF OLD PLANES ==============================================
void CsBMSReconstruction :: make_BMS_old_clusterization(){

   //few variables
   int ipl;
   unsigned int i;
   bool* used_in_cluster;
   bool is_correlated;

   vector<CsBeamCombination*>::iterator BeamCombination;

   vector<double> time_in_cluster;
   vector<double> position_y_in_cluster;
   vector<double> position_z_in_cluster;
   vector<CsDigit*> cs_hits_in_cluster;

   double time, position_y, position_z, dtime, dposition_y;

   CsBMScluster* new_cluster;
   
   for(ipl = 0; ipl < 4; ipl++) { 

      //get digits
      list<CsDigit*>digits = _planeReference[ipl]->getMyDigits(); 
      list<CsDigit*>::iterator idig1;  unsigned int int_idig1;
      list<CsDigit*>::iterator idig2;  unsigned int int_idig2;

      //histogram
      if( his_cluster_04 != NULL ) his_cluster_04->Fill(digits.size(), ipl);

      //if 0
      if(digits.size() == 0) continue;

      //create table which contains flags for used in cluster hits
      used_in_cluster = new bool[digits.size()];

      for(idig1 = digits.begin(), int_idig1 = 0; idig1 != digits.end(); idig1++, int_idig1++) {

         //histogram
         if( his_cluster_01[ipl] != NULL ) his_cluster_01[ipl]->Fill((*idig1)->getAddress(), (*idig1)->getDatum());
         if( his_cluster_02[ipl] != NULL ) his_cluster_02[ipl]->Fill(_planeDecodingMap[ipl][(*idig1)->getAddress()], (*idig1)->getDatum());
         if( his_hit_time_vs_trigger[ipl] != NULL ){
            for(int trig = 0; trig < 12; trig++){

                if( ((CsEvent::Instance()->getTriggerMask()) & (1<<trig)) == 0 ) continue;
                his_hit_time_vs_trigger[ipl]->Fill((*idig1)->getDatum(), trig);
            }
         }
            
         //check if can be correlated with any BT track
         for(BeamCombination = _BeamCombinations.begin(), is_correlated = false; BeamCombination != _BeamCombinations.end(); BeamCombination++){

            if( (*BeamCombination)->check_time((*idig1)->getDatum()) ){
               is_correlated = true;
               break;
            }
         }

         //check hit range
         if( !check_BMS_channel_address(ipl, (*idig1)->getAddress()) ){
            is_correlated = false;
         } 

         //clear flag if everything is ok, if not mark as already used
         if(is_correlated){
            used_in_cluster[int_idig1] = false;
         }else used_in_cluster[int_idig1] = true;
      }

      //loops over hits
      for(idig1 = digits.begin(), int_idig1 = 0; idig1 != digits.end(); idig1++, int_idig1++) {

         //check if used
         if( used_in_cluster[int_idig1] ) continue;

         //clear vectors for new cluster
         time_in_cluster.clear();
         position_y_in_cluster.clear();
         position_z_in_cluster.clear();
         if( _save_cs_clusters ) cs_hits_in_cluster.clear();

         //add first hit before second loop
         time_in_cluster.push_back((*idig1)->getDatum());
         position_y_in_cluster.push_back(_planeDecodingMap[ipl][(*idig1)->getAddress()]);
         position_z_in_cluster.push_back(_planeReference[ipl]->getZcm() + _planeZcorrection[(*idig1)->getAddress()][ipl]); 
         if( _save_cs_clusters ) cs_hits_in_cluster.push_back(*idig1);

         //mark as used
         used_in_cluster[int_idig1] = true;

         //second loop over hits
         for(idig2 = digits.begin(), int_idig2 = 0; idig2 != digits.end(); idig2++, int_idig2++) {

            //check if used
            if( used_in_cluster[int_idig2] ) continue;

            //check conditions
            for(i=0, is_correlated = true; i < time_in_cluster.size(); i++){
               if( 
                  !( fabs( time_in_cluster[i] - (*idig2)->getDatum() ) < (_clusterMaxTimeMult * _timeResolutionBMS[ipl]) ) || 
                  !( fabs( position_y_in_cluster[i] - _planeDecodingMap[ipl][(*idig2)->getAddress()] ) < (_clusterMaxDistMult * _spaceResolutionBMS[ipl]) ) 
               ){
                  is_correlated = false;
                  break;
               }
            }

            //hit is used
            if( is_correlated ){

               //add hit
               time_in_cluster.push_back((*idig2)->getDatum());
               position_y_in_cluster.push_back(_planeDecodingMap[ipl][(*idig2)->getAddress()]);
               position_z_in_cluster.push_back(_planeReference[ipl]->getZcm() + _planeZcorrection[(*idig2)->getAddress()][ipl]);
               if( _save_cs_clusters ) cs_hits_in_cluster.push_back(*idig2);

               //mark as used
               used_in_cluster[int_idig2] = true;
            }
         }

         //calculate mean values
         for(i=0, time=0; i < time_in_cluster.size(); i++) time += time_in_cluster[i];    
            time /= (double)(time_in_cluster.size());
         for(i=0, position_y=0; i < position_y_in_cluster.size(); i++) position_y += position_y_in_cluster[i];   
            position_y /= (double)(position_y_in_cluster.size());
         for(i=0, position_z=0; i < position_z_in_cluster.size(); i++) position_z += position_z_in_cluster[i];   
            position_z /= (double)(position_z_in_cluster.size());

         dtime = _timeResolutionBMS[ipl];
         dposition_y = _spaceResolutionBMS[ipl];

         //histogram
         if( his_cluster_03[ipl] != NULL ) his_cluster_03[ipl]->Fill(position_y, time);
         if( his_cluster_06 != NULL ) his_cluster_06->Fill(time_in_cluster.size(), ipl);

         //save cluster in memory
         new_cluster = new CsBMScluster(time, dtime, position_y, dposition_y, position_z);

         if( _save_cs_clusters ) new_cluster->set_cs_cluster(_planeReference[ipl], cs_hits_in_cluster);

         if( ipl == 0 ){ _clustersBM01.push_back(new_cluster);}
         if( ipl == 1 ){ _clustersBM02.push_back(new_cluster);}
         if( ipl == 2 ){ _clustersBM03.push_back(new_cluster);}
         if( ipl == 3 ){ _clustersBM04.push_back(new_cluster);}

         //add cluster to tracking 
         for(BeamCombination = _BeamCombinations.begin(); BeamCombination != _BeamCombinations.end(); BeamCombination++){
            (*BeamCombination)->add_cluster(ipl, new_cluster);
         }

         new_cluster = NULL;
      }

      //remove table
      delete used_in_cluster;
   }

   //histogram
   if( his_cluster_05 != NULL ) his_cluster_05->Fill(_clustersBM01.size(), 0);
   if( his_cluster_05 != NULL ) his_cluster_05->Fill(_clustersBM02.size(), 1);
   if( his_cluster_05 != NULL ) his_cluster_05->Fill(_clustersBM03.size(), 2);
   if( his_cluster_05 != NULL ) his_cluster_05->Fill(_clustersBM04.size(), 3);
}


//=== GET CLUSTERS FROM NEW PLANES ===================================================
void CsBMSReconstruction :: get_BMS_new_clusters(){

   bool is_correlated;
   int ipl;

   double dt;

   vector<CsBeamCombination*>::iterator BeamCombination;

   CsBMScluster* new_cluster;

   //loop over planes
   for(ipl = 4; ipl < 6; ipl++){

      //checked if used    
      if(ipl == 4 && (!_useBM05 )) continue;
      if(ipl == 5 && (!_useBM06 )) continue;

      //get reference to digits 
      list<CsDigit*>digits = _planeReference[ipl]->getMyDigits();
      list<CsDigit*>::iterator idig;

      //loop over digits 
      for(idig = digits.begin(); idig != digits.end(); idig++){

         //histogram
         if( his_cluster_03[ipl] != NULL ) his_cluster_03[ipl]->Fill(_planeDecodingMap[ipl][(*idig)->getAddress()], (*idig)->getDatum());

         if( his_hit_time_vs_trigger[ipl] != NULL ){
            for(int trig = 0; trig < 12; trig++){

               if( ((CsEvent::Instance()->getTriggerMask()) & (1<<trig)) == 0 ) continue;
               his_hit_time_vs_trigger[ipl]->Fill((*idig)->getDatum(), trig);
            }
         }

         //check range
         if( !check_BMS_channel_address(ipl, (*idig)->getAddress()) ) continue;
 
         //check if can be correlated with any BT track
         for(BeamCombination = _BeamCombinations.begin(), is_correlated = false; BeamCombination != _BeamCombinations.end(); BeamCombination++){
            if( (*BeamCombination)->check_time((*idig)->getDatum()) ){
               is_correlated = true;
               break;
            }
         }

         if(!is_correlated) continue;

         //save cluster in memory (errors acording to NIM paper)
         new_cluster = new CsBMScluster((*idig)->getDatum(), _timeResolutionBMS[ipl], _planeDecodingMap[ipl][(*idig)->getAddress()], _spaceResolutionBMS[ipl], _planeReference[ipl]->getZcm());

         if( _save_cs_clusters ){
            vector<CsDigit*> cs_hits_in_cluster;
            cs_hits_in_cluster.push_back(*idig);
            new_cluster->set_cs_cluster(_planeReference[ipl], cs_hits_in_cluster);
         }

         if( ipl == 4 ){ _clustersBM05.push_back(new_cluster); }
         if( ipl == 5 ){ _clustersBM06.push_back(new_cluster); }

         //add cluster to tracking 
         for(BeamCombination = _BeamCombinations.begin(); BeamCombination != _BeamCombinations.end(); BeamCombination++){
            (*BeamCombination)->add_cluster(ipl, new_cluster);
         }

         new_cluster = NULL;

      }

   }

   //histogram
   if( his_cluster_05 != NULL ) his_cluster_05->Fill(_clustersBM05.size(), 4);
   if( his_cluster_05 != NULL ) his_cluster_05->Fill(_clustersBM06.size(), 5);
}

//=== SEARCH BMS TRACK =========================================================
void CsBMSReconstruction :: search_BMS_tracks(){

   //itherators
   vector<CsBMScluster*>::iterator clus1, clus2, clus3, clus4, clus5, clus6;
   vector<CsBMScluster*>::iterator clus; 

   vector<CsBeamCombination*>::iterator BeamCombination;
   vector<CsBMSslice*>::iterator BMSslice;
   vector<CsBMStrack*>::iterator old_BMStrack;

   //few veriables
   int n_planes, i_plane;
   double dt;

   double mom, d_mom;

   double LH;

   CsBMStrack* new_BMStrack;

   //loop over beam combinations
   for(BeamCombination = _BeamCombinations.begin(); BeamCombination != _BeamCombinations.end(); BeamCombination++){

      //loop over slices
      vector<CsBMSslice*> BMSslices = (*BeamCombination)->get_BMSslices();
      for(BMSslice = BMSslices.begin(); BMSslice != BMSslices.end(); BMSslice++){

         vector<CsBMScluster*> clustersBM01 = (*BMSslice)->get_clusters(0);
         vector<CsBMScluster*> clustersBM02 = (*BMSslice)->get_clusters(1);
         vector<CsBMScluster*> clustersBM03 = (*BMSslice)->get_clusters(2);
         vector<CsBMScluster*> clustersBM04 = (*BMSslice)->get_clusters(3);
         vector<CsBMScluster*> clustersBM05 = (*BMSslice)->get_clusters(4);
         vector<CsBMScluster*> clustersBM06 = (*BMSslice)->get_clusters(5);
   
         //loops overs hits/clusters with checking of time correlation between hits
         for(clus1 = clustersBM01.begin(); clus1 != clustersBM01.end(); clus1++){
            if( (!(*clus1)->isEmpty()) && (fabs((*clus1)->get_t()-((*BeamCombination)->get_BTtrack())->getMeanTime()-_BTvsBMSCorrection[0]) >
                 _BTvsBMSTimeMult*_BTvsBMSResolution[0] ) ) continue;
         for(clus2 = clustersBM02.begin(); clus2 != clustersBM02.end(); clus2++){
            if( (!(*clus2)->isEmpty()) && (fabs((*clus2)->get_t()-((*BeamCombination)->get_BTtrack())->getMeanTime()-_BTvsBMSCorrection[1]) >
                 _BTvsBMSTimeMult*_BTvsBMSResolution[1] ) ) continue;
         for(clus3 = clustersBM03.begin(); clus3 != clustersBM03.end(); clus3++){
            if( (!(*clus3)->isEmpty()) && (fabs((*clus3)->get_t()-((*BeamCombination)->get_BTtrack())->getMeanTime()-_BTvsBMSCorrection[2]) >
                 _BTvsBMSTimeMult*_BTvsBMSResolution[2] ) ) continue;
         for(clus4 = clustersBM04.begin(); clus4 != clustersBM04.end(); clus4++){
            if( (!(*clus4)->isEmpty()) && (fabs((*clus4)->get_t()-((*BeamCombination)->get_BTtrack())->getMeanTime()-_BTvsBMSCorrection[3]) >
                 _BTvsBMSTimeMult*_BTvsBMSResolution[3] ) ) continue;
         for(clus5 = clustersBM05.begin(); clus5 != clustersBM05.end(); clus5++){
            if( (!(*clus5)->isEmpty()) && (fabs((*clus5)->get_t()-((*BeamCombination)->get_BTtrack())->getMeanTime()-_BTvsBMSCorrection[4]) >
                 _BTvsBMSTimeMult*_BTvsBMSResolution[4] ) ) continue;
         for(clus6 = clustersBM06.begin(); clus6 != clustersBM06.end(); clus6++){
            if( (!(*clus6)->isEmpty()) && (fabs((*clus6)->get_t()-((*BeamCombination)->get_BTtrack())->getMeanTime()-_BTvsBMSCorrection[5]) >
                 _BTvsBMSTimeMult*_BTvsBMSResolution[5] ) ) continue;
   
            //number of used planes
            n_planes = 0;
            if(!(*clus1)->isEmpty()) n_planes++;
            if(!(*clus2)->isEmpty()) n_planes++;
            if(!(*clus3)->isEmpty()) n_planes++;
            if(!(*clus4)->isEmpty()) n_planes++;
            if(!(*clus5)->isEmpty()) n_planes++;
            if(!(*clus6)->isEmpty()) n_planes++;
   
            //remove 0, 1 and 2 hit tracks
            if( n_planes < 3 ) continue;
   
            //check if track has hit in upstream part of the BMS
            if( (*clus1)->isEmpty() && (*clus2)->isEmpty() && (*clus5)->isEmpty() ) continue;
   
            //check if track has hit in downstream part of the BMS
            if( (*clus3)->isEmpty() && (*clus4)->isEmpty() && (*clus6)->isEmpty() ) continue;
   
            //create BMS track
            new_BMStrack = new CsBMStrack(n_planes, *clus1, *clus2, *clus3, *clus4, *clus5, *clus6);
   
            //calculate time
            new_BMStrack->calculate_time(((*BeamCombination)->get_BTtrack())->getMeanTime(), _BTvsBMSResolution, _BTvsBMSCorrection);
   
            //calculate momentum
            BeamMomCalculation->calculate_momentum(new_BMStrack, &mom, &d_mom);
            new_BMStrack->set_mom(mom);
            new_BMStrack->set_dmom(d_mom);
   
            //calculate schi2 (from back propagation)
            BeamBackPropagation->make_back_propagation((*BeamCombination)->get_lastBThelix(), new_BMStrack->get_mom());
            new_BMStrack->set_schi2(BeamBackPropagation->calculate_chi2(new_BMStrack));
   
            //cuts on prob of tchi2 and schi2
            if(_LHCalculationMode == 3 ){
               if( ( TMath::Prob(new_BMStrack->get_tchi2(), new_BMStrack->get_nPlanes()-1) < _t_LH_cut ) || ( TMath::Prob(new_BMStrack->get_schi2(), new_BMStrack->get_nPlanes()-1) < _s_LH_cut ) ){
   
                  new_BMStrack->set_is_wrong();
               }
            }
   
            //calculate total chi2
            if(_LHCalculationMode == 0){
               new_BMStrack->set_LH(-1.);
            }
            else if(_LHCalculationMode == 1){
               new_BMStrack->set_LH(TMath::Prob(new_BMStrack->get_tchi2(), new_BMStrack->get_nPlanes()-1));
            }
            else if(_LHCalculationMode == 2){
               new_BMStrack->set_LH(TMath::Prob(new_BMStrack->get_schi2(), new_BMStrack->get_nPlanes()-1));
            }
            else if(_LHCalculationMode == 3){
               new_BMStrack->set_LH(TMath::Prob(new_BMStrack->get_tchi2()+new_BMStrack->get_schi2(), 2*new_BMStrack->get_nPlanes()-2));
            }
   
            //histogram
            if( his_track_01[0] != NULL ) his_track_01[0]->Fill(new_BMStrack->get_nPlanes());
            if( his_track_02[0] != NULL ) his_track_02[0]->Fill(new_BMStrack->get_cluster_mask(), new_BMStrack->get_nPlanes());
            if( his_track_03[0] != NULL ) his_track_03[0]->Fill(new_BMStrack->get_t(), new_BMStrack->get_nPlanes());
            if( his_track_04[0] != NULL ) his_track_04[0]->Fill(new_BMStrack->get_mom(), new_BMStrack->get_nPlanes());
            if( his_track_05[0] != NULL ) his_track_05[0]->Fill(new_BMStrack->get_t()-((*BeamCombination)->get_BTtrack())->getMeanTime(), new_BMStrack->get_nPlanes());
            if( his_track_06[0] != NULL ) his_track_06[0]->Fill(TMath::Prob(new_BMStrack->get_tchi2(), new_BMStrack->get_nPlanes()-1), new_BMStrack->get_nPlanes());
            if( his_track_07[0] != NULL ) his_track_07[0]->Fill(TMath::Prob(new_BMStrack->get_schi2(), new_BMStrack->get_nPlanes()-1), new_BMStrack->get_nPlanes());
            if( his_track_08[0] != NULL ) his_track_08[0]->Fill(new_BMStrack->get_LH(), new_BMStrack->get_nPlanes());
   
            //best one
            vector<CsBMStrack*> old_BMStracks = (*BeamCombination)->get_BMStracks();
            for(old_BMStrack = old_BMStracks.begin(); old_BMStrack != old_BMStracks.end(); old_BMStrack++){
   
               //check if already marked as wrong
               if( (*old_BMStrack)->is_wrong() ) continue;
   
               //check chi2
               if( new_BMStrack->get_nPlanes() > (*old_BMStrack)->get_nPlanes() ){
                  if( ( new_BMStrack->get_LH() > (*old_BMStrack)->get_LH() ) || ( new_BMStrack->get_LH() > _reasonable_LH ) ){
                     (*old_BMStrack)->set_is_wrong();
                  }else{
                     new_BMStrack->set_is_wrong();
                     break;
                  }
               }else{
                  if( ( (*old_BMStrack)->get_LH() > new_BMStrack->get_LH() ) || ( ( new_BMStrack->get_nPlanes() != (*old_BMStrack)->get_nPlanes() ) && ( (*old_BMStrack)->get_LH() > _reasonable_LH ) ) ){
                     new_BMStrack->set_is_wrong();
                     break;
                  }else{
                     (*old_BMStrack)->set_is_wrong();
                  }
               }

            }
   
            //add track
            if(new_BMStrack->is_wrong()){
               delete new_BMStrack; 
               new_BMStrack = NULL;
            }else{
               (*BeamCombination)->add_bms_track(new_BMStrack);
               new_BMStrack = NULL;
            }
   
         } 
         }
         }
         }
         }
         }
      }

      //histogram
      new_BMStrack = (*BeamCombination)->get_best_BMStrack(&LH); 

      if(new_BMStrack != NULL){
         if( his_track_01[1] != NULL ) his_track_01[1]->Fill(new_BMStrack->get_nPlanes());
         if( his_track_02[1] != NULL ) his_track_02[1]->Fill(new_BMStrack->get_cluster_mask(), new_BMStrack->get_nPlanes());
         if( his_track_03[1] != NULL ) his_track_03[1]->Fill(new_BMStrack->get_t(), new_BMStrack->get_nPlanes());
         if( his_track_04[1] != NULL ) his_track_04[1]->Fill(new_BMStrack->get_mom(), new_BMStrack->get_nPlanes());
         if( his_track_05[1] != NULL ) his_track_05[1]->Fill(new_BMStrack->get_t()-((*BeamCombination)->get_BTtrack())->getMeanTime(), new_BMStrack->get_nPlanes());
         if( his_track_06[1] != NULL ) his_track_06[1]->Fill(TMath::Prob(new_BMStrack->get_tchi2(), new_BMStrack->get_nPlanes()-1), new_BMStrack->get_nPlanes());
         if( his_track_07[1] != NULL ) his_track_07[1]->Fill(TMath::Prob(new_BMStrack->get_schi2(), new_BMStrack->get_nPlanes()-1), new_BMStrack->get_nPlanes());
         if( his_track_08[1] != NULL ) his_track_08[1]->Fill(new_BMStrack->get_LH(), new_BMStrack->get_nPlanes());

         if( _his_mode >= 1 ){

            vector<CsBMScluster*> clusters = new_BMStrack->get_clusters();

            for(clus = clusters.begin(), i_plane = 0; clus != clusters.end(); clus++, i_plane++){

               if((*clus)->isEmpty()) continue;

               if( his_BMS_vs_BT0 != NULL ) his_BMS_vs_BT0->Fill((*clus)->get_t() - ((*BeamCombination)->get_BTtrack())->getMeanTime(), i_plane);
               if( his_BMS_vs_BT1 != NULL && TMath::Prob(new_BMStrack->get_tchi2(), new_BMStrack->get_nPlanes()-1) > 0.01 && TMath::Prob(new_BMStrack->get_schi2(), new_BMStrack->get_nPlanes()-1) > 0.01) his_BMS_vs_BT1->Fill((*clus)->get_t() - ((*BeamCombination)->get_BTtrack())->getMeanTime(), i_plane);
            }

         }

         if(new_BMStrack->get_nPlanes() == 6) histogram_residuals(new_BMStrack);

      }

   }

}

//=== HISTOGRAM RESIDUALS =======================================================
void CsBMSReconstruction :: histogram_residuals(CsBMStrack* BMStrack){

   if(BMStrack->get_nPlanes() != 6) return;

   int i, j, k;
   double a, b; 

   for(i = 0, j = -1, k = -1; i < 6; i++){
                     
      if(i == 0) { j = 4; k = 1; }
      if(i == 1) { j = 0; k = 4; }
      if(i == 2) { j = 5; k = 3; }
      if(i == 3) { j = 2; k = 5; } 
      if(i == 4) { j = 0; k = 1; }
      if(i == 5) { j = 2; k = 3; }

      a = ( (BMStrack->get_clusters()[k])->get_y() - (BMStrack->get_clusters()[j])->get_y() )/( _planeReference[k]->getZcm() - _planeReference[j]->getZcm());
      b =  (BMStrack->get_clusters()[k])->get_y() - a*_planeReference[k]->getZcm();

      if( his_residual_01[i] != NULL ) his_residual_01[i]->Fill((a*_planeReference[i]->getZcm() + b) - (BMStrack->get_clusters()[i])->get_y());
      if( his_residual_02[i] != NULL ) his_residual_02[i]->Fill( (BMStrack->get_clusters()[i])->get_y(), (a*_planeReference[i]->getZcm() + b) - (BMStrack->get_clusters()[i])->get_y());
      if( his_residual_03[i] != NULL ) his_residual_03[i]->Fill(1000*a, (a*_planeReference[i]->getZcm() + b) - (BMStrack->get_clusters()[i])->get_y());
   }

   return;
}

//=== SET SIZE AND DECODING MAP ======================================================
void CsBMSReconstruction :: set_BMS_decoding_map(){

   int ipl, ll, chn;

   //loop over planes
   for(ipl=0; ipl < 6; ipl++){

      //set size
      _planeSize[ipl]=64;
      if(ipl == 5) _planeSize[ipl] = 128;

      //create table which contains y position of each canal 
      _planeDecodingMap[ipl] = new double[_planeSize[ipl]];

      //fill table
      for(ll=1; ll <= _planeSize[ipl]; ll++){

         chn = ll-1;
         switch(ipl){
            case 0:  if    (ll<=4)              _planeDecodingMap[ipl][chn] = -87.5+5.*(ll-1)         + _planeReference[ipl]->getYcm(); 
                     else if((ll>4)&&(ll<=60))  _planeDecodingMap[ipl][chn] = -67.5+5.*((ll-5)/2)     + _planeReference[ipl]->getYcm();
                     else                       _planeDecodingMap[ipl][chn] =  72.5+5.*(ll-61)        + _planeReference[ipl]->getYcm();
                     break;
            case 1:  if    (ll<=8)              _planeDecodingMap[ipl][chn] = -42.5+5.*((ll-1)/2)     + _planeReference[ipl]->getYcm(); 
                     else if((ll>8 )&&(ll<=20)) _planeDecodingMap[ipl][chn] = -22.5+5.*((ll-9)/4)     + _planeReference[ipl]->getYcm();
                     else if((ll>20)&&(ll<=44)) _planeDecodingMap[ipl][chn] =  -7.5+5.*((ll-21)/6)    + _planeReference[ipl]->getYcm();
                     else if((ll>44)&&(ll<=56)) _planeDecodingMap[ipl][chn] =  12.5+5.*((ll-45)/4)    + _planeReference[ipl]->getYcm();
                     else                       _planeDecodingMap[ipl][chn] =  27.5+5.*((ll-57)/2)    + _planeReference[ipl]->getYcm();
                     break;
            case 2:  if    (ll<=12)             _planeDecodingMap[ipl][chn] = -47.5+5.*((ll-1)/2)     + _planeReference[ipl]->getYcm(); 
                     else if((ll>12)&&(ll<=20)) _planeDecodingMap[ipl][chn] = -17.5+5.*((ll-13)/4)    + _planeReference[ipl]->getYcm();
                     else if((ll>20)&&(ll<=44)) _planeDecodingMap[ipl][chn] =  -7.5+5.*((ll-21)/6)    + _planeReference[ipl]->getYcm();
                     else if((ll>44)&&(ll<=54)) _planeDecodingMap[ipl][chn] =  12.5+5.*((ll-45)/4)    + _planeReference[ipl]->getYcm();
                     else                       _planeDecodingMap[ipl][chn] =  27.5+5.*((ll-55)/2)    + _planeReference[ipl]->getYcm();
                     break;
            case 3:  if    (ll<=14)             _planeDecodingMap[ipl][chn] =-112.5+5.*(ll-1)         + _planeReference[ipl]->getYcm(); 
		     else if((ll>14)&&(ll<=50)) _planeDecodingMap[ipl][chn] = -42.5+5.*((ll-15)/2)    + _planeReference[ipl]->getYcm();
                     else                       _planeDecodingMap[ipl][chn] =  47.5+5.*(ll-51)        + _planeReference[ipl]->getYcm();
                     break;
	    case 4:  _planeDecodingMap[ipl][chn] = (chn - 27.5) * 2.5                                 + _planeReference[ipl]->getYcm();
                     break;
	    case 5:  _planeDecodingMap[ipl][chn] = (chn - 62.5) * 1.25                                + _planeReference[ipl]->getYcm();
                     break;
         }
      }
   }
   
   //position of slabs inside BMS01-04 crates along beam direction w.r.t plane center in mm
   //according to 
   //[channel: 1-64][plane: BM01-04]
   double planeZcorrection[64][4] = {
   
      {	   -50.0,	   -90.0,	   110.0,	   -50.0	},
      {	   -10.0,	    90.0,	  -110.0,	    90.0	},
      {	    50.0,	   -50.0,	   -90.0,	    30.0	},
      {	    10.0,	    50.0,	    90.0,	    10.0	},
      {	    30.0,	   -10.0,	   -50.0,	   110.0	},
      {	   -30.0,	    10.0,	    50.0,	   -70.0	},
      {	    70.0,	    30.0,	   -10.0,	   -90.0	},
      {	   -70.0,	   -30.0,	    10.0,	    50.0	},
      {	   110.0,	    70.0,	    30.0,	   -10.0	},
      {	  -110.0,	   -50.0,	   -30.0,	   -30.0	},
      {	   -90.0,	    50.0,	    70.0,	    70.0	},
      {	    90.0,	   -70.0,	   -70.0,	  -110.0	},
      {	   -50.0,	   110.0,	   -90.0,	   -90.0	},
      {	    50.0,	   -10.0,	   -10.0,	    50.0	},
      {	   -10.0,	    10.0,	    10.0,	   -10.0	},
      {	    10.0,	  -110.0,	    90.0,	    10.0	},
      {	    30.0,	   -50.0,	   -50.0,	    30.0	},
      {	   -30.0,	    30.0,	    30.0,	   -30.0	},
      {	    70.0,	   -30.0,	   -30.0,	    70.0	},
      {	   -70.0,	    50.0,	    50.0,	   -70.0	},
      {	   110.0,	   -10.0,	   -10.0,	   110.0	},
      {	  -110.0,	    70.0,	    70.0,	  -110.0	},
      {	   -90.0,	   -90.0,	   -90.0,	   -90.0	},
      {	    90.0,	    90.0,	    90.0,	    90.0	},
      {	   -50.0,	   -70.0,	   -70.0,	   -50.0	},
      {	    50.0,	    10.0,	    10.0,	    50.0	},
      {	   -10.0,	    30.0,	    30.0,	   -10.0	},
      {	    10.0,	   110.0,	   110.0,	    10.0	},
      {	    30.0,	   -50.0,	   -50.0,	    30.0	},
      {    -30.0,	    50.0,	    50.0,	   -30.0	},
      {    110.0,	  -110.0,	  -110.0,	   110.0	},
      {	  -110.0,	   -30.0,	   -30.0,	  -110.0	},
      {	    70.0,	   -10.0,	   -10.0,	    70.0	},
      {	   -70.0,	    70.0,	    70.0,	   -70.0	},
      {	    30.0,	   -90.0,	   -90.0,	    30.0	},
      {	   -30.0,	    90.0,	    90.0,	   -30.0	},
      {	   -10.0,	   -70.0,	   -70.0,	   -10.0	},
      {	    10.0,	    10.0,	    10.0,	    10.0	},
      {	   -50.0,	   -50.0,	   -50.0,	   -50.0	},
      {	    50.0,	    30.0,	    30.0,	    50.0	},
      {	   -90.0,	   110.0,	   110.0,	   -90.0	},
      {	    90.0,	  -110.0,	  -110.0,	    90.0	},
      {	   110.0,	   -30.0,	   -30.0,	   110.0	},
      {	  -110.0,	    50.0,	    50.0,	  -110.0	},
      {	    70.0,	   -90.0,	   -90.0,	    70.0	},
      {	   -70.0,	    70.0,	    70.0,	   -70.0	},
      {	    30.0,	   -70.0,	   -70.0,	    30.0	},
      {	   -30.0,	    90.0,	    90.0,	   -30.0	},
      {	   -10.0,	   110.0,	   -10.0,	   -10.0	},
      {	    10.0,	    30.0,	    30.0,	    10.0	},
      {	   -50.0,	   -30.0,	   -30.0,	    50.0	},
      {	    50.0,	  -110.0,	    10.0,	   -90.0	},
      {	   -90.0,	    70.0,	    70.0,	  -110.0	},
      {	    90.0,	   -10.0,	   -70.0,	    70.0	},
      {	   110.0,	    10.0,	    30.0,	   -30.0	},
      {	  -110.0,	   -70.0,	   -30.0,	   -10.0	},
      {	    70.0,	    30.0,	   -10.0,	    50.0	},
      {	   -70.0,	   -30.0,	    10.0,	   -90.0	},
      {	    30.0,	   -10.0,	   -50.0,	   -70.0	},
      {	   -30.0,	    10.0,	    50.0,	   110.0	},
      {	    10.0,	   -50.0,	   -90.0,	    10.0	},
      {	    50.0,	    50.0,	    90.0,	    30.0	},
      {	   -10.0,	   -90.0,	   110.0,	    90.0	},
      {	   -50.0,	    90.0,	  -110.0,	   -50.0	}
   };

   //copy values
   for(ipl = 0; ipl < 4; ipl++){
      for(chn = 0; chn < 64; chn++){
         _planeZcorrection[chn][ipl] = planeZcorrection[chn][ipl];
      }
   }
   
}

//=== SET HISTOGRAMS ============================================================
void CsBMSReconstruction :: set_histograms(){

   //"first" marker
   static bool first = true;
   char title[100], name[100];

   if(first){

      for(int i = 0; i < 6; i++) his_cluster_01[i] = NULL;
      for(int i = 0; i < 6; i++) his_cluster_02[i] = NULL;
      for(int i = 0; i < 6; i++) his_cluster_03[i] = NULL;
      his_cluster_04 = NULL;
      his_cluster_05 = NULL;
      his_cluster_06 = NULL;

      his_BMS_vs_BT0 = NULL;
      his_BMS_vs_BT1 = NULL;

      for(int i = 0; i < 6; i++) his_residual_01[i] = NULL;
      for(int i = 0; i < 6; i++) his_residual_02[i] = NULL;
      for(int i = 0; i < 6; i++) his_residual_03[i] = NULL;
      for(int i = 0; i < 6; i++) his_hit_time_vs_trigger[i] = NULL;

      for(int i = 0; i < 2; i++) his_track_01[i] = NULL;
      for(int i = 0; i < 2; i++) his_track_02[i] = NULL;
      for(int i = 0; i < 2; i++) his_track_03[i] = NULL;
      for(int i = 0; i < 2; i++) his_track_04[i] = NULL;
      for(int i = 0; i < 2; i++) his_track_05[i] = NULL;
      for(int i = 0; i < 2; i++) his_track_06[i] = NULL;
      for(int i = 0; i < 2; i++) his_track_07[i] = NULL;
      for(int i = 0; i < 2; i++) his_track_08[i] = NULL;

      string histograms_path =  "/BeamReconstruction/BMSReconstruction";
      CsHistograms::SetCurrentPath(histograms_path);

      for(int ipl = 0; ipl < 6; ipl++){

         //hit channel vs. time
         sprintf(title, "BM0%i, hit channel vs. time BMS", ipl+1);
         sprintf(name,  "cluster_01_BM0%i", ipl+1);
         if( _his_mode > 1 ) his_cluster_01[ipl] = new CsHist2D(name, title, _planeSize[ipl], 0, _planeSize[ipl], 60, -15, 15);

         //hit position vs. time
         sprintf(title, "BM0%i, hit position vs. time BMS", ipl+1);
         sprintf(name,  "cluster_02_BM0%i", ipl+1);
         if( _his_mode > 1 ) his_cluster_02[ipl] = new CsHist2D(name, title, 100, _planeDecodingMap[ipl][0], _planeDecodingMap[ipl][_planeSize[ipl]-1], 60, -15, 15);

         //cluster position vs. mean time
         sprintf(title, "BM0%i, cluster position vs. mean time BMS", ipl+1);
         sprintf(name,  "cluster_03_BM0%i", ipl+1);
         if( _his_mode > 1 ) his_cluster_03[ipl] = new CsHist2D(name, title, 100, _planeDecodingMap[ipl][0], _planeDecodingMap[ipl][_planeSize[ipl]-1], 60, -15, 15);

         //residual
         sprintf(title, "BM0%i, residuals - final tracks", ipl+1);
         sprintf(name,  "residual_01_BM0%i", ipl+1);
         if( _his_mode > 0 ) his_residual_01[ipl] = new CsHist1D(name, title, 100, -15, 15);

         //residual
         sprintf(title, "BM0%i, residuals vs. y - final tracks", ipl+1);
         sprintf(name,  "residual_02_BM0%i", ipl+1);
         if( _his_mode > 0 ) his_residual_02[ipl] = new CsHist2D(name, title, 100, _planeDecodingMap[ipl][0], _planeDecodingMap[ipl][_planeSize[ipl]-1], 100, -15, 15);

         //residual
         sprintf(title, "BM0%i, residuals vs. y - final tracks", ipl+1);
         sprintf(name,  "residual_03_BM0%i", ipl+1);
         if( _his_mode > 0 ) his_residual_03[ipl] = new CsHist2D(name, title, 100, -5., 5., 100, -15, 15);

         //hit time vs. trigger
         sprintf(title, "BM0%i, hit time vs. trigger", ipl+1);
         sprintf(name,  "hit_time_vs_trigger_BM0%i", ipl+1);
         if( _his_mode > 1 ) his_hit_time_vs_trigger[ipl] = new CsHist2D(name, title, 600, -30, 30, 12, 0, 12);
      }

      //BMS clusters vs. BT track clusters
      if( _his_mode > 0 ) his_BMS_vs_BT0 = new CsHist2D("his_BMS_vs_BT0", "channel time vs. BT track", 200, -5, 5, 6, 0, 6);

      //BMS clusters vs. BT track clusters
      if( _his_mode > 1 ) his_BMS_vs_BT1 = new CsHist2D("his_BMS_vs_BT1", "channel time vs. BT track, tLH, sLH > 0.01", 200, -5, 5, 6, 0, 6);

      //number of hits/event vs. plane
      if( _his_mode > 0 ) his_cluster_04 = new CsHist2D("cluster_04", "number of hits/event vs. plane", 100, 0, 100, 6, 0, 6);

      //number of clusters/event vs. plane
      if( _his_mode > 1 ) his_cluster_05 = new CsHist2D("cluster_05", "number of clusters/event vs. plane", 100, 0, 100, 6, 0, 6);

      //number of hits/cluster vs. plane
      if( _his_mode > 1 ) his_cluster_06 = new CsHist2D("cluster_06", "number of hits/cluster vs. plane", 50, 0, 50, 6, 0, 6);


      //number of fired planes/track
      if( _his_mode > 1 ) his_track_01[0] = new CsHist1D("track_01_all",     "number of fired planes/track (all tracks)",    7, 0, 7);
      if( _his_mode > 0 ) his_track_01[1] = new CsHist1D("track_01_final",   "number of fired planes/track (final tracks)",  7, 0, 7);

      //fired planes vs. track class
      if( _his_mode > 1 ) his_track_02[0] = new CsHist2D("track_02_all",     "fired planes (binary) vs. track class (all tracks)",    65, 0, 65, 7, 0, 7);
      if( _his_mode > 0 ) his_track_02[1] = new CsHist2D("track_02_final",   "fired planes (binary) vs. track class (final tracks)",  65, 0, 65, 7, 0, 7);

      //time vs. track class
      if( _his_mode > 1 ) his_track_03[0] = new CsHist2D("track_03_all",     "time vs. track class (all tracks)",   60, -15, 15, 7, 0, 7);
      if( _his_mode > 0 ) his_track_03[1] = new CsHist2D("track_03_final",   "time vs. track class (final tracks)", 60, -15, 15, 7, 0, 7);

      //momentum vs. track class
      if( _his_mode > 1 ) his_track_04[0] = new CsHist2D("track_04_all",     "momentum vs. track class (all tracks)",     80, 140, 180, 7, 0, 7);
      if( _his_mode > 0 ) his_track_04[1] = new CsHist2D("track_04_final",   "momentum vs. track class (final tracks)",   80, 140, 180, 7, 0, 7);

      //dt vs. track class
      if( _his_mode > 1 ) his_track_05[0] = new CsHist2D("track_05_all",     "dt (tBMS-tBT) vs. track class (all tracks)",   60, -15, 15, 7, 0, 7);
      if( _his_mode > 0 ) his_track_05[1] = new CsHist2D("track_05_final",   "dt (tBMS-tBT) vs. track class (final tracks)", 60, -15, 15, 7, 0, 7);

      //tchi2 vs. track class
      if( _his_mode > 1 ) his_track_06[0] = new CsHist2D("track_06_all",     "LH(time chi2, ndf) vs. track class (all tracks)",    100, 0, 1, 7, 0, 7);
      if( _his_mode > 0 ) his_track_06[1] = new CsHist2D("track_06_final",   "LH(time chi2, ndf) vs. track class (final tracks)",  100, 0, 1, 7, 0, 7);

      //schi2 vs. track class
      if( _his_mode > 1 ) his_track_07[0] = new CsHist2D("track_07_all",     "LH(spatial chi2, ndf) vs. track class (all tracks)",    100, 0, 1, 7, 0, 7);
      if( _his_mode > 0 ) his_track_07[1] = new CsHist2D("track_07_final",   "LH(spatial chi2, ndf) vs. track class (final tracks)",  100, 0, 1, 7, 0, 7);

      //chi2 vs. track class
      if( _his_mode > 1 ) his_track_08[0] = new CsHist2D("track_08_all",     "LH(total chi2, ndf) vs. track class (all tracks)",   100, 0, 1, 7, 0, 7);
      if( _his_mode > 0 ) his_track_08[1] = new CsHist2D("track_08_final",   "LH(total chi2, ndf) vs. track class (final tracks)", 100, 0, 1, 7, 0, 7);

      first = false;
   }
}

//=== INSTANCE ====================================================================
CsBMSReconstruction* CsBMSReconstruction :: instance_ = NULL;

CsBMSReconstruction* CsBMSReconstruction :: Instance() {
   if( instance_ == NULL ){
      instance_ = new CsBMSReconstruction();
   }
       
   return( instance_ );
}

