#include "CsBeamReconstructionSupport.h"

using namespace std;

void read_configuration_key(string tag, string key, double& value){

   CsOpt* opt = CsOpt::Instance();

   if(! opt->getOpt(tag, key, value) ){
      cout << "BEAM RECONSTRUCTION: option not set, tag '" << tag << "', key '" << key << "'" << endl;
      cout << "BEAM RECONSTRUCTION: are you using beam option file appropriate for this reconstruction model?" << endl;
      CsErrLog::Instance()->mes(elFatal, "Exiting ...");
   }
}

void read_configuration_key(string tag, string key, int& value){

   CsOpt* opt = CsOpt::Instance();

   if(! opt->getOpt(tag, key, value) ){
      cout << "BEAM RECONSTRUCTION: option not set, tag '" << tag << "', key '" << key << "'" << endl;
      cout << "BEAM RECONSTRUCTION: are you using beam option file appropriate for this reconstruction model?" << endl;
      CsErrLog::Instance()->mes(elFatal, "Exiting ...");
   }
}

void read_configuration_key(string tag, string key, bool& value){

   CsOpt* opt = CsOpt::Instance();

   if(! opt->getOpt(tag, key, value) ){
      cout << "BEAM RECONSTRUCTION: option not set, tag '" << tag << "', key '" << key << "'" << endl;
      cout << "BEAM RECONSTRUCTION: are you using beam option file appropriate for this reconstruction model?" << endl;
      CsErrLog::Instance()->mes(elFatal, "Exiting ...");
   }
}

void read_configuration_key(string tag, string key, string& value){

   CsOpt* opt = CsOpt::Instance();

   if(! opt->getOpt(tag, key, value) ){
      cout << "BEAM RECONSTRUCTION: option not set, tag '" << tag << "', key '" << key << "'" << endl;
      cout << "BEAM RECONSTRUCTION: are you using beam option file appropriate for this reconstruction model?" << endl;
      CsErrLog::Instance()->mes(elFatal, "Exiting ...");
   }
}

void CsBMScluster::set_cs_cluster(CsDetector* det, std::vector<CsDigit*> digits){

   double u = _y/10.;
   double v = 0.;
   double w = _z/10.;

   CLHEP::HepMatrix cov(3,3,0);
   cov(1, 1) = pow(_dy, 2);
   cov(2, 2) = 1.E6; 

   _cs_cluster = new CsCluster(u, v, w, cov);
   _cs_cluster->setTime(_t);
   _cs_cluster->setTimeError(_dt);

   _cs_cluster->addDet(*det);

   vector<CsDigit*>::iterator digit;
   for(digit = digits.begin(); digit != digits.end(); digit++) _cs_cluster->addDigit(*(*digit));

   CsEvent::Instance()->addCluster(*_cs_cluster);

   return;
}

CsBMScluster :: CsBMScluster(){
   _is_empty = true;
   _t       = 0.;
   _dt      = 0.;
   _y       = 0.;
   _dy      = 0.;
   _z       = 0.;
   _cs_cluster = NULL;
}

CsBMScluster :: CsBMScluster(double t, double dt, double y, double dy, double z){
   _is_empty = false;
   _t       = t;
   _dt      = dt;
   _y       = y;
   _dy      = dy;
   _z       = z;
   _cs_cluster = NULL;
}

CsBMScluster :: ~CsBMScluster(){
}

void CsBMStrack :: calculate_time(double BT_time, double* BTvsBMSResolution, double* BTvsBMSCorrection){

   int i;
   double weight;
   double sum_of_weights;

   //time 
   for(i = 0, _t = 0., sum_of_weights = 0.; i < 6; i++){

      if( !(*_clusters[i]).isEmpty() ){
         weight = 1./pow(BTvsBMSResolution[i]/2., 2);

         _t += (( (*_clusters[i]).get_t() + BT_time )/2.) * weight;
         sum_of_weights += weight;
      }
   }
   _t /= sum_of_weights;

   //error of the mean value
   _dt = 1./sqrt(sum_of_weights);

   //time chi2
   for(i = 0, _tchi2 = 0.; i < 6; i++){
      if( !(*_clusters[i]).isEmpty() ) _tchi2 += pow(((*_clusters[i]).get_t() - BT_time - BTvsBMSCorrection[i]) / BTvsBMSResolution[i], 2.);
   }

}

void CsBMStrack :: show_clusters(){

      cout << "show_clusters" << "\t";
      cout << CsEvent::Instance()->getRunNumber() << "\t";
      cout << CsEvent::Instance()->getBurstNumber() << "\t";
      cout << CsEvent::Instance()->getEventNumberInBurst() << "\t";
      cout << "BM01: "; if((_clusters[0])->isEmpty()){cout << "empty";} else cout << (_clusters[0])->get_y(); cout << "\t";
      cout << "BM05: "; if((_clusters[4])->isEmpty()){cout << "empty";} else cout << (_clusters[4])->get_y(); cout << "\t";
      cout << "BM02: "; if((_clusters[1])->isEmpty()){cout << "empty";} else cout << (_clusters[1])->get_y(); cout << "\t";
      cout << "BM03: "; if((_clusters[2])->isEmpty()){cout << "empty";} else cout << (_clusters[2])->get_y(); cout << "\t";
      cout << "BM06: "; if((_clusters[5])->isEmpty()){cout << "empty";} else cout << (_clusters[5])->get_y(); cout << "\t";
      cout << "BM04: "; if((_clusters[3])->isEmpty()){cout << "empty";} else cout << (_clusters[3])->get_y(); cout << endl;

      return;
}

int CsBMStrack :: get_cluster_mask(){
   
   int word = 0;

   for(int i_plane = 0; i_plane < 6; i_plane++){
      if(!((_clusters[i_plane])->isEmpty())) word = (word | 1<<i_plane);
   }

   return word;
}

CsBMStrack :: CsBMStrack(int nPlanes, CsBMScluster* hit1, CsBMScluster* hit2, CsBMScluster* hit3, 
                                      CsBMScluster* hit4, CsBMScluster* hit5, CsBMScluster* hit6){

   _nPlanes = nPlanes;

   _clusters.clear();
   _clusters.push_back(hit1);
   _clusters.push_back(hit2);
   _clusters.push_back(hit3);
   _clusters.push_back(hit4);
   _clusters.push_back(hit5);
   _clusters.push_back(hit6);

   _t = _dt = _mom = _tchi2 = _schi2 = 0.;
   _is_wrong = false;
}

CsBMStrack :: ~CsBMStrack(){};

void CsBMSslice :: add_cluster(int ipl, CsBMScluster* cluster){

   if(ipl == 0){  _clustersBM01.push_back(cluster);   return;  }
   if(ipl == 1){  _clustersBM02.push_back(cluster);   return;  }
   if(ipl == 2){  _clustersBM03.push_back(cluster);   return;  }
   if(ipl == 3){  _clustersBM04.push_back(cluster);   return;  }
   if(ipl == 4){  _clustersBM05.push_back(cluster);   return;  }
   if(ipl == 5){  _clustersBM06.push_back(cluster);   return;  }

   cout << "CsBMSslice::add_clusters() Error: Wrong index of plane" << endl;
   cout << "CsBMSslice::add_clusters() Index: " << ipl << endl;
   CsErrLog::Instance()->mes(elFatal, "Exiting ...");

   return;
}

const vector<CsBMScluster*> CsBMSslice :: get_clusters(int ipl){

   if(ipl == 0){ return _clustersBM01; }
   if(ipl == 1){ return _clustersBM02; }
   if(ipl == 2){ return _clustersBM03; }
   if(ipl == 3){ return _clustersBM04; }
   if(ipl == 4){ return _clustersBM05; }
   if(ipl == 5){ return _clustersBM06; }

   cout << "CsBMSslice::get_clusters() Error: Wrong index of plane" << endl;
   cout << "CsBMSslice::get_clusters() Index: " << ipl << endl;
   CsErrLog::Instance()->mes(elFatal, "Exiting ...");

   return _clustersBM01;
}

CsBMSslice :: CsBMSslice(double* y_min, double* y_max){

   _clustersBM01.clear();
   _clustersBM02.clear();
   _clustersBM03.clear();
   _clustersBM04.clear();
   _clustersBM05.clear();
   _clustersBM06.clear();

   for(int i_plane = 0; i_plane < 6; i_plane++){
      _y_min[i_plane] = y_min[i_plane];
      _y_max[i_plane] = y_max[i_plane];
   }
}

CsBMSslice :: ~CsBMSslice(){};

CsBMStrack* CsBeamCombination :: get_best_BMStrack(double* best_LH){

   CsBMStrack* best_BMStrack = NULL;
   *best_LH = 0.;

   vector<CsBMStrack*>::iterator BMStrack;

   //loop over reconstructed BMS tracks
   for(BMStrack = _BMStracks.begin(); BMStrack != _BMStracks.end(); BMStrack++){

      //check if wrong
      if( (*BMStrack)->is_wrong() ) continue;

      //check if the best
      if((*BMStrack)->get_LH() > *best_LH){   best_BMStrack = *BMStrack; *best_LH = (*BMStrack)->get_LH(); }
   }

   return best_BMStrack;
}

bool CsBeamCombination :: check_time(double time){

   if(time > _time_min && time < _time_max) return true;
   return false;
}

void CsBeamCombination :: add_cluster(int i_plane, CsBMScluster* cluster){

   if(! check_time(cluster->get_t()) ) return;

   bool is_added = false;
   vector<CsBMSslice*>::iterator BMSslice;

   for(BMSslice = _BMSslices.begin(); BMSslice != _BMSslices.end(); BMSslice++){

      if( 
         ( cluster->get_y() < (*BMSslice)->get_y_min(i_plane) ) ||
         ( cluster->get_y() > (*BMSslice)->get_y_max(i_plane) )  
      ) continue;

      (*BMSslice)->add_cluster(i_plane, cluster);
      is_added = true;
   }

   if( ! is_added ){
      CsErrLog::Instance()->msg(elWarning, __FILE__, __LINE__, "No BMS slice for given BMS cluster"); 
   }

   return;
}

CsBeamCombination :: CsBeamCombination(const CsTrack* BTtrack, double time_gate, double slice_mom, int slice_n, double slice_overlap, CsBMScluster* empty_cluster){

   _BTtrack = BTtrack;
   _BMStracks.clear();

   _BMSslices.clear();

   _time_min = _BTtrack->getMeanTime() - time_gate;
   _time_max = _BTtrack->getMeanTime() + time_gate;

   _slice_mom = slice_mom;
   _slice_n = slice_n;
   _slice_overlap = slice_overlap;

   BeamBackPropagation = CsBeamBackPropagation :: Instance();

   //find helix used in back propagation
   double last_helix_position;

   vector<CsHelix>helices = _BTtrack->getHelices();
   vector<CsHelix>::iterator helix;
   CsHelix *last_helix;

   for(helix = helices.begin(), last_helix_position = -1.E6, last_helix = NULL; helix != helices.end(); helix++){

      if( (*helix).getZ() > last_helix_position ){
         last_helix_position = (*helix).getZ();
         last_helix = &(*helix);
      }
   }

   if( last_helix == NULL ){
      CsErrLog::Instance()->msg(elFatal, __FILE__, __LINE__, "Beam telescope track does't have any helices"); 
   }else{
      _lastBThelix = new CsHelix(*last_helix);
   }

   //create slices
   int i_slice, i_plane;
   double mom;
   double y_min[6];
   double y_max[6];
   double y_min_final[6];
   double y_max_final[6];
   CsBMSslice *new_BMSslice;

   // < min
   mom = 120.;
   BeamBackPropagation->make_back_propagation(_lastBThelix, mom);
   for(i_plane = 0; i_plane < 6; i_plane++) y_min[i_plane] = (BeamBackPropagation->get_y(i_plane) + BeamBackPropagation->get_correctionFactor(i_plane));

   mom = 160.-_slice_mom;
   BeamBackPropagation->make_back_propagation(_lastBThelix, mom);
   for(i_plane = 0; i_plane < 6; i_plane++) y_max[i_plane] = (BeamBackPropagation->get_y(i_plane) + BeamBackPropagation->get_correctionFactor(i_plane));

   for(i_plane = 0; i_plane < 6; i_plane++){
      y_min_final[i_plane] = y_min[i_plane] - _slice_overlap*BeamBackPropagation->get_spaceResolutionBMS(i_plane);
      y_max_final[i_plane] = y_max[i_plane] + _slice_overlap*BeamBackPropagation->get_spaceResolutionBMS(i_plane);
   }

   new_BMSslice = new CsBMSslice(y_min_final, y_max_final);
   _BMSslices.push_back(new_BMSslice);
   new_BMSslice = NULL; 

   // min - max
   for(i_slice = 1; i_slice <= _slice_n; i_slice++){
      
      for(i_plane = 0; i_plane < 6; i_plane++) y_min[i_plane] = y_max[i_plane];

      mom = (160.-_slice_mom) + i_slice*(2.*_slice_mom/(double)_slice_n);
      BeamBackPropagation->make_back_propagation(_lastBThelix, mom);
      for(i_plane = 0; i_plane < 6; i_plane++) y_max[i_plane] = (BeamBackPropagation->get_y(i_plane) + BeamBackPropagation->get_correctionFactor(i_plane));

      for(i_plane = 0; i_plane < 6; i_plane++){
         y_min_final[i_plane] = y_min[i_plane] - _slice_overlap*BeamBackPropagation->get_spaceResolutionBMS(i_plane);
         y_max_final[i_plane] = y_max[i_plane] + _slice_overlap*BeamBackPropagation->get_spaceResolutionBMS(i_plane);
      }

      new_BMSslice = new CsBMSslice(y_min_final, y_max_final);
      _BMSslices.push_back(new_BMSslice);
      new_BMSslice = NULL; 
   }

   // > max
   for(i_plane = 0; i_plane < 6; i_plane++) y_min[i_plane] = y_max[i_plane];

   mom = 200.;
   BeamBackPropagation->make_back_propagation(_lastBThelix, mom);
   for(i_plane = 0; i_plane < 6; i_plane++) y_max[i_plane] = (BeamBackPropagation->get_y(i_plane) + BeamBackPropagation->get_correctionFactor(i_plane));

   for(i_plane = 0; i_plane < 6; i_plane++){
      y_min_final[i_plane] = y_min[i_plane] - _slice_overlap*BeamBackPropagation->get_spaceResolutionBMS(i_plane);
      y_max_final[i_plane] = y_max[i_plane] + _slice_overlap*BeamBackPropagation->get_spaceResolutionBMS(i_plane);
   }

   new_BMSslice = new CsBMSslice(y_min_final, y_max_final);
   _BMSslices.push_back(new_BMSslice);
   new_BMSslice = NULL;

   //add empty cluster 
   vector<CsBMSslice*>::iterator BMSslice;

   for(BMSslice = _BMSslices.begin(); BMSslice != _BMSslices.end(); BMSslice++){
      for(i_plane = 0; i_plane < 6; i_plane++) (*BMSslice)->add_cluster(i_plane, empty_cluster);
   }
}

CsBeamCombination :: ~CsBeamCombination(){

   delete _lastBThelix;
   _lastBThelix = NULL;

   vector<CsBMSslice*>::iterator BMSslice;
   for(BMSslice = _BMSslices.begin(); BMSslice != _BMSslices.end(); BMSslice++){ delete *BMSslice; *BMSslice = NULL;}
   _BMSslices.clear();

   vector<CsBMStrack*>::iterator BMStrack;
   for(BMStrack = _BMStracks.begin(); BMStrack != _BMStracks.end(); BMStrack++){ delete *BMStrack; *BMStrack = NULL;}
   _BMStracks.clear();
}
