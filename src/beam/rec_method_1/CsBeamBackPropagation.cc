#include "CsBeamBackPropagation.h"

//=== CONSTRUCTOR ====================================================================
CsBeamBackPropagation :: CsBeamBackPropagation(){

   //few variables
   int i;

   //read configuration file 
   read_configuration_key("beam", "histogram_level", _his_mode);

   read_configuration_key("beam_bms", "BM01_name", _planeName[0]);
   read_configuration_key("beam_bms", "BM02_name", _planeName[1]);
   read_configuration_key("beam_bms", "BM03_name", _planeName[2]);
   read_configuration_key("beam_bms", "BM04_name", _planeName[3]);
   read_configuration_key("beam_bms", "BM05_name", _planeName[4]);
   read_configuration_key("beam_bms", "BM06_name", _planeName[5]);
 
   //beam back propagation coefficients
   _readFromAligFile = false;
   CsOpt::Instance()->getOpt("beam_backprop", "ReadFromAligFile", _readFromAligFile);

   if( _readFromAligFile ) { 
      std::cout << "CsBeamBackPropagation.cc: Info: Values of the beam back propagation algorithm will be taken from alignment file" << std::endl;

      //find matched entry
      std::vector<CsBeamBackPropagationCoeff> bbp = CsGeom::Instance()->getBeamBackProp();
      std::vector<CsBeamBackPropagationCoeff>::iterator ibbp;

      for(i = 0; i < 6; i++) {
         
         int howMany = 0;   

         for(ibbp = bbp.begin(); ibbp != bbp.end(); ibbp++) {
            if( ibbp->get_plane_name() == _planeName[i] ){

               _correctionFactor[i] = ibbp->get_correctionFactor();
               _spaceResolutionBMS[i] = ibbp->get_spaceResolutionBMS();
               howMany++;
            }
         }

         if( howMany == 0 )
            CsErrLog::msg(elFatal,__FILE__,__LINE__,
            "CsBeamBackPropagation.cc: Error: No entry in alignment file, BM0%d", i+1);

         if( howMany > 1 )
            CsErrLog::msg(elFatal,__FILE__,__LINE__,
            "CsBeamBackPropagation.cc: Error: Two entries with the same plane name in alignment file, BM0%d", i+1);
      }


   } else {
      std::cout << "CsBeamBackPropagation.cc: Info: Values of the beam back propagation algorithm will be taken from option files" << std::endl;

      read_configuration_key("beam_backprop", "space_resolution_BM01", _spaceResolutionBMS[0]);
      read_configuration_key("beam_backprop", "space_resolution_BM02", _spaceResolutionBMS[1]);
      read_configuration_key("beam_backprop", "space_resolution_BM03", _spaceResolutionBMS[2]);
      read_configuration_key("beam_backprop", "space_resolution_BM04", _spaceResolutionBMS[3]);
      read_configuration_key("beam_backprop", "space_resolution_BM05", _spaceResolutionBMS[4]);
      read_configuration_key("beam_backprop", "space_resolution_BM06", _spaceResolutionBMS[5]);

      read_configuration_key("beam_backprop", "space_correction_BM01", _correctionFactor[0]);
      read_configuration_key("beam_backprop", "space_correction_BM02", _correctionFactor[1]);
      read_configuration_key("beam_backprop", "space_correction_BM03", _correctionFactor[2]);
      read_configuration_key("beam_backprop", "space_correction_BM04", _correctionFactor[3]);
      read_configuration_key("beam_backprop", "space_correction_BM05", _correctionFactor[4]);
      read_configuration_key("beam_backprop", "space_correction_BM06", _correctionFactor[5]);
   }

   //print
   std::cout << "CsBeamBackPropagation.cc: Info: Used bbp coefficients" << std::endl;
   for(i = 0; i < 6; i++) {
      printf("\tBM0%d\t% 8.4f\t% 8.4f\n", i+1,  _correctionFactor[i], _spaceResolutionBMS[i]);
   }

   //Set pointers to detectors
   std::list<CsDetector*>  det = CsGeom::Instance()->getDetectors();
   std::list<CsDetector*>::iterator idet;

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
         std::cout << "CsBMSReconstruction.cc: Error: Reference to BMS plane is not existed" << std::endl;
         std::cout << "CsBMSReconstruction.cc: Plane name = " << _planeName[i] << std::endl;
         CsErrLog::Instance()->mes(elFatal, "Exiting ...");
      }
   }

   //set histograms
   set_histograms();
   
   //clear
   clear();
}

//=== DESTRUCTOR ======================================================================
CsBeamBackPropagation :: ~CsBeamBackPropagation(){
   clear();
}

//=== CALCULATE POSITIONS =============================================================
bool CsBeamBackPropagation :: make_back_propagation(CsHelix* BTHelix, double mom){

   CsHelix exHel;
   double SF_y, SF_DYDZ, SF_y_err, SF_DYDZ_err;
   double SF_x, SF_DXDZ, SF_x_err, SF_DXDZ_err;

   //Tarnsport equations are defined for x and y known at 0 of Lau coordinate system
   //Lau coordinate system is shifted by 500 cm from compass coordinate system 
   //(calculated comparing detectors.dat and printout of transport equations parameters)
   //We take helix closest to z=0 abd extrapolate from it to z_lau=0
   //Then we do some conversions (Ask Martin why...)

   if( !(BTHelix->Extrapolate(-5220, exHel)) ) return false;

   SF_y = exHel.getY();
   SF_y_err = sqrt(exHel.getCov()[2]);
   SF_DYDZ = - exHel.getDYDZ()*1000;
   SF_DYDZ_err = sqrt(exHel.getCov()[9])*1000;

   SF_x = exHel.getX();
   SF_x_err = sqrt(exHel.getCov()[0]);
   SF_DXDZ = - exHel.getDXDZ()*1000;
   SF_DXDZ_err = sqrt(exHel.getCov()[5])*1000;

   double BMS_dpp = (mom - 160.)/1.6;
   double BMS_dpp_err = (0.005 * mom)/1.6;

   double R11[4];
   double R12[4];
   double R21[4];
   double R22[4];
   double R33[4];
   double R34[4];
   double R36[4];
   double R43[4];
   double R44[4];
   double R46[4];
   //matrix elements:   R11    R12       R21     R22       R33    R34       R43      R44       R16      R26     R36     R46
   //BMS1           *   -0.741 -20.395   0.055   0.170 *    0.561 -43.498   0.044  -1.640 *    0.000   0.000   13.596   0.538
   //BMS2           *   -1.482 -22.671   0.055   0.170 *   -0.032 -21.484   0.044  -1.640 *    0.000   0.000    6.376   0.538
   //BMS3           *    0.068  50.581  -0.018   1.580 *   -0.894  16.304   0.032  -1.698 *    0.000   0.000    0.610  -0.037
   //BMS4           *    0.287  30.947  -0.018   1.580 *   -1.288  37.401   0.032  -1.698 *    0.000   0.000    1.076  -0.037
   R11[0]=-0.741;
   R11[1]=-1.482;
   R11[2]=0.068;
   R11[3]=0.287;

   R12[0]=-20.395;
   R12[1]=-22.671;
   R12[2]=50.581;
   R12[3]=30.947;

   R21[0]=0.055;
   R21[1]=0.055;
   R21[2]=-0.018;
   R21[3]=-0.018;

   R22[0]=0.170;
   R22[1]=0.170;
   R22[2]=1.580;
   R22[3]=1.580;

   R33[0]=0.561;
   R33[1]=-0.032;
   R33[2]=-0.894;
   R33[3]=-1.288;

   R34[0]=-43.498;
   R34[1]=-21.484;
   R34[2]=16.304;
   R34[3]=37.401 ;

   R36[0]=13.596;
   R36[1]=6.376;
   R36[2]=0.610;
   R36[3]=1.076;

   R43[0]=0.044;
   R43[1]=0.044;
   R43[2]=0.032;
   R43[3]=0.032;

   R44[0]=-1.640;
   R44[1]=-1.640;
   R44[2]=-1.698;
   R44[3]=-1.698;

   R46[0]=0.538;
   R46[1]=0.538;
   R46[2]=-0.037;
   R46[3]=-0.037;

   for(int i = 0; i < 4; i++){
      _y[i] = R33[i] * SF_y + R34[i] * SF_DYDZ + R36[i] * BMS_dpp;
      _y_err[i][0] = fabs(R33[i] * SF_y_err);
      _y_err[i][1] = fabs(R34[i] * SF_DYDZ_err);
      _y_err[i][2] = fabs(R36[i] * BMS_dpp_err);
      _y_err[i][3] = sqrt( pow( R33[i] * SF_y_err , 2)+ pow( R34[i] * SF_DYDZ_err ,2) + pow( R36[i] * BMS_dpp_err ,2)+ pow(5/sqrt(12.),2));
      if( i == 0 || i == 1 ) _dy[i] = R43[i] * SF_y +R44[i] * SF_DYDZ +R46[i] * BMS_dpp;
      if( i == 2 || i == 3 ) _dy[i] = R43[i] * SF_y +R44[i] * SF_DYDZ +R46[i] * BMS_dpp;

      /* it is not used now
      _x[i] = R11[i] * SF_x + R12[i] * SF_DXDZ;
      _x_err[i][0] = fabs(R11[i] * SF_x_err);
      _x_err[i][1] = fabs(R12[i] * SF_DXDZ_err);
      _x_err[i][2] = sqrt( pow( R11[i] * SF_x_err , 2)+ pow( R12[i] * SF_DXDZ_err ,2) + pow(5/sqrt(12.),2));
      if( i == 0 || i == 1 ) _dx[i]= R21[i] * SF_x + R22[i] * SF_DXDZ;
      if( i == 2 || i == 3 ) _dx[i]= R21[i] * SF_x + R22[i] * SF_DXDZ;
      */
   }

   //calculate for new planes
   _y[4] = _y[1] - ((_y[1]-_y[0])/(_planeReference[1]->getZcm()-_planeReference[0]->getZcm())) * (_planeReference[1]->getZcm()-_planeReference[4]->getZcm());
   //_x[4] = _x[1] - ((_x[1]-_x[0])/(_planeReference[1]->getZcm()-_planeReference[0]-getZcm())) * (_planeReference[1]->getZcm()-_planeReference[4]-getZcm());

   _y[5] = _y[3] - ((_y[3]-_y[2])/(_planeReference[3]->getZcm()-_planeReference[2]->getZcm())) * (_planeReference[3]->getZcm()-_planeReference[5]->getZcm());
   //_x[5] = _x[3] - ((_x[3]-_x[2])/(_planeReference[3]->getZcm()-_planeReference[2]-getZcm())) * (_planeReference[3]->getZcm()-_planeReference[5]-getZcm());
    
   return true;
}

//=== EXTRAPOLATE TO CLUSTER ==========================================================
double CsBeamBackPropagation :: extrapolate_to_cluster(int ipl, CsBMScluster* cluster){

   //check if empty
   if( cluster->isEmpty() ){

      CsErrLog::msg(elFatal,__FILE__,__LINE__,
      "CsBeamBackPropagation.cc: Error: extrapolate_to_cluster() cluster is empty");
   }

   //check if position of cluster equals position of station
   if( _planeReference[ipl]->getZcm() == cluster->get_z() ) return _y[ipl];

   //if different
   double a = 1.E6;

   if     ( ipl == 0 || ipl == 1 || ipl == 4 ){    //upstream part

      a =  (_y[1]-_y[0])/(_planeReference[1]->getZcm()-_planeReference[0]->getZcm()); 

   }
   else if( ipl == 2 || ipl == 3 || ipl == 5 ){    //downstream part

      a =  (_y[3]-_y[2])/(_planeReference[3]->getZcm()-_planeReference[2]->getZcm()); 

   }
   else{

      CsErrLog::msg(elFatal,__FILE__,__LINE__,
      "CsBeamBackPropagation.cc: Error: extrapolate_to_cluster() wrong plane number");

      return 1.E6;
   }
   
   //extrapolate
   return _y[ipl] + a * (cluster->get_z() - _planeReference[ipl]->getZcm());

}

//=== CALCULATE LH ====================================================================
double CsBeamBackPropagation :: calculate_chi2(CsBMStrack* BMStrack){

   int i;

   double chi2, y_corr;

   for(i = 0, chi2 = 0.; i < 6; i++){

      if( !((BMStrack->get_clusters()[i])->isEmpty()) ){

         y_corr = extrapolate_to_cluster(i, BMStrack->get_clusters()[i]); 

         chi2 += pow( ((BMStrack->get_clusters()[i])->get_y() - y_corr - _correctionFactor[i]) / _spaceResolutionBMS[i], 2 );
      }

   }

   return chi2;
}

//=== HISTOGRAMS ====================================================================
void CsBeamBackPropagation :: histogram_residuals(CsHelix* BTHelix, double mom, CsBMStrack* BMStrack){

   make_back_propagation(BTHelix, mom);

   int i;
   double y_corr;

   for(i = 0; i < 6; i++){

      if( !((BMStrack->get_clusters()[i])->isEmpty()) ){

         y_corr = extrapolate_to_cluster(i, BMStrack->get_clusters()[i]); 

         if( his_backprop_01[i] != NULL ) his_backprop_01[i]->Fill( (BMStrack->get_clusters()[i])->get_y() - y_corr, BMStrack->get_nPlanes());
         if( his_backprop_02[i] != NULL ) his_backprop_02[i]->Fill( (BMStrack->get_clusters()[i])->get_y() - y_corr, (BMStrack->get_clusters()[i])->get_y());
      }
   }
}

void CsBeamBackPropagation :: set_histograms(){

   static bool first = true;
   char title[100], name[100];

   if(first){

      for(int i = 0; i < 6; i++) his_backprop_01[i] = NULL;
      for(int i = 0; i < 6; i++) his_backprop_02[i] = NULL;

      std::string histograms_path =  "/BeamReconstruction/BackProp";
      CsHistograms::SetCurrentPath(histograms_path);

      for(int ipl = 0; ipl < 6; ipl++){

         //residuals for finall BMS tracks
         sprintf(title, "BM0%i, Y_{meas} - Y_{BackProp} finall tracks", ipl+1);
         sprintf(name,  "back_prop_01_BM0%i", ipl+1);
         if( _his_mode > 0 ) his_backprop_01[ipl] = new CsHist2D(name, title, 200, -50., 50., 7, 0, 7);

         //residuals for final BMS tracks vs y, all cases
         sprintf(title, "BM0%i, Y_{meas} - Y_{BackProp} finall tracks", ipl+1);
         sprintf(name,  "back_prop_02_BM0%i", ipl+1);
         if( _his_mode > 1 ) his_backprop_02[ipl] = new CsHist2D(name, title, 200, -50., 50., 300, -150., 150.);
      }

      first = false;
   }

}

//=== CLEAR ===========================================================================
void CsBeamBackPropagation :: clear(){

   int i;

   for( i = 0; i < 6; i++){
      _x[i] = 0.;
      _x_err[i][0] = 0.;   _x_err[i][1] = 0.;   _x_err[i][2] = 0.;
      _dx[i] = 0.;

      _y[i] = 0.;
      _y_err[i][0] = 0.;   _y_err[i][1] = 0.;   _y_err[i][2] = 0.;   _y_err[i][3] = 0.;
      _dy[i] = 0.;
   }

}

//=== INSTANCE ========================================================================
CsBeamBackPropagation* CsBeamBackPropagation :: instance_ = NULL;

CsBeamBackPropagation* CsBeamBackPropagation :: Instance() {
   if( instance_ == NULL ){
      instance_ = new CsBeamBackPropagation();
   }

   return( instance_ );
}

