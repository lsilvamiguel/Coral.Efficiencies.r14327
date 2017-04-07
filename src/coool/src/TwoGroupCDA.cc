#include "TVector3.h"

int TwoGroupCDA( std::map<std::string, Group*> fGroupMap  ){

  typedef std::map<std::string,Group*>::iterator GI;
  
  GroupSiliTrack* beforeTarget = 0;
  GroupSiliTrack* afterTarget  = 0;
  
  for(GI gi=fGroupMap.begin();
      gi!=fGroupMap.end(); gi++ ) {
    if (TString(gi->second->GetName()).Contains("SI01")&&
	TString(gi->second->GetName()).Contains("SI02"))
      beforeTarget=(GroupSiliTrack*)gi->second;
    
    if (TString(gi->second->GetName()).Contains("SI03")&&
	TString(gi->second->GetName()).Contains("SI04")&&
	TString(gi->second->GetName()).Contains("SI05"))
      afterTarget=(GroupSiliTrack*)gi->second;
  }

  if(!(beforeTarget&&afterTarget)) return -1;

  static TTree* tprim;
  static int prima_ntrack_before;
  static int prima_ntrack_after;
  static float prima_cda;
  static float prima_angle;
  static float prima_zpos;
  static float prima_xpos;
  static float prima_ypos;
  static float prima_cda_2;
  static float prima_angle_2;
  static float prima_zpos_2;
  static float prima_xpos_2;
  static float prima_ypos_2;
  static float prima_cda_3;
  static float prima_angle_3;
  static float prima_zpos_3;
  static float prima_xpos_3;
  static float prima_ypos_3;
  static bool first(true);
      
  prima_cda     = 100;
  prima_cda_2   = 100;
  prima_cda_3   = 100;
  prima_angle   = -0.1;
  prima_angle_2 = -0.1;
  prima_angle_3 = -0.1;
  prima_cda_2  = -10000;
  prima_zpos_2 = -10000;
  prima_xpos_2 = -10000;
  prima_ypos_2 = -10000;
  prima_cda_3  = -10000;
  prima_zpos_3 = -10000;
  prima_xpos_3 = -10000;
  prima_ypos_3 = -10000;

  if(first){ 
    tprim = new TTree("prima","Primakoff"); 
    tprim->Branch("ntrb", &prima_ntrack_before, "/I");
    tprim->Branch("ntra", &prima_ntrack_after, "/I");
    tprim->Branch("cda", &prima_cda, "/F");
    tprim->Branch("angle", &prima_angle, "/F");
    tprim->Branch("zpos", &prima_zpos, "/F");
    tprim->Branch("xpos", &prima_xpos, "/F");
    tprim->Branch("ypos", &prima_ypos, "/F");
    tprim->Branch("cda_2", &prima_cda_2, "/F");
    tprim->Branch("angle_2", &prima_angle_2, "/F");
    tprim->Branch("zpos_2", &prima_zpos_2, "/F");
    tprim->Branch("xpos_2", &prima_xpos_2, "/F");
    tprim->Branch("ypos_2", &prima_ypos_2, "/F");
    tprim->Branch("cda_3", &prima_cda_3, "/F");
    tprim->Branch("angle_3", &prima_angle_3, "/F");
    tprim->Branch("zpos_3", &prima_zpos_3, "/F");
    tprim->Branch("xpos_3", &prima_xpos_3, "/F");
    tprim->Branch("ypos_3", &prima_ypos_3, "/F");
    first=false;
  }
      
  if (beforeTarget && afterTarget) {

    prima_ntrack_before = beforeTarget->ntrack;
    prima_ntrack_after  = afterTarget->ntrack;
      
    if (beforeTarget->ntrack>=1 && afterTarget->ntrack>=1) {
      TVector3 P1(beforeTarget->crossing_point[0],
		  beforeTarget->crossing_point[1],
		  beforeTarget->crossing_point[2]);
      TVector3 P2(afterTarget->crossing_point[0],
		  afterTarget->crossing_point[1],
		  afterTarget->crossing_point[2]);
      TVector3 V1(beforeTarget->slopes[0],
		  beforeTarget->slopes[1], 1.0);
      TVector3 V2(afterTarget->slopes[0],
		  afterTarget->slopes[1], 1.0);

      TVector3 nn = (V1.Cross(V2)).Unit();
      prima_cda = TVector3((P1-P2) * nn).Mag();
      prima_angle =V1.Angle(V2);
	  
      TVector3 vv = P2 - P1 - prima_cda * nn;
      double tr_beta = V2.Angle(vv);
      double tr_c = vv.Mag() * sin(TMath::Pi()-tr_beta-prima_angle) /
	sin(prima_angle);
      double tr_b = vv.Mag() * sin(tr_beta) / sin(prima_angle);
      TVector3 CDA1 = P1 - tr_b * V1.Unit();	  
      TVector3 CDA2 = P2 - tr_c * V2.Unit();	  

      prima_zpos = CDA1[2];
      prima_xpos = CDA1[0] ;
      prima_ypos = CDA1[1];

      /* std::cout << "check CDA: " << (CDA1-CDA2).Mag() << std::endl;
	 std::cout << "CDA:   " << prima_cda << std::endl;
	 std::cout << "Angle: " << prima_angle << std::endl;
	 std::cout << " zpos: " << prima_zpos 
	 << " xpos: " << prima_xpos
	 << " ypos: " << prima_ypos << std::endl;
	 std::cout << "P1 " <<P1[0]<<" " <<P1[1]<<" "<<P1[2]<< std::endl;
	 std::cout << "V1 " <<V1[0]<<" " <<V1[1]<<" "<<V1[2]<< std::endl;

	 std::cout << "P2 " <<P2[0]<<" " <<P2[1]<<" "<<P2[2]<< std::endl;
	 std::cout << "V2 " <<V2[0]<<" " <<V2[1]<<" "<<V2[2]<< std::endl; 
      */

      if (afterTarget->ntrack>=2) {
	TVector3 P2_2(afterTarget->crossing_point2[0],
		      afterTarget->crossing_point2[1],
		      afterTarget->crossing_point2[2]);
	TVector3 V2_2(afterTarget->slopes2[0],
		      afterTarget->slopes2[1], 1.0);

	TVector3 nn_2 = (V1.Cross(V2_2)).Unit();
	prima_cda_2   = TVector3((P1-P2_2) * nn_2).Mag();
	prima_angle_2 =V1.Angle(V2_2);
	    
	TVector3 vv_2 = P2_2 - P1 - prima_cda_2 * nn_2;
	double tr_beta_2 = V2_2.Angle(vv);
	double tr_c_2 = vv_2.Mag() * 
	  sin(TMath::Pi()-tr_beta_2-prima_angle_2) / sin(prima_angle_2);
	double tr_b_2 = vv_2.Mag() * sin(tr_beta_2) / sin(prima_angle_2);
	TVector3 CDA1_2 = P1   - tr_b_2 * V1.Unit();	  
	TVector3 CDA2_2 = P2_2 - tr_c_2 * V2_2.Unit();	  
	    
	prima_zpos_2 = CDA1_2[2];
	prima_xpos_2 = CDA1_2[0] ;
	prima_ypos_2 = CDA1_2[1];

	TVector3 nn_3 = (V2.Cross(V2_2)).Unit();
	prima_cda_3   = TVector3((P2-P2_2) * nn_3).Mag();
	prima_angle_3 = V2.Angle(V2_2);
	    
	TVector3 vv_3 = P2_2 - P2 - prima_cda_3 * nn_3;
	double tr_beta_3 = V2_2.Angle(vv_3);
	double tr_c_3 = vv_3.Mag() * 
	  sin(TMath::Pi()-tr_beta_3-prima_angle_3) / sin(prima_angle_3);
	double tr_b_3 = vv_3.Mag() * sin(tr_beta_3) / sin(prima_angle_3);
	TVector3 CDA1_3 = P2   - tr_b_3 * V2.Unit();	  
	TVector3 CDA2_3 = P2_2 - tr_c_3 * V2_2.Unit();	  
	    
	prima_zpos_3 = CDA1_3[2];
	prima_xpos_3 = CDA1_3[0] ;
	prima_ypos_3 = CDA1_3[1];
	    
      }
	  
      tprim->Fill();
    }
    /*std::cout << "Event "<< fNevent << ": " 
      << beforeTarget->ntrack << ", " 
      << beforeTarget->crossing_point[0] << ", " 
      << beforeTarget->crossing_point[1] << ", " 
      << beforeTarget->crossing_point[2] << ", " 
      << beforeTarget->ntrack << ", " 
      << afterTarget->ntrack  << std::endl;
    */
	
  }  

  return 0;    
}
