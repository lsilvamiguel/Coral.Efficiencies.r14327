#include "coral_config.h"
#include "CsOpt.h"
#include "CsVertex.h"
#include "CsEvent.h"
#include "CsGeom.h"
#include "CsMCHit.h"
#include "CsMCDigit.h"
#include "CsMCVertex.h"
#include "CsDriftChamberDetector.h"
#include "CsDigit.h"
#include "CsCluster.h"
#include "CsMCTrkHit.h"
#include "CsTrack.h"
#include "CsECAL1.h"
#include "CsRwChargeRecons.h"
#include <list>
#include "CsHistograms.h"
#include "Reco/CalorimeterParticle.h"
#include "Reco/Cell.h"
#include <CLHEP/Matrix/Matrix.h>

#include "CsRichWallDetector.h"

using namespace std;

CsRwChargeRecons::CsRwChargeRecons(void) {
  
  string tag, key, str;
  int n;
  CsOpt* opt = CsOpt::Instance();
  
  if( CsInit::Instance()->IsAMonteCarloJob()) RW_MC=true;
  else  RW_MC=false;
  
  tag="RW_Charge_Recons";
  key="HistogramLevel";
  if( opt->getOpt( tag, key, RW_histoLevel) ) ; else RW_histoLevel=0;
  cout<<"CsRwChargeRecons: histoLevel "<<RW_histoLevel<<endl;

  readcal=false;
  key="RWoperdata";
  if( opt->getOpt( tag, key, readcal)&&readcal==true ) {

    vector<double> vec;

    key="RW_Summ_Range";
    if( opt->getOpt( tag, key, vec ) && (int)vec.size()==3) {  for(int i=0;i<3;i++) RW_Summ_Range[i]=vec[i]; }
    else {readcal=false;  cout<< " !!!! == Attention: No/bad RW_Charge_Recons  RW_Summ_Range  data"<<endl;}

    // ECal1 Clusters overlap parameters
    
    vec.clear();

    key="RW_Space_Overlap";
    if( opt->getOpt( tag, key, vec ) && (int)vec.size()==3) {  for(int i=0;i<3;i++) deltaSpace_overlap[i]=vec[i]; }
    else {readcal=false;  cout<< " !!!! == Attention: No/bad RW_Charge_Recons  RW_Space_Overlap  data"<<endl;}

    vec.clear();

    key="RW_Projection_Overlap";
    if( opt->getOpt( tag, key, vec ) && (int)vec.size()==3) {  for(int i=0;i<3;i++) deltaPR_overlap[i]=vec[i]; }
    else {readcal=false;  cout<< " !!!! == Attention: No/bad RW_Charge_Recons  RW_Projection_Overlap  data"<<endl;}

  }

  if(!readcal) {
    RW_Summ_Range[0]=50.; RW_Summ_Range[1]=60.; RW_Summ_Range[2]=70.;
    
    deltaPR_overlap[0]=75.; deltaPR_overlap[1]=75.; deltaPR_overlap[2]=75.;
    deltaSpace_overlap[0]=200.; deltaSpace_overlap[1]=200.; deltaSpace_overlap[2]=250.;
  }

  //    Get Detectors
  
//   idet1 =CsDet::FindDetector ( "EC01P1__" );
//   calorimeter= (CsECAL1 *)idet1;
 
  const std::string RWnames[nplanes]={"DR01X1__","DR01X2__","DR02X1__","DR02X2__",
				      "DR01Y1__","DR01Y2__","DR02Y1__","DR02Y2__"};

  for (int i=0;i<nplanes;i++) {
    RW_Plane_Z[i]= 0.;
    RW_Plane_X[i]= 0.;
    RW_Plane_Y[i]= 0.;

    CsDet *idet =CsDet::FindDetector ( RWnames[i] );
    if(idet==NULL) { 
      cout<< "RW plane "<<RWnames[i]<<", wasn't found"<<endl;
      continue;
    }
    RWdet[i]= (CsDetector*)idet; 
    RW_Plane_Z[i]=RWdet[i]->getZcm();
    RW_Plane_X[i]=RWdet[i]->getXcm();
    RW_Plane_Y[i]=RWdet[i]->getYcm();
  }
  
  char name[20],title[100];
  if(RW_histoLevel>0){ 
    string pathname =  "/RwChargeRecons";
    CsHistograms::SetCurrentPath(pathname);
    hi1D[1]=new CsHist1D("Nevents","events",100,0.,100.);
    hi2D[0]=new CsHist2D("XY_undef1","XY of ECAL1 cluster undef(far overlap)",400,-2000.,2000.,300,-1500.,1500.);
    hi2D[1]=new CsHist2D("XY_0","XY of ECAL1 cluster prob=0",400,-2000.,2000.,300,-1500.,1500.);
    hi2D[2]=new CsHist2D("XY_05","XY of ECAL1 cluster prob=0.5",400,-2000.,2000.,300,-1500.,1500.);
    hi2D[3]=new CsHist2D("XY_07","XY of ECAL1 cluster prob=0.7",400,-2000.,2000.,300,-1500.,1500.);
    hi2D[4]=new CsHist2D("XY_1","XY of ECAL1 cluster prob=1",400,-2000.,2000.,300,-1500.,1500.);
    hi2D[5]=new CsHist2D("XY_undef2","XY of ECAL1 cluster undef(close overlap)",400,-2000.,2000.,300,-1500.,1500.);
    pathname =  "/";
    CsHistograms::SetCurrentPath(pathname);
  }
}

void CsRwChargeRecons::RwChargeRecons(void){

  //==================================================
  if(RW_MC) RW_MC=false;         // make MC procedure the same as for real data
  //==================================================

  if(RW_histoLevel>2) {  cout << setprecision(5);
    cout<< " =============   NEW EVENT ==================="<<endl;  }

  CsEvent* event  = CsEvent::Instance();
  vector<CsParticle*> pat=event->getParticles();
  if(RW_histoLevel>0) hi1D[1]->Fill(1.);
  if(pat.size()>300) return;         // to reject bad events
  list<CsDigit*>digit=event->getDigits();
  list<CsCluster*>cluster=event->getClusters();

  double xtarg=0., ytarg=0., ztarg=0.;

  //  const std::vector<Reco::Cell> & cells = calorimeter->GetCells();
//
//    take ECAL1 digits
//

  std::vector<ec1gmm> ec1;   
  ec1.clear();
  
  double Eg=0.,Xg=0.,Yg=0.,Zg=0., Mom=0.;

  int ngamma=0, eczone=-1;
  if(pat.size()){
    int size=pat.size();
    for(int j=0;j<size;j++) {
      CsParticle*  ip = pat[j]; 
      vector<Reco::CalorimeterParticle*> cal = ip->getCalObjects();
      int ng= cal.size();

      int ng_track = 0;
 
      for (int i=0; i<ng;i++) {

	//	cal[i]->Set_RW_chg_Prob(-1);

	if(cal[i]->GetCalorimeterName()!="EC01P1__") continue;

	Eg = cal[i]->GetE();
	if(Eg<0.2) continue;

	Xg = cal[i]->GetX();
	Yg = cal[i]->GetY();
	Zg = cal[i]->GetZ();

	if(!gmchk(Xg, Yg, 38.3)) continue;
	if(fabs(Xg)<842.6) { 
	  if(fabs(Yg)<459.6){                                       //     Gams
	    eczone=0;
	  } else {                                                   // Mainz
	    eczone=1;
	  } 
	} else {                                                     //      Olga
	  eczone=2;
	}

	ng_track++;
	ec1gmm Phot;
	Phot.zone=eczone;
	Phot.X=Xg; Phot.Y=Yg; Phot.Zx=Zg; Phot.Zy=Zg; Phot.E=Eg;
	Phot.Mom=Mom;
	Phot.nhits=0;
	Phot.nfirst=9;
	Phot.npart=j;
	Phot.calobj=i;
	ec1.push_back(Phot);
      }
      ngamma+=ng_track;
    }         
  }

  if(ngamma>50) return;
  int ngsel=ec1.size();
  if(ngsel<1) return;

  int x_pl_hit = 0;
  int y_pl_hit = 0;

  for(int ig=0; ig<ngsel; ig++) {
    ec1gmm G=ec1[ig];          // get ECAl1 photon
    int    nz=G.zone;
    double X=G.X;
    double Y=G.Y;
    double Zx=G.Zx;
    double Zy=G.Zy;
    double E=G.E;
    double Mom=G.Mom;
    double Xg4=X, Zg4x=Zx, Yg4=Y, Zg4y=Zy;
    int corcase=3;

    int pat_id      = G.npart;
    int Calo_pat_id = G.calobj;

    vector<Reco::CalorimeterParticle*> tmp_cal_pat = (pat[pat_id])->getCalObjects();    

    x_pl_hit = 0;
    y_pl_hit = 0;
       
    double  Dmin, Dx,  Dy,  x,  y, zx, x1, y1, zy1;
    gm_gm_distance(ec1,  ig,  Dmin, Dx,  Dy,  x, y, zx,  x1, y1, zy1); 

    double x_hits[nplanes]={0.,0.,0.,0.,0.,0.,0.,0.};
    double y_hits[nplanes]={0.,0.,0.,0.,0.,0.,0.,0.};
    double nhits_tot=0.;
    
    if(Dmin>=deltaSpace_overlap[nz]) {

      if(Dx>=deltaPR_overlap[nz]&&Dy>=deltaPR_overlap[nz]){

	sum_digits(0,X,Y,Zx,RW_Summ_Range[nz],x_hits,x_pl_hit);
	sum_digits(4,X,Y,Zy,RW_Summ_Range[nz],y_hits,y_pl_hit);

	if(x_pl_hit > 2 && y_pl_hit > 2){
	  tmp_cal_pat[Calo_pat_id]->SetMiscInfo(Reco::CalorimeterParticle::RW_CHRG_PROB,1);

	  if(RW_histoLevel>0) hi2D[4]->Fill(X,Y);

	}else if(y_pl_hit > 2){
	  tmp_cal_pat[Calo_pat_id]->SetMiscInfo(Reco::CalorimeterParticle::RW_CHRG_PROB,0.7);

	  if(RW_histoLevel>0) hi2D[3]->Fill(X,Y);
	  
	}else if(x_pl_hit > 2){

	  tmp_cal_pat[Calo_pat_id]->SetMiscInfo(Reco::CalorimeterParticle::RW_CHRG_PROB,0.6);

	}else if(x_pl_hit == 0 && y_pl_hit == 0){
	  tmp_cal_pat[Calo_pat_id]->SetMiscInfo(Reco::CalorimeterParticle::RW_CHRG_PROB,0);

	  if(RW_histoLevel>0) hi2D[1]->Fill(X,Y);

	}else{
	  tmp_cal_pat[Calo_pat_id]->SetMiscInfo(Reco::CalorimeterParticle::RW_CHRG_PROB,0.5);

	  if(RW_histoLevel>0) hi2D[2]->Fill(X,Y);

	}

      }else if(Dx>=deltaPR_overlap[nz]&&Dy<deltaPR_overlap[nz]){

	sum_digits(0,X,Y,Zx,RW_Summ_Range[nz],x_hits,x_pl_hit);

	if(x_pl_hit > 2){
	  tmp_cal_pat[Calo_pat_id]->SetMiscInfo(Reco::CalorimeterParticle::RW_CHRG_PROB,1);

	  if(RW_histoLevel>0) hi2D[4]->Fill(X,Y);

	}else if(x_pl_hit == 0){
	  tmp_cal_pat[Calo_pat_id]->SetMiscInfo(Reco::CalorimeterParticle::RW_CHRG_PROB,0);

	  if(RW_histoLevel>0) hi2D[1]->Fill(X,Y);

	}else{
	  tmp_cal_pat[Calo_pat_id]->SetMiscInfo(Reco::CalorimeterParticle::RW_CHRG_PROB,0.5);

	  if(RW_histoLevel>0) hi2D[2]->Fill(X,Y);

	}

      }else if(Dx<deltaPR_overlap[nz]&&Dy>=deltaPR_overlap[nz]){

	sum_digits(4,X,Y,Zy,
		   RW_Summ_Range[nz],y_hits,y_pl_hit);
	
	if(y_pl_hit > 2){
	  tmp_cal_pat[Calo_pat_id]->SetMiscInfo(Reco::CalorimeterParticle::RW_CHRG_PROB,1);
	  
	  if(RW_histoLevel>0) hi2D[4]->Fill(X,Y);
	  
	}else if(y_pl_hit == 0){
	  tmp_cal_pat[Calo_pat_id]->SetMiscInfo(Reco::CalorimeterParticle::RW_CHRG_PROB,0);
	  
	  if(RW_histoLevel>0) hi2D[1]->Fill(X,Y);
	  
	}else{
	  tmp_cal_pat[Calo_pat_id]->SetMiscInfo(Reco::CalorimeterParticle::RW_CHRG_PROB,0.5);

	  if(RW_histoLevel>0) hi2D[2]->Fill(X,Y);
	  
	}

      }else{
	if(RW_histoLevel>0) hi2D[0]->Fill(X,Y);
      }

    }else{
      if(RW_histoLevel>0) hi2D[5]->Fill(X,Y);
    }
	 
  }
  
}

bool CsRwChargeRecons::gmchk(double Xg, double Yg, double delta) {
 static double Xhl=14*38.3, Yhl=8*38.3;
 static double Xmax=22*38.3+8*143., Ymax1=38*38.3, Ymax2=1420.;
 static double Xolga=22*38.3;
 if(fabs(Xg)<Xhl+delta&&fabs(Yg)<Yhl+delta) return false;
 if(fabs(Xg)>Xmax-delta) return false;
 if(fabs(Xg)>Xolga) {
            if(fabs(Yg)>Ymax2-delta) return false;}
 else  { 
            if(fabs(Yg)>Ymax1-delta) return false;}
 return true;
}

void CsRwChargeRecons::gm_gm_distance(std::vector<ec1gmm> & ec1, int ng , double & Dmin, double& Dx, double& Dy, 
				      double & x, double& y,double& zx, double & x1, double & y1,double& zy1)  {
  Dx=1000000.; Dy=1000000., Dmin=1000000.;
  int ngsel=(int)ec1.size();
  if(ngsel<2) return;
  ec1gmm G=ec1[ng];
  double X=G.X;
  double Y=G.Y;
  for(  int i=0; i<ngsel; i++ ) {
    if(i==ng) continue;
    ec1gmm G1=ec1[i];
    if(G1.E<0.8) continue;
    double X1=G1.X;
    double Y1=G1.Y;
    double dx=fabs(X-X1);
    double dy=fabs(Y-Y1);
    double D=sqrt(dx*dx+dy*dy);
    
    if(D<1.) continue;
    if(D<Dmin) { Dmin=D;  }
    if(dx<Dx){ Dx=dx; x=X1;y=Y1; zx=G1.Zx;}
    if(dy<Dy){ Dy=dy; x1=X1;y1=Y1;zy1=G1.Zy;}
  }
}

void CsRwChargeRecons::sum_digits(int ibg, double ph_X,double ph_Y, double ph_Z,
				  double Summ_Range,double hits[], int &pl_fired)  {

  // RWdetector[i]->Wire2Pos(current_wire) + RWdetector[i]->getWirD() + RWdetector[i]->getXcm() "X"
  // RWdetector[i]->Wire2Pos(current_wire) + RWdetector[i]->getWirD() + RWdetector[i]->getYcm() "Y"
  
  int RW_hits = 0;
  int wire = 0;
  double RW_XY_cm = 0;
  double RW_XY = 0;

  std::list<CsDigit*> RWdigits;
  std::list<CsDigit*>::iterator location;

  for(int i=ibg;i<ibg+4;i++) {
    
    //     list<CsCluster*>::iterator mcl;
    //     CsDetector* myrwdet = RWdet[i];
    //     std::list<CsCluster*> cluster=myrwdet->getMyClusters();
    
    hits[i]=0;
    double x= 0;
    double xyphot= 0;

    if(i<4) {x=ph_X; xyphot=ph_X+ph_X/ph_Z*(RW_Plane_Z[i]-ph_Z); RW_XY_cm=RWdet[i]->getXcm();}
    else    {x=ph_Y; xyphot=ph_Y+ph_Y/ph_Z*(RW_Plane_Z[i]-ph_Z); RW_XY_cm=RWdet[i]->getYcm();}
    

    RWdigits= RWdet[i]->getMyDigits(); 
    
    if( !RWdigits.empty() ){
      for( location = RWdigits.begin(); location != RWdigits.end(); ++location){

	wire = (*location)->getAddress();

	int nz= ((CsRichWallDetector*) RWdet[i])->DetectorSection(wire);
	if(!htchk(i, ph_X, ph_Y, nz,Summ_Range)) continue; 
		
	RW_XY = RWdet[i]->Wire2Pos(wire) + RWdet[i]->getWirD() + RW_XY_cm;

	if(fabs(RW_XY-xyphot)>Summ_Range) continue;

	hits[i-ibg]++;
	
      }
      
    }

  }

  for(int pl_chk=0; pl_chk<4 ; pl_chk++){
    
    if(hits[pl_chk] >0){
      pl_fired++;
    }
    
  }

  

}


bool CsRwChargeRecons::htchk(int ipl, double Xg, double Yg, int nz, double delta) {
 static double Xmin=-14*38.3, Xmax=14*38.3;
 static double Ymin=-8*38.3, Ymax=8*38.3;
 if(ipl<4) {
    if(Yg<(Ymin+delta)&&nz==3) return false;
    if(Yg>(Ymax-delta)&&nz==1) return false;
    if(Yg>(Ymin+delta)&&Yg<(Ymax-delta)) {
                    if(Xg<Xmin&&nz!=0) return false;
                    if(Xg>Xmax&&nz!=2) return false;
    }
    return true;
 }  else {
    if(Xg<(Xmin+delta)&&nz==3) return false;
    if(Xg>(Xmax-delta)&&nz==1) return false;
    if(Xg>(Xmin+delta)&&Xg<(Xmax-delta)) {
                    if(Yg<Ymin&&nz!=0) return false;
                    if(Yg>Ymax&&nz!=2) return false;
    }
    return true;
 }
}

//
CsRwChargeRecons* CsRwChargeRecons::instance_ = 0;
//
CsRwChargeRecons* CsRwChargeRecons::Instance() {
  if( instance_ == 0 ) instance_ = new CsRwChargeRecons();
  return( instance_ );
}
