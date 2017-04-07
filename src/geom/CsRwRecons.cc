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
#include "CsRwRecons.h"
#include <list>
#include "CsHistograms.h"
#include "Reco/CalorimeterParticle.h"
#include "Reco/Cell.h"
#include <CLHEP/Matrix/Matrix.h>

using namespace std;

//
CsRwRecons::CsRwRecons(void) {
        string tag, key, str;
        int n;
        CsOpt* opt = CsOpt::Instance();
//
        if( CsInit::Instance()->IsAMonteCarloJob()) RW_MC=true;
        else  RW_MC=false;
//
        tag="RW_Recons";
	key="HistogramLevel";
	if( opt->getOpt( tag, key, RW_histoLevel) ) ; else RW_histoLevel=0;
	cout<<"CsRwRecons: histoLevel "<<RW_histoLevel<<endl;
//
	key="EnergyCorrection";
	if( opt->getOpt( tag, key,EnergyCorrection ) ) ; else EnergyCorrection=0;
	cout<<"CsRwRecons: EnergyCorrection "<<EnergyCorrection<<endl;
//
	key="CoordCorrection";
	if( opt->getOpt( tag, key,CoordCorrection ) ) ; else CoordCorrection=0;
	cout<<"CsRwRecons: CoordCorrection "<<CoordCorrection<<endl;
//
	key="PhotonSelection";
	if( opt->getOpt( tag, key,PhotonSelection ) ) ; else PhotonSelection=0;
	cout<<"CsRwRecons: PhotonSelection "<<PhotonSelection<<endl;
//
        readcal=false;
	key="RWcorrectiondata";
       if( opt->getOpt( tag, key, readcal)&&readcal==true ) {
          vector<double> vec;
	  key="EnergyPoints";
	  if( opt->getOpt( tag, key, vec ) ) ; else  {readcal=false; 
                        cout<< " !!!! == Attention: No/bad RW_Recons  HitsSlope data"<<endl;}
          Ncorpoints=vec.size();
          for(int i=0;i<Ncorpoints;i++) Epoints[i]=vec[i];
	  key="EnergyCorr";
	  if( opt->getOpt( tag, key, vec ) ) { 
	       if(3*Ncorpoints!=(int)vec.size()) {   readcal=false; 
                            cout << " !!! ==  Attention: mismatch between EnergyPoints and EnergyNorm vector sizes"<<endl;}
               else {for(int i=0;i<3;i++) { for(int j=0;j<Ncorpoints;j++)  Ecorr_points[i][j]=vec[i*Ncorpoints+j];}}
          }
	  else  {readcal=false;  cout<< " !!!! == Attention: No/bad RW_Recons  EnergyNorm data"<<endl;}
//
          key="HitsSlopes";
	  if( opt->getOpt( tag, key, vec ) ) {
  	       if( 3*Ncorpoints!=(int) vec.size()) {   readcal=false; 
                            cout << " !!! ==  Attention: mismatch between EnergyPoints and HitsSlopes vector sizes"<<endl;}
               else {for(int i=0;i<3;i++) { for(int j=0;j<Ncorpoints;j++)  hits_slp[i][j]=vec[i*Ncorpoints+j];}}
          }
	  else {readcal=false;  cout<< " !!!! == Attention: No/bad RW_Recons  HitsSlope data"<<endl;}
//
          key="HitsDeltaEn";
	  if( opt->getOpt( tag, key, vec ) &&(int)vec.size()==3) {  for(int i=0;i<3;i++) deltaEn[i]=vec[i]; }
	  else {readcal=false;  cout<< " !!!! == Attention: No/bad RW_Recons  HitsDeltaEn data"<<endl;}
//          
          key="HitsDeltaCR";
	  if( opt->getOpt( tag, key, vec ) &&(int)vec.size()==3) {  for(int i=0;i<3;i++) deltaCRd[i]=vec[i]; }
	  else {readcal=false;  cout<< " !!!! == Attention: No/bad RW_Recons  HitsDeltaCR data"<<endl;}
//
	  if(!readcal) cout<< " !!! ==  Attention: Error during the RW correction data reading -> default values will be used"<< endl;
       } else cout<< "CsRwRecons Warning: No RW correction data in rw....opt file --> default values will be used "<<endl;
//
       if(!readcal) {
          deltaEn[0]=75.; deltaEn[1]=75.; deltaEn[2]=75.;
          deltaCRd[0]=50.; deltaCRd[1]=60.; deltaCRd[2]=70.;
       }
       deltaCR1=15.; 
       deltaPR_overlap[0]=75.; deltaPR_overlap[1]=75.; deltaPR_overlap[2]=75.;
       deltaSpace_overlap[0]=200.; deltaSpace_overlap[1]=200.; deltaSpace_overlap[2]=250.;
//
//     Default nhits_slope/ecorr values (can be used for 2006 with ECAl calibraion early 01.03.2007??)
//
       if(RW_MC) {
              hits_slope[0]=0.017;    hits_slope[1]=0.012;   hits_slope[2]=0.01;
              //   alfaDL[0]=0.013;  alfaDL[1]=0.013;  alfaDL[2]=0.012;
              Ecorr[0]=1.;      Ecorr[1]=1.;      Ecorr[2]=1.;
       }
      else {
             hits_slope[0]=0.064;   hits_slope[1]=0.0333;  hits_slope[2]=0.029;
            //   Ecorr[0]=0.8057;  Ecorr[1]=0.875;    Ecorr[2]=0.769;
             Ecorr[0]=0.738;  Ecorr[1]=0.748;    Ecorr[2]=0.8035;
      }
      if(RW_histoLevel>1) PRrecons_consts();
//
        int plmap1[nplanes]={0,1,2,3,4,5,6,7};
        int goodmap1[nplanes]={0,1,4,5,2,3,6,7};
        for(int i=0;i<nplanes;i++) {plmap[i]=plmap1[i]; goodmap[i]=goodmap1[i];}
//
//    Get Detectors
//
        idet1 =CsDet::FindDetector ( "EC01P1__" );
        calorimeter= (CsECAL1 *)idet1;
//
        const std::string RWnames[nplanes]={"DR01X1__","DR01X2__","DR02X1__","DR02X2__",
                                    "DR01Y1__","DR01Y2__","DR02Y1__","DR02Y2__"};
        for (int i=0;i<nplanes;i++) {
          z0cl[i]=0.;
          CsDet *idet =CsDet::FindDetector ( RWnames[i] );
          if(idet==NULL) { 
                      cout<< "RW plane "<<RWnames[i]<<"  didn't found"<<endl;
                      continue;}
          RWdet[i]= (CsDetector*)idet; 
          nwires[i]=RWdet[i]->getNWir();
          z0cl[i]=RWdet[i]->getZcm();
        }
//
  float emin=0.,emax=20.;
  char name[20],title[100];
if(RW_histoLevel>0){ 
  string pathname =  "/RwRecons";
  CsHistograms::SetCurrentPath(pathname);
  hi1D[1]=new CsHist1D("Nevents","events",100,0.,100.);
  hi1D[3]=new CsHist1D("P_e","Electorn momentum",200,0.,8.);
  hi1D[4]=new CsHist1D("Zel_last","Electron last mearured Z",400,0.,12000.);
  hi1D[5]=new CsHist1D("DXelectron","Electron DX",400,-200.,200.);
  hi1D[6]=new CsHist1D("DYelectron","Electron DY",400,-200.,200.);
//
  hi1D[9]=new CsHist1D("Eg_rawG","Gams, Raw photon energies",300,0.,30.);
  hi1D[10]=new CsHist1D("Eg_rawM","Mainz, Raw photon energies",200,0.,20.);
  hi1D[11]=new CsHist1D("Eg_rawO","Olga, Raw photon enegries",200,0.,20.);
//
  hi1D[12]=new CsHist1D("Ngtot","Total number of photons in ECAL1",15,0.,15.);
  hi1D[15]=new CsHist1D("Ngsel","Number of photons after selection",15,0.,15.);
  hi1D[16]=new CsHist1D("dif0","Hit-photon distance, plane 0",100,-300.,300.);
  hi1D[17]=new CsHist1D("dif1","Hit-photon distance, plane 1",100,-300.,300.);
  hi1D[18]=new CsHist1D("dif2","Hit-photon distance, plane 2",100,-300.,300.);
  hi1D[19]=new CsHist1D("dif3","Hit-photon distance, plane 3",100,-300.,300.);
  hi1D[20]=new CsHist1D("dif4","Hit-photon distance, plane 4",100,-300.,300.);
  hi1D[21]=new CsHist1D("dif5","Hit-photon distance, plane 5",100,-300.,300.);
  hi1D[22]=new CsHist1D("dif6","Hit-photon distance, plane 6",100,-300.,300.);
  hi1D[23]=new CsHist1D("dif7","Hit-photon distance, plane 7",100,-300.,300.);
  hi1D[24]=new CsHist1D("Nfirst","PS, first fired plane ",10,0.,10.);
  hi1D[25]=new CsHist1D("m2g_raw","2 gamma mass",200,0.,2.);
  hi1D[26]=new CsHist1D("m2g_cor","2 gamma corrected mass",200,0.,2.);
  hi1D[27]=new CsHist1D("m2g_corG","2 gamma corrected mass, Gams",200,0.,2.);
  hi1D[28]=new CsHist1D("m2g_corM","2 gamma corrected mass, Mainz",200,0.,2.);
  hi1D[29]=new CsHist1D("m2g_corO","2 gamma corrected mass, Olga",200,0.,2.);
  hi1D[30]=new CsHist1D("Eg_corG","  Corrected energy, Gams",100,0.,25.);
  hi1D[31]=new CsHist1D("Eg_corM","  Corrected energy, Mainz",100,0.,25.);
  hi1D[32]=new CsHist1D("Eg_corO","  Corrected energy, Olga",100,0.,25.);
  hi1D[36]=new CsHist1D("Nhits_tot","  Total number of hits in RW",100,0.,100.);
  hi1D[37]=new CsHist1D("Nhits_tot_first_0_7","Total number of RW hits,first0-7",100,0.,100.);
  hi1D[57]=new CsHist1D("m2g_rawG","2 gamma raw mass, Gams",200,0.,2.);
  hi1D[58]=new CsHist1D("m2g_rawM","2 gamma raw mass, Mainz",200,0.,2.);
  hi1D[59]=new CsHist1D("m2g_rawO","2 gamma raw mass, Olga",200,0.,2.);
  hi1D[60]=new CsHist1D("m2g_raw_nobkg","2 gamma raw mass, no bkg",200,0.,2.);
  hi1D[61]=new CsHist1D("m2g_cor_nobkg","2 gamma cor mass, no bkg",200,0.,2.);
//
  hi1D[62]=new CsHist1D("DE_elG","Electrons, relative energy difference ECAL-Pe",200,-5.,5.);
  hi1D[63]=new CsHist1D("DE_elM","Electrons, relative energy difference ECAL-Pe",200,-5.,5.);
  hi1D[64]=new CsHist1D("DE_elO","Electrons, relative energy difference ECAL-Pe",200,-5.,5.);
  hi1D[65]=new CsHist1D("DE_elG_cor","Electrons, relative Ecor difference ECAL-Pe",200,-5.,5.);
  hi1D[66]=new CsHist1D("DE_elM_cor","Electrons, relative Ecor difference ECAL-Pe",200,-5.,5.);
  hi1D[67]=new CsHist1D("DE_elO_cor","Electrons, relative Ecor difference ECAL-Pe",200,-5.,5.);
//
  hi1D[68]=new CsHist1D("Dmin_g_g","Minimal g-g distance",200,0.,600.);
  hi1D[69]=new CsHist1D("Dmin_g_g_X","Minimal g-g distance in X-projection",200,0.,600.);
  hi1D[70]=new CsHist1D("Dmin_g_g_Y","Minimal g-g distance in Y-projection",200,0.,600.);
//
  hi2D[0]=new CsHist2D("XY_noInt","ECAL1 XY, nfirst=8",50,-2000.,2000.,50,-1500.,1500.);
  hi2D[2]=new CsHist2D("Xelvs_Xg","Xel vs Xg",50,-2000.,2000.,50,-2000.,2000.);
  hi2D[3]=new CsHist2D("Yel_vs_Yg","Yel vs Yg",50,-1500.,1500.,50,-1500.,1500.);
  hi2D[4]=new CsHist2D("DXel_vs_XG","DXel vs Xg",50,-200.,200.,50,-2000.,2000.);
  hi2D[5]=new CsHist2D("DYel_vs_YG","DYel vs Yg",50,-200.,200.,50,-1500.,1500.);
  hi2D[6]=new CsHist2D("EelG_vs_Pe","Eel vs moment Gams",30,0.,15.,30,0.,15.);
  hi2D[7]=new CsHist2D("EelM_vs_Pe","Eel vs moment Mainz",30,0.,15.,30,0.,15.);
  hi2D[8]=new CsHist2D("EelO_vs_Pe","Eel vs moment Olga",30,0.,15.,30,0.,15.);
  hi2D[11]=new CsHist2D("Eel_vs_ndig","Eel vs ndig",50,0.,100.,50,emin,emax);
  hi2D[12]=new CsHist2D("XY_gamma","ECAL1 photon impackt points",100,-50*38.3,50*38.3,80,-40*38.3,40*38.3);
  hi2D[14]=new CsHist2D("m2g_raw_vs_E","raw 2 gamma mass",50,0.,25.,100,0.,2.);
  hi2D[15]=new CsHist2D("m2g_cor_vs_E_cor","Corrected 2 gamma mass",50,0.,25.,100,0.,2.);
  hi2D[16]=new CsHist2D("E_vs_nhitG","E_ECAL1 vs nhits,2.5GeV,GAMS",50,0.,100.,50,0.,10.);
  hi2D[17]=new CsHist2D("E_vs_nhitM","E_ECAL1 vs nhits,2.5GeV,Mainz",50,0.,100.,50,0.,10.);
  hi2D[18]=new CsHist2D("E_vs_nhitO","E_ECAL1 vs nhits,2.5GeV,Olga",50,0.,100.,50,0.,10.);
  hi2D[25]=new CsHist2D("Enew_vs_nhitG","E_ECAL1_new corrected vs nhits,2.5GeV,GAMS",50,0.,100.,50,0.,10.);
  hi2D[26]=new CsHist2D("Enew_vs_nhitM","E_ECAL1_new corrected vs nhits,2.5GeV,GAMS",50,0.,100.,50,0.,10.);
  hi2D[27]=new CsHist2D("Enew_vs_nhitO","E_ECAL1_new corrected vs nhits,2.5GeV,GAMS",50,0.,100.,50,0.,10.);
  hi2D[29]=new CsHist2D("XY_first0","ECAL1 XY, nfirst=0",50,-2000.,2000.,50,-1500.,1500.);
  hi2D[30]=new CsHist2D("XY_first_0_7","ECAL1 XY, hits1_7",50,-2000.,2000.,50,-1500.,1500.);
//
  hi2D[31]=new CsHist2D("EGams_raw_vs_Pe","EGams raw vs Pe, electrons",40,0.,8.,50,0.,10.);
  hi2D[32]=new CsHist2D("EMainz_raw_vs_Pe","EMainz raw vs Pe, electrons,",40,0.,8.,50,0.,10.);
  hi2D[33]=new CsHist2D("EOlga_raw_vs_Pe","EOlga raw vs Pe, electrons",40,0.,8.,50,0.,10.);
//
  hi2D[34]=new CsHist2D("DEelG_vs_Pe","DEel vs moment",40,0.,8.,50,-5.,5.);
  hi2D[35]=new CsHist2D("DEelM_vs_Pe","DEel vs moment",40,0.,8.,50,-5.,5.);
  hi2D[36]=new CsHist2D("DEelO_vs_Pe","DEel vs moment",40,0.,8.,50,-5.,5.);
  hi2D[37]=new CsHist2D("DEelG_cor_vs_Pe","DEel_cor vs moment",40,0.,8.,50,-5.,5.);
  hi2D[38]=new CsHist2D("DEelM_cor_vs_Pe","DEel_cor vs moment",40,0.,8.,50,-5.,5.);
  hi2D[39]=new CsHist2D("DEelO_cor_vs_Pe","DEel_cor vs moment",40,0.,8.,50,-5.,5.);
  pathname =  "/";
  CsHistograms::SetCurrentPath(pathname);
 }
}
//=======================================================================
CsRwRecons::~CsRwRecons() { Clear();}
//=======================================================================
void CsRwRecons::coordw(CsCluster* mycl, double &weight) {
       weight=0.;
       if(mycl->hasMirrorCluster()){ 
              CsCluster* mircl= mycl->getMirrorCluster();
              double prob=mycl->getLRProb();
              double mirprob=mircl->getLRProb();
              if(prob==mirprob) {weight=0.5; return;}
              if(prob>mirprob) weight=1.;
       } else{
               weight=1.;
       }
//  next line just for test
//       weight=mycl->getLRProb();
}
//=======================================================================
void CsRwRecons::RwRecons() {
//==================================================
  if(RW_MC) RW_MC=false;         // make MC procedure the same as for real data
//==================================================
  Clear();
  if(RW_histoLevel>2) {  cout << setprecision(5);
                cout<< " =============   NEW EVENT ==================="<<endl;  }
//
  CsEvent*  event  = CsEvent::Instance();
  vector<CsParticle*> pat=event->getParticles();
  if(RW_histoLevel>0) hi1D[1]->Fill(1.);
  if(pat.size()>300) return;         // to reject bad events from ECAl1 recontruction
  if(RW_histoLevel>0) hi1D[1]->Fill(2.);
  list<CsDigit*>digit=event->getDigits();
  list<CsCluster*>cluster=event->getClusters();
    double xtarg=0., ytarg=0., ztarg=0.;
    const std::vector<Reco::Cell> & cells = calorimeter->GetCells();
//
//    take ECAL1 digits
//
  list<CsDigit*>::iterator mdg;
  list<CsDigit*> ec1dig;  ec1dig.clear();
  for( mdg = digit.begin(); mdg != digit.end(); mdg++ ) {
       CsDigit* mydg = *mdg;
       CsDet* dgdet = mydg->getDet();
       if(dgdet == idet1) {                           // ECAL1 digits
                ec1dig.push_back(mydg);
//                int ncell=mydg->getAddress();
//		int digsize=mydg->getDataSize();
//		double amp=0.;
//		if(digsize==1) amp = mydg->getDatum();
//		else cout<<"  Bad size of the ECAL1 digit="<<digsize<<endl;
//                cout<< "calorim.digit, cell="<<ncell<<" ampl="<< amp<<endl;
        }
  }
//
//  take  Ecal1 data
//
// const bool selel = true; //false
 const bool selel = false; //false
 double Eg=0.,Xg=0.,Yg=0.,Zg=0., Mom=0.;
 double Xel_ecal=0, Yel_ecal=0;
 std::vector<ec1rec>  ec1;   ec1.clear();
 std::vector<ec1rec>  ec1_test;   ec1_test.clear();
 CsHelix hel,helec1;
 int ngamma=0, eczone=-1;
 if(pat.size()){
     int size=pat.size();
      for(int j=0;j<size;j++) {
        CsParticle*  ip = pat[j]; 
         vector<Reco::CalorimeterParticle*> cal = ip->getCalObjects(); 
         const CsTrack* trackpnt = ip->getTrack();
         Mom=0.;
         if(selel) {                                                      // for test purposes
 	      if(trackpnt ==NULL) continue; 
              if(!trackpnt->hasRich1Probs()) continue;
              int ptID=LikePid(const_cast<CsTrack*>(trackpnt));
              if(ptID!=3) continue;
              if(RW_histoLevel>0) hi1D[1]->Fill(3.);
              const std::vector<CsHelix> h = trackpnt->getHelices();
              if(h.size()<1) continue;
              int nlst=h.size()-1;
              double Zlast=h[nlst].getZ(); 
              if(RW_histoLevel>0)  hi1D[4]->Fill(Zlast);
              if(Zlast>11000.) continue;
              hel=h[nlst];
              Mom=abs(h[nlst].getCop());
              if(Mom<=0.) continue;
              Mom=1./Mom;
              if(RW_histoLevel>0)  hi1D[3]->Fill(Mom);
              if(Mom<0.8) continue;
              if(RW_histoLevel>0) hi1D[1]->Fill(4.);
     }   else {    
 	     if(trackpnt !=NULL) continue; 
         }     
	 int ng= cal.size();
//	 if(ng!=1) continue;
//         cout<< "ng="<<ng<<endl;
         int ng_track=0;
	 for (int i=0; i<ng;i++) {
                   if(cal[i]->GetCalorimeterName()!="EC01P1__") continue;
		   Eg = cal[i]->GetE();
                   if(Eg<0.2) continue;
		   Xg = cal[i]->GetX();
		   Yg = cal[i]->GetY();
		   Zg = cal[i]->GetZ();
//	   	   cout << " calorim.track:x="<< Xg<<" y="<<Yg<<" z="<<Zg<<" E="<<Eg<<endl;
                   if(RW_histoLevel>0) hi1D[1]->Fill(5.);
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
                   if(RW_histoLevel>0) {
                                hi1D[1]->Fill(16.);
                                hi1D[9+eczone]->Fill(Eg);
                                hi2D[12]->Fill(Xg,Yg); 
                               if(selel) { 
                                   hi2D[31+eczone]->Fill(Mom,Eg);
                                   hi1D[62+eczone]->Fill(Eg/Mom-1.);
                                   hi2D[34+eczone]->Fill(Mom,Eg/Mom-1.);
                                   hel.Extrapolate(Zg, helec1, false);
                                   Xel_ecal=helec1.getX();
                                   Yel_ecal=helec1.getY();
                               }
                  }
                   ng_track++;
                   ec1rec   Phot;
                   Phot.zone=eczone;
                   Phot.X=Xg; Phot.Y=Yg; Phot.Zx=Zg; Phot.Zy=Zg; Phot.E=Eg;
                   Phot.Mom=Mom;
                   Phot.nhits=0;
                   Phot.nfirst=9;
                   if(selel) Phot.hel=hel;
                   Phot.npart=j;
                   Phot.calobj=i;
                   ec1.push_back(Phot);
                   ec1_test.push_back(Phot);
                   if(RW_histoLevel>2) std::cout<<"raw gamma, nz="<<eczone<<" E="<<Eg<<"   Xg="<<Xg<<"  Yg="<<Yg<<"  Zg="<<Zg<<std::endl;
 	 }
         ngamma+=ng_track;
      }         
  }
  if(RW_histoLevel>0) hi1D[12]->Fill(ngamma);
  if(ngamma>50) return;
  int ngraw=ec1.size();
  if(RW_histoLevel>2) std::cout<<"Number of selecter raw gammas  "<<ngraw<<std::endl;
  if(ngraw<1) return;
  if(RW_histoLevel>0) hi1D[1]->Fill(6.);
//
//       kill false gammas
//
   if(PhotonSelection) {
//  bool yes;
//  kill_all_ghnosts(ec1, yes);
// cout<< " after ngamma="<<ngamma<<"  "<<ec1.size()<<"  "<<yes<<endl;
   }
  int ngsel=ec1.size();
  if(RW_histoLevel>0)  hi1D[15]->Fill(ngsel); 
  if(ngsel<1) return;
  if(RW_histoLevel>0) hi1D[1]->Fill(7.);
//
//    Photon energy/coordinate correction using RW
//
 std::vector<ec1rec>  ec1_rw;   ec1_rw.clear();
 int nplx=4, nply=4;
 for(  int ig=0; ig<ngsel; ig++) {
         ec1rec G=ec1[ig];          // get ECAl1 photon
         int nz=G.zone;
         double X=G.X;
         double Y=G.Y;
         double Zx=G.Zx;
         double Zy=G.Zy;
         double E=G.E;
         double Mom=G.Mom;
         double Xg4=X, Zg4x=Zx, Yg4=Y, Zg4y=Zy;
         int corcase=3;
         if(RW_histoLevel>2) cout<<"Selected gamma="<<ig<<"  "<<"X="<<X<<"  Y="<<Y<<" E="<<E<<"  type="<< nz<<endl;
//
         double  Dmin, Dx,  Dy,  x,  y, zx, x1, y1, zy1;
         gm_gm_distance(ec1,  ig,  Dmin, Dx,  Dy,  x, y, zx,  x1, y1, zy1); 
         if(RW_histoLevel>0&&ngsel>1) {
               hi1D[68]->Fill(Dmin);
               hi1D[69]->Fill(Dx);
               hi1D[70]->Fill(Dy);
               if(RW_histoLevel>2)   cout<<" Dmin="<<Dmin<<" Dx="<<Dx<<" Dy="<<Dy<<endl;   
         }
//
         double xmean[nplanes]={0.,0.,0.,0.,0.,0.,0.,0.};
         double wpl[nplanes]={0.,0.,0.,0.,0.,0.,0.,0.};
         double xmean4[nplanes]={0.,0.,0.,0.,0.,0.,0.,0.};
         double w4[nplanes]={0.,0.,0.,0.,0.,0.,0.,0.};
         double wtotx,wtoty,wxmn,wymn;
         double nhits_tot=0.;
         double Eec1new=0.;
         if(RW_histoLevel>0) hi1D[1]->Fill(9.);
         if(Dmin>=deltaSpace_overlap[nz]) {
                 if(RW_histoLevel>0) hi1D[1]->Fill(10.);
                 if(Dx>=deltaPR_overlap[nz]&&Dy>=deltaPR_overlap[nz]){
                      corcase=0;  
                      if(RW_histoLevel>0) hi1D[1]->Fill(11.);
                      sum_digits(0,X,Y,Zx,deltaEn[nz],deltaCRd[nz],wtotx,wpl,xmean);
                      sum_digits(4,X,Y,Zx,deltaEn[nz],deltaCRd[nz],wtoty,wpl,xmean);
                      if(CoordCorrection) {
                        sum_CRdigits(0,deltaCR1,X,Y,xmean,w4,xmean4);
                        sum_CRdigits(4,deltaCR1,X,Y,xmean,w4,xmean4);
                        CrCorrection(nplx, nz, X, Zx, w4,xmean4, z0cl,Xg4,Zg4x,wxmn);
                        CrCorrection(nply, nz, Y, Zy, &w4[4],&xmean4[4],&z0cl[4],Yg4,Zg4y,wymn);
                      }
                      if(RW_histoLevel>2) {
//   cout<<"Case 0"<<endl;
//   cout<<xmean[0]<<"  "<<xmean[1]<<"  "<<xmean[2]<<"  "<<xmean[3]<<"  "<<xmean[4]<<"  "<<xmean[5]<<"  "<<xmean[6]<<"  "<<xmean[7]<<endl;
//   cout<<xmean4[0]<<"  "<<xmean4[1]<<"  "<<xmean4[2]<<"  "<<xmean4[3]<<"  "<<xmean4[4]<<"  "<<xmean4[5]<<"  "<<xmean4[6]<<"  "<<xmean4[7]<<std::endl;
//   cout<<wpl[0]<<"  "<<wpl[1]<<"  "<<wpl[2]<<"  "<<wpl[3]<<"  "<<wpl[4]<<"  "<<wpl[5]<<"  "<<wpl[6]<<"  "<<wpl[7]<<endl;
//   std::cout<<w4[0]<<"  "<<w4[1]<<"  "<<w4[2]<<"  "<<w4[3]<<"  "<<w4[4]<<"  "<<w4[5]<<"  "<<w4[6]<<"  "<<w4[7]<<std::endl;
    cout<<"Case0, xcor="<<Xg4<<"  Zx="<<Zg4x<<" ycor="<<Yg4<<" Zy"<<Zg4y<<" wtotx="<<wtotx<<" wtoty"<<wtoty<<endl; 
                      }
                      nhits_tot=wtotx+wtoty;
                      if(RW_histoLevel>0) {
                        int nfirst=nplanes;
                        for(int i=0; i<nplanes; i++) {
                          int ii=goodmap[i];
                          if(wpl[ii]!=0) { nfirst=i; break;}
                        }
                        G.nfirst=nfirst;
                        ec1_test[ig].nfirst=nfirst;
                        hi1D[24]->Fill(nfirst);
                        if(nfirst==8) hi2D[0]->Fill(X,Y);
                        if(nfirst==0) hi2D[29]->Fill(X,Y);
                        if(nfirst<8) hi2D[30]->Fill(X,Y);
                        hi1D[36]->Fill(nhits_tot);
                        if(nfirst<8) hi1D[37]->Fill(nhits_tot);
                      }
                } 
                 else if(Dx>=deltaPR_overlap[nz]&&Dy<deltaPR_overlap[nz]){;
                      corcase=1;  
                      if(RW_histoLevel>0) hi1D[1]->Fill(12.);
//                       cout<<"Case 1"<<endl;
                       double  wtotx1=0.;
                       double wplwk[nplanes], xmeanwk[nplanes];
                       sum_digits(0, X, Y, Zx, deltaEn[nz], deltaCRd[nz],  wtotx,  wpl,  xmean ); 
                       sum_digits(0, x, y, zx, deltaEn[nz], deltaCRd[nz],  wtotx1,  wplwk,  xmeanwk ); 
                       if(CoordCorrection) {
                         sum_CRdigits(0, deltaCR1, X, Y, xmean,  w4,  xmean4 );
                         CrCorrection(nplx, nz, X, Zx, w4,xmean4, z0cl, Xg4, Zg4x, wxmn);
                       }
                       if((wtotx+wtotx1)<1.) wtoty=0.;
                       else  { 
                            double  deltaEy=deltaEn[nz]+Dy;
                            sum_digits(4, X, Y, Zy, deltaEy, deltaCRd[nz],  wtoty,  wpl, xmean );
                            wtoty=wtoty*wtotx/(wtotx+wtotx1);
                       }
                       nhits_tot=wtotx+wtoty; 
//                       cout<<E<<"  "<<nhits<<"  "<<nz<<" Case 2"  <<endl; 
                       wymn=1000000;
                }
                 else if(Dx<deltaPR_overlap[nz]&&Dy>=deltaPR_overlap[nz]){;
                      corcase=2;  
                      if(RW_histoLevel>0) hi1D[1]->Fill(13.);
//                        cout<<"Case 2"<<endl;
                       double  wtoty1=0.;
                       double wplwk[nplanes], xmeanwk[nplanes];
                       sum_digits( 4, x1, y1,zy1, deltaEn[nz], deltaCRd[nz],  wtoty1,  wplwk, xmeanwk );
                       sum_digits( 4, X, Y, Zy, deltaEn[nz], deltaCRd[nz],  wtoty,  wpl, xmean );
                       if(CoordCorrection) {
                         sum_CRdigits(4,X,Y, deltaCR1,  xmean,  w4,  xmean4 );
                         CrCorrection(nply, nz, Y, Zy, &w4[4],&xmean4[4],&z0cl[4],Yg4,Zg4y,wymn);
                       }
                       double deltaEx=deltaEn[nz]+Dx;
                       if(wtoty+wtoty1<1.) wtotx=0.;
                       else  { 
                            sum_digits(0, X, Y, Zx, deltaEx, deltaCRd[nz],  wtotx,  wpl,  xmean ); 
                            wtotx=wtotx*wtoty/(wtoty+wtoty1);
                       }
                       nhits_tot=wtotx+wtoty; 
                       wxmn=1000000;
               }   else {  // cout<<"Case 3"<<endl;
                         if(RW_histoLevel>0) hi1D[1]->Fill(14.);
                         nhits_tot=0; Xg4=X; Yg4=Y;
                         corcase=3;  
              }
//       
              Eec1new=E;
//              if(nz==0) {                 
//                   double  xcut=100., ycut=100.;
//                   ec1energy_cor(xcut,  ycut,  ec1dig,  X, Y, cells, Eec1new);
//              }  
//
       } else {  //   cout<<"Case 5"<<endl;
                   if(RW_histoLevel>0) hi1D[1]->Fill(15.);
                   nhits_tot=0.; Eec1new=E;
                   Xg4=X; Yg4=Y;
       }
//       
//    RW+ECAL1  photon energy 
//
    double a,b;
    corcoeff(E, nhits_tot, nz,  a, b);
    double  Enew=a*(E+b*(nhits_tot-8));
    if(RW_histoLevel>0) {    
      hi1D[30+nz]->Fill(Enew);
      if(fabs(Enew-2.5)<0.5){
          hi2D[16+nz]->Fill(nhits_tot,E); 
          hi2D[25+nz]->Fill(nhits_tot,Enew); }
       if(selel) {
          hi1D[65+nz]->Fill(Enew/Mom-1);
          hi2D[37+nz]->Fill(Mom,Enew/Mom-1.);
       }
    }
//   
//      save corrected photon
//
         if(Zg4x!=Zg4y){ if(Zg4x>Zg4y) { Yg4=Yg4*Zg4x/Zg4y; Zg4y=Zg4x;}
                                  else { Xg4=Xg4*Zg4y/Zg4x; Zg4x=Zg4y;}
         }
         ec1rec NewPhot=G;
         NewPhot.X=Xg4;
         NewPhot.Y=Yg4;
         NewPhot.Zx=Zg4x;
         NewPhot.Zy=Zg4y;
         NewPhot.E=Enew;
         ec1_rw.push_back(NewPhot);
         if(RW_histoLevel>2) cout<<"Correctted gamma="<<ig<<"  "<<"X="<<Xg4<<"  Y="<<Yg4
                                               <<" E="<<Enew<<"  corcase="<< corcase<<endl;
}
//
//   make new CalorimeterParticles
//
int nrecg=ec1_rw.size();
if(nrecg>0) {
     for(int i=0;i<nrecg;i++){
        int ipt=ec1_rw[i].npart;
        int icalo=ec1_rw[i].calobj;
        CsParticle*  ip = pat[ipt]; 
        vector<Reco::CalorimeterParticle*> cal = ip->getCalObjects(); 
        Reco::CalorimeterParticle* cp = new Reco::CalorimeterParticle(*cal[icalo]);
	CsEvent::Instance()->AddCalObject(cp);                    // save pointer into "standard" CORAL place 
        cp->SetFictionalCalorimeterName("RW_ECAL1");
//        std::cout << std::setprecision(5);
//        std::cout<<ec1_rw[i].E<<" "<<ec1_rw[i].X<<" "<<ec1_rw[i].Y<<"  "<<ec1_rw[i].Zx<<std::endl;
        if(EnergyCorrection){
          cp->SetE(ec1_rw[i].E,cp->GetEerr());
        }
        if(CoordCorrection){
          cp->SetX(ec1_rw[i].X,cp->GetXerr());
          cp->SetY(ec1_rw[i].Y,cp->GetYerr());
          cp->SetZ(ec1_rw[i].Zx,cp->GetZerr());
        }
        ec1out *ecout= new ec1out;
        ecout->npart=ipt; 
        ecout->calobj=cp; 
        dataout.push_back(ecout); 
     } 
}
if(RW_histoLevel>0){
    m2g(ec1_test,25,14,57,60);
    m2g(ec1_rw,26,15,27,61);
}
return;
}
//
//==================================================================================
void CsRwRecons::CrCorrection(int npl, int nz, double Xg, double Zg, double w4[],double xmean4[], double z0cl[],
                                                 double & Xg4,double & Zg4x, double &wxmn){
double Zcorrection[3]={0.,0.,0.};
double lmax[3]={5.,10.,18.};
//double Zcorrection[3]={85.,20.,90.};
  wxmn=10000.;   
  for(int i=0; i<npl; i++) {
     if(w4[i]>0.&&w4[i]<wxmn) wxmn=w4[i];
  }
  if(wxmn==1) { 
       int nfr=0; 
       double xmn4=0., zmn4=0.;
       for(int i=0; i<npl; i++) {
          if(w4[i]==1) { 
               nfr++; 
               xmn4+=xmean4[i];
               zmn4+=z0cl[i];
          }
       }
       Xg4=xmn4/nfr;
       Zg4x=zmn4/nfr;
  }
  if(wxmn>1&&wxmn<1000.) { 
       int nfr=0; 
       double xmn4=0., zmn4=0.;
       for(int i=0; i<npl; i++) {
          if(w4[i]>0.) { 
               nfr++; 
               xmn4+=xmean4[i];
               zmn4+=z0cl[i];
          }
       }
       Xg4=xmn4/nfr;
       Zg4x=zmn4/nfr;
  }
  if(wxmn>1000.) { Xg4=Xg; 
              Zg4x=Zg+Zcorrection[nz];
              return; 
  }
//  double xg=Xg4*Zg/Zg4x;
//  if(fabs(xg-Xmax)>lmax[lb]) {Xg4=Xg; Zg4x=Zg+Zcorrection[lb];}
}
//================================================================================
void  CsRwRecons::gm_gm_distance(std::vector<ec1rec> & ec1, int ng , double & Dmin, double& Dx, double& Dy, 
                 double & x, double& y,double& zx, double & x1, double & y1,double& zy1)  {
  Dx=1000000.; Dy=1000000, Dmin=1000000.;
  int ngsel=(int)ec1.size();
  if(ngsel<2) return;
         ec1rec G=ec1[ng];
         double X=G.X;
         double Y=G.Y;
         for(  int i=0; i<ngsel; i++ ) {
             if(i==ng) continue;
             ec1rec G1=ec1[i];
             if(G1.E<0.8) continue;
             double X1=G1.X;
             double Y1=G1.Y;
             double dx=fabs(X-X1);
             double dy=fabs(Y-Y1);
             double D=sqrt(dx*dx+dy*dy);
//             if(D<1.) cout<< ng<<"  "<<i<<"================================"<<endl;
             if(D<1.) continue;
             if(D<Dmin) { Dmin=D;  }
             if(dx<Dx){ Dx=dx; x=X1;y=Y1; zx=G1.Zx;}
             if(dy<Dy){ Dy=dy; x1=X1;y1=Y1;zy1=G1.Zy;}
        }
}
//====================================================================
void  CsRwRecons::kill_ghnost(std::list<ec1rec>&  ec1, bool& yes, int mod)  {
  double distcut=450., rtcut=0.3;
  double  Dmin=1000000., ecutmx=0.;
  yes=false;
  if(ec1.size()<2) return;
  //  list<ec1rec>::iterator phot, phot1, phot2, photsv=NULL, photsv1=NULL;
  list<ec1rec>::iterator phot, phot1, phot2, photsv=ec1.end(), photsv1=ec1.end();
  for(  phot = ec1.begin();  phot != ec1.end();  phot++ ) {
         ec1rec G=*phot;
         double X=G.X;
         double Y=G.Y;
         phot2=phot;
         phot2++;
        for(  phot1 = phot2;  phot1!= ec1.end();  phot1++ ) {
             ec1rec G1=*phot1;
             double X1=G1.X;
             double Y1=G1.Y;
             double D=sqrt((X-X1)*(X-X1)+(Y-Y1)*(Y-Y1));
             if(mod==0) {
	        if(D<Dmin) {  Dmin=D;  photsv=phot; photsv1=phot1; }
             } else {
	        if(D>distcut) continue;
                double ecut=rtcut*(1.-D/distcut);
                double EE=0.;
                if(G.E>G1.E) EE=G1.E/G.E; else EE=G.E/G1.E;
                if(EE>ecut) continue;
	        if(ecut>ecutmx){ Dmin=D; photsv=phot; photsv1=phot1; ecutmx=ecut;}
	    }
       }
  }
         if(photsv==ec1.end()||photsv1==ec1.end()) return;
         ec1rec G=*photsv;
         ec1rec G1=*photsv1;
         double EE=0.;
         if(G.E>G1.E) EE=G1.E/G.E; else EE=G.E/G1.E;
         int type=0;
         if(G.zone==1&&G1.zone==1) type=1;
         if(G.zone==2&&G1.zone==2) type=2;
         if(Dmin>distcut) return;
         if(mod==1) {
            double ecut=rtcut*(1.-Dmin/distcut);
            if(EE>ecut) return;
         }
         yes=true;
         double E=G.E+G1.E;
         G.X=(G.X*G.E+G1.X*G1.E)/E;
         G.Y=(G.Y*G.E+G1.Y*G1.E)/E;
         G.E=E;
         ec1.erase(photsv1);
	 double alfa=0.,beta=0., xtarg=0.,ytarg=0.,ztarg=0.;
         find_best_MCvertex(G.X, G.Y, G.Zx, G.Zy, xtarg, ytarg, ztarg, alfa, beta);
         double X0=xtarg+alfa*G.Zx;
         double Y0=ytarg+beta*G.Zy;
         if(fabs(G.X)>=842.6) {G.zone=2;return; }                                      //     OLGA
         if(fabs(G.Y)<459.6)  {G.zone=0;return; }                                      //     Gams
         G.zone=1;                                                                                     //     Mainz
}
//======================================================================
void  CsRwRecons::kill_all_ghnosts(std::list<ec1rec>&  ec1, bool& yes1)  {
     bool yes=true;
     yes1=false;
     while(yes)   {kill_ghnost(ec1,yes, 1);  if(yes) yes1=true;}
}
//======================================================================
void  CsRwRecons::sum_digits(int ibg, double X,double Y, double Z,double deltaE,
         double deltaCR,double &wtot,double wpl[],double xmean[] )  {
    wtot=0.; 
    for(int i=ibg;i<ibg+4;i++) {
         list<CsCluster*>::iterator mcl;
         CsDetector* myrwdet = RWdet[i];
         std::list<CsCluster*> cluster=myrwdet->getMyClusters();  
         wpl[i]=0.; xmean[i]=0;
         double x,xyphot;
         if(i<4) {x=X; xyphot=X+X/Z*(z0cl[i]-Z);}
         else    {x=Y; xyphot=Y+Y/Z*(z0cl[i]-Z);}
         for( mcl = cluster.begin(); mcl != cluster.end(); mcl++ ) {
              CsCluster* mycl = *mcl;
              int nz=wrchk(i, mycl);
              if(!htchk(i, X, Y, nz,deltaE)) continue; 
              double xcl=mycl->getU();
              double weight=0.;
              coordw(mycl,weight);
//              cout<<"new pl="<<i<<"   "<<xcl<<"  "<<z0cl[i]<<"   "<<xyphot<<"  "
//                  <<weight<<"  "<<nz<<"  "<<deltaCR<<"  "<<fabs(xcl-xyphot)<<endl;
              if(RW_histoLevel>0) hi1D[16+i]->Fill(xcl-xyphot);
              if(fabs(xcl-xyphot)<deltaE) wtot+=weight; 
              if(fabs(xcl-xyphot)>deltaCR) continue;
              if(!htchk(i, X, Y, nz,deltaCR)) continue; 
              wpl[i]+=weight;
              xmean[i]+=xcl*weight; 
//             cout << " ID= "<<i<<" X="<<xcl<<" Z="<<z0cl[i]<<endl;
          }
          if(wpl[i]>0.55) xmean[i]=xmean[i]/wpl[i]; 
          else{ wpl[i]=0.;   xmean[i]=x;}
    }
 }
//=======================================================================
void  CsRwRecons::sum_CRdigits( int ibg,double delta, double Xg, double Yg, 
                               double xmean[], double w4[], double xmean4[] )  {
  for(int i=ibg;i<ibg+4;i++) {
       list<CsCluster*>::iterator mcl;
       CsDetector* myrwdet = RWdet[i];
       std::list<CsCluster*> cluster=myrwdet->getMyClusters();  
       w4[i]=0.; xmean4[i]=0.;
       for( mcl = cluster.begin(); mcl != cluster.end(); mcl++ ) {
            CsCluster* mycl = *mcl;
            int nz=wrchk(i, mycl);
            if(!htchk(i, Xg, Yg, nz,delta)) continue; 
            double xcl=mycl->getU();
            double weight=0.;
            coordw(mycl,weight);
            if(fabs(xcl-xmean[i])>delta) continue;
            xmean4[i]+=xcl*weight; 
            w4[i]+=weight;
       }
       if(w4[i]>0.55) xmean4[i]=xmean4[i]/w4[i];  else {w4[i]=0.; xmean4[i]=0.;} ;
    }
  }
void  CsRwRecons::ec1energy_cor(double xcut, double ycut,  std::list<CsDigit*>& ec1dig, double X, double Y,
                                            const  std::vector<Reco::Cell> & cells, double &Enew) {
         Enew=0.;
         list<CsDigit*>::iterator mdg;
         for( mdg = ec1dig.begin(); mdg != ec1dig.end(); mdg++ ) {
                CsDigit* mydg = *mdg;
                int ncell=mydg->getAddress();
		int digsize=mydg->getDataSize();
		double amp=0.;
		if(digsize==1) amp = mydg->getDatum();
                double xx=cells[ncell].GetX();
                double yy=cells[ncell].GetY();
                if(fabs(X-xx)<xcut&&fabs(Y-yy)<ycut)  {
                         Enew+=amp;
                }
        }
   }
void CsRwRecons::find_best_MCvertex(double x, double y, double zx,  double zy, 
                                    double xtarg,double ytarg, double ztarg, double &alfa,double & beta){
  CsEvent*  event  = CsEvent::Instance();
  list <CsMCVertex*> mcvert=event->getMCVertices();
  list<CsMCVertex*>::iterator vrt;
  double X0=0.,Y0=0.,Z0=0.,rmin=100000000.;
  bool orvertfound = false;
  for( vrt = mcvert.begin(); vrt != mcvert.end(); vrt++ ) {
      CsMCVertex* myvert = *vrt;
//      if(myvert->getGnum()!=2) continue;
      Z0=myvert->getZ();
      if(Z0<9750.||Z0>9850.) continue;
      X0=myvert->getX();
      Y0=myvert->getY();
//      cout<<" MC vertex: X="<<X0<<" Y0="<<Y0<<" Z0="<<Z0<<endl;
      orvertfound = true;
     double alfam=(X0-xtarg)/(Z0-ztarg);
     double betam=(Y0-ytarg)/(Z0-ztarg);
     X0=xtarg+alfam*zx;
     Y0=ytarg+betam*zy;
     double r=sqrt((X0-x)*(X0-x)+(Y0-y)*(Y0-y));
     if(r<rmin)  {rmin=r; alfa=alfam; beta=betam;}
  }
  if(!orvertfound) { cout << " Warning: original MC vertex not found"<<endl; return;}
  }
//
// --------------------------------------------------
double  CsRwRecons::GetLike( int tag, CsTrack *track ){

  double like = -1;
  if(!track->hasRich1Probs())
    return -1;
  else 
    {
      const double *a=track->getRich1Probs();
      int len=sizeof(*a);
      len=21;
      if (tag==0)
	like = track->PionLikelihood();
      else if(tag==1)
	like = track->KaonLikelihood();
      else if(tag==2)
	like = track->ProtonLikelihood();
      else if(tag==3)
	if (len != 21) return -1;
	else       like = track->ElectronLikelihood();
      else if(tag==4)
	if(len  != 21) return -1;
	else       like = track->MuonLikelihood();
      else if(tag==5)
	like = track->BkgLikelihood();
    }
  
  return ( like );
  
}
int CsRwRecons::LikePid(CsTrack *track)
{
  
  int   id  = -1;
  double like[6];

  if(!track->hasRich1Probs())
    return -1;
  const std::vector<CsHelix> h = track->getHelices();
  if(h.size()<1) return -1;
  double Mom=abs(h[0].getCop());
  if(Mom>0.) Mom=1./Mom;
  int NRichInf =21;
  if(Mom >50.)
    return -1;

  float max  = 0;
  if(NRichInf == 15)
    max=1;
    
  
  for (int k =0; k<6; ++k)
    {
     
      if(k!=4)
	if((Mom>8. && k!=3) || ( Mom<8.))
	  {
	    like[k] = GetLike( k, track );
	    if( like[k] > max )
	      {
		max = like[k] ;
		id = k ;
	      }
	  }
    }
  //to avoid cases in which all the likelihood are the same
  int cont=0;
  for (int k =0; k<5; ++k)
    if(like[k]==like[k+1])
      ++cont;
  if (cont>2) id=-1;
  return id;
}
//
CsRwRecons* CsRwRecons::instance_ = 0;
//
CsRwRecons* CsRwRecons::Instance() {
  if( instance_ == 0 ) instance_ = new CsRwRecons();
  return( instance_ );
}
void CsRwRecons::Clear(void){
    int sz=dataout.size();
    if(sz>0){
      for(int i=0; i<sz;i++){
        Reco::CalorimeterParticle* id = dataout[i]->calobj;
        // delete id; // commented out as pointers are saved by CsEvet::AddCalObject() cleaned by CORAL
        delete dataout[i];
      }
      dataout.clear();  
    }
}
//
bool CsRwRecons::gmchk(double Xg, double Yg, double delta) {
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
int CsRwRecons::wrchk(int prj, CsCluster* hit) {
  std::list<CsDigit*> dig =hit->getDigitsList();
  std::list<CsDigit*>:: iterator idig=dig.begin();
  int nwr=(*idig)->getAddress();
  if(prj<4) {
    if(nwr<200) return 0;
    else if(nwr<296) return 1;
    else if(nwr<496) return 2;
    else return 3;
  } else {
    if(nwr<160) return 0;
    else if(nwr<208) return 1;
    else if(nwr<368) return 2;
    else return 3;
  }
}
bool CsRwRecons::htchk(int ipl, double Xg, double Yg, int nz, double delta) {
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
void CsRwRecons::m2g(std::vector<ec1rec>&  ec1,int nh1,int nh2, int nh3, int nh4){
     int ngam=ec1.size();
     if(ngam<2) return;
     for(int i=0;i<ngam-1;i++){
       double X1=ec1[i].X;    
       double Y1=ec1[i].Y;    
       double Z1=ec1[i].Zx;    
       double E1=ec1[i].E; 
       if(E1<2.) continue;  
       int nz1= ec1[i].zone; 
       for(int j=i+1;j<ngam;j++){
          double X2=ec1[j].X;    
          double Y2=ec1[j].Y;    
          double Z2=ec1[j].Zx;    
          double E2=ec1[j].E;
          if(E2<2.) continue;   
          int nz2= ec1[j].zone; 
//          if(nz1!=nz2) continue;
          double sc=(X1*X2+Y1*Y2+Z1*Z2);
          sc=sc/sqrt(X1*X1+Y1*Y1+Z1*Z1)/sqrt(X2*X2+Y2*Y2+Z2*Z2);
          double m2g=sqrt(2.*E1*E2*(1.-sc));
          hi1D[nh1]->Fill(m2g);
          if(nh2>0) hi2D[nh2]->Fill(E1+E2,m2g);
          int nz=nz1; if(E2>E1) nz =nz2;
          if(nh3>0) hi1D[nh3+nz]->Fill(m2g);
          if(nh4>0) {if(ec1[i].nfirst==8&&ec1[j].nfirst==8) hi1D[nh4]->Fill(m2g);}
       }
     }
}
//=============================================================
void CsRwRecons:: PRrecons_consts(void){
cout << setprecision(6);
if(readcal)  cout<<" RW calibration read, Ncalibpoint="<<Ncorpoints<<endl;
else cout<<" No RW calibration read, defaults value will be ussed"<<endl;
if(readcal) {
    cout<< " Eneriges points: ";
   for(int i=0;i<Ncorpoints;i++) cout<<Epoints[i]<<"   "; 
   cout<<endl<< "  Ecorrections:"<<endl;
   for(int j=0; j<3;j++){
          for(int i=0;i<Ncorpoints;i++) cout<<Ecorr_points[j][i]<<"   ";
          cout<<endl;
   } 
   cout<< "  Hits slopes:"<<endl;
   for(int j=0; j<3;j++){
          for(int i=0;i<Ncorpoints;i++) cout<<hits_slp[j][i]<<"   ";
          cout<<endl;
   } 
  }
   cout<< "  DeltaEn: ";
   for(int i=0;i<3;i++) cout<<deltaEn[i]<<"   "; 
   cout<<endl;
   cout<< "  DeltaCR: ";
   for(int i=0;i<3;i++) cout<<deltaCRd[i]<<"   "; 
   cout<<endl<< "  Defaults values: "<<endl;
   cout<<"  Ecorr: ";
   for(int i=0;i<3;i++) cout<<Ecorr[i]<<"   "; 
   cout<<endl<< "  Hits slopes: ";
   for(int i=0;i<3;i++) cout<<hits_slope[i]<<"   "; 
    cout<<endl;
}
//=========================================================
 void CsRwRecons:: corcoeff(double Eg,double wtot, int nz, double& a,double& b){
    if(readcal) {
       double aaf[3]={1.1,0.95,1.05};
       double bbf[3]={0.022,0.022,0.020};
       double E=aaf[nz]*(Eg+wtot*bbf[nz]);
       int k=1;
       if(E<Epoints[0]) {a=Ecorr_points[nz][0];  b=hits_slp[nz][0]; return;}
       else if(E>Epoints[5]) {a=Ecorr_points[nz][5];  b=hits_slp[nz][5]; return;}
       else { for(int i=1; i<6; i++) {
                  if(Epoints[i]>E)  break; 
                  k++;
               }
      }
      a=Ecorr_points[nz][k-1]+(Ecorr_points[nz][k]-Ecorr_points[nz][k-1])*(E-Epoints[k-1])/(Epoints[k]-Epoints[k-1]); 
      b=hits_slp[nz][k-1]+(hits_slp[nz][k]-hits_slp[nz][k-1])*(E-Epoints[k-1])/(Epoints[k]-Epoints[k-1]); 
  } else{
      a=Ecorr[nz];  
      b=hits_slope[nz];}
}
