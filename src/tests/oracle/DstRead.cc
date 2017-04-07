#include "CoralUser.h"
#undef  __STRICT_ANSI__
#include "CsGeom.h"
#include "CsEvent.h"
#include "CsOraStore.h"

#include <TDirectory.h>
#include <TFile.h>
#include <TH2.h>
#include <TVector3.h>
#include <TVector3.h>

using namespace CLHEP;
using namespace std;

const double M_Pi = 0.139567;
const double M_P = 0.9382723;
const double M_K = 0.493677;
const double M_mu = 0.105658357;
const double M_D0_= 1.8645;
const double M_K0_= 0.497672;
static CsOraStore* store;

inline void Set( CsHelix& , TVector3& );
inline bool RICH_D0_OR  ( double m1, double p1, double m2, double p2 );
inline bool RICH_D0_AND ( double m1, double p1, double m2, double p2 );
inline bool RICH_aD0_OR ( double m1, double p1, double m2, double p2 );
inline bool RICH_aD0_AND( double m1, double p1, double m2, double p2 );
inline bool RICH_Phi_OR ( double m1, double p1, double m2, double p2 );
inline bool RICH_Phi_AND( double m1, double p1, double m2, double p2 );
inline bool RICH_LAMBDA_OR ( double m1, double p1, double m2, double p2 );
inline bool RICH_LAMBDA_AND( double m1, double p1, double m2, double p2 );

TDirectory *dM_K0,*dM_L,*dM_Phi,*dM_JP, *dM_D0,*dM_aD0, *dM_Dp,*dM_Dm,*dGeneral;
TH1F *hM_K0[10],*hM_L[10]  ,*hM_Phi[10],*hM_JP[10];
TH1F *hM_D0[10],*hM_aD0[10],*hM_Dp[10] ,*hM_Dm[10];
TH1F* hRICH1[10],*hMonit[12];
TH1F* prof[500];
TH2F* hRICH2[10];
//TH2F* hXY[10];
TH2F* hCluMap;

int NT,NTM,NPV,NSV;

list<CsDetector*> dtt;

bool Boost( const TVector3& v0, double m0, const TVector3& v1, double m1, TVector3& v10 )
{
  double E0=sqrt(m0*m0+v0.Mag2());
  double E1=sqrt(m1*m1+v1.Mag2());
  
  double E10=(E0*E1-v0*v1)/m0;
  v10=v1-(E1+E10)/(E0+m0)*v0;
  
  //std::cout<<" KIN E0="<<E0<<" E1="<<E1<<" E10="<<E10<<std::endl;
  //std::cout<<" KIN v0: "<<v0[0]<<" "<<v0[1]<<" "<<v0[2]<<" "<<std::endl;
  //std::cout<<" KIN v1: "<<v1[0]<<" "<<v1[1]<<" "<<v1[2]<<" "<<std::endl;
  //std::cout<<" KIN v10: "<<v10[0]<<" "<<v10[1]<<" "<<v10[2]<<" "<<std::endl;
  //std::cout<<std::endl;
  
  return true;
}

// --------------------------------------------------

bool CoralUserSetup(int argc, char* argv[])
{
  return true;
}

bool CoralUserInit() {

  store = CsOraStore::Instance();
  NT=NTM=NPV=NSV=0;

  for(int i=0;i<10;i++) {
    hM_K0[i]=0; hM_L[i]=0; hRICH1[i]=0; hRICH2[i]=0; hM_Phi[i]=0; hM_JP[i]=0;
    hM_D0[i]=0; hM_aD0[i]=0; hM_Dp[i]=0; hM_Dm[i]=0; hMonit[i]=0;
    //    hM_D0[i]=0; hM_aD0[i]=0; hM_Dp[i]=0; hM_Dm[i]=0; hXY[i]=0; hMonit[i]=0;
  }
  
  TFile *f = gFile;
  
  dM_K0 = new TDirectory("dM_K0","M_K0");
  dM_K0 -> cd();
  
  hM_K0[0] = new TH1F("hM_K0_0","\\pi^{+}\\pi^{-} invariant mass", 200,0.2,0.9 );
  hM_K0[1] = new TH1F("hM_K0_1","\\pi^{+}\\pi^{-} invariant mass, z>35cm", 200,0.2,0.9 );
  hM_K0[2] = new TH1F("hM_K0_2","\\pi^{+}\\pi^{-} invariant mass, dz>3*sigma", 200,0.2,0.9 );
  hM_K0[3] = new TH1F("hM_K0_3","\\pi^{+}\\pi^{-} invariant mass, dz>10cm", 200,0.2,0.9 );
  
  f->cd();
  dM_L = new TDirectory("dM_L","M_L");
  dM_L -> cd();
  
  hM_L[0] = new TH1F("hM_L_0","p\\pi^{-} invariant mass.", 200,1,1.4 );
  hM_L[1] = new TH1F("hM_L_1","p\\pi^{-} invariant mass, z>35cm.", 200,1,1.4 );
  hM_L[2] = new TH1F("hM_L_2","p\\pi^{-} invariant mass, dz>10cm.", 200,1,1.4 );
  hM_L[3] = new TH1F("hM_L_3","p\\pi^{-} invariant mass. RICH 'OR'", 200,1,1.4 );
  hM_L[4] = new TH1F("hM_L_4","p\\pi^{-} invariant mass, z>35cm. RICH 'OR'", 200,1,1.4 );
  hM_L[5] = new TH1F("hM_L_5","p\\pi^{-} invariant mass, dz>10cm. RICH 'OR'", 200,1,1.4 );
  hM_L[6] = new TH1F("hM_L_6","p\\pi^{-} invariant mass. RICH 'AND'.", 200,1,1.4 );
  hM_L[7] = new TH1F("hM_L_7","p\\pi^{-} invariant mass, z>35cm. RICH 'AND'.", 200,1,1.4 );
  hM_L[8] = new TH1F("hM_L_8","p\\pi^{-} invariant mass, dz>10cm. RICH 'AND'.", 200,1,1.4 );
    
  f->cd();
  dM_Phi = new TDirectory("dM_Phi","M_Phi");
  dM_Phi-> cd();
  
  hM_Phi[0] = new TH1F("hM_Phi_0","K^{+}K^{-} invariant mass", 400,0.8,1.4 );
  hM_Phi[1] = new TH1F("hM_Phi_1","K^{+}K^{-} invariant mass. RICH 'OR'.", 400,0.8,1.4 );
  hM_Phi[2] = new TH1F("hM_Phi_2","K^{+}K^{-} invariant mass. RICH 'AND'.", 400,0.8,1.4 );
   
  f->cd();
  dM_JP = new TDirectory("dM_JP","M_JP");
  dM_JP-> cd();
  
  hM_JP[0]  = new TH1F("hM_JP_0","\\mu^{+}\\mu^{-} invariant mass",400,2,3.5 );
  hM_JP[1]  = new TH1F("hM_JP_1","\\mu^{+}\\mu^{-} invariant mass. p_{t}>0.8 GeV",400,2,3.5 );
  hM_JP[2]  = new TH1F("hM_JP_2","\\mu^{+}\\mu^{-} invariant mass. p_{t}>1.3 GeV",400,2,3.5 );
  
  f->cd();
  dM_D0= new TDirectory("dM_D0","M_D0");
  dM_D0-> cd();
  
  hM_D0[0]  = new TH1F("hM_D0_0","K^{-}\\pi^{+} invariant mass", 400,0.8,2.5 );
  hM_D0[1]  = new TH1F("hM_D0_1","K^{-}\\pi^{+} invariant mass. cos(\\theta^{*})<0.5", 400,0.8,2.5 );
  hM_D0[2]  = new TH1F("hM_D0_2","K^{-}\\pi^{+} invariant mass. cos(\\theta^{*})<0.7", 400,0.8,2.5 );
  hM_D0[3]  = new TH1F("hM_D0_3","K^{-}\\pi^{+} invariant mass", 400,0.8,2.5 );
  hM_D0[4]  = new TH1F("hM_D0_4","K^{-}\\pi^{+} invariant mass. cos(\\theta^{*})<0.5", 400,0.8,2.5 );
  hM_D0[5]  = new TH1F("hM_D0_5","K^{-}\\pi^{+} invariant mass. cos(\\theta^{*})<0.7", 400,0.8,2.5 );
  hM_D0[6]  = new TH1F("hM_D0_6","K^{-}\\pi^{+} invariant mass", 400,0.8,2.5 );
  hM_D0[7]  = new TH1F("hM_D0_7","K^{-}\\pi^{+} invariant mass. cos(\\theta^{*})<0.5", 400,0.8,2.5 );
  hM_D0[8]  = new TH1F("hM_D0_8","K^{-}\\pi^{+} invariant mass. cos(\\theta^{*})<0.7", 400,0.8,2.5 );
  hM_D0[9]  = new TH1F("hM_D0_9","K^{-}\\pi^{+} invariant mass", 400,0.8,2.5 );
  
  f->cd();
  dM_aD0= new TDirectory("dM_aD0","M_aD0");
  dM_aD0-> cd();
  
  hM_aD0[0] = new TH1F("hM_aD0_0","\\pi^{-}K^{+} invariant mass", 400,0.8,2.5 );
  hM_aD0[1] = new TH1F("hM_aD0_1","\\pi^{-}K^{+} invariant mass. cos(\\theta^{*})<0.5", 400,0.8,2.5 );
  hM_aD0[2] = new TH1F("hM_aD0_2","\\pi^{-}K^{+} invariant mass. cos(\\theta^{*})<0.7", 400,0.8,2.5 );
  hM_aD0[3] = new TH1F("hM_aD0_3","\\pi^{-}K^{+} invariant mass", 400,0.8,2.5 );
  hM_aD0[4] = new TH1F("hM_aD0_4","\\pi^{-}K^{+} invariant mass. cos(\\theta^{*})<0.5", 400,0.8,2.5 );
  hM_aD0[5] = new TH1F("hM_aD0_5","\\pi^{-}K^{+} invariant mass. cos(\\theta^{*})<0.7", 400,0.8,2.5 );
  hM_aD0[6] = new TH1F("hM_aD0_6","\\pi^{-}K^{+} invariant mass", 400,0.8,2.5 );
  hM_aD0[7] = new TH1F("hM_aD0_7","\\pi^{-}K^{+} invariant mass. cos(\\theta^{*})<0.5", 400,0.8,2.5 );
  hM_aD0[8] = new TH1F("hM_aD0_8","\\pi^{-}K^{+} invariant mass. cos(\\theta^{*})<0.7", 400,0.8,2.5 );
  hM_aD0[9] = new TH1F("hM_aD0_9","\\pi^{-}K^{+} invariant mass", 400,0.8,2.5 );
  
  f->cd();
  dM_Dp= new TDirectory("dM_Dp","M_Dp");
  dM_Dp-> cd();
  
  hM_Dp[0] = new TH1F("hM_Dp_0","D0\\pi+ invariant mass",  400,1.6,2.5 );
  hM_Dp[1] = new TH1F("hM_Dp_1","D*+ - D0 invariant mass",200,0,0.3 );
  hM_Dp[2] = new TH1F("hM_Dp_2","D0\\pi+ invariant mass",  400,1.6,2.5 );
  hM_Dp[3] = new TH1F("hM_Dp_3","D*+ - D0 invariant mass",200,0,0.3 );
  hM_Dp[4] = new TH1F("hM_Dp_4","D0\\pi+ invariant mass. RICH 'AND'",  400,1.6,2.5 );
  hM_Dp[5] = new TH1F("hM_Dp_5","D*+ - D0 invariant mass. RICH 'AND'",200,0,0.3 );
  hM_Dp[6] = new TH1F("hM_Dp_6","D0\\pi+ invariant mass. RICH 'AND'",  400,1.6,2.5 );
  hM_Dp[7] = new TH1F("hM_Dp_7","D*+ - D0 invariant mass. RICH 'AND'",200,0,0.3 );
  
  hM_Dp[8] = new TH1F("hM_Dp_8","D+ invariant mass",  400,0.5,2.2 );
  
  f->cd();
  dM_Dm= new TDirectory("dM_Dm","M_Dm");
  dM_Dm-> cd();
  
  hM_Dm[0] = new TH1F("hM_Dm_0","aD0\\pi- invariant mass",  400,1.6,2.5 );
  hM_Dm[1] = new TH1F("hM_Dm_1","D*- - aD0 invariant mass",200,0,0.3 );
  hM_Dm[2] = new TH1F("hM_Dm_2","aD0\\pi- invariant mass",  400,1.6,2.5 );
  hM_Dm[3] = new TH1F("hM_Dm_3","D*- - aD0 invariant mass",200,0,0.3 );
  hM_Dm[4] = new TH1F("hM_Dm_4","aD0\\pi- invariant mass. RICH 'AND'",  400,1.6,2.5 );
  hM_Dm[5] = new TH1F("hM_Dm_5","D*- - aD0 invariant mass. RICH 'AND'",200,0,0.3 );
  hM_Dm[6] = new TH1F("hM_Dm_6","aD0\\pi- invariant mass. RICH 'AND'",  400,1.6,2.5 );
  hM_Dm[7] = new TH1F("hM_Dm_7","D*- - aD0 invariant mass. RICH 'AND'",200,0,0.3 );
  
  hM_Dm[8] = new TH1F("hM_Dm_8","D- invariant mass",  400,0.5,2.2 );
  
  f->cd();
  dGeneral= new TDirectory("dGeneral","General");
  dGeneral-> cd();
  
  hRICH1[0] = new TH1F("hM1","Mass of positive particle", 200,0,1.2 );
  hRICH1[1] = new TH1F("hM2","Mass of negative particle", 200,0,1.2 );
  hRICH2[0] = new TH2F("hPB1","\\theta vs mom (+)", 120,0,40, 100,0,70 );
  hRICH2[1] = new TH2F("hPB2","\\theta vs mom (-)", 120,0,40, 100,0,70 );
  
//    hXY[0] = new TH2F("hXY_0","Profile of FI02", 100, -60, 60, 50,-30 ,30  );
//    hXY[1] = new TH2F("hXY_1","Profile of FI03", 100, -30, 30, 50,-30 ,30  );
//    hXY[2] = new TH2F("hXY_2","Profile of GM02",  50,-200,200, 50,-200,200 );
//    hXY[3] = new TH2F("hXY_3","Profile of GM10",  50, -30,370, 50,-200,200 );
//    hXY[4] = new TH2F("hXY_4","Profile of PS01",  50,-900,900, 50,-700,700 );
  
  hMonit[0] = new TH1F("hMonit_0","tracks per event"            , 300,0,30000 );
  hMonit[1] = new TH1F("hMonit_1","tracks (with mom) per event" , 300,0,30000 );
  hMonit[2] = new TH1F("hMonit_2","primary vertices per event"  , 300,0,30000 );
  hMonit[3] = new TH1F("hMonit_3","secondary vertices per event", 300,0,30000 );
  
  hMonit[4] = new TH1F("hMonit_4","Number of tracks " , 100,0,100 );
  hMonit[5] = new TH1F("hMonit_5","Number of tracks before the target" , 100,0,100 );
  hMonit[6] = new TH1F("hMonit_6","Number of tracks before SM1" , 100,0,100 );
  hMonit[7] = new TH1F("hMonit_7","Number of tracks between M1 and M2" , 100,0,100 );
  hMonit[8] = new TH1F("hMonit_8","Number of tracks between M2 and Muon Wal" , 100,0,100 );
  hMonit[9] = new TH1F("hMonit_9","Number of tracks after Muon Wall" , 100,0,100 );
  hMonit[10] = new TH1F("hMonit_10","Total number of clusters per event " , 1000,0,4000 );
  hMonit[11] = new TH1F("hMonit_11","Total number of clusters per plane" , 100,0,100 );
  
  hCluMap = new TH2F("hCluMap","ID vs. clusters size" , 100, 0, 100, 1300, 1, 1300 );

  dtt= CsGeom::Instance()->getDetectors(); 
  list<CsDetector*>::iterator id;

  unsigned int nbins;
  float minx, maxx;
  string det_name;
  string isto_title;
  unsigned int isto_number=0;
  char* isto_name=new char[10];


  for (id=dtt.begin(); id!=dtt.end(); id++){
    det_name = (*id)-> GetTBName();
    if (det_name[4] == 'X' ||  det_name[4] == 'Y' ||
	det_name[4] == 'U' ||  det_name[4] == 'V'){

      sprintf(isto_name, "prof_%d", isto_number);
      isto_title = "Profile of " + det_name;
      double dimension, center;
      nbins = (*id)->getNWir();
      dimension = (nbins - 1) * (*id)->getWirP();
      center = (*id)-> getWirD() + dimension/2; 
      minx = center - dimension/2 - (*id)->getWirP()/2;
      maxx = center + dimension/2 + (*id)->getWirP()/2;
      cout << isto_name << " " << isto_title << endl;
      prof[isto_number] = new TH1F( isto_name,isto_title.c_str(), nbins, minx, maxx );
      isto_number++;
    }
  }

  f->cd();
  return true;
}


// --------------------------------------------------
bool CoralUserEvent() {

  if(store)
    {
      unsigned eventSize = store->getRawBufferLength();
      uint8* buffer = store->rawBuffer();
      cout << "Event size: " << eventSize << endl;
      for(int i = 0; i < 10; i++)
	cout << "byte[" << i << "] = " << (int)buffer[i] << endl;
    }
  unsigned int NEV = CsEvent::Instance()->getNumberOfEvents();
  
  const list<CsVertex*> &vrts = CsEvent::Instance()->getVertices();
  const list<CsTrack* > &Trks = CsEvent::Instance()->getTracks();
  const vector<CsParticle* > &parts = CsEvent::Instance()->getParticles();
  
  
  ///////////////// LOOP OVER PARTICLES //////////////// 
  
  CsVertex *vPr(0);    // pointer to primary vertex
  vector<CsParticle*>::const_iterator ip;
  for( ip=parts.begin(); ip!=parts.end(); ip++ ) {
    if( (*ip)->getType() != CsParticle::SPECIAL ) continue;
    const CsTrack *trk = (*ip)->getTrack();
    if( trk == 0 ) {
      cout<<"ERROR: particle is special but no track!!!"<<endl;
      continue;
    }
    vector<CsHelix> v = trk->getHelices();
    if( v[0].getZ() < 0 ) continue;  // beam
    CsVertex* vrt = trk->getFirstVertex();
    if( vrt==0 || !vrt->isPrimary() ) continue;  // no primary
    vPr = vrt;                       // OK there is primary with mu/mu'
    break;
  }
  

  ///////////////// PRIMARY VERTEX ////////////////

  double nu(0);             // to make cut on z
  double X1(0),Y1(0),Z1(0); // coordinates of primary vertex
  HepMatrix* cov1(0);       // error matrix of primary vertex
  
  if( vPr != 0 ) {
    X1 = vPr->getX();
    Y1 = vPr->getY();
    Z1 = vPr->getZ();
    cov1 = vPr->getCov(0);

    const list<CsTrack*> &trks = vPr->getTracks();
    list<CsTrack*>::const_iterator it = trks.begin();
    Cs3Vector parB,parM;
    if( vPr->getPar( *   it , parB ) && 
        vPr->getPar( *(++it), parM ) ) {
      double pB = parB.getMom();
      double pM = parM.getMom();
      nu = pB - pM;
    }
  }


  ///////////////// LOOP OVER SECONDARY VERTICES ////////////////
  
  list<CsVertex*>::const_iterator iv;
  for( iv=vrts.begin(); iv!=vrts.end(); iv++ ) {
    CsVertex *vrt = (*iv);
    
    if( vrt->isPrimary() ) { NPV++; continue; }
    else                     NSV++;
    
    const list<CsTrack*> &trks = vrt->getTracks();
    
    if( trks.size() != 2 ) {
      cout<<"ERROR: number of tracks is "<<trks.size()<<endl;
      continue;
    }
    
    CsHelix h1( trks.front()->getHelices()[0] );
    CsHelix h2( trks.back() ->getHelices()[0] );
    
    if( h1(0) < 0 || h2(0) < 0 ) continue;  // beam tracks
    if( h1(5) < 0 ) {                       // the first track must be +
      cout<<"ERROR: the first track is negative!"<<endl;
      continue;
    }
    if( h2(5) > 0 ) {                       // the second track must be -
      cout<<"ERROR: the second track is positive!"<<endl;
      continue;
    }
    
    double X = vrt->getX();
    double Y = vrt->getY();
    double Z = vrt->getZ();
    
    CsHelix ht1,ht2;
    h1.Extrapolate(Z,ht1);
    h2.Extrapolate(Z,ht2);
    
    TVector3 v1,v2;
    Set( ht1, v1 );
    Set( ht2, v2 );
    
    TVector3 v0(v1+v2);
    
    
    ///////////// CUTS FOR PRIM VERTEX /////////////// 
    
    HepMatrix* cov = vrt->getCov(0);
    
    double dist(999999);
    double SIGMA(999999);
    double distR(999999);
    if( vPr != 0 ) {
      TVector3 Dist(X-X1,Y-Y1,Z-Z1);
      dist = Dist.Mag();
      SIGMA = 3*sqrt( (*cov)(3,3) + (*cov1)(3,3) );
      distR = dist * Dist.Angle(v0);
    }
    
    
    ///////////// RICH ///////////////
    
    double m1(-1);
    double m2(-1);
    
    if( trks.front()->hasRich1Probs() ) {
      double theta1 = trks.front()->Rich1Theta();
      double mom = v1.Mag();
      if( theta1 < 50.4 && theta1 > 15 && mom > 2.5 && mom < 40 ) {
        double a = 2.314;
        double b = 35.7;
        if( theta1 < a*mom+b ) {
          double n    = 1.00129;
          double beta = 1/cos(theta1*0.001)/n;
          m1=mom * sqrt(1-beta*beta) / beta;
          hRICH1[0]->Fill( m1 );
          hRICH2[0]->Fill( mom, theta1 );
        }
      }
    }
    
    if( trks.back()->hasRich1Probs() ) {
      double theta2 = trks.back() ->Rich1Theta();
      double mom = v2.Mag();
      if( theta2 < 50.4 && theta2 > 15 && mom > 2.5 && mom < 40 ) {
        double a = 2.314;
        double b = 35.7;
        if( theta2 < a*mom+b ) {
          double n    = 1.00129;
          double beta = 1/cos(theta2*0.001)/n;
          m2=mom * sqrt(1-beta*beta) / beta;
          hRICH1[1]->Fill( m2 );
          hRICH2[1]->Fill( mom, theta2 );
        }
      }
    }
    
    
    
    ///////////// INVARIANT MASS /////////////
    
    double pt  = v1.Pt( v0 );
    double pl1 = v1 * v0 / v0.Mag();
    double pl2 = v2 * v0 / v0.Mag();
    double Alpha = (pl1-pl2)/(pl1+pl2);
    
    double E1,E2,E0,M_K0,M_L,M_Phi,M_JP,M_D0,M_aD0,M_Dp,M_Dm;
    
    TVector3 vKD0,vPiD0;
    double CosKD0,CosPiD0;



    ///////////// K0 ///////////////
    
    E1 = sqrt(M_Pi*M_Pi+v1.Mag2());
    E2 = sqrt(M_Pi*M_Pi+v2.Mag2());
    E0 = E1 + E2;
    M_K0 = sqrt(E0*E0-v0.Mag2());
    
    if( pt > 0.08 ) {
      hM_K0[0]->Fill( M_K0 );
      if( Z>350 ) hM_K0[1]->Fill( M_K0 );
      if( vPr != 0 && distR < 2 ) {
        if( dist > SIGMA ) hM_K0[2]->Fill( M_K0 );
        if( dist > 100   ) hM_K0[3]->Fill( M_K0 );
      }  
    }


    ///////////// D0 /////////////// 
    
    E1 = sqrt(M_Pi*M_Pi+v1.Mag2());
    E2 = sqrt(M_K *M_K +v2.Mag2());
    E0 = E1 + E2;
    M_D0 = sqrt(E0*E0-v0.Mag2());

    Boost(v0,M_D0, v2,M_K, vKD0 );
    Boost(v0,M_D0, v1,M_Pi,vPiD0);
    CosKD0 = vKD0 .Unit()*v0.Unit();
    CosPiD0= vPiD0.Unit()*v0.Unit();

    hM_D0[0]->Fill( M_D0 );
    if( fabs( CosKD0 ) < 0.5 && fabs( CosPiD0 ) < 0.5 ) hM_D0[1]->Fill( M_D0 );
    if( fabs( CosKD0 ) < 0.7 && fabs( CosPiD0 ) < 0.7 ) hM_D0[2]->Fill( M_D0 );
    
    if( RICH_D0_OR(m1,v1.Mag(),m2,v2.Mag()) ) {
      hM_D0[3]->Fill( M_D0 );
      if( fabs( CosKD0 ) < 0.5 && fabs( CosPiD0 ) < 0.5 ) hM_D0[4]->Fill( M_D0 );
      if( fabs( CosKD0 ) < 0.7 && fabs( CosPiD0 ) < 0.7 ) hM_D0[5]->Fill( M_D0 );
    }

    if( RICH_D0_AND(m1,v1.Mag(),m2,v2.Mag()) ) {
      hM_D0[6]->Fill( M_D0 );
      if( fabs( CosKD0 ) < 0.5 && fabs( CosPiD0 ) < 0.5 ) hM_D0[7]->Fill( M_D0 );
      if( fabs( CosKD0 ) < 0.7 && fabs( CosPiD0 ) < 0.7 ) hM_D0[8]->Fill( M_D0 );
    }


    ///////////// aD0 /////////////// 
    
    E1 = sqrt(M_Pi*M_Pi+v2.Mag2());
    E2 = sqrt(M_K *M_K +v1.Mag2());
    E0 = E1 + E2;
    M_aD0 = sqrt(E0*E0-v0.Mag2());
    
    Boost(v0,M_aD0, v1,M_K, vKD0 );
    Boost(v0,M_aD0, v2,M_Pi,vPiD0);
    CosKD0 = vKD0 .Unit()*v0.Unit();
    CosPiD0= vPiD0.Unit()*v0.Unit();

    hM_aD0[0]->Fill( M_D0 );
    if( fabs( CosKD0 ) < 0.5 && fabs( CosPiD0 ) < 0.5 ) hM_aD0[1]->Fill( M_aD0 );
    if( fabs( CosKD0 ) < 0.7 && fabs( CosPiD0 ) < 0.7 ) hM_aD0[2]->Fill( M_aD0 );
    
    if( RICH_aD0_OR(m1,v1.Mag(),m2,v2.Mag()) ) {
      hM_aD0[3]->Fill( M_D0 );
      if( fabs( CosKD0 ) < 0.5 && fabs( CosPiD0 ) < 0.5 ) hM_aD0[4]->Fill( M_aD0 );
      if( fabs( CosKD0 ) < 0.7 && fabs( CosPiD0 ) < 0.7 ) hM_aD0[5]->Fill( M_aD0 );
    }

    if( RICH_aD0_AND(m1,v1.Mag(),m2,v2.Mag()) ) {
      hM_aD0[6]->Fill( M_D0 );
      if( fabs( CosKD0 ) < 0.5 && fabs( CosPiD0 ) < 0.5 ) hM_aD0[7]->Fill( M_aD0 );
      if( fabs( CosKD0 ) < 0.7 && fabs( CosPiD0 ) < 0.7 ) hM_aD0[8]->Fill( M_aD0 );
    }

    
    ///////////// D*+ D*- D+ D- ///////////////

    list<CsTrack*>::const_iterator it;
    for( it=Trks.begin(); it!=Trks.end(); it++ ) {
      if( trks.front() == (*it) || trks.back() == (*it) ) continue;
      const vector<CsHelix> &v = (*it)->getHelices();
      if( v[0].getZ() > 3500 || v[0].getZ() < 0) continue;
      if( v[0].getCop() == 0   ) continue;
      if( fabs(1/v[0].getCop()) > 20 ) continue;
      CsHelix hV;
      v[0].Extrapolate(Z,hV);
      double r = sqrt( (hV(1)-X)*(hV(1)-X) + (hV(2)-Y)*(hV(2)-Y) );
      if( r > 15 ) continue;
      
      TVector3 vPi;
      Set( hV, vPi );
      TVector3 vDst = v0 + vPi;
      
            
      if( hV.getCop() > 0 ) {
        
        ///////////// D*+ ///////////////
      
        E1 = sqrt(M_D0*M_D0 + v0. Mag2());
        E2 = sqrt(M_Pi*M_Pi + vPi.Mag2());
        E0 = E1 + E2;
        M_Dp = sqrt(E0*E0-vDst.Mag2());
      
        if( fabs( CosKD0 ) < 0.5 && fabs( CosPiD0 ) < 0.5 ) {
          hM_Dp[0]->Fill( M_Dp );
	  hM_Dp[1]->Fill( M_Dp - M_D0 );
	  if( fabs(M_D0-M_D0_) < 0.05 ) {
	    hM_Dp[2]->Fill( M_Dp );
	    hM_Dp[3]->Fill( M_Dp - M_D0 );
	  }
        }
        
        if( fabs( CosKD0 ) < 0.5 && fabs( CosPiD0 ) < 0.5 &&
            RICH_D0_AND(m1,v1.Mag(),m2,v2.Mag()) ) {
          hM_Dp[4]->Fill( M_Dp );
	  hM_Dp[5]->Fill( M_Dp - M_D0 );
	  if( fabs(M_D0-M_D0_) < 0.05 ) {
	    hM_Dp[6]->Fill( M_Dp );
	    hM_Dp[7]->Fill( M_Dp - M_D0 );
	  }
        }
	
	///////////// D+ /////////////// 
	
	if( fabs(M_K0-M_K0_) < 0.02 && pt > 0.08 ) {
	  
	  E1 = sqrt(M_K0_*M_K0_+ v0. Mag2());
          E2 = sqrt(M_Pi *M_Pi + vPi.Mag2());
          E0 = E1 + E2;
          M_Dp = sqrt(E0*E0-vDst.Mag2());
      
	  hM_Dp[8]->Fill( M_Dp );
	}
	
      }
      
      
      if( hV.getCop() < 0 ) {
      
        ///////////// D*- ///////////////
      
        E1 = sqrt(M_aD0*M_aD0 + v0. Mag2());
        E2 = sqrt(M_Pi *M_Pi  + vPi.Mag2());
        E0 = E1 + E2;
        M_Dm = sqrt(E0*E0-vDst.Mag2());
      
        if( fabs( CosKD0 ) < 0.5 && fabs( CosPiD0 ) < 0.5 ) {
          hM_Dm[0]->Fill( M_Dm );
	  hM_Dm[1]->Fill( M_Dm - M_aD0 );
	  if( fabs(M_aD0-M_D0_) < 0.05 ) {
	    hM_Dm[2]->Fill( M_Dm );
	    hM_Dm[3]->Fill( M_Dm - M_aD0 );
	  }
        }
        if( fabs( CosKD0 ) < 0.5 && fabs( CosPiD0 ) < 0.5 &&
            RICH_aD0_AND(m1,v1.Mag(),m2,v2.Mag()) ) {
          hM_Dm[4]->Fill( M_Dm );
	  hM_Dm[5]->Fill( M_Dm - M_aD0 );
	  if( fabs(M_aD0-M_D0_) < 0.05 ) {
	    hM_Dm[6]->Fill( M_Dm );
	    hM_Dm[7]->Fill( M_Dm - M_aD0 );
	  }
        }
        
	///////////// D- ///////////////
	
	if( fabs(M_K0-M_K0_) < 0.02 && pt > 0.08 ) {
	  E1 = sqrt(M_K0_*M_K0_+ v0. Mag2());
          E2 = sqrt(M_Pi *M_Pi + vPi.Mag2());
          E0 = E1 + E2;
          M_Dm = sqrt(E0*E0-vDst.Mag2());
      
	  hM_Dm[8]->Fill( M_Dm );
	}
	
      }
      
    }

    
    
    ///////////// LAMBDA /////////////// 
    
    E1 = sqrt(M_P *M_P +v1.Mag2());
    E2 = sqrt(M_Pi*M_Pi+v2.Mag2());
    E0 = E1 + E2;
    
    M_L = sqrt(E0*E0-v0.Mag2());
    
    if( pt > 0.03 ) {
      
      hM_L[0]->Fill( M_L );
      if( pt < 0.2 && Alpha > 0.3 ) {
        if( Z > 350 ) hM_L[1]->Fill( M_L );
	if( vPr != 0 && distR < 5 && dist > 100 ) hM_L[2]->Fill( M_L );	
      }
      
      if( RICH_LAMBDA_OR(m1,v1.Mag(),m2,v2.Mag()) ) {
        hM_L[3]->Fill( M_L );
        if( pt < 0.2 && Alpha > 0.3 ) {
          if( Z > 350 ) hM_L[4]->Fill( M_L );
	  if( vPr != 0 && distR < 5 && dist > 100 ) hM_L[5]->Fill( M_L );	
        }
      }
      
      if( RICH_LAMBDA_AND(m1,v1.Mag(),m2,v2.Mag()) ) {
        hM_L[6]->Fill( M_L );
        if( pt < 0.2 && Alpha > 0.3 ) {
          if( Z > 350 ) hM_L[7]->Fill( M_L );
	  if( vPr != 0 && distR < 5 && dist > 100 ) hM_L[8]->Fill( M_L );	
        }
      }
      
    }
    
    ///////////// Phi //////////////
    
    E1 = sqrt(M_K*M_K+v1.Mag2());
    E2 = sqrt(M_K*M_K+v2.Mag2());
    E0 = E1 + E2;
    M_Phi = sqrt(E0*E0-v0.Mag2());
    
    if( pt > 0.02 ) hM_Phi[0]->Fill( M_Phi );
    
    if( RICH_Phi_OR (m1,v1.Mag(),m2,v2.Mag()) && pt > 0.02 ) hM_Phi[1]->Fill( M_Phi );
    if( RICH_Phi_AND(m1,v1.Mag(),m2,v2.Mag()) && pt > 0.02 ) hM_Phi[2]->Fill( M_Phi );

    
    ///////////// JP ////////////// 
    
    E1 = sqrt(M_mu*M_mu+v1.Mag2());
    E2 = sqrt(M_mu*M_mu+v2.Mag2());
    E0 = E1 + E2;
    
    M_JP = sqrt(E0*E0-v0.Mag2());
    
    hM_JP[0]->Fill( M_JP );
    if( pt > 0.8 ) hM_JP[1]->Fill( M_JP );
    if( pt > 1.2 ) hM_JP[2]->Fill( M_JP );
  }
  
  
  
  ///////////////// LOOP OVER TRACKS //////////////// 
  
  int NZ1(0),NZ2(0),NZ3(0),NZ4(0),NZ5(0);

  list<CsTrack*>::const_iterator it;
  for( it=Trks.begin(); it!=Trks.end(); it++ ) {
    CsTrack *trk = (*it);
    
    const vector<CsHelix> v = trk->getHelices();
    
    if( v.size()!=0 && v[0].getCop()!=0 ) NTM++;   // number of tracks with momentum
    
    ////////// number of tracks per zone ////////
    
    const list<CsZone*> zones = trk->getZones ();
    list<CsZone*>::const_iterator iz;
    for( iz=zones.begin(); iz!=zones.end(); iz++ ) {
      const string &name = (*iz)->getName();
      if( name == "before the target"        ) NZ1++;
      if( name == "before M1"                ) NZ2++;
      if( name == "between M1 and M2"        ) NZ3++;
      if( name == "between M2 and Muon Wall" ) NZ4++;
      if( name == "after Muon Wall"          ) NZ5++;
    }
    
  }
  
  
  /////////////// MONITORING //////////////////
  
  hMonit[4]->Fill( Trks.size() );
  hMonit[5]->Fill( NZ1 );
  hMonit[6]->Fill( NZ2 );
  hMonit[7]->Fill( NZ3 );
  hMonit[8]->Fill( NZ4 );
  hMonit[9]->Fill( NZ5 );
 
  const int NE = 100;
  NT += Trks.size();
  
  if( NEV%NE == 0 ) {
    hMonit[0]->Fill( NEV-0.1, NT /float(NE) );
    hMonit[1]->Fill( NEV-0.1, NTM/float(NE) );
    hMonit[2]->Fill( NEV-0.1, NPV/float(NE) );
    hMonit[3]->Fill( NEV-0.1, NSV/float(NE) );
    NT=NTM=NPV=NSV=0;
  }


  if( NEV%25 == 0 ) {
    unsigned int total_number_of_cluster_in_event=0;
    unsigned int total_number_of_cluster_per_plane=0;
    
    unsigned int isto_number=0;
    for( list<CsDetector*>::iterator k=dtt.begin();
	 k!=dtt.end(); k++ ) {
     
      list<CsCluster*> clus = (*k)->getMyClusters();
      list<CsCluster*>::iterator clu;
      string det_name = (*k)-> GetTBName();
      if (det_name[4] == 'X' ||  det_name[4] == 'Y' ||
	  det_name[4] == 'U' ||  det_name[4] == 'V'){
	
	for (clu=clus.begin(); clu!=clus.end(); clu++){
	  prof[isto_number]->Fill((*clu)->getU());
	}
	
	isto_number++;
      }
      total_number_of_cluster_in_event += clus.size();
      total_number_of_cluster_per_plane = clus.size();
      hMonit[11]->Fill( total_number_of_cluster_per_plane );
      hCluMap->Fill( total_number_of_cluster_per_plane, (*k)->GetID().GetNumber() );
    }
    hMonit[10]->Fill( total_number_of_cluster_in_event );
  }
 
  
  return true;
  
}


// --------------------------------------------------
// This fuction is called bevore Coral end. Put here your
// final lines of code...
bool CoralUserEnd() {
  return true;
}



//////////////////////////////////////

inline void Set( CsHelix &ht, TVector3& v )
{
 const double Cop = ht.getCop();
 const double dXdZ= ht.getDXDZ();
 const double dYdZ= ht.getDYDZ();
 
 v[2] = 1/fabs(Cop)/sqrt(1+dXdZ*dXdZ+dYdZ*dYdZ);
 v[0] = v[2] * dXdZ;
 v[1] = v[2] * dYdZ;
}


////////////// RICH CHECKS //////////////

const double CUTPi = 0.1;
const double CUTK  = 0.08;
const double CUTP  = 0.1;

const double THPi = 2.5;
const double THK  = 11;
const double THP  = 19;

/////// D0 ///////

inline bool RICH_D0_OR( double m1, double p1, double m2, double p2 )
{
  if( m1 == -1 && m2 == -1 ) return false;

  if( m1 != -1 && p1>THPi ) {  // pi+ 
    if( fabs(m1-M_Pi) < CUTPi  ) return true;
  }

  if( m2 != -1 && p2>THPi ) {  // K-
    if( p2>THK ) { if( fabs(m2-M_K ) < CUTK  ) return true; }
    else         { if( fabs(m2-M_Pi) > CUTPi ) return true; }
  }

  return false;
}

inline bool RICH_D0_AND( double m1, double p1, double m2, double p2 )
{
  if( m1 == -1  || m2 == -1  ) return false;
  if( p1 < THPi || p2 < THPi ) return false;

  // pi+ 
  if( fabs(m1-M_Pi) > CUTPi  ) return false;

  // K-
  if( p2>THK ) { if( fabs(m2-M_K ) > CUTK  ) return false; }
  else         { if( fabs(m2-M_Pi) < CUTPi ) return false; }
  
  return true;
}


/////// aD0 ///////

inline bool RICH_aD0_OR( double m1, double p1, double m2, double p2 )
{
  if( m1 == -1 && m2 == -1 ) return false;

  if( m1 != -1 && p1>THPi ) {  // K+
    if( p1>THK ) { if( fabs(m1-M_K ) < CUTK  ) return true; }
    else         { if( fabs(m1-M_Pi) > CUTPi ) return true; }
  }

  if( m2 != -1 && p2>THPi ) {  // pi- 
    if( fabs(m2-M_Pi) < CUTPi  ) return true;
  }


  return false;
}

inline bool RICH_aD0_AND( double m1, double p1, double m2, double p2 )
{
  if( m1 == -1  || m2 == -1  ) return false;
  if( p1 < THPi || p2 < THPi ) return false;

  // K+
  if( p1>THK ) { if( fabs(m1-M_K ) > CUTK  ) return false; }
  else         { if( fabs(m1-M_Pi) < CUTPi ) return false; }

  // pi-
  if( fabs(m2-M_Pi) > CUTPi  ) return false;
  
  return true;
}


/////// Phi ///////

inline bool RICH_Phi_OR( double m1, double p1, double m2, double p2 )
{
  if( m1 == -1 && m2 == -1 ) return false;

  if( m1 != -1 && p1>THPi ) {  // K+
    if( p1>THK ) { if( fabs(m1-M_K ) < CUTK  ) return true; }
    else         { if( fabs(m1-M_Pi) > CUTPi ) return true; }
  }

  if( m2 != -1 && p2>THPi ) {  // K-
    if( p2>THK ) { if( fabs(m2-M_K ) < CUTK  ) return true; }
    else         { if( fabs(m2-M_Pi) > CUTPi ) return true; }
  }

  return false;
}

inline bool RICH_Phi_AND( double m1, double p1, double m2, double p2 )
{
  if( m1 == -1  || m2 == -1  ) return false;
  if( p1 < THPi || p2 < THPi ) return false;

  // K+
  if( p1>THK ) { if( fabs(m1-M_K ) > CUTK  ) return false; }
  else         { if( fabs(m1-M_Pi) < CUTPi ) return false; }
  
  // K-
  if( p2>THK ) { if( fabs(m2-M_K ) > CUTK  ) return false; }
  else         { if( fabs(m2-M_Pi) < CUTPi ) return false; }
  
  return true;
}


/////// Lambda ///////

inline bool RICH_LAMBDA_OR( double m1, double p1, double m2, double p2 )
{
  if( m1 == -1 && m2 == -1 ) return false;

  if( m1 != -1 && p1>THPi ) {  // P
    if( p1>THP ) { if( fabs(m1-M_P ) < CUTP  ) return true; }
    if( p1>THK ) { if( fabs(m1-M_K ) > CUTK && 
                       fabs(m1-M_Pi) > CUTPi ) return true; }
    else         { if( fabs(m1-M_Pi) > CUTPi ) return true; }
  }
  
  if( m2 != -1 && p2>THPi ) {  // pi-
    if( fabs(m2-M_Pi) < CUTPi  ) return true;
  }

  return false;
}

inline bool RICH_LAMBDA_AND( double m1, double p1, double m2, double p2 )
{
  if( m1 == -1 && m2 == -1 ) return false;
  if( p1 < THPi || p2 < THPi ) return false;
  
  // P
  if( p1>THP ) { if( fabs(m1-M_P ) > CUTP  ) return false; }
  if( p1>THK ) { if( fabs(m1-M_K ) < CUTK || 
                     fabs(m1-M_Pi) < CUTPi ) return false; }
  else         { if( fabs(m1-M_Pi) < CUTPi ) return false; }
  
  // pi-
  if( fabs(m2-M_Pi) > CUTPi  ) return false;

  return true;
}
