/*!
   \file    CsKalmanFitting.cc
   \brief   Compass Vertex Parameters Fitting Class.
   \author  Alexandre Korzenev
   \version $Revision: 1.42 $
   \date    $Date: 2010/10/26 12:32:17 $ 

*/

#include "CsOpt.h"
#include "CsRegistry.h" 
#include "CsErrLog.h"
#include "CsInit.h"
#include "CsGeom.h"
#include "CsGeant3.h"
#include "CsEvent.h"
#include "CsKalmanFitting.h"

using namespace std;

CsKalmanFitting::CsKalmanFitting()
{

  CsRegistry reg; reg.EOJRegistration( this ); // Register for "end of job" call

  for (int i = 0; i<20; i++ ) statistics_[i] = 0;

  Schema_ = 0;     // IKF
  string msg, mcexact, cuts;
  CsOpt* opt = CsOpt::Instance();
  int NUMB;

  PrimReduce_ = 0;
  opt->getOpt("CsKalmanFitting", "PrimReduce", PrimReduce_);

  hist_ = 0;
  opt->getOpt("CsKalmanFitting", "Hist", hist_);

  Agm_ = true; // Augmenting track's error to make for ms
  opt->getOpt("CsKalmanFitting", "Agm", NUMB); Agm_ = (bool)NUMB;

  Specials_ = true;
  opt->getOpt("CsKalmanFitting", "Specials", NUMB); Specials_ = (bool)NUMB;

  Covar_ = 2;
  opt->getOpt( "CsKalmanFitting", "Covar", Covar_);
  if( Covar_ == 0 ) msg = "Covariation matrices will not be filled.";
  if( Covar_ == 1 ) msg = "Only coordinate correlations are filled.";
  if( Covar_ == 2 ) msg = "Coordinate and momentum correlations of every track are filled.";
  if( Covar_ == 3 ) msg = "Coordinate covariation, momenta covariation and correlations between them for every track are filled.";
  if( Covar_ == 4 ) msg = "All correlations are filled.";
  CsErrLog::mes( elInfo, msg);

  // ********** REFERENCE PLANES **********
  //   There are 2: "RefPlane", for spectrometer tracks, and "RefBeam". An
  // extra quantity, viz. "RefMargin", modifies the referencing. Cf.
  // "CsKalmanFitting::setTrks".
  //                                    ***** DEFAULT: NO REFERENCE
  //   Slightly more than 1e6 cm so that one can later require the RefXXX_'s
  // to be < 1e6 to enable the referencing. "RefMargin" defaults to 0,
  // meaning it's disregarded.
  RefPlane_=RefBeam_ = 10*1.1e6; RefMargin_ = 0;
  if (opt->getOpt("CsKalmanFitting","RefPlane",RefPlane_)) {
    //                                  ***** OPTION: SET REFERENCE
    CsErrLog::msg(elInfo,__FILE__,__LINE__,"Reference Plane @ z = %.2f mm",
		  RefPlane_);
    if (opt->getOpt("CsKalmanFitting","RefMargin",RefMargin_)) {
      CsErrLog::msg(elInfo,__FILE__,__LINE__,"Reference Margin = %.2f mm",
		    RefMargin_);
    }
  }
  if (opt->getOpt("CsKalmanFitting","RefBeam",RefBeam_)) {
    //                                  ***** OPTION: SET REFERENCE for BEAM
    CsErrLog::msg(elInfo,__FILE__,__LINE__,"Beam Reference @ z = %.2f mm",
		  RefBeam_);
  }
  RefPlane_ /= 10; RefBeam_ /= 10; RefMargin_ /= 10;    // ***** Coral->Traffic

  ELcorr_ = 1.0;
  opt->getOpt("CsKalmanFitting", "ELcorr", ELcorr_);
  
  for(int i=0; i<6; i++) Print_[i] = 0;
  vector<int> pr;
  opt->getOpt( "CsKalmanFitting", "Print", pr);
  for(unsigned int i=0;i<6&&i<pr.size();i++) Print_[i]= pr[i];

  list<string> Cuts;

  DirChi2Trk_=5;
  InvChi2Trk_=15;

  if( opt->getOpt("CsKalmanFitting","CUTS",Cuts) ){
    list<string>::iterator Is;
    int i = 0;
    for( Is=Cuts.begin(); Is!=Cuts.end(); Is++, i++ ) {
      if( i == 0 ) {
        istringstream( (*Is) ) >> DirChi2Trk_;
      }
      else if( i == 1 ) {
        istringstream( (*Is) ) >> InvChi2Trk_;
      }
    }
  }
  CsErrLog::msg(elInfo,__FILE__,__LINE__,
		"Following cuts were applied: %.2f %.2f",
		DirChi2Trk_,InvChi2Trk_);
  
  ScalMS_ = 1;
  opt->getOpt("CsKalmanFitting","ScalMS",  ScalMS_);
  
  Chi2VCut_ = 10;
  opt->getOpt("CsKalmanFitting","Chi2VCut",Chi2VCut_);
  
  DistInit_ = 30;
  opt->getOpt("CsKalmanFitting","DistInit",DistInit_);
  DistInit_ /= 10;

  double angle = 2;
  opt->getOpt("CsKalmanFitting","AnglePrim",angle);
  CosPrimMM_ = cos( angle * 3.141593 / 180 );
  CsErrLog::msg(elInfo,__FILE__,__LINE__,
		"Angle between beam and scattered mu < %.2f => Inverse KF for primary Vertex Reco",angle);

  // Best primary Vertex
  BpVDY_ = false;
  opt->getOpt("CsKalmanFitting","BpVDY",NUMB); BpVDY_ = (bool)NUMB;
  BpVTimeW_ = 0;
  opt->getOpt("CsKalmanFitting","BpVTimeW",BpVTimeW_);

  // Cut on bestVertex beam track time conditioning the request for re-tracking.
  rebuildTracksCut_ = .5;
  opt->getOpt("CsKalmanFitting","Retrack",rebuildTracksCut_);

  for (int i = 0; i<5; i++) { hPullsMC[i] = 0; hPullsSm[i] = 0; }

  BeamP0_ = 160 /* GeV */; // Only for histo'ing purposes: need no CsErrLog'ing
  opt->getOpt("CsKalmanFitting","BeamP0",BeamP0_);

  if (hist_) {         // ***** HISTOGRAM INSTANTIATION *****
    const double M_phi = 1.019456;

    if (CsInit::Instance()->IsAMonteCarloJob()) {
      CsHistograms::SetCurrentPath("/CsKalmanFitting/PullsMC");
      hPullsMC[0] = new CsHist1D("hPullX"   ,"MC Pulls: x"    ,100,-10,10);
      hPullsMC[1] = new CsHist1D("hPullY"   ,"MC Pulls: y"    ,100,-10,10);
      hPullsMC[2] = new CsHist1D("hPulldXdZ","MC Pulls: dX/dZ",100,-10,10);
      hPullsMC[3] = new CsHist1D("hPulldYdZ","MC Pulls: dX/dZ",100,-10,10);
      hPullsMC[4] = new CsHist1D("hPullCop" ,"MC Pulls: q/p"  ,100,-10,10);
      
      CsHistograms::SetCurrentPath("/CsKalmanFitting/PullsMCvrt");
      hPullsMCvrt[0] = new CsHist1D("hPullX","Vertex MC Pulls: x",100,-10,10);
      hPullsMCvrt[1] = new CsHist1D("hPullY","Vertex MC Pulls: y",100,-10,10);
      hPullsMCvrt[2] = new CsHist1D("hPullZ","Vertex MC Pulls: z",100,-10,10);

      CsHistograms::SetCurrentPath("/CsKalmanFitting/Delta");
      hDeltaKalman[0] = new CsHist1D("DKX","Kalman: DeltaX (MC - rec)",100,-1.5,1.5);
      hDeltaKalman[1] = new CsHist1D("DKY","Kalman: DeltaY (MC - rec)",100,-1.5,1.5);
      hDeltaKalman[2] = new CsHist1D("DKZ","Kalman: DeltaZ (MC - rec)",100,-100,100);
      
      CsHistograms::SetCurrentPath("/CsKalmanFitting/Secondaries");
      hSecDeltaKalman[0] = new CsHist1D("DKX","Kalman: DeltaX (MC - rec)",100,-0.8,0.8);
      hSecDeltaKalman[1] = new CsHist1D("DKY","Kalman: DeltaY (MC - rec)",100,-0.8,0.8);
      hSecDeltaKalman[2] = new CsHist1D("DKZ","Kalman: DeltaZ (MC - rec)",200,-16 ,16 );

      hSecPullsMCvrt[0] = new CsHist1D("hPullXmu","Vertex MC Pulls: x (scattered muon is found)",100,-10,10);
      hSecPullsMCvrt[1] = new CsHist1D("hPullYmu","Vertex MC Pulls: y (scattered muon is found)",100,-10,10);
      hSecPullsMCvrt[2] = new CsHist1D("hPullZmu","Vertex MC Pulls: z (scattered muon is found)",100,-10,10);
      
    }

    CsHistograms::SetCurrentPath("/CsKalmanFitting/Physics");
    //char hTitle[] = "K#pi  -  -2.8<D^{o}#pi-D^{*}<3.2MeV  ";
    char hTitle[] =
      /* */ "#bar{p}#pi^{+}  -  pT>100MeV,Z>100cm,vv/#deltavv>16,c#theta>.99990  ";
    sprintf(hTitle,
	    "p#pi^{-}  -  pT>%.0fMeV,vv/#deltavv>%d,c#theta>%.5f",
	    V0pTCut*1000,V0vvCut,V0cthCut);
    hm_ppi = new CsHist1D("hm_ppi", hTitle,100,M_Lam-dm_Lam,M_Lam+dm_Lam);
    sprintf(hTitle,
	    "#bar{p}#pi^{+}  -  pT>%.0fMeV,vv/#deltavv>%d,c#theta>%.5f",
	    V0pTCut*1000,V0vvCut,V0cthCut);
    hm_pip =  new CsHist1D("hm_appi",hTitle,100,M_Lam-dm_Lam,M_Lam+dm_Lam);
    sprintf(hTitle,
	    "#pi^{+}#pi^{-}  -  pT>%.0fMeV,vv/#deltavv>%d,c#theta>%.5f",
	    V0pTCut*1000,V0vvCut,V0cthCut);
    hm_pipi = new CsHist1D("hm_pipi",hTitle,100,M_K0 -dm_K0 ,M_K0 +dm_K0);

    hdE =     new CsHist1D("hdE","dE  -  (#mup-#muK^{+}K^{-}-M_{p})/2M_{p})",
			 100,-5,10);
    hm_KK  =  new CsHist1D("hm_KK","K^{+}K^{-}",100,2*M_K-.02,M_phi+dm_phi);
    sprintf(hTitle,"#mu^{+}#mu^{-}  -  #muID");
    hm_mumu = new CsHist1D("hm_mumu",hTitle,100,1.5,6.5);
    sprintf(hTitle,"K#pi  -  %.1f<D^{o}#pi-D^{*}<%.1fMeV",
	    DSLowCut*1000,DSUpCut*1000);
    hm_Kpi =  new CsHist1D("hm_Kpi",hTitle,80,M_D0+D0LowCut,M_D0+D0UpCut);
    hm_pi3m = new CsHist1D("hm_pi3m","3#pi missing M^2",160,-15,17);
    hm_pi3e = new CsHist1D("hm_pi3e","3#pi missing E",   80, -8, 8);
    hm_pi5m = new CsHist1D("hm_pi5m","5#pi missing M^2",160,-15,17);
    hm_pi5e = new CsHist1D("hm_pi5e","5#pi missing E",   80, -8, 8);
    CsHistograms::SetCurrentPath("/");
    
  }

}




bool CsKalmanFitting::end()
{
  unsigned int N = CsEvent::Instance()->getNumberOfProcessedEvents();
  cout<<setprecision(1)<<endl;
  cout<<"/------------/ Vertex Kalman Filter statistics \\---------------\\"<<endl;
  cout<<"| The number of events (def)                   "<<setw(7)<<      N       <<" "<<setw(6)<<100<<"% |\n";
  cout<<"| Was called                                   "<<setw(7)<<statistics_[0]<<" "<<setw(6)<< 100*float(statistics_[0]) /float(N)<<"% |\n";

  cout<<"| Primary was found                            "<<setw(7)<<statistics_[10]<<" "<<setw(6)<<100*float(statistics_[10])/float(N)<<"% |\n";
  cout<<"| Primary with mu'                             "<<setw(7)<<statistics_[12]<<" "<<setw(6)<<100*float(statistics_[12])/float(N)<<"% |\n";
  cout<<"| Primary with BMS                             "<<setw(7)<<statistics_[1]<<" "<<setw(6)<< 100*float(statistics_[1]) /float(N)<<"% |\n";
  cout<<"| Primary with mu' and BMS                     "<<setw(7)<<statistics_[2]<<" "<<setw(6)<< 100*float(statistics_[2]) /float(N)<<"% |\n";
  cout<<"| Primary with mu' and BMS OK                  "<<setw(7)<<statistics_[8]<<" "<<setw(6)<< 100*float(statistics_[8]) /float(N)<<"% |\n";

  cout<<"|--------------------------------------------------------------|\n";
  cout<<"| Less than 2 tracks left in vertex            "<<setw(7)<<statistics_[3]<<" "<<setw(6)<<100*float(statistics_[3])/float(N)<<"% |\n";
  cout<<"| Rejection by Chi2/ndf cut                    "<<setw(7)<<statistics_[4]<<" "<<setw(6)<<100*float(statistics_[4])/float(N)<<"% |\n";
  cout<<"|--------------------------------------------------------------|\n";
  cout<<"| # of events with #prim vert > 1              "<<setw(7)<<statistics_[13]<<" "<<setw(6)<<100*float(statistics_[13])/float(N)<<"% |\n";
  cout<<"| Events which profit from rescu procedure     "<<setw(7)<<statistics_[14]<<" "<<setw(6)<<100*float(statistics_[14])/float(N)<<"% |\n";
  cout<<"|--------------------------------------------------------------|\n";
  cout<<"| Average # of tracks in primary vertex        "<<setw(7)<<float(statistics_[6])/float(statistics_[10])<<"         |\n";
  cout<<setprecision(1);
  cout<<"| Number of secondary vertices                 "<<setw(7)<<statistics_[5]<<" "<<setw(6)<<100*float(statistics_[5])/float(N)<<"% |\n";
  cout<<"|--------------------------------------------------------------|\n";
  double TT = float(statistics_[15]);
  cout<<"| Total time spent per event                   "<<setw(7)<<setprecision(3)<<float(statistics_[15])/1000.<<" "<<setw(6)<<setprecision(0)<<100*float(statistics_[15])/TT<<"% |\n";
  cout<<"\\-----------\\______________________________________/-----------/"<<endl;
  cout<<endl;
  
  CsHistograms::SetCurrentPath("/CsKalmanFitting/");
  CsHist1S *hStat = new CsHist1S( "hStat","Statistics", 20,0,20 );
  CsHistograms::SetCurrentPath("/");
  
  hStat->SetBinContent(0,N);
  for(int i=1;i<20;i++) hStat->SetBinContent(i,statistics_[i-1]);
  
  return true;
}
