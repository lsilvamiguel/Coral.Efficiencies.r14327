/*!
   \file    CsAverPattern.cc
   \brief   Compass Vertex Pattern Class.
   \author  Alexandre Korzenev
   \version $Revision: 1.39 $
   \date    $Date: 2010/11/11 01:09:30 $ 

*/

#include "CsOpt.h"
#include "CsRegistry.h"
#include "CsErrLog.h"
#include "CsMCUtils.h"
#include "CsGeant3.h"
#include "CsEvent.h"
#include "CsAverPattern.h"

using namespace std;

//Constructor
CsAverPattern::CsAverPattern() 
{ 
  // Register for "end of job" call
  CsRegistry reg;
  reg.EOJRegistration( this );

  for( int i=0; i<20; i++ ) statistics_[i]=0;

  string msg,mcexact,cuts;
  CsOpt* opt = CsOpt::Instance();
  int NUMB;

  hist_ = 0;
  msg = "Histogram regime is ";
  opt->getOpt("CsAverPattern", "Hist", hist_);
  if( hist_ == 0 ) msg += "OFF";
  else             msg += "ON";
  CsErrLog::Instance()->mes( elInfo, msg);

  findPrim_ = 1;
  if (opt->getOpt("CsAverPattern", "findPrim", NUMB)) findPrim_ = (bool)NUMB;

  // Uncertainty on helix Distance in "FindPrimary":
  // Do (default) or do not take material (be it from map or ROOTG) into account
  DDistUseMMap_ = true;
  if (opt->getOpt("CsKalmanFitting", "DDistUseMMap", NUMB)) {
    DDistUseMMap_ = (bool)NUMB;
    CsErrLog::mes(elInfo,
		  "Take material (MatMap or ROOTG) into account in helix Distance uncertainty");
  }

  findSec_ = 0;
  msg = "Secondaries regime is ";
  if (opt->getOpt("CsAverPattern", "findSec", NUMB)) findSec_ = (bool)NUMB;
  if( findSec_ == 0 ) msg += "OFF";
  else             msg += "ON";
  CsErrLog::Instance()->mes( elInfo, msg);

  refitTracks_ = 0;
  if (opt->getOpt("CsAverPattern","Refit",NUMB)) refitTracks_ = (bool)NUMB;
  CsErrLog::msg(elInfo,__FILE__,__LINE__,
		"Refit tracks is %s",refitTracks_ ? "ON" : "OFF");

  rebuildTracks_ = 0;
  if (opt->getOpt("CsAverPattern","Retrack",rebuildTracksCut_) &&
      rebuildTracksCut_) {
    rebuildTracks_ = true;
    CsErrLog::msg(elInfo,__FILE__,__LINE__,
  "(requesting a) Re-tracking is ON w/ cut on beam track time: |t|<%f",
		  rebuildTracksCut_);
  }
  else
    CsErrLog::mes(elInfo,"(requesting a) Re-tracking is OFF");

  SecDist_ = 1;
  if( opt->getOpt("CsAverPattern", "SecDist", SecDist_ ) ) 
    SecDist_ /= 10;

  SecXX0_ = 70;
  opt->getOpt("CsAverPattern", "SecXX0", SecXX0_ );
  
  TimePrimCut_ = 100000;
  opt->getOpt("CsAverPattern", "TimePrimCut", TimePrimCut_ );
  
  TimeSecCut_ = 10000;
  opt->getOpt("CsAverPattern", "TimeSecCut", TimeSecCut_ );

  mode_ = 0;
  opt->getOpt( "CsAverPattern", "Mode", mode_);
  if( mode_ == 0 ) msg = "Standard mode";
  if( mode_ == 1 ) msg = "Mode when only correct beam track is involved to Kalman fit";
  if( mode_ == 2 ) msg = "Mode when all involved tracks are correct";
  msg += " was chosen for vertex reconstruction.";
  CsErrLog::Instance()->mes( elInfo, msg);

  NUMB = 1;
  opt->getOpt( "CsAverPattern", "NSpec", NUMB); NSpec_ = (bool)NUMB;

  for(int i=0;i<6;i++) Print_[i]=false;
  vector<int> pr;
  opt->getOpt( "CsAverPattern", "Print", pr);
  for(unsigned int i=0;i<6&&i<pr.size();i++) if(pr[i]!=0) Print_[i]=true;
  
  list<string> Cuts;
  LinDist_=100;
  HelDist_=1;
 
  if( opt->getOpt("CsAverPattern","CUTS",Cuts) ){
    list<string>::iterator Is;
    int i = 0;
    for( Is=Cuts.begin(); Is!=Cuts.end(); Is++, i++ ) {
      if( i == 0 ) {
        istringstream( (*Is) ) >> LinDist_;
	LinDist_ /= 10;
      }
      else if( i == 1 ) {
        istringstream( (*Is) ) >> HelDist_;
	//HelDist_ /= 10;
      }
      
    }
  }
  cuts="Following cuts were applied: ";
  char CUT[10];
  sprintf(CUT, "%.2f",LinDist_);
  cuts+=CUT; cuts+=", ";
  sprintf(CUT, "%.2f",HelDist_);
  cuts+=CUT; cuts+=" .";
  CsErrLog::Instance()->mes( elInfo, cuts);

  CnSigmas_[0] = CnSigmas_[1] = 16; // in mm
  CnSigmas_[2] = 2500; // in mm
  if( opt->getOpt("CsAverPattern","CnSigmas",Cuts) ){
    list<string>::iterator Is;
    int i = 0;
    for( Is=Cuts.begin(); Is!=Cuts.end(); Is++, i++ ) {
      if( i == 0 ) {
        istringstream( (*Is) ) >> CnSigmas_[0];
	CnSigmas_[0] = CnSigmas_[0] * CnSigmas_[0]; // in mm
      }
      else if( i == 1 ) {
        istringstream( (*Is) ) >> CnSigmas_[1];
	CnSigmas_[1] = CnSigmas_[1] * CnSigmas_[1]; // in mm
      }
      else if( i == 2 ) {
        istringstream( (*Is) ) >> CnSigmas_[2];
	CnSigmas_[2] = CnSigmas_[2] * CnSigmas_[2]; // in mm
      }
      
    }
  }

  if (!opt->getOpt("CsAverPattern","MomCut",MomCut_)) MomCut_ = 0.1;
  CsErrLog::msg(elInfo,__FILE__,__LINE__,
		"Momentum cut: P > %.2f GeV",MomCut_);
  
  AcceptTill_ = 90;
  opt->getOpt("CsAverPattern","AcceptTill",AcceptTill_);  
  msg="The track is accepted for prefilter if its momentum is smaller then ";
  sprintf(CUT, "%.2f",AcceptTill_);
  msg+=CUT;
  msg+="%";
  CsErrLog::Instance()->mes( elInfo, msg);
  AcceptTill_ /= 100.;

  double angle = .1; // Default 0.1 deg -> ~1.7 mrd
  opt->getOpt("CsAverPattern", "AnglePrim", angle );
  CosPrimMM_ = cos( angle * 3.141593 / 180 );
  CsErrLog::msg(elInfo,__FILE__,__LINE__,
		"AnglePrim: If mu' such that its angle w.r.t. beam < %.2f => it sets initial guess for Z of vertex",angle);
  angle = .002; // In mrd
  CosPrim_ = cos(angle);

  if( !opt->getOpt("CsAverPattern","BeamCand",BeamCand_) )
    BeamCand_=0.4;
  msg="do loop over all beam candidates inside of interval ";
  sprintf(CUT, "%.2f",BeamCand_);
  msg+=CUT;
  CsErrLog::Instance()->mes( elInfo, msg);
  
  if( hist_ ) {

    const double m_LAM = 1.115684;
    const double m_K0  = 0.497672;

    if( opt->getOpt("","Monte Carlo job") ) {
      
      if( findPrim_ ) {
	CsHistograms::SetCurrentPath("/CsAverPattern/Delta");
	hDeltaPref[0] = new CsHist1D("DPX","Prefilter: DeltaX (MC - rec)",100,-1.5,1.5);
	hDeltaPref[1] = new CsHist1D("DPY","Prefilter: DeltaY (MC - rec)",100,-1.5,1.5);
	hDeltaPref[2] = new CsHist1D("DPZ","Prefilter: DeltaZ (MC - rec)",100,-100,100);
      }

      if( findSec_ ) {
	CsHistograms::SetCurrentPath("/CsAverPattern/Secondaries");
	hSecDeltaPref[0] = new CsHist1D("DPX","Prefilter: DeltaX (MC - rec)",100,-0.8,0.8);
	hSecDeltaPref[1] = new CsHist1D("DPY","Prefilter: DeltaY (MC - rec)",100,-0.8,0.8);
	hSecDeltaPref[2] = new CsHist1D("DPZ","Prefilter: DeltaZ (MC - rec)",100,-16 ,16 );
      }

    }

    CsHistograms::SetCurrentPath("/");

  }

  return;
}


//Destructor 
CsAverPattern::~CsAverPattern() {
}

bool CsAverPattern::end()
{
  unsigned int N = CsEvent::Instance()->getNumberOfProcessedEvents();
  cout<<setprecision(1)<<endl;
  cout<<"/----------/ Vertex preliminary filter statistics \\------------\\"<<endl;
  cout<<"| The number of events (def)                   "<<setw(7)<<      N       <<" "<<setw(6)<<100<<"% |"<<endl;
  cout<<"| Was called                                   "<<setw(7)<<statistics_[0] <<" "<<setw(6)<<100*float(statistics_[0]) /float(N)<<"% |"<<endl;
  cout<<"| At least one Primary was found               "<<setw(7)<<statistics_[10]<<" "<<setw(6)<<100*float(statistics_[10])/float(N)<<"% |"<<endl;
  cout<<"| At least one Primary with mu'                "<<setw(7)<<statistics_[12]<<" "<<setw(6)<<100*float(statistics_[12])/float(N)<<"% |"<<endl;
  cout<<"|--------------------------------------------------------------|"<<endl;
  cout<<"| Beam was absent                              "<<setw(7)<<statistics_[2] <<" "<<setw(6)<<100*float(statistics_[2]) /float(N)<<"% |"<<endl;
  if( NSpec_ != 0 )
    cout<<"| Not enough special tracks                    "<<setw(7)<<statistics_[3]<<" "<<setw(6)<<100*float(statistics_[3])/float(N)<<"% |"<<endl;
  cout<<"| Only one track remains after dist cuts       "<<setw(7)<<statistics_[4] <<" "<<setw(6)<<100*float(statistics_[4]) /float(N)<<"% |"<<endl;
  cout<<"|--------------------------------------------------------------|"<<endl;
  cout<<"| # of events with #beams > 1                  "<<setw(7)<<statistics_[13]<<" "<<setw(6)<<100*float(statistics_[13])/float(N)<<"% |"<<endl;
  cout<<"| # of events with #mu' > 1                    "<<setw(7)<<statistics_[15]<<" "<<setw(6)<<100*float(statistics_[15])/float(N)<<"% |"<<endl;
  cout<<"| # of events with #prim vert > 1              "<<setw(7)<<statistics_[16]<<" "<<setw(6)<<100*float(statistics_[16])/float(N)<<"% |"<<endl;
  cout<<"|--------------------------------------------------------------|"<<endl;
  cout<<"| Average # of tracks in primary vertex        "<<setw(7)<<float(statistics_[8])/float(statistics_[1])<<"         |"<<endl;
  cout<<"| # of events where primary vertex refit       "<<setw(7)<<statistics_[9] <<" "<<setw(6)<<100*float(statistics_[9]) /float(N)<<"% |"<<endl;
  cout<<setprecision(1);
  cout<<"| Number of secondary vertices                 "<<setw(7)<<statistics_[7] <<" "<<setw(6)<<100*float(statistics_[7]) /float(N)<<"% |"<<endl;
  cout<<"|--------------------------------------------------------------|"<<endl;
  double TT = float(statistics_[14]);
  cout<<"| Total time spent per event                   "<<setw(7)<<setprecision(3)<<float(statistics_[14])/1000.<<" "<<setw(6)<<setprecision(0)<<100*float(statistics_[14])/TT<<"% |"<<endl;
  cout<<"| time for Primary Vertices                    "<<setw(7)<<setprecision(3)<<float(statistics_[17])/1000.<<" "<<setw(6)<<setprecision(0)<<100*float(statistics_[17])/TT<<"% |"<<endl;
  cout<<"| time for Secondary Vertices                  "<<setw(7)<<setprecision(3)<<float(statistics_[18])/1000.<<" "<<setw(6)<<setprecision(0)<<100*float(statistics_[18])/TT<<"% |"<<endl;
  cout<<"\\----------\\______________________________________/------------/"<<endl;
  cout<<endl;
  
  CsHistograms::SetCurrentPath("/CsAverPattern/");
  CsHist1S *hStat = new CsHist1S( "hStat","Statistics", 20,0,20 );
  CsHistograms::SetCurrentPath("/");
  
  hStat->SetBinContent(0,N);
  for(int i=1;i<20;i++) hStat->SetBinContent(i,statistics_[i-1]);
  
  return true;
}
