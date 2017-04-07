// $Id: CsPixelGEMDetector.cc 14069 2015-09-17 20:44:46Z lsilva $

/*!
   \file    CsPixelGEMDetector.cc
   \brief   Compass detector Class for pixel GEM.
   \author  Yann.Bedfer@cern.ch
   \version $Revision: 14069 $
   \date    $Date: 2015-09-17 22:44:46 +0200 (Thu, 17 Sep 2015) $
*/

#include "CsGEMDetector.h"
#include "CsPixelGEMDetector.h"

//-----------------------------------------------------------------------------
#include <math.h>
#include <string.h>
#include <strings.h>
#include <functional>
#include "CsZone.h"
#include "CsDigit.h"
#include "CsOpt.h"
#include "CsEvent.h"
#include "CsDigit.h"
#include "CsMCDigit.h"
#include "CsComgNtCommons.h"
#include "CsGeom.h"
#include "CsInit.h"
#include <cstdlib>
#include "CsMCTrkHit.h"
#include "CsRandom.h"
#include "CsGeant3.h"
#include "CsMCTrack.h"
#include "CDB.h"

#include "DaqDataDecoding/ChipAPV.h"

#include "CsGEMPlane.h"
#include "CsPixelGEMPlane.h"

using namespace std;
using namespace CLHEP;
using CS::DetID;

extern QhitType Qhit;

//=============================================================================

CsPixelGEMDetector::CsPixelGEMDetector (const int    row,
					const int    id,    const char* name, const char *TBname,
					const int    unit,  const int    type,
					const double rdLen, const double xsiz,  
					const double ysiz,  const double zsiz,
					const double xcm,   const double ycm,   
					const double zcm,   const HepMatrix rotDRS,
					const HepMatrix rotWRS,
					const double wirD,  const double ang,   
					const int    nWir,  const double wirP, 
					const double eff,   const double bkg,
					const double tGate, const double spSig,
					const double eGain, const double eGSig, 
					const double sWidth,const double tRes ) :
  CsDetector( row, id, name, TBname, unit, type, rdLen, xsiz, ysiz, zsiz,
	      xcm, ycm, zcm, rotDRS, rotWRS, wirD, ang,   
	      nWir, wirP, eff, bkg, tGate ),
  isMaster_(false), associateDet_(NULL) ,
  spSig_(spSig),
  eGain_(eGain),
  eGSig_(eGSig),
  sWidth_(sWidth),
  tRes_(tRes),
  stripplane(NULL), pixelplane(NULL) ,
  apv_chan_cals_mapped(false)
{
  //    ***** # of CHANNELS IN 2ND DIMENSION: *****
  // =0 means CsPG is of the non-pixelised kind. It's the default. It's to be
  // reset by "CsPG::AddSubDetector" for a pixelised CsPG.
  nWirV_ = 0;
  // Time resolution:
  // - Is here built-in = 12 ns, in keeping w/ the value actually retained for
  //  plain GEM's, although "tRes_" data member is not explicitly reset there.
  //  (Note: Have had to declare "CsPixelGEM" a friend class of "CsGEM" in order
  //  to enforce this setting.)
  // - This is not very statisfying: would be better to set "tRes_" from the
  //  detector table. This would pose 2 problems:
  //  i) Redefinition of the field (of the "det" entry in "detectors.dat") used
  //    by "CsGeom::readDetTable" to set the "tRes" argument of the constructor:
  //    it's presently "tslice" (obviously not a t resolution) and is in fact
  //    effectively used by the sole drift detector classes.
  // ii) Modification of all "detectors.dat's".
  tRes_ = 12;

  // Default 
  this->do_clustering=1; // 0: primitive clustering, 1: full clustering
                         //may be superseded by option

  //
  // Options 
  // 
  string TB = GetTBName().substr(0,2);  // First 2 letter of TB name

  // Which detectors to decode
  decodeCard_ = false;
  string tag = ""; 
  string key = "make decoding";
  bool status = CsOpt::Instance()->getOpt( tag, key );
  if( status ) {
    list<string> options;
    list<string>::iterator Is;
    status = CsOpt::Instance()->getOpt( tag, key, options );
    if( status ) {
      for( Is=options.begin(); Is!=options.end(); Is++ ) {
	if( *Is == TB || *Is == "PixelGEM" || *Is == "all" ) {
	  decodeCard_ = true;
	}
      }
    }
    else {
      decodeCard_ = true;
    }
  }

  key = "decode latch all";
  decodeLatch_ = false;
  if (CsOpt::Instance()->getOpt(TBname,key) ||
      CsOpt::Instance()->getOpt(TB,key)) {
    decodeLatch_ = true;
  }
  if (decodeLatch_)
    CsErrLog::msg(elWarning, __FILE__, __LINE__,
            "%s: decoding of latch-all data requested!",TBname);

  // Thresholds for single strip and cluster amplitudes
  // This cut can be made dependent upon detector's TB name, e.g:
  // GP Threshold [0-1] 3. 5.
  // GP01V1__ Threshold [0-1] 2. 5.
  key = "Threshold";
  fThresholdHit_ = 3.;                 // Default is 3 sigma for single hit
  fThresholdClu_ = 5.;                 // Default is 5 sigma for cluster
  vector<double> v;
  if (CsOpt::Instance()->getOpt(TBname,key,v) ||
      CsOpt::Instance()->getOpt(TB,key,v)) {
    if (v.size()!=2) {
      CsErrLog::msg(elFatal,__FILE__,__LINE__,
		    "%s: Syntax error in specifying thresholds!",TBname);
    }
    else {
      fThresholdHit_ = v[0]; fThresholdClu_ = v[1];
    }
  }
  CsErrLog::msg(elInfo,__FILE__,__LINE__,
		"%s: Strip amplitude >=  %.3f sigma, Cluster amplitude >= %.3f sigma",
		TBname,fThresholdHit_,fThresholdClu_);
  
  // Cut on amplitude ratios, i.e. ratios of 1st/3rd and 2nd/3rd amplitudes
  // This cut can be made dependent upon detector's TB name, e.g:
  // GP AmplitudeRatio       [0-7] 0. 0. 1.2 0. 1.2 1. 0. 1.
  // GP01V1__ AmplitudeRatio [0-7] 0. 0. 1.2 0. 1.2 1. 0. 1.
  key = "AmplitudeRatio";
  if (CsOpt::Instance()->getOpt(TBname,key,v) ||
      CsOpt::Instance()->getOpt(TB,key,v)) {
    if (v.size()%2!=0) {
      CsErrLog::msg(elFatal,__FILE__,__LINE__,
		    "%s: Syntax error in specifying cut on amplitude ratio!",
		    TBname); 
    }
    else {
      //      cout << "Amplitude ratio cut:" << endl;
      int np = v.size()/2;
      string mes = "%s: Polygon for amplitude ratio cut\n";
      fAmpRatio23_.clear();
      fAmpRatio13_.clear();
      for (int i=0; i<np; i++){
	fAmpRatio23_.push_back(v[2*i]);   // x coordinates
	fAmpRatio13_.push_back(v[2*i+1]); // y coordinates
	//	cout << v[2*i] << " " << v[2*i+1] << endl;
	char coord[100];
	sprintf(coord,"(%.3f,%.3f)\n",v[2*i],v[2*i+1]);
	mes += coord;
      }
      // Close polygon
      fAmpRatio23_.push_back(v[0]);
      fAmpRatio13_.push_back(v[1]);
      //      cout << v[0] << " " << v[1] << endl;
      CsErrLog::msg(elInfo,__FILE__,__LINE__,
		    mes.c_str(), TBname);
    }
  }

  // Cut on cluster multiplicity 
  // This cut can be made dependent upon detector's TB name, e.g:
  // GP Multiplicity       [0-1] 0 1000
  // GP01V1__ Multiplicity [0-1] 0 1
  key = "Multiplicity";
  fUppMult_ = 1000.;
  fLowMult_ = 0.;
  if (CsOpt::Instance()->getOpt(TBname,key,v) ||
      CsOpt::Instance()->getOpt(TB,key,v)) {
    if (v.size()!=2) {
      CsErrLog::msg(elFatal,__FILE__,__LINE__,
		    "%s: Syntax error in specifying multiplicity cut!",TBname);
    }
    else {
      fLowMult_ = v[0];
      fUppMult_ = v[1];
    }
  }
  CsErrLog::msg(elInfo,__FILE__,__LINE__,
		"%s: %.3f <= Multiplicity <= %.3f",
		TBname,fLowMult_,fUppMult_);

  // Switch for clustering procedure: 
  // "0" for primitive clustering (1 hit = 1 cluster)
  // "1" for full clustering (default, if not specified)
  // This switch can be made dependent upon detector's TB name, e.g:
  // GP Clustering       1
  // GP01V1__ Clustering 1 
  key = "Clustering";
  int valuei;
  do_clustering = 1;
  if (CsOpt::Instance()->getOpt(TBname,key,valuei) ||
      CsOpt::Instance()->getOpt(TB,key,valuei)) {
    do_clustering = valuei;
  }
  CsErrLog::msg(elInfo,__FILE__,__LINE__,
		"%s: Clustering %d",
		TBname,do_clustering);

  // Config mask for PixelGEM clustering procedure: 
  // bit 0 (lsb): enable hit division between adjacent clusters
  // bit 1      : enable diagonal clustering
  // bit 2      : enable time use of time information for clustering
  // bit 3      : enable size splitting
  // bit 4      : enable/disable cluster position size splitting
  // bit 5      : switch between r (1) condition and x&y (0) condition in pos-clust
  // bit 6+7    : enable/disable hit sharing (0 : no sharing; 1 : to all; 2 : fraction)
  // This switch can be made dependent upon detector's TB name, e.g:
  // GP ClusterConfig        0xb3
  // GP01V1__ ClusterConfig  0xb3 
  key = "ClusterConfigMask";
  if (CsOpt::Instance()->getOpt(TBname,key,valuei) ||
      CsOpt::Instance()->getOpt(TB,key,valuei)) {
    fClusterConfigMask = valuei;
    CsErrLog::msg(elInfo,__FILE__,__LINE__,
		  "%s: ClusterConfigMask %d",
		  TBname,fClusterConfigMask);
  }
  else {
    fClusterConfigMask = 0xb3;
    CsErrLog::msg(elInfo,__FILE__,__LINE__,
		  "%s: No ClusterConfigMask specified: Using standard mask %d",
		  TBname, fClusterConfigMask);
  }
  // ***** MONTE CARLO: CUT ON Hit Time *****

  fMCHitTMin = -tGate_/2; fMCHitTMax = tGate_/2;
  // This cut can be made dependent upon detector's TB name, cf. supra:
  key = "MCHitTime";
  if (CsOpt::Instance()->getOpt(TBname,key,v) ||
      CsOpt::Instance()->getOpt(TB,key,v)) {
    if (v.size()!=2)
      CsErrLog::msg(elFatal,__FILE__,__LINE__,
		    "%s: Syntax error in specifying MC hit time cut!",TBname);
    fMCHitTMin = v[0]; fMCHitTMax = v[1];
    CsErrLog::msg(elInfo,__FILE__,__LINE__,"%s: %f <= MC hit time <= %f",
		  TBname,fMCHitTMin,fMCHitTMax);
    if (fMCHitTMin<-tGate_/2. || fMCHitTMax>tGate_/2.)
      CsErrLog::msg(elWarning,__FILE__,__LINE__,
		    "%s: MC hit time window [%f,%f] > tGate %f",
		    TBname,fMCHitTMin,fMCHitTMax,tGate_);
  }

  // Debug flag for CsPixelGEMDetector clusterization, etc.
  // This flag can be made dependent upon detector's TB name, as above
  fDebugLevel = 0; // default: no debugging
  key = "debug level";
  if (CsOpt::Instance()->getOpt(TBname,key,fDebugLevel) ||
      CsOpt::Instance()->getOpt(TB,key,fDebugLevel)) {
  }

  //=== Check if amplitudes simulation in MC decoding should be used ===

  amplitudeMCDecoding_ = false; 
  if (CsOpt::Instance()->getOpt(TBname,"amplitudeMCDecoding") ||
      CsOpt::Instance()->getOpt(TB,"amplitudeMCDecoding")) {
    amplitudeMCDecoding_ = true;
    CsErrLog::msg(elWarning,__FILE__,__LINE__,
		  "%s: amplitudes simulation in MC decoding will be used",TBname);
  }

  // ***** MONTE CARLO: Amlitudes simulation parameters *****
  //
  // spSig_ detector space resolution (mm) 
  // eGain_ effective gain - should be tuned to reproduce cluster amlitudes   
  // eGSig_ gain sigma (a.u.) for amplitude correlation, for example.
  // sWidth_ signal width (mm) (should be tuned to have correct number of strips/cluster)
  // tRes_ time resolution (ns)
  // These pars can be made dependent upon detector's TB name, cf. supra:
  // GP ampParsMC [0-4] 0.05  2000. 20.  0.30  12.
  // GP01X1__ ampParsMC [0-4] 0.05  2000. 20. 0.30  12.
  key = "ampParsMC";
  if (CsOpt::Instance()->getOpt(TBname,key,v) ||
      CsOpt::Instance()->getOpt(TB,key,v)) {

    if (v.size()!=5)
      CsErrLog::msg(elFatal,__FILE__,__LINE__,
		    "%s: Syntax error in specifying MC amplitude parameters!",TBname);
    spSig_ = v[0]; eGain_ = v[1]; eGSig_ = v[2]; sWidth_ = v[3]; tRes_ = v[4];
    CsErrLog::msg(elInfo,__FILE__,__LINE__,"%s: space res. %f mm, gain %f, "
                  "sigma gain %f , signal width %f mm, time res. %f ns", 
                  TBname,spSig_,eGain_,eGSig_,sWidth_,tRes_);
  }

  if(amplitudeMCDecoding_) {
    if( sWidth_ <= 0. || eGain_ <=0. ){
     CsErrLog::msg(elFatal,__FILE__,__LINE__,
                    "%s: amplitude simulation impossible, wrong " 
                    "width=%f or eGain=%f  ",
                     TBname,sWidth_,eGain_);
    } 
  }

  //=== Check if amplitudes correlation in MC should be used ===

  string vs;
  ampCorrelationMC_ = false; 
  isMaster_ = false;
  if (CsOpt::Instance()->getOpt(TBname,"ampCorrelationMC") ||
      CsOpt::Instance()->getOpt(TB,"ampCorrelationMC")) {
    ampCorrelationMC_ = true;
    CsErrLog::msg(elWarning,__FILE__,__LINE__,
		  "%s: amplitudes correlation in MC decoding will be used",TBname);
    key = "Master";
    if ( CsOpt::Instance()->getOpt(TB,key,vs) ) {
      if (vs.size()==2) {
	if(GetTBName()[4]==vs[0] || GetTBName()[4]==vs[1] ) {
	  isMaster_ = true;
	}  
      }
    }  
    key = "Master";
    if(CsOpt::Instance()->getOpt(TBname,key)){
      isMaster_ = true;    
    }
    key = "Slave";
    if ( CsOpt::Instance()->getOpt(TBname,key)) {
      isMaster_ = false;
    }
    CsErrLog::msg(elInfo,__FILE__,__LINE__,
		  "%s: detector is master if non-zero: %i ",TBname,int(isMaster_) );
  }

  fAmpCorr[0] = 0; fAmpCorr[1] = 1; fAmpCorr[2] = 0;
  if (CsOpt::Instance()->getOpt(TBname,"AmplitudeCorr",v)) {
    if (v.size()!=3)
      CsErrLog::msg(elFatal,__FILE__,__LINE__,
		    "%s: Syntax error in specifying AmplitudeCorr!",TBname);
    // ***** Correction to the amplitude in order to match that of counterpart
    fAmpCorr[0] = v[0]; fAmpCorr[1] = v[1]; fAmpCorr[2] = v[2];
    CsErrLog::msg(elInfo,__FILE__,__LINE__,"%s: Amplitude Corr %f %f %f",
		  TBname,fAmpCorr[0],fAmpCorr[1],fAmpCorr[2]);
  }
  
  // read parameters for cross talk suppression algorithm
  // standard parameters according to M. Kraemers diploma thesis for the new read out design
  fCrossTalkParams[0] = 0.30; fCrossTalkParams[1] = 0.55;
  if (CsOpt::Instance()->getOpt(TBname,"CrossTalkParams",v)  ||
      CsOpt::Instance()->getOpt(TB,"CrossTalkParams",v)) {
      if (v.size()!=2)
          CsErrLog::msg(elFatal,__FILE__,__LINE__,
                        "%s: Syntax error in specifying CrossTalkParams!",TBname);
      fCrossTalkParams[0] = v[0]; fCrossTalkParams[1] = v[1];
      CsErrLog::msg(elInfo,__FILE__,__LINE__,"%s: CrossTalkParams %f %f",
                    TBname,fCrossTalkParams[0],fCrossTalkParams[1]);
  }

  double valued;
  fTimeCrossTalkParam = .0;
  if (CsOpt::Instance()->getOpt(TBname,"TimeCrossTalkParam",valued)  ||
      CsOpt::Instance()->getOpt(TB,"TimeCrossTalkParam",valued)) {
          fTimeCrossTalkParam = valued;
          CsErrLog::msg(elInfo,__FILE__,__LINE__,"%s: TimeCrossTalkParams %f",
                        TBname,fTimeCrossTalkParam);
  }

  // read parameters for position correction of pixel clusters
  std::vector<float> fPosCorrs;
  for (unsigned int i=0; i<48; i++) fPosCorrs.push_back(0.);
  if (CsOpt::Instance()->getOpt(TBname,"PosCorrs",v)  ||
      CsOpt::Instance()->getOpt(TB,"PosCorrs",v)) {
    fPosCorrs.clear();
    for (unsigned int i=0; i<v.size(); i++) fPosCorrs.push_back(v[i]);
  }

  // Create CsGEMPlane (to be used in CsPG::clusterize)
  if (TBname[4]!='P') {
    ChipsPerPlane = ChipsPerStripPlane;

    stripplane = new CsGEMPlane(GetTBName().c_str());
    stripplane->GetPar()->SetClusMeth(do_clustering);
    stripplane->GetPar()->SetThrClus(fThresholdClu_);
    stripplane->GetPar()->SetThrHit(fThresholdHit_);
    stripplane->GetPar()->SetTimeCrossTalkHitRatio(fTimeCrossTalkParam);
    stripplane->GetPar()->SetSample(2); // sample used for clustering

    // read cluster-size specific spatial resolution
    std::vector<float> clusSizeRes;
    if (CsOpt::Instance()->getOpt(TBname,"StripClusSizeRes", clusSizeRes) ||
        CsOpt::Instance()->getOpt(TB,"StripClusSizeRes", clusSizeRes)) {
      std::ostringstream msg;
      msg << TBname << ": cluster-size specific spatial resolution for strip plane:";
      if (clusSizeRes.size()>0) {
        for (size_t i=1; i<clusSizeRes.size(); i++)
          msg << " " << i << "->" << clusSizeRes[i-1];
        msg << " >=" << clusSizeRes.size() << "->" << clusSizeRes[clusSizeRes.size()-1];
        CsErrLog::msg(elInfo, __FILE__, __LINE__, msg.str().c_str());
        stripplane->GetPar()->SetClusSizeRes(clusSizeRes);
      }
    }
  }
  else {
    ChipsPerPlane = ChipsPerPixelPlane;

    pixelplane = new CsPixelGEMPlane(GetTBName().c_str()); // pixel part
    pixelplane->GetPar()->SetClusConfigMask(fClusterConfigMask);
    pixelplane->GetPar()->SetClusMeth(do_clustering);
    pixelplane->GetPar()->SetThrClus(fThresholdClu_);
    pixelplane->GetPar()->SetThrHit(fThresholdHit_);
    pixelplane->GetPar()->SetCrossTalkHitRatio(fCrossTalkParams[0]);
    pixelplane->GetPar()->SetTimeCrossTalkHitRatio(fTimeCrossTalkParam);
    pixelplane->GetPar()->SetSample(2); // sample used for clustering

    // read cluster-size specific spatial resolution
    std::vector<float> clusSizeRes;
    if (CsOpt::Instance()->getOpt(TBname,"PixelClusSizeRes", clusSizeRes) ||
        CsOpt::Instance()->getOpt(TB,"PixelClusSizeRes", clusSizeRes)) {
      std::ostringstream msg;
      msg << TBname << ": cluster-size specific spatial resolution for pixel plane:";
      if (clusSizeRes.size()>0) {
        for (size_t i=1; i<clusSizeRes.size(); i++)
          msg << " " << i << "->" << clusSizeRes[i-1];
        msg << " >=" << clusSizeRes.size() << "->" << clusSizeRes[clusSizeRes.size()-1];
        CsErrLog::msg(elInfo, __FILE__, __LINE__, msg.str().c_str());
        pixelplane->GetPar()->SetClusSizeRes(clusSizeRes);
      }
    }

    if (!pixelplane->GetPar()->SetPosCorrs(fPosCorrs))
      CsErrLog::msg(elError,__FILE__,__LINE__,"%s: Wrong format of position corrections", TBname);
  }
}

CsPixelGEMDetector::~CsPixelGEMDetector() {
  delete stripplane;
  delete pixelplane;
}

void CsPixelGEMDetector::AddSubDetector(const int    row,
					const int    id,    const char* name, const char *TBname,
					const int    unit,  const int    type,
					const double rdLen, const double xsiz,  
					const double ysiz,  const double zsiz,
					const double xcm,   const double ycm,
					const double zcm,   const HepMatrix rotDRS,
					const HepMatrix rotWRS,
					const double wirD,  const double ang,
					const int    nWir,  const double wirP, 
					const double eff,   const double bkg,
					const double tGate)
{
  // ***** 2ND PART OF THE CONSTRUCTION OF A CsPG OF THE PIXELISED KIND *****

  CsErrLog::msg(elInfo,__FILE__,__LINE__,
		"Adding subdetector to %s",GetTBName().c_str());

  if (fDebugLevel>9) {
    //    cout << "CsPixelGEMDetector::AddSubDetector : Adding subdetector to "
    //         << GetTBName().c_str() << endl;
  }

  // The following parameters must be the same in all subdetectors.

  if(strcmp(TBname,GetTBName().c_str()))
    CsErrLog::msg(elFatal,__FILE__,__LINE__,"AddSubDetector: %s +=%s: Dissimilar TBnames",GetTBName().c_str(),TBname);
  if(type!=getType())
    CsErrLog::msg(elFatal,__FILE__,__LINE__,"AddSubDetector: %s(%d) +=%s(%d): Dissimilar types",GetTBName().c_str(),getType(),TBname,type);
  if(rdLen-rdLen_)
    CsErrLog::msg(elFatal,__FILE__,__LINE__,"AddSubDetector: %s, #%d(%f) +=#%d(%f): Dissimilar rad.lengths",GetTBName().c_str(),GetID().GetNumber(),rdLen_,id,rdLen);
  if (!(rotDRS==rotDRS_))
    CsErrLog::msg(elFatal,__FILE__,__LINE__,"AddSubDetector: %s, #%d +=#%d: Added piece's MRS->DRS != this CsPixelGEM's",
		  GetTBName().c_str(),GetID().GetNumber(),id);
  HepMatrix x2y(3,3,0); x2y[1][0] = -1; x2y[0][1] = 1; x2y[2][2] = 1;
  HepMatrix M = rotWRS*x2y; M -= rotWRS_;
  double diff; int i, j; for (i = 0, diff = 0; i<3; i++) for (j = 0; j<3; j++)
    if (fabs(M[i][j])>diff) diff = fabs(M[i][j]);
  if (diff>1e-6)
    CsErrLog::msg(elFatal,__FILE__,__LINE__,"AddSubDetector: %s, #%d +=#%d: Added piece's MRS->WRS != this CsPixelGEM's + pi/2",
		  GetTBName().c_str(),GetID().GetNumber(),id);

  // Cancel amplitude MC decoding: not yet implemented for the pixelised (as
  // opposed to stripped) piece
  amplitudeMCDecoding_ = false;

  // Setting parameters specific to 2nd dimension 

  wirDV_ = wirD;		// Wires distance    in 2nd dimension
  angV_  = ang;			// Wires angle       in 2nd dimension
  sinAngV_ = sin(angV_*M_PI/180); cosAngV_ = cos(angV_*M_PI/180);
  nWirV_ = nWir;		// Number of wires   in 2nd dimension
  wirPV_ = wirP;		// Wire pitch        in 2nd dimension
}

////////////////////////////////////////////////////////////////////////////////
void CsPixelGEMDetector::BookHistograms( void )
{
  
  // Check if histograms are to be booked
  CsDetector::ReadHistLevel();
  if( hLevel_ == None ) return;
  
  string tbn  = GetTBName(); 
  CsHistograms::SetCurrentPath("/PixelGEMs");
  
  // Normal level
  if( hLevel_ >= Normal ) {
    mH1[tbn+"_cA1T"]=new CsHist1D(tbn+"_cA1T",tbn+" cluster amplitude 1 time cut",256,-2.,1022.);
    mH1[tbn+"_cA2T"]=new CsHist1D(tbn+"_cA2T",tbn+" cluster amplitude 2 time cut",256,-2.,1022.);
    mH1[tbn+"_cA3T"]=new CsHist1D(tbn+"_cA3T",tbn+" cluster amplitude 3 time cut",256,-2.,1022.);
    mH1[tbn+"_cSizeT"]=new CsHist1D(tbn+"_cSizeT",tbn+" cluster size time cut",32,-0.5,31.5);
    mH1[tbn+"_nClusT"]=new CsHist1D(tbn+"_nClusT",tbn+" number of clusters time cut",32,-0.5,31.5);
    mH1[tbn+"_cTimeT"]=new CsHist1D(tbn+"_cTimeT",tbn+" cluster time time cut",100,-45.5,50.5);
    mH2[tbn+"_cAmpRatioT"]=new CsHist2D(tbn+"_cAmpRatioT",tbn+" cluster amplitude ratio time cut",25,-0.1,2.4,25,-0.1,2.4);
    mH1[tbn+"_cCTProbC"]=new CsHist1D(tbn+"_cCTProbC",tbn+" cluster cross talk probability after cut",100,0.,1.001);

    // Pixel region, designated by a "P" in the projection field of TBname
    if (tbn.find("P", 1, 1) != string::npos) { 
      mH2[tbn+"_cMapT"]=new CsHist2D(tbn+"_cMapT",tbn+" cluster map time cut",nWirU_,-0.5,nWirU_-0.5,nWirV_,-0.5,nWirV_-0.5);
    } else { // strip region
      mH1[tbn+"_cPosT"]=new CsHist1D(tbn+"_cPosT",tbn+" cluster position time cut",2*nWir_+1,-0.5-nWir_,nWir_+0.5);
      mH1[tbn+"_hPosT"]=new CsHist1D(tbn+"_hPosT",tbn+" hit position time cut",2*nWir_+1,-0.5-nWir_,nWir_+0.5);
    }
    
  }
  
  // High level
  if( hLevel_ >= High ) {
    mH1[tbn+"_cA1"]=new CsHist1D(tbn+"_cA1",tbn+" cluster amplitude 1 all",256,-2.,1022.);
    mH1[tbn+"_cA2"]=new CsHist1D(tbn+"_cA2",tbn+" cluster amplitude 2 all",256,-2.,1022.);
    mH1[tbn+"_cA3"]=new CsHist1D(tbn+"_cA3",tbn+" cluster amplitude 3 all",256,-2.,1022.);
    mH1[tbn+"_nHit"]=new CsHist1D(tbn+"_nHit",tbn+" number of hit strips all",128,-0.5,127.5);
    mH1[tbn+"_nClus"]=new CsHist1D(tbn+"_nClus",tbn+" number of clusters all",32,-0.5,31.5);
    mH2[tbn+"_cAmpRatio"]=new CsHist2D(tbn+"_cAmpRatio",tbn+" cluster amplitude ratio all",25,-0.1,2.4,25,-0.1,2.4);
    mH1[tbn+"_cCTProbT"]=new CsHist1D(tbn+"_cCTProbT",tbn+" cluster cross talk probability time cut",100,0.,1.001);

    // Pixel region, designated by a "P" in the projection field of TBname
    if (tbn[4] == 'P') { 
      mH2[tbn+"_hMap"]=new CsHist2D(tbn+"_hMap",tbn+" hit map all",nWirU_,-0.5,nWirU_-0.5,nWirV_,-0.5,nWirV_-0.5);
      mH2[tbn+"_cMap"]=new CsHist2D(tbn+"_cMap",tbn+" cluster map all",nWirU_,-0.5,nWirU_-0.5,nWirV_,-0.5,nWirV_-0.5);
    } else { // strip region
      mH1[tbn+"_cPos"]=new CsHist1D(tbn+"_cPos",tbn+" cluster position all",2*nWir_+1,-0.5-nWir_,nWir_+0.5);
      mH1[tbn+"_hPos"]=new CsHist1D(tbn+"_hPos",tbn+" hit position all",2*nWir_+1,-0.5-nWir_,nWir_+0.5);
    }
    
    if( CsInit::Instance()->IsAMonteCarloJob() )  {  
      mH1[tbn+"_MCdA1"]=new CsHist1D(tbn+"_MCdA1",tbn+" MC digit amplitude 1 all",256,-2.,1022.);
      mH1[tbn+"_MCdA2"]=new CsHist1D(tbn+"_MCdA2",tbn+" MC digit amplitude 2 all",256,-2.,1022.);
      mH1[tbn+"_MCdA3"]=new CsHist1D(tbn+"_MCdA3",tbn+" MC digit amplitude 3 all",256,-2.,1022.);
      mH1[tbn+"_MCdnHit"] = new CsHist1D(tbn+"_MCdnHit",tbn+" number of MC hits per digit",32,-0.5,31.5);
      mH2[tbn+"_MCdAmpRatio"]=new CsHist2D(tbn+"_MCdAmpRatio",tbn+" MC digit amplitude ratio all",25,-0.1,2.4,25,-0.1,2.4); 
      mH2[tbn+"_MCnDvsT"]=new CsHist2D(tbn+"_MCnDvsT",tbn+" number of MC digits per hit vs time",
				       40,-220,100,32,-0.5,31.5);       
    }
  }  

  CsHistograms::SetCurrentPath("/");
}   

///////////////////////////////////////////////////////////////////////////////
class RawInfo {
public:
  RawInfo() : 
    fDigit(0) {
  }
  RawInfo(const CS::ChipAPV::Digit* digit) : 
    fDigit(digit) {
    for (size_t i=0; i<3; i++)
      fAmps[i] = 1023-digit->GetAmplitude()[i];
  }

  const CS::ChipAPV::Digit* getDigit()  const {
    return fDigit;
  }

  const int*                getAmps() const {
    return fAmps;
  }

  void                      subtract(int val) {
    for (size_t i=0; i<3; i++) {
      if (fAmps[i] > val)
        fAmps[i] -= val;
      else
        fAmps[i] = 0;
    }
  }

  void                      subtract(int* val) {
    for (size_t i=0; i<3; i++)
      if (fAmps[i] > val[i])
        fAmps[i] -= val[i];
      else
        fAmps[i] = 0;
  }

private:
  const CS::ChipAPV::Digit* fDigit;

  int                       fAmps[3];
};

void CsPixelGEMDetector::DecodeChipDigits(const CS::Chip::Digits &digits )
{
  if (!decodeLatch_)
    return CsDet::DecodeChipDigits(digits);

  if (fDebugLevel>9)
    std::cout << "CsPixelGEMDetector::DecodeChipDigits: " <<  GetTBName() << std::endl;
  
  // Apply mapping to channel calibrations if not done yet
  // cannot be done in CsGEMDetector::readCalibration because 
  // CS::Chip::Maps are not yet available
  if (!apv_chan_cals_mapped) apv_chan_cals_mapped = mapChannelCal();

  typedef std::multimap<CS::DetID, CS::Chip::Digit*>::const_iterator m_it; // iterator type
  // get all digits for the detector
  std::pair<m_it, m_it> m_range = digits.equal_range(GetTBName());

  std::map<std::pair<uint16, uint16>, std::vector<RawInfo> > data;

  // loop on all found digits
  for (m_it d_it=m_range.first; d_it!=m_range.second; d_it++) {
    // Check digit
    const CS::ChipAPV::Digit* d = dynamic_cast<const CS::ChipAPV::Digit*>(d_it->second);

    if (d==NULL)
      throw CS::Exception("CsPixelGEMDetector::DecodeChipDigits(): Wrong digit type.");

    if (d->IsSparsifed())
      CsErrLog::msg(elFatal, __FILE__, __LINE__,
              "%s: Cannot decode already sparsified data, when decoding latch-all data is requested.", GetTBName().c_str());

    const CS::ChipAPV::DataID& dataId = reinterpret_cast<const CS::ChipAPV::DataID&>(d->GetDataID());

    // in mapChannelCal it has been checked that this detector is only connected to one source ID,
    // each APV chip can be uniquely identified via the combination of ADC ID and chip ID
    std::pair<uint16, uint16> id(dataId.u.s.adc_id, d->GetChip());

    std::map<std::pair<uint16, uint16>, std::vector<RawInfo> >::iterator data_it = data.find(id);
    if (data_it==data.end()) {
      std::vector<RawInfo> content;
      content.resize(128);
      data.insert(std::pair<std::pair<uint16, uint16>, std::vector<RawInfo> >(id, content));
      data_it = data.find(id);
    }

    assert(data_it!=data.end());
    assert(d->GetChipChannel()<128);
    // do not check this for the moment and skip the second digit
    // this currently happens with the CMC word in latch-all mode
    //assert(data_it->second[d->GetChipChannel()].getDigit()==NULL);
    if (data_it->second[d->GetChipChannel()].getDigit()==NULL) {
      data_it->second[d->GetChipChannel()] = RawInfo(d);
    }
  }

  for (std::map<std::pair<uint16, uint16>, std::vector<RawInfo> >::iterator d_it=data.begin(); d_it!=data.end(); d_it++) {
    std::vector<RawInfo>& channels = d_it->second;
    assert(channels.size()==128);

    const uint16 adcId = d_it->first.first;
    const uint16 apvId = d_it->first.second;

    std::vector<APVCal>::const_iterator it=apv_chan_cals.begin();
    for (; it!=apv_chan_cals.end(); it++)
      if (it->adc_id==adcId && it->chip_id==apvId)
        break;

    if (it==apv_chan_cals.end())
      CsErrLog::msg(elFatal, __FILE__, __LINE__,
              "%s: Could not find pedestals for adcId=%d, chipId=%d.", GetTBName().c_str(), d_it->first.first, d_it->first.second);

    const std::vector<APVCal::Channel>& cal = it->channels;
    assert(cal.size()==128);

    int smallest(1023);
    for (size_t i=0; i<128; i++)
      if ((1024-static_cast<int>(cal[i].pedestal_mean+0.5)) < smallest)
        smallest = 1024-static_cast<int>(cal[i].pedestal_mean+0.5);

    // 1. do the pedestal subtraction
    // 2. check that this chip has been received completely
    unsigned int count(0);
    int means[3] = {0, 0, 0};
    for (std::vector<RawInfo>::iterator c_it=channels.begin(); c_it!=channels.end(); c_it++) {
      if (c_it->getDigit()==0)
        break;

      c_it->subtract(1024-static_cast<int>(cal[count].pedestal_mean+0.5)-smallest);

      for (size_t i=0; i<3; i++)
        means[i] += c_it->getAmps()[i];

      count++;
    }

    // if not all 128 channels have been present
    if (count!=128)
      CsErrLog::msg(elFatal, __FILE__, __LINE__,
              "%s: Not all 128 channels of adcId=%d, chipId=%d have been received. This should be caught by the decoding.", GetTBName().c_str(), d_it->first.first, d_it->first.second);

    for (size_t i=0; i<3; i++) {
      means[i] /= 128;
      means[i] -= 32;
    }

    int hists[3][32]; for (size_t i=0; i<3; i++) for (size_t j=0; j<32; j++) hists[i][j] = 0;
    for (std::vector<RawInfo>::iterator c_it=channels.begin(); c_it!=channels.end(); c_it++) {
      c_it->subtract(means);

      for (size_t i=0; i<3; i++) {
        for (size_t j=0; j<32; j++) {
          if (2*((int)j) > c_it->getAmps()[i])
            break;
          hists[i][j]++;
        }
      }
    }

    size_t meansIdx[3] = {0, 0, 0};
    int cmc[3] = {0, 0, 0};
    for (size_t i=0; i<3; i++) {
      while (meansIdx[i]<32 && hists[i][meansIdx[i]] >= 64)
        meansIdx[i]++;
      if (meansIdx[i]>0)
        meansIdx[i]--;  
      cmc[i] = 2*meansIdx[i];
      means[i] += cmc[i];
    }

    // subtract the remaining part of the common mode correction
    for (std::vector<RawInfo>::iterator c_it=channels.begin(); c_it!=channels.end(); c_it++)
      c_it->subtract(cmc);

    for (std::vector<RawInfo>::iterator c_it=channels.begin(); c_it!=channels.end(); c_it++) {
      const CS::ChipAPV::Digit* d = c_it->getDigit();

      // Check if pixel digit
      const CS::ChipAPV::DigitPixel* dp = dynamic_cast<const CS::ChipAPV::DigitPixel *>(d);

      CS::ChipAPV::Digit* digit(0);

      // Strip  digit 
      if (dp==NULL)
        digit = new CS::ChipAPV::Digit(*d);
      else
        digit = new CS::ChipAPV::DigitPixel(*dp);
      digit->SetAmplitude(c_it->getAmps()[0], c_it->getAmps()[1], c_it->getAmps()[2]);
      digit->SetCoNo(means[0], means[1], means[2]);
      digit->SetSparsifed(true);

      DecodeChipDigit(*digit);
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
void CsPixelGEMDetector::DecodeChipDigit(const CS::Chip::Digit &digit )
{

  string tbn = GetTBName();
  char mess[500];
  
  if (fDebugLevel>9) {
    std::cout << "CsPixelGEMDetector::DecodeChipDigit : " <<  tbn << std::endl;
  }
  
  // Apply mapping to channel calibrations if not done yet
  // cannot be done in CsGEMDetector::readCalibration because 
  // CS::Chip::Maps are not yet available
  if (!apv_chan_cals_mapped) apv_chan_cals_mapped = mapChannelCal();

  // Check digit
  const CS::ChipAPV::Digit *d = dynamic_cast<const CS::ChipAPV::Digit *>(&digit);
  if( d==NULL )
    throw CS::Exception("CsPixelGEMDetector::DecodeChipDigit(): "
			"Wrong digit type");
  
  // Check if pixel digit
  const CS::ChipAPV::DigitPixel *dp = dynamic_cast<const CS::ChipAPV::DigitPixel *>(d);

  // Strip  digit 
  if( dp==NULL ) {
    
    if (d->GetChannel()>=int(getNWir())) {
      d->Print(cerr, "BAD DIGIT:  ");
      throw CS::Exception("CsPixelGEMDetector::DecodeChipDigit(Strip) "
                          "==> Unexpected strip number %d for PixelGEM %d %s",
			  d->GetChannel(), GetID().GetNumber(), GetTBName().c_str());
    }
    
    if( d->GetChipChannel()>=ChipChannels ){
      d->Print(cerr, "BAD DIGIT:  ");
      throw CS::Exception("CsPixelGEMDetector::DecodeChipDigit(Strip): "
                          "bad channel number! chip=%d, chan=%d(max=%d)",
                          d->GetChip(), d->GetChipChannel(), ChipChannels-1);
    }
    
    if( d->IsSingleFrame() ) {
      d->Print(cerr, "BAD DIGIT:  ");
      throw CS::Exception("CsPixelGEMDetector::DecodeChipDigit(Strip): single frame is not supported yet!");
    }
    
    if( !d->IsSparsifed() ) {
      d->Print(cerr, "BAD DIGIT:  ");
      throw CS::Exception("CsPixelGEMDetector::DecodeChipDigit(Strip): "
                          "found digits without sparsification. "
                          "Latch all events not yet supported!");
    }
    
    int chan_position = d->GetChanPos();
    if (chan_position>1 || chan_position<-1)
      throw CS::Exception("CsPixelGEMDetector::DecodeChipDigit(Strip) "
                          "==> Unexpected channel position %d for PixelGEM %d %s",
			  d->GetChanPos(), GetID().GetNumber(), GetTBName().c_str());
    
    int iwir = d->GetChannel();
    
    // Get data
    vector<double> data;
    data.push_back(d->GetAmplitude()[0]);
    data.push_back(d->GetAmplitude()[1]);
    data.push_back(d->GetAmplitude()[2]);
    data.push_back(d->GetChanPos());
    data.push_back(d->GetChipChannel());
    
    // Fill histograms 
    if (hLevel_>=High) {   
      if(mH1[tbn+"_hPos"]!=NULL)  mH1[tbn+"_hPos"]->Fill(iwir); // Hit position
    }
    
    //
    // Create CORAL digit
    //
    myDigits_.push_back( new CsDigit(*this,iwir,&data[0],data.size()) );
    
    // Log
    if (fDebugLevel>9) {
      sprintf(mess,
              "CsPixelGEMDetector::DecodeChipDigit(Strip) : %s \n"
              "Created digit at %d, hem %d, amplitudes %d, %d, %d",
	      tbn.c_str(),
	      iwir,
	      chan_position,
	      d->GetAmplitude()[0], 
	      d->GetAmplitude()[1], 
	      d->GetAmplitude()[2]);
      std::cout << mess << std::endl;
    }
  }

  // Pixel digit
  else {

    if (dp->GetPixelX()>=nWirU_ || dp->GetPixelY()>=nWirV_ ||
	dp->GetPixelX()<0 || dp->GetPixelY()<0) {
      dp->Print(cerr, "BAD DIGIT:  ");
      throw CS::Exception("CsPixelGEMDetector::DecodeChipDigit(Pixel) "
                          "==> Unexpected pixel number %d, %d for PixelGEM %d %s",
			  dp->GetPixelX(), dp->GetPixelY(),
			  GetID().GetNumber(), GetTBName().c_str());
    }
    
    if (dp->GetChannel()>=int(getNWir())*int(getNWirV())) {
      dp->Print(cerr, "BAD DIGIT:  ");
      throw CS::Exception("CsPixelGEMDetector::DecodeChipDigit(Pixel) "
                          "==> Unexpected channel number %d for PixelGEM %d %s",
			  dp->GetChannel(), GetID().GetNumber(), GetTBName().c_str());
    }

    if( dp->GetChipChannel()>=ChipChannels ){
      dp->Print(cerr, "BAD DIGIT:  ");
      throw CS::Exception("CsPixelGEMDetector::DecodeChipDigit(Pixel): "
                          "bad channel number! chip=%d, chan=%d(max=%d)",
                          dp->GetChip(), dp->GetChipChannel(), ChipChannels-1);
    }
    
    if( dp->IsSingleFrame() ) {
      dp->Print(cerr, "BAD DIGIT:  ");
      throw CS::Exception("CsPixelGEMDetector::DecodeChipDigit(Pixel): "
                          "single frame is not supported yet!");
    }
    
    if( !dp->IsSparsifed() ) {
      dp->Print(cerr, "BAD DIGIT:  ");
      throw CS::Exception("CsPixelGEMDetector::DecodeChipDigit(Pixel): "
                          "found digits without sparsification. "
                          "Latch all events not yet supported!");
    }
    
    int iwir = d->GetChannel();
    int pixX = dp->GetPixelX();
    int pixY = dp->GetPixelY();

    // Get data
    vector<double> data;
    data.push_back((double)dp->GetAmplitude()[0]);
    data.push_back((double)dp->GetAmplitude()[1]);
    data.push_back((double)dp->GetAmplitude()[2]);
    data.push_back((double)(pixX+(pixY<<10)));
    data.push_back((double)dp->GetChipChannel());

    // Fill histograms 
    if (hLevel_>=High) {   
      if(mH2[tbn+"_hMap"]!=NULL)  mH2[tbn+"_hMap"]->Fill(pixX,pixY); // Hit map
    }
    
    //
    // Create CORAL digit, encode pixel address
    //
    myDigits_.push_back( new CsDigit(*this,iwir,&data[0],data.size()) );

    // Log
    if (fDebugLevel>9) {
      sprintf(mess,
              "CsPixelGEMDetector::DecodeChipDigit(Pixel) : %s \n"
              "Created digit at %d, %d, %d, amplitudes %d, %d, %d",
	      tbn.c_str(),iwir,pixX,pixY,
	      dp->GetAmplitude()[0], 
	      dp->GetAmplitude()[1], 
	      dp->GetAmplitude()[2]);
      cout << mess << endl;
    }
  }
  return;
}
  
////////////////////////////////////////////////////////////////////////////
void CsPixelGEMDetector::clusterize() {
    // Clear list of CORAL clusters
    clearClusterList();

    if (pixelplane) { // Pixel region
        clusterizePixels();
    } else if (stripplane) { // Strip region
        clusterizeStrips();
    }
}

//-----------------------------------------------------------------------------
// Clusterization for pixel part of PixelGEM detector
//-----------------------------------------------------------------------------
void CsPixelGEMDetector::clusterizePixels() {
  
  string tbn  = GetTBName();     

  // CORAL digits
  list<CsDigit*>::iterator Id;
  list<CsDigit*> digits = getMyDigits();
  vector<list<CsDigit*>::iterator> iterators;
  iterators.clear();
  
  // Protection
  if( digits.empty() ) return;
  
  HepMatrix iRotM = getRotWRSInv();
  double wireDCorrU = wirDU_ + iRotM(1,1) * _deltax + iRotM(1,2) *_deltay; 
  double wireDCorrV = wirDV_ + iRotM(2,1) * _deltax + iRotM(2,2) *_deltay; 
  
  // Map to store CORAL digits
  map<int,CsDigit*> digitmap;
  
  // Loop over digits and fill them as hits into the plane
  for( Id=digits.begin(); Id!=digits.end(); Id++ ) {
    
    // Make sure that digit belongs to this detector
    assert( (*Id)->getDet() == this );

    // Skip if already clusterized
    if( (*Id)->isClusterized() ) continue;
    
    // Get amplitudes
    const float a0 = (float)(*Id)->getData()[0];
    const float a1 = (float)(*Id)->getData()[1];
    const float a2 = (float)(*Id)->getData()[2];
    
    // Detector channel
    const int channel = (int)(*Id)->getAddress();
    
    // Decode pixel address
    int pixU = ((int)(*Id)->getData()[3]) & 0x3ff;
    int pixV = (((int)(*Id)->getData()[3])>>10) & 0x3ff;

    int debug(0);
    
    if ( pixU < 0 || pixU > 31 || pixV < 0 || pixU > 31 ) {
      if ( debug )
	cout << "--------------------------------------------------\n"
	     << "PixelGEM: Hit with pixU = " << pixU << " , pixV = " << pixV << endl
	     << "Channel = " << (int)(*Id)->getAddress() << endl
	     << "Chip channel = "  << (int)(*Id)->getData()[4] << endl
	     << "Hit ignored!\n"
	     << "--------------------------------------------------\n";
      continue;
    }

    // Add hit to plane
    pixelplane->AddHit(channel, pixU, pixV, a0, a1, a2);
    
    // Add reference to CsDigit
    digitmap[channel] = *Id;
    
  } // end of loop over digits
  
  // Fill histogram with number of hits
  if (hLevel_>=High){
    if(mH1[tbn+"_nHit"]!=NULL)  mH1[GetTBName()+"_nHit"]->Fill(pixelplane->GetNhit());
  }
  
  // Clustering
  pixelplane->Clusterize();
  
  // Display and print clusters
  if (fDebugLevel>9) {
    pixelplane->PrintClusters();
    pixelplane->Display();
  }
  
  // Interface to CORAL clusters
  int nclus = pixelplane->GetNcluster();
  int ngoodclus = 0;

  double triggerOffset(0.); // MC: Time offset of the trigger w.r.t. the event
  const CsEvent *event = CsEvent::Instance();
  bool isMC = event->isAMonteCarloEvent();
  if (isMC) triggerOffset = event->getTriggerMCOffset();

  // Get found clusters
  list<CsPixelGEMCluster*> clusters = pixelplane->GetClusters();

  // Loop over clusters to fill CORAL clusters
  list<CsPixelGEMCluster*>::iterator itclus;
  for (itclus=clusters.begin(); itclus!=clusters.end(); itclus++) {
    
    float pixU = (*itclus)->GetPositionX();
    float pixV = (*itclus)->GetPositionY();
    
    // Get cluster amplitudes
    std::vector<Float_t> Amp = (*itclus)->GetAmp();
    float Amp13 = (Amp[2]==0) ? 0 : 
      (Amp[0]/Amp[2]);
    float Amp23 = (Amp[2]==0) ? 0 :
      (Amp[1]/Amp[2]);
    
    // Fill histograms before time cut
    if (hLevel_>=High){
      if(mH1[tbn+"_cA1"]!=NULL)  mH1[tbn+"_cA1"]->Fill(Amp[0]);
      if(mH1[tbn+"_cA2"]!=NULL)  mH1[tbn+"_cA2"]->Fill(Amp[1]);
      if(mH1[tbn+"_cA3"]!=NULL)  mH1[tbn+"_cA3"]->Fill(Amp[2]);
      if(mH2[tbn+"_cMap"]!=NULL)  mH2[tbn+"_cMap"]->Fill(pixU, pixV);
      if(mH2[tbn+"_cAmpRatio"]!=NULL)  mH2[tbn+"_cAmpRatio"]->Fill(Amp23,Amp13);
    }
    
    // Cut on amplitude ratios of clusters
    if (fAmpRatio13_.size()>0 && fAmpRatio23_.size()>0) {
      if (!CsGEM::IsInside(Amp23,Amp13,fAmpRatio23_,fAmpRatio13_)) continue;
    }
    
    // Fill histograms after time cut
    if (hLevel_>=Normal){
      if(mH1[tbn+"_cA1T"]!=NULL)  mH1[tbn+"_cA1T"]->Fill(Amp[0]);
      if(mH1[tbn+"_cA2T"]!=NULL)  mH1[tbn+"_cA2T"]->Fill(Amp[1]);
      if(mH1[tbn+"_cA3T"]!=NULL)  mH1[tbn+"_cA3T"]->Fill(Amp[2]);
      if(mH2[tbn+"_cMapT"]!=NULL)  mH2[tbn+"_cMapT"]->Fill(pixU, pixV);
      if(mH1[tbn+"_cSizeT"]!=NULL) mH1[tbn+"_cSizeT"]->Fill((*itclus)->GetSize());
      if(mH2[tbn+"_cAmpRatioT"]!=NULL)  mH2[tbn+"_cAmpRatioT"]->Fill(Amp23,Amp13);
    }

    // fill histogram before cut on crosstalk probability
    if (hLevel_>=High)
      if(mH1[tbn+"_cCTProbT"]!=NULL) mH1[tbn+"_cCTProbT"]->Fill( (*itclus)->GetXTalk() );

    // Cut on cross talk probability
    if ( (*itclus)->GetXTalk() > fCrossTalkParams[1] )
      continue;

    // fill histogram after cut on crosstalk probability
    if (hLevel_>=Normal)
      if(mH1[tbn+"_cCTProbC"]!=NULL) mH1[tbn+"_cCTProbC"]->Fill( (*itclus)->GetXTalk() );

    ngoodclus++;

    // Set coordinates in detector reference system
    const double u(pixU*wirPU_+wireDCorrU); // column
    const double v(pixV*wirPV_+wireDCorrV); // rows
    double w = zcm_;

    // Set errors
    HepMatrix cov(3,3,0);  // Zero matrix
    if((*itclus)->GetPositionErr() < 0 ) { // if cluster error is unknown
      cov(1,1) = pow( wirPU_ / sqrt(12.0), 2 );
      cov(2,2) = pow( wirPV_ / sqrt(12.0), 2 );
    } else {
      cov(1,1) = pow( (*itclus)->GetPositionXErr()*wirPU_, 2 );
      cov(2,2) = pow( (*itclus)->GetPositionYErr()*wirPV_, 2 );
    }
    cov(3,3) = 1.; // assume 1 mm resolution in Z
    
    // Create CsCluster 
    CsCluster* cluster = new CsCluster( u, v, w, cov );

    // Save CsDigits the CsCluster is made of
    std::list<CsGEMHit*> hits = (*itclus)->GetHits();
    std::list<CsGEMHit*>::iterator ithit;
    for (ithit = hits.begin(); ithit != hits.end(); ithit++) {
      cluster->addDigit( *(digitmap[(*ithit)->GetChan()->GetId()->GetDetectorChannel()]) );
    }

    // Store three amplitudes
    cluster->addAnalogData( (*itclus)->GetAmp()[0], (*itclus)->GetNoise() );
    cluster->addAnalogData( (*itclus)->GetAmp()[1], (*itclus)->GetNoise() );
    cluster->addAnalogData( (*itclus)->GetAmp()[2], (*itclus)->GetNoise() );
    
    // Store x and y pos
    cluster->addAnalogData( pixU );
    cluster->addAnalogData( pixV );

    // store crosstalk probability
    cluster->addAnalogData( (*itclus)->GetXTalk() );

    // Set detector
    cluster->addDet( *this );
    
    // Store cluster time (TCS phase corrected)
    double time, etime;
    if ( (*itclus)->GetTime(time, etime) ) {
      time = event->getTCSPhaseTime() - time;
      if (isMC) time -= triggerOffset;   // MC: Correct for trigger offset
      cluster->setTime(time, etime);
    } else
      cluster->setTime(0., 1.e9);

    // Fill histograms after time cut
    if (hLevel_>=Normal){ 
      if(mH1[tbn+"_cTimeT"]!=NULL) mH1[tbn+"_cTimeT"]->Fill(time);      
    }
    
    // Add cluster
    addCluster( *cluster );
    
    // Print out clusters
    if (fDebugLevel>9) {
      // 	  CsErrLog::msg(elFatal,__FILE__,__LINE__,
      // 		      "%s: Created cluster at %d, %d, amplitudes %d, %d, %d",
      // 		      tbn.c_str(),u,v,
      // 		      (*Id)->getData()[0],
      // 		      (*Id)->getData()[1],
      // 		      (*Id)->getData()[2]);
    }
  } // end of loop over found clusters

  // Fill histogram for cluster mutliplicities before and after time cut
  if (hLevel_>=Normal){  
    if(mH1[tbn+"_nClusT"]!=NULL)  mH1[tbn+"_nClusT"]->Fill(ngoodclus);
  }
  if (hLevel_>=High){  
    if(mH1[tbn+"_nClus"]!=NULL)  mH1[tbn+"_nClus"]->Fill(nclus);
  }
  
  sortClusters();
  setClusteringDone();
  
  // Clean up
  //  delete plane;
  pixelplane->Clear();

  // Retrieve clusters
  //  list<CsCluster*>cls;
  //  list<CsCluster*>::iterator itcls;
  //  cls = getMyClusters();
  //  for (itcls=cls.begin(); itcls!=cls.end(); itcls++) {
  //    cout << "CsCluster " << (*itcls)->getU()
  //	 << " " << (*itcls)->getV()
  //	 << " " << (*itcls)->getW() << endl;
  //  }  
}

//-----------------------------------------------------------------------------
// Clusterization for strip part of PixelGEM detector
//-----------------------------------------------------------------------------
void CsPixelGEMDetector::clusterizeStrips() {
  
  string tbn  = GetTBName();     
  
  // CORAL digits
  list<CsDigit*>::iterator Id;
  list<CsDigit*> digits = getMyDigits();
  vector<list<CsDigit*>::iterator> iterators;
  iterators.clear();
  
  // Protection
  if( digits.empty() ) return;
//   if ( gem_chan_cals.GetSize()!=ChipsPerPlane*ChipChannels ) {
//     CsErrLog::msg(elFatal,__FILE__,__LINE__,
// 		  "%s: Wrong calibrations size: %i",
// 		  GetTBName().c_str(),gem_chan_cals.GetSize());
//   }
//   if ( time_cals.GetSize()==0 ){
//     CsErrLog::msg(elFatal,__FILE__,__LINE__,
// 		  "%s: Timing calibrations missing!",
// 		  GetTBName().c_str());
//   }
    
  int err;
  HepMatrix iRotM(3,3); iRotM = rotWRS_.inverse( err );
  double wireDCorr = wirD_ + iRotM(1,1) * _deltax + iRotM(1,2) *_deltay; 
  
  // Map to store CORAL digits
  map<pair<int,int>,CsDigit*> digitmap;
  
  // Loop over digits and fill them as hits into the plane
  for( Id=digits.begin(); Id!=digits.end(); Id++ ) {
    
    // Make sure that digit belongs to this detector
    assert( (*Id)->getDet() == this );
    
    // Skip if already clusterized
    if( (*Id)->isClusterized() ) continue;
    (*Id)->setClusterized();
    
    // Get amplitudes
    const float a0 = (float)(*Id)->getData()[0];
    const float a1 = (float)(*Id)->getData()[1];
    const float a2 = (float)(*Id)->getData()[2];
    
    // Detector channel
    const int channel = (int)(*Id)->getAddress();
    
    // Hemisphere
    int hem = (int)(*Id)->getData()[3];
    
    // Add hit to plane
    stripplane->AddHit(channel, hem, a0, a1, a2); 
    
    // Add reference to CsDigit
    pair<int,int> id(channel,hem);
    digitmap[id] = *Id;
    
  } // end of interface
  
  // Fill histogram with number of hits
  if (hLevel_>=High){ 
    if(mH1[tbn+"_nHit"]!=NULL)  mH1[tbn+"_nHit"]->Fill(stripplane->GetNhit());
  }
  
  // Clustering
  stripplane->Clusterize();
  
  // Display and dump
  if (fDebugLevel>9) {
    stripplane->PrintClusters();
    stripplane->Display();
  }

  // Interface to CORAL clusters
  int nclus = stripplane->GetNcluster();
  int ngoodclus = 0;

  double triggerOffset(0.); // MC: Time offset of the trigger w.r.t. the event
  const CsEvent *event = CsEvent::Instance();
  bool isMC = event->isAMonteCarloEvent();
  if (isMC) triggerOffset = event->getTriggerMCOffset();

  // Get found clusters
  list<CsGEMCluster*> clusters = stripplane->GetClusters();
  
  // Loop over clusters to fill CORAL clusters
  list<CsGEMCluster*>::iterator itclus;
  for (itclus=clusters.begin(); itclus!=clusters.end(); itclus++) {
    
    // Get cluster amplitudes
    std::vector<Float_t> Amp = (*itclus)->GetAmp();
    float Amp13 = (Amp[2]==0) ? 0 : 
      (Amp[0]/Amp[2]);
    float Amp23 = (Amp[2]==0) ? 0 :
      (Amp[1]/Amp[2]);

    // Fill histograms before time cut
    if (hLevel_>=High){ 
      if(mH1[tbn+"_cA1"]!=NULL)  mH1[tbn+"_cA1"]->Fill(Amp[0]);
      if(mH1[tbn+"_cA2"]!=NULL)  mH1[tbn+"_cA2"]->Fill(Amp[1]);
      if(mH1[tbn+"_cA3"]!=NULL)  mH1[tbn+"_cA3"]->Fill(Amp[2]);
      if(mH1[tbn+"_cPos"]!=NULL)  mH1[tbn+"_cPos"]->Fill((*itclus)->GetPosition());
      if(mH2[tbn+"_cAmpRatio"]!=NULL)  mH2[tbn+"_cAmpRatio"]->Fill(Amp23,Amp13);
    }
    
    // Cut on amplitude ratios of clusters
    if (fAmpRatio13_.size()>0 && fAmpRatio23_.size()>0) {
      if (!CsGEM::IsInside(Amp23,Amp13,fAmpRatio23_,fAmpRatio13_)) continue;
    }
    
    // Fill histograms after time cut
    if (hLevel_>=Normal){ 
      if(mH1[tbn+"_cA1T"]!=NULL)  mH1[tbn+"_cA1T"]->Fill(Amp[0]);
      if(mH1[tbn+"_cA2T"]!=NULL)  mH1[tbn+"_cA2T"]->Fill(Amp[1]);
      if(mH1[tbn+"_cA3T"]!=NULL)  mH1[tbn+"_cA3T"]->Fill(Amp[2]);
      if(mH1[tbn+"_cPosT"]!=NULL)  mH1[tbn+"_cPosT"]->Fill((*itclus)->GetPosition());
      if(mH1[tbn+"_cSizeT"]!=NULL) mH1[tbn+"_cSizeT"]->Fill((*itclus)->GetSize());
      if(mH2[tbn+"_cAmpRatioT"]!=NULL)  mH2[tbn+"_cAmpRatioT"]->Fill(Amp23,Amp13);
    }

    // fill histogram before cut on crosstalk probability
    if (hLevel_>=High)
        if(mH1[tbn+"_cCTProbT"]!=NULL) mH1[tbn+"_cCTProbT"]->Fill( (*itclus)->GetXTalk() );

    // Cut on cross talk probability
    if ( (*itclus)->GetXTalk() > fCrossTalkParams[1] )
        continue;

    // fill histogram after cut on crosstalk probability
    if (hLevel_>=Normal)
        if(mH1[tbn+"_cCTProbC"]!=NULL) mH1[tbn+"_cCTProbC"]->Fill( (*itclus)->GetXTalk() );

    ngoodclus++;

    // Set coordinates in detector reference system
    double u;
    if ( IsVarPitch() )
      u = wireDCorr + Wire2Pos( (*itclus)->GetPosition() );  // perp. to wire
    else
      u = wireDCorr + (*itclus)->GetPosition() * wirP_;  // perp. to wire
    double v = 0;                                      // parallel to wire
    double w = zcm_;                                   // z-position of det.
    
    // Set errors
    HepMatrix cov(3,3,0);  // Zero matrix
    if((*itclus)->GetPositionErr() < 0 ) { // if cluster error is unknown
      if ( IsVarPitch() ) cov(1,1) = pow( Pitch((*itclus)->GetPosition()) / sqrt(12.0), 2 );
      else                cov(1,1) = pow( wirP_ / sqrt(12.0), 2 );
    } else {
      if ( IsVarPitch() ) cov(1,1) = pow( (*itclus)->GetPositionErr()*Pitch((*itclus)->GetPosition()), 2 );
      else                cov(1,1) = pow( (*itclus)->GetPositionErr()*wirP_, 2 );
    }
    cov(2,2) = pow(max(getXsiz(), getYsiz()),2); // just very big error
    cov(3,3) = 10.;
    
    // Create CsCluster
    CsCluster* cluster = new CsCluster( u, v, w, cov );
    
    // Save CsDigits the CsCluster is made of
    list<CsGEMHit*> hits = (*itclus)->GetHits();
    list<CsGEMHit*>::iterator ithit;
    for (ithit=hits.begin(); ithit!=hits.end(); ithit++) {
      pair<int,int> id((*ithit)->GetChan()->GetId()->GetDetectorChannel(),
		       (*ithit)->GetChan()->GetId()->GetHemisphere());
      //      CsDigit* dig = digitmap[id];
      //      cout << "Saving digit " << dig->getAddress() 
      //	   << " " << dig->getData()[0]
      //	   << " " << dig->getData()[1]
      //	   << " " << dig->getData()[2] << endl;
      cluster->addDigit( *(digitmap[id]) );
    }
    
    // Store three amplitudes
    cluster->addAnalogData( (*itclus)->GetAmp()[0], (*itclus)->GetNoise() );
    cluster->addAnalogData( (*itclus)->GetAmp()[1], (*itclus)->GetNoise() );
    cluster->addAnalogData( (*itclus)->GetAmp()[2], (*itclus)->GetNoise() );

    // Store wire position
    cluster->addAnalogData( (*itclus)->GetPosition(), (*itclus)->GetHemisphere() );

    // store crosstalk probability
    cluster->addAnalogData( (*itclus)->GetXTalk() );

    cluster->addDet( *this );

    // Store cluster time (TCS phase corrected)
    double time, etime;
    if ( (*itclus)->GetTime(time, etime) ) {
      time = event->getTCSPhaseTime() - time;
      if (isMC) time -= triggerOffset;   // MC: Correct for trigger offset
      cluster->setTime(time, etime);
    } else
      cluster->setTime(0., 1.e9);

    // Fill histograms after time cut
    if (hLevel_>=Normal){ 
      if(mH1[tbn+"_cTimeT"]!=NULL) mH1[tbn+"_cTimeT"]->Fill(time);      
    }
    
    // Add cluster
    addCluster( *cluster );
    
  } // end of loop over found clusters
  
  // Fill histogram for cluster mutliplicities before and after time cut
  if (hLevel_>=Normal){  
    if(mH1[tbn+"_nClusT"]!=NULL)  mH1[tbn+"_nClusT"]->Fill(nclus);
  }
  if (hLevel_>=High){  
    if(mH1[tbn+"_nClus"]!=NULL)  mH1[tbn+"_nClus"]->Fill(ngoodclus);
  }
  
  sortClusters();
  setClusteringDone();
  
  // Clean up
  stripplane->Clear();
}

//     
// Local class for detector response simulation
// (Used in CsPixelGEMDetector::makeMCDecoding())
//

// tempopary digits
class pGMdig
{
public:
  int    wire;  // wire #
  int    hemi;  // hemisphere (-1 for lower half, +1 for upper half)
  double amp1;  // amp1
  double amp2;  // amp2
  double amp3;  // amp3
  CsMCHit* ref; // reference to MCHit
  bool operator < (const pGMdig& gd) const
  {
    return (wire==gd.wire ? hemi<gd.hemi : wire<gd.wire);
  };
}; 

void getAmp(const double& t, double amps[] );

void CsPixelGEMDetector::makeMCDecoding() {

  // Apply mapping to channel calibrations if not done yet
  // cannot be done in CsGEMDetector::readCalibration because 
  // CS::Chip::Maps are not yet available
  if (!apv_chan_cals_mapped) apv_chan_cals_mapped = mapMCChannelCal();

  if (!decode_ && !decodeCard_) return;   // Should I proceed?
  if (decodingDone_) return;              // Already done?
  myDigits_.clear();                      // Clear

  if (CsGeant3::Instance()->isAnNtFile()) {
    // COMGeant NTuples (an obsolete type of output used in earlier COMPASS
    // times) are not supported by the CsPG class 
    CsErrLog::msg(elFatal,__FILE__,__LINE__,
		  "%s: No MC decoding COMGeant NTuple",GetTBName().c_str());
  }

  if (nWirV_) { //                 Pixelised CsPixelGEMDetector
    if (!amplitudeMCDecoding()) // Use simplistics decoding
      makeMCPixelsSimple();
    else  //                       Use amplitude decoding
      makeMCPixelsAmplitude();
  } else { //                      CsPixelGEMDetector of strip kind
    if (!amplitudeMCDecoding())
      makeMCStripsSimple();
    else
      makeMCStripsAmplitude();
  }

  decodingDone_ = true;
}

void CsPixelGEMDetector::makeMCPixelsSimple() {
  list<CsMCHit*>::iterator Ih;
  for (Ih = myMCHits_.begin(); Ih!=myMCHits_.end(); Ih++ ) { // loop on hits
	if ((((*Ih)->getMCTrack())->getParticle())->getCharge()==0 && (*Ih)->getOrigin()==0) 
      continue;

	// ***** LOOP ON HITS from CHARGED PARTICLES or CHARGED PRODUCTS *****

	CsMCTrkHit *thit = dynamic_cast<CsMCTrkHit*>(*Ih);
	if (thit==0) return;

	double t  = thit->getDTime();  // Delay time (ns)
	double ui = thit->getUin();    // Hit in point (DRS)
	double vi = thit->getVin();
	double wi = thit->getWin();
	double uo = thit->getUout();   // Hit out point (DRS)
	double vo = thit->getVout();
	double wo = thit->getWout();

	//    ***** TIMING: SMEARING, TIME CUT 
	double tdc = t+tRes_*CsRandom::gauss();
	if (fMCHitTMin>tdc || tdc>fMCHitTMax) 
      continue;

	HepMatrix iRotM(3,3); iRotM = getRotWRSInv();
	double xi = rotDRS_(1,1)*ui+rotDRS_(1,2)*vi+rotDRS_(1,3)*wi+ xcm_ - _xshift;
	double yi = rotDRS_(2,1)*ui+rotDRS_(2,2)*vi+rotDRS_(2,3)*wi+ ycm_ - _yshift;
	double zi = rotDRS_(3,1)*ui+rotDRS_(3,2)*vi+rotDRS_(3,3)*wi+ zcm_;
	double xo = rotDRS_(1,1)*uo+rotDRS_(1,2)*vo+rotDRS_(1,3)*wo+ xcm_ - _xshift;
	double yo = rotDRS_(2,1)*uo+rotDRS_(2,2)*vo+rotDRS_(2,3)*wo+ ycm_ - _yshift;
	double zo = rotDRS_(3,1)*uo+rotDRS_(3,2)*vo+rotDRS_(3,3)*wo+ zcm_;
	double Ui = iRotM(1,1)*xi+iRotM(1,2)*yi+iRotM(1,3)*zi; // WRS
	double Vi = iRotM(2,1)*xi+iRotM(2,2)*yi+iRotM(2,3)*zi;
	double Wi = iRotM(3,1)*xi+iRotM(3,2)*yi+iRotM(3,3)*zi;
	double Uo = iRotM(1,1)*xo+iRotM(1,2)*yo+iRotM(1,3)*zo; // WRS
	double Vo = iRotM(2,1)*xo+iRotM(2,2)*yo+iRotM(2,3)*zo;
	double Wo = iRotM(3,1)*xo+iRotM(3,2)*yo+iRotM(3,3)*zo;

	int iu, iv, pixUF, pixVF;  // First X/Y raw for this hit
	if ((Ui-wirDU_)/wirPU_<0) pixUF = int((Ui-wirDU_)/wirPU_-0.5);
	else                      pixUF = int((Ui-wirDU_)/wirPU_+0.5);
	if ((Vi-wirDV_)/wirPV_<0) pixVF = int((Vi-wirDV_)/wirPV_-0.5);
	else                      pixVF = int((Vi-wirDV_)/wirPV_+0.5);
	int pixUL, pixVL;  // Last  X/Y raw for this hit
	if ((Uo-wirDU_)/wirPU_<0) pixUL = int((Uo-wirDU_)/wirPU_-0.5);
	else                      pixUL = int((Uo-wirDU_)/wirPU_+0.5);
	if ((Vo-wirDV_)/wirPV_<0) pixVL = int((Vo-wirDV_)/wirPV_-0.5);
	else                      pixVL = int((Vo-wirDV_)/wirPV_+0.5);
	if (pixUL<pixUF) {
	  int tmp = pixUL; pixUL = pixUF; pixUF = tmp;
	}
	if (pixVL<pixVF) {
	  int tmp = pixVL; pixVL = pixVF; pixVF = tmp;
	}

	for (iu = pixUF; iu<=pixUL; iu++) for (iv = pixVF; iv<=pixVL; iv++) {

	  // ********** DOES A CsDigit ALREADY EXIST ON THIS PIXEL? **********

	  list<CsDigit*>::iterator Id; bool found;
	  for (Id = myDigits_.begin(), found = false; Id!=myDigits_.end(); Id++) {
	    int raw = (*Id)->getAddress(), col = raw%nWirU_; raw /= nWirU_; 
	    if (iu==col && iv==raw) {        // Here it is....
	      found = true;	 // ***** ...ADD THIS HIT TO EXISTING CsDigit...
	      dynamic_cast<CsMCDigit*>(*Id)->addHit(*(*Ih));
          double info[3];
          getAmp(tdc, info);
          for (unsigned int infoidx=0; infoidx<3; infoidx++)
            (*Id)->getData()[infoidx]+=info[infoidx]*20000.;
	      break;
	    }
	  }
	  if (!found) { // No digits found with these wire and detector
	    if (iu<0 || nWirU_<=iu || iv<0 || nWirV_<=iv) {
	      CsErrLog::msg(elAnomaly,__FILE__,__LINE__,
			    "%s: Pixel(=%d,%d) outside range [0,%d[x[0,%d[",
			    GetTBName().c_str(),iu,iv,nWirU_,nWirV_); continue;
	    } else {                                        // ***** NEW CsDigit
	      // Lump the 2 pixel coord's into a single CsDigit's address
          double info[4];
          getAmp(tdc, info);
          for (unsigned int infoidx=0; infoidx<3; infoidx++)
            info[infoidx] *= 20000.;
          info[3] = (iv<<10)+iu;
	      CsDigit *digit = new CsMCDigit(*this,iu+iv*nWirU_,info, 4);
	      dynamic_cast<CsMCDigit*>(digit)->addHit(*(*Ih));
	      myDigits_.push_back( digit );            // ***** ADD IT TO LIST
	    }
	  }	  	  
	}  // End loop on hits associated to current MC hit
  }  // End of loop on MC hits

  // fill histograms
  if ( hLevel_ >= High ) {
    list<CsDigit*>::iterator Id;
    for (Id = myDigits_.begin(); Id!=myDigits_.end(); Id++) {
      double amp1 = (*Id)->getData()[0];
      double amp2 = (*Id)->getData()[1];
      double amp3 = (*Id)->getData()[2];
      string tbn  = GetTBName();
      if(mH1[tbn+"_MCdA1"]!=NULL)  mH1[tbn+"_MCdA1"]->Fill(amp1);
      if(mH1[tbn+"_MCdA2"]!=NULL)  mH1[tbn+"_MCdA2"]->Fill(amp2);
      if(mH1[tbn+"_MCdA3"]!=NULL)  mH1[tbn+"_MCdA3"]->Fill(amp3);
      if(mH2[tbn+"_MCdAmpRatio"]!=NULL)  mH2[tbn+"_MCdAmpRatio"]->Fill(amp2/amp3, amp1/amp3);
    }
  }
} 

void CsPixelGEMDetector::makeMCPixelsAmplitude() {
  CsErrLog::msg(elFatal, __FILE__, __LINE__, "Amplitude decoding for pixels not yet included");
}

void CsPixelGEMDetector::makeMCStripsSimple() {
  list<pGMdig> ld; // list of temporary digits 

  list<CsMCHit*>::iterator Ih;
  for( Ih=myMCHits_.begin(); Ih!=myMCHits_.end(); Ih++ ) { // loop on hits

	//only charged particles or charged products
	if(  (((*Ih)->getMCTrack())->getParticle())->getCharge() || (*Ih)->getOrigin()  ) {
	  	  
	  CsMCTrkHit* thit = dynamic_cast<CsMCTrkHit*>(*Ih);
	  if( thit == 0 ) return;
	  
	  double t  = thit->getDTime();  // Delay time (ns)
	  double ui = thit->getUin();    // Hit in point (DRS)
	  double vi = thit->getVin();
	  double wi = thit->getWin();
	  double uo = thit->getUout();   // Hit out point (DRS)
	  double vo = thit->getVout();
	  double wo = thit->getWout();
	  
	  // ***** TIME SMEARING:
	  // Detector resolution is hard-coded for the time being (03/02)
	  //double tdc = t + tRes_ * CsRandom::gauss();
	  double tdc = t + 12. * CsRandom::gauss();
	  if (fMCHitTMin>tdc || tdc>fMCHitTMax) continue; // ***** REQUIRE w/in T GATE
	  // Check if the hit is inside the time gate
	  //	if ( (-tGate_/2) > t ||  t > (tGate_/2) ) continue;

	  int err;
	  HepMatrix iRotM(3,3); iRotM = rotWRS_.inverse( err );
	  double xi = rotDRS_(1,1)*ui+rotDRS_(1,2)*vi+rotDRS_(1,3)*wi + xcm_ - _xshift;
	  double yi = rotDRS_(2,1)*ui+rotDRS_(2,2)*vi+rotDRS_(2,3)*wi + ycm_ - _yshift;
	  double zi = rotDRS_(3,1)*ui+rotDRS_(3,2)*vi+rotDRS_(3,3)*wi + zcm_;
	  double xo = rotDRS_(1,1)*uo+rotDRS_(1,2)*vo+rotDRS_(1,3)*wo + xcm_ - _xshift;
	  double yo = rotDRS_(2,1)*uo+rotDRS_(2,2)*vo+rotDRS_(2,3)*wo + ycm_ - _yshift;
	  double zo = rotDRS_(3,1)*uo+rotDRS_(3,2)*vo+rotDRS_(3,3)*wo + zcm_;
	  double Ui = iRotM(1,1)*xi+iRotM(1,2)*yi+iRotM(1,3)*zi; // WRS
	  double Vi = iRotM(2,1)*xi+iRotM(2,2)*yi+iRotM(2,3)*zi;
	  double Wi = iRotM(3,1)*xi+iRotM(3,2)*yi+iRotM(3,3)*zi;
	  double Uo = iRotM(1,1)*xo+iRotM(1,2)*yo+iRotM(1,3)*zo; // WRS
	  double Vo = iRotM(2,1)*xo+iRotM(2,2)*yo+iRotM(2,3)*zo;
	  double Wo = iRotM(3,1)*xo+iRotM(3,2)*yo+iRotM(3,3)*zo;
	  
	  int wireF;  // first wire for this hit
	  if( (Ui-wirD_)/wirP_ < 0 ) {
	    wireF = int( (Ui-wirD_)/wirP_-0.5 );
	  } else {
	    wireF = int( (Ui-wirD_)/wirP_+0.5 );
	  }
	  int wireL; // last  wire for this hit
	  if( (Uo-wirD_)/wirP_ < 0 ) {
	    wireL = int( (Uo-wirD_)/wirP_-0.5 );
	  } else {
	    wireL = int( (Uo-wirD_)/wirP_+0.5 );
	  } 
	  if( wireL < wireF ) {
	    int tmp = wireL;
	    wireL = wireF;
	    wireF = tmp;
	  }

      // find the hemisphere in which the hit was found
      unsigned int hemneg(0), hempos(0);
      if (GetTBName()[4]=='U') {          // remember the left-handed
        if (vi < 0.) hempos=1;            // coordinate system
        else         hemneg=1;            // for UV strips
        if (vo < 0.) hempos=1;
        else         hemneg=1;
      }
      if (GetTBName()[4]=='V') {
        if (ui > 0.) hempos=1;
        else         hemneg=1;
        if (uo > 0.) hempos=1;
        else         hemneg=1;
      }
      if (GetTBName()[4]=='X') {
        if (vi > 0.) hempos=1;
        else         hemneg=1;
        if (vo > 0.) hempos=1;
        else         hemneg=1;
      }
      if (GetTBName()[4]=='Y') {
        if (ui > 0.) hempos=1;
        else         hemneg=1;
        if (uo > 0.) hempos=1;
        else         hemneg=1;
      }
	  
	  for( int i=wireF; i<=wireL; i++ ) { 
	    if( i<0 || i>=nWir_ ) {
	      ostringstream ost;
          ost << "Unreliable wire number: " << i << " (0," << nWir_ <<"), "
		      << " detector : " << GetID() << " " << unit_ << " " << type_ << ".";
          CsErrLog::Instance()->mes( elAnomaly, ost.str() );
        } else {
          double amps[3];
          getAmp(tdc, amps);
          if (hemneg) {
            pGMdig gd;
		    gd.wire = i;     
            gd.hemi = -1;
		    gd.amp1 = 20000.*amps[0] / (hemneg+hempos);
		    gd.amp2 = 20000.*amps[1] / (hemneg+hempos);
		    gd.amp3 = 20000.*amps[2] / (hemneg+hempos); 
		    gd.ref  =  (*Ih); // pointer to CsMCHit   
		    ld.push_back(gd); // save GM digit
          }
          if (hempos) {
            pGMdig gd;
		    gd.wire = i;     
            gd.hemi = 1;
		    gd.amp1 = 20000.*amps[0] / (hemneg+hempos);
		    gd.amp2 = 20000.*amps[1] / (hemneg+hempos);
		    gd.amp3 = 20000.*amps[2] / (hemneg+hempos); 
		    gd.ref  =  (*Ih); // pointer to CsMCHit   
		    ld.push_back(gd); // save GM digit
          }
        }
      }
    }
  }      //end of MCHits loop

  list<pGMdig>::iterator id,idnext;
  bool newdig=false;
  vector<CsMCHit*> dighits;
  double ampl[4]; int iwir;
  for(id = ld.begin(); id != ld.end(); id++) { // loop over digits
	
	idnext=id; idnext++;          // look on the previous digit
    if(! newdig ) {
	  dighits.clear();
      iwir=-1;
      ampl[0]=0.; ampl[1]=0.; ampl[2]=0.;
      newdig = true; 
	}
        
    iwir=(*id).wire;
    ampl[0]+=(*id).amp1;
    ampl[1]+=(*id).amp2;
    ampl[2]+=(*id).amp3;
    ampl[3] =(*id).hemi;
    dighits.push_back((*id).ref);
	
	if(idnext == ld.end()) {
      newdig=false; 
	} else if ( (*id).wire != (*idnext).wire || (*id).hemi != (*idnext).hemi ) {
	  newdig=false;
	} else {
	  continue;
	}

    if( ! newdig ){
      CsDigit* digit = new CsMCDigit(*this, iwir, ampl, 4) ;
      vector<CsMCHit*>::iterator ih;
      for(ih = dighits.begin(); ih != dighits.end(); ih++){
        dynamic_cast<CsMCDigit*>(digit)->addHit( *(*ih) );
      }
      myDigits_.push_back( digit );
	}
  }  
}

void CsPixelGEMDetector::makeMCStripsAmplitude() {
  pGMdig gd;          
  list<pGMdig> ld; // list of temporary digits 

  list<CsMCHit*>::iterator Ih;
  for( Ih=myMCHits_.begin(); Ih!=myMCHits_.end(); Ih++ ) { // loop on hits

	//only charged particles or charged products
	if(  (((*Ih)->getMCTrack())->getParticle())->getCharge() || (*Ih)->getOrigin()  ) {
	  	  
	  CsMCTrkHit* thit = dynamic_cast<CsMCTrkHit*>(*Ih);
	  if( thit == 0 ) return;
	  
      double t  = thit->getDTime();  // Delay time (ns)
      int ndig=0;

      // ***** TIME SMEARING:
      // Detector resolution done in CsGeant3::readGeantHits()
      // if amplitude correlations are ON
      double tdc=t;
      if(!doAmpCorrelationMC()) {
	    tdc = t + tRes_ * CsRandom::gauss();
      } 

      // ***** REQUIRE w/in fiducial T Gate, reasonable: -300+100ns
	  if (fMCHitTMax<tdc || tdc<fMCHitTMin) continue; 

	  double ui = thit->getUin();    // Hit in point (DRS)
	  double vi = thit->getVin();
	  double wi = thit->getWin();
	  double uo = thit->getUout();   // Hit out point (DRS)
	  double vo = thit->getVout();
	  double wo = thit->getWout();
      double eloss = thit->getELos()*1.E6; // dE in detector, keV


      if( eloss <= 0. ) continue;

	  int err;
	  HepMatrix iRotM(3,3); iRotM = rotWRS_.inverse( err );
	  double xi = rotDRS_(1,1)*ui+rotDRS_(1,2)*vi+rotDRS_(1,3)*wi+ xcm_ - _xshift;
	  double yi = rotDRS_(2,1)*ui+rotDRS_(2,2)*vi+rotDRS_(2,3)*wi+ ycm_ - _yshift;
	  double zi = rotDRS_(3,1)*ui+rotDRS_(3,2)*vi+rotDRS_(3,3)*wi+ zcm_;
	  double xo = rotDRS_(1,1)*uo+rotDRS_(1,2)*vo+rotDRS_(1,3)*wo+ xcm_ - _xshift;
	  double yo = rotDRS_(2,1)*uo+rotDRS_(2,2)*vo+rotDRS_(2,3)*wo+ ycm_ - _yshift;
	  double zo = rotDRS_(3,1)*uo+rotDRS_(3,2)*vo+rotDRS_(3,3)*wo+ zcm_;
	  double Ui = iRotM(1,1)*xi+iRotM(1,2)*yi+iRotM(1,3)*zi; // WRS
	  double Vi = iRotM(2,1)*xi+iRotM(2,2)*yi+iRotM(2,3)*zi;
	  double Wi = iRotM(3,1)*xi+iRotM(3,2)*yi+iRotM(3,3)*zi;
	  double Uo = iRotM(1,1)*xo+iRotM(1,2)*yo+iRotM(1,3)*zo; // WRS
	  double Vo = iRotM(2,1)*xo+iRotM(2,2)*yo+iRotM(2,3)*zo;
	  double Wo = iRotM(3,1)*xo+iRotM(3,2)*yo+iRotM(3,3)*zo;

      const double s2pi=2.5066283;   // sqrt(2*pi)
      const double nfactor=eGain_;    // some normalization factor "gain"
      const double sigspace = sWidth_;  // signal width, mm (tuned to have 3.2 strips/cluster)=0.3
      const double sig2   = sigspace*sigspace*2.;
    
      double ltrack   = fabs(Uo-Ui);
      double aN;                                     // integral
      int dummyCounter = 0;
	  // randomize it
	  do{
	    aN = eloss*(nfactor+eGSig_*CsRandom::gauss()) ; 
	  } while( aN<0. && dummyCounter++<10 );
      if(aN<=0.) aN = eloss*nfactor;                  // need to keep it 
      double amp      = aN/sigspace/s2pi;             // amplitude
      double norm     = aN+amp*ltrack;                // new integral
      amp = amp*aN/norm;  

      double aUo=0.;
      double aUi=0.;
      double Utmin=0.;
      double Utmax=0.;
	  double resSpace = spSig_*CsRandom::gauss() ;  // space resolution

	  if(  Uo-Ui >0 ) {      // take into account diffusion, sort aUo ">" aUi   
	    aUo = Uo+3.*sigspace+resSpace ;
        aUi = Ui-3.*sigspace+resSpace ;
        Utmin=Ui+resSpace;
        Utmax=Uo+resSpace; 
	  } else {
        aUo = Ui+3.*sigspace+resSpace ;
        aUi = Uo-3.*sigspace+resSpace ;
        Utmin=Uo+resSpace;
        Utmax=Ui+resSpace;             
	  } 

	  int wireF;  // first wire for this hit
	  if( (aUi-wirD_)/wirP_ < 0 ) {
	    wireF = int( (aUi-wirD_)/wirP_-0.5 );
	  } else {
	    wireF = int( (aUi-wirD_)/wirP_+0.5 );
	  }
	  int wireL; // last  wire for this hit
	  if( (aUo-wirD_)/wirP_ < 0 ) {
	    wireL = int( (aUo-wirD_)/wirP_-0.5 );
	  } else {
	    wireL = int( (aUo-wirD_)/wirP_+0.5 );
	  } 
	  if( wireL < wireF ) {
	    int tmp = wireL;
	    wireL = wireF;
	    wireF = tmp;
	  }

      // find the hemisphere in which the hit was found
      unsigned int hemneg(0), hempos(0);
      if (GetTBName()[4]=='U') {          // remember the left-handed
        if (vi < 0.) hempos=1;            // coordinate system
        else         hemneg=1;            // for UV strips
        if (vo < 0.) hempos=1;
        else         hemneg=1;
      }
      if (GetTBName()[4]=='V') {
        if (ui > 0.) hempos=1;
        else         hemneg=1;
        if (uo > 0.) hempos=1;
        else         hemneg=1;
      }
      if (GetTBName()[4]=='X') {
        if (vi > 0.) hempos=1;
        else         hemneg=1;
        if (vo > 0.) hempos=1;
        else         hemneg=1;
      }
      if (GetTBName()[4]=='Y') {
        if (ui > 0.) hempos=1;
        else         hemneg=1;
        if (uo > 0.) hempos=1;
        else         hemneg=1;
      }

	  for( int i=wireF; i<=wireL; i++ ) { 
	  
	    if( i<0 || i>=nWir_ ) {
	      ostringstream ost;
	      ost << "Unreliable wire number: " << i << " (0," << nWir_ <<"), "
		      << " detector : " << GetID()  << " " << unit_ << " " << type_ <<".";
	      CsErrLog::Instance()->mes( elAnomaly, ost.str() );
	    } else {

	      double ampI; // total ampliutde for the strip "i" 
	      double Uwire =  i*wirP_+wirD_;           

	      if( Uwire <= Utmin ) {
	        ampI=amp*exp(-(Uwire-Utmin)*(Uwire-Utmin)/sig2);
	      } else if ( Utmin < Uwire && Uwire < Utmax ) {
		    ampI=amp;
	      } else {
		    ampI=amp*exp(-(Uwire-Utmax)*(Uwire-Utmax)/sig2); 
	      }

          double amps[3];
	      getAmp( tdc, amps );      
	      const float acut=10.;  // to remove meaningless hits (to be tuned)
          if( (amps[0]+amps[1]+amps[2])*ampI > acut ) {
            // create the hits in the correct hemisphere
            if (hemneg) {
              ndig++;
		      gd.wire = i;     
              gd.hemi = -1;
		      gd.amp1 = ampI*amps[0] / (hemneg+hempos);
		      gd.amp2 = ampI*amps[1] / (hemneg+hempos);
		      gd.amp3 = ampI*amps[2] / (hemneg+hempos); 
		      gd.ref  =  (*Ih); // pointer to CsMCHit   
		      ld.push_back(gd); // save GM digit
            }
            if (hempos) {
              ndig++;
		      gd.wire = i;     
              gd.hemi = 1;
		      gd.amp1 = ampI*amps[0] / (hemneg+hempos);
		      gd.amp2 = ampI*amps[1] / (hemneg+hempos);
		      gd.amp3 = ampI*amps[2] / (hemneg+hempos); 
		      gd.ref  =  (*Ih); // pointer to CsMCHit   
		      ld.push_back(gd); // save GM digit
            }
	      }
	    }
	  }	  

	  if( hLevel_ >= High )  {
	    string tbn  = GetTBName();
	    if(mH2[tbn+"_MCnDvsT"]!=NULL)  mH2[tbn+"_MCnDvsT"]->Fill(tdc,float(ndig));
	  }                   
	}    // end of charged particles
  }      // end of MCHits loop

  ld.sort();            // sort by wire# 

  // Here we've all what's needed for build digits...
      
  list<pGMdig>::iterator id,idnext;
  bool newdig=false;
  vector<CsMCHit*> dighits;
  double ampl[4]; int iwir;     
  for(id = ld.begin(); id != ld.end(); id++) { // loop over digits
	
	idnext=id; idnext++;          // look on the previous digit
    if(! newdig ) {
	  dighits.clear();
      iwir=-1;
      ampl[0]=0.; ampl[1]=0.; ampl[2]=0.;
      newdig = true; 
	}
        
    iwir=(*id).wire;
    ampl[0]+=(*id).amp1;
    ampl[1]+=(*id).amp2;
    ampl[2]+=(*id).amp3;
    ampl[3] =(*id).hemi;
    dighits.push_back((*id).ref);
	
	if(idnext == ld.end()) {
      newdig=false; 
	} else if ( (*id).wire != (*idnext).wire || (*id).hemi != (*idnext).hemi ) {
	  newdig=false;
	} else {
	  continue;
	}

    if( ! newdig ){
      CsDigit* digit = new CsMCDigit(*this, iwir, ampl, 4) ;
      vector<CsMCHit*>::iterator ih;
      for(ih = dighits.begin(); ih != dighits.end(); ih++){
        dynamic_cast<CsMCDigit*>(digit)->addHit( *(*ih) );
      }
      myDigits_.push_back( digit );

      if( hLevel_ >= High )  {
	    list<CsMCHit*> dhit = dynamic_cast<CsMCDigit*>(digit)->getHits(); 
        string tbn  = GetTBName();
	    float Amp13 = (ampl[2]==0) ? 0 : ampl[0]/ampl[2];
        float Amp23 = ( ampl[2]==0) ? 0 : ampl[1]/ampl[2];  
	    if(mH1[tbn+"_MCdA1"]!=NULL)  mH1[tbn+"_MCdA1"]->Fill(ampl[0]);
	    if(mH1[tbn+"_MCdA2"]!=NULL)  mH1[tbn+"_MCdA2"]->Fill(ampl[1]);
	    if(mH1[tbn+"_MCdA3"]!=NULL)  mH1[tbn+"_MCdA3"]->Fill(ampl[2]);
	    if(mH1[tbn+"_MCdnHit"]!=NULL)  mH1[tbn+"_MCdnHit"]->Fill(float(dhit.size()));
	    if(mH2[tbn+"_MCdAmpRatio"]!=NULL)  mH2[tbn+"_MCdAmpRatio"]->Fill(Amp23,Amp13);
	  }
	}
  }  
}

////////////////////////////////////////////////////////////////////////////////

istream & operator >> (istream &in,vector<CsPixelGEMDetector::APVCal> &c)
{
  c.clear();
//   for( size_t i=0; i<CsPixelGEMDetector::ChipsPerPlane; i++ )
//     c.push_back(CsPixelGEMDetector::APVCal());

  try
  {
    unsigned chip=9999999;
    while(1)
    {
      char s[111];
      in.getline(s,sizeof(s)); // read line with information about chips
      if( !in )
        break; // end of data stream

      if( s[0]=='#' )
      {
        char s2[sizeof(s)];
        int a, src_id, adc_id, chip_id;
        if( 5!=sscanf(s,"%s %d %d %d %d",s2,&a,&src_id,&adc_id,&chip_id) )
          throw CS::Exception("CsPixelGEMDetector::APVCal::operator>>: bad string \"%s\".",s);
        chip--;

	//       if( chip>=c.size() )
	//         throw CS::Exception("CsPixelGEMDetector::APVCal::operator>>:  Too big chip number %d, max=%d",
	//                              chip,c.size());

	CsPixelGEMDetector::APVCal cal;
	cal.src_id=src_id;
	cal.adc_id=adc_id;
	cal.chip_id=chip_id;
	
	for( int i=0; i<CsPixelGEMDetector::ChipChannels; i++ )
	  {
	    in.getline(s,sizeof(s)); // read line with information about chips
	    if( !in )
	      throw CS::Exception("CsPixelGEMDetector::APVCal::operator>>: Error in reading data.");
	    cal.channels.push_back(CsPixelGEMDetector::APVCal::Channel(s));
	  }
	
	c.push_back(cal);
      }
      else
        throw CS::Exception("CsPixelGEMDetector::APVCal::operator>>: format error.");
    }
  }
  catch(...)
  {
    cerr << "CsPixelGEMDetector::APVCal::operator>> error!\n";
    throw;
  }
    
  in.clear(); // Why should I do this? I don't know...

  return in;
}

///////////////////////////////////////////////////////////////////////////////

void CsPixelGEMDetector::readCalibration(time_t timePoint){
  
  // Clear all existing calibrations
  if (stripplane) stripplane->ClearChannels();
  else if(pixelplane) pixelplane->ClearChannels();
  apv_chan_cals.clear();
  CsGEMTimeCals time_cals; time_cals.Clear();
  apv_chan_cals_mapped = false;
  
  CDB::Time tp(timePoint,0);
  tm *t = localtime(&tp.first);

  // Read calibration data for each strip
  string strdata("");
  cdb_->read(GetTBName(), strdata, tp);
  istringstream istrdata(strdata);
  istrdata >> apv_chan_cals;
  
  // Check if calibrations for correct number of chips was found
  bool apv_chan_cals_ok = true;
  if((int)apv_chan_cals.size() == ChipsPerPlane){
    CsErrLog::msg(elInfo,__FILE__,__LINE__,
		  "%s calibration data found in CDB for local time %s",
		  GetTBName().c_str(),asctime(t));

    // Check if correct number of channels for each chip
    vector<CsPixelGEMDetector::APVCal>::iterator itcal;
    for (itcal = apv_chan_cals.begin(); itcal != apv_chan_cals.end(); itcal++) {
      vector<CsPixelGEMDetector::APVCal::Channel>& vchn = (*itcal).channels;
      if (vchn.size() != ChipChannels) {
	CsErrLog::msg(elError,__FILE__,__LINE__,
		      "%s wrong calibration data size (# channels) in CDB for local time %s",GetTBName().c_str(),asctime(t));
	apv_chan_cals_ok = false;
      }
    }

    // Debugging output
    if (getenv("COND_DB_DEBUG")!=0 || getenv("COND_DB_DEBUG_GEM")!=0) {
      cout << GetTBName() <<endl;
      cout << apv_chan_cals.size() << endl;
      vector<CsPixelGEMDetector::APVCal>::iterator itcal;
      vector<CsPixelGEMDetector::APVCal::Channel>::iterator itchn;
      for (itcal = apv_chan_cals.begin(); itcal != apv_chan_cals.end(); itcal++) {
        cout<<"next APVCal structure:\n";
        vector<CsPixelGEMDetector::APVCal::Channel>& vchn = (*itcal).channels;
        for (itchn = vchn.begin(); itchn != vchn.end(); itchn++) {
          (*itchn).Print();
        }
      }
      cout<<endl;
    }
  }

  // Calibration data not found for all chips
  else {
    CsErrLog::msg(elWarning,__FILE__,__LINE__,
		  "%s wrong calibration data size (# chips) in CDB for local time %s",
		  GetTBName().c_str(),asctime(t));
    apv_chan_cals_ok = false;
  }

  // Calibrations not found for every strip of plane
  if (!apv_chan_cals_ok) apv_chan_cals.clear();

  // Read calibration data for conversion of amplitude->time for each plane
  string strdata2("");
  cdb_->read(GetTBName(), strdata2, tp, "timing");
  istringstream istrdata2(strdata2);
  istrdata2 >> time_cals;
  
  if (time_cals.IsValid()) {
    CsErrLog::msg(elInfo,__FILE__,__LINE__,
                  "%s time calibration data found in CDB for local time %s",
		  GetTBName().c_str(),asctime(t));
  } 
  else {
    CsErrLog::msg(elWarning,__FILE__,__LINE__,
                  "%s no time calibration data found in CDB for local time %s",
		  GetTBName().c_str(),asctime(t));
    CsErrLog::msg(elWarning,__FILE__,__LINE__,
		  "Using default values.");
    
    // Default values for timing calibrations
    // please do not forget to update getAmp() if you change this !
    time_cals.Clear();
    CsGEMTimeCalOld *cal0 = new CsGEMTimeCalOld(0.2,1.1,-40.,22.,1.65);
    CsGEMTimeCalOld *cal1 = new CsGEMTimeCalOld(0.3,1.0,0.,27.,1.23);
    time_cals.Set(cal0,cal1);
  }
  if (stripplane) stripplane->SetTimeCals(time_cals);
  else if (pixelplane) pixelplane->SetTimeCals(time_cals); 
  
  // Just debug printouts
  if (getenv("COND_DB_DEBUG")!=0 || getenv("COND_DB_DEBUG_GEM")!=0) {
    cout<<endl<<"Timing calibrations read in for "<<GetTBName()<<endl;
    if (stripplane) {
      const CsGEMTimeCals *cals = stripplane->GetTimeCals();
      cals->Print();
    }
    else if (pixelplane) {
      const CsGEMTimeCals *cals = pixelplane->GetTimeCals();
      cals->Print();
    }
  }
}

//-----------------------------------------------------------------------------
// Map channel calibrations
//-----------------------------------------------------------------------------
bool CsPixelGEMDetector::mapChannelCal() {
  // All decoding maps
  const CS::Chip::Maps &daq_maps = 
    CsInit::Instance()->getDaqMaps(); 

  // Create a short name for iterator through map
  typedef CS::Chip::Maps::const_iterator m_it; 

  if ( daq_maps.size() != 0 ) {

    // Apply mapping to calibrations for this plane
    unsigned  gem_chan_cals_ok = 0;
    if (apv_chan_cals.size()!=0) {

      // Loop over all chips and get calibration data for each chip
      for(unsigned chip=0;chip<apv_chan_cals.size();chip++) {
        APVCal &c = apv_chan_cals[chip];
      
        // Get Source ID and ADC ID and Chip ID of chip
        int src_id=c.src_id, adc_id=c.adc_id, chip_id=c.chip_id; 
        int wire=-9999; int hem=0;
        int pixX=-1; int pixY=-1;
        //	cout << "Looking for source Id " << src_id 
        //	     << ", ADC Id " << adc_id 
        //	     << ", chip Id " << chip_id << endl;
      
	    // Loop over all channels of chip  
        for(unsigned chip_chan=0 ;chip_chan<c.channels.size();chip_chan++) {
	  
          // All maps with given data ID
	      CS::ChipAPV::DataID data_id(src_id,adc_id,chip_id,chip_chan);
	  
	      //	  cout << "Data Id " << data_id << endl;
	      const pair<m_it,m_it> m_range = 
	        daq_maps.equal_range(data_id); 
	  
	      // Read map
	      bool mapped=false;
	      for( m_it cc=m_range.first; cc!=m_range.second; cc++ ) {
	        const CS::ChipAPV::Digit *digit1 = 
	          dynamic_cast<CS::ChipAPV::Digit*>(cc->second); 
	    
	        if( digit1==NULL ){
	          CsErrLog::msg(elWarning,__FILE__,__LINE__,
			                "%s: ChipAPV wrong map!",GetTBName().c_str());
	          continue;
	        }
	    
	        if( GetTBName()!=digit1->GetDetID().GetName() ){
	          CsErrLog::msg(elError,__FILE__,__LINE__,
			                "%s: Inconsistency between mapping file and calibration file!",GetTBName().c_str());
	          CsErrLog::msg(elError,__FILE__,__LINE__,
			                "%s from calibration file does not match TBName",
			                digit1->GetDetID().GetName().c_str());
	          continue;
	        }
	    
	        // Map found
	    
	        // Check if pixel digit
	        const CS::ChipAPV::DigitPixel *digit1p = 
	          dynamic_cast<const CS::ChipAPV::DigitPixel*>(digit1);
	    
	        // Strip digit
	        if( digit1p==NULL ) {
	      
	          // Apply mapping
	          try { 
		        wire=digit1->GetChannel();
		        hem=digit1->GetChanPos();
	          }
	          catch(...) {
		        CsErrLog::msg(elError,__FILE__,__LINE__,
			                  "%s: Could not apply mapping!",
			                  GetTBName().c_str());
		        continue;
	          }
	          mapped=true;
	          break;
	        }
	    
	        // Pixel digit
	        else {
	      
	          // Apply mapping
	          try { 
		        wire=digit1p->GetChannel();
	          }
	          catch(...) {
		        CsErrLog::msg(elError,__FILE__,__LINE__,
			                  "%s: Could not apply mapping!",
			                  GetTBName().c_str());
		        continue;
	          }
	      
	          // Get pixel coordinates
	          pair<int,int> xy=pixgem::detch2xy(wire);
	          pixX = xy.first;
	          pixY = xy.second;
	      
	          // Pixel coordinates invalid
	          if (pixX==-1 && pixY==-1) 
		      {
		        CsErrLog::msg(elError,__FILE__,__LINE__,
				              "%s: Could not apply mapping!",
				              GetTBName().c_str());
		        continue;
		      }
	      
	          // Invert X pads if detector orientation negative (det facing upstream)
	          if( digit1p->GetDetOrientation() < 0 )
		      {
		        pixX = 31 - pixX;
		      }
	          mapped=true;
	          break;
	        }
	    
	      } // Loop over maps 
	  
	      // Calibrations mapped
	      if( mapped ){
	      /*
	        cout << "The source Id is " 
	             << src_id << " " << adc_id << " " << chip_id << " " << endl; 
	        cout << "  wire=" << wire << " ";
	        cout << c.channels[chip_chan].flag << " ";
	        cout << c.channels[chip_chan].pedestal_mean << " ";
	        cout << c.channels[chip_chan].pedestal_sigma << " ";
	        cout << c.channels[chip_chan].calibration_mean << " ";
	        cout << c.channels[chip_chan].calibration_sigma <<endl;
	      */
            gem_chan_cals_ok++;
	    
	        // Create strip channel with calibrations
	        if (stripplane) {
	          stripplane->AddChan(wire, hem, c.channels[chip_chan].flag,
			      	                 c.channels[chip_chan].pedestal_mean,
				                 c.channels[chip_chan].pedestal_sigma,
                                                 c.chip_id, chip_chan);
	        }

	        // Create pixel channel with calibrations
	        else if (pixelplane) {
	          pixelplane->AddChan(wire, pixX, pixY, c.channels[chip_chan].flag,
			  	                        c.channels[chip_chan].pedestal_mean,
			    	                        c.channels[chip_chan].pedestal_sigma,
                                                        c.chip_id, chip_chan);
	        }

	      }
	      else {
	        CsErrLog::msg(elError,__FILE__,__LINE__,
		   	              "%s: Mapping failed for ChipId: %i, ChipChannel: %i!",
                          GetTBName().c_str(), chip_id, chip_chan);
	      }	    
	    } // Loop over chip channels
      } // Loop over chips
    } // Mapping of APV calibrations data

    // Apply default calibrations
    if ( gem_chan_cals_ok != daq_maps.GetWires( DetID(GetTBName()) ) ) {
      CsErrLog::msg(elError,__FILE__,__LINE__,
                    "%s: Using default calibrations!",GetTBName().c_str());
    
      // Delete any existing calibrations
      if (stripplane) {
        if (stripplane->GetNchannel()!=0) {
          stripplane->ClearChannels();
        }
      }
      // Delete any existing calibrations
      if (pixelplane) {
        if (pixelplane->GetNchannel()!=0) {
	      pixelplane->ClearChannels();
        }
      }

      set<uint16> srcIDs;
      daq_maps.GetSrcIDs( DetID(GetTBName()), srcIDs );
      if ( srcIDs.size() != 1 ) {
        CsErrLog::msg(elFatal,__FILE__, __LINE__,
                      "%s: One plane must be connected to only one source ID!",
                      GetTBName().c_str() );
      }
      uint16 srcID = *srcIDs.begin();
    
      // All maps with given data ID
      CS::ChipAPV::DataID start_id (srcID,   0, 0, 0);
      CS::ChipAPV::DataID finish_id(srcID+1, 0, 0, 0);

      const pair<m_it, m_it> m_start  = daq_maps.equal_range(start_id);
      const pair<m_it, m_it> m_finish = daq_maps.equal_range(finish_id);

      // Fill calibrations for all wires of plane
      for (m_it cc=m_start.first; cc!=m_finish.second; cc++) {
        const CS::ChipAPV::Digit *digit1 =
          dynamic_cast<CS::ChipAPV::Digit*>(cc->second);

        if (digit1==NULL) {
          CsErrLog::msg(elWarning, __FILE__, __LINE__,
                        "%s: ChipAPV wrong map!", GetTBName().c_str() );\
          continue;
        }

        if ( GetTBName() != digit1->GetDetID().GetName() )
          continue;

        int flag = 1;
        float ped = 750.;
        float sigma = 4.;
        int wire, hem, PX, PY;
        if (stripplane) {
          try {
            wire = digit1->GetChannel();
            hem  = digit1->GetChanPos();
          }
          catch (...) {
            CsErrLog::msg(elError,__FILE__,__LINE__,
		                  "%s: Could not apply mapping!",GetTBName().c_str());
            continue;
          }
          stripplane->AddChan(wire,hem,flag,ped,sigma);
        }
        else if (pixelplane) {
          const CS::ChipAPV::DigitPixel *digit1p =
            dynamic_cast<const CS::ChipAPV::DigitPixel*>(digit1);

          if (digit1p == NULL) {
            CsErrLog::msg(elWarning,__FILE__,__LINE__,
		                  "%s: ChipAPV wrong map!",GetTBName().c_str());
            continue;
          }

          // Apply mapping
          try {
            wire = digit1p->GetChannel();
          }
          catch (...) {
            CsErrLog::msg(elError,__FILE__,__LINE__,
		                  "%s: Could not apply mapping!",GetTBName().c_str());
            continue;
          }
 
          // Get pixel coordinates
          pair<int,int> xy=pixgem::detch2xy(wire);
	      PX = xy.first;
	      PY = xy.second;

          // Pixel coordinates invalid
          if (PX == -1 && PY == -1) {
		    CsErrLog::msg(elError,__FILE__,__LINE__,
			              "%s: Could not apply mapping!",
			              GetTBName().c_str());
		    continue;
		  }
	      
	      // Invert X pads if detector orientation negative (det facing upstream)
	      if( digit1p->GetDetOrientation() < 0 ) {
            PX = 31 - PX;
	      }

	      pixelplane->AddChan(wire,PX,PY,flag,ped,sigma);
        }
      }
    }
  } // Daq maps found
      
  // No daq maps found
  else {
    CsErrLog::msg(elFatal,__FILE__,__LINE__,
                  "%s: Daq maps not found!",GetTBName().c_str());
    return false;
  }

  // Debugging output
  if (getenv("COND_DB_DEBUG")!=0 || getenv("COND_DB_DEBUG_GEM")!=0) {
    if (stripplane) stripplane->PrintChannels();
    else if (pixelplane) pixelplane->PrintChannels();
  }

  return true;
}

bool CsPixelGEMDetector::mapMCChannelCal() {
  if (!stripplane && !pixelplane) {
    CsErrLog::msg(elFatal,__FILE__,__LINE__,
                  "%s: mapMCChannelCal() was called, but plane objects do not exist!", 
                  GetTBName().c_str());
  
    return false;
  }

  int   flag  =   1;
  float ped   = 750.;
  float sigma =   4.;
  int   hem   =   1;

  if (stripplane) {
    for (int wire=0; wire<nWir_; wire++) {
      stripplane->AddChan(wire, hem,flag,ped,sigma);
      stripplane->AddChan(wire,-hem,flag,ped,sigma);
    }

    // set the standard time calibration, this time calibration is used
    // in MC to simulate the amplitude of the three samples, though either
    // change all occurences or no!!!
    CsGEMTimeCals time_cals; time_cals.Clear();
    CsGEMTimeCalOld *cal0 = new CsGEMTimeCalOld(0.2,1.1,-40.,22.,1.65);
    CsGEMTimeCalOld *cal1 = new CsGEMTimeCalOld(0.3,1.0,0.,27.,1.23);
    time_cals.Set(cal0,cal1);
    stripplane->SetTimeCals(time_cals);
  } else { // pixelplane
    for (int pixU=0; pixU<nWir_; pixU++)
      for (int pixV=0; pixV<nWirV_; pixV++)
        pixelplane->AddChan(pixU + pixV*nWir_, pixU, pixV, flag, ped, sigma);

    // set the standard time calibration, this time calibration is used
    // in MC to simulate the amplitude of the three samples, though either
    // change all occurences or no!!!
    CsGEMTimeCals time_cals; time_cals.Clear();
    CsGEMTimeCalOld *cal0 = new CsGEMTimeCalOld(0.2,1.1,-40.,22.,1.65);
    CsGEMTimeCalOld *cal1 = new CsGEMTimeCalOld(0.3,1.0,0.,27.,1.23);
    time_cals.Set(cal0,cal1);
    pixelplane->SetTimeCals(time_cals);
  }

  return true;
}

