/*!
  \file    doFitting.cc
  \brief   Vertex fitting
  \author  Alexandre Korzenev
  \version $Revision: 14094 $
  \date    $Date: 2015-11-06 16:28:48 +0100 (Fri, 06 Nov 2015) $ 

  \param \e vrts = Reference to list of CsVertex objects to be fitted.
  \param \e specials = Reference to map of scattered-mu ID.
  \param \e reTrackingON: true if expecting to be in the context of a 2nd pass of vertexing, after a re-tracking has taken place.
  \param \e T0: If not a null pointer, on input it holds the T0 used in tracking: if this turns out to deviate too much from eventual best vertex' time, and \e reTrackingON is non true, then doFitting aborts and returns false so that a re-tracking be requested. On output, \e T0 holds the Best Vertex time, which is the recommended reference time to use as T0 in the re-tracking.  

  The method "CsKalmanFitting::doFitting" sets the Best Primary Vertex:
(For the Drell-Yan case (enabled by option "CsKalmanFitting BpVDY"), see infra
method "getBestVertexDY".)
  - Basic definition:
    - largest #tracks first,
   and, in order to raise ambiguities:
    - interacting pVertex wins over non interacting one, w/ interaction being
     decided on view of the momentum transfer to the leading 2ndary track and
     scattering angle,
    - else, best chi2.
    The problematics: select one* true interaction among otherwise fake vertices
   built on a non-interacting (or merely elastically scattered) beam track,
   continued into the spectrometer and possibly picking up some 2ndaries from
   the true interaction.
   (*): The probability for two interactions to be accidentally coincident (w/in
   typical time gate) is negligible in the muon beam case. Can happen in the
   hadron case: whether one could define a criterion to select the best one
   among these two hasn't been considered. 
  - mu' (whether it's attached to the vertex or not) is not considered (except
   special cases, cf. infra). For it can be a mere non-interacting halo track
   continued into the spectrometer (provided an actual interaction is not
   requested, the coincidence rate is sizable) instead of the interaction that
   pulled the trigger.
    (Note: That halo tracks may have supplied the hodo component of the trigger,
   while an interaction outside hodo acceptance (Saleve side, low-Q^2, ...)
   would have supplied the calo component, if required.)
  - Special case I: Can happen that 2 vertices correspond to exactly the same
   set of tracks, when there are 2 mu' candidates. 2 vertices are built: on each
   mu'. One vertex sees its candidate loose its mu'-ID, being rejected from the
   list of tracks (cf. cut on "d12" in "CsAP::FindPrimary"), while being later
   recovered ("CsKF::Refinding"). In that case, let's give precedence to the
   vertex w/ mu'. Btw, it has no real impact on the reconstruction (same event
   time, e.g.). Only that it allows to include this genuine mu/mu' vertex into
   the statistics and histogramming.
  - Special case II: 2 pVertices that turn out to be non interacting.
  - The singling out the best pVertex can become tricky when one is considering
   channels w/ neutrals in the final state. One has then often to arbitrate
   between a non-interacting candidate and a more interesting one, which both
   have the same (small, e.g.: two) # of 2ndary tracks. Examples: DIS, or
   photo-production, w/ neutrals (DVCS, pi0, K0S, etc...), diffractive w/ FS =
   diffracted pi + neutrals (e.g., pi eta eta, pi K0S K0S). It is then important
   to select it as best vertex already at the reconstruction time (as opposed to
   physics analysis time), in order to have the correct event time. Event
   time is (can be) set by the best vertex' beam track time: it then serves as a
   reference for cutting on ECAL cell time or for tracking pi's from K0S decay
   in drift detectors. (Note that the present selection of best vertex is not
   necessarily ported to PHAST, which may set its own Best Primary Vertex
   according to its own rules.)
*/

#include <cassert>
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "CsInit.h"
#include "CsOpt.h"
#include "CsStopwatch.h"
#include "CsErrLog.h"
#include "CsEvent.h"
#include "CsGeant3.h"
#include "DaqDataDecoding/DaqOption.h"
#include "CsKalmanFitting.h"
#include "CsBeam.h"
#include "CsGeom.h"

using namespace std;
using namespace CLHEP;

extern "C" float prob_(float*,int*);
extern double Q2(double,double,double);
extern double xbj(double,double,double);
int ReFinding(TMtx &Xn,TMtx &Cn, list<CsVTrack>& vtrks, list<CsTrack*>& trks);

const double M_p   = 0.9382723,    M2_p  = M_p*M_p;
const double M_mu  = 0.1056583568, M2_mu = M_mu*M_mu;
const double M_pi  = 0.13957018,   M2_pi = M_pi*M_pi; 

bool CsKalmanFitting::doFitting(list<CsVertex*> &vrts,
				map<CsTrack*,bool> &specials,
				bool reTrackingON, double *T0)
{
  static CsStopwatch _chronometer(1);
  static int _chrono = _chronometer.start();
  static int _nfittingTime = 0; static double _totTime = 0;
  double dt = _chronometer.inter(_chrono);

  double chi2; TMtx Cn(3,3), Xn(3); CsTrack *trk;
  list<CsTrack*> trks;  list<CsTrack*>::iterator it;
  list<CsVTrack> vtrks; list<CsVTrack>::iterator ivt;
  list<CsVertex*> vrt_to_erase;

  // STREAMER
  static unsigned int doStream = 0; static int iStream;
  const unsigned int pStream = 0x1, sStream = 0x2, phStream = 0x4;

  static double zT, zTMn, zTMx, zTR2;// Target centre, min, max, radius^2 (mm)
  bool hadronRun = CsInit::Instance()->IsAHadronJob();

  //    ******************** HISTOGRAMS ********************
  static TH1F     *OM_xbj,*OM_Q2,*OM_xbjQ1;
  static TH1F     *OM_xbjMT,*OM_xbjLT,*OM_xbjOT,*OM_xbjCT,*OM_xbjLAST,*OM_xbjIMT;
  static CsHist1F *OM_ybjMT,*OM_ybjLT,*OM_ybjOT,*OM_ybjCT,*OM_ybjLAST,*OM_ybjIMT;
  static bool first = true;
  static TH1D *hBeamT, *hTrkT, *hTrkTdT, *hBeamP;
  static TH1D *hVrtChi2,*hVrtProb, *hTrkVChi2, *hTrkVProb, *hpiNChi2;
  static TH1D *hVrtX, *hVrtY, *hVrtZ;
  static CsHist1F *hPullX,*hPullY,*hPullZ;
  static CsHist1F *hSigmaX,*hSigmaY,*hSigmaZ;
  static CsHist2F *hSigmaX_N,*hSigmaY_N,*hSigmaZ_N;
  static CsHist2F *hSigmaZ_X,*hSigmaZ_Y,*hSigmaZ_Z;
  static CsHist2F *hProb_T,*hProb_Mom;
  static TH2D *hKZ_T, *hKZ_Q2, *hKZ_Q21h;
  static TH2D *hTvsVProb;
  static TH1D *hNTrk,  *hNTrk_SI,  *hNTrk_C,  *hNTrk_O;
  static TH1D *hNTrkM, *hNTrkM_SI, *hNTrkM_C, *hNTrkM_O;
  static TH1D *hNTrkG, *hNTrkG_SI, *hNTrkG_C, *hNTrkG_O;
  static TH1D *hNVtx,  *hNVtx_SI,  *hNVtx_C,  *hNVtx_O;
  static TH1D *hNVtxM, *hNVtxM_SI, *hNVtxM_C, *hNVtxM_O;
  static TH1D *hNVtxG, *hNVtxG_SI, *hNVtxG_C, *hNVtxG_O;
  static TH2D *hBVQ2,  *hBVQ2_SI,  *hBVQ2_C,  *hBVQ2_O;
  static TH2D *hBVQ2h, *hBVQ2h_SI, *hBVQ2h_C, *hBVQ2h_O;
  //     ********** SINGLE OUT SPECIAL TYPES of TRIGGERS **********
  // - Special types are:
  //   - "stdTrigs" = standard triggers.
  //   - "specialTrigs", to be inclusively histo'd,
  //   - "pureTrigs", to be histo'd only when pure (typically their overlap w/
  //    the rest of the triggers is large).
  //  These 3 types allow to single out interesting events for histo'ing,
  //  whatever the data taking.
  // - Muon data taking: "specialTrigs" is assigned OT. "pureTrigs" is assigned
  //  pure CT (w/ allowed combination w/ highQ2T). And "stdTrigs" are IMLT.
  // - Hadron data taking: "specialTrigs" is assigned NIT, i.e. non interacting
  //  (interesting for it allows to monitor the reconstruction of paraxial
  //  tracks, much akin to that of the scattered hadron, but which is
  //  unfortunately not enabled in physics data taking). "pureTrigs" are
  //  beam-based triggers, BT|aBT (which, w/ their large overlap w/ the rest of
  //  the triggers, are interesting to consider only when pure). And "stdTrigs"
  //  are RPD-based triggers.
  // - DVCS running: "specialTrigs" is assigned "DVCS" in 2009, and "Tiger" in 
  //  2012. And "stdTrigs" are MLO[L]T. "pureTrigs" are BT|aBT.
  // - Primakoff: "stdTrigs" is assigned Primakoff1|2. "pureTrigs" are BT|aBT.
  //  "specialTrigs" is RPD.
  // - Low-t: "stdTrigs" is assigned lowt1|2. "pureTrigs" are BT|aBT.
  //  "specialTrigs" is RPD.
  // - DY: "stdTrigs" is assigned MTLT|OTLT|LAST2, "pureTrigs" is CT,
  //  "specialTrigs" is "MT|OT|LT".
  static unsigned int pureTrigs, specialTrigs , stdTrigs = 0x7;

  if (hist_ && first) {
    const double Zint = 5000; double SL, SZ;
    if (CsEvent::Instance()->isAMonteCarloEvent()) { SL = 0.20; SZ = 150; }
    else                                           { SL = 2;    SZ = 400; }

    bool hadronRun = CsInit::Instance()->IsAHadronJob();
    bool primakoffRun, lowtRun = false, dyRun; int  dvcsRun;

    const unsigned int caloTrig = 0x10;
    unsigned int caloBasedTrigs = caloTrig, beamTrigs = 0;
    const CS::DaqOption &opts = CsInit::Instance()->getDaqEventsManager().GetDaqOptions();
    int bit = -1;// Bit 0x200 doesn't necessarily mean high-Q2 trigger: check it
    try { bit = opts.GetTriggerBit("LargeQ2Trigger"); } catch ( ... ) { }
    if (bit>=0) caloBasedTrigs |= 1<<bit;
    bit = -1;
    try { bit = opts.GetTriggerBit("BeamTrigger"); } catch ( ... ) { }
    if (bit>=0) beamTrigs |= 1<<bit;
    bit = -1;
    try { bit = opts.GetTriggerBit("altBeamTrigger"); } catch ( ... ) { }
    if (bit>=0) beamTrigs |= 1<<bit;
    pureTrigs =  hadronRun ? beamTrigs : caloBasedTrigs;
    unsigned int niTrig =// As of 2008/12, NIT's not described in the mapping...
      0x10; // ...=> Built-in (although NIT is not present throughout)
    specialTrigs = hadronRun ? niTrig : 0x8;
    // ***** DVCS data taking
    bit = -1;
    try { bit = opts.GetTriggerBit("DVCS"); } catch ( ... ) { }
    if (bit>=0) {
      specialTrigs = 1<<bit; pureTrigs = beamTrigs; stdTrigs = 0xe;
      dvcsRun = 2009;
    }
    else {
      dvcsRun = 0;
    }
    try { bit = opts.GetTriggerBit("Tiger"); } catch ( ... ) { }
    if (bit>=0) {
      specialTrigs = 1<<bit; pureTrigs = beamTrigs; stdTrigs = 0x20e;
      dvcsRun = 2012;
    }
    // ***** DY data taking
    bit = -1;
    try { bit = opts.GetTriggerBit("LAST2"); } catch ( ... ) { }
    if (bit>=0) {
      specialTrigs = 0x20a; pureTrigs = 0x10; stdTrigs = 0x105;
      dyRun = true;
    }
    else {
      dyRun = false;
    }
    // ***** Primakoff data taking
    bit = -1;
    try { bit = opts.GetTriggerBit("Primakoff1"); } catch ( ... ) { }
    if (bit<0) {
      try { bit = opts.GetTriggerBit("Prim1"); } catch ( ... ) { }
    }
    if (bit>=0) {
      stdTrigs = 1<<bit; primakoffRun = true;
    }
    else primakoffRun = false;
    bit = -1;
    try { bit = opts.GetTriggerBit("Primakoff2"); } catch ( ... ) { }
    if (bit<0) {
      try { bit = opts.GetTriggerBit("Prim2"); } catch ( ... ) { }
    }
    if (bit>=0) {
      if (primakoffRun) stdTrigs |= 1<<bit;
      else              stdTrigs = 1<<bit;
      primakoffRun = true;
    }
    if (primakoffRun) {
      pureTrigs = beamTrigs; specialTrigs = 0x1;
    }
    else {
      // ***** Low-t data taking
      // (Note: This is conditioned by !=Primakoff, because lowt triggers are
      // still described in the mapping (DAQ.xml) for the 2009 Primakoff data
      // taking, although they, probably, were not active.)
      bit = -1;
      try { bit = opts.GetTriggerBit("lowt1"); } catch ( ... ) { }
      if (bit>=0) {
	stdTrigs = 1<<bit; lowtRun = true;
      }
      else lowtRun = false;
      bit = -1;
      try { bit = opts.GetTriggerBit("lowt2"); } catch ( ... ) { }
      if (bit>=0) {
	if (lowtRun) stdTrigs |= 1<<bit;
	else         stdTrigs = 1<<bit;
	lowtRun = true;
      }
      if (lowtRun) {
	pureTrigs = beamTrigs; specialTrigs = 0x1;
      }
    }

    CsHistograms::SetCurrentPath("/CsKalmanFitting/");
    hBeamP = new TH1D("hBeamP","Beam momentum",160,BeamP0_-40,BeamP0_+40);
    hVrtChi2 = new TH1D("hVrtChi2","Best pVertex #chi^{2}/NDF", 160,0,20);
    hVrtProb = new TH1D("hVrtProb","Best pVertex #chi^{2}Proba",100,0,1);
    hTrkVChi2 = new TH1D("hTrkVChi2","Track #chi^{2}/NDF in fit to pVertex",
			 160,0,20);
    hTrkVProb = new TH1D("hTrkVProb","Track #chi^{2}Proba in fit to pVertex",
			 100,0,1);
    hpiNChi2 = new TH1D("hpiNChi2","Track #chi^{2}/NDF for N#pi excl. events",
			100,0,20);

    zT = CsGeom::Instance()->getTargetCenter();
    double Q2bins[] = {.001,.01,.1,1,10,100};
    int nQ2bins = sizeof(Q2bins)/sizeof(double)-1; 
    hKZ_Q2   = new TH2D("KZ_Q2",  "Q2 vs. Z",     250,zT-Zint,zT+Zint,
			nQ2bins,Q2bins);
    hKZ_Q21h = new TH2D("KZ_Q21h","Q2 vs. Z - 1H",250,zT-Zint,zT+Zint,
			nQ2bins,Q2bins);

    char hN[]   = "hNBwBMS_SI";
    //char hN[] = "hBVQ2h_SI";
    char hT[] =
      "Q^{2}+h: all,.05<y<.95,-100<Z<+100/R<1.6,BMS,#chi^{2}OK - PrimakoffT ";
    //"#Tracks/#mu#mu'Vtx w/ BMS - PrimakoffT ";
    const char *tabs[4] = {"","_SI","_C","_O"};
    const char *tags[7][4]= { {""," - ILMT",      " - CT",   " - OT"},
			      {""," - RPDT",      " - (a)BT"," - NIT"},
			      {""," - LMOT",      " - (a)BT"," - DVCST"},
			      {""," - PrimakoffT"," - (a)BT"," - DT0"},
			      {""," - lowtT",     " - (a)BT"," - DT0"},
			      {""," - LMOLT",     " - (a)BT"," - TigerT"},
			      {""," - 2muT",      " - CT",   " - 1muT"} };
    int i = 0; if (hadronRun) i = 1;
    if (dvcsRun==2009) i = 2; if (dvcsRun==2012) i = 5; if (dyRun) i = 6;
    if (primakoffRun) i = 3; if (lowtRun) i = 4;
    TH1D **hs[4][6] =
      { {&hNTrk,   &hNTrkM,   &hNTrkG,   &hNVtx,   &hNVtxM,   &hNVtxG},
	{&hNTrk_SI,&hNTrkM_SI,&hNTrkG_SI,&hNVtx_SI,&hNVtxM_SI,&hNVtxG_SI},
	{&hNTrk_C, &hNTrkM_C, &hNTrkG_C, &hNVtx_C, &hNVtxM_C, &hNVtxG_C},
	{&hNTrk_O, &hNTrkM_O, &hNTrkG_O, &hNVtx_O, &hNVtxM_O, &hNVtxG_O} };
    TH2D **h2s[4][2] = { {&hBVQ2,  &hBVQ2h},  {&hBVQ2_SI,&hBVQ2h_SI},
			 {&hBVQ2_C,&hBVQ2h_C},{&hBVQ2_O, &hBVQ2h_O} };
    // Determine the target zone extension as mid way between target ends and
    // closest upstream (resp. downstream) tracking detector. This in order to
    // histogram the vertexing efficiency:
    // - w/in the zone of highest interest, viz. the target's,
    // - excepting edge effects, hence away from w/ target ends, so as to be
    //  independent of vertex resolution,
    // - excepting the intricacies associated w/ vertices orginating in tracking
    //  detectors.
    const vector<CsTargetCell> &cells = CsGeom::Instance()->getTargetCells();
    double zTR; int nCells = (int)cells.size(); if (nCells) {
      int ic; for (ic = 0, zTMn=zTMx = zT, zTR = 0; ic<nCells; ic++) {
	double z = cells[ic].x[0], l = cells[ic].length;
	if (z-l/2<zTMn) zTMn = z-l/2; if (z+l/2>zTMx) zTMx = z+l/2;
	if (cells[ic].radius>zTR) zTR = cells[ic].radius;
      }
      const list<CsDetector*> &dets = CsGeom::Instance()->getDetectors();
      list<CsDetector*>::const_iterator id;
      double zDUp, zDDown; for (id = dets.begin(), zDUp = zTMn, zDDown = zTMx;
			       id!=dets.end(); id++) {
	double zD = (*id)->getZcm(); if (zD<zTMn) zDUp = zD;
	else if (zD>zTMx) {                       zDDown = zD; break; }
      }
      zTMn = (zDUp+zTMn)/2; zTMx = (zTMx+zDDown)/2; 
    }
    else { // Default, corresponding the [2002,2004] case.
      zTMn = zT-1000; zTMx = zT+1000; zTR = 15;
    }
    zTMn = 10.*int(zTMn/10+.5); zTMx = 10.*int(zTMx/10+.5);
    zTR = int(zTR*1.33+.5); zTR2 = zTR*zTR;
    for (int j = 0; j<4; j++) {
      // # of tracks per pVertex and # of pVertices
      sprintf(hN,"hNTrk%s", tabs[j]);
      sprintf(hT,"#Tracks/pVertex%s",          tags[i][j]);
      *hs[j][0] = new TH1D(hN,hT,15,0,15);
      sprintf(hN,"hNTrkM%s",tabs[j]);
      sprintf(hT,"#Tracks/#mu#mu'Vtx%s",       tags[i][j]);
      *hs[j][1] = new TH1D(hN,hT,15,0,15);
      sprintf(hN,"hNTrkG%s",tabs[j]);
      sprintf(hT,"#Tracks/#mu#mu'Vtx w/ BMS%s",tags[i][j]);
      *hs[j][2] = new TH1D(hN,hT,15,0,15);
      sprintf(hN,"hNVtx%s", tabs[j]);
      sprintf(hT,"#pVertices%s",               tags[i][j]);
      *hs[j][3] = new TH1D(hN,hT, 8,0, 8);
      sprintf(hN,"hNVtxM%s",tabs[j]);
      sprintf(hT,"##mu#mu'Vertices%s",         tags[i][j]);
      *hs[j][4] = new TH1D(hN,hT, 8,0, 8);
      sprintf(hN,"hNVtxG%s",tabs[j]);
      sprintf(hT,"##mu#mu'Vertices w/ BMS%s",  tags[i][j]);
      *hs[j][5] = new TH1D(hN,hT, 8,0, 8);
      sprintf(hN,"hBVQ2%s", tabs[j]);
      sprintf(hT,"Q^{2}: all,.05<y<.95,%.0f<Z<%.0f/R<%.1f,BMS,#chi^{2}OK%s",
	      zTMn/10,zTMx/10,zTR/10,          tags[i][j]);
      *h2s[j][0] = new TH2D(hN,hT,nQ2bins,Q2bins,5,-.5,4.5);
      sprintf(hN,"hBVQ2h%s", tabs[j]);
      sprintf(hT,"Q^{2}+h: all,.05<y<.95,%.0f<Z<%.0f/R<%.1f,BMS,#chi^{2}OK%s",
	      zTMn/10,zTMx/10,zTR/10,          tags[i][j]);
      *h2s[j][1] = new TH2D(hN,hT,nQ2bins,Q2bins,5,-.5,4.5);
    }

    hProb_T = new CsHist2F("P_T","Probability vs Theta", 50,0,1, 50,0,0.03);
    hProb_Mom = new CsHist2F("P_M","Probability vs Momentum", 50,0,1, 100,0,160);
    hTvsVProb = new TH2D("hTvsVProb","Track #chi^{2}Proba in fit to pV vs. Vertex #chi^{2}Proba",50,0,1,50,0,1);

    hBeamT  = new TH1D("hBeamT", "Best pVertex Beam Time",           60,-15,15);
    hTrkT   = new TH1D("hTrkT",  "TrackT-beamT in best pVertex",     50,-25,25);
    hTrkTdT = new TH1D("hTrkTdT","(TrackT-beamT)/dT in best pVertex",50,-5,5);

    CsHistograms::SetCurrentPath("/CsKalmanFitting/Delta");
    hVrtX = new TH1D("hVrtX","Best pVertex X",100,-40,40);
    hVrtY = new TH1D("hVrtY","Best pVertex Y",100,-40,40);
    hVrtZ = new TH1D("hVrtZ","Best pvertex Z",500,zT-Zint/2,zT+Zint);

    hSigmaX = new CsHist1F("KSX","Sigma X",100,0,SL);
    hSigmaY = new CsHist1F("KSY","Sigma Y",100,0,SL);
    hSigmaZ = new CsHist1F("KSZ","Sigma Z",500,0,SZ);

    hSigmaX_N = new CsHist2F("KSX_N","Sigma X vs N tracks", 10,0,10, 100,0,SL);
    hSigmaY_N = new CsHist2F("KSY_N","Sigma Y vs N tracks", 10,0,10, 100,0,SL);
    hSigmaZ_N = new CsHist2F("KSZ_N","Sigma Z vs N tracks", 10,0,10, 100,0,SZ);

    hSigmaZ_X = new CsHist2F("KZ_X","Z vs Sigma X", 250,zT-Zint,zT+Zint, 200,0,SL );
    hSigmaZ_Y = new CsHist2F("KZ_Y","Z vs Sigma Y", 250,zT-Zint,zT+Zint, 200,0,SL );
    hSigmaZ_Z = new CsHist2F("KZ_Z","Z vs Sigma Z", 550,zT-Zint,zT+Zint, 200,0,SZ );

    hKZ_T  = new TH2D("KZ_T","#thetaMX vs. Z",250,zT-Zint,zT+Zint,50,0,0.05);

    CsHistograms::SetCurrentPath("/CsKalmanFitting/PullsSm");
    hPullX = new CsHist1F("hPullvX","Pulls: X",100,-10,10);
    hPullY = new CsHist1F("hPullvY","Pulls: Y",100,-10,10);
    hPullZ = new CsHist1F("hPullvZ","Pulls: Z",100,-10,10);

    CsHistograms::SetCurrentPath("/CsKalmanFitting/OnlineMonitor");
    OM_Zvtx  = new CsHist1F("OM_Zvtx"  ,"OM - Zvtx distribution", 100, -1000., 1000.);
    OM_ZvsX  = new CsHist2F("OM_ZvsX"  ,"OM - Xvtx vs Zvtx",       50, -1000., 1000., 50, -50., 50.);
    OM_ZvsY  = new CsHist2F("OM_ZvsY"  ,"OM - Yvtx vs Zvtx",       50, -1000., 1000., 50, -50., 50.);
    OM_pipi  = new CsHist1F("OM_pipi"  ,"OM - K invariant mass",50,0.45,0.55);
    OM_NTrkM = new CsHist1F("OM_NTrkM" ,"OM - #Tracks/#mu#mu'BestV",15,0,15);
    OM_BeamP = new CsHist1F("OM_BeamP" ,"OM - beam momentum  ",80,120,200);
    OM_BeamPz= new CsHist1F("OM_BeamPz","OM - beam P_{z}     ",80,120,200);
    OM_BeamPx= new CsHist1F("OM_BeamPx","OM - beam P_{x}     ",50, -1,  1);
    OM_BeamPy= new CsHist1F("OM_BeamPy","OM - beam P_{y}     ",50, -1,  1);
    OM_BeamTx= new CsHist1F("OM_BeamTx","OM - beam #theta_{x}",50,-.005,0.005);
    OM_BeamTy= new CsHist1F("OM_BeamTy","OM - beam #theta_{y}",50,-.005,0.005);
    OM_ybj   = new CsHist1F("OM_ybj"   ,"OM - y all triggers ",50,0,1);
    OM_zhad  = new CsHist1F("OM_zhad"  ,"OM - z hadrons      ",50,0,1);
    OM_mult  = new CsHist1F("OM_mult"  ,"OM - had multiplicity",20,0,20);


    const int nbinx=100;
    const int ndimx=101;
    const int nbiny=100;
    const int ndimy=101;
    
    double xmin=1.e-4;
    double xmax=1.   ;
    float  xbins[ndimx];
    
    double ymin=1.e-1;
    double ymax=1.e+2;
    float  ybins[ndimy];
    
    double xstep=(log10(xmax)-log10(xmin))/nbinx;
    double ystep=(log10(ymax)-log10(ymin))/nbiny;
    
    for(int ii=0;ii<=nbinx;ii++){
      xbins[ii]=float(pow(10.,log10(xmin)+xstep*ii));
      ybins[ii]=float(pow(10.,log10(ymin)+ystep*ii));
    }
    
    
    OM_xbj     = new     TH1F("OM_xbj"    ,"OM - xbj all triggers"  ,nbinx,xbins);
    OM_xbjMT   = new     TH1F("OM_xbjMT"  ,"OM - xbj middle trg  "  ,nbinx,xbins);
    OM_xbjLT   = new     TH1F("OM_xbjLT"  ,"OM - xbj ladder trg  "  ,nbinx,xbins);
    OM_xbjOT   = new     TH1F("OM_xbjOT"  ,"OM - xbj outer  trg  "  ,nbinx,xbins);
    OM_xbjCT   = new     TH1F("OM_xbjCT"  ,"OM - xbj calo   trg  "  ,nbinx,xbins);
    OM_xbjLAST = new     TH1F("OM_xbjLAST","OM - xbj LAST   trg  "  ,nbinx,xbins);
    OM_xbjIMT  = new     TH1F("OM_xbjIMT" ,"OM - xbj incM   trg  "  ,nbinx,xbins);
    OM_xbjQ1   = new     TH1F("OM_xbjQ1"  ,"OM - xbj AT Q^{2}>1  "  ,nbinx,xbins);
    OM_Q2      = new     TH1F("OM_Q2"     ,"OM - Q^{2} all triggers",nbiny,ybins);
    OM_ybjMT   = new CsHist1F("OM_ybjMT"  ,"OM - y middle trg   ",50,0,1);
    OM_ybjLT   = new CsHist1F("OM_ybjLT"  ,"OM - y ladder trg   ",50,0,1);
    OM_ybjOT   = new CsHist1F("OM_ybjOT"  ,"OM - y outer  trg   ",50,0,1);
    OM_ybjCT   = new CsHist1F("OM_ybjCT"  ,"OM - y calo   trg   ",50,0,1);
    OM_ybjLAST = new CsHist1F("OM_ybjLAST","OM - y LAST   trg   ",50,0,1);
    OM_ybjIMT  = new CsHist1F("OM_ybjIMT" ,"OM - y incM   trg   ",50,0,1);

    CsHistograms::SetCurrentPath("/");
  }

  if (first) {
    list<string> streams;        // ********** OUTPUT STREAMS **********
    // Syntax: CsKalmanfitting WriteBack <file_name> [Primary] [Secondary] [Physics]
    if (CsOpt::Instance()->getOpt("CsKalmanFitting","WriteBack",streams)) {
      list<string>::iterator is = streams.begin();
      CsEvent *event = CsEvent::Instance(); iStream = -1;
      //            ***** FIRST FIELD = FILE NAME *****
      if (event->getNOutStreams()) {// If a write-back stream already opened,...
	// ...we would like to have the possibility to use it also for the
	// present CsKalmanFitting write-back of interesting events. To that
	// end, we have to check that the file name associated w/ this other
	// stream is the same as the one we are having here. Unfortunately I
	// don't know how to retrieve this file name. The best I can do is to
	// consider the case of a mass production job where we have 2 kinds of
	// write-back: CsKalmanFitting's and CsEvent's and the opening CsEvent's
	// one must already have taken place, given coral's flowchart. =>
	//  - Let's check that a "CsEvent WriteBack" has been requested, w/ same
	//   associated file name.
	//  - Let's require that there is a unique stream opened so far.
	string csEvtFile; if (CsOpt::Instance()->getOpt("CsEvent","WriteBack",
							csEvtFile)) {
	  if (*is==csEvtFile) {
	    if (event->getNOutStreams()>1)
	      CsErrLog::msg(elFatal,__FILE__,__LINE__,"3 write-back requested (\"CsKalmanFitting %s\" + \"CsEvent %s\" + else) => don't know how to sort them out!",(*is).c_str(),csEvtFile.c_str());
	    iStream = 0;
	  }
	}
      }
      if (iStream<0) {
	if ((iStream = CsEvent::Instance()->openRawEventOutputStream(*is))<0)
	  CsErrLog::msg(elFatal,__FILE__,__LINE__,
			"Cannot open WriteBack stream \"%s\"",(*is).c_str());
      }
      is++; // ***** NEXT FIELDS: SPECIF OF KINDS of EVENTS TO WRITE-BACK *****
      for (; is!=streams.end(); is++) {
	if      (*is=="Primary")   doStream |= pStream;
	else if (*is=="Secondary") doStream |= sStream;
	else if (*is=="Physics")   doStream |= phStream;
	else
	  CsErrLog::msg(elFatal,__FILE__,__LINE__,
			"Unknown WriteBack selection \"%s\"",(*is).c_str());
      }
      if (doStream==0)
	CsErrLog::mes(elFatal,"WriteBack: No selection specified!");
    }

    first = false;
  }            // ***** END of INIT BLOCK *****

  CsEvent *cse = CsEvent::Instance();// Not const: cf. "outputRawEventToStream".
  const unsigned int allTrigs = 0xffff; // Although TCS can handle 0x7fffff, it was decided to limit max. #triggers to 16 and reserve for an other, not related to trigger, use the MS bits of the trigger TDC module.
  unsigned int evTrig = cse->getTriggerMask();
  evTrig &= allTrigs;        // Cut away trailing end bits, i.e. online filter

  BestPVertex_ = 0; // CsKalmanFitting's pointer to Best Primary Vertex

  list<CsVertex*>::iterator iv;

  static unsigned int prvENum = 0, eNum; // Statistics: # of calls to fitting
  if ((eNum = cse->getEventNumberInRun())!=prvENum) {// Avoid double counting...
    // ...those events where the 1st pass of vertex fitting requests a
    // re-tracking (cf. infra), and then a 2nd pass of vertex fitting.
    prvENum = eNum; 
    for (iv = vrts.begin(); iv!=vrts.end(); iv++)
      if ((*iv)->isPrimary()) {// At least one pVertex identified in the PR step
	statistics_[0]++; break;                   // for statistics
      }
  }
  int nPVrts, nPVrtsM /* mu/mu' */, nPVrtsG /* Golden: mu/mu' w/ BMS */, jv;
  // Best pVertex = largest #tracks. If ambiguity, best chi2 and whether the
  // vertex is an interacting one (criteria being significant momentum tranfer
  // to the target, scattering angle, and momentum transfer to leading 2ndary
  // track).
  CsVertex *bestV; int nonInteractingBV;
  static bool hasMupBV, hasBMSBV; static int piNBV;
  static int nGoodTrksBV, nMuswTBV; // Caracteristics of BpV in the DY case
  static double tChi2BV;
  unsigned int rejection;
  for (iv = vrts.begin(), nPVrts=nPVrtsM=nPVrtsG=jv = 0, bestV = 0,
	 rejection = 0, nonInteractingBV = 0; iv!=vrts.end(); iv++) {

    // ********************************************************************
    // ******************** LOOP on CANDIDATE VERTICES ********************
    // ********************************************************************

    CsVertex *vrt = *iv;
    if (Print_[0]) printf("\nCsKF::doFitting==> ***** VERTEX %d *****\n",jv++);

    // Coral -> Traffic
    Xn(1) = vrt->getZ()/10; Xn(2) = vrt->getX()/10; Xn(3) = vrt->getY()/10;
    HepMatrix Ctmp = *(vrt->getCov(0));
    TMtx X0(Xn); HepMatrix C0(Ctmp); // Remember X0,C0, in view of pseudo-Pulls

    trks = vrt->getTracks();
    
    // ********** KALMAN FIT AND FILTER **********
    int NIter = 0, NTrkReFind = 0, NTrkReFindIn = 0;
#ifdef doFitting_RESTRICTED_RESCUE
    // Disabled by default ('cause not significantly effective)
    double chi2CCut_bak = Chi2VCut_;  // Remember initial v chi2 cut
#endif
    do {                         // ********* LOOP on 2 ITERATIONS **********
      if (Print_[0])
	cout<<"CsKF::doFitting==> ***** ITERATION #="<<NIter<<" *****\n\n";
      Cn(1,1) = Ctmp[2][2]/100; Cn(1,2) = Ctmp[2][0]/100; Cn(1,3) = Ctmp[2][1]/100;
      Cn(2,1) = Ctmp[0][2]/100; Cn(2,2) = Ctmp[0][0]/100; Cn(2,3) = Ctmp[0][1]/100;
      Cn(3,1) = Ctmp[1][2]/100; Cn(3,2) = Ctmp[1][0]/100; Cn(3,3) = Ctmp[1][1]/100;
      // Built list of CsVTracks "vtrks" according to Z of "Xn"
      vtrks.clear(); setTrks(Xn,trks,vrt->isPrimary(),specials,vtrks);

      Kalman(Xn,Cn,vtrks,chi2);    // ********** VERTEX FIT **********

      // Special case of mu w/o BMS: restore sigma(1/p)^2 of beam track.
      // (In that case, p and sigma(1/p) have been assigned values taken
      // directly from option, w/ sigma(1/p) large (expressing the large spread
      // in p of the mu beam and hence the large p uncertainty arising when p is
      // not actually measured). In order for this large uncertainty not to
      // affect the vertex fit (by allowing beam's p to vary greatly), and also
      // in order to leave the p value, which is expected to be a round number
      // and hence a good signature of the case, unchanged, signa(1/p) was set
      // at a low value inside "CskalmanFitting::Kalman".)
      CsVTrack &vBeam = vtrks.front();
      const CsTrack *beam = vBeam.getAssociatedTrk();
      const list<CsZone*> &zones = beam->getZones();
      list<CsZone*>::const_iterator iz; int hasBMS;
      for (iz = zones.begin(), hasBMS = 0; iz!=zones.end(); iz++)
	if ((*iz)->getName()=="BMSzone") { hasBMS = 1; break; }
      if (!hasBMS) {
	const CsHelix &trackingHelix = beam->getHelices().front();
	double r = trackingHelix.getCop(), sr2 =  trackingHelix.getCov()[14];
	if (fabs(sqrt(sr2)/r)>.02) {// E.g. P=160+/-15/sqrt(12) GeV yields ~2.7%
	  vBeam.Dnk(3,3) = sr2;
	  for (int k = 1; k<=2; k++) vBeam.Dnk(k,3)=vBeam.Dnk(3,k) = 0;
	  for (int k = 1; k<=3; k++) vBeam.Enk(k,3) = 0;
	}
      }

      if (NIter==0) {                           // ***** RESCUE...*****
	NTrkReFindIn = vtrks.size();
	NTrkReFind = ReFinding( Xn,Cn, vtrks,trks );
#ifdef doFitting_RESTRICTED_RESCUE
	// Accept rescued tracks only if don't degrade chi2
	if (chi2*2<Chi2VCut_) Chi2VCut_ = chi2*2;
#endif
      }
      else if (NTrkReFind) { // What's the use of this?!
	trks.clear();
	for( ivt=vtrks.begin(); ivt!=vtrks.end(); ivt++ ) 
	  trks.push_back( (*ivt).getAssociatedTrk() );
      }
      NIter++;
    } while (NTrkReFind && NIter<2);
#ifdef doFitting_RESTRICTED_RESCUE
    Chi2VCut_ = chi2CCut_bak; // Restore initial v chi2 cut
#endif

    if (vtrks.size()<2) {                     // ***** REQUIRE AT LEAST 2 TRACKS
      vrt_to_erase.push_back(vrt); if (vrt->isPrimary()) rejection |= 0x1;
      continue;
    }

    if (chi2/(2*vtrks.size()-3)>Chi2VCut_) {            // ***** VERTEX CHI2 CUT
      if( vrt->isPrimary() ) {
	CsErrLog::msg(elError,__FILE__,__LINE__,
		      "Primary vertex (Nt=%d) chi2=%.1f > %.1f => Removed.",
		      vtrks.size(),chi2/(2*vtrks.size()-3),Chi2VCut_);
        rejection |= 0x2;
      }
      vrt_to_erase.push_back( vrt ); continue;
    }

    ELossCorrection( Xn, vtrks );           // ***** CORRECTIONS FOR ENERGY LOSS

    //                                 ***** UPDATE LIST OF VERTICES IN CsTracks
    for (it = trks.begin(); it!=trks.end(); it++) (*it)->addVertex( vrt );

    //                                                          ***** SET VERTEX
    vrt->setVertex(10*Xn(2),10*Xn(3),10*Xn(1)); vrt->setChi2( chi2 );

    //                                      ***** SET TRACK PARAMETERS AT VERTEX
    for( ivt=vtrks.begin(); ivt!=vtrks.end(); ivt++ ) {
      Cs3Vector vec( (*ivt).Qnk(1), (*ivt).Qnk(2), (*ivt).Qnk(3) );
      trk = (*ivt).getAssociatedTrk();
      vrt->addTrackAtVertex( trk, vec );
    }


    vrt->clearCov();    //                 ***** SET COVARIANCE MATRIX OF VERTEX
    // As an example: case of a vertex of 2 tracks...
    //       x y z   dxdz dydz q/p   dxdz dydz q/p  
    // x 
    // y       0           3               4
    // z
    //
    // dxdz
    // dydz                1               5
    // q/p
    //
    // dxdz
    // dydz                                2
    // q/p
   if( Covar_ >= 1 ) {
      // TRAFFIC -> CORAL 
      Ctmp[0][0] = Cn(2,2)*100;  Ctmp[0][1] = Cn(2,3)*100;   Ctmp[0][2] = Cn(2,1)*100;
      Ctmp[1][0] = Cn(3,2)*100;  Ctmp[1][1] = Cn(3,3)*100;   Ctmp[1][2] = Cn(3,1)*100;
      Ctmp[2][0] = Cn(1,2)*100;  Ctmp[2][1] = Cn(1,3)*100;   Ctmp[2][2] = Cn(1,1)*100;
      vrt->addCov( new HepMatrix( Ctmp ) );
    }

    // SET COVAR MATRICES OF MOMENTA OF TRACKS IN VERTEX 
    if( Covar_ >= 2 ) {
      for( ivt=vtrks.begin(); ivt!=vtrks.end(); ivt++ ) {
        HepMatrix &hm = (* new HepMatrix(3,3) );
        hm[0][0] = (*ivt).Dnk(1,1); hm[0][1] = (*ivt).Dnk(1,2); hm[0][2] = (*ivt).Dnk(1,3);
        hm[1][0] = (*ivt).Dnk(2,1); hm[1][1] = (*ivt).Dnk(2,2); hm[1][2] = (*ivt).Dnk(2,3);
        hm[2][0] = (*ivt).Dnk(3,1); hm[2][1] = (*ivt).Dnk(3,2); hm[2][2] = (*ivt).Dnk(3,3);
        vrt->addCov( &hm );
      }
    }

    // SET COORDINATE-MOMENTUM CORRELATIONS OF TRACKS IN VERTEX 
    for( ivt=vtrks.begin(); ivt!=vtrks.end(); ivt++ ) {
      if( Covar_ >= 3 ) {
        HepMatrix &hm = (* new HepMatrix(3,3) );
        hm[0][0] = (*ivt).Enk(1,1); hm[0][1] = (*ivt).Enk(1,2); hm[0][2] = (*ivt).Enk(1,3);
        hm[1][0] = (*ivt).Enk(2,1); hm[1][1] = (*ivt).Enk(2,2); hm[1][2] = (*ivt).Enk(2,3);
        hm[2][0] = (*ivt).Enk(3,1); hm[2][1] = (*ivt).Enk(3,2); hm[2][2] = (*ivt).Enk(3,3);
        vrt->addCov( &hm );
      }
    }

    if( Covar_ >= 4 ) {
      // SET CROSS-MOMENTUM CORRELATIONS OF TRACKS IN VERTEX
      // Dummy values for the time being (03/02)
      for( ivt=vtrks.begin(); ivt!=vtrks.end(); ivt++ ) {
	list<CsVTrack>::iterator jvt = ivt;
	for (jvt++; jvt!=vtrks.end(); jvt++) {
	  HepMatrix &hm = (* new HepMatrix(3,3) );
	  vrt->addCov( &hm );
	}
      }
    }

    vrt->setTracks(trks); int nTrks = vrt->getNTracks();


    //     ********** OUTPUT STREAM and STATISTICS **********
    bool hasMup = false;    // "hasMup" (used infra) is set by the same token.
    bool hasBMS = false;    // Same thing for "hasBMS". 
    if (vrt->isPrimary()) {
      nPVrts++;
      int phiCandidate = nTrks==4 ? 0x1 : 0;  // Init phi search
      int DSCandidate  = nTrks>=5 ? 0x1 : 0;  // Init D* search
      int piNCandidate = 0;                   // Init 3/5pi search
      if      (nTrks==4) piNCandidate = 0x1;
      else if (nTrks==6) piNCandidate = 0x3;
      else if (nTrks==8) piNCandidate = 0x7;
      list<CsVTrack>::iterator imu, iKp, iKm;
      ivt = vtrks.begin(); 
      const CsTrack *beam = (*ivt).getAssociatedTrk();         // ***** HAS BMS?
      const list<CsZone*> &zones = beam->getZones();
      list<CsZone*>::const_iterator iz;
      for (iz = zones.begin(); iz!=zones.end(); iz++)
	if((*iz)->getName()=="BMSzone") { hasBMS = true; break; }
      TLorentzVector *missing = 0; static int piNQ; // 3/5pi missing E-P, charge
      for (ivt++; ivt!=vtrks.end(); ivt++) {

	//     ********** LOOP ON SPECTROMETER TRACKS **********

	if ((*ivt).getStatus()==5) {                           // ***** HAS mu'?
	  nPVrtsM++; hasMup = true; if (hasBMS) nPVrtsG++;
	  phiCandidate |= 0x2; DSCandidate |= 0x2; imu = ivt;
	  //#define doFitting_DIS_KINEMATICS 1
#ifdef doFitting_DIS_KINEMATICS
	  const double M_mu  = 0.1056583568, M2_mu = M_mu*M_mu;
	  CsVTrack &mui = vtrks.front(), &mus = *ivt;
	  TLorentzVector p(0,0,0,M_p);
	  // w = (1+tgx^2+tgy^2)^1/2; Px = P.tgx/w;  Py = P.tgy/ w; Pz = P/ w
	  // dPx/dtgx = Pz-Px.tgx/w^2; dPx/dtgy = -Py.tgx/w^2; dPx/d1/P = -Px.P
	  // dPz/dtgx =   -Px    /w^2; dPz/dtgy = -Py    /w^2; dPz/d1/P = -Pz.P
	  // E = (1/r^2+m^2) dE/dr = -1/E * P^3
	  // mu'
	  double tx = mus.Qnk(1), ty = mus.Qnk(2);
	  double w2 = 1+tx*tx+ty*ty, w = sqrt(w2);
	  double Ks = 1/fabs(mus.Qnk(3)), pz = Ks/w, px = tx*pz, py = ty*pz;
	  double Es = sqrt(Ks*Ks+M2_mu), dEsdr = -1/Es*Ks*Ks*Ks;
	  TLorentzVector ks(TVector3(px,py,pz),Es);
	  TLorentzVector dksdtx(TVector3(pz-px*tx/w2,  -px*ty/w2,-px/w2)  ,0);
	  TLorentzVector dksdty(TVector3(  -py*tx/w2,pz-py*ty/w2,-py/w2)  ,0);
	  TLorentzVector dksdr (TVector3(  -px*Ks ,    -py*Ks   ,-pz*Ks),dEsdr);
	  // mu
	  tx = mui.Qnk(1); ty = mui.Qnk(2); w2 = 1+tx*tx+ty*ty; w = sqrt(w2);
	  double Ki = 1/fabs(mui.Qnk(3)); pz = Ki/w; px = tx*pz; py = ty*pz;
	  double Ei = sqrt(Ki*Ki+M2_mu), dEidr = -1/Ei*Ki*Ki*Ki;
	  TLorentzVector ki(TVector3(px,py,pz),Ei);
	  TLorentzVector dkidtx(TVector3(pz-px*tx/w2,  -px*ty/w2,-px/w2) ,0);
	  TLorentzVector dkidty(TVector3(  -py*tx/w2,pz-py*ty/w2,-py/w2)  ,0);
	  TLorentzVector dkidr (TVector3(  -px*Ki ,    -py*Ki   ,-pz*Ki),dEidr);
	  TLorentzVector q = ki-ks; double Q2 = -q.M2();
	  double pq = q.Dot(p), pk = ki.Dot(p);
	  double xB = Q2/2/pq, y = pq/pk;
	  double dQ2ds[3] = { 2*dksdtx.Dot(q), 2*dksdty.Dot(q), 2*dksdr.Dot(q)};
	  double dQ2di[3] = {-2*dkidtx.Dot(q),-2*dkidty.Dot(q),-2*dkidr.Dot(q)};
	  double dpqds[3] = {-  dksdtx.Dot(p),-  dksdty.Dot(p),-  dksdr.Dot(p)};
	  double dpqdi[3] = {   dkidtx.Dot(p),   dkidty.Dot(p),   dkidr.Dot(p)};
	  double dxBds[3], dxBdi[3], dyds[3], dydi[3]; int i;
	  for (i = 0; i<3; i++) {
	    dxBds[i] = (dQ2ds[i]*pq-Q2*dpqds[i])/pq/pq/2;
	    dxBdi[i] = (dQ2di[i]*pq-Q2*dpqdi[i])/pq/pq/2;
	    dyds[i] = dpqds[i]/pk; dydi[i] = (dpqdi[i]*pk-pq*dpqdi[i])/pk/pk; 
	  }
	  double dQ2s, dQ2i, dxBs, dxBi, dys, dyi;
	  for (i = 1, dQ2s=dQ2i=dxBs=dxBi=dys=dyi = 0; i<=3; i++) {
	    double csii = mus.Dnk(i,i), ciii = mui.Dnk(i,i);
	    double dQ2si = dQ2ds[i-1], dQ2ii = dQ2di[i-1];
	    dQ2s += dQ2si*dQ2si*csii; dQ2i += dQ2ii*dQ2ii*ciii;
	    double dxBsi = dxBds[i-1], dxBii = dxBdi[i-1];
	    dxBs += dxBsi*dxBsi*csii; dxBi += dxBii*dxBii*ciii;
	    double dysi = dyds[i-1], dyii = dydi[i-1];
	    dys += dysi*dysi*csii; dyi += dyii*dyii*ciii;
	    for (int j = i+1; j<=3; j++) {
	      double csij = mus.Dnk(i,j), ciij = mui.Dnk(i,j);
	      double dQ2sj = dQ2ds[j-1], dQ2ij = dQ2di[j-1];
	      dQ2s += 2*dQ2si*dQ2sj*csij; dQ2i += 2*dQ2ii*dQ2ij*ciij; 
	      double dxBsj = dxBds[j-1], dxBij = dxBdi[j-1];
	      dxBs += 2*dxBsi*dxBsj*csij; dxBi += 2*dxBii*dxBij*ciij; 
	      double dysj = dxBds[j-1], dyij = dydi[j-1];
	      dys += 2*dysi*dysj*csij; dyi += 2*dyii*dyij*ciij; 
	    } 
	  }
	  double dQ2 = sqrt(dQ2s+dQ2i), dxB = sqrt(dxBs+dxBi), dy = sqrt(dys+dyi);
	  CsTrack *scmu = mus.getAssociatedTrk();
#  if doFitting_DIS_KINEMATICS > 1
	  double bChi2 = beam->getChi2()/(beam->getClusters().size()-4);
	  int sNHs = scmu->getClusters().size(), sNDWs = 0;
	  static int iwDW = -2; static unsigned int patDWs[2]; if (iwDW<-1) {
	    iwDW = -1;
	    list<CsDetector*> dets = CsGeom::Instance()->getDetectors();
	    list<CsDetector*>::iterator id; int i;
	    for (id = dets.begin(), i = 0, patDWs[0]=patDWs[1] = 0;
		 id!=dets.end(); id++, i++) {
	      if ((*id)->GetTBName().find("DW")==0) {
		int iw = i/32, ib = i%32; if (iwDW<0) iwDW = iw; iw -= iwDW;
		if (iw>1)
		  printf("doFitting_DIS_KINEMATICS:\a DWs don't fit w[%d:%d]\n",
			 iwDW,iwDW+1);
		else patDWs[iw] |= 1<<ib;
	      }
	    }
	  }
	  if (iwDW>=0) {
	    for (int iw = iwDW; iw<=(iwDW<(CSTRACK_MAPSIZE-1)?iwDW+1:(CSTRACK_MAPSIZE-1)); iw++) {
	      unsigned int fired = scmu->getFiredDetsBitmap()[iw];
	      unsigned int pat = patDWs[iw-iwDW]&fired, b;
	      for (b = 1; b!=0x80000000; b <<= 1) if (pat&b) sNDWs++;
	    }
	  }
	  double sChi2 = scmu->getChi2()/(sNHs-5);
#  endif
	  double bT = beam->hasMeanTime()?beam->getMeanTime():1000;
	  double sT = scmu->hasMeanTime()?scmu->getMeanTime():1000;
	  double bP = 1/beam->getHelices()[0].getCop();
	  double sP = 1/scmu->getHelices()[0].getCop();
	  const CsHelix &hL = scmu->getHelices().back();
	  double sL = const_cast<CsHelix&>(hL)(0)/1000;
	  printf("==================== DIS KINEMATICS ====================\n");
#  if doFitting_DIS_KINEMATICS > 1
	  printf("Evt %d#%9d #trks %d",
		 cse->getRunNumber(),cse->getEventNumberInRun(),vtrks.size());
	  printf("  BEAM %4.1fNDF %4.1fns %6.2f/%6.2f+/-%4.2fGeV",
		 bChi2,           bT,bP,Ki,sqrt(mui.Dnk(2,2))*Ki*Ki);
	  printf("  MU' %4.1fNDF #hits %d(%2d DWs) %4.1fm %4.1fns %6.2f/%6.2f+/-%4.2fGeV",
		 sChi2,sNHs,sNDWs,sL,sT,sP,Ks,sqrt(mus.Dnk(2,2))*Ks*Ks);
#  else
	  printf("Evt %d#%9d",cse->getRunNumber(),cse->getEventNumberInRun());
	  printf("  BEAM %4.1fns %6.2f/%6.2f+/-%4.2fGeV",
		 bT,bP,Ki,sqrt(mui.Dnk(2,2))*Ki*Ki);
	  printf("  MU' %4.1fns %6.2f/%6.2f+/-%4.2fGeV",
		 sT,sP,Ks,sqrt(mus.Dnk(2,2))*Ks*Ks);
#  endif
	  printf(" Z %5.1fcm Q2,x,y %.4g+/-%.3g %.3g+/-%.2g %.3f+/-%.3f\n",
		 Xn(1),Q2,dQ2,xB,dxB,y,dy);
	  printf("========================================================\n");
#  if doFitting_DIS_KINEMATICS > 1
	  static FILE *fp = 0; if (!fp) {
	    string fName("~/DIS_KINEMATICS.1.61.log"); CsOpt::expand(fName);
	    fp = fopen(fName.c_str(),"a");
	  }
	  fprintf(fp,"Evt %d#%9d #trks %d",
		 cse->getRunNumber(),cse->getEventNumberInRun(),vtrks.size());
	  fprintf(fp,"  BEAM %4.1fNDF %4.1fns %6.2f/%6.2f+/-%4.2fGeV",
		  bChi2,           bT,bP,Ki,sqrt(mui.Dnk(2,2))*Ki*Ki);
	  fprintf(fp,"  MU' %4.1fNDF #hits %d(%2d DWs) %4.1fm %4.1fns %6.2f/%6.2f+/-%4.2fGeV",
		 sChi2,sNHs,sNDWs,sL,sT,sP,Ks,sqrt(mus.Dnk(2,2))*Ks*Ks);
	  fprintf(fp," Z %5.1fcm Q2,x,y %.4g+/-%.3g %.3g+/-%.2g %.3f+/-%.3f\n",
		  Xn(1),Q2,dQ2,xB,dxB,y,dy);
#  endif
#  define doFitting_DIS_WRITEBACK
#  ifdef doFitting_DIS_WRITEBACK
	  if (Q2>.95) cse->outputRawEventToStream(iStream);
#  endif
#endif
	} // End mu'
	else {
	  const CsTrack *t = (*ivt).getAssociatedTrk();
	  if (t->getZones().size()==1) continue; // Exclude fringe-field tracks
	  double cop = (*ivt).Qnk(3);
	  if (phiCandidate&0x1) {
	    if      (!(phiCandidate&0x4) && cop>0) {
	      phiCandidate |= 0x4; iKp = ivt;
	    }
	    else if (!(phiCandidate&0x8) && cop<0) {
	      phiCandidate |= 0x8; iKm = ivt;
	    }
	  }
	  if (DSCandidate&0x1) {
	    if      (cop>0) DSCandidate |= 0x4;
	    else if (cop<0) DSCandidate |= 0x8;
	  }
	}
	if (piNCandidate) {                          // ***** 3/5pi MISSING MASS
	  // Here the 3/5 pions are collected. Incrementing the "piNCandidate"
	  // flag, the missing charge and the missing E-P. Note still that
	  // fringe-field track are excluded (cf. supra). Therefore the 3/5
	  // pions may not make it all up here.
	  double tx, ty, w2, w, K, px, py, pz, Epi;
	  if (piNCandidate&0x1) { // Init w/ beam E-P
	    CsVTrack &pii = vtrks.front();
	    tx = pii.Qnk(1), ty = pii.Qnk(2), w2 = 1+tx*tx+ty*ty, w = sqrt(w2);
	    K = 1/fabs(pii.Qnk(3)), pz = K/w, px = tx*pz, py = ty*pz;
	    Epi = sqrt(K*K+M2_pi);
	    missing = new TLorentzVector(0,0,0,M_p);
	    *missing += TLorentzVector(TVector3(px,py,pz),Epi);
	    piNQ = pii.Qnk(3)>0 ? -1 : 1;
	  }
	  CsVTrack &pis = *ivt; piNCandidate <<= 1;
	  piNQ +=  pis.Qnk(3)>0 ? 1 : -1;
	  tx = pis.Qnk(1), ty = pis.Qnk(2), w2 = 1+tx*tx+ty*ty, w = sqrt(w2);
	  K = 1/fabs(pis.Qnk(3)), pz = K/w, px = tx*pz, py = ty*pz;
	  Epi = sqrt(K*K+M2_pi);
	  *missing += TLorentzVector(TVector3(-px,-py,-pz),-Epi);
	  if      (piNCandidate==0x8 && piNQ==0) {
	    hm_pi3m->Fill(missing->M2()); hm_pi3e->Fill(missing->E());
	  }
	  else if (piNCandidate==0x60 && piNQ==0) {
	    hm_pi5m->Fill(missing->M2()); hm_pi5e->Fill(missing->E());
	    if (fabs(missing->E()-M_p)<5 && (doStream&phStream))
	      // ***** WRITE-BACK 5pi CANDIDATES w/ |missing E| < 5 GeV
	      cse->outputRawEventToStream(iStream);
	  }
	  else if (piNCandidate==0x380 && piNQ==0) {
	    if (fabs(missing->E()-M_p)<5 && (doStream&phStream))
	      // ***** WRITE-BACK 7pi CANDIDATES w/ |missing E| < 5 GeV
	      cse->outputRawEventToStream(iStream);
	  }
	  //#define doFitting_DEBUG_piN
#ifdef doFitting_DEBUG_piN
	  printf("3/5pi missing mass: 0x%04x %5.1f GeV  => %+d (%5.1f %4.1f %4.1f %5.1f GeV) %.2f GeV^2\n",
		 piNCandidate,Epi,
		 piNQ,missing->Px(),missing->Py(),missing->Pz(),missing->E(),missing->M2());
#endif
	}
      } // End loop on tracks in pVertex
      if (missing) delete missing; // Cleaning 3-pi candidate attempt
      if (doStream&pStream)              // ***** WRITE-BACK ALL PRIMARIES *****
	cse->outputRawEventToStream(iStream);
      //                                    ***** WRITE-BACK phi'S... *****
      if (phiCandidate==0xf) {           // mu, K+, K- found
	static HepLorentzVector lvp(0,0,0,M_p);
	const double M_mu  = 0.1056583568;
	const HepLorentzVector &lvbeam = (*vtrks.begin()).LzVec(M_mu);
	const HepLorentzVector &lvmu = (*imu).LzVec(M_mu);
	const HepLorentzVector &lvKp = (*iKp).LzVec(M_K);
	const HepLorentzVector &lvKm = (*iKm).LzVec(M_K);
	HepLorentzVector lvV0 = lvKp+lvKm;                  // V0
	HepLorentzVector lvRp = lvp+lvbeam-lvmu-lvV0;       // Recoil proton
	double dEK = (lvRp.m2()-M_p*M_p)/2/M_p;             // Inelasticity
	if (hist_) // Possible double counting! Histo filled for all pVertices
	  hdE->Fill(dEK);
	if (dEK<2.5) {                                    // Cut on inelasticity
	  double m_KK = lvV0.m();
	  if (m_KK<1.08) {
	    if (doStream&phStream) cse->outputRawEventToStream(iStream);
	    //#define doFitting_DEBUG_WRITEBACK
#ifdef doFitting_DEBUG_WRITEBACK
	    printf("WriteBack: phi KK %f\n",m_KK);
#endif
	    if (hist_) // Possible double counting! Histo filled for all pVert's
	      hm_KK->Fill(m_KK);
	  }
	}
      }
      //                                    ***** WRITE-BACK D*'s *****
      if (DSCandidate==0xf) {           // mu, +, - found
	list<CsVTrack>::iterator ivt1, ivt2, ivt3, ivt0;
	ivt1 = vtrks.begin(); ivt1++; ivt1++ /* skip beam, mu' */; ivt0 = ivt1;
	for (; ivt1!=vtrks.end(); ivt1++) { // Loop on spectro tracks for K
	  CsVTrack &K = *ivt1;
	  const HepLorentzVector &lvK = K.LzVec(M_K);
	  for (ivt2 = ivt0; ivt2!=vtrks.end(); ivt2++) { // Loop on for pi
	    CsVTrack &pi = *ivt2;
	    if (K.Qnk(3)*pi.Qnk(3)>0) continue; // pi's charge opposite K's
	    const HepLorentzVector &lvpi = pi.LzVec(M_pi);
	    HepLorentzVector lvKpi = lvK+lvpi;
	    const double M_D02 = M_D0*M_D0; double m_Kpi = lvKpi.m();
	    if (m_Kpi-M_D0<D0LowCut || D0UpCut<m_Kpi-M_D0) continue;
	    Hep3Vector pKpi = lvK.v(); pKpi += lvpi.v(); 
	    HepLorentzVector lvD0(pKpi,sqrt(pKpi*pKpi+M_D02));
	    for (ivt3 = ivt0; ivt3!=vtrks.end(); ivt3++) { // Loop on for piS
	      if (ivt3==ivt1 || ivt3==ivt2) continue;
	      if (K.Qnk(3)*(*ivt3).Qnk(3)>0) continue;
	      HepLorentzVector lvD0pi = (*ivt3).LzVec(M_pi); lvD0pi += lvD0;
	      double m_D0pi = lvD0pi.m();
	      const double M_DS  = 2.0100;
	      if (DSLowCut<m_D0pi-M_DS && m_D0pi-M_DS<DSUpCut) {
		if (doStream&phStream) cse->outputRawEventToStream(iStream);
		if (hist_)
		  // Possible double counting! Histo filled for all pVert's
		  hm_Kpi->Fill(m_Kpi);
#ifdef doFitting_DEBUG_WRITEBACK
		printf("WriteBack: D* Kpi %f D0pi %f\n",m_Kpi,m_D0pi);
#endif
	      }
	    }
	  }
	}
      }

      //         ********** BEST PRIMARY VERTEX **********
      // (Cf. definition in this file's header.)
      if (BpVDY_)
	// (Note the in the DY case, "hasMup|BMS" and "piNCandidate" are
	// meaningless. Yet, we update them, to be on the safe side...)
	getBestVertexDY(vrt,bestV,nGoodTrksBV,nMuswTBV,tChi2BV,
			hasMup,      hasMupBV,
			hasBMS,      hasBMSBV,
			piNCandidate,piNBV);
      else
	getBestVertex  (vrt,bestV,nonInteractingBV,
			hasMup,      hasMupBV,
			hasBMS,      hasBMSBV,
			piNCandidate,piNBV);

      //#define doFitting_PRINT 2
#ifdef doFitting_PRINT      //             ********** PRINT-OUT **********
      if (evTrig || cse->isAMonteCarloEvent()) {
	// Printed are statistics for current event and cumulated statistics.
	// For the latter the aim is to consider only the best vertex, which is
	// not yet known => Show the cumulated statistics obtained by adding
	// the current vertex to the statistics of previous event. (And
	// statistics themselves are updated later on.)
	int npvrts = statistics_[10]+1, ntrks = statistics_[6]+nTrks;
	int mup = hasMup?1:0, mupBMS = (hasBMS&hasMup)?1:0, ntmup = hasMup?nTrks:0;
# if doFitting_PRINT > 2
	const list<CsTrack*> &trks = vrt->getTracks();
	Cs3Vector parb, part; bool vtxOK = true; double theta;
	vtxOK &= vrt->getPar(trks.front(),parb);
	vtxOK &= vrt->getPar(trks.back(), part);
	if (vtxOK) {
	  double tx, ty, a, tDCZ, tDCX, tDCY, bDCZ, bDCX, bDCY;
	  tx = parb  (0); ty = parb  (1); a = sqrt(1+tx*tx+ty*ty);
	  bDCZ = 1/a, bDCX = tx/a, bDCY = ty/a;
	  tx = part  (0); ty = part  (1); a = sqrt(1+tx*tx+ty*ty);
	  tDCZ = 1/a; tDCX = tx/a; tDCY = ty/a;
	  theta   = acos(tDCX*bDCX+tDCY*bDCY+tDCZ*bDCZ)*1000;
	}
	else theta = 0;
	printf("pVtx Z,th,chi2,#T  #V,#T,T/V,mu'..: 0x%04x %5.0f %6.2f %.1f %2d   %4d %5d %.3f %4d %4d %5d\n",
	       evTrig,vrt->getZ()/10,theta,vrt->getChi2()/(2*nTrks-3),nTrks,npvrts,ntrks,
	       (float)ntrks/npvrts,statistics_[12]+mup,statistics_[2]+mupBMS,
	       statistics_[7]+ntmup);
#else
	printf("pVtx chi2,#T  #V,#T,T/V,mu'...: 0x%04x %.1f %2d   %4d %5d %.3f %4d %4d %5d\n",
	       evTrig,vrt->getChi2()/(2*nTrks-3),nTrks,npvrts,ntrks,
	       (float)ntrks/npvrts,statistics_[12]+mup,statistics_[2]+mupBMS,
	       statistics_[7]+ntmup);
#  endif
#  if doFitting_PRINT > 1
	ivt = vtrks.begin();
	CsTrack *t = (*ivt).getAssociatedTrk();
	printf("tracks(t,p): %.1f %.3f  ",t->getMeanTime(),
	       1/t->getHelices()[0].getCop());
	for (ivt++;  ivt!=vtrks.end(); ivt++) {
	  CsVTrack &vt = *ivt; CsTrack *t = vt.getAssociatedTrk();
	  printf("  %.1f %.3f",t->hasMeanTime()?t->getMeanTime():1000,
		 1/t->getHelices()[0].getCop());
	}
	printf("\n");
#  endif
      }
#endif
      if (hist_) {
	//        *************** HISTOGRAMS ***************
	// Note: these histograms are filled for all vertices. Would be better
	// to fill them w/ best vertex only, given that extra vertices usually
	// come, at least in the muon case, from non interacting incident
	// particles continued into the spectrometers, which constitutes a very
	// specific, and uninteresting case. BTW, some histos are already filled
	// w/ best vertex only, cf. infra.
	Pulls( vtrks, Cn, Xn );

	//                                          ********** TRACKS **********
	double bDCX(0), bDCY(0), bDCZ(0);
	float vChi2 = vrt->getChi2(); int vNDF = 2*nTrks-3;
	float vProb = prob_(&vChi2,&vNDF);
        for (ivt = vtrks.begin(); ivt!=vtrks.end(); ivt++) {
          CsVTrack &vt = *ivt; CsTrack *t = vt.getAssociatedTrk();
          int tNDF = 2; float tChi2 = vt.getChi2(); hTrkVChi2->Fill(tChi2/tNDF);
	  float tProb = prob_(&tChi2,&tNDF );       hTrkVProb->Fill(tProb);
	  hTvsVProb->Fill(vProb,tProb);
	  if (ivt==vtrks.begin()) {                                // ***** BEAM
            double a = sqrt(1+vt.getDXDZ()*vt.getDXDZ()+vt.getDYDZ()*vt.getDYDZ());
            bDCZ = 1/a; bDCX = vt.getDXDZ()/a; bDCY = vt.getDYDZ()/a;
          }
	  else {                                         // ***** SPECTRO TRACKS
            double a = sqrt(1+vt.getDXDZ()*vt.getDXDZ()+vt.getDYDZ()*vt.getDYDZ());
            double DCZ = 1/a, DCX = vt.getDXDZ()/a, DCY = vt.getDYDZ()/a;
	    double theta = acos(DCX*bDCX+DCY*bDCY+DCZ*bDCZ);
	    hProb_T->Fill(tProb,theta); hProb_Mom->Fill(tProb,vt.getMom());
          }
        } // End loop on CsVTrack's in pVertex


	hSigmaX->Fill( 10 * sqrt( Cn(2,2) ) );
        hSigmaY->Fill( 10 * sqrt( Cn(3,3) ) );
        hSigmaZ->Fill( 10 * sqrt( Cn(1,1) ) );
	
	hSigmaX_N->Fill( nTrks, 10 * sqrt( Cn(2,2) ) );
        hSigmaY_N->Fill( nTrks, 10 * sqrt( Cn(3,3) ) );
        hSigmaZ_N->Fill( nTrks, 10 * sqrt( Cn(1,1) ) );
	
        hSigmaZ_X->Fill( vrt->getZ(), 10 * sqrt( Cn(2,2) ) );
        hSigmaZ_Y->Fill( vrt->getZ(), 10 * sqrt( Cn(3,3) ) );
        hSigmaZ_Z->Fill( vrt->getZ(), 10 * sqrt( Cn(1,1) ) );
	
        hPullX->Fill( ( X0(2) - vrt->getX() ) / sqrt( fabs( Cn(2,2) - C0[0][0] ) ) );
        hPullY->Fill( ( X0(3) - vrt->getY() ) / sqrt( fabs( Cn(3,3) - C0[1][1] ) ) );
        hPullZ->Fill( ( X0(1) - vrt->getZ() ) / sqrt( fabs( Cn(1,1) - C0[2][2] ) ) );
	
      }
    }
    else { // ************************* SECONDARIES *************************
      if (trks.size()!=2)
	CsErrLog::mes(elFatal,"Number of tracks in secondary < 2");
      if (!reTrackingON)  // If not re-vertexing after retracking...
	statistics_[5]++;
      // ...else the statistics have already been updated (maybe not accurately
      // though: re-tracking may yield fewer/more vertices than first tracking.)
      if (doStream&sStream)  // Output all secondaries
	cse->outputRawEventToStream(iStream);

      //     *************** SECONDARIES' STREAMS ***************
      double Zs = vrt->getZ();                          // ***** V0 SEARCH *****
      CsTrack* trk1 = trks.front(), *trk2 = trks.back();
      if (trk1->getZones().size()==1 || trk2->getZones().size()==1)
	continue;  // Exclude fring-field tracks
      const Hep3Vector &v1 = vtrks.front().Vec(), &v2 = vtrks.back().Vec();
      Hep3Vector v0(v1); v0 += v2;
      double pT = v1.perp(v0); if (pT<V0pTCut) continue;
      double p02 = v0*v0;
      double p12 = v1*v1, e1p = sqrt(p12+M2_p), e1pi = sqrt(p12+M2_pi);
      double p22 = v2*v2, e2p = sqrt(p22+M2_p), e2pi = sqrt(p22+M2_pi);
      double m_ppi  = sqrt((e1p+e2pi)*(e1p+e2pi)-p02);
      double m_pip  = sqrt((e1pi+e2p)*(e1pi+e2p)-p02);
      double m_pipi = sqrt((e1pi+e2pi)*(e1pi+e2pi)-p02);
      bool isLambda  = fabs(m_ppi-M_Lam)<dm_Lam;
      bool isALambda = fabs(m_pip-M_Lam)<dm_Lam;
      bool isK0 = fabs(m_pipi-M_K0)<dm_K0;
      double e1mu = sqrt(p12+M2_mu) , e2mu = sqrt(p22+M2_mu);
      double m_mumu = sqrt((e1mu+e2mu)*(e1mu+e2mu)-p02);
      if (m_mumu>1.5) {
	// I.e.: typically the lower limit of the fitting range of the J/psi
	// peak in DY w/ hadron absorber
	unsigned int isJpsi = 0; if (trk1->getXX0()>15) isJpsi |= 0x1;
	else if                     (trk2->getXX0()>15) isJpsi |= 0x2;
	if (isJpsi) {
	  if (doStream&phStream) cse->outputRawEventToStream(iStream);
	  if (isJpsi==0x3) hm_mumu->Fill(m_mumu);
	}
      }
      if (!(isLambda || isALambda || isK0)) continue;
      double Xs = vrt->getX(), Ys = vrt->getY();
      list<CsVertex*>::iterator ivp;
      for (ivp = vrts.begin(); ivp!=vrts.end(); ivp++) {
	CsVertex *vrtp = *ivp; if (vrtp->isPrimary()) { // Loop on primaries
	  double Zp = vrtp->getZ(); if (Zp>Zs) continue;
	  double Xp = vrtp->getX(), Yp = vrtp->getY();
	  HepVector vv(3); vv[0] = Xs-Xp; vv[1] = Ys-Yp; vv[2] = Zs-Zp;
	  HepVector tmps(vv), tmpp(vv);
	  tmps = *(vrt->getCov(0))*tmps; tmpp = *(vrtp->getCov(0))*tmpp;
	  double dist = vv.norm(), ddist = dot(tmps,vv)+dot(tmpp,vv);
	  ddist = sqrt(ddist)/fabs(dist); if (dist<V0vvCut*ddist) continue;
	  Hep3Vector vv3(Xs-Xp,Ys-Yp,Zs-Zp);
	  double ctheta = vv3*v0/sqrt(p02)/dist;
	  if (ctheta>V0cthCut) {
	    if (doStream&phStream &&
		(isK0 ||
		 hadronRun && // Hadron run ...
		 // ..., at least the early low intensity part of 2008, is clean
		 // enough that the (A)Lambda mass range can be written back.
		 (isLambda || isALambda)))
	      cse->outputRawEventToStream(iStream);
	    if (hist_ && vrtp==bestV) {
	      if (isLambda)  hm_ppi->Fill(m_ppi);
	      if (isALambda) hm_pip->Fill(m_pip);
	      if (isK0)      hm_pipi->Fill(m_pipi);
	    }
#ifdef doFitting_DEBUG_WRITEBACK
	    if (hadronRun) {
	      if (isLambda)  printf("WriteBack: Lambda ppi %f\n",m_ppi);
	      if (isALambda) printf("WriteBack: ALambda pip %f\n",m_pip);
	    }
	    if (isK0)      printf("WriteBack: K0 pipi %f\n",m_pipi);
#endif
	  }
	}
      }
    }  // End primary else secondary
  } // End loop on candidate vertices

  //                            ***** ERASE VERTICES WHICH DID NOT PASS CRITERIA
  for( iv=vrt_to_erase.begin(); iv!=vrt_to_erase.end(); iv++ ) {
    list<CsVertex*>::iterator iv_erase = find( vrts.begin(), vrts.end(), *iv );
    delete *iv;
    vrts.erase( iv_erase );
  }
  vrt_to_erase.clear();

  if (bestV) {
    // ***** CHECK BEST pVERTEX is COMPATIBLE W/ T0 USED in TRACKING
    const CsTrack *bestVBeam = bestV->getTracks().front();
    if (T0 && bestVBeam->hasMeanTime()) {
      double bestVTime = bestVBeam->getMeanTime();
      if (fabs(bestVTime-*T0)>rebuildTracksCut_) {
	if (!nonInteractingBV) {
	  if (reTrackingON)
	    CsErrLog::msg(elError, __FILE__, __LINE__,
 "Wrong (not best pVertex') T0 used in track fitting: %.2f instead of %.2f ns!",
			  *T0,bestVTime);
	  else {
	    CsErrLog::msg(elError, __FILE__, __LINE__,
 "Wrong (not best pVertex') T0 used in track fitting: %.2f instead of %.2f ns"
			" => Re-tracking requested",
			  *T0,bestVTime);
	    *T0 = bestVTime; return false;
	  }
	}
	else if (Print_[0] || Print_[3]) {
	  printf("Best pVertex: T0(=%.2f ns) != track fitting's(=%.2f ns). "
		 "It's nonInteracting, though.\n",bestVTime,*T0);
	  if (!reTrackingON) printf("   => Re-tracking NOT requested\n");
	}
      }
    }
    bestV->setBestVertex(true);
  }
  else { // No best primary vertex, meaning vertexing failed:...
    // ...update statistics about the reason for failure.
    if      (rejection&0x2) statistics_[4]++;
    else if (rejection&0x1) statistics_[3]++;
  }

  BestPVertex_ = bestV; // Set CsKalmanFitting's pointer to Best Primary Vertex

  if (bestV && !hasMupBV &&       // ***** BEST pV W/O MU': LET'S CHECK AGAIN...
      specials.size()>1) { // ...in case there are several mu' candidates
    // (Indeed, several mu' may lead to the situation where the pV that wins the
    // best pV contest turns out to have been built on a some choice for the mu'
    // (so-called "special", cf. "CsAverPattern::FindPrimary") that was later on
    // rejected, while retaining another candidate mu' (not "special"): it would
    // only have had to get a slightly better chi2 than the pV w/ exactly the
    // same set of 2ndary tracks but built on the correct choice of mu'. Some 
    // action has already been taken to rescue (meaning get it w/ a "hasMup"
    // flag) such a pV, cf. "(vChi2<chi2BV+X && same && hasMup && !hasMupBV"
    // selection criterion supra, but it may not suffice if "X" is kept low so
    // as to avoid selecting a best pV w/ too bad a chi2 (so that in turn, e.g.,
    // it does not get overridden but yet another pV). Here, we try another 
    // work-around, by re-evaluating the "hasMup" flag.)
    const list<CsTrack*> trks = bestV->getTracks();
    list<CsTrack*>::const_iterator it = trks.begin();
    for (it++; it!=trks.end(); it++) { // Loop on 2ndary tracks
      if (specials.find(*it)!=specials.end()) { hasMupBV = true; break; }
    }
  }

  if (PrimReduce_==1 && bestV) {

    //   *************** (UPON OPTION) RETAIN ONLY BEST VERTEX ***************

    iv = vrts.begin(); while (iv!=vrts.end()) {
      CsVertex *vrt = *iv;
      if (vrt->isPrimary() && vrt!=bestV) {
	// Erase all other vertices. Update list of vertices in tracks
	list<CsTrack*> trks = vrt->getTracks();
	for (it = trks.begin(); it!=trks.end(); it++) (*it)->rmVertex(vrt);
	delete vrt;
	vrts.erase(iv++);
      }
      else iv++;
    }
    nPVrts = 1; nPVrtsM = hasMupBV?1:0; nPVrtsG = hasBMSBV?1:0;
  }

  //              *************** HISTOGRAMMING ***************
  //  Histograms and statistics of #tracks and other vertex characteristics are
  // made to correspond to the best vertex. It's the only meaningful choice in
  // the muon case, where accidentally coincident interactions are unlikely and
  // all but best vertex are expected to correspond to non interacting incident
  // muons, continued into the spectrometer, that can at most highjack tracks
  // from the main interaction.
  bool hasBokBV = // Best pVertex w/ BMS fully OK, including back-propagation
    false;
  if (hist_) {

    //                                                           ***** #VERTICES
    hNVtx->Fill(nPVrts);      hNVtxM->Fill(nPVrtsM);    hNVtxG->Fill(nPVrtsG);
    if (evTrig&stdTrigs) {                // Single out ``standard'' triggers...
      hNVtx_SI->Fill(nPVrts); hNVtxM_SI->Fill(nPVrtsM); hNVtxG_SI->Fill(nPVrtsG);
    }
    if ((evTrig&pureTrigs)==evTrig) {     // ...and CT (+? highQ2) or beam
      hNVtx_C->Fill(nPVrts);  hNVtxM_C->Fill(nPVrtsM);  hNVtxG_C->Fill(nPVrtsG);
    }
    if (evTrig&specialTrigs) {            // ...and OT or NIT or DVCS
      hNVtx_O->Fill(nPVrts);  hNVtxM_O->Fill(nPVrtsM);  hNVtxG_O->Fill(nPVrtsG);
    }

    if (bestV) {                 // ***** HISTOGRAM BEST VERTEX STATISTICS *****
      // Note: Histograms from phi candidate et al., are filled for all
      // vertices, cf. supra. So are pseuo-pulls.
      int nTrksBV = bestV->getNTracks();
      hNTrk->Fill(nTrksBV); if (hasMupBV) {                     // ***** #TRACKS
	hNTrkM->Fill(nTrksBV);      if (hasBMSBV) hNTrkG->Fill(nTrksBV);
      }
      if (evTrig&stdTrigs) {              // Single out ``standard'' triggers...
	hNTrk_SI->Fill(nTrksBV); if (hasMupBV) {
	  hNTrkM_SI->Fill(nTrksBV); if (hasBMSBV) hNTrkG_SI->Fill(nTrksBV);
	}
      }
      if ((evTrig&pureTrigs)==evTrig) {   // ...CT (+? highQ2) or beam
	hNTrk_C->Fill(nTrksBV);  if (hasMupBV) {
	  hNTrkM_C->Fill(nTrksBV);  if (hasBMSBV) hNTrkG_C->Fill(nTrksBV);
	}
      }
      if (evTrig&specialTrigs) {          // ...OT or NIT or DVCS
	hNTrk_O->Fill(nTrksBV);  if (hasMupBV) {
	  hNTrkM_O->Fill(nTrksBV);  if (hasBMSBV) hNTrkG_O->Fill(nTrksBV);
	}
      }

      int vNDF = 2*nTrksBV-3; float vChi2 = bestV->getChi2();      // ***** CHI2
      hVrtProb->Fill(prob_(&vChi2,&vNDF ));
      double vChi2BV = vChi2/vNDF; hVrtChi2->Fill(vChi2BV);

      // ********** HISTOGRAM (some of) BEST VERTEX CHARACTERISTICS **********
      // Note: these histograms are filled w/ best vertex only, given that extra
      // vertices usually come, at least in the muon case, from non interacting
      // incident particles continued into the spectrometers, which constitutes
      // a very specific, and uninteresting case. Some other histos are still
      // filled w/ all vertices only, cf. supra.

      //                                                    ***** VERTEX COORD'S
      double xV = bestV->getX(), yV =  bestV->getY(), zV = bestV->getZ();
      hVrtX->Fill(xV); hVrtY->Fill(yV); hVrtZ->Fill(zV);

      if (cse->isAMonteCarloEvent()) {                        // ***** MC: PULLS
	double XXMC,XYMC,XZMC;
	CsGeant3* geant3 = CsGeant3::Instance();
	const CsMCVertex *vrtMC = CsGeant3::Instance()->getMCVertices().front();
	double xVMC = vrtMC->getX(), yVMC = vrtMC->getY(), zVMC = vrtMC->getZ();
	const HepMatrix &cov = *(bestV->getCov(0));
	hDeltaKalman[0]->Fill(xV-xVMC);
	hDeltaKalman[1]->Fill(yV-yVMC);
	hDeltaKalman[2]->Fill(zV-zVMC);
	hPullsMCvrt[0]-> Fill((xV-xVMC)/sqrt(cov(2,2)));
	hPullsMCvrt[1]-> Fill((yV-yVMC)/sqrt(cov(3,3)));
	hPullsMCvrt[2]-> Fill((zV-zVMC)/sqrt(cov(1,1)));
      }

      //                                                       ********** TRACKS
      list<CsTrack*> trks = bestV->getTracks();
      list<CsTrack*>::const_iterator it = trks.begin();
      CsTrack *beam = *it;                                         // ***** BEAM
      CsBeam *csbeam = dynamic_cast<CsBeam*>(beam);
      if (csbeam) hasBokBV = !csbeam->getChi2CutFlag();
      Cs3Vector par;
      if (!bestV->getPar(beam,par)) CsErrLog::msg(elFatal,__FILE__,__LINE__,
 "Evt #%d: Vertex has track w/o vertex parameters",cse->getEventNumberInRun());
      double tx = par(0), ty = par(1), a = sqrt(1+tx*tx+ty*ty);
      double bDCz = 1/a, bDCx = tx*bDCz, bDCy = ty*bDCz;
      double bP = 1/fabs(par(2)); hBeamP->Fill(bP);
      double bT, bdT; if (beam->hasMeanTime()) {
	bT =  beam->getMeanTime(); hBeamT->Fill(bT);
	bdT = beam->getMeanTimeError();
	if (Print_[0])
	  printf("\nCsKF::doFitting==> Best Vertex %.2f,%.2f,%.2f cm,%.2f ns\n",
	       xV,yV,zV,bT);
      }
      else {
	bT = 0; bdT = 0;
	CsErrLog::msg(elError,__FILE__,__LINE__,
 "Best Vertex %.2f,%.2f,%.2f cm, has no Timing\n",xV,yV,zV);
      }
      const CsTrack *mus; double musP, Q2 = 0, xB = 0, y = 0, thetaMx;
      for (it++, musP=thetaMx = 0, mus = 0; it!=trks.end(); it++) {
	CsTrack *t = *it;                                // ***** SPECTRO TRACKS
	if (!bestV->getPar(t,par)) CsErrLog::msg(elFatal,__FILE__,__LINE__,
 "Evt #%d: Vertex has track w/o vertex parameters",cse->getEventNumberInRun());
	tx = par(0); ty = par(1); a = sqrt(1+tx*tx+ty*ty);
	// a = (1+tgx^2+tgy^2)^1/2; DCz = 1/a; DCx = tgx/a;  DCy = tgy/a
	double tDCz = 1/a, tDCx = tx*bDCz, tDCy = ty*bDCz, tP = 1/fabs(par(2));
	double theta = acos(bDCx*tDCx+bDCy*tDCy+bDCz*tDCz);
	if (theta>thetaMx) thetaMx = theta;
	if (t->hasMeanTime()) {
	  double T = t->getMeanTime(); hTrkT  ->Fill(T-bT);
	  double tdT = t->getMeanTimeError(), dT = sqrt(tdT*tdT+bdT*bdT);
	  hTrkTdT->Fill((T-bT)/dT);
	}
	if (piNBV) hpiNChi2->Fill(t->getChi2()/(t->getNDFs()-5));
	if (specials.find(*it)!=specials.end() &&          // ***** SCATTERED MU
	    (!mus || tP>musP)) { // If several mu's, retain highest momentum.
	  mus = t; musP = tP;
	  // Pz = P.DCz = P/a; Px = tgx.P/a; Py = tgy.P/a; E = (P^2+m^2)
	  // mu'
	  double pz = tDCz*tP, px = tx*pz, py = ty*pz, Es = sqrt(tP*tP+M2_mu);
	  TLorentzVector ks(TVector3(px,py,pz),Es);
	  // mu
	  pz = bDCz*bP; px = bDCx*bP; py = bDCy*bP; double Ei = sqrt(bP*bP+M2_mu);
	  TLorentzVector ki(TVector3(px,py,pz),Ei);
	  TLorentzVector p(0,0,0,M_p), q = ki-ks; Q2 = -q.M2();
	  double pq = q.Dot(p), pk = ki.Dot(p); xB = Q2/2/pq; y = pq/pk;
	}
      } // End loop on CsVTrack's in pVertex
      if (mus) {
#ifdef doFitting_PRINT   //  ***** MONITORING of PERFS: PRINT-OUT and WRITE-BACK
	if (Q2>.1) {
	  printf("CsKF: Evt %d Q2,xB,y %6.2f,%6.4f,%5.2f #tracks %d\n",
		 eNum,Q2,xB,y,(int)trks.size());
	}
#endif	  
	hKZ_Q2->Fill(zV,Q2); if (trks.size()>2) hKZ_Q21h->Fill(zV,Q2);
	//              ***** HISTOGRAM Q2 W/ INCREASINGLY DEMANDING CRITERIA...
	int binMx = 0; if (.05<y && y<.95) {  // ...Physically reachable y range
	  binMx = 1;
	  if (zTMn<zV && zV<zTMx && xV*xV+yV*yV<zTR2) {        // ...W/in target
	    binMx = 2; if (hasBMSBV) binMx = hasBokBV ? 4 : 3; // ...W/ BMS[&OK]
	  }
	}
	for (int bin = 0; bin<=binMx; bin++) {
	  hBVQ2->Fill(Q2,bin);                           // All triggers
	  if (evTrig&stdTrigs) hBVQ2_SI->Fill(Q2,bin);   // Standard triggers
	  if ((evTrig&pureTrigs)==evTrig)
	    hBVQ2_C->Fill(Q2,bin);                       // CT[+highQ2T] or beam
	  if (evTrig&specialTrigs) hBVQ2_O->Fill(Q2,bin);// OT or NIT or DVCS
	  if (nTrksBV<3) continue;          // Requiring now one hadron in FS...
	  hBVQ2h->Fill(Q2,bin);
	  if (evTrig&stdTrigs) hBVQ2h_SI->Fill(Q2,bin);
	  if ((evTrig&pureTrigs)==evTrig)
	    hBVQ2h_C->Fill(Q2,bin);
	  if (evTrig&specialTrigs) hBVQ2h_O->Fill(Q2,bin);
	}
      }
      hKZ_T->Fill(zV,thetaMx);
    }
  }


  //              *************** HISTOGRAMMING ONLINE MONITORING ***************
  if (hist_) {

    if (bestV) {                 // ***** HISTOGRAM BEST VERTEX STATISTICS *****
      int nTrksBV = bestV->getNTracks();
      double xV = bestV->getX(), yV =  bestV->getY(), zV = bestV->getZ();
      //                                                       ********** TRACKS
      list<CsTrack*> trks = bestV->getTracks();
      list<CsTrack*>::const_iterator it = trks.begin();
      CsTrack *beam = *it;                                         // ***** BEAM
      Cs3Vector par;

      if (!bestV->getPar(beam,par))
	CsErrLog::msg(elFatal,__FILE__,__LINE__,
		      "Evt #%d: BEAM1 Vertex has track w/o vertex parameters",
		      cse->getEventNumberInRun());
      double tx = par(0), ty = par(1), a = sqrt(1+tx*tx+ty*ty);
      double bDCz = 1/a, bDCx = tx*bDCz, bDCy = ty*bDCz;
      double bP = 1/fabs(par(2)); 
      double bT, bdT; 
      
      const CsTrack *mus; 
      double musP, Q2 = 0, xB = 0, y = 0, nu=0;
      TLorentzVector q(0,0,0,0);
      list<CsTrack*>::const_iterator qqivt  = trks.begin();
      //      beam = (*ivt).getAssociatedTrk();
      if (!bestV->getPar(beam,par)) CsErrLog::msg(elFatal,__FILE__,__LINE__,
 "Evt #%d: BEAM2 Vertex has track w/o vertex parameters",cse->getEventNumberInRun());
      double tbx = par(0), tby = par(1), ba = sqrt(1+tbx*tbx+tby*tby);
      double bbDCz = 1/ba, bbDCx = tbx*bbDCz, bbDCy = tby*bbDCz;
      
      double bbP = 1/fabs(par(2))/ba;
      double bbPz= 1/fabs(par(2))   ;
      double bbPx= bbPz * tbx;
      double bbPy= bbPz * tby;
      double bbTx=        tbx;
      double bbTy=        tby;

      
      bool mumuprim=false;

      for (it++, musP = 0, mus = 0; it!=trks.end(); it++) {
	CsTrack *t = *it;                                // ***** SPECTRO TRACKS

      if (!bestV->getPar(t,par))
	CsErrLog::msg(elFatal,__FILE__,__LINE__,
		      "Evt #%d: HADRONS Vertex has track w/o vertex parameters",
		      cse->getEventNumberInRun());
	tx = par(0); ty = par(1); a = sqrt(1+tx*tx+ty*ty);
	double tDCz = 1/a, tDCx = tx*bDCz, tDCy = ty*bDCz, tP = 1/fabs(par(2));
	double theta = acos(bDCx*tDCx+bDCy*tDCy+bDCz*tDCz);
	if (specials.find(*it)!=specials.end() &&          // ***** SCATTERED MU
	    (!mus || tP>musP)) { // If several mu's, retain highest momentum.
	  mus = t; musP = tP;
	  // Pz = P.DCz = P/a; Px = tgx.P/a; Py = tgy.P/a; E = (P^2+m^2)
	  // mu'
	  double pz = tDCz*tP, px = tx*pz, py = ty*pz, Es = sqrt(tP*tP+M2_mu);
	  TLorentzVector ks(TVector3(px,py,pz),Es);
	  // mu
	  pz = bDCz*bP; px = bDCx*bP; py = bDCy*bP; double Ei = sqrt(bP*bP+M2_mu);
	  TLorentzVector ki(TVector3(px,py,pz),Ei);
	  TLorentzVector p(0,0,0,M_p); q = ki-ks; Q2 = -q.M2(); nu = Ei-Es;
	  double pq = q.Dot(p), pk = ki.Dot(p); xB = Q2/2/pq; y = pq/pk;
	  
	  if(Q2>.7 && xB>0 && xB<1 && y>0 && y<1){

	    mumuprim=true;

	    OM_Q2 ->Fill(float(Q2));
	    OM_xbj->Fill(float(xB));
	    OM_ybj->Fill(float(y));
	    if(Q2>1) OM_xbjQ1->Fill(float(xB));

	    if(evTrig&0x02)  OM_xbjMT  -> Fill(float(xB));
	    if(evTrig&0x04)  OM_xbjLT  -> Fill(float(xB));
	    if(evTrig&0x08)  OM_xbjOT  -> Fill(float(xB));
	    if(evTrig&0x10)  OM_xbjCT  -> Fill(float(xB));
	    if(evTrig&0x100) OM_xbjIMT -> Fill(float(xB));
	    if(evTrig&0x200) OM_xbjLAST-> Fill(float(xB));
	    if(evTrig&0x02)  OM_ybjMT  -> Fill(float(y ));
	    if(evTrig&0x04)  OM_ybjLT  -> Fill(float(y ));
	    if(evTrig&0x08)  OM_ybjOT  -> Fill(float(y ));
	    if(evTrig&0x10)  OM_ybjCT  -> Fill(float(y ));
	    if(evTrig&0x100) OM_ybjIMT -> Fill(float(y ));
	    if(evTrig&0x200) OM_ybjLAST-> Fill(float(y ));
	    
	    OM_BeamP ->Fill(float(bbP ));
	    OM_BeamPz->Fill(float(bbPz));
	    OM_BeamPx->Fill(float(bbPx));
	    OM_BeamPy->Fill(float(bbPy));
	    OM_BeamTx->Fill(float(bbTx));
	    OM_BeamTy->Fill(float(bbTy));

	    OM_Zvtx->Fill(float(zV));
	    OM_ZvsX->Fill(float(zV),float(xV));
	    OM_ZvsY->Fill(float(zV),float(yV));
	    if (hasMupBV) {OM_NTrkM->Fill(nTrksBV);}

	    for (iv = vrts.begin(); iv!=vrts.end(); iv++) {

	      //	      cout << "looping on secondaries " << vrts.size() << endl;

	      if ((*iv)->isPrimary()) continue;
	      CsVertex *vrt = *iv;
	      double Zs = vrt->getZ();                          // ***** V0 SEARCH *****
	      CsTrack* trk1 = trks.front(), *trk2 = trks.back();
	      if (trk1->getZones().size()==1 || trk2->getZones().size()==1)
		continue;  // Exclude fring-field tracks
	      const Hep3Vector &v1 = vtrks.front().Vec(), &v2 = vtrks.back().Vec();
	      Hep3Vector v0(v1); v0 += v2;
	      double pT = v1.perp(v0); if (pT<V0pTCut) continue;
	      double p02 = v0*v0;
	      double p12 = v1*v1, e1p = sqrt(p12+M2_p), e1pi = sqrt(p12+M2_pi);
	      double p22 = v2*v2, e2p = sqrt(p22+M2_p), e2pi = sqrt(p22+M2_pi);
	      double m_ppi  = sqrt((e1p+e2pi)*(e1p+e2pi)-p02);
	      double m_pip  = sqrt((e1pi+e2p)*(e1pi+e2p)-p02);
	      double m_pipi = sqrt((e1pi+e2pi)*(e1pi+e2pi)-p02);
	      bool isLambda  = fabs(m_ppi-M_Lam)<dm_Lam;
	      bool isALambda = fabs(m_pip-M_Lam)<dm_Lam;
	      bool isK0 = fabs(m_pipi-M_K0)<0.1;
	      //bool isK0 = fabs(m_pipi-M_K0)<dm_K0;
	      //	      cout << " m_pipi " << m_pipi << " dm " << dm_K0 << endl;
	      if (!(isLambda || isALambda || isK0)) continue;
	      double Xs = vrt->getX(), Ys = vrt->getY();
	      
	      if (zV>Zs) continue;
	      
	      HepVector vv(3); vv[0] = Xs-xV; vv[1] = Ys-yV; vv[2] = Zs-zV;
	      HepVector tmps(vv), tmpp(vv);
	      tmps = *(vrt->getCov(0))*tmps; tmpp = *(bestV->getCov(0))*tmpp;
	      double dist = vv.norm(), ddist = dot(tmps,vv)+dot(tmpp,vv);
	      ddist = sqrt(ddist)/fabs(dist); if (dist<V0vvCut*ddist) continue;
	      Hep3Vector vv3(Xs-xV,Ys-yV,Zs-zV);
	      double ctheta = vv3*v0/sqrt(p02)/dist;
	      if (ctheta>V0cthCut) {
		//	      if (isLambda)  hm_ppi->Fill(m_ppi);
		//	      if (isALambda) hm_ppi->Fill(m_pip);
		if (isK0)      OM_pipi->Fill(float(m_pipi));
	      }
	    }
	  }
	}
      }


      int mult=0;
      for (it++, musP = 0, mus = 0; it!=trks.end(); it++) {
	CsTrack *t = *it;       
      if (!bestV->getPar(t,par))
	CsErrLog::msg(elFatal,__FILE__,__LINE__,
		      "Evt #%d: Vertex has track w/o vertex parameters",
		      cse->getEventNumberInRun());
	tx = par(0); ty = par(1); a = sqrt(1+tx*tx+ty*ty);
	double tDCz = 1/a, tDCx = tx*bDCz, tDCy = ty*bDCz, tP = 1/fabs(par(2));
	double theta = acos(bDCx*tDCx+bDCy*tDCy+bDCz*tDCz);
	if ( specials.find(*it)!=specials.end() ) {
	} else {
	  double pz = tDCz*tP, px = tx*pz, py = ty*pz, Es = sqrt(tP*tP+M2_pi);
	  TLorentzVector ks(TVector3(px,py,pz),Es);
	  // hadrons
	  double xzh=Es/nu;
	  if(Q2>.1 && xB>0 && xB<1 && y>0 && y<1){
	    OM_zhad->Fill(float(xzh));
	    mult++;
	  }
	}
      }
      if(mult>0) OM_mult->Fill(float(mult+0.5));
    } 
  }
  
  
  //        ******************** STATISTICS ********************
  int nTrksBV = bestV ? bestV->getNTracks() : 0;
  statistics_[6] += nTrksBV; statistics_[7] += hasMupBV?nTrksBV:0;
  statistics_[9] += nPVrts; statistics_[11] += nPVrtsM;
  if (nPVrts) {
    statistics_[10]++; if (hasMupBV) statistics_[12]++;
    if (hasBMSBV) {
      statistics_[1]++; if (hasMupBV) {
	statistics_[2]++; if (hasBokBV) statistics_[8]++;
      }
    }
    if (nPVrts>1) statistics_[13]++;
  }

  dt -= _chronometer.inter(_chrono);    // time statistics
  if (dt<0) {
    _nfittingTime++; _totTime -= dt;
    statistics_[15] = (int)(1000 *_totTime / _nfittingTime);
  }
  return true;
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~        ReFinding        ~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int ReFinding(TMtx &Xn,TMtx &Cn, list<CsVTrack>& vtrks, list<CsTrack*>& trks)
{
  int NTrkReFind = 0;
  list<CsTrack*> Trks;
  
  //set tracks which were already proved by Kalman
  list<CsVTrack>::iterator ivt;
  for( ivt=vtrks.begin(); ivt!=vtrks.end(); ivt++ ) 
    Trks.push_back( (*ivt).getAssociatedTrk() );
  
  // finding of track which are close to the vertex
  list<CsTrack*>::iterator it;
  for( it=trks.begin(); it!=trks.end(); it++ ) {    
    CsTrack* trk = (*it);
    bool exist = false;
    for( ivt=vtrks.begin(); ivt!=vtrks.end(); ivt++ )
      if( trk == (*ivt).getAssociatedTrk() ) exist=true;
    
    if( !exist ) {
      const vector<CsHelix> &tpar = trk->getHelices();
      CsHelix hV;
      tpar[0].Extrapolate(10*Xn(1),hV);

      double x1=Xn(2), x2=hV(1)/10;
      double y1=Xn(3), y2=hV(2)/10;
      
      double cx1=Cn(2,2),    cy1=Cn(3,3),    cxy1=Cn(2,3);
      double cx2=hV(1,1)/10, cy2=hV(2,2)/10, cxy2=hV(1,2)/10;
      
      double dist = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) );
      
      // error on distance between two helices
      double sig = 1 / dist * sqrt( (x1-x2)*(x1-x2)*(cx1 +cx2 ) +
				    (y1-y2)*(y1-y2)*(cy1 +cy2 ) +
				  2*(x1-x2)*(y1-y2)*(cxy1+cxy2) );
      
      if( dist < 2*sig ) {
	NTrkReFind++; Trks.push_back( trk );
      }
    }
    
  }
  
  trks.clear();
  trks = Trks;
  return NTrkReFind;
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~      getBestVertex      ~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void CsKalmanFitting::getBestVertex(CsVertex *vrt, CsVertex *&bestV,
				    int &nonInteractingBV,
				    bool hasMup,      bool &hasMupBV,
				    bool hasBMS,      bool &hasBMSBV,
				    int piNCandidate, int &piNBV)
{
  bool printBV = Print_[0] || Print_[3];

  int nTrks = vrt->getNTracks(), vNDF = 2*nTrks-3;
  double vChi2 = vrt->getChi2(); vChi2 /= vNDF;

  static int nTrksBV; static double chi2BV;
  if (bestV) {
    nTrksBV = bestV->getNTracks(); chi2BV = bestV->getChi2()/(2*nTrksBV-3);
  }

  if (!bestV ||
      // Apply best vertex' basic definition: Largest #Tracks first and chi2
      // as 2nd criterion, to raise ambiguities. Plus some liability for
      // special cases.
      (nTrks>nTrksBV || nTrks==nTrksBV && (vChi2<chi2BV+.1 || nTrks<=5))) {
    bool same = bestV && nTrks==nTrksBV && fabs(vChi2-chi2BV)<.11;
    int whichInteract = bestV && nTrks==nTrksBV && nTrks<=5 ? -1 : 0;
    if (whichInteract) {          // Suspecting non-interacting track
      const list<CsTrack*> &trksBV = bestV->getTracks();
      Cs3Vector parbBV, partBV, parb, part; bool ok;
      ok = bestV->getPar(trksBV.front(),parbBV);
      int q = parbBV(2)>0 ? 1 : -1;
      ok &= bestV->getPar(trksBV.back(),partBV);
      const list<CsTrack*> &trks = vrt->getTracks();
      ok &= vrt->getPar(trks.front(),   parb);
      ok &= vrt->getPar(trks.back() ,   part);
      int idxBV = nTrks, idx = nTrks; // Indices of scattered particle
      if (ok && nTrks>=3) { // 2 2ndaries: retain fastest, w/ sign = beam's
	list<CsTrack*>::const_iterator jt; Cs3Vector par;
	jt = trksBV.begin(); int j; for (jt++, j = 2; j<nTrks; jt++, j++) {
	  CsTrack *t = *jt; if (bestV->getPar(t,par)) {
	    double r = par(2)*q, rp = partBV(2)*q;
	    if (r>0 && (rp<0 || r<rp)) { partBV = par; idxBV = j; }
	  }
	  else { ok = false; break; }
	}
	jt = trks.begin();          for (jt++, j = 2; j<nTrks; jt++, j++) {
	  CsTrack *t = *jt; if (vrt->getPar(t,par)) {
	    double r = par(2)*q, rp = part(2)*q;
	    if (r>0 && (rp<0 || r<rp)) { part = par;   idx = j; }
	  }
	  else { ok = false; break; }
	}
      }
      if (!ok) CsErrLog::msg(elFatal,__FILE__,__LINE__,
			     "Evt #%d: Vertex has track w/o vertex parameters",
			     CsEvent::Instance()->getEventNumberInRun());
      HepMatrix *covbBV, *covtBV, *covb, *covt;
      ok &= (covbBV = bestV->getCov(1))!=0;
      ok &= (covtBV = bestV->getCov(idxBV))!=0;
      ok &= (covb   = vrt  ->getCov(1))!=0;
      ok &= (covt   = vrt  ->getCov(idx))!=0;
      if (!ok) CsErrLog::msg(elFatal,__FILE__,__LINE__,
			     "Evt #%d: Vertex has track w/o covariant matrix",
			     CsEvent::Instance()->getEventNumberInRun());
      double drBV = q*(partBV(2)-parbBV(2)), d2rb, d2rt;
      d2rb = (*covbBV)(3,3), d2rt = (*covtBV)(3,3);
      double ddrBV = sqrt(d2rb+d2rt); // Neglecting correlations
      double dr =   q*(part  (2)-parb  (2));
      d2rb = (*covb  )(3,3), d2rt = (*covt  )(3,3);
      double ddr =   sqrt(d2rb+d2rt); // Neglecting correlations
      drBV /= ddrBV; dr /= ddr; // "r"(=c/P) conservation in units of sigma
      double thetaBV = 0, theta = 0;
      if      (drBV<4 && 4<dr) whichInteract = 1;     // Interact'g = current V
      else if (dr<4 && 4<drBV) whichInteract = 2;     //        ... = bestV
      else if (8<dr && 8<drBV) whichInteract = 3;     //        ... = both
      else {                                          //        ... = none
	// Checking the conservation of c/P should suffice to identify the
	// non-interacting beam track. But in practice, it is conditioned by
	// the quality of the BMS reco, which may have gone wrong. In that
	// case, the event is lost for a physics analysis point of view. But
	// it may still be interesting to identify the true interaction,
	// e.g. in order to evidence the bad behaviour of the BMS reco.
	// => Therefore, in case of ambiguity (e.g., in the case of 2 beam
	// tracks, one interacting and one not, the error in BMS reco being
	// affecting either the interacting beam (disguising it into a non-
	// interacting one), or vice-versa being affecting the non-
	// interacting one (disguising it into an interacting one): let's
	// check the scattering angle.
	double tx, ty, a, bDCX, bDCY, bDCZ, tDCX, tDCY, tDCZ;
	tx = parbBV(0); ty = parbBV(1); a = sqrt(1+tx*tx+ty*ty);
	bDCZ = 1/a; bDCX = tx/a; bDCY = ty/a;
	tx = partBV(0); ty = partBV(1); a = sqrt(1+tx*tx+ty*ty);
	tDCZ = 1/a; tDCX = tx/a; tDCY = ty/a;
	thetaBV = acos(tDCX*bDCX+ tDCY*bDCY+tDCZ*bDCZ);
	tx = parb  (0); ty = parb  (1); a = sqrt(1+tx*tx+ty*ty);
	bDCZ = 1/a; bDCX = tx/a; bDCY = ty/a;
	tx = part  (0); ty = part  (1); a = sqrt(1+tx*tx+ty*ty);
	tDCZ = 1/a; tDCX = tx/a; tDCY = ty/a;
	theta   = acos(tDCX*bDCX+tDCY*bDCY+tDCZ*bDCZ);
	if      (thetaBV<.002 && .002<theta)
	  whichInteract = 1;                          // Interact'g = currentV
	else if (.002<thetaBV && theta<.002)
	  whichInteract = 2;                          //        ... = bestV
	else if (thetaBV<.002 && theta<.002)
	  whichInteract = 0;                          //        ... = none
	else
	  whichInteract = 3;                          //        ... = both
      }
      if (printBV && whichInteract!=3) {
	printf("CsKF::doFitting==> 1/P transfer beam -> leading track: "
	       "best V(%.1fNDF) %+.2f sigma, "
	       "current V(%.1fNDF) %+.2f sigma\n",chi2BV,drBV,vChi2,dr);
	if (thetaBV) printf("CsKF::doFitting==> Scattering angle: "
			    "best V %.3f mrd, current V %.3f mrd\n",
			    thetaBV,theta);
	if (whichInteract==1) printf("  => Retain the latter.\n");
      }
      if (whichInteract==3 && vChi2<2 && chi2BV<2) {
	// Two, seemingly interacting, very good chi2 vertices: best chi2 is
	// then not so meaningful =>  Have to look for other criteria...
	// Let's try to reject a non interacting particle that would not
	// have been associated to its continuation in the spectrometer but
	// would have stolen instead a low energy track, e.g. a V0 decay
	// from the true interaction.
	double z = parb(2)/part(2), zBV = parbBV(2)/partBV(2);
	if (fabs(zBV)<.2 && .2<fabs(z) && z>0) {
	  whichInteract = 1;
	  if (printBV)
	    printf("CsKF::doFitting==> Momentum fraction of leading track: "
		   "best V(%.1fNDF) z=%+.2f, "
		   "current V(%.1fNDF) z=%+.2f => Retain the latter.\n",
		   chi2BV,zBV,vChi2,z);
	}
      }
    }
    if (whichInteract==3 && same) { // Suspecting same set of tracks
      const list<CsTrack*> trks = vrt->getTracks();
      const list<CsTrack*> trksBV = bestV->getTracks();
      list<CsTrack*>::const_iterator it, jt;
      for (it = trks.begin(); it!=trks.end(); it++) {
	CsTrack *t = *it; bool match;
	for (jt = trksBV.begin(), match = false; jt!=trksBV.end(); jt++) {
	  if (*jt==t) { match = true; break; }
	}
	if (!match) { same = false; break; }
      }
    }
    bool trivialCase = same ||// Same or 2 non-interact'g, good chi2 vertices...
      // The one w/ mu', if any, will be retained, because it can help to
      // have so an accounting of the trigger decision, while it does not
      // impact physics, and hence cannot bias it.
      whichInteract==0 && vChi2<2 && chi2BV<2;
    if (!bestV ||
	nTrks>nTrksBV ||
	// At this point, implicitly, nTrks==nTrksBV ambiguous case, cf.
	// supra. Then...
	// ...Replace former non-interacting, if lone, prioritarily.
	whichInteract==1 ||
	// ...Disregard current if, lone, non-interacting.
	whichInteract!=2 &&
	// ...Same track or otherwise trivial case: gave precedence to mu'.
	(trivialCase && hasMup && !hasMupBV || 
	 // ...Best chi2 raises ambiguity left over.
	 !trivialCase && vChi2<chi2BV)) {
      bestV = vrt;
      hasMupBV = hasMup; hasBMSBV = hasBMS; piNBV = piNCandidate;
      // => "whichInteract" !=0 only if newly declared "bestV" does interact
      whichInteract &= 1;
    }
    else
      // => "whichInteract" !=0 only if old "bestV" does interact
      whichInteract &= 2;
    if (whichInteract==0) nonInteractingBV = 1;
  }
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~     getBestVertexDY     ~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void CsKalmanFitting::getBestVertexDY(CsVertex *vrt, CsVertex *&bestV,
				      int &nGoodTrksBV, int &nMuswTBV,
				      double &tChi2BV, // Defined upon return only if nMuswTBV>0
				      bool hasMup,      bool &hasMupBV,
				      bool hasBMS,      bool &hasBMSBV,
				      int piNCandidate, int &piNBV)
{
  // - The selection of the BestPrimayVertex done by coral is only indicative:
  //  (almost) all the info needed to make a different, yet well documented,
  //  choice is availbale at PHAST level. In addition, in the Drell-Yan case,
  //  our COMPASS foreseen strategy is to move all of the vertexing to PHAST.
  // - Reason why, I take here the liberty no to put too much thinking in this
  //  DY-dedicated BestVertex algorithm. Which is the following.
  //   I) Largest # of tracks, excluding obvious beam muons.
  //     (Note that it would be better to exclude those tracks at en earlier
  //     step.)
  //  II) If ambiguity: largest # of muon tracks w/ good timing (I have not
  //     retained the #muons criterion as my first choice, in order no to
  //     penalize deficient muID. Good timing is requested so as to make a
  //     smooth transition to criterion #3.)
  // III) If ambiguity, best overall space+time chi2, the relative contribution
  //     of timing being specified by option "CsKalmanFitting BpVTimeW".
  //      Time is restricted to muon tracks.
  //      Good timing is defined as "track w/ hodo hit assocated". The reason
  //     being that:
  //     i) Hodo provide good timing precision,
  //    ii) Detector groups have always tended to neglect timing. =>
  //     Better restrict ourselves to a small sets of detectors (hodos are such)
  //     which timing response can be easily be monitored and calibrated by
  //     end users.

  static int first = true; static unsigned int hodoPats[CSTRACK_MAPSIZE];
  if (first) {
    // To ease to checking that a given track has hodo hits associated
    first = false;
    for (int word = 0; word<CSTRACK_MAPSIZE; word++) hodoPats[word] = 0;
    list<CsDetector*> dets = CsGeom::Instance()->getDetectors();
    list<CsDetector*>::iterator Id; int idet;
    for (Id = dets.begin(), idet = 0; Id!=dets.end(); Id++, idet++) {
      if ((*Id)->GetTBName()[0]=='H') {
	int word = idet/32, bit = idet%32; hodoPats[word] |= 1<<bit;
      }
    }
  }

  int nGoodTrks, nMuswT; double tChi2;
  list<CsTrack*> trks =vrt->getTracks(); list<CsTrack*>::iterator it;
  it = trks.begin(); const CsTrack *trk = *it;
  if (!trk->getMeanTime()) {
    CsErrLog::msg(elError,__FILE__,__LINE__,
		  "Evt #%d: Vertex has beam track w/o timing",
		  CsEvent::Instance()->getEventNumberInRun());
    nGoodTrks=nMuswT = 0; tChi2 = -1;
  }
  else {
    const CsHelix &Hb = trk->getHelices().front();
    double Pb = fabs(1/Hb.getCop());
    double tb = trk->getMeanTime(), dTb = trk->getMeanTimeError();
    double s1, st;
    for (it++, nGoodTrks=nMuswT = 0, s1=st = 0; it!=trks.end(); it++) {
      const CsTrack *trk = *it;
      // Best would be here to retrieve muID performed by "PIDdoMuonID" and
      // stored in the CsParticle to which the CsTrack is associated, so that we
      // have a unique muID procedure. But this would require inverting the
      // the CsParticle->CsTrack association
      if (trk->getXX0()>15) nGoodTrks++;
      else {
	// Let's exclude beam muons, defined as P > Pbeam-20% (Pbeam is meant
	// to be the beamline and is taken from the beam attached to current
	// vertex for simplicity's sake), and paraxial (which is taken as
	// angle in x and y < 5mrd (noting that when it comes to cut on the
	// former, comes into play the incidence of the beam on target due to
	// the chicane when the target dipole is off: we neglect the effect)).
	const CsHelix &Ht = trk->getHelices().front();
	double Pt = fabs(1/Ht.getCop());
	if (fabs(Ht.getDXDZ())>.005 || fabs(Ht.getDYDZ())>.005 ||
	    Pt<Pb*.80) nGoodTrks++;
	continue;
      }
      if (!trk->hasMeanTime()) continue;
      const unsigned int *firedDets = trk->getFiredDetsBitmap();
      int word, hasHodo; for (word=hasHodo = 0; word<CSTRACK_MAPSIZE; word++) {
	if (firedDets[word]&hodoPats[word]) { hasHodo = 1; break; }
      }
      if (!hasHodo) continue;
      nMuswT++;
      double tt = trk->getMeanTime(), dTt = trk->getMeanTimeError();
      double w = 1/dTt/dTt; s1 += w; st += tt*w;
    }
    if (nMuswT) {
      double ts = st/s1, dTs2 = 1/s1;
      tChi2 = (ts-tb)*(ts-tb)/(dTb*dTb+dTs2);;
    }
    else tChi2 = -1;
  }

  if (bestV && nGoodTrks<nGoodTrksBV) return; // ***** GIVE UP FOR CRITERION (I)

  if (bestV) {
    if (nGoodTrks==nGoodTrksBV) {
      if      (nMuswT<nMuswTBV) return;      // ***** GIVE UP FOR CRITERION (II)
      else if (nMuswT==nMuswTBV) {
	int nTrks = vrt->getNTracks(), nTrksBV = bestV->getNTracks();
	int ndf = 2*nTrks-3, ndfBV = 2*nTrksBV-3;
	double vChi2 = vrt->getChi2()/ndf, vChi2BV = bestV->getChi2()/ndfBV;
	if (nMuswT) {
	  vChi2 += BpVTimeW_*tChi2/nMuswT; vChi2BV += BpVTimeW_*tChi2BV/nMuswTBV;
	}
	if (vChi2>vChi2BV) return;          // ***** GIVE UP FOR CRITERION (III)
      }
    }
  }

  bestV = vrt;

  // Port attributes of current vertex to "BV"-type attributes
  // (In principle, this is not usefule, since in the DY case, mu' is not
  // defined, BMS is not in operation and pi multiplicity meaningless. Yet,
  // to be on the safe side...)
  nGoodTrksBV = nGoodTrks; nMuswTBV = nMuswT; tChi2BV = tChi2;
  hasMupBV = hasMup; hasBMSBV = hasBMS; piNBV = piNCandidate;
 
}
