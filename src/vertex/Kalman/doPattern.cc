/*!
  \file    doPattern.cc
  \brief   Vertex pattern recognition. May perform a re-fitting (of the tracks) or request a re-tracking.
  \author  Alexandre Korzenev
  \version $Revision: 14094 $
  \date    $Date: 2015-11-06 16:28:48 +0100 (Fri, 06 Nov 2015) $

  Selection:
  - w/ momentum P (can be granted by fringe field).
  - zFirst upstream of SM2 (to remove tracks assigned a P by the bridging over
   the muonFilter #2).
  - P > Pmin (from option).
  - Special tracks (mu' and beam) are exempted from these requirements: if
   required, have to be checked by the method granting the "special" attribute.

  Re-fitting and re-tracking:
  - They may be required if the interaction time departs significantly from
  the trigger time (because the T0 used in the tracking for drift detectors is
  then significantly wrong). This can happen if trigger time resolution is bad,
  e.g. CT or RPD trigger.
  - Interaction time can be reasonably well determined by the beam track's
  time (because of scifis and/or silicons). If there are several accidentally
  coincident beam tracks, the ambiguity can only raised by selecting a best
  vertex. Which can be done, w/ some reasonable reliability, at the PR level.
  - Re-fitting is executed from w/in the "doPattern" method. Re-tracking is
  only requested upon return to the calling routine.
*/
#include "TH1.h"
#include "CsInit.h"
#include "CsStopwatch.h"
#include "CsErrLog.h"
#include "CsGeom.h"
#include "CsEvent.h"
#include "CsMCUtils.h"
#include "DaqDataDecoding/DaqOption.h"
#include "THlx.h"
#include "CsTrafficRefitting.h"
#include "CsAverPattern.h"
#include "CsBeam.h"

using namespace std;

bool CsAverPattern::doPattern(vector<CsParticle*> &parts, // Input particles. 
			      double *reTrackT0) // If non null on input, "reTrackT0" allows "doPattern" to request a re-tracking, by returning on output the event time to be used in the re-tracking. The request is granted when the beam track time of the best vertex turns out to be large (w.r.t. cut suplied by option). Also required is that the # of beam tracks be > 1. Conversely, a null input means that one is precisely re-vertexing after a re-tracking has been performed (which is also flagged by CsEvent's attribute "reTrackingON"); this is taken into account in the flowchart: avoid double counting in the statistics book-keeping and the histogramming (but then only the status of the initial tracking and vertexing is recorded) and abstain from refitting.
{
  static CsStopwatch _chronometer(3);
  static int _chronoTOT = _chronometer.start();
  static int _chronoPR  = _chronometer.start();
  static int _chronoSC  = _chronometer.start();
  static int    _nTotTime = 0,   _nPrTime  = 0,  _nScTime  = 0;
  static double _totTimeTOT = 0, _totTimePR = 0, _totTimeSC = 0;
  double dtTOT = _chronometer.inter(_chronoTOT);

  bool hadronRun = CsInit::Instance()->IsAHadronJob();

  static TH1D *hNEVrt,*hNEVrtM;
  static TH1D *hNBeams, *hNBeams_SI, *hNBeams_C, *hNBeams_O;
  static TH1D *hNBwBMS, *hNBwBMS_SI, *hNBwBMS_C, *hNBwBMS_O;
  static double zTarget, zSM2;  // in mm
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
  static unsigned int pureTrigs, specialTrigs, stdTrigs = 0x7;

  static bool first = true; if (first) {
    //  **********************************************************************
    //  *************************   INITIALIZATIONS  *************************
    //  **********************************************************************
    first = false;

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

    //     ********** GET SETUP FROM CsGeom **********
    CsGeom *g = CsGeom::Instance();
    CsField *f = g->getCsField(); CsMagInfo *mag = f->getMagInfo();
    int nMags = f->getNumOfMags();
    if (nMags!=3) CsErrLog::msg(elFatal,__FILE__,__LINE__,
				"# of magnets = %d, while expected = 3",nMags);
    zSM2 = mag[2].zcm; // SM2 = 3rd magnet, whatever the setup, hadron or muon.
    zTarget = g->getTargetCenter();

    if (hist_) {
      //       *************** BOOK HISTOGRAMS ***************
      CsHistograms::SetCurrentPath("/CsAverPattern/Distribution");
      hNEVrt = new TH1D("hNEVrt", "N of Evts w/ pVertex before fit",  8,0,4 );
      hNEVrtM= new TH1D("hNEVrtM","N of Evts w/ #mu/#mu' before fit", 8,0,4 );
      char hN[]  = "hNBwBMS_SI"; char hT[] = "#Beams w/ BMS - PrimakoffT   ";
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
      TH1D **hs[4][2] = { {&hNBeams,  &hNBwBMS},   {&hNBeams_SI,&hNBwBMS_SI},
			  {&hNBeams_C,&hNBwBMS_C}, {&hNBeams_O, &hNBwBMS_O} };
      for (int j = 0; j<4; j++) {
	sprintf(hN,"hNBeams%s",tabs[j]); sprintf(hT,"#Beams%s",tags[i][j]);
	*hs[j][0] = new TH1D(hN,hT,8,0,8);
	sprintf(hN,"hNBwBMS%s",tabs[j]); sprintf(hT,"#Beams w/ BMS%s",tags[i][j]);
	*hs[j][1] = new TH1D(hN,hT,8,0,8);
      }
      CsHistograms::SetCurrentPath("/");
    }
  }

  //         *************** INIT VERTEX PATTERN ***************
  TMtx Xn(3),Cn(3,3);
  CsEvent *event  = CsEvent::Instance();

  bool reVertex =// Will bypass (some) statistics book-keeping and histogramming
    event->reTrackingON();
  if (reTrackT0) { // ***** REQUESTING a RE-TRACKING is ENABLED *****
    // (Cf. supra the commenting the list of arguments of "doPattern".)
    *reTrackT0 = 0; // Set the value to be return on out
    if (reVertex) CsErrLog::mes(elFatal,
      "Requesting re-tracking while re-tracking is already going on!");
  }

  const int NTRK_MAX = 32;
  int vIact[NTRK_MAX]; CsTrack* TrkRef[NTRK_MAX]; THlx* HTrMom[NTRK_MAX];
  for (int i = 0; i<NTRK_MAX; i++) {
    HTrMom[i] = 0; TrkRef[i] = 0;
    vIact[i] = 0; // all tracks are accepted
  }

  vrts_.clear();  // Clean vertex list from previous call

  list<CsTrack*> Tracks, Beams, Specials; int nBwBMS;
  vector<CsParticle*>::iterator ip;
  for (ip =  parts.begin(), specials_.clear(), nBwBMS = 0; ip!=parts.end();
       ip++ ) {

    // ********** LOOP ON CsParticles -> "Beams","Specials","Tracks" ***********

    CsTrack *trk  = const_cast<CsTrack*>( (*ip)->getTrack() );
    if (trk==0) continue;
    const vector<CsHelix> &tpar = trk->getHelices();
    
    if ((*ip)->getType()== CsParticle::SPECIAL) {     // ***** BEAM or MU' *****
      if (tpar[0].getCop()!= 0) {
	specials_[ trk ] = true;
	if (tpar[0].getZ()<zTarget ) {                          // -> "Beams"
	  Beams.push_back(trk);
	  const list<CsZone*> &zones = trk->getZones();
	  list<CsZone*>::const_iterator iz;
	  for (iz = zones.begin(); iz!=zones.end(); iz++)
	    if((*iz)->getName()=="BMSzone") { nBwBMS++; break; }
	}
	else                          Specials.push_back(trk);  // -> "Specials"
      }
      else if(Print_[1])
	cout<<"CsAP::doPattern==>ERROR: Special track q/p="
	    <<tpar[0].getCop()<<". Z="<<tpar[0].getZ()<<endl;
    }
    else {                                         // ***** ORDINARY TRACK *****
      double zFirst = tpar[0].getZ();
      if (zFirst>zTarget &&// Reject beams.
	  //                  Reject tracks w/ P yet starting downstream of SM2.
	  // (Tracks bridging over muFilter #2 can be assigned a P (= beam's P,
	  // e.g. in "TEv::TracksFit2").)
	  zFirst<zSM2)
	Tracks.push_back(trk);                                  // -> "Tracks"
    }
  }
  if (Print_[0]) cout<<"CsAP::doPattern==> #beams="<<
		   Beams.size()<<" #mu'="<<Specials.size()<<endl;

  if (Tracks.size()+Specials.size()==0) {
    if(Print_[0]) cout<<"CsAP::doPattern==> No reco'ed spectrometer tracks\n";
    return 0;
  }
  int N_Specials = Specials.size();  // #muons (incident|scattered) in current
  if (!reVertex) {
    if (Beams.size()>1) statistics_[13]++;   // Stats: # of events with #beams>1
    if (N_Specials  >1) statistics_[15]++;   // Stats: # of events with #mu'>1
  }

  int tn;       // *************** TRACK INDEX: 1ST IS BEAM ***************
  if   (findPrim_) tn=1;
  else             tn=0;
  list<CsTrack*>::iterator it;
  for (it = Tracks.begin(); it!=Tracks.end(); it++) {
    if (tn>=NTRK_MAX) break;

    //       *********** LOOP ON ORDINARY TRACKS: ... **********
    //       ...SELECTION: REQUIRE MOMENUM P, P > P_min.
    //       => FILL WORKING ARRAYS
    const vector<CsHelix> &tpar = (*it)->getHelices();
    if (tpar[0].getCop()!=0 && fabs(1/tpar[0].getCop())>MomCut_) {
      if (Print_[0])
	printf("CsAP::doPattern==> Added: Mom =%#7.4g, z=%#7.4g cm\n",
	       1/tpar[0].getCop(),tpar[0].getZ()/10);
      HTrMom[tn] = new THlx();
      HTrMom[tn]->ImportHelix(tpar[0]); TrkRef[tn] = *it;
      tn++;
    }
  }

  T0SettingVertex_ = 0; // CsAverpattern's pointer to, preliminary, Best pVertex
  if (findPrim_) {
    double dtPR = _chronometer.inter(_chronoPR);

    // *********************** PRIMARY VERTEX ***********************

    //                                  ***** STATISTICS and HISTOGRAMMING *****
    bool pVertexWMu = false;
    if (!reVertex) {
      statistics_[0]++;                      // Stats: (vertexing) Was called
      if (Beams.size()==0) statistics_[2]++; // Stats: Beam was absent
    }
    const unsigned int allTrigs = 0xffff; // Although TCS can handle 0x7fffff, it was decided to limit max. #triggers to 16. And take advantage of the fact to make use of the other bits of the trigger TDC module.
    unsigned int evTrig = event->getTriggerMask();
    evTrig &= allTrigs;        // Cut away trailing end bits, i.e. online filter
    //#define doPattern_PRINT
#ifdef doPattern_PRINT // Some print-out for debugging purposes
    static int NBeams = 0, NBwBMS = 0, NEvts_Beam = 0, NEvts_BwBMS = 0;
    if (Beams.size() && !reVertex) {
      NBeams += Beams.size(); NBwBMS += nBwBMS;
      NEvts_Beam++; if (nBwBMS) NEvts_BwBMS++;
      printf("#B %d(%d)  #B,#E %d(%d) %d(%d)\n",
	     (int)Beams.size(),nBwBMS,NBeams,NBwBMS,NEvts_Beam,NEvts_BwBMS);
    }
#endif
    if (!reVertex) {
      hNBeams->Fill(Beams.size()); hNBwBMS->Fill(nBwBMS);
      if      (evTrig&stdTrigs) {        // Single out ``standard'' triggers...
	hNBeams_SI->Fill(Beams.size()); hNBwBMS_SI->Fill(nBwBMS);
      }
      if ((evTrig&pureTrigs)==evTrig) {  // ...and CT (+? highQ2) or beam
	hNBeams_C->Fill(Beams.size());  hNBwBMS_C->Fill(nBwBMS);
      }
      if (evTrig&specialTrigs) {         // ...and OT or NIT or DVCS
	hNBeams_O->Fill(Beams.size());  hNBwBMS_O->Fill(nBwBMS);
      }
      if (NSpec_ && N_Specials == 0)
	statistics_[3]++;                    // Stats: Not enough special tracks
    }

    if (!NSpec_ ||         // If Special (e.g. mu') not mandatory...
	N_Specials!= 0) {  // ... or available
      list<CsTrack*>::const_iterator iB; int nB;
      for (iB = Beams.begin(), nB = 0; iB!=Beams.end(); iB++, nB++) {
      
	// *************** LOOP OVER BEAMS for FINDING PRIMARY ***************

	//                   ***** INSERT CURRENT BEAM into WORKING ARRAYS *****
	const vector<CsHelix> &tpar= (*iB)->getHelices();
	if (Print_[0]) {
	  const CsBeam *csb = dynamic_cast<CsBeam*>(*iB);
	  char flag = csb && csb->getChi2CutFlag() ? '?' : '\0';
	  printf("\nCsAP::doPattern==> Added: Mom =%#7.4g, z=%#7.4g cm, t=%#7.4g ns   (beam%c)\n",
		 1/tpar[0].getCop(),tpar[0].getZ()/10,(*iB)->hasMeanTime()?(*iB)->getMeanTime():1000,flag);
	}
	if (HTrMom[0]!=0) delete HTrMom[0]; HTrMom[0] = new THlx();
	HTrMom[0]->ImportHelix(tpar[0]); TrkRef[0]= *iB; vIact[0] = 5;
	if (TrkRef[0]==0)
	  // This condition should never be fulfilled. Nevertheless, the check
	  // had been implemented, for some unfortunately now forgotten, reason:
	  //if (Print_[1]) cout<<"\nCsAP::doPattern==> No beam track.\n\n";
	  //continue;
	  // => Let's keep it, but making it now a fatal error, in order to be
	  //   able to verify that it indeed never happens.
	  CsErrLog::mes(elFatal,"Inconsistency: Iterator on beam tracks points to NULL");

	// ***** STRATEGY for HANDLING MU'
	// - Mu' is given privilege in the vertex PR (so that pVertex is kind of
	//  built upon mu').
	//  I) All candidate mu's are first tried.
	// II) If none's successful (possibly because no candidate mu' at all,
	//    in e.g., the hadron case (which latter case is technically singled
	//    out by an explicit "CsKalmanFitting Specials" option), one tries
	//    "FindPrimary" w/o any mu'.
	list<CsTrack*>::const_iterator iS, jS; int nS, nVs;
	for (iS = Specials.begin(), nS=nVs = 0; iS!=Specials.end();
	     iS++, nS++) {
	  if (Print_[1])
	    printf("\nCsAP::doPattern==> beam %d/%zu, mu' #%d/%d\n",
		   nB,Beams.size(),nS,N_Specials);

	  // ********** LOOP ON MU's: FIND PRIMARY w/ CURRENT MU' **********

	  int TN = tn;       // ***** TRACK INDEX FOR CURRENT (beam,mu') *****

	  //                                ***** CURRENT MU' IS SPECIAL *****
	  //   ***** ALL OTHER MU's -> WORKING ARRAYS AS ORDINARY Tracks *****
	  for (jS = Specials.begin(); jS!=Specials.end(); jS++) {
	    if (jS==iS) continue; // Skip current mu'
	    const vector<CsHelix> &tpar = (*jS)->getHelices();
	    //       ...SELECTION: REQUIRE MOMENUM P, P > P_min (as supra).
	    if (tpar[0].getCop()!=0 && fabs(1/tpar[0].getCop())>MomCut_) {
	      if (Print_[0])
		cout<<setprecision(4)<<"CsAP::doPattern==> Added: Mom = "
		    <<setw(6)<<1/tpar[0].getCop()
		    <<", z="<<tpar[0].getZ()/10<<" cm  (mu)\n";
	      if (TN>=NTRK_MAX-1) { // Reserve one component for mu'
		CsErrLog::mes(elError,
			      "Primary Vertex: too many tracks => Skip special!");
		break;
	      }
	      if (HTrMom[TN]!=0) delete HTrMom[TN];
	      HTrMom[TN] = new THlx(); HTrMom[TN]->ImportHelix(tpar[0]);
	      TrkRef[TN] = *jS; vIact[TN] = 1; TN++;
	    }
	  }

	  //        ***** CURRENT MU' -> WORKING ARRAYS AS SPECIAL TRACK *****
	  if (TN>=NTRK_MAX) {
	    CsErrLog::mes(elError,
	      "Primary Vertex: too many tracks => Skip mu' & give up!");
	    continue;
	  }
	  const vector<CsHelix> &tpar = (*iS)->getHelices();
	  if(Print_[0])
	    printf("CsAP::doPattern==> Added: Mom =%#7.4g, z=%#7.4g cm, t=%#7.4g ns   (mu')\n",
		   1/tpar[0].getCop(),tpar[0].getZ()/10,(*iS)->hasMeanTime()?(*iS)->getMeanTime():1000);
	  if (HTrMom[TN]!=0) delete HTrMom[TN]; HTrMom[TN] = new THlx();
	  HTrMom[TN]->ImportHelix(tpar[0]); TrkRef[TN] = *iS; vIact[TN] = 5;
	  TN++;

	  if (!FindPrimary(HTrMom,vIact,TrkRef,TN)) {
	    if (Print_[1])
	      printf("CsAP::doPattern==> beam %d/%zu, mu' %d/%d: "
		     "No pVertex.\n\n",nB,Beams.size(),nS,N_Specials);
	  }
	  else {
	    if (!reVertex) statistics_[1]++; // Stats: total # of pVertices
	    pVertexWMu = true; nVs++;
	    if (Print_[1])
	      printf("\nCsAP::doPattern==> beam %d/%zu, mu' %d/%d: "
		     "pVertex FOUND!\n\n",nB,Beams.size(),nS,N_Specials);
	  }
	}  // End loop on mu's

	if (!nVs) {// ********** NO pVERTEX w/ MU' (be it no mu' at all)...
	  //                          => FIND PRIMARY w/o MU'   **********
	  if (Print_[1]) cout<<"CsAP::doPattern==> NO mu'.\n";
	  if (!FindPrimary(HTrMom,vIact,TrkRef,tn)) {
	    if (Print_[1])
	      printf("CsAP::doPattern==> beam %d/%zu: No pVertex.\n",
		     nB,Beams.size());
	  }
	  else {
	    if (!reVertex) statistics_[1]++; // Stats: total # of pVertices
	    if (Print_[1]) printf("CsAP::doPattern==> beam %d/%zu:"
	      " pVertex is FOUND!\n",nB,Beams.size());
	  }
	}  // End find primary w/o mu'
      }  // End loop on beams
    }
    else if (!NSpec_ && Print_[1])
      cout<<"CsAP::doPattern==> NO special tracks.\n\n";
    
    //   ************* HAVE TO RE-FIT or RE-BUILD TRACKS? **********
    // - The reference time used for the processing of the current event (which
    //  is typically initalised to be trigger time) may turn out to be badly
    //  off. This, one can determine at this point of the vertex reco, where
    //  a definition of the best vertex, which embodies the interaction at the
    //  origin of the event, is available.
    // - Best definition of best vertex is vertex w/ largest number of tracks.
    //  The chi2 value (a crude preliminary one at this stage) is used to raise
    //  ambiguities.
    // - The time of the interaction (and hence most appropraite reference time
    //  for the processing) is taken from the best vertex' beam track time. For
    //  beam track time is (a priori) precisely determined. (One could think of
    //  using the average time of all tracks in the best vertex in a later
    //  version of the software.)
    if (refitTracks_ ||    // Refitting tracks allowed (by option)
	rebuildTracks_ &&  // Requesting a re-tracking allowed by option...
	reTrackT0) {       // ...and via argument for this call to "doPattern"

      //       ********** DETERMINE BEST PRIMARY VERTEX **********
      CsVertex *bestVertex; list<CsVertex*>::iterator iv;
      static int nTsBest; static double chi2Best;
      for (iv = vrts_.begin(), bestVertex = 0; iv!=vrts_.end(); iv++) {
	CsVertex *v = *iv;
	if (!v->isPrimary()) CsErrLog::mes(elFatal,
	  "Secondary vertex found while search for secondaries still pending!");
	CsTrack *beam = v->getTracks().front();
	int nTs = v->getNTracks(); double chi2 = v->getChi2()/nTs;
	if (chi2<0) {// Current pVertex is tagged as violating P conservation...
	  // ...cf. "CsAverPattern::FindPrimary", which means its set of tracks
	  // includes unduly associated tracks and its #tracks is overestimated.
	  // => Let's decrease that # by one unit. Which may still leave it
	  // overestimated. But must solves most cases. Indeed, that a fake
	  // pVertex has a large #tracks can only be due to the fact that its
	  // incident particle is close to the true pVertex's and picks up true
	  // secondaries. And that it has an even larger #tracks than the true
	  // pVertex (or equal # w/ better chi2) can only be due to either
	  // tertiaries that turn out to better match the wrong incident track
	  // or a bad reco (due in turn to, e.g., a bad T0). In any case, the
	  // pool of available tracks is the set of tracks of the true pVertex,
	  // and it must be very rare that the # of them picked up by the fake
	  // pVertex be larger by more than one unit than the # picked up by the
	  // true pVertex.
	  nTs -= 1; chi2 = fabs(chi2); // Restore the true (estimate of) chi2
	}
	int socId = beam->getAssociatedTrack();
	if (socId>=0) {  // Suspecting a non-interacting beam...
	  // ...The associate is presumed to be a mere continuation of the beam.
	  bool continuation = hadronRun; if (!hadronRun) {
	    // Association was set in tracking. There, in the muon case, beam
	    // has no precise momentum yet. And the check that we're not facing
	    // a Q^2=0 event could not be accurately check. Let's do it now.
	    const list<CsZone*> &zones = beam->getZones();
	    list<CsZone*>::const_iterator iz;
	    for (iz = zones.begin(); iz!=zones.end(); iz++) {
	      if ((*iz)->getName()!="BMSzone") continue;
	      const CsHelix &h0 = beam->getHelices()[0];
	      double r0 = h0.getCop(), dr20 = h0.getCov()[14];
	      const CsTrack *soc; list<CsTrack*>::iterator it;
	      list<CsTrack*> &tracks = event->getTracks();
	      for (it = tracks.begin(), soc = 0; it!=tracks.end(); it++) {
		const CsTrack *cst = *it;
		if (cst->getId()==socId) { soc = cst; break; }
	      }
	      if (!soc) CsErrLog::msg(elFatal,__FILE__,__LINE__,
  "CsTrack #%d's associate(=%d) not in event's list!",beam->getId(),socId);
	      const CsHelix &hs = soc->getHelices()[0];
	      double rs = hs.getCop(), dr2s = hs.getCov()[14];
	      if (fabs(r0-rs)<5*sqrt(dr20+dr2s)) continuation = true;
	      break;
	    }
	  }
	  if (continuation)   // If non-interacting indeed...
	    // ...Several actions could be taken. Most legitimate would be to
	    // require the continuation to be included in the vertex: too late.
	    // Here we decide to downgrade the corresponding vertex. A decrement
	    // of 1 should be enough, given that a non-interacting beam taking
	    // precedence over the genuine one (most probably) corresponds to
	    // cases where the 2 beams ``intersect'' and have then same #tracks.
	    // Yet, cf. #245 in 2006.03 D*.
	    nTs -= 3;
	}
	if (bestVertex==0 || nTs>nTsBest ||   // ***** BEST = LARGEST #TRACKS...
	    nTs==nTsBest && chi2<chi2Best) {
	  bestVertex = v; nTsBest = nTs; chi2Best = chi2;
	}
      }
      if (bestVertex) {
	CsTrack *beam = bestVertex->getTracks().front();
	if (event->isAMonteCarloEvent()) {
	  const CsMCTrack *mcb; if ((mcb = beam->getAssociatedMCTrack())) {
	    if (mcb!=event->getMCTracks().front()) {
	      cout<<"CsAP::doPattern: Best vertex beam != genuine MC beam\n";
	    }
	  }
	}
	if (rebuildTracks_ && reTrackT0) {// Requesting a re-tracking allowed...
	  if (Beams.size()>1 && // ...Require #beams >1, otherwise the true event time must have already been properly set in TraFDic
	      beam->hasMeanTime()) {
	    double t = beam->getMeanTime(), dt = beam->getMeanTimeError();
	    // ...Check beam track time against re-tracking cut.
	    // (We assume here that the time reference is zero...)
	    if (t+3*dt<-rebuildTracksCut_ || rebuildTracksCut_<t-3*dt ||
		// Rescue: the uncertainty on the beam time can be large in the
		// hadron case, where there may be no minimum requirement on the
		// # of scifi hits (and particularly in 2009 where, because of
		// the early triggering bug of the RPD-based trigger, the
		// interaction time may fall well outside the scifi time gate)
		// => Let's rescue those cases where the beam time fulfills the
		// cut +/- 5ns.
		t+5<-rebuildTracksCut_ || rebuildTracksCut_<t-5) {
	      *reTrackT0 = t; return 0; // ...W/in cut? exit current "doPattern"
	    }
	  }
	}
	//if (!reTrackT0): means that we are w/in the re-tracking phase...
	// ...(this is checked supra, cf. CsEvent's "reTrackingON"), in which
	// the event time is known from the start, and hence we need not refit.
	if (reTrackT0) {
	  static CsTrafficRefitting *refit = 0; // Instantiate refitting method
	  if (!refit) refit = new CsTrafficRefitting;
	  T0SettingVertex_ = bestVertex;
	  list<CsTrack*> tracksToBeDeleted; 
	  bool modified = refit->doRefitting(bestVertex,tracksToBeDeleted);
	  if (modified) { // Refitting had a non-null effect

	    // ********** UPDATE WORKING ARRAYS AND...               **********
	    // *****  ... RE-BUILD ALL VERTICES w/ BEST VERTEX' BEAM **********
	    if (Print_[0])
	      cout<<"\n\nCsAP::doPattern: ====== Rebuilding best vertex =====\n\n";
	    statistics_[9]++;                // # of events where pVertex refit

	    //                ***** INSERT BEST'S BEAM into WORKING ARRAYS *****
	    const vector<CsHelix> &tpar = beam->getHelices();
	    if (Print_[0]) {
	      CsBeam *csb = dynamic_cast<CsBeam*>(beam);
	      char flag = csb && csb->getChi2CutFlag() ? '?' : '\0';
	      printf("CsAP::doPattern==> Added: Mom =%#7.4g, z=%#7.4g cm, t=%#7.4g ns   (beam%c)\n",
		     1/tpar[0].getCop(),tpar[0].getZ()/10,beam->hasMeanTime()?beam->getMeanTime():1000,flag);
	    }
	    if (HTrMom[0]!=0) delete HTrMom[0]; HTrMom[0] = new THlx();
	    HTrMom[0]->ImportHelix(tpar[0]);
	    int bId = beam->getId(); list<CsTrack*>::iterator iB;
	    for (iB = Beams.begin(); iB!=Beams.end(); iB++)
	      if ((*iB)->getId()==bId) break;
	    if (iB==Beams.end()) CsErrLog::msg(elFatal,__FILE__,__LINE__,
	"Rebuiding best vertex: beam track (Id=%d) not in \"Beams\" list",bId);
	    TrkRef[0]= *iB; vIact[0] = 5;

	    //                  ***** DELETE ALL VERTICES WITH BEST'S BEAM *****
	    iv = vrts_.begin(); while (iv!=vrts_.end()) {
	      if ((*iv)->getTracks().front()->getId()==bId && *iv != bestVertex) {
		delete *iv;
		vrts_.erase(iv++);
		statistics_[1]--;            // Stats: total # of pVertices
	      }
	      else iv++;
	    }

	    if (tracksToBeDeleted.size()!=0) {
	      // ********** SOME TRACKS ARE TO BE DELETED **********

	      //       ***** DELETE ALL VERTICES WITH TRACKS TO BE DELETED *****
	      iv = vrts_.begin(); while (iv!=vrts_.end()) {
		CsVertex *v = *iv; list<CsTrack*> vTracks = v->getTracks();
		bool toBeDeleted;
		for (it = vTracks.begin(), toBeDeleted = false;
		     it!=vTracks.end(); it++) {
		  int tId = (*it)->getId(); list<CsTrack*>::iterator itp;
		  for (itp = tracksToBeDeleted.begin();
		       itp!=tracksToBeDeleted.end(); itp++) {
		    if ((*itp)->getId()==tId) { toBeDeleted = true; break;}
		  }
		  if (toBeDeleted) break;
		}
		if (toBeDeleted) {
		  if (*iv == bestVertex) {
		    CsErrLog::msg(elError, __FILE__, __LINE__,
		                  "Best primary vertex removed after refitting of tracks. Different t0 used in vertex (0.0) and track (%.1f) fit.",
		                  beam->hasMeanTime()?beam->getMeanTime():0);

		    bestVertex = NULL;
		    T0SettingVertex_ = NULL;
		  }
		  delete *iv;
		  vrts_.erase(iv++);
		  statistics_[1]--;          // Stats: total # of pVertices
		}
		else iv++;
	      }

	      specials_.clear(); Specials.clear(); // ***** REASSESS SPECIALS...
	      Tracks.clear();                    // ...AND ORDINARY TRACKS *****
	      if (Print_[0])
		printf("\nCsAP::doPattern: Erasing CsParticles w/ CsTrack deleted by refit, viz. %zu out of %zu\n",tracksToBeDeleted.size(),parts.size());
	      // "parts" is a vector: be careful when erasing components:...
	      int jp = 0; // ...make use of an integer counter
	      ip = parts.begin();  while (jp<(int)parts.size()) {
		CsTrack *trk  = const_cast<CsTrack*>((*ip)->getTrack());
		if (trk==0) { ip++; jp++; continue; }
		int tId = trk->getId(); bool toBeDeleted;
		for (it = tracksToBeDeleted.begin(), toBeDeleted = false;
		     it!=tracksToBeDeleted.end(); it++) {
		  if ((*it)->getId()==tId) { toBeDeleted = true; break; }
		}
		if (toBeDeleted) {               // ***** ERASE CsParticle *****
		  if (Print_[0]) printf("CsAP::doPattern: Deleting CsParticle type %d, born to CsTrack #%d\n",(*ip)->getType(),trk->getId());
		  parts.erase(ip); continue;
		}
		const vector<CsHelix> &tpar = trk->getHelices();
		if ((*ip)->getType()==CsParticle::SPECIAL) { // Beam or mu'
		  if (tpar[0].getCop()!=0) {
		    specials_[ trk ] = true;
		    if (tpar[0].getZ()>zTarget) 
		      Specials.push_back(trk);                  // -> "Specials"
		  }
		  else if (Print_[1])
		    cout<<"CsAP::doPattern: ERROR: Special track q/p="
			<<tpar[0].getCop()<<". Z="<<tpar[0].getZ()<<endl;
		}
		else if (tpar[0].getZ()>zTarget)            // Ordinary track
		  Tracks.push_back(trk);                          // -> "Tracks"
		ip++; jp++;
	      }  // End loop on CsParticle's
	      if (Print_[0]) cout<<"CsAP::doPattern: #beams="<<Beams.size()<<
			       " #mu'="<<Specials.size()<<endl;
	  
	      list<CsTrack*> eTracks = event->getTracks();
	      if (Print_[0])
		printf("\nCsAP::doPattern: Erasing CsTracks deleted by refit, viz. %zu out of %zu\n",tracksToBeDeleted.size(),eTracks.size());
	      it = eTracks.begin(); while (it!=eTracks.end()) {
		CsTrack *trk  = const_cast<CsTrack*>(*it);
		int tId = trk->getId(); bool toBeDeleted;
		list<CsTrack*>::iterator itp;
		for (itp = tracksToBeDeleted.begin(), toBeDeleted = false;
		     itp!=tracksToBeDeleted.end(); itp++) {
		  if ((*itp)->getId()==tId) {
		    toBeDeleted = true; tracksToBeDeleted.erase(itp); break;
		  }
		}
		if (toBeDeleted) {                  // ***** ERASE CsTrack *****
		  if (Print_[0]) printf("CsAP::doPattern: Deleting CsTrack #%d\n",
					trk->getId());
		  eTracks.erase(it++); continue;
		}
		else it++;
	      }

	      if (Tracks.size()+Specials.size()==0) {
		if(Print_[0])
		  cout<<"CsAP::doPattern==> No reco'ed spectrometer tracks\n";
		return 0;
	      }
	      int oldN = N_Specials; N_Specials = Specials.size();
	      if (N_Specials==0) {
		if (NSpec_) statistics_[3]++;// Stats: not enough special tracks
		if (oldN)  statistics_[15]--;// Stats: # of events with #mu'>1
	      }
	    }  // End "some tracks are to be deleted"

	    tn = 1; // Since tn==0 has been booked by beam
	    for (it = Tracks.begin(); it!=Tracks.end(); it++) {
	      if (tn>=NTRK_MAX) break;
	      //     ****** LOOP ON ORDINARY TRACKS: UPDATE WORKING ARRAYS *****
	      const vector<CsHelix> &tpar = (*it)->getHelices();
	      if (tpar[0].getCop()!=0 && fabs(1/tpar[0].getCop())>MomCut_) {
		if (Print_[0])
		  cout<<setprecision(4)<<"CsAP::doPattern==> Added: Mom = "
		      <<setw(6)<<1/tpar[0].getCop()
		      <<", z="<<tpar[0].getZ()/10<<" cm\n";
		if (HTrMom[tn]!=0) delete HTrMom[tn]; HTrMom[tn] = new THlx();
		HTrMom[tn]->ImportHelix(tpar[0]); TrkRef[tn] = *it;
		tn++;
	      }
	    }

	    // ***** STRATEGY for HANDLING MU' (cf. comment supra)
	    list<CsTrack*>::const_iterator iS, jS; int nS, nVs;
	    for (iS = Specials.begin(), nS=nVs = 0; iS!=Specials.end();
		 iS++, nS++) {
	      if (Print_[1])
		printf("\nCsAP::doPattern==> mu' #%d/%d\n",nS,N_Specials);

	      // ********** LOOP ON MU's: FINDING PRIMARY **********

	      int TN = tn;   // ***** TRACK INDEX FOR CURRENT (beam,mu') *****

	      //                            ***** CURRENT MU' IS SPECIAL *****
	      //                 ***** ALL OTHER MU's -> ORDINARY Tracks *****
	      for (jS = Specials.begin(); jS!=Specials.end(); jS++) {
		if (jS==iS) continue; // Skip current mu'
		const vector<CsHelix> &tpar = (*jS)->getHelices();
		if (tpar[0].getCop()!=0 &&
		    fabs(1/tpar[0].getCop())>MomCut_) {
		  if (Print_[0])
		    cout<<setprecision(4)<<"CsAP::doPattern==> Added: Mom = "
			<<setw(6)<<1/tpar[0].getCop()
			<<", z="<<tpar[0].getZ()/10<<" cm  (fomer special)\n";
		  if (TN>=NTRK_MAX-1) { // Reserve one component for mu'
		    CsErrLog::mes(elError,
		      "Primary Vertex: too many tracks => Skip special!");
		    break;
		  }
		  if (HTrMom[TN]!=0) delete HTrMom[TN];
		  HTrMom[TN] = new THlx(); HTrMom[TN]->ImportHelix(tpar[0]);
		  TrkRef[TN] = *jS; vIact[TN] = 1; TN++;
		}
	      }

	      //    ***** CURRENT MU' -> WORKING ARRAYS AS SPECIAL TRACK *****
	      if (TN>=NTRK_MAX) {
		CsErrLog::mes(elError,
		  "Primary Vertex: too many tracks => Skip mu' & give up!");
		continue;
	      }
	      const vector<CsHelix> &tpar = (*iS)->getHelices();
	      if(Print_[0])
		printf("\nCsAP::doPattern==> Added: Mom =%#7.4g, z=%#7.4g cm, t=%#7.4g ns   (mu')\n",
		       1/tpar[0].getCop(),tpar[0].getZ()/10,(*iS)->hasMeanTime()?(*iS)->getMeanTime():1000);
	      if (HTrMom[TN]!=0) delete HTrMom[TN]; HTrMom[TN] = new THlx();
	      HTrMom[TN]->ImportHelix(tpar[0]); TrkRef[TN] = *iS;
	      vIact[TN] = 5;
	      TN++;

	      if (!FindPrimary(HTrMom,vIact,TrkRef,TN)) {
		if (Print_[1])
		  printf("CsAP::doPattern==> mu' %d/%d: No pVertex!\n\n",
			 nS,N_Specials);
	      }
	      else {
		statistics_[1]++;        // Stats: total # of pVertices
		pVertexWMu = true; nVs++;
		if (Print_[1])
		  printf("\nCsAP::doPattern==> mu' %d/%d: pVertex FOUND!\n\n",
			 nS,N_Specials);
	      }
	    }  // End loop on mu's

	    if (!nVs) {// ********** NO pVERTEX w/ MU' (be it no mu' at all)...
	      //                           => FIND PRIMARY w/o MU'   **********
	      if (Print_[1]) cout<<"\nCsAP::doPattern==> NO mu'.\n\n";
	      if (!FindPrimary(HTrMom,vIact,TrkRef,tn)) {
		if (Print_[1])
		  cout<<"CsAP::doPattern==> No primary vertex.\n\n";
	      }
	      else {
		statistics_[1]++;            // Stats: total # of pVertices
		if (Print_[1])
		  cout<<"CsAP::doPattern==> pVertex is FOUND AGAIN!\n\n";
	      }
	    }  // End find primary w/o mu'
	  }  // End modified after refit
	}  // End refit
      }  // End "if pVertex"
    }  // End "refitTracks"

    if (vrts_.size()!=0) {                      // ***** STATISTICS *****
      statistics_[10]++;                     // Stats: # of events w/ pVertex
      if (pVertexWMu) statistics_[12]++;     // Stats: # of evts w/ mu/mu'
      if (vrts_.size()>1) statistics_[16]++; // Stats: # of evts w/ #pVertices>1
    }
    if (hist_) {                                // ***** HISTOGRAMMING *****
      hNEVrt->Fill(vrts_.size()?1:0); hNEVrtM->Fill(pVertexWMu?1:0);
    }
    dtPR -= _chronometer.inter(_chronoPR);      // ***** TIME STATISTICS *****
    if (dtPR<=0) {
      _nPrTime++; _totTimePR -= dtPR;
      statistics_[17] = (int)(1000 *_totTimePR/_nPrTime);
    }
    
  } ///////////// End of primary /////////////////

  if (findSec_) {
    double dtSC = _chronometer.inter(_chronoSC);
    
    // *********************** SECONDARY VERTICES ***********************

    //          ***** ADD MU'S TRACKS -> WORKING ARRAYS AS ORDINARY TRACKS *****
    list<CsTrack*>::iterator is;
    for (is = Specials.begin(); is!=Specials.end(); is++) {
      const vector<CsHelix> &tpar = (*is)->getHelices();
      if (tpar[0].getCop()!=0 && fabs(1/tpar[0].getCop())>MomCut_) {
	if (Print_[0]) cout<<setprecision(4)<<"CsAP::doPattern==> Added: Mom = "
			   <<setw(6)<<1/tpar[0].getCop()
			   <<", z="<<tpar[0].getZ()/10<<" cm  (fomer special, for secondary)"<<endl;

	if (tn>=NTRK_MAX) {
	  CsErrLog::
	    mes(elError,"2ndary Vertex: too many tracks => Skip special!");
	  break;
	}
	if (HTrMom[tn]!=0) delete HTrMom[tn]; HTrMom[tn] = new THlx();
	HTrMom[tn]->ImportHelix(tpar[0]); TrkRef[tn] = *is; vIact[tn] = 1;
	tn++;
      }
    }

    FindSecondaries( HTrMom, vIact, TrkRef, tn );

    dtSC -= _chronometer.inter(_chronoSC);      // time statistics
    if( dtSC <= 0 ) {
      _nScTime++; _totTimeSC -= dtSC;
      statistics_[18] = (int)(1000 *_totTimeSC / _nScTime);
    }
  } ///////////// end of secondary /////////////////

  for (int i = 0; i<NTRK_MAX; i++) if (HTrMom[i]!=0) delete HTrMom[i];

  dtTOT -= _chronometer.inter(_chronoTOT);    // time statistics
  if (dtTOT<=0) {
    _nTotTime++; _totTimeTOT -= dtTOT;
    statistics_[14] = (int)(1000 *_totTimeTOT / _nTotTime);
  }
  return true;
}
