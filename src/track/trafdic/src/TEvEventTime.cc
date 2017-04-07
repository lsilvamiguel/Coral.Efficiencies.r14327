// $Id: TEvEventTime.cc 13354 2012-04-14 23:56:53Z ybedfer $

/*! 
   "SetEventTime": Derive "eventTime":
    I) From alternative master trigger time, when:
      i) defined,
     ii) consistent w/ at least one beam track time.
   II) From beam track time, if it can be uniquely defined (unique track).
  III) From beam track time determined in a previous pass of the track + vertex
      reconstruction (this latter case occuring if the event time turns out to
      differ very much from 0, cf. "CsAverPattern::doPattern").
  Event time will be used later on to update drift THits, provided it deviates
  from trigger time by more than "TOpt::dCut[71]",
   - W/ no other condition in cases (I) and (III).
   - Else depending upon trigger pattern matching "TOpt::ReMode[35,36,37]" (
    cf. also "TOpt::ReMode[40]" for the special case of high-Q2 trigger in
    earlier 2007).

   "UpdateDrifts":
   - Update CsClusters from drift-like detectors w/ "eventTime", prior to track
    finding in the spectrometer, and hence while the propagation time along the
    anode wires is not yet known. Update takes place only under certain
    requirements, cf. "SetEventTime" supra.
   - "eventTime" is stored into "eventTRef" and reset = 0.
   - Note: The propagation time is corrected for later on, cf. "FitSegments",
    "TracksFit2" and "TracksRefit":
     - The correction bear on THit's, as opposed to CsCluster's.
     - It's ported to CsCluster's depending upon option "ReMode[29]&0x4".

  \author  Yann.Bedfer@cern.ch
*/

#include "TH1.h"
#include "TH2.h"
#include "CsErrLog.h"
#include "CsInit.h"
#include "DaqDataDecoding/DaqOption.h"
#include "CsDetector.h"
#include "CsEvent.h"
#include "TEv.h"
#include "THit.h"

using namespace std;

/*
  All methods are specific to "lattice" alternative of "traffic".
 */

void TEv::SetEventTime()
{
  // ******************** (UPON 1ST CALL) INITIALISATIONS ********************

  static unsigned int hodoBasedTrigs = 0x0;  // Typically = IMLO|Incl.M, if they exist (DIS[+photoprod] data taking, cf. init block infra
  static unsigned int notCaloTrigs;
  static TH2D *hBeamTime = 0, *hBeamUniT = 0;
  static TH2D *hBeamMult = 0, *hBeamAltT = 0, *hBeamAltB = 0;
  static bool first = true; if (first) {     // ***** BOOK HISTOS *****
    first = false;
    if (TOpt::Hist[8]) {
      //             ***** HISTOGRAM BOOKING *****
      CsHistograms::SetCurrentPath("/Traffic");
      int nbinsY = 1;   // Bin 0 = inclusive
      for (int bit = 0; bit<16; bit++) if (1<<bit&TOpt::Hist[8]) nbinsY++;
      char hTitle[] = "Unique Beam Time vs. Trigger(0xffffffff)  ";
      sprintf(hTitle,"Beam Time vs. Trigger(0x%x)",TOpt::Hist[8]);
      hBeamTime = new TH2D("hBeamTime",hTitle,60,-15,15,nbinsY,-.5,nbinsY-.5);
      sprintf(hTitle,"Unique Beam Time vs. Trigger(0x%x)",TOpt::Hist[8]);
      hBeamUniT = new TH2D("hBeamUniT",hTitle,60,-15,15,nbinsY,-.5,nbinsY-.5);
      sprintf(hTitle,"Beam Mult vs. Trigger(0x%x)",TOpt::Hist[8]);
      hBeamMult = new TH2D("hBeamMult",hTitle,8,-.5,7.5, nbinsY,-.5,nbinsY-.5);
      hBeamAltT = new TH2D("hBeamAltT","Alter master trigger Time vs. beam reco",
			   40,-10,10,2,-.5,1.5);
      hBeamAltB = new TH2D("hBeamAltB","Alter master trigger vs. beam reco",
			   16,-0.5,15.5,2,-.5,1.5);
      CsHistograms::SetCurrentPath("/");

      // ***** GET TRIGGER BIT# of special interest for histogramming *****
      static unsigned int caloBasedTrigs = 0x0;
      if (!isMC) {
	const CS::DaqOption &opts = CsInit::Instance()->getDaqEventsManager().GetDaqOptions();
	int bit = -1;// Get Calo bit
	try { bit = opts.GetTriggerBit("CalorimeterTrigger"); } catch ( ... ) { }
	if (bit>=0) caloBasedTrigs |= 1<<bit;
	bit = -1;    // Get high-Q2T bit
	try { bit = opts.GetTriggerBit("LargeQ2Trigger"); } catch ( ... ) { }
	if (bit>=0) caloBasedTrigs |= 1<<bit;
	bit = -1;    // Get IMLO...
	try { bit = opts.GetTriggerBit("InnerTrigger"); } catch ( ... ) { }
	if (bit>=0) hodoBasedTrigs = 1<<bit;
	bit = -1;    // ...
	try { bit = opts.GetTriggerBit("MiddleTrigger"); } catch ( ... ) { }
	if (bit>=0) hodoBasedTrigs = 1<<bit;
	bit = -1;    // ...
	try { bit = opts.GetTriggerBit("LadderTrigger"); } catch ( ... ) { }
	if (bit>=0) hodoBasedTrigs = 1<<bit;
	bit = -1;    // ...
	try { bit = opts.GetTriggerBit("OuterTrigger"); } catch ( ... ) { }
	if (bit>=0) hodoBasedTrigs = 1<<bit;
	bit = -1;    // Inclusive Middle
	try { bit = opts.GetTriggerBit("InclMiddleTrigger"); } catch ( ... ) { }
	if (bit>=0) hodoBasedTrigs = 1<<bit;
	bit = -1;    // There may be an Incl.OT, at times,...
	try { bit = opts.GetTriggerBit("InclOuterTrigger"); } catch ( ... ) { }
	if (bit>=0) hodoBasedTrigs = 1<<bit;
      }
      notCaloTrigs = ~caloBasedTrigs;
    }
  }

  if (ptrEvt()->reTrackingON()) { // ***** 2ND ITERATION of TRACKING... *****
    eventTime =
      ptrEvt()->getReTrackT0(); // ***** ...EVENT'S TIME retrieved from 1ST ITER
    // The event time derived from previous, 1st iteration takes precedence
    // over all others. And histos already filled (in 1st iter). Therefore...
    return; // ...exit now.
  }

  //                       ***** GET TRIGGER *****
  const unsigned int allTrigs = 0xffff; // Although TCS can handle 0x7fffff, it was decided to limit max. #triggers to 16. And take advantage of the fact to make use of the other bits of the trigger TDC module.
  unsigned int evtTrig = trig_mask;
  evtTrig &= allTrigs;  // Cut away trailing end bits (i.e. online filter...)

  int altMasterTrig; static double altT0, extraTime;
  if (!isMC) {                        // ***** ALTERNATIVE MASTER TRIGGER *****
    // The master trigger (i.e. the trigger selected from the trigger pattern by
    // the decoding library) may not be optimum. In that case an alternative
    // trigger may be defined and used for serving as T0. Which means:
    // - Defining time gates: this is achieved in "TEv::ImportCluster"
    // - Defining T0 of drift detectors, which will be achieved be setting here
    //  the event's time equal to the alternative trigger's (and relying on
    //  "TEv::UpdateDrifts" to re-calculate the drift clusters).
    altMasterTrig = ptrEvt()->getAlterMasterTrigger();
    altT0 = ptrEvt()->getTriggerTimeCorr();
    extraTime = // Extra time set by "CsEvent::_clusterize" to make for alternative trigger (if exists) jitter
      ptrEvt()->getExtraTimeWidth();
  }
  else altMasterTrig = -1;

  //                                          ***** GET BEAM TIME *****
  //  I) Get time of beam track if unique
  // II) Check that a beam track (among possibly several ones) compatible w/
  //    alternative master trigger exists.
  int nBeams = listTrack.size(); double beamTime; int altBeam;
  list<TTrack>::iterator it; for (it = listTrack.begin(), altBeam = 0,
				    beamTime = 0; it!=listTrack.end(); it++) {
    TTrack &t = *it; t.UseHitTime();    if (t.SigmaTime<0) continue;
    double bT = t.MeanTime;
    if (nBeams==1) beamTime = bT;
    if (altMasterTrig>=0 && fabs(bT-altT0)<4*t.SigmaTime+extraTime)
      // Found at least one beam track compatible w/ alternative master trigger
      altBeam = 1;
    if (TOpt::Hist[8]) {
      hBeamTime->Fill(bT,0);
      unsigned int bit; int bin; for (bit=bin = 1; bit<0x10000;
				      bit <<= 1, bin++) {
	if (!(bit&TOpt::Hist[8])) continue; if (!(bit&evtTrig)) continue;
	if (bit==evtTrig ||
	    (bit&hodoBasedTrigs) &&// Hodo triggers. Allow them to have calo bit
	    bit==(evtTrig&notCaloTrigs)) hBeamTime->Fill(bT,bin);
      }
    }
  }

  if (TOpt::Hist[8]) {                    // ***** FILL GLOBAL HISTOS *****
    hBeamMult->Fill(nBeams,0); if (nBeams==1) hBeamUniT->Fill(beamTime,0);
    unsigned int bit; int bin; for (bit=bin = 1; bit<0x10000;
				    bit <<= 1, bin++) {
      if (!(bit&TOpt::Hist[8])) continue; if (!(bit&evtTrig)) continue;
      if (bit==evtTrig ||
	  (bit&hodoBasedTrigs) && // Hodo triggers. Allow them to have calo bit
	  bit==(evtTrig&notCaloTrigs)) {
	hBeamMult->Fill(nBeams,bin);
	if (nBeams==1) hBeamUniT->Fill(beamTime,bin);
      }
    }
    hBeamAltB->Fill(altMasterTrig,altBeam);
    if (altMasterTrig>=0) hBeamAltT->Fill(altT0,altBeam);
  }

  //       ********** INIT EVENT TIME STRICTLY = 0. **********
  // This is expected to prevent any attempt at updating drift THits.
  eventTime = 0;

  //       ********** (CONDITIONAL) SET EVENT TIME **********
  if (isMC && TOpt::ReMode[4])       // ***** IDEAL PR
    // This can be overwritten infra: makes it a little less ideal...
    eventTime = CsEvent::Instance()->getTriggerMCOffset();
  if (nBeams==1)                     // ***** UNIQUE BEAM: EVENT'S TIME = BEAM'S
    eventTime = beamTime;
  if (altMasterTrig>=0)       // ***** ALTER MASTER TRIG: EVENT'S TIME = ALTER'S
    eventTime = altT0;
}

// ***************************************************************************
// ****************************** UpdateDrifts *******************************
// ***************************************************************************
void TEv::UpdateDrifts()
{
  const TSetup& setup = TSetup::Ref();

  //                       ***** GET TRIGGER *****
  const unsigned int allTrigs = 0xffff; // Although TCS can handle 0x7fffff, it was decided to limit max. #triggers to 16. And take advantage of the fact to make use of the other bits of the trigger TDC module.
  unsigned int evtTrig = trig_mask;
  evtTrig &= allTrigs;  // Cut away trailing end bits (i.e. online filter...)

  //     ***** REQUIRE THAT TRIGGER MATCHES ReMode[35]... *****
  unsigned int inclMask =// Trigger pattern for which bits we require only an...
    // ...inclusive match. This can be used, e.g., to allow in 2007 combinations
    // of CT (0x10) and highQ2T (0x200) and require only CT to be exclusive, w/
    // the following: ReMode[35] = 0x210 and ReMode[40] = 0x10 ( which could be
    // useful for those periods where highQ2T was not delayed).
    TOpt::ReMode[35]&(~TOpt::ReMode[40]);
  if ((( evtTrig&TOpt::ReMode[35])!=evtTrig && // Strict match 
       !(evtTrig&inclMask) ||                  // Inclusive match
       // The following can happen in MC: event did not pass the trigger...
       evtTrig==0) && // => Exclude the case. No possible trigger offset anyway.
      // ...OR THAT ALTERNATIVE MASTER TRIGGER is defined
      ptrEvt()->getAlterMasterTrigger()<0 &&
      // ...OR 2ND ITERATION of TRACKING GOING ON
      !ptrEvt()->reTrackingON() ||
      //    ***** ...AND THAT "eventTime" LARGE ENOUGH *****
      fabs(eventTime)<TOpt::dCut[71]) return;

  const vector<TDetect> &dets = setup.vDetect();
  const vector<TPlane> &planes = setup.vPlane();
  for (int i = 0; i<(int)setup.vDetect().size(); i++) {

    //                ********** LOOP OVER ALL DETS **********

    if (!planes[i].IFlag) continue;
    const CsDetector *csd = dets[i].PtrDet();
    if (!csd->hasDrift()) continue;
    const_cast<CsDetector*>(csd)->updateClusters(eventTime);
  }
  eventTRef = eventTime; eventTime = 0;
}
