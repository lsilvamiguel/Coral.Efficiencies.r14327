// $Id: CsTrafficPrepattern.cc 14069 2015-09-17 20:44:46Z lsilva $

/*!
   \file    CsTrafficPrepattern.cc
   \brief   Traffic Prepattern Derived Class.
   (Version specific to "lattice" alternative)
   \author  ...
   \version $Revision: 14069 $
   \date    $Date: 2015-09-17 22:44:46 +0200 (Thu, 17 Sep 2015) $ 

*/

/*
  Changes w/ respect to "traffic/CsTrafficPrepattern":
  i) The entire flowchart is modified: split into 2 cases:
    I) Standard 5 zones PR of physics data:
     - Import coral clusters in 2 steps: 1st, beam telescope.
     - Get event time from CsEvent or, by default, from beam-tracks's time if
      uniquely defined.
     - Upon option, update drift CsCluster's w/ event time.
     - Upon option, perfrom LR ambiguity raising.
     - Then only, import CsCluster's from the rest of the spectrometer.
   II) All other cases: import CsCluster's in one go.
  i) Call "DumpEvent"
  ii) "ImportClusters" w/ LR selection.
*/

#include "CsInit.h"
#include "CsErrLog.h"
#include "CsEvent.h"
#include "CsEventUtils.h"
#include "CsTrafficPrepattern.h"
#include "Traffic.h"
#include "TOpt.h"
#include "TWatches.h"
#include "TDisplay.h"
#include "TSetup.h"
#include "TEv.h"

using namespace std;

//Constructor
CsTrafficPrepattern::CsTrafficPrepattern() {
  new Traffic; // create Traffic package object
}

//Destructor (empty)
CsTrafficPrepattern::~CsTrafficPrepattern() {}

//Pattern recognition method
bool CsTrafficPrepattern::doPrepattern( const list<CsCluster*> &clusters, 
					const list<CsZone*> &zones) {

  if (TOpt::ReMode[0]>0) return(true); // ALL tracking is OFF
  if (TOpt::ReMode[2]>0) return(true); // Prepattern is OFF 

  Traffic::Ref().Stopwatch.Start(1);
#ifdef TraFDic_HISTO_STAT
  if (TOpt::Hist[1])  // So that TEv::Monitor be called and Stopwatch.Stop(8)
    Traffic::Ref().Stopwatch.Start(8);
#endif

  TSetup::Ptr()->Update();  //   *************** UPDATE SETUP ***************


  TEv *ev_p = TEv::Ptr(); bool reTrack = ev_p!=NULL; double bMSBackup = 0;
  if (reTrack) {           // ***** ALREADY EXISTING TEv: ReTRACKING *****
    if (ev_p->IsMC())// MC: backup BMS smearing value found in initial tracking.
      bMSBackup = ev_p->GetBMSSmearing();
    delete ev_p;
  }
  new TEv; TEv &ev = TEv::Ref();  // ***** INSTANTIATE TEv OBJECT *****
  if (reTrack) {
    ev.FlagReTracking();
    if (ev.IsMC())
      ev.SetBMSSmearing(bMSBackup); // Re-instate initial BMS smearing value.
  }

  if (TOpt::Print[0]&0x2) {   // ***** (OPTIONAL) PRINT EVENT NUMBER *****
    int ievt = ev.ptrEvt()->getEventNumberInRun();
    int jevt = ev.ptrEvt()->getEventNumberInBurst();
    printf("Evt %d %d\n",ievt,jevt);
  }
  bool trigger_ok = true;        // ***** TRIGGER SELECTION *****
  if (TOpt::iCut[0]) {
    const unsigned int allTrigs = 0xffff; // Although TCS can handle 0x7fffff, it was decided to limit max. #triggers to 16. And take advantage of the fact to make use of the other bits of the trigger TDC module.
    unsigned int evTrig = ev.TrigMask();
    evTrig &= allTrigs;  // Cut away trailing end bits (i.e. online filter...)
    trigger_ok =
      (evTrig&TOpt::iCut[0]) && // Require trigger to be included in selection...
      (evTrig&(~TOpt::iCut[0]))==0;  // ...AND strictly included
  }
  if (trigger_ok) {

    // *************** PrePattern IN SPECTROMETER'S ZONES ***************

    if (TOpt::ReMode[11]==2 &&
	(zones.size()==5 ||                     // Standard COMPASS configuration
	 zones.size()==6 && TOpt::ReMode[49])) {// Special Drell-Yan config. w/ a 6th zone for the vectex detector
      //      ********** TraFDic ON PHYSICS DATA **********
      //   ***** PROCESS BEAM TELESCOPE (==last zone) FIRST *****
      ev.ImportClusters(clusters,1); 

      ev.PrePattern2(TSetup::Ref().Group2Zone(4));
      //     ********** EVENT TIME ("TEv::eventTime/eventTRef") **********
      // - Get it from CsEvent or, by default, from beam-tracks's time if
      //  uniquely defined.                                      (SetEventTime)
      // - Disregard if smaller than "dCut[71]".                 (SetEventTime)
      // - If trigger matches "ReMode[35]", update drift CsClusters (i.e. the U
      //  coordinate of CsClusters from drift-like detectors).   (UpdateDrifts)
      //   If indeed, "eventTime" is copied to "eventTRef" and reset to 0.
      //  This availability of 2 different event times is not taken advantage of
      //  so far (as of 2008/02): might be in the future.
      // - The fact that the U coordinate of CsClusters has been updated means
      //  that the event time will be accounted for when it comes to correct
      //  drift hits for propagation time. The latter is done depending upon
      //  "TOpt::ReMode[36]".                           (FitSegments,TrackFit2)
      // - The event time may be available somewhat later in TraFDiC flowchart,
      //  viz. @ "TrackRefit" time. It may then be assigned to "TEv::eventTime",
      //  and taken into account when correcting drift hits for time
      //  propagation, depending upon "TOpt::ReMode[37]".         (TracksRefit)
      //  (Note: above matches are incl/exclusive depending upon "ReMode[40]")
      // - This last update of drift hits will be ported to CsClusters if
      //  "ReMode[29]&0x4".
      ev.SetEventTime();
      ev.UpdateDrifts();
      if (TOpt::ReMode[28]&0x2)  // Upon option, perform LR association on
	CsEventUtils::associateClusters();
      ev.ImportClusters(clusters,-1);
      if (TOpt::Print[0]&0x4) ev.DumpEvent();
 
      list<CsZone*>::const_iterator iZ;
      for (iZ = zones.begin(); iZ!=zones.end(); iZ++) { // ***** LOOP OVER ZONES
	// Do track segments reco in the requested zones...
	if (*iZ==TSetup::Ref().Group2Zone(4)) continue; // ...but beam telescope
	if (zones.size()==6 &&
	    *iZ==TSetup::Ref().Group2Zone(5)) continue; // ...and DY 6th zone
	ev.PrePattern2(*iZ);
      }
    } // End PR on physics data w/ TraFDic
    else {
      //  ********** OTHER KINDS of DATA, OTHER PR OPTIONS **********
      ev.ImportClusters(clusters,0);
      if (TOpt::Print[0]&0x4) ev.DumpEvent();
      list<CsZone*>::const_iterator iZ;
      for (iZ = zones.begin(); iZ!=zones.end(); iZ++) { // ***** LOOP OVER ZONES
	// Do track segments reco in the requested zones
	switch(TOpt::ReMode[11]) {
	case 0: ev.PrePattern (*iZ); break;  
	case 1: ev.PrePattern1(*iZ); break; // ... alternative PrePattern method
	case 2: ev.PrePattern2(*iZ); break;
	case 3:
	  if (TSetup::Ref().Group2Zone(0)==*iZ)
	    ev.PrePattern3(*iZ);   // Works only for 1st group
	  else ev.PrePattern(*iZ);
	  break;
	default: 
	  CsErrLog::msg(elFatal,__FILE__,__LINE__,
			"Unexpected ReMode[11] value = %d", TOpt::ReMode[11]);
	}
      }
    } // End PR
  } // End trigger OK

  Traffic::Ref().Stopwatch.Stop(1);

  return(true);
}

// Returns found track segments
bool CsTrafficPrepattern::getPatterns(list<CsTrack*>& tracks) {
  if(TOpt::ReMode[0] > 0) return(true); // do nothing
  if(TOpt::ReMode[1] < 2 ) TEv::Ref().ExportTracks(tracks);
  return(true);
}

//Accessor to unused clusters
list<CsCluster*> CsTrafficPrepattern::getUnusedClusters(void) const {
  list<CsCluster*> clusters;
  if(TOpt::ReMode[0] > 0) return(clusters); // do nothing
  if(TOpt::ReMode[1]!=2 ) clusters = TEv::Ref().ExportClusters("unused");
  return( clusters );
}
