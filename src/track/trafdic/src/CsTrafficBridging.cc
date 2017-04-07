// $Id: CsTrafficBridging.cc 14069 2015-09-17 20:44:46Z lsilva $

#include "CsTrafficBridging.h"
#include "CsErrLog.h"
#include "Traffic.h"
#include "TDisplay.h"
#include "TEv.h"
#include "TOpt.h"
#include "TWatches.h"

using namespace std;

//Constructor
CsTrafficBridging::CsTrafficBridging() {
  if(Traffic::Ptr() == NULL) new Traffic; // create Traffic package object 
}

//Destructor
CsTrafficBridging::~CsTrafficBridging() {}

//Bridging method
bool CsTrafficBridging::doBridging( list<CsTrack*>& tracks, list<CsCluster*>& clusters ) {
  if(TOpt::ReMode[0] > 0) return(true); // ALL tracking is OFF
  
  Traffic::Ref().Stopwatch.Start(2);

  if(TOpt::ReMode[1] == 0 && TEv::Ptr() != NULL) delete TEv::Ptr(); 
  if(TEv::Ptr() == NULL) new TEv; // create TEv object instance 
  TEv& ev = TEv::Ref();

  Traffic::Ref().Stopwatch.Start(13);
  if(TOpt::ReMode[1] < 2 ) {
    ev.ImportTracks(tracks, "mv");
    ev.ImportClusters(clusters);
  }
  Traffic::Ref().Stopwatch.Stop(13);
  

  // Fit track segments befor bridging
  ev.FitSegments();
  
  // Have the Pre-Pattern step of the ``step-by-step'' graphics
  // here: in order to benefit from "FitSegments".
  if(TOpt::Graph[0] > 0 && TOpt::Graph[6] > 0 ) {
    cout<<"Pre-Pattern is done"<<endl;
    TDisplay::Ref().Draw(3);  // draw tracks
    TDisplay::Ref().Draw(4);  // menu
  }

  // do the bridging through the magnets
  switch(TOpt::ReMode[12]) {
  case 0:
    ev.BridgeSegments (); // default
    break;
  case 1:
    ev.BridgeSegments1(); // ... alternative
    break;
  case 2:
    ev.BridgeSegments2();
    break;
  case 3:
    ev.BridgeSegments3();
    break;
  default:
    cout<<"CsTrafficBridging::doBridging() ==> unexpected ReMode[12] switch value : "
	<<TOpt::ReMode[12]<<endl;
    assert(false);
  }

  // do the bridging through the muon wall
  if(TOpt::ReMode[16] == 0) ev.BridgeMuons();


  if(TOpt::Graph[0] > 0 && TOpt::Graph[6] > 0 ) {
    cout<<"Bridging is done"<<endl;
    TDisplay::Ref().Draw(3);  // draw tracks
    TDisplay::Ref().Draw(4);  // draw menu
  }
  
  // export bridged tracks back to CORAL
  Traffic::Ref().Stopwatch.Start(14);
  if(TOpt::ReMode[1] < 2 ) {
    ev.ExportTracks(tracks);
    if (TOpt::ReMode[1]==0) clusters = ev.ExportClusters("unused");
  }
  Traffic::Ref().Stopwatch.Stop(14);

  Traffic::Ref().Stopwatch.Stop(2);

  return( true );
}










