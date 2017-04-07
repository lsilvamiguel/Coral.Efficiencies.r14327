// $Id: CsTrafficFitting.cc 13148 2011-12-28 16:55:25Z kbicker $

#include "CsTrafficFitting.h"
#include "Traffic.h"
#include "TEv.h"
#include "TOpt.h"
#include "TWatches.h"

using namespace std;

// Constructor
CsTrafficFitting::CsTrafficFitting() {
  if(Traffic::Ptr() == NULL) new Traffic; // create Traffic package object 
}

//Destructor (empty)
CsTrafficFitting::~CsTrafficFitting() {}

// Fitting method
bool CsTrafficFitting::doFitting( list<CsTrack*>& tracks, const list<CsCluster*>& clusters) {
  if(TOpt::ReMode[0] > 0) return(true); // Traffic is OFF


  if(TOpt::ReMode[1] == 0 && TEv::Ptr() != NULL) delete TEv::Ptr();
  if(TEv::Ptr() == NULL) new TEv; // create TEv object instance
  TEv& ev = TEv::Ref();

  Traffic::Ref().Stopwatch.Start(3);

  Traffic::Ref().Stopwatch.Start(15);
  if(TOpt::ReMode[1] < 2 ) {
    ev.ImportTracks(tracks, "mv");
    ev.ImportClusters(clusters);
  }
  Traffic::Ref().Stopwatch.Stop(15);


  // do the fit (with possible clusters re-assignment)
  switch(TOpt::ReMode[13]){
  case 0:
    ev.TracksFit ();  // default
    break;
  case 1:
    ev.TracksFit1();  // ... alternative
    break;
  case 2:
    ev.TracksFit2();
    break;
  case 3:
    ev.TracksFit3();
    break;
  default:
    cout<<"CsTrafficFitting::doFitting ==> unexpected ReMode[13] switch value : "
	<<TOpt::ReMode[13]<<endl;
    assert(false);
  }

  // save fitted tracks
  ev.ExportTracks(tracks);

  Traffic::Ref().Stopwatch.Stop(3);

  if(TOpt::Graph[0] > 0 && TOpt::Graph[6] > 0 ) {
    //TDisplay::Ref().Draw(3);  // draw fited tracks
    cout<<"Fitting is done"<<endl;
  }

  return( true );
}
















