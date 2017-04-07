#include "GroupMumega.h"
#include "TThread.h"
#include "GroupMumegaPanel.h"
#include "ChipF1.h"

#include "GeomPlane.h"

ClassImp(GroupMumega);

const int GroupMumega::fGATE = 300;

#define TOTMIN  850
#define MULTMAX  10
#define DIFFXMAX 10


GroupMumega::GroupMumega(const char* name)
  : Group(name) {}


void GroupMumega::Init() {

#if USE_TRACK == 1
  fTrackMaker = new TrackMaker(GetName(),2,0.1);
  if (fDoTracking) {
    for(unsigned i=0; i<fPlanes.size(); i++) {
      const GeomPlane *geom = fPlanes[i]->GetGeometry();
      int id = geom->GetID();
      const char *name = geom->GetName();
      int nwires = geom->GetNWires();
      float x = geom->GetX();
      float y = geom->GetY();
      float z = geom->GetZ();
      float dx = geom->GetDX();
      float dy = geom->GetDY();
      float dz = geom->GetDZ();
      float angle = geom->GetAngle();
      float inpitch = .036;
      float outpitch = .042;
      float ctmin = -500;
      float ctmax =  500;
      float ctotmin = 600;
      fTrackMaker->AddTracker(id,name,nwires,x,y,z,dx,dy,dz,angle,inpitch,outpitch,
			      ctmin,ctmax,ctotmin);
    }
  }

  fTrackMaker->Init();

  // Adds tracking histogram to this's histogram list
  std::vector<TH1*>& trackhists=fTrackMaker->GetHistograms();
  for(unsigned i=0; i<trackhists.size();i++)
    fHistList.push_back(trackhists[i]);
#endif
}


void GroupMumega::ResetHistograms() {
  Group::ResetHistograms();

#if USE_TRACK == 1
  if (thr_flag) TThread::Lock();
  fTrackMaker->ResetHistos();
  if (thr_flag) TThread::UnLock();
#endif
}


void GroupMumega::Tracking() {

#if USE_TRACK == 1
  if (thr_flag) TThread::Lock();

  // reset tracking
  fTrackMaker->Reset();

  // passing clusters to tracking
  for(unsigned i=0; i<fPlanes.size(); i++) {
    if(const GeomPlane *geom = fPlanes[i]->GetGeometry()) {
      const std::set<CCluster1>& clusters = geom->GetClusters();
      for(std::set<CCluster1>::iterator ic=clusters.begin(); ic!=clusters.end(); ic++)
	fTrackMaker->AddClusterWRS(geom->GetID(), ic->pos, ic->size, ic->res);
    }
  }

  for(unsigned i=0; i<fPlanes.size(); i++) {
    if(const GeomPlane *geom = fPlanes[i]->GetGeometry()) {
      Tracker *ct=fTrackMaker->GetTracker(geom->GetID());
      if(ct) ct->Activate(0);
      else {
	std::cerr<<"Cannot find Tracker with ID "<<geom->GetID()<<std::endl;
	exit(1);
      }
      fTrackMaker->Clusterize(1);
      fTrackMaker->Candidates();
      if(ct) ct->Activate(1);
    }
  }
  if (thr_flag) TThread::UnLock();
#endif
}


void GroupMumega::ControlPanel(const TGWindow *p, const TGWindow *main) {

  if (!fControlPanel) fControlPanel = new GroupMumegaPanel(p, main, 100, 100, this);
}













